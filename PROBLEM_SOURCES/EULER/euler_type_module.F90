MODULE euler_type_MODULE
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE def_type_periodic
   USE euler_bc_arrays
   USE Butcher_tableau
   USE euler_matrices_module
   USE space_dim
   IMPLICIT NONE

   ABSTRACT INTERFACE
      FUNCTION function_template_pressure(rho, ie) RESULT(vv)
         REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, ie
         REAL(KIND = 8), DIMENSION(SIZE(rho, 1)) :: vv
      END FUNCTION function_template_pressure
   END INTERFACE

   ABSTRACT INTERFACE
      SUBROUTINE function_template_impose_bc(un, euler_bc, mesh, time)
         USE def_type_mesh
         USE euler_bc_arrays
         REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
         TYPE(mesh_type) :: mesh
         TYPE(euler_bc_type) :: euler_bc
         REAL(KIND = 8) :: time
      END SUBROUTINE function_template_impose_bc
   END INTERFACE

   TYPE euler_type
      MPI_Comm :: communicator
      TYPE(mesh_type), POINTER :: mesh
      TYPE(petsc_csr_LA), POINTER :: LA
      TYPE(periodic_type), POINTER :: per
      PROCEDURE(function_template_pressure), NOPASS, POINTER :: pressure
      PROCEDURE(function_template_impose_bc), NOPASS, POINTER :: impose_bc
      TYPE(BT), PUBLIC :: ERK
      TYPE(euler_bc_type) :: euler_bc
      TYPE(euler_matrices_type) :: matrices
      REAL(KIND = 8) :: dt, time, in_tol
      LOGICAL :: no_iter
      INTEGER :: syst_dim = k_dim + 2

      Vec, PRIVATE :: x1vec, x2vec, x3vec, x2_ghost, vec_loc
      INTEGER, DIMENSION(:), POINTER :: tab
   CONTAINS
      PROCEDURE, PUBLIC :: init => init_euler
      PROCEDURE, PUBLIC :: update
      PROCEDURE, PRIVATE :: compute_dij

   END TYPE euler_type

CONTAINS
   SUBROUTINE init_euler(this, communicator, mesh, LA, per, pressure, erk_sv, impose_bc, time_init)
      USE st_matrix
      CLASS(euler_type), INTENT(INOUT) :: this
      MPI_Comm, INTENT(IN) :: communicator
      TYPE(mesh_type), TARGET, INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
      TYPE(periodic_type), TARGET, INTENT(IN) :: per
      INTEGER :: erk_sv, ierr, n
      REAL(KIND = 8) :: time_init
      PROCEDURE(function_template_pressure) :: pressure
      PROCEDURE(function_template_impose_bc) :: impose_bc
      INTEGER, POINTER, DIMENSION(:) :: ifrom

      this%mesh => mesh
      this%communicator = communicator
      this%LA => LA
      this%per => per
      this%pressure => pressure
      this%impose_bc => impose_bc
      this%euler_bc%syst_dim = this%syst_dim
      this%time = time_init

      this%in_tol = 1.d-2
      this%no_iter = .true.

      CALL this%ERK%init(erk_sv)
      CALL this%euler_bc%construct_euler_bc(this%mesh)
      CALL this%matrices%construct(this%communicator, this%mesh, this%LA)

      CALL create_my_ghost(this%mesh, this%LA, ifrom)
      CALL VecCreateGhost(this%communicator, this%mesh%dom_np, &
           PETSC_DETERMINE, SIZE(ifrom), ifrom, this%x1vec, ierr)

      CALL VecDuplicate(this%x1vec, this%x2vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x3vec, ierr)

      CALL VecGhostGetLocalForm(this%x2vec, this%x2_ghost, ierr)

      CALL VecCreateSeq(PETSC_COMM_SELF, this%mesh%dom_np, this%vec_loc, ierr)
      ALLOCATE(this%tab(this%mesh%dom_np))
      DO n = 1, this%mesh%dom_np
         this%tab(n) = n - 1
      END DO

   END SUBROUTINE init_euler

   SUBROUTINE update(this, un)
      USE petsc_tools
      USE euler_flux
      USE st_matrix
      CLASS(euler_type) :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim) :: ff
      REAL(KIND = 8), DIMENSION(this%mesh%np) :: rk
      REAL(KIND = 8), DIMENSION(this%mesh%dom_np) :: dij_diag
      REAL(KIND = 8) :: dt_min_loc, dt_min_glob
      INTEGER :: k, comp, ierr

      !===compute dij
      CALL this%compute_dij(un)

      CALL MatGetDiagonal(this%matrices%dij, this%vec_loc, ierr)
      CALL VecGetValues(this%vec_loc, this%mesh%dom_np, this%tab, dij_diag, ierr)


      dij_diag = this%matrices%lumped_mass(1:this%mesh%dom_np) / ABS(dij_diag)
      dt_min_loc = MINVAL(dij_diag)


      CALL MPI_ALLREDUCE(dt_min_loc, dt_min_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PETSC_COMM_WORLD, ierr)
      this%dt = dt_min_glob
            IF (this%mesh%rank == 0) write(*, *) this%time, this%dt

      this%dt = 2.d-1 / real(SUM(this%mesh%domnp))
      this%time = this%time + this%dt

      DO comp = 1, this%syst_dim
         ff = flux(comp, un)

         CALL VecSet(this%x3vec, 0.d0, ierr)
         DO k = 1, k_dim
            !=== set flux_k in x1vec
            CALL array_to_petsc_vec(ff(:, k), this%x1vec, this%mesh, this%LA, 'insert')
            !=== compute sum_j cij_k * fluxj_k in x2vec
            CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
            !=== construct sum_k sum_j cij_k flux_k into x3vec
            CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
         END DO

         !=== set un(comp) in x1vec
         CALL array_to_petsc_vec(un(:, comp), this%x1vec, this%mesh, this%LA, 'insert')

         !=== add dij un(comp)to x3vec in x2vec
         CALL MatMultAdd(this%matrices%dij, this%x1vec, this%x3vec, this%x2vec, ierr)

         CALL VecGhostGetLocalForm(this%x2vec, this%x2_ghost, ierr)
         CALL VecGhostUpdateBegin(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL extract(this%x2_ghost, 1, 1, this%LA, rk)

         rk = rk * this%dt / this%matrices%lumped_mass

         un(:, comp) = un(:, comp) + rk

         CALL this%impose_bc(un, this%euler_bc, this%mesh, this%time)
      END DO

   END SUBROUTINE update


   SUBROUTINE compute_dij(this, un)
      USE space_dim
      USE petsc
      USE my_util
      USE def_type_mesh
      USE arbitrary_eos_lambda_module
      IMPLICIT NONE
      CLASS(euler_type) :: this
      TYPE(mesh_type), POINTER :: mesh
      TYPE(petsc_csr_LA), POINTER :: LA
      REAL(KIND = 8), DIMENSION(:, :) :: un
      INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge
      INTEGER, DIMENSION(1) :: i_t, j_t, idx, jdx
      REAL(KIND = 8), DIMENSION(1, k_dim) :: nij_c
      REAL(KIND = 8), DIMENSION(1) :: norm_c, dij_c
      REAL(KIND = 8), DIMENSION(2) :: u, rho, ie, p, lambda_max
      REAL(KIND = 8) :: pstar
      LOGICAL, DIMENSION(this%mesh%medge) :: virgin_edge

      CALL MatZeroEntries(this%matrices%dij, ierr)

      mesh => this%mesh
      LA => this%LA

      virgin_edge = .true.
      nw = mesh%gauss%n_w

      DO m = 1, mesh%me
         DO n = 1, mesh%gauss%n_e
            IF (mesh%attr_e(mesh%jce(n, m))) THEN
               edge = mesh%jce_loc(n, m)
               IF (.not. virgin_edge(edge)) CYCLE
               virgin_edge(edge) = .false.

               ni = MOD(n, nw) + 1
               nj = MOD(n + 1, nw) + 1
               i = mesh%jj(ni, m)
               j = mesh%jj(nj, m)
               i_t = i
               j_t = j

               DO k = 1, k_dim
                  CALL MatGetValues(this%matrices%nij_loc(k), 1, i_t - 1, 1, j_t - 1, nij_c(:, k), ierr)
               END DO

               rho(1) = un(i, 1)
               rho(2) = un(j, 1)

               u(1) = SUM(un(i, 2:1 + k_dim) * nij_c(1, :)) / rho(1)
               u(2) = SUM(un(j, 2:1 + k_dim) * nij_c(1, :)) / rho(2)

               ie(1) = un(i, k_dim + 2) / rho(1) - 0.5d0 * u(1) * u(1)
               ie(2) = un(j, k_dim + 2) / rho(2) - 0.5d0 * u(2) * u(2)

               p = this%pressure(rho, ie)

               CALL lambda_arbitrary_eos(rho, u, ie, p, this%in_tol, this%no_iter, lambda_max, pstar)

               CALL MatGetValues(this%matrices%cij_norm_loc, 1, i_t - 1, 1, j_t - 1, norm_c, ierr)

               dij_c = MAXVAL(lambda_max) * norm_c

               IF (mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j

                  DO k = 1, k_dim
                     CALL MatGetValues(this%matrices%nij_loc(k), 1, j_t - 1, 1, i_t - 1, nij_c(:, k), ierr)
                  END DO

                  u(1) = SUM(un(i, 2:1 + k_dim) * nij_c(1, :)) / rho(1)
                  u(2) = SUM(un(j, 2:1 + k_dim) * nij_c(1, :)) / rho(2)

                  rho = (/rho(2), rho(1)/)
                  ie = (/ie(2), ie(1)/)
                  p = (/p(2), p(1)/)

                  CALL lambda_arbitrary_eos(rho, u, ie, p, this%in_tol, this%no_iter, lambda_max, pstar)

                  dij_c = MAX(dij_c, MAXVAL(lambda_max) * norm_c)

               END IF

               idx = LA%loc_to_glob(1, i) - 1
               jdx = LA%loc_to_glob(1, j) - 1

               CALL MatSetValues(this%matrices%dij, 1, idx, 1, jdx, dij_c, ADD_VALUES, ierr)
               CALL MatSetValues(this%matrices%dij, 1, jdx, 1, idx, dij_c, ADD_VALUES, ierr)
               CALL MatSetValues(this%matrices%dij, 1, idx, 1, idx, -dij_c, ADD_VALUES, ierr)
               CALL MatSetValues(this%matrices%dij, 1, jdx, 1, jdx, -dij_c, ADD_VALUES, ierr)
            END IF

         END DO

      END DO

      CALL MatAssemblyBegin(this%matrices%dij, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd  (this%matrices%dij, MAT_FINAL_ASSEMBLY, ierr)

   END SUBROUTINE compute_dij


END MODULE euler_type_MODULE
