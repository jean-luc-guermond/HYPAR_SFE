MODULE euler_type_MODULE
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE euler_bc_arrays
   USE Butcher_tableau
   USE euler_matrices_module
   USE periodic_data_module
   USE cell_limiting_engine_parallel_module
   USE read_inputs_module
   IMPLICIT NONE

   ABSTRACT INTERFACE
      FUNCTION function_template_pressure(rho, ie) RESULT(vv)
         REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, ie
         REAL(KIND = 8), DIMENSION(SIZE(rho, 1))  :: vv
      END FUNCTION function_template_pressure
   END INTERFACE

   ABSTRACT INTERFACE
      SUBROUTINE function_template_impose_bc(un, euler_bc, mesh, time)
         USE def_type_mesh
         USE euler_bc_arrays
         REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
         TYPE(mesh_type)     :: mesh
         TYPE(euler_bc_type) :: euler_bc
         REAL(KIND = 8)      :: time
      END SUBROUTINE function_template_impose_bc
   END INTERFACE

   TYPE argument_euler_type
      CHARACTER(LEN=rec_length) :: CFL       = '=== CFL ? ==='
      CHARACTER(LEN=rec_length) :: method = '=== Which method to solve Euler (viscous, high) ? ==='
      CHARACTER(LEN=rec_length) :: erk_sv
      CHARACTER(LEN=rec_length) :: eos_param
   END TYPE argument_euler_type

   TYPE euler_type
      !===Parameters read from data
      REAL(KIND=8)                 :: CFL       = 0.5d0
      CHARACTER(LEN=rec_length)    :: method = 'viscous'
      REAL(KIND = 8), DIMENSION(1) :: eos_param = 0.d0
      INTEGER                      :: erk_sv    = -21
      !===Parameters built along way
      MPI_Comm :: communicator
      Vec :: x1vec, x2vec, x3vec, x4vec, x5vec, x2_ghost, vec_loc
      CHARACTER(100) :: name
      LOGICAL :: no_iter
      INTEGER :: syst_dim
      REAL(KIND = 8) :: dt, time, final_time, in_tol
      INTEGER, DIMENSION(:), POINTER :: tab
      TYPE(mesh_type),     POINTER :: mesh
      TYPE(petsc_csr_LA),  POINTER :: LA
      TYPE(BT),             PUBLIC :: ERK
      TYPE(euler_bc_type)          :: euler_bc
      TYPE(euler_matrices_type)    :: matrices
      TYPE(limiting_type)          :: limiting
      PROCEDURE(function_template_pressure),  NOPASS, POINTER :: pressure
      PROCEDURE(function_template_impose_bc), NOPASS, POINTER :: impose_bc
   CONTAINS
      PROCEDURE, PUBLIC  :: init => init_euler
      PROCEDURE, PUBLIC  :: read_euler_data
      PROCEDURE, PUBLIC  :: update
      PROCEDURE, PRIVATE :: compute_dij, compute_dt
      PROCEDURE, PRIVATE :: flux
      PROCEDURE, PRIVATE :: compute_dK, compute_dt_from_dK
   END TYPE euler_type

CONTAINS
   SUBROUTINE init_euler(this, communicator, name, mesh, LA, pressure, impose_bc, times)
   ! SUBROUTINE init_euler(this, communicator, name, mesh, LA, per, pressure, impose_bc, times)
      USE space_dim
      USE st_matrix
      USE periodic_data_module
      IMPLICIT NONE
      CLASS(euler_type), INTENT(INOUT) :: this
      MPI_Comm, INTENT(IN) :: communicator
      CHARACTER(100) :: name
      TYPE(mesh_type), TARGET, INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
      ! TYPE(periodic_type), TARGET, INTENT(IN) :: per
      INTEGER :: ierr, n
      REAL(KIND = 8), DIMENSION(2) :: times
      PROCEDURE(function_template_pressure) :: pressure
      PROCEDURE(function_template_impose_bc) :: impose_bc
      INTEGER, POINTER, DIMENSION(:) :: ifrom

      this%syst_dim = k_dim + 2

      this%name = name
      this%mesh => mesh
      this%communicator = communicator
      this%LA => LA
      this%pressure => pressure
      this%impose_bc => impose_bc
      this%euler_bc%syst_dim = this%syst_dim
      this%time = times(1) !<==initial_time
      this%final_time = times(2) !<==final_time

      !===Parameters for lambda_arbitrary_eos
      this%in_tol = 1.d-2
      this%no_iter = .TRUE.

      CALL this%read_euler_data("EULER PARAMETERS")

      !=== new Butcher module
      this%ERK%sv = this%erk_sv
      CALL this%ERK%init()
      !=== end new Butcher module

      this%matrices%method = this%method !<==transfer this%method to this%matrices

      !===Boundary conditions
      CALL this%euler_bc%construct_euler_bc(this%mesh)

      !===Periodic boundary if any
      CALL this%matrices%construct(this%communicator, this%mesh, this%LA)

      !===Goshting structures
      CALL create_my_ghost(this%mesh, this%LA, ifrom)
      CALL VecCreateGhost(this%communicator, this%mesh%dom_np, &
           PETSC_DETERMINE, SIZE(ifrom), ifrom, this%x1vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x2vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x3vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x4vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x5vec, ierr)
      CALL VecGhostGetLocalForm(this%x2vec, this%x2_ghost, ierr)

      CALL VecCreateSeq(PETSC_COMM_SELF, this%mesh%dom_np, this%vec_loc, ierr)
      ALLOCATE(this%tab(this%mesh%dom_np))
      DO n = 1, this%mesh%dom_np
         this%tab(n) = n - 1
      END DO

      !===Limiting
      CALL this%limiting%init(this%communicator, this%name, this%mesh, this%LA)

   END SUBROUTINE init_euler

   SUBROUTINE read_euler_data(this, section_name)
     IMPLICIT NONE
     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name

     CLASS(euler_type), INTENT(INOUT) :: this
     TYPE(argument_euler_type)        :: argument_data


     !===Initialize data arguments (depends on the name)
     argument_data%eos_param = '===' // TRIM(ADJUSTL(this%name)) // ': b_covolume ?==='
     argument_data%erk_sv = '===' // TRIM(ADJUSTL(this%name)) // ': ERK ?==='

     !================
     !=== MANDATORY Reading all data file
     !================
     IF (PRESENT(section_name)) THEN
        CALL read_data_init_list(section_name)
     ELSE
        CALL read_data_init_list()
     END IF

     !================
     !=== We now find the relevant information for this specific Euler data
     !================
     !===CFL
     CALL read_data(argument_data%CFL, this%CFL)

     !===b_covolume
     CALL read_data(argument_data%eos_param, this%eos_param(1))

     !===ERK
     CALL read_data(argument_data%erk_sv, this%erk_sv)

     !===Method order
     CALL read_data(argument_data%method, this%method)

     !================
     !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
     !================

     CALL finalize_rewrite_data
   END SUBROUTINE read_euler_data

   SUBROUTINE update(this, un)
     USE space_dim
     USE petsc_tools
     USE st_matrix
     USE my_util
     USE cell_limiting_engine_module
     USE sub_plot
     CLASS(euler_type)                                                     :: this
     REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
     REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim)                :: un_temp
     REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim) :: ff
     REAL(KIND = 8), DIMENSION(this%mesh%np)                       :: rk
     REAL(KIND = 8), DIMENSION(this%mesh%np,2)                     :: bounds
     INTEGER :: comp, k, ierr, it, code
     REAL(KIND = 8) :: norm
     INTEGER, PARAMETER :: limit_max = 2 !<<FIXME
     un_temp = un
     SELECT CASE(this%method)
     CASE('viscous')
        !===compute dijL and dt
        CALL this%compute_dij(un_temp, bounds)
        CALL this%compute_dt

        this%time = this%time + this%dt

        DO comp = 1, this%syst_dim
           ff = 0.d0
           ff = this%flux(comp, un_temp)

           CALL VecSet(this%x3vec, 0.d0, ierr)
           DO k = 1, k_dim
              !=== set flux_k in x1vec
              CALL array_to_petsc_vec(ff(:, k), this%x1vec, this%mesh, this%LA, 'insert')
              !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
              CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
              !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
              CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
           END DO

           !=== set un(comp) in x1vec
           CALL array_to_petsc_vec(un_temp(:, comp), this%x1vec, this%mesh, this%LA, 'insert')
           !=== add dij un(comp)to x3vec in x2vec
           CALL MatMultAdd(this%matrices%dijL, this%x1vec, this%x3vec, this%x2vec, ierr)
           CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x2vec, this%LA)
           CALL VecGhostGetLocalForm(this%x2vec, this%x2_ghost, ierr)
           CALL VecGhostUpdateBegin(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
           CALL VecGhostUpdateEnd(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
           CALL extract(this%x2_ghost, 1, 1, this%LA, rk)

           rk = rk * this%dt / this%matrices%lumped_mass

           un(:, comp) = un_temp(:, comp) + rk

           DO k = 1, this%mesh%per%nb_bords
              un(this%mesh%per%list(k)%DIL, comp) = un(this%mesh%per%perlist(k)%DIL, comp)
           END DO

           CALL this%impose_bc(un, this%euler_bc, this%mesh, this%time)
        END DO
     CASE('high')
        !===compute dijL and dt
        CALL this%compute_dij(un_temp, bounds)
        CALL this%compute_dt

        this%time = this%time + this%dt

        DO comp = 1, this%syst_dim
           ff = this%flux(comp, un_temp)

           CALL VecSet(this%x3vec, 0.d0, ierr)
           DO k = 1, k_dim
              !=== set flux_k in x1vec
              CALL array_to_petsc_vec(ff(:, k), this%x1vec, this%mesh, this%LA, 'insert')
              !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
              CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
              !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
              CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
           END DO

           !=== set un(comp) in x1vec
           CALL array_to_petsc_vec(un_temp(:, comp), this%x1vec, this%mesh, this%LA, 'insert')
           !TEST Low order
!!$           IF (comp==1) THEN
!!$              CALL MatMultAdd(this%matrices%dijL, this%x1vec, this%x3vec, this%x4vec, ierr)
!!$              CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x4vec, this%LA)
!!$              CALL VecGhostGetLocalForm(this%x4vec, this%x2_ghost, ierr)
!!$              CALL VecGhostUpdateBegin(this%x4vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
!!$              CALL VecGhostUpdateEnd(this%x4vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
!!$              CALL extract(this%x2_ghost, 1, 1, this%LA, rk)
!!$              rk = rk * this%dt / this%matrices%lumped_mass
!!$              un(:, comp) = un_temp(:, comp) + rk
!!$              bounds(:,2) = un(:, 1)
!!$           END IF
           !END TEST
           !=== add dij un(comp)to x3vec in x2vec
           CALL MatMultAdd(this%matrices%dijH, this%x1vec, this%x3vec, this%x2vec, ierr)
           CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x2vec, this%LA)
           CALL VecGhostGetLocalForm(this%x2vec, this%x2_ghost, ierr)
           CALL VecGhostUpdateBegin(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
           CALL VecGhostUpdateEnd(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
           CALL extract(this%x2_ghost, 1, 1, this%LA, rk)

           !rk = rk * this%dt / this%matrices%lumped_mass
           !===FIXME
           rk = rk*this%dt
           CALL divide_by_mass(this,rk)
           !===FIXME

           un(:, comp) = un_temp(:, comp) + rk
        END DO
        !===Limiting
        IF (this%limiting%if_limiting) THEN
           DO it = 1, limit_max
               CALL this%limiting%iterative_cell_limiting_procedure(un,bounds(:,1),&
                   psi_rho_min,zero_of_psi_rho_min,un_temp)
               un(:,:) = un_temp(:,:)
           END DO
           DO it = 1, limit_max
               CALL this%limiting%iterative_cell_limiting_procedure(un,bounds(:,2),&
                   psi_rho_max,zero_of_psi_rho_max,un_temp)
               un(:,:) = un_temp(:,:)
           END DO
        END IF
        !===Periodicity
        DO comp = 1, this%syst_dim 
           DO k = 1, this%mesh%per%nb_bords
              un(this%mesh%per%list(k)%DIL, comp) = un(this%mesh%per%perlist(k)%DIL, comp)
           END DO
        END DO

        !===Boundary conditions
        CALL this%impose_bc(un, this%euler_bc, this%mesh, this%time)
     CASE DEFAULT
        CALL error_petsc('Wrong method in euler update')
     END SELECT

   CONTAINS
     FUNCTION psi_rho_min(x,psi_m) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:) :: x
       REAL(KIND=8) :: psi_m, v
       v = x(1)-psi_m
     END FUNCTION psi_rho_min

     FUNCTION zero_of_psi_rho_min(psi_m,u0,P) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:) :: u0, P
       REAL(KIND=8) :: psi_m, v
       v = (psi_m-u0(1))/P(1)
     END FUNCTION zero_of_psi_rho_min

     FUNCTION psi_rho_max(x,psi_m) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:) :: x
       REAL(KIND=8) :: psi_m, v
       v = psi_m-x(1)
     END FUNCTION psi_rho_max

     FUNCTION zero_of_psi_rho_max(psi_m,u0,P) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:) :: u0, P
       REAL(KIND=8) :: psi_m, v
       v = (psi_m-u0(1))/P(1)
     END FUNCTION zero_of_psi_rho_max

   END SUBROUTINE update

   FUNCTION flux(this, comp, un) RESULT(vv)  
      USE space_dim
      USE my_util
      IMPLICIT NONE
      CLASS(euler_type)                           :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
      INTEGER,                         INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(SIZE(un, 1), k_dim) :: vv
      REAL(KIND = 8), DIMENSION(SIZE(un, 1))                       :: H, u, ie
      INTEGER :: k

      SELECT CASE(comp)
      CASE(1)
         DO k = 1, k_dim
            vv(:, k) = un(:, k + 1)
         END DO
      CASE(2:k_dim + 1)
         u = un(:, comp) / un(:, 1)
         DO k = 1, k_dim
            vv(:, k) = un(:, k + 1) * u
         END DO
         ie = un(:, k_dim + 2) / un(:, 1)
         DO k = 1, k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO
         vv(:, comp - 1) = vv(:, comp - 1) + this%pressure(un(:, 1), ie)
      CASE(k_dim + 2)
         ie = un(:, k_dim + 2) / un(:, 1)
         DO k = 1, k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO

         H = un(:, comp) + this%pressure(un(:, 1), ie)
         DO k = 1, k_dim
            vv(:, k) = (un(:, k + 1) / un(:, 1)) * H
         END DO
      CASE DEFAULT
         CALL error_petsc(' BUG in flux, wrong comp = '//itoa(comp)//" with k_dim="//itoa(k_dim))
      END SELECT
    END FUNCTION flux

    FUNCTION pressure_from_state(this, un) RESULT(vv)  
      USE space_dim
      IMPLICIT NONE
      CLASS(euler_type)                           :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
      REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: e, vv
      INTEGER :: k
      e = 0.d0
      DO k = 1, k_dim
         e = e + un(:,k+1)**2
      END DO
      e = un(:,k_dim+2)/un(:,1) - 0.5d0*e/un(:,1)**2
      vv = this%pressure(un(:,1),e)
    END FUNCTION pressure_from_state

!========================================================
!========== PRIVATE PROCEDURES ==========================
!========================================================


   SUBROUTINE compute_dt(this)
      IMPLICIT NONE
      CLASS(euler_type) :: this
      REAL(KIND = 8), DIMENSION(this%mesh%dom_np) :: dijL_diag
      REAL(KIND = 8) :: dt_min_loc, dt_min_glob
      INTEGER :: ierr

      CALL MatGetDiagonal(this%matrices%dijL, this%vec_loc, ierr)
      CALL VecGetValues(this%vec_loc, this%mesh%dom_np, this%tab, dijL_diag, ierr)
      dijL_diag = this%matrices%lumped_mass(1:this%mesh%dom_np) / ABS(dijL_diag)

      dt_min_loc = MINVAL(dijL_diag) / 2.d0

      CALL MPI_ALLREDUCE(dt_min_loc, dt_min_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PETSC_COMM_WORLD, ierr)
      this%dt = this%CFL * dt_min_glob
   END SUBROUTINE compute_dt

   SUBROUTINE compute_dij(this, un, bounds)
     USE space_dim
     USE petsc
     USE my_util
     USE def_type_mesh
     USE arbitrary_eos_lambda_module
     USE compute_periodic
     USE petsc_tools, ONLY : array_to_petsc_vec
     IMPLICIT NONE
     CLASS(euler_type) :: this
     TYPE(mesh_type), POINTER :: mesh
     TYPE(petsc_csr_LA), POINTER :: LA
     REAL(KIND = 8), DIMENSION(:, :) :: un, bounds
     INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge, comp, rank
     INTEGER, DIMENSION(1) :: i_t, j_t, idx, jdx
     REAL(KIND = 8), DIMENSION(1, k_dim) :: nij_c
     REAL(KIND = 8), DIMENSION(1) :: norm_c, dijL_c
     REAL(KIND = 8), DIMENSION(1) :: dijH_c, extrem_value
     REAL(KIND = 8), DIMENSION(2) :: u, rho, ie, p, lambda_max
     REAL(KIND = 8) :: pstar, max_lambda, uijbar, norm_petsc
     LOGICAL, DIMENSION(this%mesh%medge) :: virgin_edge
     REAL(KIND = 8), DIMENSION(this%mesh%np)  :: alpha !<==commutator in (0,1)

     !===Compute commutator if needed
     IF (this%method=='high') THEN
        CALL commutator(this, un, alpha)
     END IF

     !===Compute dijL
     CALL MatZeroEntries(this%matrices%dijL, ierr)
     IF (this%method=='high') THEN
        CALL MatZeroEntries(this%matrices%dijH, ierr)
     END IF
     mesh => this%mesh
     LA => this%LA

     virgin_edge = .TRUE.
     nw = mesh%gauss%n_w
     bounds(:,1) = un(:,1)
     bounds(:,2) = un(:,1)
     DO m = 1, mesh%me
        DO n = 1, mesh%gauss%n_e
           IF (mesh%attr_e(mesh%jce(n, m))) THEN
              edge = mesh%jce_loc(n, m)
              IF (.NOT. virgin_edge(edge)) CYCLE
              virgin_edge(edge) = .FALSE.

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

              CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, lambda_max, pstar)
              CALL MatGetValues(this%matrices%cij_norm_loc, 1, i_t - 1, 1, j_t - 1, norm_c, ierr)

              max_lambda = MAXVAL(lambda_max)
              dijL_c = max_lambda * norm_c

              IF (mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j

                 DO k = 1, k_dim
                    CALL MatGetValues(this%matrices%nij_loc(k), 1, j_t - 1, 1, i_t - 1, nij_c(:, k), ierr)
                 END DO

                 u(1) = SUM(un(i, 2:1 + k_dim) * nij_c(1, :)) / rho(1)
                 u(2) = SUM(un(j, 2:1 + k_dim) * nij_c(1, :)) / rho(2)

                 rho = (/rho(2), rho(1)/)
                 ie = (/ie(2), ie(1)/)
                 p = (/p(2), p(1)/)

                 CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, lambda_max, pstar)


                 dijL_c = MAX(dijL_c, MAXVAL(lambda_max) * norm_c)
                 max_lambda = MAX(max_lambda,MAXVAL(lambda_max))

              END IF

              idx = LA%loc_to_glob(1, i) - 1
              jdx = LA%loc_to_glob(1, j) - 1

              CALL MatSetValues(this%matrices%dijL, 1, idx, 1, jdx, dijL_c, ADD_VALUES, ierr)
              CALL MatSetValues(this%matrices%dijL, 1, jdx, 1, idx, dijL_c, ADD_VALUES, ierr)

              CALL MatSetValues(this%matrices%dijL, 1, idx, 1, idx, -dijL_c, ADD_VALUES, ierr) !===add value on diagonal
              CALL MatSetValues(this%matrices%dijL, 1, jdx, 1, jdx, -dijL_c, ADD_VALUES, ierr) !===add value on diagonal
              IF (this%method=='high') THEN
                 dijH_c = dijL_c*(alpha(i)+alpha(j))/2
                 CALL MatSetValues(this%matrices%dijH, 1, idx, 1, jdx, dijH_c, ADD_VALUES, ierr)
                 CALL MatSetValues(this%matrices%dijH, 1, jdx, 1, idx, dijH_c, ADD_VALUES, ierr)
                 CALL MatSetValues(this%matrices%dijH, 1, idx, 1, idx, -dijH_c, ADD_VALUES, ierr) !===add value on diagonal
                 CALL MatSetValues(this%matrices%dijH, 1, jdx, 1, jdx, -dijH_c, ADD_VALUES, ierr) !===add value on diagonal
                 !===Compute low-order update to estimate bounds
                 uijbar = (un(i, 1)+un(j, 1))/2 &
                      - SUM((un(j, 2:k_dim+1) - un(i, 2:k_dim+1))*nij_c(1, :))/(2*max_lambda)
                 bounds(i,1) = MIN(bounds(i,1),uijbar)
                 bounds(i,2) = MAX(bounds(i,2),uijbar)
                 bounds(j,1) = MIN(bounds(j,1),uijbar)
                 bounds(j,2) = MAX(bounds(j,2),uijbar)
                 !===End compute low-order update to estimate bounds
              END IF
           END IF
        END DO

     END DO

     CALL MatAssemblyBegin(this%matrices%dijL, MAT_FINAL_ASSEMBLY, ierr)
     CALL MatAssemblyEnd  (this%matrices%dijL, MAT_FINAL_ASSEMBLY, ierr)

     IF (this%method=='high') THEN
        CALL MatAssemblyBegin(this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)
        CALL MatAssemblyEnd  (this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)
        !VB: CORRECTION ==> compatibility with several procs
        CALL array_to_petsc_vec(bounds(:,1), this%x1vec, this%mesh, this%LA, 'insert')
        CALL VecGhostGetLocalForm(this%x1vec, this%x2_ghost, ierr)
        CALL VecGhostUpdateBegin(this%x1vec, MIN_VALUES, SCATTER_FORWARD, ierr)
        CALL VecGhostUpdateEnd(this%x1vec, MIN_VALUES, SCATTER_FORWARD, ierr)
        CALL extract(this%x2_ghost, 1, 1, this%LA, bounds(:, 1))

        CALL array_to_petsc_vec(bounds(:,2), this%x1vec, this%mesh, this%LA, 'insert')
        CALL VecGhostGetLocalForm(this%x1vec, this%x2_ghost, ierr)
        CALL VecGhostUpdateBegin(this%x1vec, MAX_VALUES, SCATTER_FORWARD, ierr)
        CALL VecGhostUpdateEnd(this%x1vec, MAX_VALUES, SCATTER_FORWARD, ierr)
        CALL extract(this%x2_ghost, 1, 1, this%LA, bounds(:, 2))
        !VB: CORRECTION ==> compatibility with several procs
     END IF
   END SUBROUTINE compute_dij

   SUBROUTINE divide_by_mass(this,rk)
     USE petsc_tools
     IMPLICIT NONE
     CLASS(euler_type) :: this

     REAL(KIND = 8), DIMENSION(:) :: rk
     REAL(KIND = 8), DIMENSION(SIZE(rk)) :: rk_cp
     INTEGER :: ierr
     rk = rk/this%matrices%lumped_mass
     CALL array_to_petsc_vec(rk, this%x1vec, this%mesh, this%LA, 'insert')
     CALL MatMult(this%matrices%mass, this%x1vec, this%x2vec, ierr)
     CALL VecGhostGetLocalForm(this%x2vec, this%x2_ghost, ierr)
     CALL VecGhostUpdateBegin(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     CALL VecGhostUpdateEnd(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     CALL extract(this%x2_ghost, 1, 1, this%LA, rk_cp)
     rk = 2*rk - rk_cp/this%matrices%lumped_mass
   END SUBROUTINE divide_by_mass

   SUBROUTINE commutator(this, un, alpha)
     USE space_dim
     USE petsc_tools
     USE st_matrix
     USE sub_plot
     IMPLICIT NONE
     CLASS(euler_type) :: this
     REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: un
     REAL(KIND = 8), DIMENSION(:), INTENT(OUT):: alpha
     REAL(KIND = 8), DIMENSION(this%mesh%np)  :: rk, rk_norm, eta, logeta
     INTEGER :: k, ierr, np_tot
     REAL(KIND = 8) :: norm_diff, norm_log, pi=ACOS(-1.d0)
     CHARACTER(5) :: char
     PetscReal :: norm
     CALL VecGetSize(this%x5vec, np_tot, ierr)
     !===
     CALL VecSet(this%x4vec, 0.d0, ierr)
     CALL VecSet(this%x5vec, 0.d0, ierr)
     !eta = pressure_from_state(this, un)/un(:,1)**1.4
     !eta = pressure_from_state(this, un)
     eta = un(:,1)
     logeta = log(abs(eta))
     norm_diff = 0.d0
     norm_log = 0.d0

     DO k = 1, k_dim
        CALL array_to_petsc_vec(logeta, this%x1vec, this%mesh, this%LA, 'insert') !<==v1 = log(eta)
        CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr) !<==v2 = dk(log(eta))
        
        CALL array_to_petsc_vec(eta,    this%x1vec, this%mesh, this%LA, 'insert') !<==v1 = eta
        CALL VecPointwiseMult(this%x3vec,  this%x1vec, this%x2vec, ierr) !<==v3 = eta*dk(log(eta))

        CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)          !<==v2 = dk(eta))
        CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr) !<==v3 = eta*dk(log(eta)) - dk(eta)
  
        CALL VecNorm(this%x3vec, Norm_1, norm, ierr)
        norm_diff = norm_diff + norm

        CALL VecNorm(this%x2vec, Norm_1, norm, ierr)
        norm_log = norm_log + norm
        
        CALL VecAbs(this%x3vec,ierr)
        CALL VecAXPY(this%x4vec, 1.d0, this%x3vec, ierr) !<==v4 = sum_k |dk(eta)-eta*dk(log(eta))|

        CALL VecAbs(this%x2vec,ierr)
        CALL VecAXPY(this%x5vec, 1.d0, this%x2vec, ierr) !<==v5 = sum_k |dk(eta)
     END DO
     !IF (this%mesh%rank==0) write(*,*) 'error', norm_diff/norm_log, norm_diff, norm_log
     
     CALL VecGhostGetLocalForm(this%x4vec, this%x2_ghost, ierr)
     CALL VecGhostUpdateBegin(this%x4vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     CALL VecGhostUpdateEnd(this%x4vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     CALL extract(this%x2_ghost, 1, 1, this%LA, rk)

     CALL VecGhostGetLocalForm(this%x5vec, this%x2_ghost, ierr)
     CALL VecGhostUpdateBegin(this%x5vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     CALL VecGhostUpdateEnd(this%x5vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
     CALL extract(this%x2_ghost, 1, 1, this%LA, rk_norm)
    
     norm_log = norm_log/np_tot
  
     rk = abs(rk)/max(abs(rk_norm),1.d-1*norm_log)
     alpha = MIN(10*rk,1.d0)
     alpha = threshold(alpha)
    
     !IF (this%time+1.1*this%dt>this%final_time .AND. stage==this%ERK%s+1) THEN
     IF (this%time+1.5*this%dt>this%final_time) THEN
        WRITE(char, '(I5)') this%mesh%rank
        CALL plot_scalar_field(this%mesh%jj, this%mesh%rr, alpha, 'a'//trim(adjustl(char))//'.plt')
        CALL plot_scalar_field(this%mesh%jj, this%mesh%rr, eta, 'eta'//trim(adjustl(char))//'.plt')
     END IF
   END SUBROUTINE commutator

 FUNCTION threshold(x) RESULT(g)
    IMPLICIT NONE
    INTEGER, PARAMETER :: exp=3
    REAL(KIND=8), DIMENSION(:)  :: x
    REAL(KIND=8), DIMENSION(SIZE(x))  :: z, t, zp, relu, f, g
    REAL(KIND=8), PARAMETER :: x0 = .5d0
    REAL(KIND=8), PARAMETER :: x1=SQRT(3.d0)*x0
    SELECT CASE(exp)
    CASE(2)
       !===Quadratic threshold
       z = x-x0
       zp = x-2*x0
       relu = (zp+ABS(zp))/2
       f = -z*(z**2-x1**2)  + relu*(z-x0)*(z+2*x0)
       g = (f + 2*x0**3)/(4*x0**3)
    CASE(3)
       !===Cubic threshold
       relu = ((x-2*x0)+ABS(x-2*x0))/2
       t = x/(2*x0)
       g = t**3*(10-15*t+6*t**2) - relu*(t-1)**2*(6*t**2+3*t+1)/(2*x0)
    END SELECT
    RETURN
  END FUNCTION threshold



  


  !GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE
  !GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE
  !GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE
  SUBROUTINE compute_dk (this, un)
    USE arbitrary_eos_lambda_module
    IMPLICIT NONE
    CLASS(euler_type) :: this
    REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
    INTEGER, DIMENSION(1) :: i_t, j_t
    REAL(KIND = 8), DIMENSION(1, this%mesh%gauss%k_d) :: nij_c
    REAL(KIND = 8), DIMENSION(1) :: norm_c, dijL_c
    REAL(KIND = 8), DIMENSION(2) :: u, rho, ie, p, lambda_max
    LOGICAL, DIMENSION(this%mesh%medge) :: virgin_edge
    REAL(KIND = 8) :: pstar
    LOGICAL :: bug
    INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge, divider, nb_shared_cell
    nw = this%mesh%gauss%n_w

    bug = .FALSE.
    SELECT CASE(this%mesh%gauss%k_d)
    CASE(1)
       nb_shared_cell = 1
       IF (this%mesh%gauss%n_w/=2) bug=.TRUE.
    CASE(2)
       nb_shared_cell = 2
       IF (this%mesh%gauss%n_w/=3) bug=.TRUE.
    END SELECT
    IF (bug) THEN
       CALL error_petsc('Wrong polynomial degree for low-order viscosity')
    END IF

    DO m = 1, this%mesh%dom_me
       DO n = 1, this%mesh%gauss%n_e
          IF (this%mesh%attr_e(this%mesh%jce(n, m))) THEN
             edge = this%mesh%jce_loc(n, m)
             IF (.NOT. virgin_edge(edge)) CYCLE
             virgin_edge(edge) = .FALSE.
             ni = MOD(n, nw) + 1
             nj = MOD(n + 1, nw) + 1
             i = this%mesh%jj(ni, m)
             j = this%mesh%jj(nj, m)
             i_t = i
             j_t = j
             DO k = 1, this%mesh%gauss%k_d
                CALL MatGetValues(this%matrices%nij_loc(k), 1, i_t - 1, 1, j_t - 1, &
                     nij_c(:, k), ierr)
             END DO
             rho(1) = un(i, 1)
             rho(2) = un(j, 1)
             u(1) = SUM(un(i, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(1)
             u(2) = SUM(un(j, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(2)
             ie(1) = un(i, this%mesh%gauss%k_d + 2) / rho(1) - 0.5d0 * u(1) * u(1)
             ie(2) = un(j, this%mesh%gauss%k_d + 2) / rho(2) - 0.5d0 * u(2) * u(2)
             p = this%pressure(rho, ie)
             CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, &
                  lambda_max, pstar)
             dijL_c = MAXVAL(lambda_max) * norm_c
             divider = nb_shared_cell

             IF (this%mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j
                DO k = 1, this%mesh%gauss%k_d
                   CALL MatGetValues(this%matrices%nij_loc(k), 1, j_t - 1, 1, i_t - 1, &
                        nij_c(:, k), ierr)
                END DO
                u(1) = SUM(un(i, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(1)
                u(2) = SUM(un(j, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(2)
                rho = (/rho(2), rho(1)/)
                ie = (/ie(2), ie(1)/)
                p = (/p(2), p(1)/)
                CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, &
                     lambda_max, pstar)
                dijL_c = MAX(dijL_c, MAXVAL(lambda_max) * norm_c)
                divider = 1
             END IF

             this%matrices%dK(m) = MAX(this%matrices%dK(m),dijL_c(1)/divider)
          END IF
       END DO
    END DO
  END SUBROUTINE compute_dk

   SUBROUTINE compute_dt_from_dK(this)
     IMPLICIT NONE
     CLASS(euler_type) :: this
     REAL(KIND = 8), DIMENSION(this%mesh%dom_np) :: dijL_diag
     REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w) :: v_loc
     INTEGER, DIMENSION(this%mesh%gauss%n_w) :: idxm
     INTEGER :: i, m, ni, iglob
     REAL(KIND = 8) :: dt_min_loc, dt_min_glob
     Vec                                         :: vect
     PetscErrorCode                              :: ierr
     CALL VecSet(vect, 0.d0, ierr)

     WRITE(*,*) "VB: WARNING ==> this subroutine does not call any ghost points???"
     STOP

     DO m = 1, this%mesh%me
        v_loc = 0.d0
        DO ni = 1, this%mesh%gauss%n_w
           i = this%mesh%jj(ni, m)
           iglob = this%LA%loc_to_glob(1, i)
           idxm(ni) = iglob - 1
           v_loc(ni) = v_loc(ni) + this%matrices%dK(m)
        ENDDO
        CALL VecSetValues(vect, this%mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
     ENDDO
     CALL VecAssemblyBegin(vect, ierr)
     CALL VecAssemblyEnd(vect, ierr)

     CALL VecGetValues(this%vec_loc, this%mesh%dom_np, this%tab, dijL_diag, ierr)
     dijL_diag = this%matrices%lumped_mass(1:this%mesh%dom_np) / ABS(dijL_diag)

     dt_min_loc = MINVAL(dijL_diag) / 2.d0

     CALL MPI_ALLREDUCE(dt_min_loc, dt_min_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PETSC_COMM_WORLD, ierr)
     this%dt = this%CFL * dt_min_glob
   END SUBROUTINE compute_dt_from_dK

   SUBROUTINE compute_flux(this, ff, Vect)
      USE space_dim
      IMPLICIT NONE
      CLASS(euler_type) :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim) :: ff
      REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w) :: v_loc
      REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w, k_dim) :: f_loc
      REAL(KIND = 8), DIMENSION(this%mesh%np) :: v_glb
      REAL(KIND = 8), DIMENSION(this%mesh%me) :: volK
      REAL(KIND = 8), DIMENSION(this%mesh%gauss%l_G) :: wwrj
      INTEGER, DIMENSION(this%mesh%gauss%n_w) :: idxm, jj_loc
      REAL(KIND = 8) :: x
      INTEGER :: i, k, m, ni, nj, iglob
      Vec                                         :: vect
      PetscErrorCode                              :: ierr
      CALL VecSet(vect, 0.d0, ierr)
      v_glb = 0.d0
      DO m = 1, this%mesh%dom_me
         jj_loc = this%mesh%jj(:, m)
         f_loc = ff(jj_loc,:)
         !<==recompute cij on the fly
         DO ni = 1, this%mesh%gauss%n_w
            !wwrj = this%mesh%gauss%ww(ni,:)*this%mesh%gauss%rj(:,m)
            x = 0.d0
            DO k = 1, this%mesh%gauss%k_d
               DO nj = 1, this%mesh%gauss%n_w
                  x = x + f_loc(nj,k)* &
                        !SUM(this%mesh%gauss%dw(k,nj,:,m)*wwrj)
                     SUM(this%mesh%gauss%dw(k,nj,:,m)*this%mesh%gauss%ww(ni,:)*this%mesh%gauss%rj(:,m))
               ENDDO
            ENDDO
            v_loc(ni) = x
         ENDDO
         idxm = this%LA%loc_to_glob(1, jj_loc) -1
         v_loc = -v_loc
         CALL VecSetValues(vect, this%mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
   !!$        v_glb(jj_loc) = v_glb(jj_loc) - v_loc
      ENDDO
   !!$     CALL VecSetValues(vect, this%mesh%np, this%LA%loc_to_glob(1,:)-1, v_glb, INSERT_VALUES, ierr)
      CALL VecAssemblyBegin(vect, ierr)
      CALL VecAssemblyEnd(vect, ierr)
   END SUBROUTINE compute_flux


 END MODULE euler_type_MODULE
