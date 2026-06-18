MODULE hyperbolic_matrices_module
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE solver_petsc
   USE periodic_data_module
   USE compute_periodic
   USE fem_petsc_matrix_factory_module, &
              ONLY : construct_lumped_mass_vector, construct_cij
   TYPE hyperbolic_matrices_type
      INTEGER                          :: method, which_mass
      REAL(KIND=8), DIMENSION(:,:,:),   POINTER :: cij_norm_loc_array
      REAL(KIND=8), DIMENSION(:,:,:,:), POINTER :: nij_loc_array

      Mat :: mass, dijL, stiffL, cij_norm_loc, Al_mass
      Mat :: dijH
      Mat, DIMENSION(:),   ALLOCATABLE :: cij, nij_loc
      Mat, DIMENSION(:,:), ALLOCATABLE :: cij_loc
      Vec                              :: lump_mass_vec
      KSP                              :: ksp_consistent_mass
   CONTAINS
      PROCEDURE, PUBLIC :: construct => construct_hyperbolic_matrices
      PROCEDURE, PRIVATE :: construct_loc_nij, construct_consistent_mass_solver, construct_Ar
   END TYPE hyperbolic_matrices_type

   INTEGER, PRIVATE, PARAMETER :: METHOD_VISCOUS=1, METHOD_HIGH=2
   INTEGER, PRIVATE, PARAMETER :: LUMPED_MASS=1, QUASI_CONSISTENT_MASS=2, CONSISTENT_MASS=3

CONTAINS

   SUBROUTINE construct_hyperbolic_matrices(this, communicator, mesh, LA)
      USE space_dim
      USE fem_M
      USE st_matrix, ONLY: create_my_ghost
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      TYPE(mesh_type),    INTENT(IN) :: mesh
      type(petsc_csr_LA), INTENT(IN) :: LA
      INTEGER                        :: k, ierr
      INTEGER, DIMENSION(:), POINTER :: ifrom
      MPI_Comm                       :: communicator
      IS,      DIMENSION(1)          :: is

      !===Init global vectors
      IF (.NOT. ALLOCATED(this%cij)) THEN
         ALLOCATE(this%cij(k_dim))
         ALLOCATE(this%nij_loc(k_dim))
         ALLOCATE(this%cij_loc(1, k_dim))
      END IF

      !===Mat allocations
      CALL create_local_petsc_matrix(communicator, LA, this%mass, clean = .FALSE.)
      CALL MatDuplicate(this%mass, MAT_DO_NOT_COPY_VALUES, this%dijL, ierr)
      DO k = 1, k_dim
         CALL create_local_petsc_matrix(communicator, LA, this%cij(k), clean = .FALSE.)
      END DO

      !===Mat construction
      CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, this%mass)
      CALL periodic_matrix_petsc(mesh%per, LA, this%mass)
      
      CALL create_my_ghost(mesh, LA, ifrom)
      CALL VecCreateGhost(communicator, mesh%dom_np, &
           PETSC_DETERMINE, SIZE(ifrom), ifrom, this%lump_mass_vec, ierr)
      CALL construct_lumped_mass_vector(mesh, LA, this%mass, this%lump_mass_vec)
      

      CALL construct_cij(mesh, LA, this%cij)

      CALL ISCreateGeneral(communicator, mesh%np, LA%loc_to_glob(1, :) - 1, PETSC_COPY_VALUES, is(1), ierr)
      DO k = 1, k_dim
         CALL MatCreateSubMatrices(this%cij(k), 1, is, is, MAT_INITIAL_MATRIX, this%cij_loc(:, k), ierr)
         CALL MatDuplicate(this%cij_loc(1, k), MAT_DO_NOT_COPY_VALUES, this%nij_loc(k), ierr)
      END DO
      CALL MatDuplicate(this%cij_loc(1, 1), MAT_DO_NOT_COPY_VALUES, this%cij_norm_loc, ierr)
      ALLOCATE(this%cij_norm_loc_array(mesh%gauss%n_w, mesh%gauss%n_w, mesh%me))
      ALLOCATE(this%nij_loc_array(mesh%gauss%n_w, mesh%gauss%n_w, k_dim, mesh%me))
      CALL this%construct_loc_nij(mesh)

      IF (this%method==METHOD_HIGH) THEN
         CALL MatDuplicate(this%mass, MAT_DO_NOT_COPY_VALUES, this%dijH, ierr)
         CALL MatDuplicate(this%mass, MAT_DO_NOT_COPY_VALUES, this%stiffL, ierr)
         CALL qs_mass_diff_M (mesh, 0.d0, 1.d0, LA, this%stiffL)
      END IF

      IF (this%which_mass==CONSISTENT_MASS) THEN
         CALL this%construct_consistent_mass_solver(communicator)
      ELSE IF (this%which_mass==QUASI_CONSISTENT_MASS) THEN
         CALL this%construct_Ar
      END IF

   END SUBROUTINE construct_hyperbolic_matrices

   SUBROUTINE construct_loc_nij(this, mesh)
      USE space_dim
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      TYPE(mesh_type), INTENT(IN) :: mesh
      REAL(KIND = 8), DIMENSION(1, k_dim) :: cij_c
      REAL(KIND = 8), DIMENSION(1, 1) :: nij_c, norm
      REAL(KIND =8) :: xx
      INTEGER, DIMENSION(1) :: i, j
      LOGICAL, DIMENSION(mesh%medge) :: virgin_edge
      INTEGER :: k, m, n, ni, nj, ierr, nw, edge

      nw = mesh%gauss%n_w
      virgin_edge = .true.
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

               !=== fill blocks (i, j)
               norm = 0.d0
               DO k = 1, k_dim
                  CALL MatGetValues(this%cij_loc(1, k), 1, i - 1, 1, j - 1, cij_c(:, k), ierr)
                  norm = norm + cij_c(1, k)**2
               END DO
               norm = SQRT(norm)

               CALL MatSetValues(this%cij_norm_loc, 1, i - 1, 1, j - 1, norm, ADD_VALUES, ierr)
               DO k = 1, k_dim
                  nij_c = cij_c(1, k) / norm
                  CALL MatSetValues(this%nij_loc(k), 1, i - 1, 1, j - 1, nij_c, ADD_VALUES, ierr)
               END DO

               !=== fill blocks (j, i)
               norm = 0.d0
               DO k = 1, k_dim
                  CALL MatGetValues(this%cij_loc(1, k), 1, j - 1, 1, i - 1, cij_c(:, k), ierr)
                  norm = norm + cij_c(1, k)**2
               END DO
               norm = SQRT(norm)

               CALL MatSetValues(this%cij_norm_loc, 1, j - 1, 1, i - 1, norm, ADD_VALUES, ierr)
               DO k = 1, k_dim
                  nij_c = cij_c(1, k) / norm
                  CALL MatSetValues(this%nij_loc(k), 1, j - 1, 1, i - 1, nij_c, ADD_VALUES, ierr)
               END DO

            END IF
         END DO
      END DO

      DO k = 1, k_dim
         CALL MatAssemblyBegin(this%nij_loc(k), MAT_FINAL_ASSEMBLY, ierr)
         CALL MatAssemblyEnd  (this%nij_loc(k), MAT_FINAL_ASSEMBLY, ierr)
      END DO

      CALL MatAssemblyBegin(this%cij_norm_loc, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd  (this%cij_norm_loc, MAT_FINAL_ASSEMBLY, ierr)


      this%cij_norm_loc_array = 0.d0
      this%nij_loc_array = 0.d0
      DO m = 1, mesh%me
         DO ni = 1, mesh%gauss%n_w
            DO nj=1, mesh%gauss%n_w
               IF (ni==nj) CYCLE
               i = mesh%jj(ni, m)
               j = mesh%jj(nj, m)

               !=== fill blocks (i, j)
               xx = 0.d0
               DO k = 1, k_dim
                  CALL MatGetValues(this%cij_loc(1, k), 1, i - 1, 1, j - 1, cij_c(:, k), ierr)
                  xx = xx + cij_c(1, k)**2
               END DO
               xx = SQRT(xx)
               this%cij_norm_loc_array(ni, nj, m) = xx

               DO k = 1, k_dim
                  this%nij_loc_array(ni, nj, k, m) = cij_c(1, k)/xx
               END DO
               ! !=== fill blocks (i, j)
               ! norm = 0.d0
               ! DO k = 1, k_dim
               !    CALL MatGetValues(this%cij_loc(1, k), 1, i - 1, 1, j - 1, cij_c(:, k), ierr)
               !    norm = norm + cij_c(1, k)**2
               ! END DO
               ! norm = SQRT(norm)

               ! CALL MatSetValues(this%cij_norm_loc, 1, i - 1, 1, j - 1, norm, ADD_VALUES, ierr)
               ! DO k = 1, k_dim
               !    nij_c = cij_c(1, k) / norm
               !    CALL MatSetValues(this%nij_loc(k), 1, i - 1, 1, j - 1, nij_c, ADD_VALUES, ierr)
               ! END DO
            END DO
         END DO
      END DO


   END SUBROUTINE construct_loc_nij

   SUBROUTINE construct_consistent_mass_solver(this, communicator)
      USE solver_petsc
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      TYPE(solver_param)              :: my_par_solver
      MPI_Comm                       :: communicator

      !===Create ksp solver
      my_par_solver%it_max  = 5000
      my_par_solver%rel_tol = 1.d-10
      my_par_solver%abs_tol = 1.d-18
      my_par_solver%verbose = .FALSE.
      CALL init_solver(my_par_solver, this%ksp_consistent_mass, this%mass, communicator, &
      solver = "GMRES", precond = 'MUMPS')
   END SUBROUTINE construct_consistent_mass_solver

   SUBROUTINE construct_Ar(this)
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      INTEGER :: ierr
      Vec     :: xx
      CALL VecDuplicate(this%lump_mass_vec, xx, ierr)
      CALL VecSet(xx, 1.d0, ierr)
      CALL MatDuplicate(this%mass, MAT_COPY_VALUES, this%Al_mass, ierr)
      CALL VecPointWiseDivide(xx, xx, this%lump_mass_vec, ierr)
      CALL MatDiagonalScale(this%Al_mass, PETSC_NULL_VEC, xx, ierr)
      CALL MatScale(this%Al_mass, -1.d0, ierr)
      CALL MatShift(this%Al_mass, 1.d0, ierr)
      CALL VecDestroy(xx, ierr)
   END SUBROUTINE construct_Ar

END MODULE hyperbolic_matrices_module
