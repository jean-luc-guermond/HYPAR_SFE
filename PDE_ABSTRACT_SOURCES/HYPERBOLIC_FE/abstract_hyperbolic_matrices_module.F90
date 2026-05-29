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
      CHARACTER(LEN=200) :: method
      Mat :: mass, dijL, stiffL, cij_norm_loc
      Mat :: dijH
      Mat, DIMENSION(:),   ALLOCATABLE :: cij, nij_loc
      Mat, DIMENSION(:,:), ALLOCATABLE :: cij_loc
      Vec                              :: lump_mass_vec
   CONTAINS
      PROCEDURE, PUBLIC :: construct => construct_hyperbolic_matrices
      PROCEDURE, PRIVATE :: construct_loc_nij
   END TYPE hyperbolic_matrices_type

CONTAINS

   SUBROUTINE construct_hyperbolic_matrices(this, communicator, mesh, LA)
      USE space_dim
      USE fem_M
      USE st_matrix, ONLY: create_my_ghost
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      TYPE(mesh_type), INTENT(IN) :: mesh
      type(petsc_csr_LA), INTENT(IN) :: LA
      INTEGER :: k, ierr
      INTEGER, DIMENSION(:), POINTER :: ifrom
      real(kind=8) :: norm
      MPI_Comm       :: communicator
      IS, DIMENSION(1) :: is

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
      !>>> because this subroutine is designed to solve AX=Y but does not technically periodize A
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
      CALL this%construct_loc_nij(mesh)

      IF (TRIM(ADJUSTL(this%method))=='high') THEN
         CALL MatDuplicate(this%mass, MAT_DO_NOT_COPY_VALUES, this%dijH, ierr)
         CALL MatDuplicate(this%mass, MAT_DO_NOT_COPY_VALUES, this%stiffL, ierr)
         CALL qs_mass_diff_M (mesh, 0.d0, 1.d0, LA, this%stiffL)
      END IF

   END SUBROUTINE construct_hyperbolic_matrices

   SUBROUTINE construct_loc_nij(this, mesh)
      USE space_dim
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      TYPE(mesh_type), INTENT(IN) :: mesh
      REAL(KIND = 8), DIMENSION(1, k_dim) :: cij_c
      REAL(KIND = 8), DIMENSION(1, 1) :: norm, nij_c
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

   END SUBROUTINE construct_loc_nij

END MODULE hyperbolic_matrices_module
