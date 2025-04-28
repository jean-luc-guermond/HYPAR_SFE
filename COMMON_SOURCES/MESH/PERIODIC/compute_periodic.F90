MODULE compute_periodic
   USE def_type_mesh
   USE periodic_data_module
   IMPLICIT NONE

   PUBLIC :: periodic_matrix_petsc, periodic_rhs_petsc
   PRIVATE

CONTAINS

   SUBROUTINE periodic_matrix_petsc(periodic, LA, matrix)
      USE dyn_line_type
      USE def_type_mesh
      USE my_util
#include "petsc/finclude/petsc.h"
      USE petsc
      IMPLICIT NONE
      TYPE(periodic_type), INTENT(IN) :: periodic
      INTEGER :: n_bord
      TYPE(dyn_int_line), DIMENSION(:), POINTER :: list, perlist
      TYPE(petsc_csr_la), INTENT(IN) :: LA
      INTEGER, PARAMETER :: nmaxcols = 300
      INTEGER :: ncols
      INTEGER, DIMENSION(nmaxcols) :: cols
      REAL(KIND = 8), DIMENSION(nmaxcols) :: vals
      INTEGER, DIMENSION(:), ALLOCATABLE :: n_cols_i
      INTEGER, DIMENSION(1) :: idxn
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: jdxn
      REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: vals_pi
      INTEGER :: n, l, i, pi, n_D, k
      Mat                                          :: matrix
      PetscErrorCode                               :: ierr

      CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      CALL MatSetOption (matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)

      n_bord = periodic%n_bord
      list = periodic%list
      perlist = periodic%perlist

      DO k = 1, SIZE(LA%loc_to_glob, 1)
         DO n = 1, n_bord
            n_D = SIZE(list(n)%DIL)
            IF (n_D /=0) THEN
               ALLOCATE(jdxn(n_D, nmaxcols), vals_pi(n_D, nmaxcols), n_cols_i(n_D))
               jdxn = 0
               vals_pi = 0.d0
               n_cols_i = 0

               DO l = 1, SIZE(list(n)%DIL)
                  idxn(1) = LA%loc_to_glob(k, list(n)%DIL(l))
                  CALL MatGetRow(matrix, idxn(1) - 1, ncols, cols, vals, ierr)
                  n_cols_i(l) = ncols
                  jdxn(l, 1:ncols) = cols(1:ncols)
                  vals_pi(l, 1:ncols) = vals(1:ncols)
                  CALL MatRestoreRow(matrix, idxn(1) - 1, ncols, cols, vals, ierr)
               END DO

               DO l = 1, n_D
                  idxn(1) = LA%loc_to_glob(k, perlist(n)%DIL(l)) - 1
                  CALL MatSetValues(matrix, 1, idxn, n_cols_i(l), jdxn(l, 1:n_cols_i(l)), &
                       vals_pi(l:l, 1:n_cols_i(l)), ADD_VALUES, ierr)
               END DO
               DEALLOCATE(jdxn, vals_pi, n_cols_i)

            END IF
         END DO
         CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
         CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

         DO n = 1, n_bord
            n_D = SIZE(list(n)%DIL)
            !CALL MatZeroRows(matrix, n_D, LA%loc_to_glob(k,list(n)%DIL(:))-1, 1.d0, &
            !     PETSC_NULL_OBJECT, PETSC_NULL_OBJECT, ierr) !petsc.3.7.4
            CALL MatZeroRows(matrix, n_D, LA%loc_to_glob(k, list(n)%DIL(:)) - 1, 1.d0, &
                 PETSC_NULL_VEC, PETSC_NULL_VEC, ierr) !(JLG) Feb 20, 2019, petsc.3.8.4
         END DO
         !!$       CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
         !!$       CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

         DO n = 1, n_bord
            DO l = 1, SIZE(list(n)%DIL)
               i = LA%loc_to_glob(k, list(n)%DIL(l))
               pi = LA%loc_to_glob(k, perlist(n)%DIL(l))
               CALL MatSetValue(matrix, i - 1, pi - 1, -1.d0, INSERT_VALUES, ierr)
            END DO
         END DO
         CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
         CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

      END DO

   END SUBROUTINE periodic_matrix_petsc

   SUBROUTINE periodic_rhs_petsc(n_bord, list, perlist, v_rhs, LA)
      USE dyn_line_type
      USE def_type_mesh
#include "petsc/finclude/petsc.h"
      USE petsc
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n_bord
      TYPE(dyn_int_line), DIMENSION(:), INTENT(IN) :: list, perlist
      TYPE(petsc_csr_la), INTENT(IN) :: LA
      INTEGER, DIMENSION(:), ALLOCATABLE :: idxn, jdxn
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: vals, bs
      INTEGER :: n, k, n_D
      Vec                                          :: v_rhs
      PetscErrorCode                               :: ierr


      DO k = 1, SIZE(LA%loc_to_glob, 1)
         DO n = 1, n_bord
            n_D = SIZE(list(n)%DIL)
            ALLOCATE(idxn(n_D), vals(n_D), jdxn(n_D), bs(n_D))
            idxn = LA%loc_to_glob(k, list(n)%DIL(:)) - 1
            jdxn = LA%loc_to_glob(k, perlist(n)%DIL(:)) - 1
            CALL VecGetValues(v_rhs, n_D, idxn, vals, ierr)
            CALL VecAssemblyBegin(v_rhs, ierr)
            CALL VecAssemblyEnd(v_rhs, ierr)

            bs = 0.d0
            CALL VecSetValues(v_rhs, n_D, jdxn, vals, ADD_VALUES, ierr)
            CALL VecAssemblyBegin(v_rhs, ierr)
            CALL VecAssemblyEnd(v_rhs, ierr)
            CALL VecSetValues(v_rhs, n_D, idxn, bs, INSERT_VALUES, ierr)
            CALL VecAssemblyBegin(v_rhs, ierr)
            CALL VecAssemblyEnd(v_rhs, ierr)

            IF (ALLOCATED(idxn)) DEALLOCATE(idxn, jdxn, vals, bs)
         END DO
      END DO

   END SUBROUTINE periodic_rhs_petsc

END MODULE compute_periodic