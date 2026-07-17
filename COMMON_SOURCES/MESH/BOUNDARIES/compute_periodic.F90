
MODULE compute_periodic
   USE def_type_mesh
   USE periodic_data_module
   IMPLICIT NONE

   PUBLIC :: periodic_matrix_petsc, periodic_rhs_petsc, periodic_vector_petsc, periodic_add_vector_petsc
   PRIVATE

CONTAINS

  SUBROUTINE periodic_matrix_petsc(periodic, LA, matrix)
    USE dyn_line_type
    USE def_type_mesh
    USE petsc_csr_LA_module
    USE my_util
    USE petscmat
    IMPLICIT NONE
    TYPE(periodic_type), INTENT(IN) :: periodic
    INTEGER :: nb_per_edges
    TYPE(petsc_csr_la), INTENT(IN) :: LA
    INTEGER, PARAMETER :: nmaxcols = 300
    INTEGER :: ncols
#if (PETSC_VERSION_MINOR < 23)
    INTEGER, DIMENSION(nmaxcols) :: cols
    REAL(KIND = 8), DIMENSION(nmaxcols) :: vals
#else
    INTEGER,        DIMENSION(:), POINTER :: cols
    REAL(KIND = 8), DIMENSION(:), POINTER :: vals
#endif
    INTEGER, DIMENSION(:), ALLOCATABLE :: n_cols_i
    INTEGER, DIMENSION(1) :: idxn
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: jdxn
    REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: vals_pi
    INTEGER :: n, l, i, pi, n_D, k
    TYPE(tMat)                                          :: matrix
    INTEGER                               :: ierr

   !  CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_TRUE, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
    CALL MatSetOption (matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)

    nb_per_edges = periodic%nb_bords

    DO k = 1, SIZE(LA%loc_to_glob, 1)
   !  DO k = 1, SIZE(LA%loc_to_glob, 1)
       DO n = 1, nb_per_edges
         !  CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
         !  CALL MatSetOption (matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
          n_D = SIZE(periodic%list(n)%DIL)
         !  IF (n_D /=0) THEN
             ALLOCATE(jdxn(n_D, nmaxcols))
             ALLOCATE(vals_pi(n_D, nmaxcols), SOURCE=0.d0)
             ALLOCATE(n_cols_i(n_D), SOURCE=0)
            !  vals_pi = 0.d0
            !  n_cols_i = 0

             DO l = 1, n_D!SIZE(periodic%list(n)%DIL)
                idxn(1) = LA%loc_to_glob(k, periodic%list(n)%DIL(l)) - 1
                CALL MatGetRow(matrix, idxn(1), ncols, cols, vals, ierr)
                n_cols_i(l) = ncols
                jdxn(l, 1:ncols) = cols(1:ncols)
                vals_pi(l, 1:ncols) = vals(1:ncols)
                CALL MatRestoreRow(matrix, idxn(1), ncols, cols, vals, ierr)
             END DO

             DO l = 1, n_D
                idxn(1) = LA%loc_to_glob(k, periodic%perlist(n)%DIL(l)) - 1
#if (23 <= PETSC_VERSION_MINOR) && (PETSC_VERSION_MINOR < 25)
                CALL MatSetValues(matrix, 1, idxn, n_cols_i(l), jdxn(l, 1:n_cols_i(l)), &
                     vals_pi(l, 1:n_cols_i(l)), ADD_VALUES, ierr)
#else
                CALL MatSetValues(matrix, 1, idxn, n_cols_i(l), jdxn(l, 1:n_cols_i(l)), &
                     vals_pi(l:l, 1:n_cols_i(l)), ADD_VALUES, ierr)
#endif
             END DO
             DEALLOCATE(jdxn, vals_pi, n_cols_i)
! if (k==2) stop

         !  END IF
          CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
          CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
       END DO

       DO n = 1, nb_per_edges
          n_D = SIZE(periodic%list(n)%DIL)
          CALL MatZeroRows(matrix, n_D, LA%loc_to_glob(k, periodic%list(n)%DIL(:)) - 1, 1.d0, &
               PETSC_NULL_VEC, PETSC_NULL_VEC, ierr) !(JLG) Feb 20, 2019, petsc.3.8.4
       END DO

       DO n = 1, nb_per_edges
          DO l = 1, SIZE(periodic%list(n)%DIL)
             i = LA%loc_to_glob(k, periodic%list(n)%DIL(l))
             pi = LA%loc_to_glob(k, periodic%perlist(n)%DIL(l))
             CALL MatSetValue(matrix, i - 1, pi - 1, -1.d0, INSERT_VALUES, ierr)
          END DO
       END DO
       CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
       CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

    END DO
  END SUBROUTINE periodic_matrix_petsc

   SUBROUTINE periodic_rhs_petsc(list, perlist, v_rhs, LA)
      USE dyn_line_type
      USE def_type_mesh
      USE petsc_csr_LA_module

      USE petscvec
      IMPLICIT NONE
      TYPE(dyn_int_line), DIMENSION(:), INTENT(IN) :: list, perlist
      TYPE(petsc_csr_la),               INTENT(IN) :: LA
      INTEGER                            :: nb_per_edges
      INTEGER, DIMENSION(:), ALLOCATABLE :: idxn, jdxn
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: vals, bs
      INTEGER :: n, k, n_D
      TYPE(tVec)                                          :: v_rhs
      INTEGER                               :: ierr

      nb_per_edges = SIZE(list)

      DO k = 1, SIZE(LA%loc_to_glob, 1)
         DO n = 1, nb_per_edges
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

   SUBROUTINE periodic_vector_petsc(list, perlist, vec_in, LA)
      USE dyn_line_type
      USE def_type_mesh
      USE petsc_csr_LA_module

      USE petscvec
      IMPLICIT NONE
      TYPE(dyn_int_line), DIMENSION(:), INTENT(IN) :: list, perlist
      TYPE(petsc_csr_la),               INTENT(IN) :: LA
      INTEGER                            :: nb_per_edges
      INTEGER, DIMENSION(:), ALLOCATABLE :: idxn, jdxn
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: vals, bs
      INTEGER :: n, k, n_D
      TYPE(tVec)                                          :: vec_in
      INTEGER                               :: ierr

      nb_per_edges = SIZE(list)

      DO k = 1, SIZE(LA%loc_to_glob, 1)
         DO n = 1, nb_per_edges
            n_D = SIZE(list(n)%DIL)
            ALLOCATE(idxn(n_D), vals(n_D), jdxn(n_D), bs(n_D))
            idxn = LA%loc_to_glob(k, list(n)%DIL(:)) - 1
            jdxn = LA%loc_to_glob(k, perlist(n)%DIL(:)) - 1
            CALL VecGetValues(vec_in, n_D, jdxn, vals, ierr)
            CALL VecAssemblyBegin(vec_in, ierr)
            CALL VecAssemblyEnd(vec_in, ierr)

            CALL VecSetValues(vec_in, n_D, idxn, vals, INSERT_VALUES, ierr)
            CALL VecAssemblyBegin(vec_in, ierr)
            CALL VecAssemblyEnd(vec_in, ierr)

            IF (ALLOCATED(idxn)) DEALLOCATE(idxn, jdxn, vals, bs)
         END DO
      END DO

   END SUBROUTINE periodic_vector_petsc

   SUBROUTINE periodic_add_vector_petsc(list, perlist, vec_in, LA)
      USE dyn_line_type
      USE def_type_mesh
      USE petsc_csr_LA_module
      USE petscvec
      IMPLICIT NONE
      TYPE(dyn_int_line), DIMENSION(:), INTENT(IN) :: list, perlist
      TYPE(petsc_csr_la),               INTENT(IN) :: LA
      INTEGER                            :: nb_per_edges
      INTEGER, DIMENSION(:), ALLOCATABLE :: idxn, jdxn
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: vals_list, vals_perlist, bs
      INTEGER :: n, k, n_D
      TYPE(tVec)                                          :: vec_in
      INTEGER                               :: ierr

      nb_per_edges = SIZE(list)

      DO k = 1, SIZE(LA%loc_to_glob, 1)
         DO n = 1, nb_per_edges
            n_D = SIZE(list(n)%DIL)
            ALLOCATE(idxn(n_D), vals_list(n_D), vals_perlist(n_D), jdxn(n_D), bs(n_D))
            idxn = LA%loc_to_glob(k, list(n)%DIL(:)) - 1
            jdxn = LA%loc_to_glob(k, perlist(n)%DIL(:)) - 1
            CALL VecGetValues(vec_in, n_D, idxn, vals_list, ierr)
            CALL VecGetValues(vec_in, n_D, jdxn, vals_perlist, ierr)
            CALL VecAssemblyBegin(vec_in, ierr)
            CALL VecAssemblyEnd(vec_in, ierr)

            CALL VecSetValues(vec_in, n_D, idxn, vals_perlist, ADD_VALUES, ierr)
            CALL VecSetValues(vec_in, n_D, jdxn, vals_list, ADD_VALUES, ierr)
            CALL VecAssemblyBegin(vec_in, ierr)
            CALL VecAssemblyEnd(vec_in, ierr)

            IF (ALLOCATED(idxn)) DEALLOCATE(idxn, jdxn, vals_list, vals_perlist, bs)
         END DO
      END DO

   END SUBROUTINE periodic_add_vector_petsc

END MODULE compute_periodic
