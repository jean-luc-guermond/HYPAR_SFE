MODULE petsc_tools
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh

CONTAINS

   SUBROUTINE array_to_petsc_vec(uu, xx, LA, operation)
      !> convert uu(HYPAR vector) to xx(petsc vec)
      !! LA
      !! operation = 'insert'/'add'/'min'/'max' (for ghost values)
      USE my_util
      IMPLICIT NONE
      TYPE(petsc_csr_LA),           INTENT(IN) :: LA
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: uu
      REAL(KIND = 8), DIMENSION(:), POINTER    :: x_loc
      CHARACTER(LEN = *),           INTENT(IN) :: operation
      INTEGER, DIMENSION(SIZE(uu)) :: idxm
      INTEGER :: i, np, k, kmax, ierr
      Vec     :: xx, xx_ghost

      kmax = SIZE(LA%loc_to_glob,1)
      np = SIZE(uu)/kmax

      IF (np /= SIZE(LA%loc_to_glob, 2)) THEN
         CALL error_petsc("BUG in array to petsc, wrong np")
      END IF

      DO k=1, kmax
         DO i = 1, np
            idxm(i+(k-1)*np) = LA%loc_to_glob(k, i) - 1
         END DO
      END DO
      SELECT CASE (operation)
      CASE('insert')
         CALL VecSetValues(xx, np*kmax, idxm, uu, INSERT_VALUES, ierr)
      CASE('add')
         CALL VecSetValues(xx, np*kmax, idxm, uu, ADD_VALUES, ierr)
      CASE('min')
         !=== Define ghost vectors LOCALLY
         CALL VecGhostGetLocalForm(xx, xx_ghost, ierr)
         CALL VecGetArrayF90(xx_ghost, x_loc, ierr)
         x_loc(:) = uu(:)
         CALL VecRestoreArrayF90(xx_ghost, x_loc, ierr)
         CALL VecGhostRestoreLocalForm(xx, xx_ghost, ierr)

         !=== At global level, computing min values on Ghost points across processes and:
            !=== 1) Update dom_np+1:np (compute min between overlapping values, update ghost values)
         CALL VecGhostUpdateBegin(xx, MIN_VALUES, SCATTER_REVERSE, ierr)
         CALL VecGhostUpdateEnd(xx, MIN_VALUES, SCATTER_REVERSE, ierr)
            !=== 2) Update 1:dom_np (update owned dofs with ghost values computed with min)
         CALL VecGhostUpdateBegin(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)

      CASE('max')
         !=== Define ghost vectors LOCALLY
         CALL VecGhostGetLocalForm(xx, xx_ghost, ierr)
         CALL VecGetArrayF90(xx_ghost, x_loc, ierr)
         x_loc(:) = uu(:)
         CALL VecRestoreArrayF90(xx_ghost, x_loc, ierr)
         CALL VecGhostRestoreLocalForm(xx, xx_ghost, ierr)
         !=== At global level, computing max values on Ghost points across processes and:
            !=== 1) Update dom_np+1:np (compute max between overlapping values, update ghost values)
         CALL VecGhostUpdateBegin(xx, MAX_VALUES, SCATTER_REVERSE, ierr)
         CALL VecGhostUpdateEnd(xx, MAX_VALUES, SCATTER_REVERSE, ierr)
            !=== 2) Update 1:dom_np (update owned dofs with ghost values computed with max)
         CALL VecGhostUpdateBegin(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)


      CASE DEFAULT
         CALL error_petsc(&
         &'Wrong option operation '//operation//' in array_to_petsc_vec ==> should be: add/insert/min/max')
      END SELECT

      CALL VecAssemblyBegin(xx, ierr)
      CALL VecAssemblyEnd(xx, ierr)
   END SUBROUTINE array_to_petsc_vec

END MODULE petsc_tools