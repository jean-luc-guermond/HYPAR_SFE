MODULE petsc_tools
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh

CONTAINS
   SUBROUTINE array_to_petsc_vec(uu, xx, mesh, LA, operation)
      USE my_util
      IMPLICIT NONE
      TYPE(petsc_csr_LA), INTENT(IN) :: LA
      TYPE(mesh_type), INTENT(IN) :: mesh
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: uu
      CHARACTER(LEN = 6), INTENT(IN) :: operation
      INTEGER, DIMENSION(SIZE(uu)) :: idxm
      INTEGER :: i, ierr
      Vec     :: xx
      IF (mesh%np.NE.SIZE(uu)) THEN
         CALL error_Petsc('Bug: array_to_petsc_vec, mesh%np>SIZE(uu)')
      END IF
      DO i = 1, mesh%np
         idxm(i) = LA%loc_to_glob(1, i) - 1
      END DO
      SELECT CASE (operation)
      CASE('insert')
         CALL VecSetValues(xx, mesh%np, idxm, uu, INSERT_VALUES, ierr)
      CASE('add')
         CALL VecSetValues(xx, mesh%np, idxm, uu, ADD_VALUES, ierr)
      CASE DEFAULT
         CALL error_petsc('Wrong option in array_to_petsc_vec for operation.')
      END SELECT

      CALL VecAssemblyBegin(xx, ierr)
      CALL VecAssemblyEnd(xx, ierr)
   END SUBROUTINE array_to_petsc_vec

END MODULE petsc_tools

