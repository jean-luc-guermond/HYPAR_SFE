MODULE fem_rhs
#include "petsc/finclude/petsc.h"
   USE petsc
 CONTAINS

   SUBROUTINE qs_00 (mesh, LA, ff, vect)
      !=================================
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type), TARGET :: mesh
      type(petsc_csr_LA) :: LA
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: ff
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w) :: ff_loc
      INTEGER, DIMENSION(mesh%gauss%n_w) :: jj_loc
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w) :: v_loc
      INTEGER, DIMENSION(mesh%gauss%n_w) :: idxm
      INTEGER :: i, m, l, ni, iglob
      REAL(KIND = 8) :: fl
      Vec                                         :: vect
      PetscErrorCode                              :: ierr
      CALL VecSet(vect, 0.d0, ierr)

      DO m = 1, mesh%dom_me
         jj_loc = mesh%jj(:, m)
         ff_loc = ff(jj_loc)
         DO ni = 1, mesh%gauss%n_w
            i = mesh%jj(ni, m)
            iglob = LA%loc_to_glob(1, i)
            idxm(ni) = iglob - 1
         END DO

         v_loc = 0.d0
         DO l = 1, mesh%gauss%l_G
            fl = SUM(ff_loc * mesh%gauss%ww(:, l)) * mesh%gauss%rj(l, m)
            DO ni = 1, mesh%gauss%n_w
               v_loc(ni) = v_loc(ni) + mesh%gauss%ww(ni, l) * fl
            END DO
         ENDDO
         CALL VecSetValues(vect, mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
      ENDDO
      CALL VecAssemblyBegin(vect, ierr)
      CALL VecAssemblyEnd(vect, ierr)
   END SUBROUTINE qs_00

END MODULE fem_rhs
