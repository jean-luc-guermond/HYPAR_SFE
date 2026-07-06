MODULE fem_rhs
#include "petsc/finclude/petsc.h"
   USE petsc
   CONTAINS
   
   SUBROUTINE qs_00 (mesh, LA, ff, vect)
      !> transfer of  vect(PETSc projected on test functions) <== ff(HYPAR)
      !! mesh(hypar)
      !! LA(petsc_csr_MLA)
      !=================================
      USE def_type_mesh
      USE petsc_csr_LA_module
      IMPLICIT NONE
      TYPE(mesh_type)    :: mesh
      type(petsc_csr_LA) :: LA
      REAL(KIND = 8), DIMENSION(:), INTENT(IN)  :: ff
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w) :: ff_loc
      INTEGER, DIMENSION(mesh%gauss%n_w)        :: jj_loc
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w) :: v_loc
      INTEGER, DIMENSION(mesh%gauss%n_w)        :: idxm
      INTEGER        :: i, m, l, ni, iglob
      REAL(KIND = 8) :: fl
      Vec                                         :: vect
      PetscErrorCode                              :: ierr
      
      CALL VecZeroEntries(vect, ierr)
      
      DO m = 1, mesh%me
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
   
   SUBROUTINE qs_00_block (mesh, LA, ff, vect)
      !=================================
      USE def_type_mesh
      USE petsc_csr_LA_module
      IMPLICIT NONE
      TYPE(mesh_type), TARGET :: mesh
      type(petsc_csr_LA) :: LA
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: ff
      INTEGER, DIMENSION(mesh%gauss%n_w) :: jj_loc
      REAL(KIND = 8), DIMENSION(SIZE(LA%loc_to_glob, 1) * mesh%gauss%n_w) :: ff_loc
      REAL(KIND = 8), DIMENSION(SIZE(LA%loc_to_glob, 1) * mesh%gauss%n_w) :: v_loc
      INTEGER, DIMENSION(SIZE(LA%loc_to_glob, 1) * mesh%gauss%n_w) :: idxm
      INTEGER :: i, m, l, ni, iglob, nw, ki, ix, kmax
      REAL(KIND = 8) :: fl
      Vec                                         :: vect
      PetscErrorCode                              :: ierr
      CALL VecSet(vect, 0.d0, ierr)

      kmax = SIZE(LA%loc_to_glob, 1)
      nw = mesh%gauss%n_w

      DO m = 1, mesh%me
         jj_loc = mesh%jj(:, m)
         DO ki = 1, kmax
            !ff_loc((ki - 1) * nw + 1:ki * nw) = ff(ki, jj_loc)
            ff_loc((ki - 1) * nw + 1:ki * nw) = ff(jj_loc,ki)
            DO ni = 1, nw
               i = jj_loc(ni)
               iglob = LA%loc_to_glob(ki, i)
               ix = (ki - 1) * nw + ni
               idxm(ix) = iglob - 1
            END DO
         END DO

         v_loc = 0.d0
         DO l = 1, mesh%gauss%l_G
            DO ki = 1, kmax
               fl = SUM(ff_loc((ki - 1) * nw + 1:ki * nw) * mesh%gauss%ww(:, l)) * mesh%gauss%rj(l, m)
               DO ni = 1, nw
                  ix = (ki - 1) * nw + ni
                  v_loc(ix) = v_loc(ix) + mesh%gauss%ww(ni, l) * fl
               END DO
            END DO
         ENDDO

         CALL VecSetValues(vect, kmax * nw, idxm, v_loc, ADD_VALUES, ierr)
      ENDDO
      CALL VecAssemblyBegin(vect, ierr)
      CALL VecAssemblyEnd(vect, ierr)
   END SUBROUTINE qs_00_block

   SUBROUTINE qs_11 (mesh, LA, ff, vect)
      !> transfer of vect(PETSc projected on test functions) <== grad(ff)(HYPAR)
      !! mesh(hypar)
      !! LA(petsc_csr_MLA)
      !=================================
      USE def_type_mesh
      USE petsc_csr_LA_module
      IMPLICIT NONE
      TYPE(mesh_type)         :: mesh
      type(petsc_csr_LA)      :: LA
      REAL(KIND = 8), DIMENSION(:), INTENT(IN)  :: ff
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w) :: ff_loc
      INTEGER, DIMENSION(mesh%gauss%n_w)        :: jj_loc
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w) :: v_loc
      INTEGER, DIMENSION(mesh%gauss%n_w)        :: idxm
      INTEGER        :: i, k, m, l, ni, iglob
      REAL(KIND = 8) :: dfl(mesh%gauss%k_d)
      Vec                                         :: vect
      PetscErrorCode                              :: ierr
      
      CALL VecZeroEntries(vect, ierr)
      
      DO m = 1, mesh%me
         jj_loc = mesh%jj(:, m)
         ff_loc = ff(jj_loc)
         DO ni = 1, mesh%gauss%n_w
            i = mesh%jj(ni, m)
            iglob = LA%loc_to_glob(1, i)
            idxm(ni) = iglob - 1
         END DO
         
         v_loc = 0.d0
         DO l = 1, mesh%gauss%l_G
            DO k = 1, mesh%gauss%k_d
               dfl(k) = SUM(ff_loc * mesh%gauss%dw(k, :, l, m))
            END DO
            dfl = dfl * mesh%gauss%rj(l, m)
            DO ni = 1, mesh%gauss%n_w
               v_loc(ni) = v_loc(ni) + SUM(mesh%gauss%dw(:, ni, l, m) * dfl)
            END DO
         ENDDO
         CALL VecSetValues(vect, mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
      ENDDO
      CALL VecAssemblyBegin(vect, ierr)
      CALL VecAssemblyEnd(vect, ierr)
   END SUBROUTINE qs_11
   
END MODULE fem_rhs
