#include "petsc/finclude/petsc.h"
MODULE fem_rhs
   USE petsc
   CONTAINS
   
   SUBROUTINE qs_00 (mesh, LA, ff, vect)
      !> Projection of ff(Fortran array) onto test functions phi_i
      !! insert into vect(PETSc)
      !! (IN): mesh, LA
      !! (IN): ff(mesh%np)
      !! (OUT): vect(PETSc)
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

   SUBROUTINE qs_00_dot (mesh, LA, ff, gg, vect)
      !> Projection of ff(Fortran array).gg(Fortran array) onto test functions phi_i
      !! insert into vect(PETSc)
      !! (IN): mesh, LA
      !! (IN): ff(mesh%np), g(mesh%np)
      !! (OUT): vect(PETSc)
      !=================================
      USE def_type_mesh
      USE petsc_csr_LA_module
      IMPLICIT NONE
      TYPE(mesh_type)    :: mesh
      type(petsc_csr_LA) :: LA
      REAL(KIND = 8), DIMENSION(:,:), INTENT(IN)  :: ff, gg
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w,SIZE(ff,2)) :: ff_loc, gg_loc
      INTEGER, DIMENSION(mesh%gauss%n_w)        :: jj_loc
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w) :: v_loc
      INTEGER, DIMENSION(mesh%gauss%n_w)        :: idxm
      INTEGER        :: i, m, l, ni, iglob, k, k_max
      REAL(KIND = 8) :: fl
      Vec                                         :: vect
      PetscErrorCode                              :: ierr
      
      CALL VecZeroEntries(vect, ierr)
      k_max = SIZE(ff, 2)
      DO m = 1, mesh%me
         jj_loc = mesh%jj(:, m)
         ff_loc = ff(jj_loc,:)
         gg_loc = gg(jj_loc,:)
         DO ni = 1, mesh%gauss%n_w
            i = mesh%jj(ni, m)
            iglob = LA%loc_to_glob(1, i)
            idxm(ni) = iglob - 1
         END DO
         
         v_loc = 0.d0
         DO l = 1, mesh%gauss%l_G
            fl = 0.d0
            DO k = 1, k_max
               fl = fl + SUM(ff_loc(:,k) * mesh%gauss%ww(:, l)) * SUM(gg_loc(:,k) * mesh%gauss%ww(:, l)) * mesh%gauss%rj(l, m)
            END DO
            DO ni = 1, mesh%gauss%n_w
               v_loc(ni) = v_loc(ni) + mesh%gauss%ww(ni, l) * fl
            END DO
         ENDDO
         CALL VecSetValues(vect, mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
      ENDDO
      CALL VecAssemblyBegin(vect, ierr)
      CALL VecAssemblyEnd(vect, ierr)
   END SUBROUTINE qs_00_dot

   SUBROUTINE qs_00_block (mesh, LA, ff, vect)
      !> Projection of ff(Fortran array) onto test functions
      !! insert into vect(PETSc)
      !! (IN): mesh, LA (with LA%loc_to_glob(k_max, mesh%np))
      !! (IN): ff(mesh%np, k_max)
      !! (OUT): vect(PETSc)
      !=================================
      USE def_type_mesh
      USE petsc_csr_LA_module
      USE petscvec
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
      TYPE(tVec) :: vect
      INTEGER :: ierr
      ! Vec                                         :: vect
      ! PetscErrorCode                              :: ierr

      CALL VecZeroEntries(vect, ierr)

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
      !> Projection of grad(ff)(Fortran array) onto test functions grad(phi_i)
      !! insert into vect(PETSc)
      !! (IN): mesh, LA
      !! (IN): ff(mesh%np)
      !! (OUT): vect(PETSc)
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
