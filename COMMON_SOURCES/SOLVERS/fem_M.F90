#include "petsc/finclude/petsc.h"
MODULE fem_M
   USE petsc_csr_LA_module
   USE def_type_mesh
   USE petsc

   PUBLIC :: qs_mass_diff_M, qs_mass_block_M, qs_var_mass_block_M
   PRIVATE
CONTAINS

   SUBROUTINE qs_mass_diff_M (mesh, mass, visco, LA, matrix)
      !=================================================
      IMPLICIT NONE
      TYPE(mesh_type), TARGET :: mesh
      REAL(KIND = 8), INTENT(IN) :: mass, visco
      type(petsc_csr_LA) :: LA
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w, mesh%gauss%n_w) :: mat_loc

      INTEGER, DIMENSION(mesh%gauss%n_w) :: idxn
      INTEGER :: m, ni, nj, k
      REAL(KIND = 8), DIMENSION(mesh%gauss%l_G) :: al, bl
      Mat            :: matrix
      PetscErrorCode :: ierr
      CALL MatZeroEntries (matrix, ierr)

      DO m = 1, mesh%me
         idxn = LA%loc_to_glob(1, mesh%jj(:, m)) - 1

         al = visco * mesh%gauss%rj(:, m)
         bl = mass * mesh%gauss%rj(:, m)
         DO nj = 1, mesh%gauss%n_w;
            DO ni = 1, mesh%gauss%n_w;
               mat_loc(nj, ni) = SUM(mesh%gauss%ww(ni, :) * mesh%gauss%ww(nj, :) * bl)
               DO k = 1, mesh%gauss%k_d
                  mat_loc(nj, ni) = mat_loc(nj, ni) + SUM(mesh%gauss%dw(k, nj, :, m) * mesh%gauss%dw(k, ni, :, m)  * al)
               END DO
            ENDDO
         ENDDO
         CALL MatSetValues(matrix, mesh%gauss%n_w, idxn, mesh%gauss%n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO

      CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
   END SUBROUTINE qs_mass_diff_M

   SUBROUTINE qs_mass_block_M (mesh, mass, LA, matrix)
      !> simple mass matrix by blocks of k_dim
      !! (IN): mesh, LA
      !! (IN): mass
      !! (OUT): matrix
      !! WARNING ON DIMENSION:
      !!     - LA%loc_to_glob(k_max, mesh%np) where k_max is arbitrary (say 1, k_dim, 3, etc...)
      !!     - mass: scalar
      !!     - matrix : (k_max*np) x (k_max*np) is filled by blocks of mass matrices along its diagonal
      !=================================================
      USE petscmat
      USE space_dim
      IMPLICIT NONE
      TYPE(mesh_type)                       :: mesh
      TYPE(petsc_csr_LA) :: LA
      REAL(KIND=8) :: mass
      INTEGER, DIMENSION(mesh%gauss%n_w) :: jj_loc
      INTEGER, DIMENSION(k_dim*mesh%gauss%n_w) :: idxm, idxn
      INTEGER :: m, mi, i, ki, iglob, ix, nj, j, kj, jx, l, n_w, ni, k1, jglob
      REAL(KIND=8) :: x, y
      REAL(KIND=8), DIMENSION(k_dim*mesh%gauss%n_w, k_dim*mesh%gauss%n_w) :: mat_loc
      REAL(KIND=8), DIMENSION(mesh%gauss%l_G) :: bl
      INTEGER    :: ierr
      TYPE(tMat) :: matrix

      CALL MatZeroEntries (matrix, ierr)
      CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      CALL MatSetOption (matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
      n_w = mesh%gauss%n_w
      DO m=1, mesh%me
         jj_loc = mesh%jj(:, m)
         mat_loc = 0.d0
         bl(:) =  mass * mesh%gauss%rj(:, m)
         DO ni = 1, n_w
               i = jj_loc(ni)
               DO ki = 1, k_dim
                  iglob = LA%loc_to_glob(ki, i)
                  ix = (ki - 1) * n_w + ni
                  idxm(ix) = iglob - 1
                  DO nj = 1, n_w
                     j = jj_loc(nj)
                     kj = ki
                     jglob = LA%loc_to_glob(kj, j)
                     jx = (kj - 1) * n_w + nj
                     idxn(jx) = jglob - 1
                     x = SUM(mesh%gauss%ww(ni,:)*mesh%gauss%ww(nj,:)*bl)
                     mat_loc(ix,jx) = x
                  END DO
               END DO
         END DO
         CALL MatSetValues(matrix, k_dim * n_w, idxm, k_dim * n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO

      CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

   END SUBROUTINE qs_mass_block_M

   SUBROUTINE qs_var_mass_block_M (matrix, mesh, LA, mass, rho)
      !> generalization of mass_vel_M: builds mass matrix with mass (constant) and rho (space-dependant)
      !! (IN): mesh, LA
      !! (IN): mass, rho
      !! (OUT): matrix
      !! WARNING ON DIMENSION:
      !!     - LA%loc_to_glob(k_max, mesh%np) where k_max is arbitrary (say 1, k_dim, 3, etc...)
      !!     - mass: scalar
      !!     - rho(mesh%np)
      !!     - matrix : (k_max*np) x (k_max*np) is filled by blocks of mass matrices along its diagonal
      !=================================================
      USE petscmat
      IMPLICIT NONE
      TYPE(mesh_type),    INTENT(IN)    :: mesh
      TYPE(petsc_csr_LA), INTENT(IN)    :: LA
      REAL(KIND = 8),               INTENT(IN) :: mass
      REAL(KIND = 8), DIMENSION(mesh%np), INTENT(IN) :: rho
      REAL(KIND = 8), DIMENSION(:,:), ALLOCATABLE :: mat_loc
      INTEGER, DIMENSION(:), ALLOCATABLE :: idxn, idxm

      INTEGER, DIMENSION(mesh%gauss%n_w) :: jj_loc
      INTEGER :: m, ni, nj, ki, kj, l, k_max, n_w
      REAL(KIND = 8), DIMENSION(mesh%gauss%l_G) :: bl
      TYPE(tMat) :: matrix
      INTEGER :: ierr
      INTEGER :: i, iglob, ix, j, jglob, jx

      CALL MatZeroEntries (matrix, ierr)
      k_max = SIZE(LA%loc_to_glob, 1)
      n_w = mesh%gauss%n_w
      ALLOCATE(mat_loc(k_max*n_w, k_max*n_w))
      ALLOCATE(idxn(k_max*n_w))
      ALLOCATE(idxm(k_max*n_w))

      DO m = 1, mesh%me
         jj_loc = mesh%jj(:, m)
         mat_loc = 0.d0
         DO l = 1, mesh%gauss%l_G
               bl(l) =  SUM(rho(mesh%jj(:,m))*mesh%gauss%ww(:,l)) * mass * mesh%gauss%rj(l, m)
         END DO
         DO ni=1, n_w
               i = jj_loc(ni)
               DO ki=1, k_max
                  iglob = LA%loc_to_glob(ki, i)
                  ix = (ki - 1) * n_w + ni
                  idxm(ix) = iglob - 1
                  DO nj = 1, n_w
                     j = jj_loc(nj)
                     kj = ki
                     jglob = LA%loc_to_glob(kj, j)
                     jx = (kj - 1) * n_w + nj
                     idxn(jx) = jglob - 1
                     mat_loc(ix,jx) = SUM(bl*mesh%gauss%ww(ni,:)*mesh%gauss%ww(nj,:))
                  END DO
               END DO
         END DO
         CALL MatSetValues(matrix, k_max*n_w, idxm, k_max*n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO

      CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
   END SUBROUTINE qs_var_mass_block_M

   SUBROUTINE inverse_mat(a, n)       ! Invert matrix by Gauss method
      ! --------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: n
      REAL(8) :: a(n, n)

      ! - - - Local Variables - - -
      REAL(8) :: b(n, n), c, d, temp(n)
      INTEGER :: i, j, k, m, imax(1), ipvt(n)
      ! - - - - - - - - - - - - - -

      b = a
      ipvt = (/ (i, i = 1, n) /)

      DO k = 1, n
         imax = MAXLOC(ABS(b(k:n, k)))
         m = k - 1 + imax(1)

         IF (m /= k) THEN
            ipvt((/m, k/)) = ipvt((/k, m/))
            b((/m, k/), :) = b((/k, m/), :)
         END IF
         d = 1 / b(k, k)

         temp = b(:, k)
         DO j = 1, n
            c = b(k, j) * d
            b(:, j) = b(:, j) - temp * c
            b(k, j) = c
         END DO
         b(:, k) = temp * (-d)
         b(k, k) = d
      END DO

      a(:, ipvt) = b

   END SUBROUTINE inverse_mat
END MODULE fem_M
