MODULE fem_M
   USE def_type_mesh
#include "petsc/finclude/petsc.h"
   USE petsc
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

      DO m = 1, mesh%dom_me
         idxn = LA%loc_to_glob(1, mesh%jj(:, m)) - 1

         al = visco * mesh%gauss%rj(:, m)
         bl = mass * mesh%gauss%rj(:, m)
         DO nj = 1, mesh%gauss%n_w;
            DO ni = 1, mesh%gauss%n_w;
               mat_loc(nj, ni) = mesh%gauss%ww(ni, :) * mesh%gauss%ww(nj, :) * bl
               DO k = 1, mesh%gauss%k_d
                  mat_loc(nj, ni) = mat_loc(nj, ni) + SUM(mesh%gauss%dw(k, nj, :, m) * mesh%gauss%dw(k, ni, :, m)) * al
               END DO
            ENDDO
         ENDDO

         CALL MatSetValues(matrix, mesh%gauss%n_w, idxn, mesh%gauss%n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO

      CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
   END SUBROUTINE qs_mass_diff_M

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
