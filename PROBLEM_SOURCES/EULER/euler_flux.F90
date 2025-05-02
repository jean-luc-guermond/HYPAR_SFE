MODULE euler_flux
   USE space_dim
   USE eos !<===FIX ME
   PUBLIC :: flux
CONTAINS
   FUNCTION flux(comp, un) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(SIZE(un, 1), k_dim) :: vv
      REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: H, u, ie
      INTEGER :: k
      SELECT CASE(comp)
      CASE(1)
         DO k = 1, k_dim
            vv(:, k) = un(:, k + 1)
         END DO
      CASE(2:k_dim + 1)
         u = un(:, comp) / un(:, 1)
         DO k = 1, k_dim
            vv(:, k) = un(:, k + 1) * u
         END DO
         ie = un(:, k_dim + 2) / un(:, 1)
         DO k = 1, k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO
         vv(:, comp - 1) = vv(:, comp - 1) + pressure(un(:, 1), ie)
      CASE(k_dim + 2)
         ie = un(:, k_dim + 2) / un(:, 1)
         DO k = 1, k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO

         H = un(:, comp) + pressure(un(:, 1), ie)
         DO k = 1, k_dim
            vv(:, k) = (un(:, k + 1) / un(:, 1)) * H
         END DO
      CASE DEFAULT
         WRITE(*, *) ' BUG in flux'
         STOP
      END SELECT
   END FUNCTION flux
END MODULE euler_flux
