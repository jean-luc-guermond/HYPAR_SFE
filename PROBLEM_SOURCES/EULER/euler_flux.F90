MODULE euler_flux
   ! USE space_dim
   USE mesh_parameters
   USE eos

   PUBLIC :: flux
CONTAINS
   FUNCTION flux(comp, un) RESULT(vv)  
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(SIZE(un, 1), mesh_data_info%k_dim) :: vv
      REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: H, u, ie
      INTEGER :: k

      ! SELECT CASE(comp)
      ! CASE(1)
      IF (comp == 1) THEN
         DO k = 1, mesh_data_info%k_dim
            vv(:, k) = un(:, k + 1)
         END DO
      ! CASE(2:k_dim + 1)
      ELSE IF ((comp>=2) .AND. (comp<=mesh_data_info%k_dim+1)) THEN
         u = un(:, comp) / un(:, 1)
         DO k = 1, mesh_data_info%k_dim
            vv(:, k) = un(:, k + 1) * u
         END DO
         ie = un(:, mesh_data_info%k_dim + 2) / un(:, 1)
         DO k = 1, mesh_data_info%k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO
         vv(:, comp - 1) = vv(:, comp - 1) + pressure(un(:, 1), ie)
      ! CASE(mesh_data_info%k_dim + 2)
      ELSE IF (comp == mesh_data_info%k_dim + 2) THEN
         ie = un(:, mesh_data_info%k_dim + 2) / un(:, 1)
         DO k = 1, mesh_data_info%k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO

         H = un(:, comp) + pressure(un(:, 1), ie)
         DO k = 1, mesh_data_info%k_dim
            vv(:, k) = (un(:, k + 1) / un(:, 1)) * H
         END DO
      ! CASE DEFAULT
      ELSE
         WRITE(*, *) ' BUG in flux'
         STOP
      ! END SELECT
      END IF
   END FUNCTION flux
END MODULE euler_flux
