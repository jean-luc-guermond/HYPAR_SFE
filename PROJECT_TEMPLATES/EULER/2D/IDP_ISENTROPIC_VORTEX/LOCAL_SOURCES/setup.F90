MODULE setup
   USE mesh_parameters
   PUBLIC :: sol_anal, init, rho_anal, press_anal, mt_anal, E_anal, impose_bc_euler, pressure
   PRIVATE
   REAL(KIND=8), PARAMETER :: r0=1.d0, x0=0d0, y0=0.0d0
   REAL(KIND=8), PARAMETER :: u_infty=0.d0, rho_infty=1.d0, p_infty=1.d0, beta0=5.d0, gamma = 1.4d0
   REAL(KIND=8) :: chi, beta
   CONTAINS

!==========================================================================
!================= DEF PRESSURE FOR SETUP =================================
!==========================================================================

   FUNCTION pressure(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      !REAL(KIND = 8) :: gamma
      !gamma = 7.0 / 5.0
      vv = rho * e * (gamma - 1)
   END FUNCTION pressure

!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================

   SUBROUTINE impose_bc_euler(un, euler_bc, mesh, time)
      USE euler_bc_arrays
      USE def_type_mesh
      TYPE(mesh_type) :: mesh
      TYPE(euler_bc_type) :: euler_bc
      REAL(KIND = 8) :: time
      REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
      INTEGER :: comp
      DO comp = 1, euler_bc%syst_dim
         !  SELECT CASE(comp)
         !  CASE(1)
         IF (comp == 1) THEN
            un(euler_bc%rho_bc%jsd, comp) = rho_anal(time, mesh%rr(:, euler_bc%rho_bc%jsd))
         ELSE IF ((2<=comp) .AND. (comp<=mesh_data_info%k_dim + 1)) THEN
            !
            ! CASE(2:mesh_data_info%k_dim + 1)
            un(euler_bc%rho_bc%jsd, comp) = mt_anal(comp - 1, time, mesh%rr(:, euler_bc%rho_bc%jsd))
         ELSE IF (comp == mesh_data_info%k_dim + 2) THEN
            ! CASE(mesh_data_info%k_dim + 2)
            un(euler_bc%rho_bc%jsd, comp) = E_anal(time, mesh%rr(:, euler_bc%rho_bc%jsd))
         END IF
         ! END SELECT
      END DO
   END SUBROUTINE impose_bc_euler

   SUBROUTINE init(un, time, rr)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), mesh_data_info%k_dim + 2), INTENT(OUT) :: un
      REAL(KIND = 8), INTENT(IN) :: time
      INTEGER :: comp
      REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0)

      beta = beta0/(2*pi)
      chi=((gamma-1)/(2*gamma))*beta**2

      DO comp = 1, mesh_data_info%k_dim+2
         IF (comp == 1) THEN
            un(:, comp) = rho_anal(time, rr)
         ELSE IF ((2<=comp) .AND. (comp<=mesh_data_info%k_dim + 1)) THEN
            un(:, comp) = mt_anal(comp - 1, time, rr)
         ELSE IF (comp == mesh_data_info%k_dim + 2) THEN
            un(:, comp) = E_anal(time, rr)
         END IF
      END DO
   END SUBROUTINE init

   FUNCTION rho_anal(time,rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv, z
      REAL(KIND=8) :: rsq
      INTEGER :: n
      DO n = 1, SIZE(rr,2)
         rsq = (rr(1,n)-x0-u_infty*time)**2 + (rr(2,n)-y0)**2
         z(n) = exp(1-rsq)
      END DO
      vv = (1-chi*z)**(1.d0/(gamma-1.d0))
   END FUNCTION rho_anal

   FUNCTION press_anal(time,rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
      vv = p_infty*(rho_anal(time,rr)/rho_infty)**gamma
   END FUNCTION press_anal

   FUNCTION vit_anal(comp,time,rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER,                             INTENT(IN) :: comp
      REAL(KIND = 8),                      INTENT(IN) :: time
      REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
      REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv, z
      REAL(KIND=8) :: rsq
      INTEGER :: n

      DO n = 1, SIZE(rr,2)
         rsq = (rr(1,n)-x0-u_infty*time)**2 + (rr(2,n)-y0)**2
         z(n) = exp(0.5d0*(1-rsq))
      END DO

      IF (comp==1) THEN
         vv = u_infty - beta*z*(rr(2,:)-y0)
      ELSE
         vv = beta*z*(rr(1,:)-x0-u_infty*time)
      END IF
   END FUNCTION vit_anal

   FUNCTION E_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = press_anal(time, rr) / (gamma - 1.d0) &
      + rho_anal(time, rr) * (vit_anal(1, time, rr)**2) / 2
   END FUNCTION E_anal

   FUNCTION mt_anal(comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = rho_anal(time, rr) * vit_anal(comp, time, rr)
   END FUNCTION mt_anal

   FUNCTION sol_anal(comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      SELECT CASE(comp)
      CASE(1)
         vv = rho_anal(time, rr)
      CASE(2,3)
         vv = mt_anal(comp-1, time, rr)
      CASE(4)
         vv = E_anal(time, rr)
      CASE DEFAULT
         WRITE(*, *) ' BUG in sol_anal'
         STOP
      END SELECT
   END FUNCTION sol_anal
END MODULE setup
