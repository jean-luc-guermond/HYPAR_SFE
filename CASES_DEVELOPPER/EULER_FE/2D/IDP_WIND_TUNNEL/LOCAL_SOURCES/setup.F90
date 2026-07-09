MODULE setup
   USE space_dim, ONLY : k_dim
   USE euler_bc_arrays, ONLY: euler_bc_type, mt_anal_rho_times_vit, E_anal_ideal_gas, scal_one
   USE euler_type_module, ONLY: euler_type
   USE euler_eta_commute, ONLY: default_eta_commute

   PUBLIC :: pressure, init_state_functions

   PRIVATE
   REAL(KIND=8), PARAMETER :: r0=0.15d0, x0=0d0, y0=0.0d0
   REAL(KIND=8), PARAMETER :: u_infty=0.d0, rho_infty=1.d0, p_infty=1.d0, beta0=5.d0, gamma = 1.4d0
     
   CONTAINS
   
!==========================================================================
!================= DEF PRESSURE FOR SETUP =================================
!==========================================================================

   FUNCTION pressure(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      vv = rho * e * (gamma - 1)
   END FUNCTION pressure

!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================
   
   SUBROUTINE init_state_functions(euler)
      IMPLICIT NONE
      TYPE(euler_type), INTENT(INOUT) :: euler

      euler%bc%gamma = gamma
      euler%pressure => pressure

      euler%bc%mt_anal    => mt_anal_rho_times_vit
      euler%bc%E_anal     => E_anal_ideal_gas
      euler%bc%press_anal => scal_one

      euler%bc%rho_anal   => rho_anal_wind_tunnel
      euler%bc%vit_anal   => vit_anal_wind_tunnel
      
      euler%eta_commute => default_eta_commute
   END SUBROUTINE init_state_functions

   FUNCTION rho_anal_wind_tunnel(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8),                  INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      IF (SIZE(vv)==0) RETURN
      vv = gamma
      !===dummy to avoid warning in compilation===!
      RETURN
      vv = this%gamma; vv = time
      !===dummy to avoid warning in compilation===!
   END FUNCTION rho_anal_wind_tunnel
   
   ! FUNCTION press_anal(this, time, rr) RESULT(vv)
   !    IMPLICIT NONE
   !    CLASS(euler_bc_type), INTENT(INOUT) :: this
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
   !    REAL(KIND = 8), INTENT(IN) :: time
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
   !    IF (SIZE(vv)==0) RETURN
   !    vv = 1.d0
   ! END FUNCTION press_anal
   
   FUNCTION vit_anal_wind_tunnel(this, comp, time, rr) RESULT(vv)
      USE my_util, ONLY: error_petsc, to_str
      IMPLICIT NONE
      CLASS(euler_bc_type),         INTENT(INOUT) :: this
      INTEGER,                         INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8),                  INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      INTEGER :: n
      IF (SIZE(vv)==0) RETURN
      SELECT CASE(comp)
      CASE(1)
         IF (time<1.d-8) THEN
            vv = 3.d0
            RETURN
         END IF
         DO n = 1, SIZE(vv)
            IF (rr(1, n)<1.d-8) THEN
               vv(n) = 3.0
            ELSE
               vv(n) = 0.d0
            END IF
         END DO
      CASE(2)
         vv = 0.d0
      CASE DEFAULT
         CALL error_petsc("BUG in vit_anal_wind_tunnel: wrong component "//to_str(comp))
      END SELECT
      !===dummy to avoid warning in compilation===!
      RETURN
      vv = this%gamma
      !===dummy to avoid warning in compilation===!
   END FUNCTION vit_anal_wind_tunnel
   
END MODULE setup
