MODULE setup
   USE space_dim, ONLY : k_dim
   USE euler_bc_arrays, ONLY: euler_bc_type, mt_anal_rho_times_vit, E_anal_ideal_gas, scal_one

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
   
   SUBROUTINE init_state_functions(bc)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: bc

      bc%gamma = gamma

      bc%mt_anal    => mt_anal_rho_times_vit
      bc%E_anal     => E_anal_ideal_gas
      bc%press_anal => scal_one

      bc%rho_anal   => rho_anal_wind_tunnel
      bc%vit_anal   => vit_anal_wind_tunnel
   END SUBROUTINE init_state_functions

   FUNCTION rho_anal_wind_tunnel(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      IF (SIZE(vv)==0) RETURN
      vv = gamma
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
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      IF (SIZE(vv)==0) RETURN
      IF (comp==1) THEN
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
      ELSE IF (comp==2) THEN
         vv = 0.d0
      ELSE
         WRITE(*, *) ' BUG '
         STOP
      END IF
   END FUNCTION vit_anal_wind_tunnel
   
   ! FUNCTION E_anal(this, time, rr) RESULT(vv)
   !    IMPLICIT NONE
   !    CLASS(euler_bc_type), INTENT(INOUT) :: this
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
   !    REAL(KIND = 8), INTENT(IN) :: time
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
   !    vv = this%press_anal(time, rr) / (gamma - 1.d0) &
   !    + this%rho_anal(time, rr) * (this%vit_anal(1, time, rr)**2 + this%vit_anal(2, time, rr)**2) / 2
   ! END FUNCTION E_anal
   
   ! FUNCTION mt_anal(this, comp, time, rr) RESULT(vv)
   !    IMPLICIT NONE
   !    CLASS(euler_bc_type), INTENT(INOUT) :: this
   !    INTEGER, INTENT(IN) :: comp
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
   !    REAL(KIND = 8), INTENT(IN) :: time
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
   !    vv = this%rho_anal(time, rr) * this%vit_anal(comp, time, rr)
   ! END FUNCTION mt_anal
   
END MODULE setup
