MODULE setup
   USE space_dim, ONLY : k_dim
   USE euler_bc_arrays, ONLY: euler_bc_type, mt_anal_rho_times_vit, E_anal_ideal_gas
   USE euler_type_module, ONLY: euler_type
   USE euler_eta_commute, ONLY: default_eta_commute

   PUBLIC :: pressure, init_state_functions

   PRIVATE
   REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0)
   REAL(KIND=8), PARAMETER :: r0=0.15d0, x0=0d0, y0=0.0d0
   REAL(KIND=8), PARAMETER :: u_infty=0.d0, rho_infty=1.d0, p_infty=1.d0, beta0=5.d0, gamma = 1.4d0
   REAL(KIND=8), PARAMETER :: beta=beta0/(2*pi), chi=((gamma-1)/(2*gamma))*beta**2

   CONTAINS
   
!==========================================================================
!================= DEF PRESSURE FOR SETUP =================================
!==========================================================================

   FUNCTION pressure(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      vv = rho * e * (gamma - 1.d0)
   END FUNCTION pressure
   
   FUNCTION eta_commute(un) RESULT(eta)
      USE space_dim
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:, :), INTENT(IN)    :: un
      REAL(KIND=8), DIMENSION(SIZE(un,1))         :: eta, e_tot
   
      e_tot = un(:, k_dim+2) - 0.5d0*SUM(un(:, 2:k_dim+1)**2, DIM=2)/un(:, 1)
      eta = (e_tot * (gamma - 1.d0)) / (un(:, 1))**gamma

   END FUNCTION eta_commute
!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================
      

   SUBROUTINE init_state_functions(euler)
      IMPLICIT NONE
      TYPE(euler_type), INTENT(INOUT) :: euler

      euler%bc%gamma = gamma
      euler%pressure => pressure

      euler%bc%mt_anal      => mt_anal_rho_times_vit
      euler%bc%E_anal       => E_anal_ideal_gas

      euler%bc%rho_anal     => rho_anal_isentropic
      euler%bc%vit_anal     => vit_anal_isentropic
      euler%bc%press_anal   => press_anal_isentropic
      euler%eta_commute => eta_commute
      ! euler%eta_commute => default_eta_commute

   END SUBROUTINE init_state_functions

   FUNCTION rho_anal_isentropic(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv, z
      REAL(kind=8) :: length=2.d0, x_drift
      REAL(KIND=8) :: rsq
      INTEGER :: n, k_x
      DO n = 1, SIZE(rr,2)
         k_x = FLOOR(((rr(1,n)-x0-u_infty*time) + 1)/length)
         x_drift = rr(1,n)-x0-u_infty*time - k_x*length
         rsq = (x_drift)**2 + (rr(2,n)-y0)**2
         ! rsq = (rr(1,n)-x0-u_infty*time)**2 + (rr(2,n)-y0)**2
         z(n) = exp(1-rsq/(r0**2))
      END DO
      vv = (1-chi*z)**(1.d0/(gamma-1.d0))
      !===dummy to avoid warning in compilation===!
      RETURN
      z = this%gamma
      !===dummy to avoid warning in compilation===!
   END FUNCTION rho_anal_isentropic
   
   FUNCTION press_anal_isentropic(this, time,rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
      vv = p_infty*(this%rho_anal(time,rr)/rho_infty)**gamma
   END FUNCTION press_anal_isentropic
   
   FUNCTION vit_anal_isentropic(this, comp,time,rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      INTEGER,                             INTENT(IN) :: comp
      REAL(KIND = 8),                      INTENT(IN) :: time
      REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
      REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
      REAL(KIND=8) :: rsq, z
      REAL(kind=8) :: length=2.d0, x_drift
      INTEGER :: n, k_x
      
      DO n = 1, SIZE(rr,2)
         k_x = FLOOR(((rr(1,n)-x0-u_infty*time) + 1)/length)
         x_drift = rr(1,n)-x0-u_infty*time - k_x*length
         rsq = (x_drift)**2 + (rr(2,n)-y0)**2
         ! rsq = (rr(1,n)-x0-u_infty*time)**2 + (rr(2,n)-y0)**2
         z = exp(0.5d0*(1-rsq/(r0**2)))
      
         IF (comp==1) THEN
            vv(n) = u_infty - beta*z*(rr(2,n)-y0)/r0
         ELSE
            vv(n) = beta*z*(x_drift)/r0
         END IF
      
      END DO
      
      !===dummy to avoid warning in compilation===!
      RETURN
      z = this%gamma
      !===dummy to avoid warning in compilation===!
   END FUNCTION vit_anal_isentropic
   
END MODULE setup
