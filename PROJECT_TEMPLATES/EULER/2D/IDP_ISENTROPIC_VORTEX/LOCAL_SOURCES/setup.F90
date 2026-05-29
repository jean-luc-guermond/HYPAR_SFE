MODULE setup
   USE space_dim, ONLY : k_dim
   USE euler_bc_arrays, ONLY: euler_bc_type, mt_anal_rho_times_vit, E_anal_ideal_gas

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
      vv = rho * e * (gamma - 1)
   END FUNCTION pressure
   
!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================
      

   SUBROUTINE init_state_functions(bc)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: bc

      bc%gamma = gamma

      bc%mt_anal      => mt_anal_rho_times_vit
      bc%E_anal       => E_anal_ideal_gas

      bc%rho_anal     => rho_anal_isentropic
      bc%vit_anal     => vit_anal_isentropic
      bc%press_anal   => press_anal_isentropic

   END SUBROUTINE init_state_functions

   FUNCTION rho_anal_isentropic(this, time,rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv, z
      REAL(KIND=8) :: rsq
      INTEGER :: n
      DO n = 1, SIZE(rr,2)
         rsq = (rr(1,n)-x0-u_infty*time)**2 + (rr(2,n)-y0)**2
         z(n) = exp(1-rsq/(r0**2))
      END DO
      vv = (1-chi*z)**(1.d0/(gamma-1.d0))
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
      REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv, z
      REAL(KIND=8) :: rsq
      INTEGER :: n
      
      DO n = 1, SIZE(rr,2)
         rsq = (rr(1,n)-x0-u_infty*time)**2 + (rr(2,n)-y0)**2
         z(n) = exp(0.5d0*(1-rsq/(r0**2)))
      END DO
      
      IF (comp==1) THEN
         vv = u_infty - beta*z*(rr(2,:)-y0)/r0
      ELSE
         vv = beta*z*(rr(1,:)-x0-u_infty*time)/r0
      END IF
   END FUNCTION vit_anal_isentropic
   
END MODULE setup
