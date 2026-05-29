MODULE setup
   USE space_dim, ONLY : k_dim
   USE euler_type_module, ONLY: euler_type
   USE euler_bc_arrays, ONLY: euler_bc_type
   USE euler_bc_arrays, ONLY: mt_anal_rho_times_vit, E_anal_ideal_gas, scal_one, vect_one
   USE euler_eta_commute, ONLY: default_eta_commute
   PUBLIC :: pressure, init_state_functions
   PRIVATE
   REAL(KIND = 8), PARAMETER, PRIVATE :: x0=0.1d0, x1=0.3d0, gamma=1.4d0

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
      euler%bc%vit_anal   => vect_one

      euler%bc%rho_anal   => rho_anal_smooth

      euler%eta_commute => default_eta_commute
   END SUBROUTINE init_state_functions

   FUNCTION rho_anal_smooth(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF ((rr(1, n) - time)<x0 .OR. (rr(1, n) - time)>x1) THEN
            vv(n) = 1.d0
         ELSE
            vv(n) = 1 + (2 / (x1 - x0))**6 * (rr(1, n) - time - x0)**3 * (x1 - rr(1, n) + time)**3
         END IF
      END DO
   END FUNCTION rho_anal_smooth

END MODULE setup
