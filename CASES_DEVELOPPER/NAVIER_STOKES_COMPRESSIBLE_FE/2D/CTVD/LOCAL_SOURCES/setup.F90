MODULE setup
   USE space_dim, ONLY : k_dim
   USE euler_bc_arrays, ONLY: euler_bc_type
   USE euler_bc_arrays, ONLY: mt_anal_rho_times_vit, E_anal_ideal_gas
   USE euler_type_module, ONLY: euler_type
   USE euler_eta_commute, ONLY: default_eta_commute
   USE stokes_parabolic_module
   USE stokes_bc_arrays
   PUBLIC :: pressure, init_state_functions

   PRIVATE
   !===Ref J. Fluid Mech. (2013), vol. 726, R4
   REAL(KIND = 8), PARAMETER :: gamma = 1.4d0, cv=1.d0/(gamma-1.d0)
   REAL(KIND=8) :: x0=0.50000000001d0
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

   FUNCTION temperature(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      vv = cv * e
   END FUNCTION temperature 

!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================
   
   SUBROUTINE init_state_functions(euler, stokes, thermal_diffusivity)
      IMPLICIT NONE
      TYPE(euler_type), INTENT(INOUT) :: euler
      TYPE(stokes_parabolic_type), INTENT(INOUT) :: stokes
      REAL(KIND = 8),   INTENT(IN)    :: thermal_diffusivity

      !=== Init Euler
      euler%bc%gamma = gamma
      euler%pressure => pressure
      
      euler%bc%mt_anal    => mt_anal_rho_times_vit
      euler%bc%E_anal     => E_anal_ideal_gas

      euler%bc%rho_anal   => rho_anal
      euler%bc%press_anal => press_anal
      euler%bc%vit_anal   => vit_anal_euler

      !euler%eta_commute => default_eta_commute
      euler%eta_commute => eta_commute
      !=== Init Euler

      !=== Init Stokes
      stokes%temperature  => temperature !<=== FIXME
      stokes%bc%vit_anal  => vit_anal_stokes
      stokes%bc%temp_anal => temp_anal_stokes
      !=== Init Stokes

   END SUBROUTINE init_state_functions

   FUNCTION eta_commute(un) RESULT(eta)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:, :), INTENT(IN)    :: un
      REAL(KIND=8), DIMENSION(SIZE(un,1))         :: eta
      eta = (un(:,k_dim+2) - 0.5d0*SUM(un(:,2:k_dim+1)**2,dim=2)/un(:,1))/un(:,1)**gamma
   END FUNCTION eta_commute

   FUNCTION temp_anal_stokes(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      
      vv = 0.d0
      IF (SIZE(rr)>0) THEN
         CALL local_error_petsc('Bug in stokes bc, temp_anal: not defined')
      END IF

   END FUNCTION temp_anal_stokes


   FUNCTION rho_anal(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (rr(1,n).LE.x0) THEN
            vv(n) = 120.d0
         ELSE
            vv(n) = 1.2d0
         END IF
      END DO
      !===dummy to avoid warning in compilation===!
      RETURN
      xi = this%gamma 
      !===dummy to avoid warning in compilation===!
   END FUNCTION rho_anal

   FUNCTION press_anal(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: zz
      INTEGER :: n
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (rr(1,n).LE.x0) THEN
            vv(n) = 120.d0/gamma
         ELSE
            vv(n) = 1.2d0/gamma
         END IF
      END DO
      !===dummy to avoid warning in compilation===!
      RETURN
      xi = this%gamma 
      !===dummy to avoid warning in compilation===!
   END FUNCTION press_anal

   FUNCTION vit_anal_euler(this, comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = 0.d0
   END FUNCTION vit_anal_euler
   
   FUNCTION vit_anal_stokes(comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = 0.d0
   END FUNCTION vit_anal_stokes

END MODULE setup