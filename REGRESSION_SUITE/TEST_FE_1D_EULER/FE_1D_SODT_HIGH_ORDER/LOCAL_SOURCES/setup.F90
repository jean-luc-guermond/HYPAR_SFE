MODULE setup
   USE space_dim, ONLY : k_dim
   USE euler_bc_arrays, ONLY: euler_bc_type, mt_anal_rho_times_vit, E_anal_ideal_gas

   PUBLIC :: pressure, init_state_functions

   PRIVATE
   REAL(KIND = 8), PARAMETER :: long=1.d0
   REAL(KIND = 8), PARAMETER :: x0=0.5d0*long

   REAL(KIND = 8), PARAMETER :: gamma = 1.4d0
   REAL(KIND = 8), PARAMETER :: rhoL=1.d0, rhor=0.125d0, pl=1.d0, pr=0.1d0, ul=0.d0, ur=0.d0,&
        l1m=-1.183215956619923d0, l1p=-0.07027281256118334d0, &
        l3=1.7521557320301779, ustar=0.92745262004894991d0, rhoLstar=0.4263194281784952d0, &
        rhoRstar=0.26557371170530708d0, pstar=0.3031301780506468, cL=SQRT(gamma*pL/rhoL) 

CONTAINS

!==========================================================================
!================= DEF PRESSURE FOR SETUP =================================
!==========================================================================

   FUNCTION pressure(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      REAL(KIND = 8) :: gamma
      gamma = 7.0 / 5.0
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

      bc%rho_anal   => rho_anal_sodt
      bc%press_anal => press_anal_sodt
      bc%vit_anal   => vit_anal_sodt
   END SUBROUTINE init_state_functions

   FUNCTION rho_anal_sodt(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (time.LE.1.d-12) THEN
            IF (rr(1, n)<x0) THEN
               vv(n) = rhol
            ELSE
               vv(n) = rhor
            END IF
         ELSE
            xi = (rr(1,n)-x0)/time
            IF (xi.LE.l1m) THEN
               vv(n) = rhoL
            ELSE IF (xi.LE.l1p) THEN
               vv(n) = rhoL*(2/(gamma+1) + (uL-xi)*(gamma-1)/((gamma+1)*cL))**(2/(gamma-1))
            ELSE IF (xi.LE.ustar) THEN
               vv(n) = rhoLstar
            ELSE IF (xi.LE.l3) THEN
               vv(n) = rhoRstar
            ELSE
               vv(n) = rhoR
            END IF
         END IF
      END DO
   END FUNCTION rho_anal_sodt

   FUNCTION press_anal_sodt(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (time.LE.1.d-12) THEN
            IF (rr(1, n)<x0) THEN
               vv(n) = pl
            ELSE
               vv(n) = pr
            END IF
         ELSE
            xi = (rr(1,n)-x0)/time
            IF (xi.LE.l1m) THEN
               vv(n) = pL
            ELSE IF (xi.LE.l1p) THEN
               vv(n) = pL*(2/(gamma+1) + (uL-xi)*(gamma-1)/((gamma+1)*cL))**(2*gamma/(gamma-1))
            ELSE IF (xi.LE.l3) THEN
               vv(n) = pstar
            ELSE
               vv(n) = pR
            END IF
         END IF
      END DO
   END FUNCTION press_anal_sodt

   FUNCTION vit_anal_sodt(this, comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (time.LE.1.d-12) THEN
            IF (rr(1, n)<x0) THEN
               vv(n) = ul
            ELSE
               vv(n) = ur
            END IF
         ELSE
            xi = (rr(1,n)-x0)/time
            IF (xi.LE.l1m) THEN
               vv(n) = uL
            ELSE IF (xi.LE.l1p) THEN
               vv(n) = (2/(gamma+1))*(cL + uL*(gamma-1)/2+xi)
            ELSE IF (xi.LE.l3) THEN
               vv(n) = ustar
            ELSE
               vv(n) = uR
            END IF
         END IF
      END DO
   END FUNCTION vit_anal_sodt

END MODULE setup
