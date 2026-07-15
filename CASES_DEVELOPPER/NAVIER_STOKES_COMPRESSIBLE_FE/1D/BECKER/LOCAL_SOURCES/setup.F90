MODULE becker
   REAL(KIND=8), PUBLIC :: becker_vmin, becker_vmax, becker_const, &
         becker_expomax, becker_expomin, becker_x, becker_vorigin
   PUBLIC :: set_becker_parameters, psi, dpsi
   PRIVATE
CONTAINS
   SUBROUTINE set_becker_parameters
      becker_expomax = becker_vmax/(becker_vmax-becker_vmin)
      becker_expomin = -becker_vmin/(becker_vmax-becker_vmin)
   END SUBROUTINE set_becker_parameters
   FUNCTION psi(v) RESULT(vout)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: v
      REAL(KIND=8) :: vout, ratio
      !ratio = ((becker_vmax-v)/(becker_vmax-becker_vorigin))**becker_expomax &
      !     * ((v-becker_vmin)/(becker_vorigin-becker_vmin))**becker_expomin
      !vout = becker_const*LOG(ratio)
      vout = becker_const*(becker_expomax*LOG(becker_vmax-v)              + becker_expomin*LOG(v-becker_vmin) &
                        - becker_expomax*LOG(becker_vmax-becker_vorigin) - becker_expomin*LOG(becker_vorigin-becker_vmin)) - becker_x
   END FUNCTION  psi
   FUNCTION dpsi(v) RESULT(vout)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: v
      REAL(KIND=8) :: vout
      vout = becker_const*(-becker_expomax/(becker_vmax-v) + becker_expomin/(v-becker_vmin))
   END FUNCTION  dpsi
END MODULE becker

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
   REAL(KIND = 8), PARAMETER :: long=1.d0
   REAL(KIND = 8), PARAMETER :: x0=0.5d0*long

   !===Ref J. Fluid Mech. (2013), vol. 726, R4
   REAL(KIND = 8), PARAMETER :: gamma = 1.4d0
   REAL(KIND = 8), PARAMETER :: Mach = 3.d0, rho0=1.d0, v0=1.d0, m0=rho0*v0, cv=1.d0/(gamma-1.d0), cp=gamma*cv, &
                                Pr=0.75d0, Rinfty=(gamma+1)/(gamma-1), uninfty = 0.2d0

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
      USE becker
      IMPLICIT NONE
      TYPE(euler_type), INTENT(INOUT) :: euler
      TYPE(stokes_parabolic_type), INTENT(INOUT) :: stokes
      REAL(KIND = 8),   INTENT(IN)    :: thermal_diffusivity

      !=== Init Euler
      euler%bc%gamma = gamma
      euler%pressure => pressure
      
      euler%bc%mt_anal    => mt_anal_rho_times_vit
      euler%bc%E_anal     => E_anal_ideal_gas

      euler%bc%rho_anal   => rho_anal_becker
      euler%bc%press_anal => press_anal_becker
      euler%bc%vit_anal   => vit_anal_becker_euler

      !euler%eta_commute => default_eta_commute
      euler%eta_commute => eta_commute
      !=== Init Euler

      !=== Init Stokes
      stokes%temperature  => temperature !<=== FIXME
      stokes%bc%vit_anal  => vit_anal_becker
      stokes%bc%temp_anal => temp_anal_becker
      !=== Init Stokes

      !=== Init Becker
      becker_vmax=v0
      becker_vmin=becker_vmax*(gamma-1+2/Mach**2)/(gamma+1) !===(2.10)
      becker_vorigin = sqrt(becker_vmax*becker_vmin) !===Origin at adiabatic sonic point ((3.6) to (3.7))
      becker_const=(2.d0/(gamma+1))*thermal_diffusivity/(m0*cv) !===(3.4) and (3.6)
      CALL set_becker_parameters
      !=== Init Becker


   END SUBROUTINE init_state_functions

   FUNCTION eta_commute(un) RESULT(eta)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:, :), INTENT(IN)    :: un
      REAL(KIND=8), DIMENSION(SIZE(un,1))         :: eta
      eta = (un(:,k_dim+2) - 0.5d0*SUM(un(:,2:k_dim+1)**2,dim=2)/un(:,1))/un(:,1)**gamma
   END FUNCTION eta_commute

   FUNCTION temp_anal_becker(time, rr) RESULT(vv)
      USE my_util, ONLY: local_error_petsc
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      
      vv = 0.d0
      IF (SIZE(rr)>0) THEN
         CALL local_error_petsc('Bug in stokes bc, temp_anal: not defined')
      END IF

   END FUNCTION temp_anal_becker


   FUNCTION rho_anal_becker(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      vv = m0/(this%vit_anal(1,time,rr)-uninfty)
      !===dummy to avoid warning in compilation===!
      RETURN
      xi = this%gamma 
      !===dummy to avoid warning in compilation===!
   END FUNCTION rho_anal_becker

   FUNCTION press_anal_becker(this, time, rr) RESULT(vv)
      USE becker
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: zz
      INTEGER :: n
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      zz = this%vit_anal(1,time,rr)-uninfty
      vv = (m0/zz)*(Rinfty*becker_vmin*becker_vmax - zz**2)/(2*cp) !===(3.7)
      !===dummy to avoid warning in compilation===!
      RETURN
      xi = this%gamma 
      !===dummy to avoid warning in compilation===!
   END FUNCTION press_anal_becker

   FUNCTION vit_anal_becker_euler(this, comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = vit_anal_becker(comp, time, rr)
   END FUNCTION vit_anal_becker_euler
   
   FUNCTION vit_anal_becker(comp, time, rr) RESULT(vv)
      USE becker
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      REAL(KIND = 8) :: xi, vold
      REAL(KIND = 8), PARAMETER  :: tol=1.d-10
      IF (SIZE(vv)==0) RETURN
      vv(1) = 0.99*becker_vmax+0.01*becker_vmin
      DO n = 1, SIZE(rr,2)
         IF (n>1) THEN
            vv(n) = vold
         END IF
         becker_x = rr(1,n)-uninfty*time
         CALL newton_iter(vv(n),becker_vmin,becker_vmax,psi,dpsi,tol)
         vold = vv(n)
         vv(n) = vv(n)+uninfty
      END DO
      !===dummy to avoid warning in compilation===!
      RETURN
      n = comp
      !===dummy to avoid warning in compilation===!
   CONTAINS
      SUBROUTINE newton_iter(v,vmin,vmax,psi_func,dpsi_func,tol)
         IMPLICIT NONE
         INTERFACE
            FUNCTION psi_func(vin) RESULT(vout)
               IMPLICIT NONE
               REAL(KIND=8), INTENT(IN) :: vin
               REAL(KIND=8) :: vout
            END FUNCTION psi_func
            FUNCTION dpsi_func(vin) RESULT(vout)
               IMPLICIT NONE
               REAL(KIND=8), INTENT(IN) :: vin
               REAL(KIND=8) :: vout
            END FUNCTION dpsi_func
         END INTERFACE
         REAL(KIND=8), INTENT(IN) :: vmin, vmax, tol
         REAL(KIND=8) :: v
         REAL(KIND=8) :: norm, psi, vnext, dp
         INTEGER :: iter
         norm = MAX(ABS(psi_func(0.25*vmin+0.75*vmax)), ABS(psi_func(0.75*vmin+0.25*vmax)))
         !===Test trivial cases
         IF (psi_func(vmin+tol).LE.tol) THEN
            v = vmin+tol
            RETURN
         ELSE IF (psi_func(vmax-tol).GE.tol) THEN
            v=vmax-tol
            RETURN
         END IF
         
         iter = 0
         psi =  psi_func(v)
         DO WHILE (ABS(psi).GE.tol*norm)
            !dp = (psi_func(v+sqrt(tol)) - psi_func(v-sqrt(tol)))/(2*sqrt(tol))
            dp = dpsi_func(v)
            vnext = v - psi/dp
            IF (ABS(vnext-v).LT.tol*(vmax+vmin)/2) THEN
               RETURN
            ELSE IF (vnext<vmin) THEN
               v = (vmin+v)/2
            ELSE IF (vnext>vmax) THEN
               v = (vmax+v)/2
            ELSE
               v = vnext
            END IF
            psi = psi_func(v)
            iter = iter+1
         END DO
      END SUBROUTINE newton_iter
   END FUNCTION vit_anal_becker

END MODULE setup