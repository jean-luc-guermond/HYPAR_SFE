MODULE setup_module
  USE eos_module
  PUBLIC ::  init, rho_r, press_r, mt_r, E_r
  REAL(KIND = 8), PARAMETER, PUBLIC :: rhoL=3.d0,  pL=1.d0,                  cL=sqrt(gamma*pL/rhoL), ul=cL
  REAL(KIND = 8), PARAMETER, PUBLIC :: rhor=0.5d0, pR=pL*(rhoR/rhoL)**gamma, cR=sqrt(gamma*pR/rhoR),&
       ur=uL+(2/(gamma-1))*(cL-cR)
  PRIVATE
  REAL(KIND = 8) :: x0, x1
  REAL(KIND = 8), PARAMETER :: l1m=uL-cL, L1p=uR-cR
CONTAINS

  FUNCTION init(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND = 8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(fourier_param%Nmax,2,3) :: vv
    x0 = fourier_param%Length * 0.3d0
    x1 = fourier_param%Length - x0
    vv(:,:,1) = exact_sol_F(fourier_param,time,rho_r)
    vv(:,:,2) = exact_sol_F(fourier_param,time,mt_r)
    vv(:,:,3) = exact_sol_F(fourier_param,time,E_r)
  END FUNCTION  init

  FUNCTION rho_r(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND = 8) :: time
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
    INTEGER :: n
    REAL(KIND = 8) :: xi
    DO n = 1, fourier_param%Nmax_real
       IF (time.LE.1.d-12) THEN
          IF (fourier_param%rr(n)<x0) THEN
             vv(n) = rhol
          ELSE
             vv(n) = rhor
          END IF
       ELSE
          xi = (fourier_param%rr(n)-x0)/time
          IF (xi.LE.l1m) THEN
             vv(n) = rhoL
          ELSE IF (xi.LE.l1p) THEN
             vv(n) = rhoL*(2/(gamma+1) + (uL-xi)*(gamma-1)/((gamma+1)*cL))**(2/(gamma-1))
          ELSE IF (fourier_param%rr(n).LE.x1) THEN
             vv(n) = rhoR
          ELSE
             vv(n) = rhoL
          END IF
       END IF
    END DO
  END FUNCTION rho_r

  FUNCTION press_r(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND = 8) :: time
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
    INTEGER :: n
    REAL(KIND = 8) :: xi
    DO n = 1, SIZE(vv)
       IF (time.LE.1.d-12) THEN
          IF (fourier_param%rr(n)<x0) THEN
             vv(n) = pl
          ELSE
             vv(n) = pr
          END IF
       ELSE
          xi = (fourier_param%rr(n)-x0)/time
          IF (xi.LE.l1m) THEN
             vv(n) = pL
          ELSE IF (xi.LE.l1p) THEN
             vv(n) = pL*(2/(gamma+1) + (uL-xi)*(gamma-1)/((gamma+1)*cL))**(2*gamma/(gamma-1))
          ELSE IF (fourier_param%rr(n).LE.x1) THEN
             vv(n) = pR
          ELSE
             vv(n) = pL
          END IF
       END IF
    END DO
  END FUNCTION press_r

  FUNCTION vit_r(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND = 8) :: time
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
    INTEGER :: n
    REAL(KIND = 8) :: xi
    DO n = 1, SIZE(vv)
       IF (time.LE.1.d-12) THEN
          IF (fourier_param%rr(n)<x0) THEN
             vv(n) = ul
          ELSE
             vv(n) = ur
          END IF
       ELSE
          xi = (fourier_param%rr(n)-x0)/time
          IF (xi.LE.l1m) THEN
             vv(n) = uL
          ELSE IF (xi.LE.l1p) THEN
             vv(n) = (2/(gamma+1))*(cL + uL*(gamma-1)/2+xi)
          ELSE IF (fourier_param%rr(n).LE.x1) THEN
             vv(n) = uR
          ELSE
             vv(n) = uL
          END IF
       END IF
    END DO
  END FUNCTION vit_r

  FUNCTION E_r(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND = 8) :: time
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
    vv = press_r(fourier_param,time) / (gamma - 1.d0) &
         + rho_r(fourier_param,time) * vit_r(fourier_param,time)**2 / 2
  END FUNCTION E_r

  FUNCTION mt_r(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND = 8) :: time
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
    vv = rho_r(fourier_param,time) * vit_r(fourier_param,time)
  END FUNCTION mt_r

  FUNCTION exact_sol_F(fourier_param,time,exact_sol_R) RESULT(cs_v)
    USE fourier_param_module
    USE fft_1D
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(fourier_param%Nmax,2) :: cs_v
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: r_v
    REAL(KIND=8) :: time
    INTERFACE
       FUNCTION exact_sol_R(fourier_param,time) RESULT(vv)
         USE fourier_param_module
         IMPLICIT NONE
         TYPE(fourier_param_type) :: fourier_param
         REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
         REAL(KIND=8) :: time
       END FUNCTION exact_sol_R
    END INTERFACE
    r_v = exact_sol_R(fourier_param,time)
    CALL real_to_fourier(r_v,cs_v)
  END FUNCTION exact_sol_F

END MODULE setup_module
