MODULE setup_module
  !TYPE setup_type
  ! CONTAINS
  !    PROCEDURE, PUBLIC, NOPASS :: flux=>linear_flux
  !   PROCEDURE, PUBLIC, NOPASS :: lambda_max=>linear_lambda_max
  !END TYPE setup_type
CONTAINS

  FUNCTION flux(u) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: u
    REAL(KIND = 8), DIMENSION(SIZE(u)) :: vv
    vv = sin(u)
  END FUNCTION  flux

  FUNCTION lambda_max(ul,ur) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND = 8), INTENT(IN) :: ul, ur
    REAL(KIND = 8) :: vv
    vv = 1.d0
  END FUNCTION lambda_max

  FUNCTION exact_sol_step_R(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
    REAL(KIND=8) :: time, xi
    REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0), a=1.d0, b=0.d0
    INTEGER :: i
    DO i = 1, Fourier_param%NMax_real
       xi = Fourier_param%rr(i)-Fourier_param%Length/2
       IF (xi<time*COS((2+a)*pi)) THEN
          vv(i) = (2+a)*pi
       ELSE IF (xi<0) THEN
          vv(i) = 3*pi - ACOS(ABS(xi)/time)
       ELSE IF (xi<time*COS(b*pi)) THEN
          vv(i) = ACOS(xi/time)
       ELSE
          vv(i) = b*pi 
       END IF
    END DO
  END FUNCTION exact_sol_step_R

  FUNCTION exact_sol_F(fourier_param,time,exact_sol_R) RESULT(cs_v)
    USE fourier_param_module
    USE fft_1D
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(fourier_param%Nmax,2) :: cs_v
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: r_v
    REAL(KIND=8) :: time
    interface
       function exact_sol_R(fourier_param,time) RESULT(vv)
         USE fourier_param_module
         IMPLICIT NONE
         TYPE(fourier_param_type) :: fourier_param
         REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
         REAL(KIND=8) :: time
       END FUNCTION exact_sol_R
    end interface
    r_v = exact_sol_R(fourier_param,time)
    CALL real_to_fourier(r_v,cs_v)
  END FUNCTION exact_sol_F

END MODULE setup_module
