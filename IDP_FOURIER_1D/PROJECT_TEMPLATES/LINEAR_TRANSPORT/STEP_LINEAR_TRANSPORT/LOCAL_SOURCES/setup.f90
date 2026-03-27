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
    vv = u
  END FUNCTION  flux

  FUNCTION flux_prime(u) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: u
    REAL(KIND = 8), DIMENSION(SIZE(u)) :: vv
    vv = 1
  END FUNCTION  flux_prime
  
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
    REAL(KIND=8) :: time, xi, x0, x1, delta
    INTEGER :: i, k
    delta = fourier_param%length/20
    x0= fourier_param%length/2 - delta
    x1= fourier_param%length/2 + delta
    vv = 0.d0
    DO i = 1, fourier_param%Nmax_real
       xi =  fourier_param%rr(i) - time
       k = floor(xi/fourier_param%length)
       xi = xi- k*fourier_param%length
       IF (x0<xi .AND. xi<x1) THEN 
          vv(i) = 1.d0
       END IF
    END DO
  END FUNCTION exact_sol_step_R

  FUNCTION exact_sol_step_R1(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
    REAL(KIND=8) :: time
    REAL(KIND=8), PARAMETER :: pi=ACOS(-1d0)
    vv = SIN(2*pi*(fourier_param%rr - time)*fourier_param%length)
  END FUNCTION exact_sol_step_R1

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
