MODULE setup_module
CONTAINS

  FUNCTION flux(u) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: u
    REAL(KIND = 8), DIMENSION(SIZE(u)) :: vv
    vv = u**2/2
  END FUNCTION  flux
  
  FUNCTION flux_prime(u) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: u
    REAL(KIND = 8), DIMENSION(SIZE(u)) :: vv
    vv = u
  END FUNCTION  flux_prime
  
  FUNCTION lambda_max(ul,ur) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND = 8), INTENT(IN) :: ul, ur
    REAL(KIND = 8) :: vv
    IF (ul>ur) THEN
       vv = abs(ul+ur)/2
    ELSE
       vv = max(abs(ul),abs(ur))
    END IF
  END FUNCTION lambda_max

  FUNCTION exact_sol_step_R(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
    REAL(KIND=8) :: time
    REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0)
    vv = SIN(2*pi*Fourier_param%rr/Fourier_param%Length)
  END FUNCTION exact_sol_step_R

  FUNCTION exact_sol_step_R1(fourier_param,time) RESULT(vv)
    USE fourier_param_module
    IMPLICIT NONE
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(fourier_param%Nmax_real) :: vv
    REAL(KIND=8) :: time
    REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0)
    INTEGER :: i
    DO i = 1, Fourier_param%NMax_real
       IF (Fourier_param%rr(i)< Fourier_param%Length/2) THEN
          vv(i) = Fourier_param%rr(i)/(1+time)
       ELSE
          vv(i) = (Fourier_param%rr(i)-Fourier_param%Length)/(1+time)
       END IF
    END DO
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
