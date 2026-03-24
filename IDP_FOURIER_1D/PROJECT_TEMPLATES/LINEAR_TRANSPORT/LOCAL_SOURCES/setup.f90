MODULE setup_module
  TYPE setup_type
   CONTAINS
     PROCEDURE, PUBLIC :: flux=>linear_flux
     PROCEDURE, PUBLIC :: lambda_max=>linear_lambda_max  
  END TYPE setup_type
CONTAINS
  FUNCTION linear_flux(this,u) RESULT(vv)
    IMPLICIT NONE
    CLASS(setup_type)             :: this
    REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: u
    REAL(KIND = 8), DIMENSION(SIZE(u)) :: vv
    vv = u
  END FUNCTION  linear_flux
  FUNCTION linear_lambda_max(this,ul,ur) RESULT(vv)
    IMPLICIT NONE
    CLASS(setup_type)             :: this
    REAL(KIND = 8), INTENT(IN) :: ul, ur
    REAL(KIND = 8) :: vv
    vv = 1.d0
  END FUNCTION linear_lambda_max

  FUNCTION exact_sol(fourier_param,time) RESULT(vv)
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
       xi = (i-1)*fourier_param%dx - time
       k = int(xi/fourier_param%length)
       xi = xi- k*fourier_param%length
       IF (x0<xi .AND. xi<x1) THEN 
          vv(i) = 1.d0
       END IF
    END DO
  END FUNCTION exact_sol
END MODULE setup_module
