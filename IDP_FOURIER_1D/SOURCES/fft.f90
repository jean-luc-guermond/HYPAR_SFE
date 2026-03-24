MODULE fft_1D
  REAL(KIND=8), PRIVATE, PARAMETER :: pi=ACOS(-1.d0)
CONTAINS
  SUBROUTINE real_to_fourier(r_in,cs_out)
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    REAL(KIND=8), DIMENSION(:) :: r_in !===size 2N-1
    COMPLEX(KIND=8), DIMENSION(SIZE(r_in)/2+1)   :: cout !===Complex size N
    REAL(KIND=8),    DIMENSION(SIZE(r_in)/2+1,2) :: cs_out
    REAL(KIND=8) :: scaling
    INTEGER(kind=8) :: plan
    INTEGER :: N, i
    N = SIZE(r_in)
    scaling = 2.d0/N
    CALL dfftw_plan_dft_r2c_1d(plan,N,r_in,cout,FFTW_ESTIMATE)
    CALL dfftw_execute_dft_r2c(plan, r_in, cout)
    CALL dfftw_destroy_plan(plan)
    DO i = 1, (N+1)/2
       cs_out(i,1) = REAL (cout(i),KIND=8)
       cs_out(i,2) = -AIMAG (cout(i))
    END DO
    cs_out = scaling*cs_out 
    cs_out (1,:) =  cs_out(1,:)/2 !===mode 0
  END SUBROUTINE real_to_fourier

  SUBROUTINE fourier_to_real(cs_in,r_out)
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    REAL(KIND=8),    DIMENSION(:,:)            :: cs_in !===(N,2)
    COMPLEX(KIND=8), DIMENSION(SIZE(cs_in,1))  :: cin !===Complex size N
    REAL(KIND=8), DIMENSION(2*SIZE(cs_in,1)-1) :: r_out !===size 2N-1
    REAL(KIND=8) :: scaling
    INTEGER(kind=8) :: plan
    INTEGER :: N
    N = 2*SIZE(cs_in,1)-1
    scaling = 2.d0/N
    cin = 0.5d0*CMPLX(cs_in(:,1),-cs_in(:,2),KIND=8)
    cin(1) = 2*cin(1)
    CALL dfftw_plan_dft_c2r_1d(plan,N,cin,r_out,FFTW_ESTIMATE)
    CALL dfftw_execute_dft_c2r(plan, cin, r_out)
    CALL dfftw_destroy_plan(plan)
  END SUBROUTINE fourier_to_real
  
  SUBROUTINE fourier_derivative(cs_in,cs_out,Length)
    IMPLICIT NONE
    REAL(KIND=8),    DIMENSION(:,:) :: cs_in 
    REAL(KIND=8),    DIMENSION(SIZE(cs_in,1),SIZE(cs_in,2)) :: cs_out
    REAL(KIND=8) :: Length
    INTEGER :: i
    DO i = 1, SIZE(cs_in,1)
       cs_out(i,1) =  (i-1)*cs_in(i,2)
       cs_out(i,2) = -(i-1)*cs_in(i,1)
    END DO
    cs_out = 2*pi*cs_out/Length
  END SUBROUTINE fourier_derivative

  SUBROUTINE real_derivative(r_in,r_out,Length)
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    REAL(KIND=8), DIMENSION(:) :: r_in !===size 2N-1
    REAL(KIND=8), DIMENSION(SIZE(r_in)) :: r_out
    REAL(KIND=8), DIMENSION(SIZE(r_in)/2+1,2) :: cs_in, cs_out
    REAL(KIND=8) :: Length
    CALL real_to_fourier(r_in,cs_in)
    CALL fourier_derivative(cs_in,cs_out,Length)
    CALL fourier_to_real(cs_out,r_out)
  END SUBROUTINE real_derivative
END MODULE fft_1D
