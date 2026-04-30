MODULE fft_1D
  REAL(KIND=8), PRIVATE, PARAMETER :: pi=ACOS(-1.d0)
CONTAINS

!=== FFT and inverse FFT (with/without padding)

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
    !=== WE CHOOSE TO NORMALIZE FORWARD (i.e AVG(un) = un(mF=0))
    scaling = 1.d0/N
    CALL dfftw_plan_dft_r2c_1d(plan,N,r_in,cout,FFTW_ESTIMATE)
    CALL dfftw_execute_dft_r2c(plan, r_in, cout)
    CALL dfftw_destroy_plan(plan)
    DO i = 1, (N+1)/2
       cs_out(i,1) = REAL (cout(i),KIND=8)
       cs_out(i,2) = -AIMAG (cout(i))
    END DO
    cs_out = scaling*cs_out 
    ! cs_out (2:,:) =  cs_out(2:,:)*2 !===modes > 0 ==> factor 2 from exp to cs represt
    cs_out(:,:) = cs_out(:,:) * 2 !=== ==> in case there is a single Fourier mode
    cs_out(1,:) = cs_out(1,:) / 2 
  END SUBROUTINE real_to_fourier

  SUBROUTINE fourier_to_real(cs_in,r_out)
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    REAL(KIND=8),    DIMENSION(:,:)            :: cs_in !===(N,2)
    COMPLEX(KIND=8), DIMENSION(SIZE(cs_in,1))  :: cin !===Complex size N
    REAL(KIND=8), DIMENSION(2*SIZE(cs_in,1)-1) :: r_out !===size 2N-1
    INTEGER(kind=8) :: plan
    INTEGER :: N
    N = 2*SIZE(cs_in,1)-1
    cin = 0.5d0*CMPLX(cs_in(:,1),-cs_in(:,2),KIND=8)
    cin(1) = 2*cin(1)
    CALL dfftw_plan_dft_c2r_1d(plan,N,cin,r_out,FFTW_ESTIMATE)
    CALL dfftw_execute_dft_c2r(plan, cin, r_out)
    CALL dfftw_destroy_plan(plan)
  END SUBROUTINE fourier_to_real

  SUBROUTINE fourier_to_real_padded(cs_in,r_out_pad,Nmax_pad)
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    REAL(KIND=8),    DIMENSION(:,:)            :: cs_in !===(N,2)
    REAL(KIND=8), DIMENSION(Nmax_pad, 2)       :: cs_in_pad !===(N,2)
    ! COMPLEX(KIND=8), DIMENSION(SIZE(cs_in,1))  :: cin !===Complex size N
    INTEGER, INTENT(IN) :: Nmax_pad
    REAL(KIND=8), DIMENSION(2*Nmax_pad-1) :: r_out_pad !===size 2N-1
    INTEGER :: Nmax, N
    Nmax = SIZE(cs_in,1)
    N = 2*SIZE(cs_in,1)-1
    IF (Nmax_pad < Nmax) THEN
       WRITE(*,*) 'BUG in fft.F90: Nmax_pad too small in fourier_to_real_padded=>',&
        Nmax_pad, '<',Nmax
       STOP
    END IF
    !scaling = 2.d0/N
    cs_in_pad = 0.d0
    cs_in_pad(1:Nmax,:) = cs_in(:,:)
    CALL fourier_to_real(cs_in_pad,r_out_pad)
  END SUBROUTINE fourier_to_real_padded

!========= Derivatives using Fourier space

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

!========= Product using padding (or not)
  SUBROUTINE fourier_product(cs_1, cs_2, cs_prod, opt_Nmax)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),              INTENT(IN)   :: cs_1, cs_2
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE               :: r_1, r_2, r_prod
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)  :: cs_prod
    REAL(KIND=8), DIMENSION(SIZE(cs_1,1),SIZE(cs_1,2)) :: cs_out
    INTEGER, INTENT(IN), OPTIONAL :: opt_Nmax
    INTEGER :: raw_Nmax, Nmax, Nmax_real

    raw_Nmax = SIZE(cs_1,1)
    IF (SIZE(cs_1, 1) /= SIZE(cs_2,1)) THEN
       WRITE(*,*) 'BUG in fft.F90: Incompatible sizes in fourier_product=>',&
        SIZE(cs_1, 1), SIZE(cs_2, 1)
       STOP
    END IF
    IF (PRESENT(opt_Nmax)) THEN
       Nmax = opt_Nmax
       IF (Nmax > raw_Nmax) THEN
          WRITE(*,*) 'BUG in fft.F90: opt_Nmax too large in fourier_product=>',&
           opt_Nmax, '>',raw_Nmax
          STOP
       END IF
    ELSE
        Nmax = raw_Nmax
    END IF
    Nmax_real = 2*Nmax - 1

    ALLOCATE(r_1(Nmax_real), r_2(Nmax_real), r_prod(Nmax_real))
    CALL fourier_to_real_padded(cs_1, r_1, Nmax)
    CALL fourier_to_real_padded(cs_2, r_2, Nmax)
    r_prod(:) = r_1(:)*r_2(:)
    DEALLOCATE(r_1, r_2)
    CALL real_to_fourier(r_prod, cs_prod)
    cs_out(:,:) = cs_prod(1:raw_Nmax,:)
    DEALLOCATE(r_prod, cs_prod)

  END SUBROUTINE fourier_product


END MODULE fft_1D
