MODULE fft_1D_plan

  IMPLICIT NONE
  PUBLIC    :: real_to_fourier_plan, fourier_to_real_plan, fourier_derivative_plan, real_derivative_plan

  REAL(KIND=8), PRIVATE, PARAMETER :: pi=ACOS(-1.d0)
CONTAINS

  SUBROUTINE real_to_fourier_plan(r_in,cs_out)
    !>  This subroutine performs FFT of a 2D array r_in(mesh_np, Theta)
    !!  along its second dimension

    USE petsc
#include "petsc/finclude/petsc.h"
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: r_in !===size 2N-1
    COMPLEX(KIND=8), DIMENSION(SIZE(r_in,1),SIZE(r_in,2)/2+1)   :: c_out !===Complex size N
    REAL(KIND=8),    DIMENSION(SIZE(r_in,1),SIZE(r_in,2)/2+1,2), INTENT(OUT) :: cs_out
    REAL(KIND=8)          :: scaling
    INTEGER(KIND=8)       :: plan
    INTEGER, DIMENSION(1) :: N_fourier, N_theta, inemebed
    INTEGER :: i, howmany_simult_FFT, fft_dim, mesh_np
    INTEGER :: inembed,istride,idist,onembed,ostride,odist

    fft_dim = 1                        ! 1d_fft
    N_fourier = (SIZE(r_in, 2)+1)/2    ! nb fourier modes
    N_theta = 2*N_fourier-1
    howmany_simult_FFT = SIZE(r_in, 1) ! total number of FFTs performed

    mesh_np = SIZE(r_in, 1)

    !=========================================================================================
    !===== FFTW PARAMETERS IF INPUT ARRAY OF TYPE arr(mesh%np, MF) ===========================
    !=========================================================================================
    !>   exchange the values istride <--> idist as well as ostride <--> odist
    !!   if the input array is arr(MF, mesh%np)
    !! <<<< WARNING >>>> be aware that in Fortran, moving of 1 element in array amounts to 
    !!                   an increment of 1 in the last dimension (i.e move across columns, not lines)

    inembed = mesh_np ! parameter > N_theta iff FFT performed on a subset of input array
    istride = mesh_np ! distance between two consecutive arrays to be FFT-ed
    idist = 1 ! distance between two consecutive Fourier modes inside a contiguous array

    onembed = mesh_np ! parameter > N_fourier iff FFT performed on a subset of input array
    ostride = mesh_np ! distance between two consecutive arrays in output FFT-ed
    odist = 1 ! distance between two consecutive angles inside a contiguous array in output

    scaling = 1.d0/N_theta(1)

    CALL dfftw_plan_many_dft_r2c(plan, fft_dim, N_theta, howmany_simult_FFT, &
         r_in,  inembed, istride, idist, &
         c_out, onembed, ostride, odist, &
         FFTW_ESTIMATE)    

    CALL dfftw_execute_dft_r2c(plan, r_in, c_out)
    CALL dfftw_destroy_plan(plan)

    DO i = 1, (N_theta(1)+1)/2
       cs_out(:,i,1) = REAL (c_out(:, i),KIND=8)
       cs_out(:,i,2) = -AIMAG (c_out(:, i))
    END DO
    cs_out = scaling*cs_out 
    ! cs_out (2:,:) =  cs_out(2:,:)*2 !===modes > 0 ==> factor 2 from exp to cs represt
    cs_out(:,:,:) = cs_out(:,:,:) * 2 !=== ==> in case there is a single Fourier mode
    cs_out(:,1,:) = cs_out(:,1,:) / 2 
  END SUBROUTINE real_to_fourier_plan


  SUBROUTINE fourier_to_real_plan(cs_in,r_out)
    !>  This subroutine performs IFFT of a 3D array r_in(mesh_np, MF, cs)
    !!  along its second/third dimension

    USE petsc
#include "petsc/finclude/petsc.h"
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    REAL(KIND=8),    DIMENSION(:,:,:), INTENT(IN)            :: cs_in !===(N,2)
    COMPLEX(KIND=8), DIMENSION(SIZE(cs_in,1),SIZE(cs_in,2))  :: c_in !===Complex size N
    REAL(KIND=8), DIMENSION(SIZE(cs_in,1),2*SIZE(cs_in,2)-1), INTENT(OUT) :: r_out !===size 2N-1
    INTEGER(KIND=8)       :: plan
    INTEGER, DIMENSION(1) :: N_fourier, N_theta, inemebed
    INTEGER :: i, howmany_simult_FFT, fft_dim, mesh_np
    INTEGER :: inembed,istride,idist,onembed,ostride,odist

    N_theta = 2*SIZE(cs_in,2)-1
    c_in(:,:) = 0.5d0*CMPLX(cs_in(:,:,1),-cs_in(:,:,2),KIND=8)
    c_in(:,1) = 2*c_in(:,1) !=== MODE ZERO different treatment

    fft_dim = 1                        ! 1d_fft
    N_fourier = SIZE(cs_in,2)          ! nb fourier modes
    howmany_simult_FFT = SIZE(cs_in, 1) ! total number of FFTs performed

    mesh_np = SIZE(cs_in, 1)

    !=========================================================================================
    !===== FFTW PARAMETERS IF INPUT ARRAY OF TYPE arr(mesh%np, MF) ===========================
    !=========================================================================================
    !>   exchange the values istride <--> idist as well as ostride <--> odist
    !!   if the input array is arr(MF, mesh%np)
    !! <<<< WARNING >>>> be aware that in Fortran, moving of 1 element in array amounts to 
    !!                   an increment of 1 in the last dimension (i.e move across columns, not lines)

    inembed = mesh_np ! parameter > N_theta iff FFT performed on a subset of input array
    istride = mesh_np ! distance between two consecutive arrays to be FFT-ed
    idist = 1 ! distance between two consecutive Fourier modes inside a contiguous array

    onembed = mesh_np ! parameter > N_fourier iff FFT performed on a subset of input array
    ostride = mesh_np ! distance between two consecutive arrays in output FFT-ed
    odist = 1 ! distance between two consecutive angles inside a contiguous array in output

    CALL dfftw_plan_many_dft_c2r(plan, fft_dim, N_theta, howmany_simult_FFT, &
         c_in, onembed, ostride, odist, &
         r_out,  inembed, istride, idist, &
         FFTW_ESTIMATE)

    CALL dfftw_execute_dft_c2r(plan, c_in, r_out)
    CALL dfftw_destroy_plan(plan)

  END SUBROUTINE fourier_to_real_plan

!========= Derivatives using Fourier space

  SUBROUTINE fourier_derivative_plan(cs_in,cs_out,Length)
    IMPLICIT NONE
    REAL(KIND=8),    DIMENSION(:,:,:), INTENT(IN) :: cs_in 
    REAL(KIND=8), DIMENSION(SIZE(cs_in,1), SIZE(cs_in,2), SIZE(cs_in,3)), INTENT(OUT) :: cs_out
    REAL(KIND=8) :: Length
    INTEGER :: i
    DO i = 1, SIZE(cs_in,2)
       cs_out(:,i,1) =  (i-1)*cs_in(:,i,2)
       cs_out(:,i,2) = -(i-1)*cs_in(:,i,1)
    END DO
    cs_out = 2*pi*cs_out/Length

  END SUBROUTINE fourier_derivative_plan

  SUBROUTINE real_derivative_plan(r_in,r_out,Length)
    IMPLICIT NONE
    INCLUDE 'fftw3.f'
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: r_in !===size 2N-1
    REAL(KIND=8), DIMENSION(SIZE(r_in,1), SIZE(r_in,2)), INTENT(OUT) :: r_out
    REAL(KIND=8), DIMENSION(SIZE(r_in,1), (SIZE(r_in,2)+1)/2,2) :: cs_in, cs_out
    REAL(KIND=8), INTENT(IN) :: Length
    CALL real_to_fourier_plan(r_in,cs_in)
    CALL fourier_derivative_plan(cs_in,cs_out,Length)
    CALL fourier_to_real_plan(cs_out,r_out)
  END SUBROUTINE real_derivative_plan



END MODULE fft_1D_plan
