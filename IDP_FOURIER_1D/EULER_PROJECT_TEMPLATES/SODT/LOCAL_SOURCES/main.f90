PROGRAM linear_transport
  USE setup_module
  USE start_setup_MODULE
  USE my_util
  USE fft_1D
  IMPLICIT NONE
  INTEGER :: n, it
  REAL(KIND = 8) :: tps
  
  CALL start_setup

  euler%un = init(fourier_param,euler%time) 
  !CALL fourier_param%plot_1d(rho_R(fourier_param, euler%time), 'rho_init.plt')
  !CALL fourier_param%plot_1d(mt_R(fourier_param, euler%time), 'mt_init.plt')
  !CALL fourier_param%plot_1d(E_R(fourier_param, euler%time), 'E_init.plt')
  tps = user_time()
  n = 0
  DO WHILE(euler%time < setup_data%final_time)
     CALL euler%update(fourier_param)
     n = n + 1
     !write(*, *) n, euler%time, euler%dt
  END DO
  tps = user_time() - tps
  WRITE(*,*) ' Time per time step per dof times', tps/(fourier_param%Nmax_real*n*euler%ERK%s), tps, n
  CALL errors(euler)
CONTAINS
  SUBROUTINE errors(euler)
    USE euler_module
    USE setup_module
    USE fft_1D
    IMPLICIT NONE
    type(euler_type) :: euler
    REAL(KIND = 8), DIMENSION(fourier_param%Nmax_real) :: r_un
    REAL(KIND=8) :: error, norm
    INTEGER :: i

    IF (euler%time>0.145) RETURN
    CALL fourier_to_real(euler%un(:,:,1),r_un)
    DO i = 1, fourier_param%Nmax_real
       IF (fourier_param%rr(i)<0.22) r_un(i)=rhol
       IF (fourier_param%rr(i)>0.730) r_un(i)=rhor
    END DO
    CALL fourier_param%plot_1d(r_un, 'rho.plt')
    CALL fourier_param%plot_1d(rho_R(fourier_param, euler%time), 'urhoe.plt')
    error = SUM(ABS(r_un-rho_r(fourier_param, euler%time)))
    norm = SUM(ABS(r_un))
    WRITE(*,*) ' Error relative L1-norm', error/norm 
  END SUBROUTINE errors
END PROGRAM linear_transport
