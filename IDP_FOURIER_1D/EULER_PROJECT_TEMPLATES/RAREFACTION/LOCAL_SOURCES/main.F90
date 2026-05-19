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
  CALL fourier_param%plot_1d(rho_R(fourier_param, euler%time), 'rho_init.plt')
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
    REAL(KIND = 8), DIMENSION(fourier_param%Nmax_real,3) :: r_un
    REAL(KIND = 8), DIMENSION(fourier_param%Nmax_real) :: e, p, rho_e
    REAL(KIND=8) :: error, norm
    INTEGER :: i, k

    DO k = 1, 3
       CALL fourier_to_real(euler%un(:,:,k),r_un(:,k))
    END DO
    e = r_un(:,3)/r_un(:,1) - 0.5d0*(r_un(:,2)/r_un(:,1))**2
    p = euler%eos%pressure(r_un(:,1),e)
    CALL fourier_param%plot_1d(r_un(:,1), 'rhoN.plt')
    CALL fourier_param%plot_1d(p, 'pN.plt')

    rho_e = rho_r(fourier_param, euler%time)
    DO i = 1, fourier_param%Nmax_real
       IF (fourier_param%rr(i)>1.4) THEN
          rho_e(i) = rhor
          r_un(i,1)=rhor
          p(i) = pr
       END IF
    END DO
    CALL fourier_param%plot_1d(r_un(:,1), 'rho.plt')
    CALL fourier_param%plot_1d(p, 'p.plt')

    CALL fourier_param%plot_1d(rho_e, 'rho_e.plt')
    error = SUM(ABS(r_un(:,1)-rho_e))
    norm = SUM(ABS(rho_e))

    WRITE(*,*) ' Error relative L1-norm', error/norm
  END SUBROUTINE errors
END PROGRAM linear_transport
