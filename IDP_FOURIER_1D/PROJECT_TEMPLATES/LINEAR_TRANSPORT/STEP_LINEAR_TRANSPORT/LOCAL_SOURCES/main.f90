PROGRAM linear_transport
  USE setup_module
  USE start_setup_MODULE
  USE my_util
  USE fft_1D
  IMPLICIT NONE
  INTEGER :: n, it
  REAL(KIND = 8) :: tps
  
  CALL start_setup

  nl_scalar_cons%un = exact_sol_F(fourier_param, nl_scalar_cons%time,exact_sol_step_R) 
  CALL fourier_param%plot_1d(exact_sol_step_R(fourier_param, nl_scalar_cons%time), 'u_init.plt')

  tps = user_time()
  n = 0
  DO WHILE(nl_scalar_cons%time < setup_data%final_time)
     CALL nl_scalar_cons%update(fourier_param)
     n = n + 1
     !write(*, *) n, nl_scalar_cons%time, nl_scalar_cons%dt
  END DO
  tps = user_time() - tps
  WRITE(*,*) ' Time per time step per dof times', tps/(fourier_param%Nmax_real*n*nl_scalar_cons%ERK%s), tps, n
  CALL errors(nl_scalar_cons)
CONTAINS
  SUBROUTINE errors(nl_scalar_cons)
    USE nl_scalar_cons_module
    USE fft_1D
    IMPLICIT NONE
    type(nl_scalar_cons_type) :: nl_scalar_cons
    REAL(KIND = 8), DIMENSION(fourier_param%Nmax_real) :: r_un
    REAL(KIND=8) :: error, norm
    
    CALL fourier_to_real(nl_scalar_cons%un,r_un)
    CALL fourier_param%plot_1d(r_un, 'u.plt')
    CALL fourier_param%plot_1d(exact_sol_step_R(fourier_param, nl_scalar_cons%time), 'ue.plt')
    error = SUM(ABS(r_un-exact_sol_step_R(fourier_param, nl_scalar_cons%time)))
    norm = SUM(ABS(r_un))
    WRITE(*,*) ' Error relative L1-norm', error/norm 
  END SUBROUTINE errors
END PROGRAM linear_transport
