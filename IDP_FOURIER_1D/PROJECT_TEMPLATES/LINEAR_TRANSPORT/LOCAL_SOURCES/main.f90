PROGRAM linear_transport
  USE setup_module
  USE start_setup_MODULE
  USE my_util
  USE fft_1D
  IMPLICIT NONE
  INTEGER :: n
  REAL(KIND = 8) :: tps
  REAL(KIND = 8), DIMENSION(:), POINTER :: r_un
  CALL start_setup

  nl_scalar_cons%un = exact_sol_F(fourier_param, nl_scalar_cons%time,exact_sol_step_R) 
  CALL fourier_param%plot_1d(exact_sol_step_R(fourier_param, nl_scalar_cons%time), 'u_init.plt')

  tps = user_time()
  n = 0
  DO WHILE(nl_scalar_cons%time < setup_data%final_time)
     CALL nl_scalar_cons%update(fourier_param)
     n = n + 1
     write(*, *) n, nl_scalar_cons%time, nl_scalar_cons%dt
  END DO
  tps = user_time() - tps

  ALLOCATE(r_un(fourier_param%Nmax_real))
  CALL fourier_to_real(nl_scalar_cons%un,r_un)
  CALL fourier_param%plot_1d(r_un, 'u.plt')
END PROGRAM linear_transport
