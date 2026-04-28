PROGRAM nl_scalar_cons_fft_1d
  USE setup_module
  USE start_setup_MODULE
  USE my_util
  USE fft_1D
  IMPLICIT NONE
  INTEGER :: n, num_test
  REAL(KIND = 8) :: tps

!========================!
!==== INITIALIZATION ====!
!========================!

  CALL start_setup

  nl_scalar_cons%un = exact_sol_F(fourier_param, nl_scalar_cons%time,exact_sol_R)
  CALL fourier_param%plot_1d(exact_sol_R(fourier_param, nl_scalar_cons%time), 'u_init.plt')

!=====================!
!==== SOLVER LOOP ====!
!=====================!

  tps = user_time()
  n = 0
  DO WHILE(nl_scalar_cons%time < setup_data%final_time)
     CALL nl_scalar_cons%update
     n = n + 1
     !write(*, *) n, nl_scalar_cons%time, nl_scalar_cons%dt
  END DO
  tps = user_time() - tps

!=========================!
!==== POST-PROCESSING ====!
!=========================!

  WRITE(*,*) ' Time per time step per dof times', &
       tps/(fourier_param%Nmax_real*n*nl_scalar_cons%ERK%s), tps, n

!=========================!
!==== REGRESSION TEST ====!
!=========================!

    CALL errors(nl_scalar_cons)

!=====================!
!==== END PROGRAM ====!
!=====================!

CONTAINS
  SUBROUTINE errors(nl_scalar_cons)
    USE nl_scalar_cons_module
    USE fft_1D
    USE post_processing_debug_MODULE
    IMPLICIT NONE
    type(nl_scalar_cons_type) :: nl_scalar_cons
    REAL(KIND = 8), DIMENSION(fourier_param%Nmax_real) :: r_un
    REAL(KIND = 8) :: error, norm
    REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: tab_norm

    CALL fourier_to_real(nl_scalar_cons%un,r_un)
    CALL fourier_param%plot_1d(r_un, 'u_fin.plt')
    CALL fourier_param%plot_1d(exact_sol_R(fourier_param, nl_scalar_cons%time), 'u_th.plt')
    error = SUM(ABS(r_un-exact_sol_R(fourier_param, nl_scalar_cons%time)))
    norm = SUM(ABS(r_un))
    WRITE(*,*) ' Error relative L1-norm', error/norm

    ! IF (setup_data%if_regression_test) THEN
    !    ALLOCATE(tab_norm(1))
    !    tab_norm(1) = norm
    !    CALL get_num_test(num_test)
    !    CALL regression(tab_norm, opt_num_test=num_test)
    ! END IF
    IF (setup_data%if_regression_test) THEN
      ALLOCATE(tab_norm(1))
      IF (setup_data%if_analytical_ref) THEN
        error = SUM(ABS(r_un-exact_sol_R(fourier_param, nl_scalar_cons%time)))
        norm = SUM(ABS(r_un))
        tab_norm(1) = error/norm
      ELSE
        norm = SUM(ABS(r_un))
        tab_norm(1) = norm
      END IF
      CALL get_num_test(num_test)
      CALL regression(tab_norm, opt_num_test=num_test)
    END IF

  END SUBROUTINE errors


  ! SUBROUTINE errors
  !   USE fem_tn
  !   USE post_processing_debug_MODULE
  !   IMPLICIT NONE

  !   REAL(KIND=8) :: error_loc, norm_loc, norm_anal_loc, error, norm, norm_anal
  !   REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: tab_norm
  !   INTEGER :: n, code

  !   IF (setup_data%if_analytical_ref) THEN
  !      CALL ns_l1(mesh, un(:,1)-rho_anal(euler%time,mesh%rr), error_loc)
  !      CALL MPI_ALLREDUCE(error_loc,error,1,MPI_DOUBLE_PRECISION,MPI_SUM,euler%communicator,code)
  !      CALL ns_l1(mesh, rho_anal(euler%time,mesh%rr), norm_anal_loc)
  !      CALL MPI_ALLREDUCE(norm_anal_loc,norm_anal,1,MPI_DOUBLE_PRECISION,MPI_SUM,euler%communicator,code)
  !      IF(euler%mesh%rank==0) WRITE(*, '(A,g12.3)') 'Error density relative, L1-norm ', error/norm_anal
  !   END IF

  !   IF (setup_data%if_regression_test) THEN
  !     ALLOCATE(tab_norm(size(un, 2)))
  !     DO n=1, SIZE(un,2)
  !       IF (setup_data%if_analytical_ref) THEN
  !         CALL ns_l1(mesh, un(:,n)-sol_anal(n, euler%time,mesh%rr), error_loc)
  !         CALL ns_l1(mesh, sol_anal(n, euler%time,mesh%rr), norm_loc)
  !         norm_loc = error_loc/norm_loc
  !       ELSE
  !         CALL ns_l1(mesh, un(:,n), norm_loc)
  !       END IF
  !       CALL MPI_ALLREDUCE(norm_loc,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,euler%communicator,code)
  !       tab_norm(n) = norm
  !     END DO
  !     CALL get_num_test(num_test)
  !     CALL regression(tab_norm, opt_num_test=num_test)
  !   END IF
  ! END SUBROUTINE errors

END PROGRAM nl_scalar_cons_fft_1d
