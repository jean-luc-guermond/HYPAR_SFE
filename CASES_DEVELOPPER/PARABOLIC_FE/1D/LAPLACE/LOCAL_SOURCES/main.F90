PROGRAM test_matrix

  USE petsc
  ! USE my_laplace_module, ONLY: laplace, ex_sol, source
  USE fem_tn, ONLY: ns_l1_par
  USE start_setup_MODULE

  USE sub_plot
  USE post_processing_debug_MODULE
  USE setup, ONLY: source, ex_sol
  USE my_util, ONLY: user_time
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhs, exact_solution, un_out
  REAL(KIND=8) :: tps
  INTEGER      :: tot_np
  CHARACTER(5) :: char

  INTEGER :: ierr


!========================!
!==== INITIALIZATION ====!
!========================!

  CALL start_setup
  ALLOCATE(rhs(Laplace%mesh%np), un_out(Laplace%mesh%np))
  rhs = source(Laplace%mesh%rr)
  WRITE(char, '(I5)') Laplace%mesh%rank
  CALL plot_scalar_field(Laplace%mesh%jj, Laplace%mesh%rr, rhs, 'rhs'//TRIM(ADJUSTL(char))//'.plt')
  IF (setup_data%if_analytical_ref) THEN
    ALLOCATE(exact_solution(Laplace%mesh%np))
    exact_solution = Laplace%dir_bc(Laplace%mesh%rr)
    CALL plot_scalar_field(Laplace%mesh%jj, Laplace%mesh%rr, exact_solution, 'sol_exact_'//TRIM(ADJUSTL(char))//'.plt')
  END IF


!=====================!
!==== SOLVER LOOP ====!
!=====================!

  tps = user_time()
  CALL Laplace%solve(rhs, un_out)
  tps = user_time() - tps


!=========================!
!==== POST-PROCESSING ====!
!=========================!

  CALL MPI_ALLREDUCE(Laplace%mesh%dom_np,tot_np,1,MPI_INTEGER,MPI_SUM,Laplace%communicator,ierr)
  IF(Laplace%mesh%rank==0) THEN
     WRITE(*,*) ' tot_np', tot_np
     WRITE(*,*) ' Time per dof times proc', tps/tot_np, tps
  END IF
  CALL plot_scalar_field(Laplace%mesh%jj, Laplace%mesh%rr, un_out(:), 'result_' // TRIM(ADJUSTL(char)) // '.plt')
  
  IF (setup_data%if_analytical_ref) THEN
    CALL plot_scalar_field(mesh%jj, mesh%rr, un_out-exact_solution, 'error' // trim(adjustl(char)) // '.plt')
  END IF

!=========================!
!==== REGRESSION TEST ====!
!=========================!

  CALL errors

!=====================!
!==== END PROGRAM ====!
!=====================!
  CALL PetscFinalize(ierr)


CONTAINS
  SUBROUTINE errors
    USE fem_tn
    USE post_processing_debug_MODULE
    IMPLICIT NONE

    REAL(KIND=8) :: error, norm, norm_anal
    REAL(KIND=8), DIMENSION(setup_data%syst_size) :: tab_norm
    INTEGER      :: n, num_test

!==== Put final processing stuff here ====!
    DO n=1, setup_data%syst_size
      IF (setup_data%if_analytical_ref) THEN
        CALL ns_l1_PAR(mesh, un_out(:)-exact_solution, error, Laplace%communicator)
        CALL ns_l1_PAR(mesh, exact_solution, norm_anal, Laplace%communicator)
        norm = error/norm_anal
        IF(Laplace%mesh%rank==0) WRITE(*, *) 'Comp = ',n,'; Relative error, L1-norm = ', error/norm_anal
        ! IF(Laplace%mesh%rank==0) WRITE(*, '(A,I0,A,g12.3)') 'Comp = ',n,'; Relative error, L1-norm = ', error/norm_anal
      ELSE
        CALL ns_l1_PAR(mesh, un_out(:), norm, Laplace%communicator)
        IF(Laplace%mesh%rank==0) WRITE(*, *) 'Comp = ',n,'; no analytical ref, L1-norm = ', norm
      END IF
    END DO

    
!==== For regression tests ====!
    IF (setup_data%if_regression_test) THEN
      DO n=1, setup_data%syst_size
        IF (setup_data%if_analytical_ref) THEN
          CALL ns_l1_PAR(mesh, un_out(:)-exact_solution, error, Laplace%communicator)
          CALL ns_l1_PAR(mesh, exact_solution, norm_anal, Laplace%communicator)
          norm = error/norm_anal
        ELSE
          CALL ns_l1_PAR(mesh, un_out(:), norm, Laplace%communicator)
        END IF
        tab_norm(n) = norm
      END DO
      CALL get_num_test(num_test)
      CALL regression(tab_norm, opt_num_test=num_test)
    END IF

  END SUBROUTINE errors

END PROGRAM test_matrix
