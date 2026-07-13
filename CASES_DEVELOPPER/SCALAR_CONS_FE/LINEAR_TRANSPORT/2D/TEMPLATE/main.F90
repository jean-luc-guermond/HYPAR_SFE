PROGRAM prog
#include "petsc/finclude/petsc.h"
  USE petsc
  USE start_setup_MODULE
  USE setup
  USE sub_plot
  USE my_util
  IMPLICIT NONE
  REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: un
  REAL(KIND = 8) :: tps
  CHARACTER(5) :: char
  INTEGER :: n, tot_np, code, num_test

!========================!
!==== INITIALIZATION ====!
!========================!

  CALL start_setup
  ALLOCATE(un(mesh%np, linear_transport%syst_dim))
  CALL linear_transport%bc%initial_condition(un, 0.d0, linear_transport%mesh%rr)

  WRITE(char, '(I5)') linear_transport%mesh%rank
  CALL plot_scalar_field(linear_transport%mesh%jj, linear_transport%mesh%rr, un(:, 1), 'initrho'//TRIM(ADJUSTL(char))//'.plt')

!=====================!
!==== SOLVER LOOP ====!
!=====================!

  tps = user_time()
  n = 0
  DO WHILE(linear_transport%time < setup_data%final_time)
    CALL linear_transport%update(un)
    n = n + 1
    IF (MOD(n, setup_data%verbose_freq)==0) THEN
        IF (linear_transport%mesh%rank==0) write(*, *) n, linear_transport%time, linear_transport%dt
    END IF
    IF (n == setup_data%max_it) THEN
        IF (linear_transport%mesh%rank==0) WRITE(*,*) "max_it reached, exiting solver loop"
        EXIT
    END IF
  END DO
  tps = user_time() - tps

!=========================!
!==== POST-PROCESSING ====!
!=========================!

  CALL MPI_ALLREDUCE(linear_transport%mesh%dom_np,tot_np,1,MPI_INTEGER,MPI_SUM,linear_transport%communicator,code)
  IF(linear_transport%mesh%rank==0) THEN
     WRITE(*,*) ' tot_np', tot_np
     WRITE(*,*) ' Time per time step per dof times proc', tps/(tot_np*n), tps, n
  END IF
  CALL plot_scalar_field(linear_transport%mesh%jj, linear_transport%mesh%rr, un(:, 1), 'rho' // TRIM(ADJUSTL(char)) // '.plt')
  CALL plot_scalar_field(linear_transport%mesh%jj, linear_transport%mesh%rr, linear_transport%bc%sol_anal(1, linear_transport%time,mesh%rr), 'sol_exact' // TRIM(ADJUSTL(char)) // '.plt')

!=========================!
!==== REGRESSION TEST ====!
!=========================!

    CALL errors

!=====================!
!==== END PROGRAM ====!
!=====================!
    CALL PetscFinalize(code)

CONTAINS
  SUBROUTINE errors
    USE fem_tn
    USE post_processing_debug_MODULE
    USE options_module
    IMPLICIT NONE

    REAL(KIND=8) :: error, norm, norm_anal
    REAL(KIND = 8), DIMENSION(size(un, 2)) :: tab_norm
    INTEGER :: n

!==== Put final processing stuff here ====!
    DO n=1, SIZE(un,2)
      IF (setup_data%if_analytical_ref) THEN
        CALL ns_l1_PAR(mesh, un(:,n)-linear_transport%bc%sol_anal(n, linear_transport%time,mesh%rr), error, linear_transport%communicator)
        CALL ns_l1_PAR(mesh, linear_transport%bc%sol_anal(n, linear_transport%time,mesh%rr), norm_anal, linear_transport%communicator)
        norm = error/norm_anal
        IF(linear_transport%mesh%rank==0) WRITE(*, *) 'Comp ', linear_transport%name_comp(n), '; Relative error, L1-norm = ', error/norm_anal
        ! IF(linear_transport%mesh%rank==0) WRITE(*, '(A,I0,A,g12.3)') 'Comp = ',n,'; Relative error, L1-norm = ', error/norm_anal
      ELSE
        CALL ns_l1_PAR(mesh, un(:,n), norm, linear_transport%communicator)
        IF(linear_transport%mesh%rank==0) WRITE(*, *) 'Comp ', linear_transport%name_comp(n), '; no analytical ref, L1-norm = ', norm
      END IF
    END DO

    
!==== For regression tests ====!
    IF (options%if_regression) THEN
      DO n=1, SIZE(un,2)
        IF (setup_data%if_analytical_ref) THEN
          CALL ns_l1_PAR(mesh, un(:,n)-linear_transport%bc%sol_anal(n, linear_transport%time,mesh%rr), error, linear_transport%communicator)
          CALL ns_l1_PAR(mesh, linear_transport%bc%sol_anal(n, linear_transport%time,mesh%rr), norm_anal, linear_transport%communicator)
          IF (norm_anal<1d-13) THEN
            norm = error
          ELSE
            norm = error/norm_anal
          END IF
        ELSE
          CALL ns_l1_PAR(mesh, un(:,n), norm, linear_transport%communicator)
        END IF
        tab_norm(n) = norm
      END DO
      write(*,*) 'tab_norm = ', tab_norm
      CALL regression(tab_norm, opt_num_test=options%num_regex)
    END IF

  END SUBROUTINE errors
END PROGRAM prog
