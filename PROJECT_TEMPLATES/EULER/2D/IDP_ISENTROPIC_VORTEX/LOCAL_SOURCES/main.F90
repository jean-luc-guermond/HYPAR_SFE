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
  ALLOCATE(un(mesh%np, euler%syst_dim))
  CALL euler%bc%initial_condition(un, 0.d0, euler%mesh%rr)

  WRITE(char, '(I5)') euler%mesh%rank
  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:, 1), 'initrho'//TRIM(ADJUSTL(char))//'.plt')

!=====================!
!==== SOLVER LOOP ====!
!=====================!

  tps = user_time()
  n = 0
  DO WHILE(euler%time < setup_data%final_time)
    CALL euler%update(un)
    n = n + 1
    IF (MOD(n, setup_data%verbose_freq)==0) THEN
        IF (euler%mesh%rank==0) write(*, *) n, euler%time, euler%dt
    END IF
    IF (n == setup_data%max_it) THEN
        IF (euler%mesh%rank==0) WRITE(*,*) "max_it reached, exiting solver loop"
        EXIT
    END IF
  END DO
  tps = user_time() - tps

!=========================!
!==== POST-PROCESSING ====!
!=========================!

  CALL MPI_ALLREDUCE(euler%mesh%dom_np,tot_np,1,MPI_INTEGER,MPI_SUM,euler%communicator,code)
  IF(euler%mesh%rank==0) THEN
     WRITE(*,*) ' tot_np', tot_np
     WRITE(*,*) ' Time per time step per dof times proc', tps/(tot_np*n), tps, n
  END IF
  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:, 1), 'rho' // TRIM(ADJUSTL(char)) // '.plt')

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
    IMPLICIT NONE

    REAL(KIND=8) :: error_loc, norm_loc, norm_anal_loc, error, norm, norm_anal
    REAL(KIND = 8), DIMENSION(size(un, 2)) :: tab_norm
    INTEGER :: n, code

!==== Put final processing stuff here ====!
    DO n=1, SIZE(un,2)
      IF (setup_data%if_analytical_ref) THEN
        CALL ns_l1_PAR(mesh, un(:,n)-euler%bc%sol_anal(n, euler%time,mesh%rr), error, euler%communicator)
        CALL ns_l1_PAR(mesh, euler%bc%sol_anal(n, euler%time,mesh%rr), norm_anal, euler%communicator)
        norm = error/norm_anal
        IF(euler%mesh%rank==0) WRITE(*, *) 'Comp = ',n,'; Relative error, L1-norm = ', error/norm_anal, norm_anal
      ELSE
        CALL ns_l1_PAR(mesh, un(:,n), norm, euler%communicator)
        IF(euler%mesh%rank==0) WRITE(*, *) 'Comp = ',n,'; no analytical ref, L1-norm = ', norm
      END IF
    END DO

    
!==== For regression tests ====!
    IF (setup_data%if_regression_test) THEN
      DO n=1, SIZE(un,2)
        IF (setup_data%if_analytical_ref) THEN
          CALL ns_l1_PAR(mesh, un(:,n)-euler%bc%sol_anal(n, euler%time,mesh%rr), error, euler%communicator)
          CALL ns_l1_PAR(mesh, euler%bc%sol_anal(n, euler%time,mesh%rr), norm_anal, euler%communicator)
          norm = error/norm_anal
        ELSE
          CALL ns_l1_PAR(mesh, un(:,n), norm, euler%communicator)
        END IF
        tab_norm(n) = norm
      END DO
      CALL get_num_test(num_test)
      CALL regression(tab_norm, opt_num_test=num_test)
    END IF

  END SUBROUTINE errors
END PROGRAM prog
