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
  ALLOCATE(un(mesh%np, burgers%syst_dim))
  CALL burgers%bc%initial_condition(un, 0.d0, burgers%mesh%rr)

  WRITE(char, '(I5)') burgers%mesh%rank
  CALL plot_scalar_field(burgers%mesh%jj, burgers%mesh%rr, un(:, 1), 'initrho'//TRIM(ADJUSTL(char))//'.plt')

!=====================!
!==== SOLVER LOOP ====!
!=====================!

  tps = user_time()
  n = 0
  DO WHILE(burgers%time < setup_data%final_time)
    CALL burgers%update(un)
    n = n + 1
    IF (MOD(n, setup_data%verbose_freq)==0) THEN
        IF (burgers%mesh%rank==0) write(*, *) n, burgers%time, burgers%dt
    END IF
    IF (n == setup_data%max_it) THEN
        IF (burgers%mesh%rank==0) WRITE(*,*) "max_it reached, exiting solver loop"
        EXIT
    END IF
  END DO
  tps = user_time() - tps

!=========================!
!==== POST-PROCESSING ====!
!=========================!

  CALL MPI_ALLREDUCE(burgers%mesh%dom_np,tot_np,1,MPI_INTEGER,MPI_SUM,burgers%communicator,code)
  IF(burgers%mesh%rank==0) THEN
     WRITE(*,*) ' tot_np', tot_np
     WRITE(*,*) ' Time per time step per dof times proc', tps/(tot_np*n), tps, n
  END IF
  CALL plot_scalar_field(burgers%mesh%jj, burgers%mesh%rr, un(:, 1), 'rho' // TRIM(ADJUSTL(char)) // '.plt')
  CALL plot_scalar_field(burgers%mesh%jj, burgers%mesh%rr, burgers%bc%sol_anal(1, burgers%time,mesh%rr), 'sol_exact' // TRIM(ADJUSTL(char)) // '.plt')

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

    REAL(KIND=8) :: error, norm, norm_anal
    REAL(KIND = 8), DIMENSION(size(un, 2)) :: tab_norm
    INTEGER :: n

!==== Put final processing stuff here ====!
    DO n=1, SIZE(un,2)
      IF (setup_data%if_analytical_ref) THEN
        CALL ns_l1_PAR(mesh, un(:,n)-burgers%bc%sol_anal(n, burgers%time,mesh%rr), error, burgers%communicator)
        CALL ns_l1_PAR(mesh, burgers%bc%sol_anal(n, burgers%time,mesh%rr), norm_anal, burgers%communicator)
        norm = error/norm_anal
        IF(burgers%mesh%rank==0) WRITE(*, *) 'Comp ', burgers%name_comp(n), '; Relative error, L1-norm = ', error/norm_anal
        ! IF(burgers%mesh%rank==0) WRITE(*, '(A,I0,A,g12.3)') 'Comp = ',n,'; Relative error, L1-norm = ', error/norm_anal
      ELSE
        CALL ns_l1_PAR(mesh, un(:,n), norm, burgers%communicator)
        IF(burgers%mesh%rank==0) WRITE(*, *) 'Comp ', burgers%name_comp(n), '; no analytical ref, L1-norm = ', norm
      END IF
    END DO

    
!==== For regression tests ====!
    IF (setup_data%if_regression_test) THEN
      DO n=1, SIZE(un,2)
        IF (setup_data%if_analytical_ref) THEN
          CALL ns_l1_PAR(mesh, un(:,n)-burgers%bc%sol_anal(n, burgers%time,mesh%rr), error, burgers%communicator)
          CALL ns_l1_PAR(mesh, burgers%bc%sol_anal(n, burgers%time,mesh%rr), norm_anal, burgers%communicator)
          IF (norm_anal<1d-13) THEN
            norm = error
          ELSE
            norm = error/norm_anal
          END IF
        ELSE
          CALL ns_l1_PAR(mesh, un(:,n), norm, burgers%communicator)
        END IF
        tab_norm(n) = norm
      END DO
      write(*,*) 'tab_norm = ', tab_norm
      CALL get_num_test(num_test)
      CALL regression(tab_norm, opt_num_test=num_test)
    END IF

  END SUBROUTINE errors
END PROGRAM prog
