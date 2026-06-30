PROGRAM prog
#include "petsc/finclude/petsc.h"
  USE petsc
  USE start_setup_MODULE
  USE setup
  USE sub_plot
  USE ETA_module
  USE euler_post_proc_module
  USE my_util
  IMPLICIT NONE
  REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: un
  REAL(KIND = 8) :: tps
  CHARACTER(5) :: char
  INTEGER :: n, tot_np, code, num_test
  TYPE(ETA_type) :: ETA
!========================!
!==== INITIALIZATION ====!
!========================!

  CALL start_setup
  ALLOCATE(un(mesh%np, navier_stokes%euler%syst_dim))
  CALL navier_stokes%euler%bc%initial_condition(un, 0.d0, navier_stokes%euler%mesh%rr)

  WRITE(char, '(I5)') navier_stokes%euler%mesh%rank
  CALL plot_scalar_field(navier_stokes%euler%mesh%jj, navier_stokes%euler%mesh%rr, un(:, 1), 'initrho'//TRIM(ADJUSTL(char))//'.plt')

!=====================!
!==== SOLVER LOOP ====!
!=====================!
  CALL ETA%init(navier_stokes%mesh%rank, navier_stokes%time)
  tps = user_time()
  n = 0
  DO WHILE(navier_stokes%time < setup_data%final_time)
    CALL navier_stokes%update(un)
    CALL ETA%update(navier_stokes%dt)
    n = n + 1
    IF (MOD(n, setup_data%verbose_freq)==0) THEN
        CALL ETA%print(navier_stokes%time, setup_data%final_time)
        !IF (navier_stokes%euler%mesh%rank==0) write(*, *) n, navier_stokes%time, navier_stokes%dt
    END IF
    IF (n == setup_data%max_it) THEN
        IF (navier_stokes%euler%mesh%rank==0) WRITE(*,*) "max_it reached, exiting solver loop"
        EXIT
    END IF
  END DO
  tps = user_time() - tps

!=========================!
!==== POST-PROCESSING ====!
!=========================!

  CALL MPI_ALLREDUCE(navier_stokes%euler%mesh%dom_np,tot_np,1,MPI_INTEGER,MPI_SUM,navier_stokes%euler%communicator,code)
  IF(navier_stokes%euler%mesh%rank==0) THEN
     WRITE(*,*) ' tot_np', tot_np
     WRITE(*,*) ' Time per time step per dof times proc', navier_stokes%mesh%nb_proc*tps/(tot_np*n), tps, n
  END IF
  CALL plot_scalar_field(navier_stokes%euler%mesh%jj, navier_stokes%euler%mesh%rr, un(:, 1), 'rho' // TRIM(ADJUSTL(char)) // '.plt')
  CALL plot_scalar_field(navier_stokes%euler%mesh%jj, navier_stokes%euler%mesh%rr, &
  navier_stokes%euler%bc%sol_anal(1, navier_stokes%euler%time,mesh%rr), 'rhoexact' // TRIM(ADJUSTL(char)) // '.plt')
write(*,*) 'time', navier_stokes%time, navier_stokes%stokes%time, navier_stokes%euler%time
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
        CALL ns_l1_PAR(mesh, un(:,n)-navier_stokes%euler%bc%sol_anal(n, navier_stokes%euler%time,mesh%rr), error, navier_stokes%euler%communicator)
        CALL ns_l1_PAR(mesh, navier_stokes%euler%bc%sol_anal(n, navier_stokes%euler%time,mesh%rr), norm_anal, navier_stokes%euler%communicator)
        norm = error/norm_anal
        IF(navier_stokes%euler%mesh%rank==0) WRITE(*, *) 'Comp ', navier_stokes%euler%name_comp(n), '; Relative error, L1-norm = ', error/norm_anal
        ! IF(euler%mesh%rank==0) WRITE(*, '(A,I0,A,g12.3)') 'Comp = ',n,'; Relative error, L1-norm = ', error/norm_anal
      ELSE
        CALL ns_l1_PAR(mesh, un(:,n), norm, navier_stokes%euler%communicator)
        IF(navier_stokes%euler%mesh%rank==0) WRITE(*, *) 'Comp ', navier_stokes%euler%name_comp(n), '; no analytical ref, L1-norm = ', norm
      END IF
    END DO

    
!==== For regression tests ====!
    IF (setup_data%if_regression_test) THEN
      DO n=1, SIZE(un,2)
        IF (setup_data%if_analytical_ref) THEN
          CALL ns_l1_PAR(mesh, un(:,n)-navier_stokes%euler%bc%sol_anal(n, navier_stokes%euler%time,mesh%rr), error, navier_stokes%euler%communicator)
          CALL ns_l1_PAR(mesh, navier_stokes%euler%bc%sol_anal(n, navier_stokes%euler%time,mesh%rr), norm_anal, navier_stokes%euler%communicator)
          IF (norm_anal<1d-13) THEN
            norm = error
          ELSE
            norm = error/norm_anal
          END IF
        ELSE
          CALL ns_l1_PAR(mesh, un(:,n), norm, navier_stokes%euler%communicator)
        END IF
        tab_norm(n) = norm
      END DO
      write(*,*) 'tab_norm = ', tab_norm
      CALL get_num_test(num_test)
      CALL regression(tab_norm, opt_num_test=num_test)
    END IF

  END SUBROUTINE errors
END PROGRAM prog
