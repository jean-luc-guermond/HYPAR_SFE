PROGRAM prog

  USE petscmpi
  USE start_setup_MODULE
  USE setup
  USE sub_plot
  USE euler_post_proc_module
  USE plot_vtu_module
  USE my_util
  USE ETA_module
  IMPLICIT NONE
  REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: grad
  REAL(KIND = 8) :: tps, time_backup, time_snapshot
  CHARACTER(5) :: char
  INTEGER :: n, tot_np, code, num_test
  TYPE(ETA_type) :: ETA

!========================!
!==== INITIALIZATION ====!
!========================!

  CALL start_setup
  time_backup   = euler%time + setup_data%checkpointing_freq
  time_snapshot = euler%time + setup_data%snapshot_freq

  WRITE(char, '(I5)') euler%mesh%rank
  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:, 1), 'initrho'//TRIM(ADJUSTL(char))//'.plt')

!=====================!
!==== SOLVER LOOP ====!
!=====================!

  CALL ETA%init(euler%mesh%rank, euler%time)
  tps = user_time()
  n = 0
  DO WHILE(euler%time < setup_data%final_time)
    !===  Update Euler
    CALL euler%update(un)
    CALL ETA%update(euler%dt)
    n = n + 1

    !=== Post-proc
    IF (MOD(n, setup_data%verbose_freq)==0) THEN
        CALL ETA%print(euler%time, setup_data%final_time)
    END IF
    IF (euler%time > time_backup) THEN
        IF (euler%mesh%rank==0) THEN
          WRITE(*,*) 'overwriting Backup series'
        END IF
        time_backup = time_backup + setup_data%checkpointing_freq
        CALL RW%write_restart(mesh, euler%time, un, euler%name, opt_if_series=.FALSE.)
    END IF
    IF (euler%time > time_snapshot) THEN
        IF (euler%mesh%rank==0) THEN
          WRITE(*,*) 'Writing snapshot file ', RW%counter
        END IF
        time_snapshot = time_snapshot + setup_data%snapshot_freq
        CALL RW%write_restart(mesh, euler%time, un, euler%name, opt_if_series=.TRUE.)
    END IF

    !=== Stop loop
    IF (euler%dt < 1.d-12) THEN
      CALL error_petsc("BUG in main loop: timestep too small "//to_str(euler%dt))
    ELSEIF (n == setup_data%max_it) THEN
        IF (euler%mesh%rank==0) WRITE(*,*) "max_it reached, exiting solver loop"
        EXIT
    END IF
    !===
  END DO
  tps = user_time() - tps

!=========================!
!==== POST-PROCESSING ====!
!=========================!

  CALL MPI_ALLREDUCE(euler%mesh%dom_np,tot_np,1,MPI_INTEGER,MPI_SUM,euler%communicator,code)
  IF(euler%mesh%rank==0) THEN
     WRITE(*,*) ' tot_np', tot_np
     WRITE(*,*) ' Time per time step per dof times proc', euler%mesh%nb_proc*tps/(tot_np*n*(euler%ERK%s)), tps, n*(euler%ERK%s)
  END IF

  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:, 1), 'rho' // TRIM(ADJUSTL(char)) // '.plt')

  IF (setup_data%if_final_post_proc) THEN
    ALLOCATE(grad(euler%mesh%np))
    CALL schlieren(euler,un(:, 1),grad)
    CALL make_vtu_file_2D(euler%communicator, euler%mesh, 'rho', un(:, 1), 'density', 'new', opt_it=0)
    CALL make_vtu_file_2D(euler%communicator, euler%mesh, 'rho_schlieren', grad, 'density_schlieren', 'new', opt_it=0)
  END IF

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
        CALL ns_l1_PAR(mesh, un(:,n)-euler%bc%sol_anal(n, euler%time,mesh%rr), error, euler%communicator)
        CALL ns_l1_PAR(mesh, euler%bc%sol_anal(n, euler%time,mesh%rr), norm_anal, euler%communicator)
        norm = error/norm_anal
        IF(euler%mesh%rank==0) WRITE(*, *) 'Comp ', euler%name_comp(n), '; Relative error, L1-norm = ', error/norm_anal
        ! IF(euler%mesh%rank==0) WRITE(*, '(A,I0,A,g12.3)') 'Comp = ',n,'; Relative error, L1-norm = ', error/norm_anal
      ELSE
        CALL ns_l1_PAR(mesh, un(:,n), norm, euler%communicator)
        IF(euler%mesh%rank==0) WRITE(*, *) 'Comp ', euler%name_comp(n), '; no analytical ref, L1-norm = ', norm
      END IF
    END DO

    
!==== For regression tests ====!
    IF (options%if_regression) THEN
      DO n=1, SIZE(un,2)
        IF (setup_data%if_analytical_ref) THEN
          CALL ns_l1_PAR(mesh, un(:,n)-euler%bc%sol_anal(n, euler%time,mesh%rr), error, euler%communicator)
          CALL ns_l1_PAR(mesh, euler%bc%sol_anal(n, euler%time,mesh%rr), norm_anal, euler%communicator)
          IF (norm_anal<1d-13) THEN
            norm = error
          ELSE
            norm = error/norm_anal
          END IF
        ELSE
          CALL ns_l1_PAR(mesh, un(:,n), norm, euler%communicator)
        END IF
        tab_norm(n) = norm
      END DO
      write(*,*) 'tab_norm = ', tab_norm
      CALL regression(tab_norm, opt_num_test=options%num_regex)
    END IF

  END SUBROUTINE errors
END PROGRAM prog
