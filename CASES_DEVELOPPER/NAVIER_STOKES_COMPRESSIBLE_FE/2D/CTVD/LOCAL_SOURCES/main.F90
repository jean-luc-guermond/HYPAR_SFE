PROGRAM prog

  USE petscmpi
  USE start_setup_MODULE
  USE setup
  USE sub_plot
  USE my_util
  USE euler_post_proc_module
  USE plot_vtu_module
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
  time_backup   = navier_stokes%time + setup_data%checkpointing_freq
  time_snapshot = navier_stokes%time + setup_data%snapshot_freq


  WRITE(char, '(I5)') navier_stokes%euler%mesh%rank
  CALL plot_scalar_field(navier_stokes%euler%mesh%jj, navier_stokes%euler%mesh%rr, un(:, 1), 'initrho'//TRIM(ADJUSTL(char))//'.plt')

!=====================!
!==== SOLVER LOOP ====!
!=====================!
  CALL ETA%init(navier_stokes%mesh%rank, navier_stokes%time)
  tps = user_time()
  n = 0
  DO WHILE(navier_stokes%time < setup_data%final_time)
    !=== Update Navier-Stokes
    CALL navier_stokes%update(un)
    CALL ETA%update(navier_stokes%dt)
    n = n + 1
    
    !=== Post-proc
    IF (MOD(n, setup_data%verbose_freq)==0) THEN
        CALL ETA%print(navier_stokes%time, setup_data%final_time)
    END IF
    IF (navier_stokes%time > time_backup) THEN
        IF (navier_stokes%mesh%rank==0) THEN
          WRITE(*,*) 'overwriting Backup series'
        END IF
        time_backup = time_backup + setup_data%checkpointing_freq
        CALL RW%write_restart(mesh, navier_stokes%time, un, navier_stokes%name, opt_if_series=.FALSE.)
    END IF
    IF (navier_stokes%time > time_snapshot) THEN
        IF (navier_stokes%mesh%rank==0) THEN
          WRITE(*,*) 'Writing snapshot file ', RW%counter
        END IF
        time_snapshot = time_snapshot + setup_data%snapshot_freq
        CALL RW%write_restart(mesh, navier_stokes%time, un, navier_stokes%name, opt_if_series=.TRUE.)
    END IF

    !=== Stop loop
    IF (navier_stokes%dt < 1.d-12) THEN
      CALL error_petsc("BUG in main loop: timestep too small "//to_str(navier_stokes%dt))
    ELSEIF (n == setup_data%max_it) THEN
        IF (navier_stokes%mesh%rank==0) WRITE(*,*) "max_it reached, exiting solver loop"
        EXIT
    END IF
    !===
  END DO
  
  !=== final backup
  CALL RW%write_restart(mesh, navier_stokes%time, un, navier_stokes%name, opt_if_series=.FALSE.)

  tps = user_time() - tps

!=========================!
!==== POST-PROCESSING ====!
!=========================!

  CALL MPI_ALLREDUCE(navier_stokes%mesh%dom_np,tot_np,1,MPI_INTEGER,MPI_SUM,navier_stokes%euler%communicator,code)
  IF(navier_stokes%mesh%rank==0) THEN
     WRITE(*,*) ' tot_np', tot_np
     WRITE(*,*) ' Time per (time step times rk) per dof times proc', navier_stokes%mesh%nb_proc*tps/(tot_np*n*navier_stokes%euler%ERK%s), tps, n
  END IF
  ! CALL plot_scalar_field(navier_stokes%euler%mesh%jj, navier_stokes%euler%mesh%rr, un(:, 1), 'rho' // TRIM(ADJUSTL(char)) // '.plt')
  ! CALL plot_scalar_field(navier_stokes%euler%mesh%jj, navier_stokes%euler%mesh%rr, &
  ! navier_stokes%euler%bc%sol_anal(1, navier_stokes%euler%time,mesh%rr), 'rhoexact' // TRIM(ADJUSTL(char)) // '.plt')

  ALLOCATE(grad(navier_stokes%mesh%np))
  CALL schlieren(navier_stokes%euler,un(:, 1),grad)
  CALL make_vtu_file_2D(navier_stokes%communicator, navier_stokes%mesh, 'rho', un(:, 1), 'density', 'new', opt_it=0)
  CALL make_vtu_file_2D(navier_stokes%communicator, navier_stokes%mesh, 'rho_schlieren', grad, 'density_schlieren', 'new', opt_it=0)

call navier_stokes%profiler%output()
call navier_stokes%stokes%profiler%output()
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
    IF (options%if_regression) THEN
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
      CALL regression(tab_norm, opt_num_test=options%num_regex)
    END IF

  END SUBROUTINE errors
END PROGRAM prog
