PROGRAM test_matrix

  USE petsc
  USE fem_tn, ONLY: ns_l1_par
  USE fem_rhs
  USE start_setup_MODULE

  USE sub_plot
  USE post_processing_debug_MODULE
  USE setup
  USE my_util, ONLY: user_time
  USE space_dim
  USE st_matrix
  USE petsc_tools
  USE compute_periodic
  USE solver_petsc
  USE dir_nodes_petsc
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhs, exact_solution, un_out, lhs, mmt_contiguous
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ff_vel, mmt, u_out
  INTEGER,      DIMENSION(:), POINTER :: ifrom
  REAL(KIND=8) :: tps, time, dt, tps_solver=0.d0, tps_loop=0.d0, tps_init=0.d0, norm, error, tps_ratio
  REAL(KIND = 8) :: local_tps_solver_ref, tps_solver_ref, local_tps_solver_one_iter, tps_solver_one_iter
  INTEGER      :: tot_np, k, np, n, count
  CHARACTER(5) :: char

  TYPE(tVec) :: x1vec_vel, x2vec_vel, x3vec_vel, sol_vec
  INTEGER :: ierr


!========================!
!==== INITIALIZATION ====!
!========================!

  CALL start_setup
  CALL create_my_ghost(mesh, stokes_matrices%LA_vel, ifrom)
  CALL VecCreateGhost(communicator, k_dim*mesh%dom_np, &
      PETSC_DETERMINE, SIZE(ifrom), ifrom, x1vec_vel, ierr)
  CALL VecDuplicate(x1vec_vel, x2vec_vel, ierr)
  CALL VecDuplicate(x1vec_vel, x3vec_vel, ierr)
  CALL VecDuplicate(x1vec_vel, sol_vec, ierr)

  np = mesh%np
  ALLOCATE(lhs(mesh%np*k_dim))
  ALLOCATE(mmt(mesh%np,k_dim))
  ALLOCATE(mmt_contiguous(mesh%np*k_dim))
  ALLOCATE(ff_vel(mesh%np,k_dim))


  time = 0.d0
  dt = 1.d0/SQRT(1.d0*mesh%disp(mesh%nb_proc+1))

!=====================!
!==== SOLVER LOOP ====!
!=====================!
  
  tps_ratio = 1000.d0

  DO n=1, setup_data%max_it
    tps = user_time()

    time = time + dt
    !=== Build LHS operator
    DO k = 1, k_dim
        lhs((k-1)*np+1:k*np) = stokes_matrices%scal_lumped_mass*rho(mesh%rr, time)
    END DO
    CALL array_to_petsc_vec(lhs, x1vec_vel, stokes_matrices%LA_vel, 'insert')
    
    CALL MatCopy(stokes_matrices%vel_diff_mat, stokes_matrices%vel_mat, SAME_NONZERO_PATTERN, ierr)
    CALL MatScale(stokes_matrices%vel_mat, dt, ierr)
    CALL MatDiagonalSet(stokes_matrices%vel_mat, x1vec_vel, ADD_VALUES, ierr)
    !CALL periodic_matrix_petsc(mesh%per, stokes_matrices%LA_vel, stokes_matrices%vel_mat)
    DO k = 1, k_dim
      CALL Dirichlet_M_parallel(stokes_matrices%vel_mat, LA_vel%loc_to_glob(k,dir%jsd))
    END DO
    !=== Build LHS operator

    !=== Build RHS
    mmt_contiguous(1:mesh%np) =  momentum(1, mesh%rr, time) * stokes_matrices%scal_lumped_mass
    mmt_contiguous(mesh%np+1:) = momentum(2, mesh%rr, time) * stokes_matrices%scal_lumped_mass
    CALL array_to_petsc_vec(mmt_contiguous, x2vec_vel, stokes_matrices%LA_vel, 'insert')

    ff_vel = dt * source(mesh%rr,setup_data%mu_viscosity,setup_data%lambda_viscosity)
    CALL qs_00_block (mesh, LA_vel, ff_vel, x1vec_vel)
    CALL VecAXPY(x1vec_vel, 1.d0, x2vec_vel, ierr)
    DO k = 1, k_dim
      CALL dirichlet_rhs(stokes_matrices%LA_vel%loc_to_glob(k, dir%jsd)-1, velocity(k, mesh%rr(:,dir%jsd)), x1vec_vel)
    END DO
    !=== Build RHS
    IF (n>1) tps_loop = tps_loop + (user_time()-tps)
    ! !=== Solve
    ! IF (tps_ratio>2) THEN
    !   count = 0
    !   IF (mesh%rank==0) WRITE(*,*) 'reinitializing solver precond'
    !   tps = user_time()
    !   CALL MatCopy(stokes_matrices%vel_mat, stokes_matrices%precond_vel_mat, SAME_NONZERO_PATTERN, ierr)
    !   CALL KSPDestroy(stokes_matrices%vel_ksp, ierr)
    !   ! IF (once) THEN
    !   !   once = .FALSE.
    !     CALL init_solver(communicator, stokes_matrices%elasticity_solver_param, stokes_matrices%vel_ksp, &
    !     stokes_matrices%vel_mat, opt_mat_pre=stokes_matrices%precond_vel_mat)!, opt_re_init=.TRUE.)  
    !   ! END IF  
    !   tps_init = user_time() - tps
    ! ELSE
    !   CALL KSPSetInitialGuessNonzero(stokes_matrices%vel_ksp, PETSC_TRUE, ierr)
    !   CALL KSPSetReusePreconditioner(stokes_matrices%vel_ksp, PETSC_TRUE, ierr)
    ! END IF
    ! IF (n>1) tps_loop = tps_loop + (user_time()-tps)

    ! tps = user_time()
    ! CALL solver(stokes_matrices%vel_ksp, x1vec_vel, sol_vec, reinit = .FALSE., verbose = .FALSE.)
    ! count = count + 1

    ! SELECT CASE(count)
    ! CASE(1)
    !   tps_ratio = 1.d0
    ! CASE(2)
    !   local_tps_solver_ref = user_time()-tps
    !   CALL MPI_ALLREDUCE(local_tps_solver_ref, tps_solver_ref, 1, MPI_DOUBLE, MPI_MIN, communicator, ierr)
    !   local_tps_solver_one_iter = user_time()-tps
    !   CALL MPI_ALLREDUCE(local_tps_solver_one_iter, tps_solver_one_iter, 1, MPI_DOUBLE, MPI_MIN, communicator, ierr)
    !   tps_ratio = tps_solver_one_iter/tps_solver_ref
    ! CASE DEFAULT
    !   local_tps_solver_one_iter = user_time()-tps
    !   CALL MPI_ALLREDUCE(local_tps_solver_one_iter, tps_solver_one_iter, 1, MPI_DOUBLE, MPI_MIN, communicator, ierr)
    !   tps_ratio = tps_solver_one_iter/tps_solver_ref
    ! END SELECT
    tps = user_time()
    CALL iterative_LA(stokes_matrices%elasticity_solver_param, stokes_matrices%vel_mat, stokes_matrices%vel_ksp, &
    stokes_matrices%precond_vel_mat, x1vec_vel, sol_vec)

    IF (n>1) THEN
      tps_solver = tps_solver + (user_time()-tps)
    END IF

    IF (mesh%rank==0) write(*,*) 'n=', n, (user_time()-tps)/mesh%np
    !=== Solve
  END DO

!=========================!
!==== POST-PROCESSING ====!
!=========================!
  
  WRITE(char, '(I5)') mesh%rank
  IF (mesh%rank==0) THEN
    WRITE(*,*) 'np_tot = ', mesh%disp(mesh%nb_proc+1)
    WRITE(*,*) 'np_loc = ', mesh%np
    WRITE(*,*) 'init = ', tps_init, tps_init/mesh%np
    WRITE(*,*) 'solver = ', tps_solver/(setup_data%max_it-1), tps_solver/(setup_data%max_it-1)/mesh%np
    WRITE(*,*) 'loop = ', tps_loop/(setup_data%max_it-1), tps_loop/(setup_data%max_it-1)/mesh%np
    WRITE(*,*) 'HYPRE strong threshold = ', elasticity_solver_param%boomeramg_strong_threshold
  END IF

  ALLOCATE(u_out(mesh%np,k_dim))
  DO k=1, k_dim
    CALL extract_through_ghost(sol_vec, k, k, stokes_matrices%LA_vel, u_out(:, k), opt_assemble=.FALSE.)
  END DO

  CALL plot_scalar_field(mesh%jj, mesh%rr, u_out(:,1), 'ux_' // TRIM(ADJUSTL(char)) // '.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, u_out(:,2), 'uy_' // TRIM(ADJUSTL(char)) // '.plt')
  
  IF (setup_data%if_analytical_ref) THEN
    CALL plot_scalar_field(mesh%jj, mesh%rr, u_out(:,1)-velocity(1,mesh%rr), 'errorx_' // trim(adjustl(char)) // '.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr,u_out(:,2)-velocity(2,mesh%rr), 'errory_' // trim(adjustl(char)) // '.plt')
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
    REAL(KIND=8), DIMENSION(k_dim) :: tab_norm
    INTEGER      :: n, num_test, k

!==== Put final processing stuff here ====!
    DO k=1, k_dim
      IF (setup_data%if_analytical_ref) THEN
        CALL ns_l1_PAR(mesh, u_out(:,k)-velocity(k,mesh%rr), error, communicator)
        CALL ns_l1_PAR(mesh, velocity(k,mesh%rr), norm_anal, communicator)
        norm = error/norm_anal
        IF(mesh%rank==0) WRITE(*, *) 'Comp = ',k,'; Relative error, L1-norm = ', error/norm_anal
        ! IF(Laplace%mesh%rank==0) WRITE(*, '(A,I0,A,g12.3)') 'Comp = ',n,'; Relative error, L1-norm = ', error/norm_anal
      ELSE
        CALL ns_l1_PAR(mesh, u_out(:,k), norm, communicator)
        IF(mesh%rank==0) WRITE(*, *) 'Comp = ',k,'; no analytical ref, L1-norm = ', norm
      END IF
    END DO

!==== For regression tests ====!
    IF (setup_data%if_regression_test) THEN
      DO k=1, k_dim
        IF (setup_data%if_analytical_ref) THEN
          CALL ns_l1_PAR(mesh, u_out(:,k)-velocity(k,mesh%rr), error, communicator)
          CALL ns_l1_PAR(mesh, velocity(k,mesh%rr), norm_anal, communicator)
          norm = error/norm_anal
        ELSE
          CALL ns_l1_PAR(mesh, u_out(:,k), norm, communicator)
        END IF
        tab_norm(n) = norm
      END DO
      CALL get_num_test(num_test)
      CALL regression(tab_norm, opt_num_test=num_test)
    END IF

  END SUBROUTINE errors


    SUBROUTINE iterative_LA(solver_param, matrix, matrix_ksp, matrix_precond, vec_rhs, vec_sol)
        USE my_util
        USE solver_petsc
        IMPLICIT NONE
        TYPE(solver_data_type)       :: solver_param
        REAL(KIND = 8) :: tps, local_tps_solver_ref, local_tps_solver_one_iter, tps_solver_one_iter
        INTEGER :: ierr
        TYPE(tMat) :: matrix, matrix_precond
        TYPE(tKSP) :: matrix_ksp
        TYPE(tVec) :: vec_rhs, vec_sol

        !=== Solver linear system
        IF (solver_param%tps_ratio>2) THEN
            solver_param%count = 0
            ! IF (this%mesh%rank==0) WRITE(*,*) 'reinitializing solver precond'
            CALL MatCopy(matrix, matrix_precond, SAME_NONZERO_PATTERN, ierr)
            CALL KSPDestroy(matrix_ksp, ierr) !<=== FIXME
            CALL init_solver(communicator, solver_param, matrix_ksp, &
            matrix, opt_mat_pre=matrix_precond)!, opt_re_init=.TRUE.)  
        ELSE
            CALL KSPSetInitialGuessNonzero(matrix_ksp, PETSC_TRUE, ierr)
            CALL KSPSetReusePreconditioner(matrix_ksp, PETSC_TRUE, ierr)
        END IF

        tps = user_time()
        CALL solver(matrix_ksp, vec_rhs, vec_sol, reinit = .FALSE., verbose = .FALSE.)
        solver_param%count = solver_param%count + 1

        SELECT CASE(solver_param%count)
        CASE(1)
            solver_param%tps_ratio = 1.d0
        CASE(2)
            local_tps_solver_ref = user_time()-tps
            CALL MPI_ALLREDUCE(local_tps_solver_ref, solver_param%tps_solver_ref, 1, MPI_DOUBLE, MPI_MIN, communicator, ierr)
            local_tps_solver_one_iter = user_time()-tps
            CALL MPI_ALLREDUCE(local_tps_solver_one_iter, tps_solver_one_iter, 1, MPI_DOUBLE, MPI_MIN, communicator, ierr)
            solver_param%tps_ratio = tps_solver_one_iter/solver_param%tps_solver_ref
        CASE DEFAULT
            local_tps_solver_one_iter = user_time()-tps
            CALL MPI_ALLREDUCE(local_tps_solver_one_iter, tps_solver_one_iter, 1, MPI_DOUBLE, MPI_MIN, communicator, ierr)
            solver_param%tps_ratio = tps_solver_one_iter/solver_param%tps_solver_ref
        END SELECT

    END SUBROUTINE iterative_LA


END PROGRAM test_matrix
