MODULE stokes_parabolic_module
    USE Implicit_Butcher_tableau
#include "petsc/finclude/petsc.h"
    USE petsc
    USE def_type_mesh,                        ONLY: mesh_type, petsc_csr_LA
    USE read_inputs_module,                   ONLY: rec_length
    USE space_dim,                            ONLY: k_dim
    USE stokes_parabolic_matrices_module
    USE stokes_bc_arrays

    INTEGER, PRIVATE, PARAMETER :: METHOD_LUMPED=1, METHOD_FULL=2
    CHARACTER(LEN=20), DIMENSION(2), PRIVATE, PARAMETER  :: list_method = &
                [CHARACTER(LEN=20) :: 'lumped', 'full']

    TYPE argument_stokes_parabolic_type
        CHARACTER(LEN=rec_length) :: char_method          = '=== Which method to solve Stokes parabolic equation (lumped,full)? ==='
        CHARACTER(LEN=rec_length) :: ratio_reinit_precond_vel = '=== Reinit ratio for velocity preconditioner? ==='
    END TYPE argument_stokes_parabolic_type

    TYPE stokes_parabolic_type
        !===Parameters read from data
        INTEGER                       :: method              = METHOD_LUMPED
        CHARACTER(LEN=rec_length)     :: char_method         = 'lumped'
        REAL(KIND = 8)                :: thermal_diffusivity = 0.d0
        REAL(KIND = 8)                :: mu_viscosity        = 0.d0
        REAL(KIND = 8)                :: lambda_viscosity    = 0.d0
        REAL(KIND = 8)                :: cv                  = 2.5d0
        INTEGER                       :: ratio_reinit_precond_vel = 2.d0
        !===Parameters built along the way
        MPI_Comm :: communicator
        Vec          :: vel_x1vec, vel_x2vec
        Vec          :: temp_x1vec, temp_x2vec
        Vec          :: x4vec, x5vec !!!!!Conveniance vectors to be used only inside procedures!!!!
        Vec          :: sol_vel_vec, sol_temp_vec
        Vec, DIMENSION(:), ALLOCATABLE :: flux_rk_at_dof_vel, flux_rk_at_dof_temp
        CHARACTER(LEN=:), ALLOCATABLE :: name
        TYPE(mesh_type),     POINTER  :: mesh
        TYPE(petsc_csr_LA)            :: LA_vel, LA_temp
        TYPE(IBT), POINTER,   PUBLIC  :: IRK
        TYPE(stokes_parabolic_matrices_type), POINTER :: matrices
        TYPE(stokes_bc_type)          :: bc
        INTEGER                       :: irk_sv
        INTEGER                       :: syst_dim
        REAL(KIND = 8)                :: dt, time, final_time
    CONTAINS
        PROCEDURE :: read => read_stokes_parabolic_data
        PROCEDURE :: init => init_stokes_parabolic
        PROCEDURE, PUBLIC  :: one_step_IRK
        PROCEDURE, PRIVATE :: one_step_IRK_full
        PROCEDURE, PRIVATE :: init_vectors
        PROCEDURE, PRIVATE :: construct_stokes_bc
        PROCEDURE, PRIVATE :: iterative_LA
    END TYPE stokes_parabolic_type

CONTAINS
    SUBROUTINE init_stokes_parabolic(this, communicator, name, mesh, times)
        USE my_util,            ONLY: error_petsc, to_str
        USE space_dim
        USE st_matrix

        IMPLICIT NONE
        CLASS(stokes_parabolic_type), INTENT(INOUT) :: this
        MPI_Comm,                   INTENT(IN) :: communicator
        CHARACTER(100),             INTENT(IN) :: name
        TYPE(mesh_type), TARGET,    INTENT(IN) :: mesh
        REAL(KIND = 8), DIMENSION(2) :: times

        this%name = name
        this%mesh => mesh
        this%communicator = communicator
        CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, this%LA_temp)
        CALL st_aij_csr_glob_block_with_extra_layer(communicator, k_dim, mesh, this%LA_vel)

        this%time = times(1) !<==initial_time
        this%final_time = times(2) !<==final_time

        CALL this%read("STOKES PARABOLIC PARAMETERS FOR "//trim(adjustl(this%name)))

        ! === Build IRK structure
        CALL this%IRK%init(this%irk_sv)

        !===Matrices
        ALLOCATE(this%matrices)
        this%matrices%method              = this%method
        this%matrices%thermal_diffusivity = this%thermal_diffusivity
        CALL this%matrices%construct(this%communicator, this%mesh, this%LA_vel, this%LA_temp)

    END SUBROUTINE init_stokes_parabolic

    SUBROUTINE read_stokes_parabolic_data(this, section_name)
        USE read_inputs_module
        USE my_util, ONLY: get_tab_idx_char
        IMPLICIT NONE
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name

        CLASS(stokes_parabolic_type), INTENT(INOUT) :: this
        TYPE(argument_stokes_parabolic_type)        :: argument_data

        !================
        !=== MANDATORY Reading all data file
        !================
        IF (PRESENT(section_name)) THEN
            CALL read_data_init_list(section_name)
        ELSE
            CALL read_data_init_list()
        END IF

        !================
        !=== We now find the relevant information for this specific Stokes data
        !================

        !===Method
        CALL read_data(argument_data%char_method, this%char_method, opt_name=this%name)
        CALL get_tab_idx_char(this%char_method, list_method, this%method)

        !================
        !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
        !================

        CALL finalize_rewrite_data
    END SUBROUTINE read_stokes_parabolic_data

    SUBROUTINE one_step_IRK(this, stage, urk, flux_rk_at_dof_vel, flux_rk_at_dof_temp)
        USE my_util, ONLY : error_petsc, to_str
        IMPLICIT NONE
        CLASS(stokes_parabolic_type)                                                :: this
        REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%IRK%s+1)  :: urk
        REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%IRK%s)    :: flux_rk_at_dof_vel
        REAL(KIND = 8), DIMENSION(this%mesh%np, this%IRK%s)                   :: flux_rk_at_dof_temp 

        INTEGER, INTENT(IN) :: stage

        SELECT CASE(this%method)
        ! CASE(METHOD_LUMPED)
        !     CALL one_step_IRK_LUMPED(this,stage,urk,flux_rk_at_dof)
        CASE(METHOD_FULL)
            CALL one_step_IRK_FULL(this,stage,urk,flux_rk_at_dof_vel,flux_rk_at_dof_temp)
        CASE DEFAULT
            CALL error_petsc("wrong method in one_step_IRK "//to_str(this%method))
        END SELECT
    END SUBROUTINE one_step_IRK

    SUBROUTINE one_step_IRK_FULL(this,stage,urk,flux_rk_at_dof_vel,flux_rk_at_dof_temp)
        USE my_util, ONLY : error_petsc, to_str, user_time
        USE petsc_tools
        USE space_dim
        USE dir_nodes_petsc
        USE st_matrix
        IMPLICIT NONE
        CLASS(stokes_parabolic_type)                                                :: this
        REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%IRK%s+1), TARGET   :: urk
        REAL(KIND = 8), DIMENSION(:), POINTER                          :: rho
        REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim)                 :: velocity, vel_out
        REAL(KIND = 8), DIMENSION(k_dim*this%mesh%np)                  :: rhs
        REAL(KIND = 8), DIMENSION(this%mesh%np)                        :: rhs_temp, scal_temp, temp_out
        REAL(KIND = 8), DIMENSION(k_dim*this%mesh%np)                  :: lhs_mass
        REAL(KIND = 8), DIMENSION(k_dim*this%mesh%np, this%IRK%s)      :: flux_rk_at_dof_vel
        REAL(KIND = 8), DIMENSION(this%mesh%np, this%IRK%s)            :: flux_rk_at_dof_temp
        INTEGER, INTENT(IN) :: stage
        INTEGER :: k, l, np, stage_prime, ierr
        INTEGER,        POINTER :: count_vel, count_temp
        REAL(KIND = 8), POINTER :: tps_ratio_vel, tps_ratio_temp
        REAL(KIND = 8), POINTER :: tps_solver_ref_vel, tps_solver_ref_temp
        REAL(KIND = 8) :: local_tps_solver_ref, local_tps_solver_one_iter, tps_solver_one_iter, tps

        !=== Build counters
        count_vel          => this%matrices%elasticity_solver_param%count
        tps_ratio_vel      => this%matrices%elasticity_solver_param%tps_ratio
        tps_solver_ref_vel => this%matrices%elasticity_solver_param%tps_solver_ref
        !=== Build counters
        np = this%mesh%np

        !===Define density and velocity at stage (comes from Euler(stage))
        rho => urk(:, 1, stage)
        DO k = 1, k_dim
            velocity(1:np,k) = urk(:,k,stage)/rho!
        END DO
        !===end define density and velocity at stage 

        stage_prime = this%IRK%lp_of_l(stage) 

        !==========================!
        !======== VELOCITY ========!
        !==========================!

        !===Combine parabolic fluxes with IRK coefficients
        CALL VecMAXPY(this%vel_x1vec, stage-1, -this%dt*this%IRK%MatRK(stage,1:stage-1), this%flux_rk_at_dof_vel(1:stage-1), ierr) !<=== x1 receives sum of IRK fluxes
        DO k = 1, k_dim
            rhs((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*urk(:,k+1,stage_prime)
        END DO
        CALL array_to_petsc_vec(rhs, this%vel_x2vec, this%LA_vel, 'insert')
        CALL VecAXPY(this%vel_x1vec, 1.d0, this%vel_x2vec, ierr)
        !===Combine parabolic fluxes with IRK coefficients

        
        DO k = 1, k_dim
            lhs_mass((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*rho
        END DO
        CALL array_to_petsc_vec(lhs_mass, this%vel_x2vec, this%matrices%LA_vel, 'insert')
        
        IF (stage < this%IRK%s + 1) THEN
            
            !=== rhs BC
            DO k = 1, k_dim
                CALL dirichlet_rhs(this%matrices%LA_vel%loc_to_glob(k, this%bc%vel(k)%jsd)-1, &
                                this%bc%vit_anal(k, this%time, this%mesh%rr(:,this%bc%vel(k)%jsd)), this%vel_x1vec)
            END DO
            !=== rhs BC

            !=== LHS matrix construction + BC
            CALL MatCopy(this%matrices%vel_diff_mat, this%matrices%vel_mat, SAME_NONZERO_PATTERN, ierr)
            CALL MatScale(this%matrices%vel_mat, this%dt*this%IRK%MatRK(stage, stage), ierr)
            CALL MatDiagonalSet(this%matrices%vel_mat, this%vel_x2vec, ADD_VALUES, ierr)
            !CALL periodic_matrix_petsc(mesh%per, this%matrices%LA_vel, this%matrices%vel_mat)
            DO k = 1, k_dim
                CALL Dirichlet_M_parallel(this%matrices%vel_mat, this%LA_vel%loc_to_glob(k,this%bc%vel(k)%jsd))
            END DO
            !=== LHS matrix construction + BC

            !=== Solver linear system
            CALL this%iterative_LA(this%matrices%elasticity_solver_param, &
                                this%matrices%vel_mat, this%matrices%vel_ksp, this%matrices%precond_vel_mat,&
                                this%vel_x1vec, this%sol_vel_vec)
            ! IF (tps_ratio_vel>2) THEN
            !     count_vel = 0
            !     ! IF (this%mesh%rank==0) WRITE(*,*) 'reinitializing solver precond'
            !     CALL MatCopy(this%matrices%vel_mat, this%matrices%precond_vel_mat, SAME_NONZERO_PATTERN, ierr)
            !     CALL KSPDestroy(this%matrices%vel_ksp, ierr) !<=== FIXME
            !     CALL init_solver(this%communicator, this%matrices%elasticity_solver_param, this%matrices%vel_ksp, &
            !     this%matrices%vel_mat, opt_mat_pre=this%matrices%precond_vel_mat)!, opt_re_init=.TRUE.)  
            ! ELSE
            !     CALL KSPSetInitialGuessNonzero(this%matrices%vel_ksp, PETSC_TRUE, ierr)
            !     CALL KSPSetReusePreconditioner(this%matrices%vel_ksp, PETSC_TRUE, ierr)
            ! END IF

            ! tps = user_time()
            ! CALL solver(this%matrices%vel_ksp, this%vel_x1vec, this%sol_vel_vec, reinit = .FALSE., verbose = .FALSE.)
            ! count_vel = count_vel + 1

            ! SELECT CASE(count_vel)
            ! CASE(1)
            !     tps_ratio_vel = 1.d0
            ! CASE(2)
            !     local_tps_solver_ref = user_time()-tps
            !     CALL MPI_ALLREDUCE(local_tps_solver_ref, tps_solver_ref_vel, 1, MPI_DOUBLE, MPI_MIN, this%communicator, ierr)
            !     local_tps_solver_one_iter = user_time()-tps
            !     CALL MPI_ALLREDUCE(local_tps_solver_one_iter, tps_solver_one_iter, 1, MPI_DOUBLE, MPI_MIN, this%communicator, ierr)
            !     tps_ratio_vel = tps_solver_one_iter/tps_solver_ref_vel
            ! CASE DEFAULT
            !     local_tps_solver_one_iter = user_time()-tps
            !     CALL MPI_ALLREDUCE(local_tps_solver_one_iter, tps_solver_one_iter, 1, MPI_DOUBLE, MPI_MIN, this%communicator, ierr)
            !     tps_ratio_vel = tps_solver_one_iter/tps_solver_ref_vel
            ! END SELECT
            !=== Solver linear system
            CALL MatMult(this%matrices%vel_diff_mat, this%sol_vel_vec, this%flux_rk_at_dof_vel(stage), ierr)
        ELSE
            CALL VecPointWiseDivide(this%sol_vel_vec, this%vel_x1vec, this%vel_x2vec, ierr)
            CALL MatMult(this%matrices%vel_diff_mat, this%sol_vel_vec, this%flux_rk_at_dof_vel(1), ierr)
        END IF

        !====extract velocity and update momentum
        DO k=1, k_dim
            CALL extract_through_ghost(this%sol_vel_vec, k, k, this%LA_vel, vel_out(:, k), opt_assemble=.FALSE.)
            urk(:,k+1, stage) = rho*vel_out(:, k)
        END DO

        !============= FIXME COMPUTE VISCOUS PRODUCTION

        !=============================!
        !======== TEMPERATURE ========!
        !=============================!

        !===Combine parabolic fluxes with IRK coefficients
        CALL VecMAXPY(this%temp_x1vec, stage-1, -this%dt*this%IRK%MatRK(stage,1:stage-1), &
                                                this%flux_rk_at_dof_temp(1:stage-1), ierr) !<=== x1 receives sum of IRK fluxes
        scal_temp = urk(:,k_dim+2,stage_prime) - 0.5d0*rho*sum(vel_out**2,dim=2)

        rhs_temp = this%matrices%scal_lumped_mass*scal_temp
        CALL array_to_petsc_vec(rhs_temp, this%temp_x2vec, this%LA_temp, 'insert')
        CALL VecAXPY(this%temp_x1vec, 1.d0, this%temp_x2vec, ierr)
        !===Combine parabolic fluxes with IRK coefficients

        !===LHS 
        scal_temp = this%matrices%scal_lumped_mass*rho*this%cv
        CALL array_to_petsc_vec(scal_temp, this%temp_x2vec, this%matrices%LA_temp, 'insert')

        IF (stage < this%IRK%s + 1) THEN
            !=== rhs BC
            CALL dirichlet_rhs(this%matrices%LA_temp%loc_to_glob(1, this%bc%temp%jsd)-1, &
                            this%bc%temp_anal(this%time, this%mesh%rr(:,this%bc%temp%jsd)), this%temp_x1vec)
            !=== rhs BC

            !=== LHS matrix construction + BC
            CALL MatCopy(this%matrices%temp_diff_mat, this%matrices%temp_mat, SAME_NONZERO_PATTERN, ierr)
            CALL MatScale(this%matrices%temp_mat, this%dt*this%IRK%MatRK(stage, stage), ierr)
            CALL MatDiagonalSet(this%matrices%temp_mat, this%temp_x2vec, ADD_VALUES, ierr)
            !CALL periodic_matrix_petsc(mesh%per, this%matrices%LA_temp, this%matrices%temp_mat)
            CALL Dirichlet_M_parallel(this%matrices%temp_mat, this%LA_temp%loc_to_glob(1,this%bc%temp%jsd))
            !=== LHS matrix construction + BC
            CALL this%iterative_LA(this%matrices%temperature_solver_param, &
                                this%matrices%temp_mat, this%matrices%temp_ksp, this%matrices%precond_temp_mat,&
                                this%temp_x1vec, this%sol_temp_vec)
            CALL MatMult(this%matrices%temp_diff_mat, this%sol_temp_vec, this%flux_rk_at_dof_temp(stage), ierr)
        ELSE
            CALL VecPointWiseDivide(this%sol_temp_vec, this%temp_x1vec, this%temp_x2vec, ierr)
            CALL MatMult(this%matrices%temp_diff_mat, this%sol_temp_vec, this%flux_rk_at_dof_temp(1), ierr)
        END IF

        !====extract velocity and update momentum
        CALL extract_through_ghost(this%sol_temp_vec, 1, 1, this%LA_temp, temp_out, opt_assemble=.FALSE.)
        urk(:,k_dim+2, stage) = rho*this%cv*temp_out + 0.5d0*rho*sum(vel_out**2,dim=2)

    END SUBROUTINE one_step_IRK_FULL


    SUBROUTINE construct_stokes_bc(this, mesh, LA)
        USE petsc
#include "petsc/finclude/petsc.h"

        USE space_dim,           ONLY: k_dim
        IMPLICIT NONE
        CLASS(stokes_parabolic_type), INTENT(INOUT)        :: this
        TYPE(mesh_type)                            :: mesh
        TYPE(petsc_csr_LA)                         :: LA

        CALL this%bc%vel(1)%set(mesh, "ux")
        CALL this%bc%temp%set(mesh, "temperature")
        
        IF (k_dim>1) THEN
            CALL this%bc%vel(2)%set(mesh, "uy")
        END IF

    END SUBROUTINE construct_stokes_bc


   SUBROUTINE init_vectors(this)
        USE space_dim
        USE st_matrix, ONLY : create_my_ghost
        USE petsc
#include "petsc/finclude/petsc.h"

        IMPLICIT NONE
        CLASS(stokes_parabolic_type) :: this
        INTEGER, POINTER, DIMENSION(:) :: ifrom
        INTEGER :: n, ierr

        !=== Vel vectors
        CALL create_my_ghost(this%mesh, this%LA_vel, ifrom)
        CALL VecCreateGhost(this%communicator, k_dim*this%mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%vel_x1vec, ierr)
        CALL VecDuplicate(this%vel_x1vec, this%vel_x2vec, ierr)
        ALLOCATE(this%flux_rk_at_dof_vel(this%IRK%s))
        DO n=1, this%IRK%s
            CALL VecDuplicate(this%vel_x1vec, this%flux_rk_at_dof_vel(n), ierr)
        END DO
        CALL VecDuplicate(this%vel_x1vec, this%sol_vel_vec, ierr)
        !=== Vel vectors

        !=== Temperature vectors
        CALL create_my_ghost(this%mesh, this%LA_temp, ifrom)
        CALL VecCreateGhost(this%communicator, this%mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%temp_x1vec, ierr)
        CALL VecDuplicate(this%temp_x1vec, this%temp_x2vec, ierr)
        ALLOCATE(this%flux_rk_at_dof_temp(this%IRK%s))
        DO n=1, this%IRK%s
            CALL VecDuplicate(this%temp_x1vec, this%flux_rk_at_dof_temp(n), ierr)
        END DO
        CALL VecDuplicate(this%temp_x1vec, this%sol_temp_vec, ierr)
        !=== Temperature vectors


        ! CALL VecDuplicate(this%x1vec, this%x2vec, ierr)
        ! CALL VecDuplicate(this%x1vec, this%x3vec, ierr)
        ! CALL VecDuplicate(this%x1vec, this%x4vec, ierr)
        ! CALL VecDuplicate(this%x1vec, this%x5vec, ierr)

    END SUBROUTINE init_vectors

    SUBROUTINE iterative_LA(this, solver_param, matrix, matrix_ksp, matrix_precond, vec_rhs, vec_sol)
        USE my_util
        USE solver_petsc
        IMPLICIT NONE
        CLASS(stokes_parabolic_type) :: this
        TYPE(solver_data_type)       :: solver_param
        REAL(KIND = 8) :: tps, local_tps_solver_ref, local_tps_solver_one_iter, tps_solver_one_iter
        INTEGER :: ierr
        Mat :: matrix, matrix_precond
        KSP :: matrix_ksp
        Vec :: vec_rhs, vec_sol

        !=== Solver linear system
        IF (solver_param%tps_ratio>2) THEN
            solver_param%count = 0
            ! IF (this%mesh%rank==0) WRITE(*,*) 'reinitializing solver precond'
            CALL MatCopy(matrix, matrix_precond, SAME_NONZERO_PATTERN, ierr)
            CALL KSPDestroy(matrix_ksp, ierr) !<=== FIXME
            CALL init_solver(this%communicator, solver_param, matrix_ksp, &
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
            CALL MPI_ALLREDUCE(local_tps_solver_ref, solver_param%tps_solver_ref, 1, MPI_DOUBLE, MPI_MIN, this%communicator, ierr)
            local_tps_solver_one_iter = user_time()-tps
            CALL MPI_ALLREDUCE(local_tps_solver_one_iter, tps_solver_one_iter, 1, MPI_DOUBLE, MPI_MIN, this%communicator, ierr)
            solver_param%tps_ratio = tps_solver_one_iter/solver_param%tps_solver_ref
        CASE DEFAULT
            local_tps_solver_one_iter = user_time()-tps
            CALL MPI_ALLREDUCE(local_tps_solver_one_iter, tps_solver_one_iter, 1, MPI_DOUBLE, MPI_MIN, this%communicator, ierr)
            solver_param%tps_ratio = tps_solver_one_iter/solver_param%tps_solver_ref
        END SELECT

    END SUBROUTINE iterative_LA

END MODULE stokes_parabolic_module