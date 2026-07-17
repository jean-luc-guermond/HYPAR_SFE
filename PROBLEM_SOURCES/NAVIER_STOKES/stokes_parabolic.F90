
MODULE stokes_parabolic_module
    USE Implicit_Butcher_tableau
    USE petscvec
    USE def_type_mesh,                        ONLY: mesh_type
    USE petsc_csr_LA_module,                  ONLY: petsc_csr_LA

    USE read_inputs_module,                   ONLY: rec_length
    USE space_dim,                            ONLY: k_dim
    USE stokes_parabolic_matrices_module
    USE stokes_bc_arrays
USE profiler_module

    PUBLIC :: stokes_parabolic_type, function_template_temperature, template_forcing

    PRIVATE

    INTEGER, PRIVATE, PARAMETER :: METHOD_LUMPED=1, METHOD_FULL=2
    CHARACTER(LEN=20), DIMENSION(2), PRIVATE, PARAMETER  :: list_method = &
                [CHARACTER(LEN=20) :: 'lumped', 'full']

    ABSTRACT INTERFACE
        FUNCTION function_template_temperature(rho, ie) RESULT(vv)
            IMPLICIT NONE
            REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, ie
            REAL(KIND = 8), DIMENSION(SIZE(rho, 1))  :: vv
        END FUNCTION function_template_temperature

        FUNCTION template_forcing(comp, rr, time) RESULT(vv)
            IMPLICIT NONE
            REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
            REAL(KIND=8), INTENT(IN) :: time
            INTEGER,      INTENT(IN) :: comp
            REAL(KIND = 8), DIMENSION(SIZE(rr, 2))  :: vv
        END FUNCTION template_forcing
    END INTERFACE

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
        INTEGER :: communicator
        TYPE(tVec)          :: vel1_vec, vel2_vec
        TYPE(tVec)          :: temp1_vec, temp2_vec
        TYPE(tVec)          :: x4vec, x5vec !!!!!Conveniance vectors to be used only inside procedures!!!!
        TYPE(tVec)          :: sol_vel_vec, sol_temp_vec
        TYPE(tVec), DIMENSION(:), ALLOCATABLE :: vel_flux_rk_at_dof, temp_flux_rk_at_dof, forcing_rk_at_dof
        CHARACTER(LEN=:), ALLOCATABLE :: name
        TYPE(mesh_type),     POINTER  :: mesh
        TYPE(petsc_csr_LA)            :: LA_vel, LA_temp
        TYPE(IBT),            PUBLIC  :: IRK
        TYPE(stokes_parabolic_matrices_type) :: matrices
        TYPE(stokes_bc_type)          :: bc
        INTEGER                       :: irk_sv
        INTEGER                       :: syst_dim
        REAL(KIND = 8)                :: dt, time, final_time
        PROCEDURE(function_template_temperature),  NOPASS, POINTER :: temperature => NULL()
        PROCEDURE(template_forcing),               NOPASS, POINTER :: forcing => NULL()
type(profiler_type) :: profiler
    CONTAINS
        PROCEDURE :: read => read_stokes_parabolic_data
        PROCEDURE :: init => init_stokes_parabolic
        PROCEDURE, PUBLIC :: set_times =>set_times_stokes 
        PROCEDURE, PUBLIC  :: one_step_IRK, update
        PROCEDURE, PRIVATE :: one_step_IRK_lumped
        ! PROCEDURE, PRIVATE :: one_step_IRK_full
        PROCEDURE, PRIVATE :: init_vectors
        PROCEDURE, PRIVATE :: construct_stokes_bc
        PROCEDURE, PRIVATE :: iterative_LA
    END TYPE stokes_parabolic_type

CONTAINS

!===================================
!==== INITIALIZATIONS ==============
!===================================


    SUBROUTINE init_stokes_parabolic(this, communicator, name, mesh)
        USE my_util,            ONLY: error_petsc, to_str
        USE space_dim
        USE st_matrix

        IMPLICIT NONE
        CLASS(stokes_parabolic_type), INTENT(INOUT) :: this
        INTEGER,                   INTENT(IN) :: communicator
        CHARACTER(LEN=*),             INTENT(IN) :: name
        TYPE(mesh_type), TARGET,    INTENT(IN) :: mesh
character(len=10), dimension(:), allocatable :: list_profiler

        this%name = name
        this%mesh => mesh
        this%communicator = communicator
        this%syst_dim = k_dim + 2
        
        !===Init Dirichlet BC
        CALL this%construct_stokes_bc(mesh)

        CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, this%LA_temp)
        CALL st_aij_csr_glob_block_with_extra_layer(communicator, k_dim, mesh, this%LA_vel)

        CALL this%read("STOKES PARABOLIC PARAMETERS FOR "//trim(adjustl(this%name)))

        ! === Build IRK structure
        CALL this%IRK%init(this%irk_sv)

        !===Matrices/vectors
        this%matrices%method              = this%method
        this%matrices%thermal_diffusivity = this%thermal_diffusivity
        this%matrices%mu_viscosity = this%mu_viscosity
        this%matrices%lambda_viscosity = this%lambda_viscosity

        CALL this%matrices%construct(this%communicator, this%mesh, this%LA_vel, this%LA_temp)

        CALL this%init_vectors

        !===Profiling
        ALLOCATE(list_profiler(16))
        list_profiler(1) = 'stokes tot'
        list_profiler(2) = 'vel dissip'
        list_profiler(3) = 'vel flux'
        list_profiler(4) = 'lhs vel lump'
        list_profiler(5) = 'rhs Dir'
        list_profiler(6) = 'lhs vel mat'
        list_profiler(7) = 'vel solve'
        list_profiler(8) = 'vel extract'
        list_profiler(9) = 'temp dissip'
        list_profiler(10) = 'temp flux'
        list_profiler(11) = 'rhs temp lump'
        list_profiler(12) = 'lhs temp lump'
        list_profiler(13) = 'temp Dir'
        list_profiler(14) = 'lhs temp mat'
        list_profiler(15) = 'temp solve'
        list_profiler(16) = 'temp extract'
        CALL this%profiler%init(list_profiler)
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

    SUBROUTINE set_times_stokes(this, times)
        IMPLICIT NONE
        CLASS(stokes_parabolic_type)              :: this
        REAL(KIND=8), DIMENSION(2), INTENT(IN)    :: times
        this%time = times(1) !<==initial_time
        this%final_time = times(2) !<==final_time
    END SUBROUTINE set_times_stokes

    SUBROUTINE construct_stokes_bc(this, mesh)
        USE space_dim,           ONLY: k_dim
        IMPLICIT NONE
        CLASS(stokes_parabolic_type), INTENT(INOUT)        :: this
        TYPE(mesh_type)                            :: mesh

        CALL this%bc%vel(1)%set(mesh, "ux "//TRIM(ADJUSTL(this%name)), "DIRICHLET BC PARAMETERS FOR "//TRIM(ADJUSTL(this%name)))
        CALL this%bc%temp%set(mesh, "temperature "//TRIM(ADJUSTL(this%name)))
        
        IF (k_dim>1) THEN
            CALL this%bc%vel(2)%set(mesh, "uy "//TRIM(ADJUSTL(this%name)))
        END IF

    END SUBROUTINE construct_stokes_bc


   SUBROUTINE init_vectors(this)
        USE space_dim
        USE st_matrix, ONLY : create_my_ghost

        IMPLICIT NONE
        CLASS(stokes_parabolic_type) :: this
        INTEGER, POINTER, DIMENSION(:) :: ifrom
        INTEGER :: n, ierr

        !=== Vel vectors
        CALL create_my_ghost(this%mesh, this%LA_vel, ifrom)
        CALL VecCreateGhost(this%communicator, k_dim*this%mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%vel1_vec, ierr)
        CALL VecDuplicate(this%vel1_vec, this%vel2_vec, ierr)
        ALLOCATE(this%vel_flux_rk_at_dof(this%IRK%s))
        DO n=1, this%IRK%s
            CALL VecDuplicate(this%vel1_vec, this%vel_flux_rk_at_dof(n), ierr)
        END DO
        CALL VecDuplicate(this%vel1_vec, this%sol_vel_vec, ierr)
        CALL VecZeroEntries(this%sol_vel_vec, ierr)
        !=== Vel vectors

        !=== Temperature vectors
        CALL create_my_ghost(this%mesh, this%LA_temp, ifrom)
        CALL VecCreateGhost(this%communicator, this%mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%temp1_vec, ierr)
        CALL VecDuplicate(this%temp1_vec, this%temp2_vec, ierr)
        ALLOCATE(this%temp_flux_rk_at_dof(this%IRK%s))
        DO n = 1, this%IRK%s
            CALL VecDuplicate(this%temp1_vec, this%temp_flux_rk_at_dof(n), ierr)
        END DO
        CALL VecDuplicate(this%temp1_vec, this%sol_temp_vec, ierr)
        !=== Temperature vectors

        !=== Forcing vectors
        CALL create_my_ghost(this%mesh, this%LA_vel, ifrom)
        ALLOCATE(this%forcing_rk_at_dof(this%IRK%s))
        DO n = 1, this%IRK%s
            CALL VecDuplicate(this%vel1_vec, this%forcing_rk_at_dof(n), ierr)
        END DO
        !=== Forcing vectors

    END SUBROUTINE init_vectors

!===================================
!==== GENERAL UPDATES ==============
!===================================

    SUBROUTINE update(this, un_in)
        !> IMPLICIT UPDATE of Stokes
        !! This subroutine is only used when Stokes is decoupled from any other problem (e.g no Euler)
        !! User requirements:
        !!     rho
        !!     dt
        !! TODO:
        !!     test of temperature stepping
        USE space_dim
        IMPLICIT NONE
        CLASS(stokes_parabolic_type)                                       :: this
        REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT):: un_in
        REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim, this%IRK%s+1) :: urk
        INTEGER  :: stage, stage_prime, k

        urk = 0.d0
        urk(:,:,1) = un_in
        DO stage = 2, this%IRK%s+1
            stage_prime = this%IRK%lp_of_l(stage) 
            urk(:,1,stage)     = this%bc%rho_imposed(this%time + this%IRK%C(stage)*this%dt, this%mesh%rr) !=== recompute rho
            DO k=1, k_dim
                !=== rescale momentum with new density (we are solving for vel, not momentum!!!)
                urk(:,k+1,stage)   = urk(:,k+1,stage_prime) /  urk(:,1,stage_prime) * urk(:,1,stage)
            END DO
            !=== We do not care about energy, so simply copy previous step
            urk(:,k_dim+2,stage) = urk(:,k_dim+2,stage-1) 

            CALL this%one_step_IRK(stage,urk)
        END DO
        un_in = urk(:,:,this%IRK%s+1)
        this%time = this%time + this%dt
    END SUBROUTINE update

    SUBROUTINE one_step_IRK(this, stage, urk)
        USE my_util, ONLY : error_petsc, to_str
        IMPLICIT NONE
        CLASS(stokes_parabolic_type)                                                :: this
        REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%IRK%s+1)  :: urk

        INTEGER, INTENT(IN) :: stage

        SELECT CASE(this%method)
        CASE(METHOD_LUMPED)
            CALL one_step_IRK_LUMPED(this,stage,urk)
        CASE(METHOD_FULL)
            CALL error_petsc("BUG in stokes%one_step_IRK: method full not validated on temperature yet, please use lumped")
            CALL one_step_IRK_FULL(this,stage,urk)
        CASE DEFAULT
            CALL error_petsc("wrong method in one_step_IRK "//to_str(this%method))
        END SELECT
    END SUBROUTINE one_step_IRK

!===================================
!==== NUMERICAL SCHEMES ============
!===================================

    SUBROUTINE one_step_IRK_LUMPED(this,stage,urk)
        USE my_util, ONLY : error_petsc, to_str, user_time
        USE petsc_tools
        USE space_dim
        USE dir_nodes_petsc
        USE st_matrix
        USE fem_rhs, ONLY: qs_00_block
        USE sub_plot
        USE fem_tn, ONLY: ns_l1_par
        IMPLICIT NONE
        CLASS(stokes_parabolic_type)                                                :: this
        REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%IRK%s+1), TARGET   :: urk
        REAL(KIND = 8), DIMENSION(:), POINTER                          :: rho_lm1, rho_l
        REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim)                 :: velocity_lm1, vel_out, ff_vel
        REAL(KIND = 8), DIMENSION(k_dim*this%mesh%np)                  :: rhs
        REAL(KIND = 8), DIMENSION(this%mesh%np)                        :: rhs_temp, scal_temp, temp_out
        REAL(KIND = 8), DIMENSION(k_dim*this%mesh%np)                  :: lhs_mass
        INTEGER, INTENT(IN) :: stage
        INTEGER :: k, l, np, ierr!, stage_prime
        REAL(KIND = 8) :: local_tps_solver_ref, local_tps_solver_one_iter, tps_solver_one_iter, tps, time_stage, error, norm, max_vel_loc, max_vel
        !=== Build pointers
        rho_lm1 => urk(:, 1, stage-1)
        rho_l   => urk(:, 1, stage)
        !=== Build pointers
        np = this%mesh%np
        ! stage_prime = this%IRK%lp_of_l(stage) 

        IF (stage==2) THEN
            DO l=1, this%IRK%s
                CALL VecZeroEntries(this%vel_flux_rk_at_dof(l), ierr)
            END DO
        END IF

        !=== Init rhs vector for velocity problem
        CALL VecZeroEntries(this%vel1_vec, ierr)

        !=============================================!
        !======== FORCING CONTRIBUTION AT STAGE 2 ====!
        !=============================================!

        !=== Forcing ===!
        IF (ASSOCIATED(this%forcing)) THEN
            IF (stage == 2) THEN
                time_stage = this%time+this%IRK%C(stage-1)*this%dt

                !=== Lumped version
                DO k=1, k_dim
                    rhs((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*this%forcing(k, this%mesh%rr, time_stage)
                END DO
                CALL array_to_petsc_vec(-rhs, this%vel_flux_rk_at_dof(stage-1), this%LA_vel, 'add', opt_include_ghost=.FALSE.)
            END IF
        END IF
        !=== Forcing ===!

        !========================================================!
        !======== VELOCITY VISCOUS DISSIPATION at stage-1========!
        !========================================================!

        !===Define velocity at stage-1
        DO k = 1, k_dim
            velocity_lm1(:,k) = urk(:,k+1,stage-1)/rho_lm1
            rhs((k-1)*np+1:k*np) = velocity_lm1(:,k)
        END DO
        CALL array_to_petsc_vec(rhs, this%vel2_vec, this%LA_vel, 'insert')
        !=== (-1)Div(sigma(vel))
        CALL MatMultAdd(this%matrices%vel_diff_mat, this%vel2_vec, this%vel_flux_rk_at_dof(stage-1), this%vel_flux_rk_at_dof(stage-1), ierr) 

        !================================!
        !======== VELOCITY UPDATE========!
        !================================!
        !===Combine parabolic fluxes with IRK coefficients (notice (-1)*dt*IRK%MatRK)
        IF (stage < this%IRK%s + 1 .AND. ASSOCIATED(this%forcing)) THEN
            time_stage = this%time+this%IRK%C(stage)*this%dt
            !=== Lumped version
            DO k=1, k_dim
                rhs((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*this%forcing(k, this%mesh%rr, time_stage)
            END DO
            CALL array_to_petsc_vec(-rhs, this%vel_flux_rk_at_dof(stage), this%LA_vel, 'add', opt_include_ghost=.FALSE.)

            CALL VecMAXPY(this%vel1_vec, stage, -this%dt*this%IRK%MatRK(stage,1:stage), this%vel_flux_rk_at_dof(1:stage), ierr) !<=== x1 receives sum of IRK fluxes
        ELSE
            CALL VecMAXPY(this%vel1_vec, stage-1, -this%dt*this%IRK%MatRK(stage,1:stage-1), this%vel_flux_rk_at_dof(1:stage-1), ierr) !<=== x1 receives sum of IRK fluxes
        END IF


        DO k = 1, k_dim
            !=== WARNING HERE: stage_prime involved in explicit step or must be set manually if no prior explicit step 
            ! rhs((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*urk(:,k+1,stage_prime)
            rhs((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*urk(:,k+1,stage)
        END DO
        CALL array_to_petsc_vec(rhs, this%vel1_vec, this%LA_vel, 'add', opt_include_ghost=.FALSE.)

        !===RHS stored in this%vel1_vec
        !===End Combine parabolic fluxes with IRK coefficients

        !====Construct matrix
        DO k = 1, k_dim
            lhs_mass((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*rho_l
        END DO
        CALL array_to_petsc_vec(lhs_mass, this%vel2_vec, this%matrices%LA_vel, 'insert')
        !===NOTICE: rho*ML is stored in this%vel2_vec
        
        IF (stage < this%IRK%s + 1) THEN
            !=== VB 14/07/2026: rhs BC at time_stage
            time_stage = this%time+this%IRK%C(stage)*this%dt
            DO k = 1, k_dim
                CALL dirichlet_rhs(this%matrices%LA_vel%loc_to_glob(k, this%bc%vel(k)%jsd)-1, &
                                this%bc%vit_anal(k, time_stage, this%mesh%rr(:,this%bc%vel(k)%jsd)), this%vel1_vec)
            END DO

            !=== LHS matrix construction + BC
            CALL MatCopy(this%matrices%vel_diff_mat, this%matrices%vel_mat, SAME_NONZERO_PATTERN, ierr)
            CALL MatScale(this%matrices%vel_mat, this%dt*this%IRK%MatRK(stage, stage), ierr)
            CALL MatDiagonalSet(this%matrices%vel_mat, this%vel2_vec, ADD_VALUES, ierr)
            !CALL periodic_matrix_petsc(mesh%per, this%matrices%LA_vel, this%matrices%vel_mat) !<===FIXME: PERIODIC BCs NOT DONE YET
            DO k = 1, k_dim
                CALL Dirichlet_M_parallel(this%matrices%vel_mat, this%LA_vel%loc_to_glob(k,this%bc%vel(k)%jsd))
            END DO
            !=== LHS matrix construction + BC

            !=== Solver linear system
            CALL this%iterative_LA(this%matrices%elasticity_solver_param, &
                                this%matrices%vel_mat, this%matrices%vel_ksp, this%matrices%precond_vel_mat,&
                                this%vel1_vec, this%sol_vel_vec)
        ELSE
            !===Divide by rho*ML stored in this%vel2_vec
            CALL VecPointWiseDivide(this%sol_vel_vec, this%vel1_vec, this%vel2_vec, ierr)

            !=== VB 14/07/2026: explicit impose BC
            time_stage = this%time+this%IRK%C(stage)*this%dt
            DO k = 1, k_dim
                CALL dirichlet_rhs(this%matrices%LA_vel%loc_to_glob(k, this%bc%vel(k)%jsd)-1, &
                                this%bc%vit_anal(k, time_stage, this%mesh%rr(:,this%bc%vel(k)%jsd)), this%sol_vel_vec)
            END DO
            !=== explicit impose BC
        END IF
        !====extract velocity and update momentum
        DO k=1, k_dim
            CALL extract_through_ghost(this%sol_vel_vec, k, k, this%LA_vel, vel_out(:, k), opt_assemble=.FALSE.)
            urk(:,k+1, stage) = rho_l*vel_out(:, k)
        END DO

        !===================================!
        !======== TEMPERATURE UPDATE========!
        !===================================!

        CALL VecZeroEntries(this%temp1_vec, ierr)
        
        !=== (-1)Div(sigma(vel).vel)
        IF (stage==2) THEN
            CALL viscous_production (this, velocity_lm1, this%temp_flux_rk_at_dof(stage-1))
        END IF

        !===Define temperature at stage-1
        scal_temp = ((urk(:,k_dim+2,stage-1)/rho_lm1 - 0.5d0*SUM(velocity_lm1**2,DIM=2)))/this%cv !<===Temp= (E/rho - 1/2 * vel**2)/cv
        CALL array_to_petsc_vec(scal_temp, this%temp2_vec, this%LA_temp, 'insert')
        CALL MatMultAdd(this%matrices%temp_diff_mat, this%temp2_vec, this%temp_flux_rk_at_dof(stage-1), this%temp_flux_rk_at_dof(stage-1), ierr)

        IF (stage < this%IRK%s + 1) THEN
            CALL viscous_production (this, vel_out, this%temp_flux_rk_at_dof(stage))

            !===Combine parabolic fluxes with IRK coefficients (notice (-1)*dt*IRK%MatRK)
            !=== Notice: sum goes from 1 to stage, not stage - 1!! (first solve for vel, then for temp)
            CALL VecMAXPY(this%temp1_vec, stage, -this%dt*this%IRK%MatRK(stage,1:stage), & 
                                                this%temp_flux_rk_at_dof(1:stage), ierr) !<=== x1 receives sum of IRK fluxes
        ELSE
            CALL VecMAXPY(this%temp1_vec, stage-1, -this%dt*this%IRK%MatRK(stage,1:stage-1), & 
                                                this%temp_flux_rk_at_dof(1:stage-1), ierr) !<=== x1 receives sum of IRK fluxes
        END IF

        scal_temp = urk(:,k_dim+2,stage) - 0.5d0*rho_l*sum(vel_out**2,dim=2)
        !=== Multiply by lumped mass
        rhs_temp = this%matrices%scal_lumped_mass*scal_temp
        CALL array_to_petsc_vec(rhs_temp, this%temp1_vec, this%LA_temp, 'add', opt_include_ghost=.FALSE.)
        !===RHS stored in this%temp1_vec
        !===End Combine parabolic fluxes with IRK coefficients

        !===Construct Matrix  
        scal_temp = this%matrices%scal_lumped_mass*rho_l*this%cv
        CALL array_to_petsc_vec(scal_temp, this%temp2_vec, this%matrices%LA_temp, 'insert')
        !===NOTICE: rho*cv*ML is stored in this%vel2_vec
        IF (stage < this%IRK%s + 1) THEN
            !=== rhs BC
            CALL dirichlet_rhs(this%matrices%LA_temp%loc_to_glob(1, this%bc%temp%jsd)-1, &
                            this%bc%temp_anal(this%time, this%mesh%rr(:,this%bc%temp%jsd)), this%temp1_vec)
            !=== rhs BC

            !=== LHS matrix construction + BC
            CALL MatCopy(this%matrices%temp_diff_mat, this%matrices%temp_mat, SAME_NONZERO_PATTERN, ierr)
            CALL MatScale(this%matrices%temp_mat, this%dt*this%IRK%MatRK(stage, stage), ierr)
            CALL MatDiagonalSet(this%matrices%temp_mat, this%temp2_vec, ADD_VALUES, ierr)
            !CALL periodic_matrix_petsc(mesh%per, this%matrices%LA_temp, this%matrices%temp_mat) !<===FIXME: PERIODIC BCs NOT DONE YET
            CALL Dirichlet_M_parallel(this%matrices%temp_mat, this%LA_temp%loc_to_glob(1,this%bc%temp%jsd))
            !=== LHS matrix construction + BC

            CALL this%iterative_LA(this%matrices%temperature_solver_param, &
                                this%matrices%temp_mat, this%matrices%temp_ksp, this%matrices%precond_temp_mat,&
                                this%temp1_vec, this%sol_temp_vec)
        ELSE
            CALL VecPointWiseDivide(this%sol_temp_vec, this%temp1_vec, this%temp2_vec, ierr)
        END IF
        !====extract velocity and update momentum
        CALL extract_through_ghost(this%sol_temp_vec, 1, 1, this%LA_temp, temp_out, opt_assemble=.FALSE.)
        urk(:,k_dim+2, stage) = rho_l*this%cv*temp_out + 0.5d0*rho_l*sum(vel_out**2,dim=2)
    END SUBROUTINE one_step_IRK_LUMPED

    SUBROUTINE one_step_IRK_FULL(this,stage,urk)
        USE my_util, ONLY : error_petsc, to_str, user_time
        USE petsc_tools
        USE space_dim
        USE dir_nodes_petsc
        USE st_matrix
        USE fem_rhs, ONLY: qs_00_block, qs_00
        USE fem_M,   ONLY: qs_var_mass_block_M
        USE sub_plot
        USE fem_tn, ONLY: ns_l1_par
        IMPLICIT NONE
        CLASS(stokes_parabolic_type)                                                :: this
        REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%IRK%s+1), TARGET   :: urk
        REAL(KIND = 8), DIMENSION(:), POINTER                          :: rho_lm1, rho_l
        REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim)                 :: velocity_lm1, vel_out, ff_vel
        REAL(KIND = 8), DIMENSION(k_dim*this%mesh%np)                  :: rhs
        REAL(KIND = 8), DIMENSION(this%mesh%np)                        :: rhs_temp, scal_temp, temp_out
        REAL(KIND = 8), DIMENSION(k_dim*this%mesh%np)                  :: lhs_mass
        INTEGER, INTENT(IN) :: stage
        INTEGER :: k, l, np, ierr!, stage_prime
        REAL(KIND = 8) :: local_tps_solver_ref, local_tps_solver_one_iter, tps_solver_one_iter, tps, time_stage, error, norm, max_vel_loc, max_vel
        !=== Build pointers
        rho_lm1 => urk(:, 1, stage-1)
        rho_l   => urk(:, 1, stage)
        !=== Build pointers
        np = this%mesh%np
        ! stage_prime = this%IRK%lp_of_l(stage) 

        !=== Init rhs vector for velocity problem
        CALL VecZeroEntries(this%vel1_vec, ierr)
        !============================================================!
        !======== FORCING CONTRIBUTION AT STAGES stage-1 & stage ====!
        !============================================================!

        !=== Forcing ===!
        IF (ASSOCIATED(this%forcing)) THEN
            IF (stage == 2) THEN
                time_stage = this%time+this%IRK%C(stage-1)*this%dt

                !=== Consistent version
                DO k=1, k_dim
                    ff_vel(:, k) = this%forcing(k, this%mesh%rr, time_stage)
                END DO
                CALL qs_00_block (this%mesh, this%LA_vel, -ff_vel, this%forcing_rk_at_dof(stage-1))

                !=== Lumped version
                ! DO k=1, k_dim
                !     rhs((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*this%forcing(k, this%mesh%rr, time_stage)
                ! END DO
                ! CALL array_to_petsc_vec(-rhs, this%forcing_rk_at_dof(stage-1), this%LA_vel, 'insert')
            END IF
            IF (stage < this%IRK%s + 1) THEN
                time_stage = this%time+this%IRK%C(stage)*this%dt

                !=== Consistent version
                DO k=1, k_dim
                    ff_vel(:, k) = this%forcing(k, this%mesh%rr, time_stage)
                END DO
                CALL qs_00_block (this%mesh, this%LA_vel, -ff_vel, this%forcing_rk_at_dof(stage))
                
                !=== Lumped version
                ! DO k=1, k_dim
                !     rhs((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*this%forcing(k, this%mesh%rr, time_stage)
                ! END DO
                ! CALL array_to_petsc_vec(-rhs, this%forcing_rk_at_dof(stage), this%LA_vel, 'insert')

                CALL VecMAXPY(this%vel1_vec, stage, -this%dt*this%IRK%MatRK(stage,1:stage), this%forcing_rk_at_dof(1:stage), ierr) !<=== x1 receives sum of ERK forcing
            ELSE
                CALL VecMAXPY(this%vel1_vec, stage-1, -this%dt*this%IRK%MatRK(stage,1:stage-1), this%forcing_rk_at_dof(1:stage-1), ierr) !<=== x1 receives sum of ERK forcing
            END IF
        END IF
        !=== Forcing ===!

        !========================================================!
        !======== VELOCITY VISCOUS DISSIPATION at stage-1========!
        !========================================================!

        !===Define velocity at stage-1
        DO k = 1, k_dim
            velocity_lm1(:,k) = urk(:,k+1,stage-1)/rho_lm1
            rhs((k-1)*np+1:k*np) = velocity_lm1(:,k)
        END DO
        CALL array_to_petsc_vec(rhs, this%vel2_vec, this%LA_vel, 'insert')
        !=== (-1)Div(sigma(vel))
        CALL MatMult(this%matrices%vel_diff_mat, this%vel2_vec, this%vel_flux_rk_at_dof(stage-1), ierr) 

        !================================!
        !======== VELOCITY UPDATE========!
        !================================!
        !===Combine parabolic fluxes with IRK coefficients (notice (-1)*dt*IRK%MatRK)

        CALL VecMAXPY(this%vel1_vec, stage-1, -this%dt*this%IRK%MatRK(stage,1:stage-1), this%vel_flux_rk_at_dof(1:stage-1), ierr) !<=== x1 receives sum of IRK fluxes

        ! DO k = 1, k_dim
        !     !=== WARNING HERE: stage_prime involved in explicit step or must be set manually if no prior explicit step 
        !     rhs((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*urk(:,k+1,stage)
        ! END DO
        ! CALL array_to_petsc_vec(rhs, this%vel1_vec, this%LA_vel, 'add', opt_include_ghost=.FALSE.)
        CALL qs_00_block (this%mesh, this%LA_vel, urk(:,2:k_dim+1,stage), this%vel2_vec)
        CALL VecAXPY(this%vel1_vec, 1.d0, this%vel2_vec, ierr)


        !===RHS stored in this%vel1_vec
        !===End Combine parabolic fluxes with IRK coefficients

        !====Construct matrix
        ! DO k = 1, k_dim
            ! lhs_mass((k-1)*np+1:k*np) = this%matrices%scal_lumped_mass*rho_l
        ! END DO
        ! CALL array_to_petsc_vec(lhs_mass, this%vel2_vec, this%matrices%LA_vel, 'insert')
        CALL qs_var_mass_block_M (this%matrices%vel_mass_mat, this%mesh, this%LA_vel, 1.d0, rho_l)
        !===NOTICE: rho*ML is stored in this%vel2_vec

        !=== VB 14/07/2026: rhs BC at time_stage
        time_stage = this%time+this%IRK%C(stage)*this%dt
        DO k = 1, k_dim
            CALL dirichlet_rhs(this%matrices%LA_vel%loc_to_glob(k, this%bc%vel(k)%jsd)-1, &
                            this%bc%vit_anal(k, time_stage, this%mesh%rr(:,this%bc%vel(k)%jsd)), this%vel1_vec)
        END DO       

        IF (stage < this%IRK%s + 1) THEN
            !=== LHS matrix construction + BC
            CALL MatCopy(this%matrices%vel_diff_mat, this%matrices%vel_mat, SAME_NONZERO_PATTERN, ierr)
            CALL MatScale(this%matrices%vel_mat, this%dt*this%IRK%MatRK(stage, stage), ierr)
            CALL MatAXPY(this%matrices%vel_mat, 1.d0, this%matrices%vel_mass_mat, SAME_NONZERO_PATTERN, ierr)
        ELSE
            ! CALL MatZeroEntries(this%matrices%vel_mat, ierr)
            CALL MatCopy(this%matrices%vel_mass_mat, this%matrices%vel_mat, SAME_NONZERO_PATTERN, ierr)
        ENDIF
        ! CALL MatDiagonalSet(this%matrices%vel_mat, this%vel2_vec, ADD_VALUES, ierr)
        
            !CALL periodic_matrix_petsc(mesh%per, this%matrices%LA_vel, this%matrices%vel_mat) !<===FIXME: PERIODIC BCs NOT DONE YET
        
        DO k = 1, k_dim
            CALL Dirichlet_M_parallel(this%matrices%vel_mat, this%LA_vel%loc_to_glob(k,this%bc%vel(k)%jsd))
        END DO
        !=== LHS matrix construction + BC

        !=== Solver linear system
        CALL this%iterative_LA(this%matrices%elasticity_solver_param, &
                            this%matrices%vel_mat, this%matrices%vel_ksp, this%matrices%precond_vel_mat,&
                            this%vel1_vec, this%sol_vel_vec)


        !====extract velocity and update momentum
        DO k=1, k_dim
            CALL extract_through_ghost(this%sol_vel_vec, k, k, this%LA_vel, vel_out(:, k), opt_assemble=.FALSE.)
            urk(:,k+1, stage) = rho_l*vel_out(:, k)
        END DO

        !===================================!
        !======== TEMPERATURE UPDATE========!
        !===================================!

        !=== (-1)Div(sigma(vel).vel)
        IF (stage==2) THEN
            CALL viscous_production (this, velocity_lm1, this%temp_flux_rk_at_dof(stage-1))
        END IF
        !===Define temperature at stage-1
        scal_temp = ((urk(:,k_dim+2,stage-1)/rho_lm1 - 0.5d0*SUM(velocity_lm1**2,DIM=2)))/this%cv !<===Temp= (E/rho - 1/2 * vel**2)/cv
        CALL array_to_petsc_vec(scal_temp, this%temp1_vec, this%LA_temp, 'insert')
        CALL MatMult(this%matrices%temp_diff_mat, this%temp1_vec, this%temp2_vec, ierr)
        CALL VecAXPY(this%temp_flux_rk_at_dof(stage-1), 1.d0, this%temp2_vec, ierr)

        CALL VecZeroEntries(this%temp1_vec, ierr)
        IF (stage < this%IRK%s + 1) THEN
            CALL viscous_production (this, vel_out, this%temp_flux_rk_at_dof(stage))

            !===Combine parabolic fluxes with IRK coefficients (notice (-1)*dt*IRK%MatRK)
            !=== Notice: sum goes from 1 to stage, not stage - 1!! (first solve for vel, then for temp)
            CALL VecMAXPY(this%temp1_vec, stage, -this%dt*this%IRK%MatRK(stage,1:stage), & 
                                                this%temp_flux_rk_at_dof(1:stage), ierr) !<=== x1 receives sum of IRK fluxes
        ELSE
            CALL VecMAXPY(this%temp1_vec, stage-1, -this%dt*this%IRK%MatRK(stage,1:stage-1), & 
                                                this%temp_flux_rk_at_dof(1:stage-1), ierr) !<=== x1 receives sum of IRK fluxes
        END IF

        scal_temp = urk(:,k_dim+2,stage) - 0.5d0*rho_l*sum(vel_out**2,dim=2)
        !=== Multiply by lumped mass
        ! CALL qs_00(this%mesh, this%LA_temp, scal_temp, this%temp2_vec)
        ! CALL VecAXPY(this%temp1_vec, 1.d0, this%temp2_vec, ierr)
        rhs_temp = this%matrices%scal_lumped_mass*scal_temp
        CALL array_to_petsc_vec(rhs_temp, this%temp1_vec, this%LA_temp, 'add', opt_include_ghost=.FALSE.)
        !===RHS stored in this%temp1_vec
        !===End Combine parabolic fluxes with IRK coefficients

        !===Construct Matrix  
        ! CALL qs_var_mass_block_M (this%matrices%temp_mass_mat, this%mesh, this%LA_temp, this%cv, rho_l)
        scal_temp = this%matrices%scal_lumped_mass*rho_l*this%cv
        CALL array_to_petsc_vec(scal_temp, this%temp2_vec, this%matrices%LA_temp, 'insert')
        !===NOTICE: rho*cv*ML is stored in this%vel2_vec
        !=== rhs BC
        CALL dirichlet_rhs(this%matrices%LA_temp%loc_to_glob(1, this%bc%temp%jsd)-1, &
                        this%bc%temp_anal(this%time, this%mesh%rr(:,this%bc%temp%jsd)), this%temp1_vec)
        !=== rhs BC

        IF (stage < this%IRK%s + 1) THEN

            !=== LHS matrix construction + BC
            CALL MatCopy(this%matrices%temp_diff_mat, this%matrices%temp_mat, SAME_NONZERO_PATTERN, ierr)
            CALL MatScale(this%matrices%temp_mat, this%dt*this%IRK%MatRK(stage, stage), ierr)
            ! CALL MatAXPY(this%matrices%temp_mat, 1.d0, this%matrices%temp_mass_mat, SAME_NONZERO_PATTERN, ierr)
            CALL MatDiagonalSet(this%matrices%temp_mat, this%temp2_vec, ADD_VALUES, ierr)
            !CALL periodic_matrix_petsc(mesh%per, this%matrices%LA_temp, this%matrices%temp_mat) !<===FIXME: PERIODIC BCs NOT DONE YET
        ELSE
            CALL MatZeroEntries(this%matrices%temp_mat, ierr)
            CALL MatDiagonalSet(this%matrices%temp_mat, this%temp2_vec, ADD_VALUES, ierr)
            ! CALL MatCopy(this%matrices%temp_mass_mat, this%matrices%temp_mat, SAME_NONZERO_PATTERN, ierr)
        ENDIF

        CALL Dirichlet_M_parallel(this%matrices%temp_mat, this%LA_temp%loc_to_glob(1,this%bc%temp%jsd))
        !=== LHS matrix construction + BC

        CALL this%iterative_LA(this%matrices%temperature_solver_param, &
                    this%matrices%temp_mat, this%matrices%temp_ksp, this%matrices%precond_temp_mat,&
                    this%temp1_vec, this%sol_temp_vec)
        
        !====extract velocity and update momentum
        CALL extract_through_ghost(this%sol_temp_vec, 1, 1, this%LA_temp, temp_out, opt_assemble=.FALSE.)
        urk(:,k_dim+2, stage) = rho_l*this%cv*temp_out + 0.5d0*rho_l*sum(vel_out**2,dim=2)
    END SUBROUTINE one_step_IRK_FULL

    SUBROUTINE iterative_LA(this, solver_param, matrix, matrix_ksp, matrix_precond, vec_rhs, vec_sol)
        USE petscksp
        USE my_util
        USE solver_petsc
        IMPLICIT NONE
        CLASS(stokes_parabolic_type) :: this
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
            CALL init_solver(this%communicator, solver_param, matrix_ksp, &
            matrix, opt_mat_pre=matrix_precond)!, opt_re_init=.TRUE.)  
        ELSE
            CALL KSPSetInitialGuessNonzero(matrix_ksp, PETSC_TRUE, ierr)
            CALL KSPSetReusePreconditioner(matrix_ksp, PETSC_TRUE, ierr)
        END IF

        tps = user_time()
        CALL solver(matrix_ksp, vec_rhs, vec_sol, reinit = .FALSE., verbose = solver_param%if_verbose)
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

    SUBROUTINE viscous_production (this, vel, vect)
        USE space_dim
        IMPLICIT NONE
        CLASS(stokes_parabolic_type) :: this
        REAL(KIND=8), DIMENSION(this%mesh%np,k_dim), INTENT(IN) :: vel
        REAL(KIND=8), DIMENSION(k_dim,this%mesh%gauss%n_w) :: dw_loc
        REAL(KIND=8), DIMENSION(this%mesh%gauss%n_w)       :: rhs_loc
        REAL(KIND=8), DIMENSION(this%mesh%np)              :: rhs
        REAL(KIND=8), DIMENSION(k_dim,k_dim)               :: gradl, gradTl
        REAL(KIND=8), DIMENSION(k_dim)                     :: vell, xl, tensorkl
        INTEGER,      DIMENSION(this%mesh%gauss%n_w)       :: j_loc, idxm
        INTEGER,      DIMENSION(this%mesh%np)       :: idx
        REAL(KIND=8)     :: bulk_visc, divl
        INTEGER          ::  m, l, ni, k, kp, i, iglob
        TYPE(tVec)              :: vect
        INTEGER   :: ierr
        CALL VecZeroEntries(vect, ierr)

        DO i = 1, this%mesh%np
            idx(i) = this%LA_temp%loc_to_glob(1, i)-1
        END DO

        rhs = 0.d0
        bulk_visc = this%lambda_viscosity- 2*this%mu_viscosity/3
        DO m = 1, this%mesh%me
            j_loc = this%mesh%jj(:,m)

            rhs_loc = 0.d0
            DO l = 1, this%mesh%gauss%l_G
                DO k = 1, k_dim
                    vell(k) = SUM(vel(j_loc,k)*this%mesh%gauss%ww(:,l)) !<==v_k
                END DO
                dw_loc = this%mesh%gauss%dw(:,:,l,m)
                divl = 0.d0
                DO k = 1, k_dim
                    DO kp = 1, k_dim
                        gradl(k,kp) = SUM(vel(j_loc,k)*dw_loc(kp,:)) !<==d(v_k)/d(x_kp)
                        gradTl(kp,k) = gradl(k,kp)
                    END DO
                    divl = divl + gradl(k,k)
                END DO
                xl = 0.d0
                DO k = 1, k_dim
                    DO kp = 1, k_dim
                        tensorkl(kp) = this%mu_viscosity*(gradl(k,kp) + gradTl(k,kp))
                    END DO
                    tensorkl(k) = tensorkl(k) + bulk_visc*divl !<==2*mu(Grad+GradT)/2 + (lambda-2*mu/3)Div
                    xl(k) = sum(tensorkl*vell) !<== Tensor.vel
                END DO
                xl = xl*this%mesh%gauss%rj(l,m) 
                DO ni = 1, this%mesh%gauss%n_w
                    rhs_loc(ni) = rhs_loc(ni) + SUM(xl*dw_loc(:,ni)) !<== (Tensor.vel).Grad(phi_i); accumulation over Gauss points
                END DO
            END DO
            rhs(j_loc) =  rhs(j_loc) + rhs_loc
        ENDDO
        CALL VecSetValues(vect, this%mesh%np, idx, rhs, ADD_VALUES, ierr)
        CALL VecAssemblyBegin(vect, ierr)
        CALL VecAssemblyEnd(vect, ierr)
        !===Try to use on long vector of size np to do the VecAssembly

    END SUBROUTINE viscous_production

END MODULE stokes_parabolic_module