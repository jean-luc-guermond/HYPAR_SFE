MODULE start_setup_MODULE


    USE petsc
    USE def_type_mesh
    USE read_inputs_module
    USE setup
    USE stokes_parabolic_matrices_module
    USE solver_data_module
    USE dirichlet_type_module

    TYPE argument_setup_data_type
        CHARACTER(LEN=rec_length) :: if_analytical_ref    = '=== Do we compare with analytical reference? (true/false) ==='
        CHARACTER(LEN=rec_length) :: mu_viscosity         = '=== Value of mu viscosity ==='
        CHARACTER(LEN=rec_length) :: lambda_viscosity     = '=== Value of lambda viscosity ==='
        CHARACTER(LEN=rec_length) :: thermal_diffusivity  = '=== Value of thermal diffusivity ==='
        CHARACTER(LEN=rec_length) :: cv                   = '=== Value of thermal capacity at constant volume ==='
        CHARACTER(LEN=rec_length) :: max_it               = '=== Maximum number of timesteps ==='
    END TYPE argument_setup_data_type

    TYPE setup_data_type
        LOGICAL        :: if_regression_test  = .FALSE.
        LOGICAL        :: if_analytical_ref   = .FALSE.
        INTEGER        :: syst_size
        REAL(KIND = 8)                :: thermal_diffusivity = 0.d0
        REAL(KIND = 8)                :: mu_viscosity = 0.d0
        REAL(KIND = 8)                :: lambda_viscosity =0.d0
        REAL(KIND = 8)                :: cv = 0.d0
        INTEGER                       :: max_it=10
    CONTAINS
        PROCEDURE, PUBLIC :: read => read_setup_data
        PROCEDURE, PUBLIC :: init => init_setup_data
    END TYPE setup_data_type

    TYPE(mesh_type),                   PUBLIC :: mesh
    TYPE(petsc_csr_LA),                PUBLIC :: LA_temp, LA_vel
    TYPE(periodic_type), DIMENSION(1), PUBLIC :: per
    TYPE(dirichlet_bc),                PUBLIC :: dir
    TYPE(setup_data_type),             PUBLIC :: setup_data
    TYPE(stokes_parabolic_matrices_type), PUBLIC :: stokes_matrices

    TYPE(solver_data_type), PUBLIC :: elasticity_solver_param
    INTEGER, PUBLIC :: communicator
    PUBLIC :: start_setup
    PRIVATE

CONTAINS

  SUBROUTINE start_setup
    USE construct_mesh,     ONLY: get_mesh
    USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
    USE setup
    USE space_dim
    IMPLICIT NONE
    CHARACTER(100) :: name
    INTEGER        :: rank
    INTEGER :: ierr

    !===Start PETSC and MPI (mandatory)

    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    communicator = PETSC_COMM_WORLD
    CALL MPI_Comm_rank(communicator, rank, ierr)

    !===Clean data once
    CALL clean_data_once

    !===Construct mesh
    CALL get_mesh(communicator, mesh)
    
    !===Construct dirichlet parallel
    CALL dir%set_parallel(mesh, 'scalar')

    !===Construct LA
    CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA_temp)
    CALL st_aij_csr_glob_block_with_extra_layer(communicator, k_dim, mesh, LA_vel)
    
    CALL elasticity_solver_param%init('elasticity')
    ! CALL elasticity_solver_param%set(elasticity_ksp_par)
    
    !===Read
    CALL setup_data%init
    stokes_matrices%thermal_diffusivity = setup_data%thermal_diffusivity
    stokes_matrices%mu_viscosity        = setup_data%mu_viscosity
    stokes_matrices%lambda_viscosity    = setup_data%lambda_viscosity
    stokes_matrices%cv                  = setup_data%cv

    !=== Create matrices
    CALL stokes_matrices%construct(communicator, mesh, LA_vel, LA_temp)

  END SUBROUTINE start_setup

  SUBROUTINE init_setup_data(this)
    CLASS(setup_data_type), INTENT(INOUT) :: this
    CALL this%read
  END SUBROUTINE init_setup_data

  SUBROUTINE read_setup_data(this)
    IMPLICIT NONE

    CHARACTER(LEN=rec_length) :: section_name='SETUP PARAMETERS'
    CLASS(setup_data_type)             :: this
    TYPE(argument_setup_data_type)     :: argument_data
    CHARACTER(LEN=rec_length)          :: string

    !================
    !=== MANDATORY Reading all data file
    !================
    CALL read_data_init_list(section_name)

    !================
    !=== We now find the relevant information for this setup
    !================

    !===Analytical reference
    CALL read_data(argument_data%if_analytical_ref, this%if_analytical_ref)

    !===Regression test
    CALL getarg(1, string)
    IF (trim(adjustl(string))=='regression') THEN
       this%if_regression_test = .true.
    ELSE
       this%if_regression_test = .false.
    END IF

    !===mu
    CALL read_data(argument_data%mu_viscosity, this%mu_viscosity)

    !===lambda
    CALL read_data(argument_data%lambda_viscosity, this%lambda_viscosity)

    !===lambda
    CALL read_data(argument_data%thermal_diffusivity, this%thermal_diffusivity)
    
    !===cv
    CALL read_data(argument_data%cv, this%cv)

    !===Maximum number of iterations
    CALL read_data(argument_data%max_it, this%max_it)

    !================
    !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
    !================
    CALL finalize_rewrite_data

  END SUBROUTINE read_setup_data


END MODULE start_setup_MODULE
