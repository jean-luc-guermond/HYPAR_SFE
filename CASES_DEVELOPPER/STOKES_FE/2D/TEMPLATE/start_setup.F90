#include "petsc/finclude/petsc.h"
MODULE start_setup_MODULE

    USE petsc
    USE def_type_mesh
    USE read_inputs_module
    USE setup
    USE solver_data_module
    USE petsc_csr_LA_module, ONLY: petsc_csr_LA
    USE stokes_parabolic_module
    USE stokes_bc_arrays

    TYPE argument_setup_data_type
        CHARACTER(LEN=rec_length) :: if_analytical_ref    = '=== Do we compare with analytical reference? (true/false) ==='
        CHARACTER(LEN=rec_length) :: test_case            = '=== Which test case? (1:6) ==='
        CHARACTER(LEN=rec_length) :: mu_viscosity         = '=== Value of mu viscosity ==='
        CHARACTER(LEN=rec_length) :: lambda_viscosity     = '=== Value of lambda viscosity ==='
        CHARACTER(LEN=rec_length) :: thermal_diffusivity  = '=== Value of thermal diffusivity ==='
        CHARACTER(LEN=rec_length) :: cv                   = '=== Value of thermal capacity at constant volume ==='
        CHARACTER(LEN=rec_length) :: verbose_freq         = '=== Frequency for run verbose ==='
        CHARACTER(LEN=rec_length) :: dt                   = '=== timestep ==='
        CHARACTER(LEN=rec_length) :: max_it               = '=== Maximum number of timesteps ==='
        CHARACTER(LEN=rec_length) :: irk_sv               = '=== IRK ? ==='

    END TYPE argument_setup_data_type

    TYPE setup_data_type
        LOGICAL        :: if_analytical_ref   = .FALSE.
        INTEGER        :: test_case           = 4
        INTEGER        :: syst_size
        REAL(KIND = 8)                :: thermal_diffusivity = 0.d0
        REAL(KIND = 8)                :: mu_viscosity        = 0.d0
        REAL(KIND = 8)                :: lambda_viscosity    = 0.d0
        REAL(KIND = 8)                :: cv                  = 0.d0
        REAL(KIND = 8)                :: CFL                 = 0.5d0
        REAL(KIND = 8)                :: t_final             = 1.d0
        REAL(KIND = 8)                :: dt                  = 1.d-3
        INTEGER                       :: verbose_freq        = 100
        INTEGER                       :: max_it              = 10
        INTEGER                       :: irk_sv              = -21
    CONTAINS
        PROCEDURE, PUBLIC :: read => read_setup_data
        PROCEDURE, PUBLIC :: init => init_setup_data
    END TYPE setup_data_type

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: un
    TYPE(mesh_type),                   PUBLIC :: mesh
    TYPE(petsc_csr_LA),                PUBLIC :: LA_temp, LA_vel
    TYPE(setup_data_type),             PUBLIC :: setup_data
    TYPE(stokes_parabolic_type),       PUBLIC :: stokes

    TYPE(solver_data_type), PUBLIC :: elasticity_solver_param
    MPI_Comm, PUBLIC :: communicator
    PUBLIC :: start_setup
    PRIVATE

CONTAINS

  SUBROUTINE start_setup
    USE construct_mesh,     ONLY: get_mesh
    USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
    USE setup
    USE space_dim
    USE options_module
    IMPLICIT NONE
    CHARACTER(100) :: name = 'stokes'
    INTEGER        :: rank, k
    REAL(KIND=8), DIMENSION(2) :: times
    INTEGER :: ierr

    !===Start PETSC and MPI (mandatory)

    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    communicator = PETSC_COMM_WORLD
    CALL MPI_Comm_rank(communicator, rank, ierr)

    !===Read executable arguments
    CALL read_all_arguments
  
    !===Clean data once
    CALL clean_data_once

    !===Construct mesh
    CALL get_mesh(communicator, mesh)
    
    !===Construct LA
    CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA_temp)
    CALL st_aij_csr_glob_block_with_extra_layer(communicator, k_dim, mesh, LA_vel)
    
    CALL elasticity_solver_param%init('elasticity')
    
    !===Read
    CALL setup_data%init
    stokes%thermal_diffusivity = setup_data%thermal_diffusivity
    stokes%mu_viscosity        = setup_data%mu_viscosity
    stokes%lambda_viscosity    = setup_data%lambda_viscosity
    stokes%cv                  = setup_data%cv
    stokes%irk_sv = setup_data%irk_sv
    stokes%syst_dim = k_dim+2

    !=== Create Stokes object
    CALL stokes%init(communicator, name, mesh)
    CALL init_state_functions(stokes, setup_data%test_case)

    !=== Init times
    times(1) = 0.d0
    times(2) = setup_data%dt*setup_data%max_it
    CALL stokes%set_times(times)
    stokes%dt = setup_data%dt

    !=== Initial conditions
    ALLOCATE(un(mesh%np, k_dim+2))
    un(:, 1) = stokes%bc%rho_imposed(stokes%time, mesh%rr)
    DO k=1, k_dim
      un(:, k+1) = stokes%bc%vit_anal(k, stokes%time, mesh%rr)*stokes%bc%rho_imposed(stokes%time, mesh%rr)
    END DO
    un(:, k_dim+2) = Etot(mesh%rr, stokes%time) !<=== FIXME

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

    !================
    !=== MANDATORY Reading all data file
    !================
    CALL read_data_init_list(section_name)

    !================
    !=== We now find the relevant information for this setup
    !================

    !===Analytical reference
    CALL read_data(argument_data%if_analytical_ref, this%if_analytical_ref)

    !===test_case
    CALL read_data(argument_data%test_case, this%test_case)

    !===mu
    CALL read_data(argument_data%mu_viscosity, this%mu_viscosity)

    !===lambda
    CALL read_data(argument_data%lambda_viscosity, this%lambda_viscosity)

    !===thermal_diffusivity
    CALL read_data(argument_data%thermal_diffusivity, this%thermal_diffusivity)
    
    !===cv
    CALL read_data(argument_data%cv, this%cv)
    
    !===timestep
    CALL read_data(argument_data%dt, this%dt)
    
    !===Maximum number of iterations
    CALL read_data(argument_data%max_it, this%max_it)

    !===verbose_freq
    CALL read_data(argument_data%verbose_freq, this%verbose_freq)

    !===irk_sv
    CALL read_data(argument_data%irk_sv, this%irk_sv)

    !================
    !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
    !================
    CALL finalize_rewrite_data

  END SUBROUTINE read_setup_data

END MODULE start_setup_MODULE
