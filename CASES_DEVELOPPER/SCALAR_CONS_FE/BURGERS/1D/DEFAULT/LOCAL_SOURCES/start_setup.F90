MODULE start_setup_MODULE
#include "petsc/finclude/petsc.h"
  USE petsc
  USE def_type_mesh, ONLY: mesh_type
  USE burgers_type_module, ONLY: burgers_type
  USE petsc_csr_LA_module, ONLY: petsc_csr_LA
  USE read_inputs_module
  USE setup,           ONLY: init_state_functions
  USE burgers_cell_limiting_engine_parallel_module
  USE limiting_functionals_burgers_module, ONLY: psi_scalar_rho_min, zero_of_psi_scalar_rho_min, &
                                                          psi_scalar_rho_max, zero_of_psi_scalar_rho_max
  USE burgers_limiter_cell_elt
  USE burgers_uij_bar_bounds

  TYPE argument_setup_data_type
      CHARACTER(LEN=rec_length) :: if_restart         = '=== Restart (true/false) ==='
      CHARACTER(LEN=rec_length) :: checkpointing_freq = '=== Checkpointing frequency ==='
      CHARACTER(LEN=rec_length) :: verbose_freq       = '=== Frequency for run verbose ==='
      CHARACTER(LEN=rec_length) :: final_time         = '=== Final time ==='
      CHARACTER(LEN=rec_length) :: max_it             = '=== Maximum number of timesteps ==='
      CHARACTER(LEN=rec_length) :: erk_sv             = '=== ERK ? ==='
      CHARACTER(LEN=rec_length) :: if_analytical_ref  = '=== Do we compare with analytical reference? (true/false) ==='
      CHARACTER(LEN=rec_length) :: char_which_init    = '=== Which initial condition? (sine, step) ==='
  END TYPE argument_setup_data_type

  TYPE setup_data_type
     LOGICAL        :: if_restart          = .FALSE.
     REAL(KIND = 8) :: checkpointing_freq  = 1.d20
     INTEGER        :: verbose_freq        = 1000000
     REAL(KIND = 8) :: final_time          = 0.1d0
     INTEGER        :: max_it              = 1000000
     INTEGER        :: erk_sv              = -31
     LOGICAL        :: if_analytical_ref   = .FALSE.
     CHARACTER(LEN=rec_length) :: char_which_init = 'step'
     INTEGER        :: which_init
     INTEGER        :: syst_size
   CONTAINS
     PROCEDURE, PUBLIC :: read => read_setup_data
     PROCEDURE, PUBLIC :: init => init_setup_data
  END TYPE setup_data_type

  INTEGER, PRIVATE, PARAMETER :: INIT_BUMP_RHO=1, INIT_SINE_RHO=2, INIT_STEP_RHO=3
  CHARACTER(LEN=20), DIMENSION(2), PRIVATE, PARAMETER  :: list_init_rho = &
              [CHARACTER(LEN=20) :: 'sine', 'step']

  TYPE(mesh_type),                   PUBLIC :: mesh
  TYPE(burgers_type),                PUBLIC :: burgers
  TYPE(setup_data_type),             PUBLIC :: setup_data
  TYPE(periodic_type), DIMENSION(1), PUBLIC :: per
  TYPE(limiting_functional_type), DIMENSION(:), ALLOCATABLE, PRIVATE :: limiting_functionals_burgers
  MPI_Comm :: communicator
  PUBLIC :: start_setup
  PRIVATE

CONTAINS

  SUBROUTINE start_setup
    USE construct_mesh,     ONLY: get_mesh
    USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
    USE setup
    USE options_module
    IMPLICIT NONE
    PetscErrorCode :: ierr
    REAL(KIND = 8), DIMENSION(2) :: times = (/0.d0,1.d0/)
    CHARACTER(100) :: name = 'burgers'
    INTEGER :: rank

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

    !===Read
    CALL setup_data%init

    !===Start burgers
    times(2) = setup_data%final_time
    CALL burgers%init_burgers(name)

    !=== Define burgers limiting bounds (should we put this in PROBLEM_SOURCES instead?)
    burgers%compute_bounds_uijbar => burgers_compute_bounds_uijbar

    ALLOCATE(limiting_functionals_burgers(2))
    limiting_functionals_burgers(1)%psi => psi_scalar_rho_min
    limiting_functionals_burgers(1)%zero_of_psi => zero_of_psi_scalar_rho_min
    limiting_functionals_burgers(1)%name = "min"
    limiting_functionals_burgers(2)%psi => psi_scalar_rho_max
    limiting_functionals_burgers(2)%zero_of_psi => zero_of_psi_scalar_rho_max
    limiting_functionals_burgers(2)%name = "minus max"

    limiting_functionals_burgers(1)%spe_iterative_cell_limiting_procedure => iterative_cell_limiting_procedure_scalar_rho_min
    limiting_functionals_burgers(2)%spe_iterative_cell_limiting_procedure => iterative_cell_limiting_procedure_scalar_rho_max 
    !=== Define burgers limiting bounds
    
    burgers%erk_sv = setup_data%erk_sv
    CALL burgers%init_hyperbolic(communicator, name, mesh, limiting_functionals_burgers)
    CALL init_state_functions(burgers, setup_data%which_init)

    CALL burgers%set_times(times)


  END SUBROUTINE start_setup

  SUBROUTINE init_setup_data(this)
    CLASS(setup_data_type), INTENT(INOUT) :: this
    CALL this%read
  END SUBROUTINE init_setup_data

  SUBROUTINE read_setup_data(this)
    USE my_util, ONLY: get_tab_idx_char
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

    !===Restart
    CALL read_data(argument_data%if_restart, this%if_restart)

    !===Checkpointing
    CALL read_data(argument_data%checkpointing_freq, this%checkpointing_freq)

    !===Verbose frequency
    CALL read_data(argument_data%verbose_freq, this%verbose_freq)

    !===Final time
    CALL read_data(argument_data%final_time, this%final_time)

    !===Maximum number of iterations
    CALL read_data(argument_data%max_it, this%max_it)

    !===ERK Linear Transport
    CALL read_data(argument_data%erk_sv, this%erk_sv)

    !===Analytical reference
    CALL read_data(argument_data%if_analytical_ref, this%if_analytical_ref)

    !===Which init
    CALL read_data(argument_data%char_which_init, this%char_which_init)
    CALL get_tab_idx_char(this%char_which_init, list_init_rho, this%which_init)

    !================
    !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
    !================
    CALL finalize_rewrite_data

  END SUBROUTINE read_setup_data

END MODULE start_setup_MODULE
