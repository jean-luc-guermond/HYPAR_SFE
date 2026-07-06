MODULE start_setup_MODULE
#include "petsc/finclude/petsc.h"
  USE petsc
  USE def_type_mesh
  USE linear_transport_type_module

  USE read_inputs_module
  USE setup,           ONLY: init_state_functions, my_linear_transport
  USE cell_limiting_engine_parallel_module
  TYPE argument_setup_data_type
     CHARACTER(LEN=rec_length) :: if_restart         = '=== Restart (true/false) ==='
     CHARACTER(LEN=rec_length) :: checkpointing_freq = '=== Checkpointing frequency ==='
     CHARACTER(LEN=rec_length) :: verbose_freq       = '=== Frequency for run verbose ==='
     CHARACTER(LEN=rec_length) :: final_time         = '=== Final time ==='
     CHARACTER(LEN=rec_length) :: max_it             = '=== Maximum number of timesteps ==='
     CHARACTER(LEN=rec_length) :: if_analytical_ref  = '=== Do we compare with analytical reference? (true/false) ==='
     CHARACTER(LEN=rec_length) :: erk_sv             = '=== ERK ? ==='
  END TYPE argument_setup_data_type

  TYPE setup_data_type
     LOGICAL        :: if_regression_test  = .FALSE.
     LOGICAL        :: if_restart          = .FALSE.
     REAL(KIND = 8) :: checkpointing_freq  = 1.d20
     INTEGER        :: verbose_freq        = 1000000
     REAL(KIND = 8) :: final_time          = 0.1d0
     INTEGER        :: max_it              = 1000000
     LOGICAL        :: if_analytical_ref   = .FALSE.
     INTEGER        :: syst_size
     INTEGER        :: erk_sv              = -31 
   CONTAINS
     PROCEDURE, PUBLIC :: read => read_setup_data
     PROCEDURE, PUBLIC :: init => init_setup_data
  END TYPE setup_data_type

  TYPE(mesh_type),                   PUBLIC :: mesh
  TYPE(petsc_csr_LA),               PRIVATE :: LA
  TYPE(my_linear_transport),        PUBLIC :: linear_transport
  TYPE(setup_data_type),             PUBLIC :: setup_data
  TYPE(periodic_type), DIMENSION(1), PUBLIC :: per
  TYPE(limiting_functional_type), DIMENSION(:), ALLOCATABLE, PRIVATE :: limiting_functionals_linear_transport
  MPI_Comm :: communicator
  PUBLIC :: start_setup
  PRIVATE

CONTAINS

  SUBROUTINE start_setup
    USE construct_mesh,     ONLY: get_mesh
    USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
    USE setup
    USE limiting_functionals_euler_module, ONLY: psi_rho_min, zero_of_psi_rho_min, psi_rho_max, zero_of_psi_rho_max
    IMPLICIT NONE
    PetscErrorCode :: ierr
    REAL(KIND = 8), DIMENSION(2) :: times = (/0.d0,1.d0/)
    CHARACTER(100) :: name = 'linear_transport 1'
    INTEGER :: rank

    !===Start PETSC and MPI (mandatory)
    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    communicator = PETSC_COMM_WORLD
    CALL MPI_Comm_rank(communicator, rank, ierr)

    !===Clean data once
    CALL clean_data_once

    !===Construct mesh
    CALL get_mesh(communicator, mesh)

    !===Construct LA
    CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA)
    
    !===Read
    CALL setup_data%init

    !===Start linear_transport
    times(2) = setup_data%final_time
    CALL linear_transport%init_linear_transport(name)

    !=== Define linear_transport limiting bounds (should we put this in PROBLEM_SOURCES instead?)
    ! ALLOCATE(limiting_functionals_linear_transport(1))
    ALLOCATE(limiting_functionals_linear_transport(2))
    limiting_functionals_linear_transport(1)%psi => psi_rho_min
    limiting_functionals_linear_transport(1)%zero_of_psi => zero_of_psi_rho_min
    limiting_functionals_linear_transport(1)%name = "min"
    limiting_functionals_linear_transport(2)%psi => psi_rho_max
    limiting_functionals_linear_transport(2)%zero_of_psi => zero_of_psi_rho_max
    limiting_functionals_linear_transport(2)%name = "minus max"
    !=== Define linear_transport limiting bounds
    
    linear_transport%erk_sv = setup_data%erk_sv
    CALL linear_transport%init_hyperbolic(communicator, name, mesh, limiting_functionals_linear_transport)
    CALL init_state_functions(linear_transport)
    CALL linear_transport%set_times(times)

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

    !===Analytical reference
    CALL read_data(argument_data%if_analytical_ref, this%if_analytical_ref)

    !===erk_sv
    CALL read_data(argument_data%erk_sv, this%erk_sv)

    !===Regression test
    CALL getarg(1, string)
    IF (trim(adjustl(string))=='regression') THEN
       this%if_regression_test = .true.
    ELSE
       this%if_regression_test = .false.
    END IF

    !================
    !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
    !================
    CALL finalize_rewrite_data

  END SUBROUTINE read_setup_data

END MODULE start_setup_MODULE
