MODULE start_setup_MODULE
#include "petsc/finclude/petsc.h"
  USE petsc
  USE def_type_mesh
  USE navier_stokes_module
  USE limiting_functionals_euler_module, ONLY: psi_euler_rho_min, psi_euler_rho_max,&
                               zero_of_psi_euler_rho_max, zero_of_psi_euler_rho_min
  USE euler_cell_limiting_engine_parallel_module
  USE petsc_csr_LA_module, ONLY: petsc_csr_LA

  USE read_inputs_module
  USE setup,           ONLY: init_state_functions
  TYPE argument_setup_data_type
     CHARACTER(LEN=rec_length) :: if_restart         = '=== Restart (true/false) ==='
     CHARACTER(LEN=rec_length) :: checkpointing_freq = '=== Checkpointing frequency ==='
     CHARACTER(LEN=rec_length) :: verbose_freq       = '=== Frequency for run verbose ==='
     CHARACTER(LEN=rec_length) :: final_time         = '=== Final time ==='
     CHARACTER(LEN=rec_length) :: max_it             = '=== Maximum number of timesteps ==='
     CHARACTER(LEN=rec_length) :: imex_sv            = '=== IMEX rk ? ==='
     CHARACTER(LEN=rec_length) :: if_analytical_ref  = '=== Do we compare with analytical reference? (true/false) ==='
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
     INTEGER        :: imex_sv             = -31
   CONTAINS
     PROCEDURE, PUBLIC :: read => read_setup_data
     PROCEDURE, PUBLIC :: init => init_setup_data
  END TYPE setup_data_type

  TYPE(mesh_type),                   PUBLIC :: mesh
  TYPE(navier_stokes_type),          PUBLIC :: navier_stokes
  TYPE(setup_data_type),             PUBLIC :: setup_data
  TYPE(limiting_functional_type), DIMENSION(:), ALLOCATABLE, PRIVATE :: limiting_functionals_euler
  TYPE(periodic_type), DIMENSION(1), PUBLIC :: per
  MPI_Comm :: communicator
  PUBLIC :: start_setup
  PRIVATE

CONTAINS

  SUBROUTINE start_setup
    USE construct_mesh,     ONLY: get_mesh
    USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
    USE setup
    IMPLICIT NONE
    PetscErrorCode :: ierr
    REAL(KIND = 8), DIMENSION(2) :: times = (/0.d0,1.d0/)
    CHARACTER(100) :: name = 'NS'
    INTEGER :: rank

    !===Start PETSC and MPI (mandatory)

    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    communicator = PETSC_COMM_WORLD
    CALL MPI_Comm_rank(communicator, rank, ierr)

    !===Clean data once
    CALL clean_data_once

    !===Construct mesh
    CALL get_mesh(communicator, mesh)

    !===Read
    CALL setup_data%init
    times(2) = setup_data%final_time

    !=== Define Euler limiting bounds (should we put this in PROBLEM_SOURCES instead?)
    ALLOCATE(limiting_functionals_euler(2))
    limiting_functionals_euler(1)%psi => psi_euler_rho_min
    limiting_functionals_euler(1)%zero_of_psi => zero_of_psi_euler_rho_min
    limiting_functionals_euler(1)%name = 'rho min'
    limiting_functionals_euler(2)%psi => psi_euler_rho_max
    limiting_functionals_euler(2)%zero_of_psi => zero_of_psi_euler_rho_max
    limiting_functionals_euler(2)%name = 'minus rho max'
    !=== Define Euler limiting bounds

    !===Start Navier-Stokes
    navier_stokes%imex_sv = setup_data%imex_sv
    CALL navier_stokes%init(communicator, name, mesh, limiting_functionals_euler)
    CALL init_state_functions(navier_stokes%euler, navier_stokes%stokes, navier_stokes%thermal_diffusivity)
    ! CALL init_state_functions(navier_stokes%stokes)
    CALL navier_stokes%set_times(times)
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

    !===IMEX rk
    CALL read_data(argument_data%imex_sv, this%imex_sv)

    !===Analytical reference
    CALL read_data(argument_data%if_analytical_ref, this%if_analytical_ref)

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
