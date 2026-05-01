MODULE start_setup_MODULE
#include "petsc/finclude/petsc.h"
  USE petsc
  USE def_type_mesh
  USE euler_type_module
  USE read_inputs_module, ONLY : clean_data_once
  MPI_Comm        :: communicator


  TYPE argument_setup_data_type
     CHARACTER(LEN=rec_length) :: if_restart         = '=== Restart (true/false) ==='
     CHARACTER(LEN=rec_length) :: checkpointing_freq = '=== Checkpointing frequency ==='
     CHARACTER(LEN=rec_length) :: final_time         = '=== Final time ==='
     CHARACTER(LEN=rec_length) :: if_analytical_ref  = '=== Do we compare with analytical reference? (true/false) ==='
  END TYPE argument_setup_data_type

  TYPE setup_data_type
     LOGICAL        :: if_regression_test  = .FALSE.
     LOGICAL        :: if_restart          = .FALSE.
     REAL(KIND = 8) :: checkpointing_freq  = 1.d20
     REAL(KIND = 8) :: final_time          = 0.1d0
     LOGICAL        :: if_analytical_ref   = .FALSE.
     INTEGER        :: syst_size
   CONTAINS
     PROCEDURE, PUBLIC :: read => read_setup_data
     PROCEDURE, PUBLIC :: init => init_setup_data
  END TYPE setup_data_type

  TYPE(mesh_type),                   PUBLIC :: mesh
  TYPE(petsc_csr_LA),               PRIVATE :: LA
  TYPE(euler_type),                  PUBLIC :: euler
  TYPE(setup_data_type),             PUBLIC :: setup_data
  TYPE(periodic_type), DIMENSION(1), PUBLIC :: per
  PUBLIC :: start_setup
  PRIVATE

CONTAINS

  SUBROUTINE start_setup
    use periodic_data_module
    USE construct_mesh
    USE st_matrix
    USE setup
    IMPLICIT NONE
    PetscErrorCode :: ierr
    REAL(KIND = 8), DIMENSION(2) :: times = (/0.d0,1.d0/)
    CHARACTER(100) :: name = 'Euler 1'
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

    !===Start Euler
    times(2) = setup_data%final_time
    CALL euler%init(communicator, name, mesh, LA, pressure, impose_bc_euler, times)
    ! CALL euler%init(communicator, name, mesh, LA, mesh%per, pressure, impose_bc_euler, times)

    !===Read data setup
  END SUBROUTINE start_setup

  SUBROUTINE init_setup_data(this)
    CLASS(setup_data_type), INTENT(INOUT) :: this
    CALL this%read
  END SUBROUTINE init_setup_data

  SUBROUTINE read_setup_data(this)
    USE read_inputs_module
    IMPLICIT NONE

    CHARACTER(LEN=rec_length) :: section_name='SETUP PARAMETERS'

    CLASS(setup_data_type)             :: this
    TYPE(argument_setup_data_type)     :: argument_data

    CHARACTER(LEN=rec_length)     :: string

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

    !===Final time
    CALL read_data(argument_data%final_time, this%final_time)

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
