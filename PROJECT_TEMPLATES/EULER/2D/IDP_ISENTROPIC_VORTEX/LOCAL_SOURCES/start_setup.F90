MODULE start_setup_MODULE
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE eos
   USE euler_type_module
   USE character_strings, ONLY : clean_data_once
   MPI_Comm        :: communicator
   
   INTEGER, PARAMETER, PRIVATE :: rec_length = 200

   TYPE argument_setup_data_type
      CHARACTER(LEN=rec_length) :: if_restart         = '=== Restart (true/false) ==='
      CHARACTER(LEN=rec_length) :: checkpointing_freq = '=== Checkpointing frequency ==='
      CHARACTER(LEN=rec_length) :: final_time         = '=== Final time ==='
   END TYPE argument_setup_data_type
   
   TYPE setup_data_type
      LOGICAL        :: if_regression_test  = .FALSE.
      LOGICAL        :: if_restart          = .FALSE. 
      REAL(KIND = 8) :: checkpointing_freq  = 1.d20
      REAL(KIND = 8) :: final_time          = 0.1d0
      INTEGER        :: syst_size
   CONTAINS 
      PROCEDURE, PUBLIC :: read => read_setup_data
      PROCEDURE, PUBLIC :: init => init_setup_data
   END TYPE setup_data_type
   
   TYPE(mesh_type), PUBLIC :: mesh
   TYPE(petsc_csr_LA), PRIVATE :: LA
   TYPE(euler_type), PUBLIC :: euler
   TYPE(setup_data_type), PUBLIC :: setup_data
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
      REAL(KIND = 8) :: init_time = 0.d0
      CHARACTER(100) :: name = 'Euler 1'
      INTEGER :: ier, rank

      !===Start PETSC and MPI (mandatory)
      CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      communicator = PETSC_COMM_WORLD
      CALL MPI_Comm_rank(communicator, rank, ierr)
       
      !===Clean data once
      CALL clean_data_once

      !===Construct mesh
      CALL per(1)%init("global", "PERIODIC BC PARAMETERS")
      CALL get_mesh(communicator, mesh, opt_pers = per)
      CALL per(1)%set(mesh)
      !===Construct LA
      CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA, opt_per = per(1))
      !===Read
      CALL setup_data%init

      !===Start Euler
      !FIXE ME init_time too
      CALL euler%init(communicator, name, mesh, LA, per(1), pressure, impose_bc_euler, init_time)

      !===Read data setup
   END SUBROUTINE start_setup

   SUBROUTINE init_setup_data(this)
      CLASS(setup_data_type), INTENT(INOUT) :: this
      CALL this%read
   END SUBROUTINE init_setup_data

   SUBROUTINE read_setup_data(this)
      USE character_strings
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

      !===Checkpointing
    CALL read_real_data(argument_data%final_time, this%final_time)

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
