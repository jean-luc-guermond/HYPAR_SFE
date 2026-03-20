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
      INTEGER :: syst_size
   CONTAINS 
      PROCEDURE, PUBLIC :: read => read_setup_data
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
      CALL clean_data_once(rank)

      !===Construct mesh
      CALL per(1)%read("global")
      CALL get_mesh(communicator, mesh, opt_pers = per)
      CALL per(1)%set(mesh)
      !===Construct LA
      CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA, opt_per = per(1))
      !===Read
      CALL setup_data%read(rank)

      !===Start Euler
      !FIXE ME init_time too
      CALL euler%init(communicator, name, mesh, LA, per(1), pressure, impose_bc_euler, init_time)

      !===Read data setup
   END SUBROUTINE start_setup

   SUBROUTINE read_setup_data(this, rank)
      USE character_strings
      IMPLICIT NONE
      INTEGER, PARAMETER :: in_unit = 21, list_length=200, length_begin=28, length_end=26
      CHARACTER(LEN=length_begin), PARAMETER :: begin_section ='%%% BEGIN SECTION: SETUP %%%'
      CHARACTER(LEN=length_end),   PARAMETER :: end_section   ='%%% END SECTION: SETUP %%%'
      CHARACTER(LEN=length_begin), PARAMETER :: char_begin    ='%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      CHARACTER(LEN=length_end),   PARAMETER :: char_end      ='%%%%%%%%%%%%%%%%%%%%%%%%%%'
      
      CLASS(setup_data_type)             :: this
      TYPE(argument_setup_data_type)     :: argument_data

      CHARACTER(LEN=rec_length), DIMENSION(list_length) :: list, record
      CHARACTER(LEN=rec_length)     :: string_default
      LOGICAL :: okay
      INTEGER :: rank, record_size, i_list, j


      !===Initialize data to zero and false by default
      list = ''
      record = ''
      
      !===Initializing record
      CALL read_data_in_record(record_size, record, begin_section, end_section)

      !===Now we reorganize record

      i_list = 1
      WRITE(list(i_list), *) REPEAT('|', 70)
      i_list = i_list + 1
      list(i_list) = char_begin
      i_list = i_list + 1
      list(i_list) = begin_section
      i_list = i_list + 1
      list(i_list) = char_begin

      !===Restart
      WRITE(string_default,*) this%if_restart
      CALL compare_string(record, list, argument_data%if_restart, string_default, okay, i_list, j)
      IF (okay) THEN
          READ(list(i_list),*) this%if_restart
      END IF

      !===Checkpointing
      WRITE(string_default,*) this%checkpointing_freq
      CALL compare_string(record, list, argument_data%checkpointing_freq, string_default, okay, i_list, j)
      IF (okay) THEN
          READ(list(i_list),*) this%checkpointing_freq
      END IF

      !===Checkpointing
      WRITE(string_default,*) this%final_time
      CALL compare_string(record, list, argument_data%final_time, string_default, okay, i_list, j)
      IF (okay) THEN
          READ(list(i_list),*) this%final_time
      END IF

      !===Regression test
      CALL getarg(1, string_default)
      IF (trim(adjustl(string_default))=='regression') THEN
         this%if_regression_test = .true.
      ELSE
         this%if_regression_test = .false.
      END IF

      i_list = i_list + 1
      list(i_list) = char_end
      i_list = i_list + 1
      list(i_list) = end_section
      i_list = i_list + 1
      list(i_list) = char_end
      !===Closing unit
      CALL rewrite_data_from_list_record(rank, list, record, i_list, record_size)
   END SUBROUTINE read_setup_data


END MODULE start_setup_MODULE
