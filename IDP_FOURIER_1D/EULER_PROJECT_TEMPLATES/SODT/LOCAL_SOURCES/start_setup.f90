MODULE start_setup_MODULE
  USE setup_module
  USE fourier_param_module
  USE euler_module
  TYPE argument_setup_data_type
     CHARACTER(LEN=rec_length) :: if_restart         = '=== Restart (true/false) ==='
     CHARACTER(LEN=rec_length) :: checkpointing_freq = '=== Checkpointing frequency ==='
     CHARACTER(LEN=rec_length) :: final_time         = '=== Final time ==='
  END TYPE argument_setup_data_type
  TYPE setup_data_type
     LOGICAL        :: if_restart          = .FALSE. 
     REAL(KIND = 8) :: checkpointing_freq  = 1.d20
     REAL(KIND = 8) :: final_time          = 0.1d0
     INTEGER :: syst_size
   CONTAINS 
     PROCEDURE, PUBLIC :: read => read_setup_data
     PROCEDURE, PUBLIC :: init => init_setup_data
  END TYPE setup_data_type
  TYPE(setup_data_type)     :: setup_data
  TYPE(fourier_param_type)  :: fourier_param
  TYPE(eos_type)            :: eos            
  TYPE(euler_type) :: euler
CONTAINS
  SUBROUTINE start_setup
#include "petsc/finclude/petsc.h"
    USE petsc
    !USE euler_flux
    USE arbitrary_eos_lambda_module
    IMPLICIT NONE
    REAL(KIND = 8) :: init_time = 0.d0
    INTEGER :: ierr
    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    CALL fourier_param%init
    CALL setup_data%init
    !CALL euler%init(flux,pressure,lambda_arbitrary_eos,fourier_param,init_time,setup_data%final_time)
    CALL euler%init(eos,lambda_arbitrary_eos,fourier_param,init_time,setup_data%final_time)
  END SUBROUTINE start_setup

  SUBROUTINE init_setup_data(this)
    IMPLICIT NONE
    CLASS(setup_data_type), INTENT(INOUT) :: this
    INTEGER :: i
    CALL this%read()
  END SUBROUTINE init_setup_data

  SUBROUTINE read_setup_data(this)
    USE character_strings
    IMPLICIT NONE
    INTEGER, PARAMETER :: in_unit = 21, list_length=200
    CHARACTER(1), PARAMETER :: begin_section ='~'
    CHARACTER(1), PARAMETER :: end_section   ='~'
    CLASS(setup_data_type)             :: this
    TYPE(argument_setup_data_type)     :: argument_data
    CHARACTER(LEN=rec_length), DIMENSION(list_length) :: list, record
    CHARACTER(LEN=rec_length)     :: string_default
    LOGICAL :: okay
    INTEGER :: rank, record_size, i_list, j

    !===Initialize data to zero and false by default
    list = ''
    record = ''
    i_list = 1

    !===Initializing record
    CALL read_data_in_record(record_size, record, begin_section, end_section)

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

    !===Final time
    WRITE(string_default,*) this%final_time
    CALL compare_string(record, list, argument_data%final_time, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%final_time
    END IF

    !===Closing unit
    rank = 0
    CALL rewrite_data_from_list_record(rank, list, record, i_list, record_size)

  END SUBROUTINE read_setup_data

END MODULE start_setup_MODULE
