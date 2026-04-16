MODULE start_setup_MODULE
  USE setup_module
  USE fourier_param_module
  USE nl_scalar_cons_module
  INTEGER, PRIVATE, PARAMETER :: rec_length = 200, list_length=200
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
     PROCEDURE, PUBLIC :: init => init_setup_data
  END TYPE setup_data_type
  TYPE(setup_data_type)     :: setup_data
  TYPE(fourier_param_type)  :: fourier_param
  TYPE(nl_scalar_cons_type) :: nl_scalar_cons
CONTAINS
  SUBROUTINE start_setup
#include "petsc/finclude/petsc.h"
    USE petsc
    USE character_strings, ONLY : clean_data_once
    IMPLICIT NONE
    REAL(KIND = 8) :: init_time = 0.d0
    INTEGER :: ierr
    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    !===Clean data once
    CALL clean_data_once
    !===Init parameters
    CALL fourier_param%init
    CALL setup_data%init
    CALL nl_scalar_cons%init(flux,flux_prime,lambda_max,fourier_param,init_time,setup_data%final_time)
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
