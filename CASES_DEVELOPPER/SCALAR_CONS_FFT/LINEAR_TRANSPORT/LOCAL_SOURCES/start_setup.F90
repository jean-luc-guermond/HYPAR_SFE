MODULE start_setup_MODULE
  USE setup_module
  USE fourier_param_module
  USE nl_scalar_cons_module

  TYPE argument_setup_data_type
     CHARACTER(LEN=rec_length) :: if_restart         = '=== Restart (true/false) ==='
     CHARACTER(LEN=rec_length) :: checkpointing_freq = '=== Checkpointing frequency ==='
     CHARACTER(LEN=rec_length) :: final_time         = '=== Final time ==='
     CHARACTER(LEN=rec_length) :: if_analytical_ref  = '=== Do we compare with analytical reference? (true/false) ==='
     CHARACTER(LEN=rec_length) :: char_which_init    = '=== Which initial condition? (smooth, step) ==='
  END TYPE argument_setup_data_type
  TYPE setup_data_type
     LOGICAL        :: if_restart          = .FALSE.
     REAL(KIND = 8) :: checkpointing_freq  = 1.d20
     REAL(KIND = 8) :: final_time          = 0.1d0
     INTEGER :: syst_size
     LOGICAL        :: if_analytical_ref   = .FALSE.
     CHARACTER(LEN=rec_length) :: char_which_init = 'step'
     INTEGER        :: which_init
  CONTAINS
     PROCEDURE, PUBLIC :: read => read_setup_data
     PROCEDURE, PUBLIC :: init => init_setup_data
  END TYPE setup_data_type

  INTEGER, PRIVATE, PARAMETER :: INIT_STEP_RHO=1, INIT_SMOOTH_RHO=2
  CHARACTER(LEN=20), DIMENSION(2), PRIVATE, PARAMETER  :: list_init_rho = &
              [CHARACTER(LEN=20) :: 'step', 'smooth']

  TYPE(setup_data_type)     :: setup_data
  TYPE(fourier_param_type)  :: fourier_param
  TYPE(nl_scalar_cons_type) :: nl_scalar_cons
CONTAINS
  SUBROUTINE start_setup
    USE petscsysdef
    USE petscmpi
    USE petscdmda

    USE read_inputs_module, ONLY : clean_data_once
    USE my_util, ONLY: error_petsc, to_str
    USE options_module
    IMPLICIT NONE
    REAL(KIND = 8) :: init_time = 0.d0
    INTEGER :: ierr
    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    !===Read executable arguments
    CALL read_all_arguments
    !===Clean data once
    CALL clean_data_once
    !===Init parameters
    CALL fourier_param%init
    CALL setup_data%init
    
    SELECT CASE(setup_data%which_init)
    CASE(INIT_SMOOTH_RHO)
      exact_sol_R => exact_sol_smooth_R
    CASE(INIT_STEP_RHO)
      exact_sol_R => exact_sol_step_R
    CASE DEFAULT
      CALL error_petsc("BUG in start_setup => wrong init case "//to_str(setup_data%which_init))
    END SELECT

    CALL nl_scalar_cons%init(flux,flux_prime,lambda_max,fourier_param,init_time,setup_data%final_time)
  END SUBROUTINE start_setup

  SUBROUTINE init_setup_data(this)
    IMPLICIT NONE
    CLASS(setup_data_type), INTENT(INOUT) :: this
    CALL this%read()
  END SUBROUTINE init_setup_data

   SUBROUTINE read_setup_data(this)
      USE read_inputs_module
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

      !===Final time
    CALL read_data(argument_data%final_time, this%final_time)

      !===if_analytical_ref
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
