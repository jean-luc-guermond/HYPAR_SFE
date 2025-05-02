MODULE setup_data_module
   IMPLICIT NONE
   TYPE setup_data_type
    LOGICAL                        :: if_regression_test
    LOGICAL                        :: if_restart
    REAL(KIND=8)                   :: checkpointing_freq
    REAL(KIND=8)                   :: Tfinal
    INTEGER                        :: syst_size
   CONTAINS
      PROCEDURE, PUBLIC :: init
   END TYPE setup_data_type
CONTAINS
   SUBROUTINE init(a)
      CLASS(setup_data_type), INTENT(INOUT) :: a
      !===Logicals
      a%if_regression_test = .FALSE.
      a%if_restart = .FALSE.
      !===Reals
      a%checkpointing_freq=1.d20
      a%Tfinal=0.d0
      !===Characters
      !===Integers
    END SUBROUTINE init
END MODULE setup_data_module

MODULE input_setup_data
   USE setup_data_module
   IMPLICIT NONE
   PUBLIC :: read_setup_data
   TYPE(setup_data_type), PUBLIC :: setup_data
   PRIVATE
CONTAINS
  SUBROUTINE read_setup_data(data_fichier)
    USE character_strings
    IMPLICIT NONE
    INTEGER, PARAMETER :: in_unit = 21
    CHARACTER(len = *), INTENT(IN) :: data_fichier
    CHARACTER(LEN = 100) :: argument
    LOGICAL :: okay 
    !===Initialize data to zero and false by default
    CALL setup_data%init


    !===Restart
    CALL find_string(in_unit, "===Restart (true/false)===",okay)
    IF (okay) THEN
       READ (in_unit,*) setup_data%if_restart
    END IF
    CALL find_string(in_unit, "===Checkpointing frequency===",okay)
    IF (okay) THEN
       READ (in_unit,*) setup_data%checkpointing_freq
    END IF

    !===Regression test
    CALL getarg(1,argument)
    IF (trim(adjustl(argument))=='regression') THEN
       setup_data%if_regression_test = .true.
    ELSE
       setup_data%if_regression_test = .false.
    END IF

    CLOSE(in_unit)
  END SUBROUTINE read_setup_data
    
END MODULE input_setup_data
