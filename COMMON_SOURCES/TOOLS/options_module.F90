MODULE options_module
    
    PUBLIC :: read_all_arguments
    PRIVATE

    INTEGER, PRIVATE, PARAMETER :: arg_length = 30

    TYPE options_name_type
        CHARACTER(LEN=arg_length) :: data_in       = 'data_in'
        CHARACTER(LEN=arg_length) :: data_save     = 'data_save'
        CHARACTER(LEN=arg_length) :: data_out      = 'data_out'
        CHARACTER(LEN=arg_length) :: if_regression = 'if_regression'
        CHARACTER(LEN=arg_length) :: num_regex     = 'num_regex'
    END TYPE options_name_type

    TYPE options_type
        CHARACTER(LEN=arg_length) :: data_in   = 'data'
        CHARACTER(LEN=arg_length) :: data_out  = 'data'
        CHARACTER(LEN=arg_length) :: data_save = 'previous_data'
        LOGICAL :: if_regression = .FALSE.
        INTEGER :: num_regex     = 1
    END TYPE options_type

    TYPE(options_type), PUBLIC :: options

CONTAINS

   !==============================
   !=========== read arguments ===
   !==============================

    SUBROUTINE read_all_arguments
        USE my_util, ONLY: write_rank_0, to_str
        IMPLICIT NONE
        TYPE(options_name_type) :: options_name
        INTEGER :: i
        CHARACTER(LEN=:), ALLOCATABLE :: val_read

        CALL read_given_argument(options_name%data_in, val_read)
        IF (val_read /= "NONE") THEN 
            READ(val_read,*) options%data_in
            CALL write_rank_0("Reading... " // options_name%data_in // to_str(options%data_in))
        END IF

        CALL read_given_argument(options_name%data_out, val_read)
        IF (val_read /= "NONE") THEN
            READ(val_read,*) options%data_out
            CALL write_rank_0("Reading... " // options_name%data_out // to_str(options%data_out))
        END IF

        CALL read_given_argument(options_name%data_save, val_read)
        IF (val_read /= "NONE") THEN
            READ(val_read,*) options%data_save
            CALL write_rank_0("Reading... " // options_name%data_save // to_str(options%data_save))
        END IF

        CALL read_given_argument(options_name%if_regression, val_read)
        IF (val_read /= "NONE") THEN
            READ(val_read,*) options%if_regression
            CALL write_rank_0("Reading... " // options_name%if_regression // to_str(options%if_regression))
        END IF

        CALL read_given_argument(options_name%num_regex, val_read)
        IF (val_read /= "NONE") THEN
            READ(val_read,*) options%num_regex
            CALL write_rank_0("Reading... " // options_name%num_regex // to_str(options%num_regex))
        END IF

    END SUBROUTINE read_all_arguments

    SUBROUTINE read_given_argument(argument_in, val_out)
        USE my_util, ONLY: error_petsc
        IMPLICIT NONE
        CHARACTER(LEN=*),               INTENT(IN) :: argument_in
        CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: val_out
        CHARACTER(LEN=:), ALLOCATABLE :: arg_value, cur_arg
        INTEGER :: num_args, i, l, arg_len, idx_arg
        LOGICAL :: okay

        val_out = 'NONE'

        num_args = command_argument_count()
        DO i=1, num_args
            CALL get_command_argument(number=i, length=arg_len)
            
            ALLOCATE(CHARACTER(LEN=arg_len) :: arg_value)
            
            CALL get_command_argument(number=i, value=arg_value)
            okay = .FALSE.
            DO l=1, LEN(arg_value)
                IF (arg_value(l:l)=="=") THEN
                    okay = .TRUE.
                    EXIT
                END IF
            END DO
            IF (.NOT. okay) THEN
                CALL error_petsc("BUG in read_given_argument: attempting to read "//trim(adjustl(arg_value))//&
                &" but the structure is not of the form 'argument=my_value'.")
            END IF

            IF (TRIM(ADJUSTL(argument_in)) == arg_value(1:l-1)) THEN
                val_out = arg_value(l+1:LEN(arg_value))
                DEALLOCATE(arg_value)
                RETURN
            ELSE
                DEALLOCATE(arg_value)
            END IF

        END DO

    END SUBROUTINE read_given_argument

END MODULE options_module