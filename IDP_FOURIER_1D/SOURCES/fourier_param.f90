MODULE fourier_param_module
  INTEGER, PARAMETER, PRIVATE :: rec_length = 200, list_length=200
  TYPE argument_fourier_param_type
     CHARACTER(LEN=rec_length) :: Nmax       = '=== Number of Fourier modes ? ==='
     CHARACTER(LEN=rec_length) :: Length     = '=== Fourier periodic domain Length ? ==='
  END type argument_fourier_param_type
  TYPE fourier_param_type
     INTEGER :: Nmax = 10
     REAL(KIND=8) :: Length=1.d0
     INTEGER :: Nmax_real
     REAL(KIND=8) :: dx
   CONTAINS
     PROCEDURE, PUBLIC :: init => init_fourier_param
     PROCEDURE, PUBLIC :: read => read_fourier_param
  END type fourier_param_type
CONTAINS
  SUBROUTINE init_fourier_param(this)
    IMPLICIT NONE
    CLASS(fourier_param_type), INTENT(INOUT) :: this
    CALL this%read()
    this%Nmax_real = 2*this%Nmax-1
    this%dx = this%length/this%Nmax_real
  END SUBROUTINE init_fourier_param

  SUBROUTINE read_fourier_param(this)
    USE character_strings
    IMPLICIT NONE
    CHARACTER(1), PARAMETER :: begin_section ='~'
    CHARACTER(1), PARAMETER :: end_section   ='~'
    INTEGER, PARAMETER :: in_unit = 21
    CLASS(fourier_param_type), INTENT(INOUT):: this
    TYPE(argument_fourier_param_type)  :: argument_data
    CHARACTER(LEN=rec_length), DIMENSION(list_length) :: list, record
    CHARACTER(LEN=rec_length)                         :: string_default
    LOGICAL :: okay
    INTEGER :: rank, record_size, i_list, j

    !===Initialize data to zero and false by default
    list = ''
    record = ''
    i_list = 1
    
    !===Initializing record
    CALL read_data_in_record(record_size, record, begin_section, end_section)

    !===Nmax
    WRITE(string_default,*) this%Nmax
    CALL compare_string(record, list, argument_data%Nmax, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%Nmax
    END IF

    !===Length
    WRITE(string_default,*) this%Length
    CALL compare_string(record, list, argument_data%Length, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%Length
    END IF

    !===Closing unit
    rank = 0
    CALL rewrite_data_from_list_record(rank, list, record, i_list, record_size)     
  END SUBROUTINE read_fourier_param
END MODULE fourier_param_module
