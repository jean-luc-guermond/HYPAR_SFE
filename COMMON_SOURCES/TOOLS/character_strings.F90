!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE character_strings

  PUBLIC :: last_c_leng, last_of_string, start_of_string, default_data

CONTAINS

  FUNCTION last_c_leng (len_str, string) RESULT (leng)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: len_str
    CHARACTER (LEN=len_str), INTENT(IN) :: string
    INTEGER :: leng

    INTEGER :: i

    leng = len_str

    DO i=1,len_str
       IF ( string(i:i) .EQ. ' ' ) THEN
          leng = i-1; EXIT
       ENDIF
    ENDDO

  END FUNCTION last_c_leng

  !========================================================================

  FUNCTION  eval_blank(len_str, string) RESULT (leng)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: len_str
    CHARACTER (LEN=len_str), INTENT(IN) :: string
    INTEGER :: leng

    INTEGER :: i

    leng = len_str

    DO i=1,len_str
       IF ( string(i:i) .NE. ' ' ) THEN
          leng = i; EXIT
       ENDIF
    ENDDO

  END FUNCTION eval_blank

  !========================================================================

  FUNCTION start_of_string (string) RESULT (start)

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: string
    INTEGER :: start

    INTEGER :: i

    start = 1

    DO i = 1, LEN(string)
       IF ( string(i:i) .NE. ' ' ) THEN
          start = i; EXIT
       ENDIF
    ENDDO

  END FUNCTION start_of_string

  !========================================================================

  FUNCTION last_of_string (string) RESULT (last)

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: string
    INTEGER :: last

    INTEGER :: i

    last = 1

    DO i = LEN(string), 1, -1
       IF ( string(i:i) .NE. ' ' ) THEN
          last = i; EXIT
       ENDIF
    ENDDO

  END FUNCTION last_of_string
  !========================================================================

  SUBROUTINE read_until(unit, string, error)
    IMPLICIT NONE
    INTEGER, PARAMETER                 :: long_max=128
    INTEGER,                INTENT(IN) :: unit
    CHARACTER(LEN=*),       INTENT(IN) :: string
    CHARACTER(len=long_max)            :: control
    INTEGER                            :: d_end, d_start
    LOGICAL, OPTIONAL                  :: error
    IF (PRESENT(error)) error =.FALSE.
    REWIND(unit)
    DO WHILE (.TRUE.)
       READ(unit,'(64A)',ERR=11,END=22) control
       d_start = start_of_string(control)
       d_end =   last_of_string(control)
       IF (control(d_start:d_end)==string) RETURN
    END DO

    RETURN
11  WRITE(*,*) ' Error reading data file '; IF (PRESENT(error)) error=.TRUE.; RETURN
22  WRITE(*,*) ' Data string ',string,' not found '; IF (PRESENT(error)) error=.TRUE.; RETURN

  END SUBROUTINE read_until

  SUBROUTINE find_string(unit, string, okay)
    IMPLICIT NONE
    INTEGER, PARAMETER                 :: long_max=128
    INTEGER,                INTENT(IN) :: unit
    CHARACTER(LEN=*),       INTENT(IN) :: string
    CHARACTER(len=long_max)            :: control
    LOGICAL                            :: okay

    okay = .TRUE.
    REWIND(unit)
    DO WHILE (.TRUE.)
       READ(unit,'(64A)',ERR=11,END=22) control
       IF (TRIM(ADJUSTL(control))==string) RETURN
    END DO

11  WRITE(*,*) ' File reading error for: ', string; STOP
22  okay = .FALSE.; RETURN

  END SUBROUTINE find_string
  !========================================================================

  SUBROUTINE default_data(rank, in_unit, argument, opt_char)
    IMPLICIT NONE
    INTEGER :: rank, in_unit
    CHARACTER(*), INTENT(IN) :: argument
    CHARACTER(*), OPTIONAL, INTENT(IN) :: opt_char
    IF (rank==0) THEN
       WRITE(*, *)  TRIM(ADJUSTL(argument)) // ' not found.'
       IF (PRESENT(opt_char)) THEN
          WRITE(*, *) 'Set to ' // TRIM(ADJUSTL(opt_char)) // ' by default and added to data file.'
       END IF
    END IF
    CLOSE(in_unit)
    OPEN(UNIT = in_unit, FILE = "data", FORM = 'formatted', STATUS = 'unknown', position = 'append')
    IF (rank==0) THEN
       WRITE(in_unit, '(g0)') TRIM(ADJUSTL(argument))
       IF (PRESENT(opt_char)) WRITE(in_unit, '(g0)') TRIM(ADJUSTL(opt_char))
    END IF
  END SUBROUTINE default_data

  SUBROUTINE compare_string(record, list, string, string_default, okay, i_list, j)
    IMPLICIT NONE
    CHARACTER(LEN=*), DIMENSION(:) :: record, list
    CHARACTER(LEN=*)               :: string, string_default
    LOGICAL, INTENT(OUT)           :: okay
    INTEGER, INTENT(OUT)           :: j
    INTEGER, INTENT(INOUT)         :: i_list
    INTEGER                        :: i
    okay = .TRUE.
    i_list = i_list+1
    list(i_list) = string
    DO i = 1, SIZE(record)
       IF (TRIM(ADJUSTL(record(i)))==list(i_list)) THEN
          j = i
          record(j) = ''
          i_list = i_list + 1
          list(i_list) = record(j+1)
          record(j+1) = ''
          RETURN
       END IF
    END DO
    WRITE(*,*) ' File reading error for list(i_list) =   ', list(i_list)
    i_list = i_list+1
    list(i_list) = string_default
    okay = .FALSE.
    j = -1
    RETURN
  END SUBROUTINE compare_string

  SUBROUTINE read_data_in_record(record_size, record, begin_section, end_section)
    use petsc
    IMPLICIT NONE
    
    INTEGER, INTENT(OUT)                        :: record_size
    CHARACTER(LEN=*), DIMENSION(*), INTENT(OUT) :: record
    CHARACTER(LEN=*), INTENT(IN)                 :: begin_section, end_section
    
    INTEGER, PARAMETER :: in_unit=21
    CHARACTER(LEN=200) :: control, string
    CHARACTER(LEN=5)   :: fmt
    INTEGER            :: length_begin, length_end, line_begin_section, line_end_section
    INTEGER            :: code

    line_begin_section = -1
    line_end_section = -1

  !  record(:) = ""
    record_size = 0

    OPEN(UNIT = in_unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
    length_begin = LEN(TRIM(ADJUSTL(begin_section))) 
    length_end = LEN(TRIM(ADJUSTL(end_section))) 
    !===Read data file into record
    DO
       READ(in_unit,'(A)',END=100) control
       IF (TRIM(ADJUSTL(control))=='') CYCLE
       record_size = record_size+1
       record(record_size)=control
       WRITE(fmt, '("(A", I0, ")")') length_begin
       WRITE(string,fmt) TRIM(ADJUSTL(control))
       IF (string==begin_section) THEN
          line_begin_section = record_size
       END IF
       WRITE(fmt, '("(A", I0, ")")') length_end
       WRITE(string,fmt) TRIM(ADJUSTL(control))
       IF (string==end_section) THEN
          line_end_section = record_size
       END IF
    END DO
100 CONTINUE
    CLOSE(in_unit)
    IF (line_begin_section .NE. -1) THEN
        record(line_begin_section:record_size-1) = record(line_begin_section+1: record_size)
        record_size = record_size -1
        IF (line_end_section > line_begin_section) line_end_section = line_end_section -1
    ELSE
        line_begin_section = 1
    END IF
    IF (line_end_section == -1) THEN
        line_end_section = record_size
    ELSE
        record(line_end_section:record_size-1) = record(line_end_section+1: record_size)
        record_size = record_size -1
    END IF

    CALL MPI_BARRIER(PETSC_COMM_WORLD, code)

  END SUBROUTINE read_data_in_record

  SUBROUTINE rewrite_data_from_list_record(rank, list, record, i_list, record_size, opt_maybe)
    use petsc
    IMPLICIT NONE
    
    CHARACTER(LEN=*), DIMENSION(*), INTENT(IN) :: record, list
    INTEGER,                        INTENT(IN) :: rank, i_list, record_size
    LOGICAL, OPTIONAL                          :: opt_maybe
    LOGICAL                                    :: write_list
    INTEGER, PARAMETER                         :: in_unit=21
    INTEGER                                    :: j, code

!!    list(i_list+1) = ''

    OPEN(unit=in_unit,file='data',FORM='FORMATTED',STATUS='UNKNOWN')
       IF (rank == 0) THEN 
          write_list = .TRUE.
          IF (PRESENT(opt_maybe)) THEN
                  write_list = (.NOT. opt_maybe)
                  IF (write_list) WRITE(*,*) "writing only record"
          ENDIF
          IF (write_list) THEN
             DO j = 1, record_size
                IF (TRIM(ADJUSTL(record(j)))=='') CYCLE
                WRITE(in_unit,'(A)') TRIM(ADJUSTL(record(j)))
             END DO
          ENDIF 
          write_list = .TRUE.
          IF (PRESENT(opt_maybe)) THEN
                  write_list = (opt_maybe)
                  IF (write_list) WRITE(*,*) "writing only list"
          ENDIF
          IF (write_list) THEN
             DO j = 1, i_list
                WRITE(in_unit,'(A)') TRIM(ADJUSTL(list(j)))
             END DO
          ENDIF 
       END IF
    CLOSE(in_unit)

    CALL MPI_BARRIER(PETSC_COMM_WORLD, code)
  
  END SUBROUTINE
END MODULE character_strings
