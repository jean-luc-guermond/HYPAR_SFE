!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE character_strings
   USE my_util

   PUBLIC :: last_c_leng, last_of_string, start_of_string, read_data
   
   INTERFACE read_data
      MODULE PROCEDURE read_real_data, read_integer_data, read_integer_array_data, read_character_data, read_logical_data
   END INTERFACE read_data
   
   
   INTEGER, PARAMETER, PRIVATE :: rec_length=200,list_length=200
   CHARACTER(LEN = rec_length), DIMENSION(:), ALLOCATABLE, PRIVATE :: record_info_from_data, list_info_for_new_data
   INTEGER, PARAMETER, PRIVATE :: in_unit=21
   INTEGER, PRIVATE :: index_list_info_data, record_size
   LOGICAL, PRIVATE :: data_cleaned = .FALSE.
   
CONTAINS
   
   FUNCTION last_c_leng (len_str, string) RESULT (leng)    
      IMPLICIT NONE
      INTEGER, INTENT(IN)                 :: len_str
      CHARACTER (LEN=len_str), INTENT(IN) :: string
      INTEGER                             :: leng
      INTEGER                             :: i
      
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
      INTEGER, INTENT(IN)                 :: len_str
      CHARACTER (LEN=len_str), INTENT(IN) :: string
      INTEGER                             :: leng
      INTEGER                             :: i
      
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
      CHARACTER (LEN=*), INTENT(IN)       :: string
      INTEGER                             :: start
      INTEGER                             :: i
      
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
      CHARACTER (LEN=*), INTENT(IN)       :: string
      INTEGER                             :: last
      INTEGER                             :: i
      
      last = 1
      
      DO i = LEN(string), 1, -1
         IF ( string(i:i) .NE. ' ' ) THEN
            last = i; EXIT
         ENDIF
      ENDDO
      
   END FUNCTION last_of_string
   !========================================================================
   
   FUNCTION itoa(i) RESULT (str)
      INTEGER, INTENT(IN) :: i
      CHARACTER(LEN=:), ALLOCATABLE :: str
      CHARACTER(LEN=32) :: tmp

      WRITE(tmp, '(I0)') i
      str = trim(tmp)
   END FUNCTION itoa

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
   
   SUBROUTINE compare_string_to_record(argument, string_default, okay, end_idx_record)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)   :: argument
      CHARACTER(LEN=*), INTENT(IN)   :: string_default
      LOGICAL, INTENT(OUT)           :: okay
      INTEGER, INTENT(OUT)           :: end_idx_record
      INTEGER                        :: i, j
      
      okay = .TRUE.
      index_list_info_data = index_list_info_data+1
      list_info_for_new_data(index_list_info_data) = argument
      DO i = 1, SIZE(record_info_from_data)
         IF (TRIM(ADJUSTL(record_info_from_data(i)))==list_info_for_new_data(index_list_info_data)) THEN
            record_info_from_data(i) = ''
            index_list_info_data = index_list_info_data + 1
            list_info_for_new_data(index_list_info_data) = record_info_from_data(i+1)
            record_info_from_data(i+1) = ''
            end_idx_record = i+1
            RETURN
         END IF
      END DO
      WRITE(*,*) ' File reading error for list(index_list_info_data) =   ', TRIM(ADJUSTL(list_info_for_new_data(index_list_info_data)))
      index_list_info_data = index_list_info_data+1
      list_info_for_new_data(index_list_info_data) = string_default
      okay = .FALSE.
      RETURN
   END SUBROUTINE compare_string_to_record
   

   SUBROUTINE add_dummy_string_to_record(start_idx, str_added)
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: start_idx
      INTEGER                        :: j
      LOGICAL, INTENT(OUT), OPTIONAL :: str_added
      
      IF (PRESENT(str_added)) str_added = .FALSE.
      j = start_idx + 1
      DO WHILE(j < SIZE(record_info_from_data))
         IF (TRIM(ADJUSTL(record_info_from_data(j)(1:3)))=='===') THEN
            EXIT
         ELSE
            IF (TRIM(ADJUSTL(record_info_from_data(j))) /= '') THEN
               index_list_info_data = index_list_info_data + 1
               list_info_for_new_data(index_list_info_data) = record_info_from_data(j)
               IF (PRESENT(str_added)) str_added = .TRUE.
               record_info_from_data(j) = ''
            END IF
            j = j + 1
         END IF 
      END DO
   END SUBROUTINE add_dummy_string_to_record
   !===================================================================================
   !========== MANDATORY DATA CLEANING BEFORE ANYTHING
   !===================================================================================
   
   
   SUBROUTINE clean_data_once
#include "petsc/finclude/petsc.h"
      USE PETSC
      IMPLICIT NONE
      INTEGER                                            :: rank, code, record_size_clean, raw_record_size, j
      CHARACTER(LEN=rec_length), DIMENSION(list_length)  :: record, raw_record
      CHARACTER(LEN=rec_length)                          :: control 
      
      PetscErrorCode :: ierr
      
      !===Make sure data is cleaned only once
      data_cleaned = .TRUE.
      
      CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
      
      IF (rank == 0) THEN 
         raw_record_size = 0
         record_size_clean = 0
         
         !===Cleaning data
         OPEN(UNIT = in_unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
         DO
            READ(in_unit,'(A)',END=100) control
            !===Will be used to rewrite the previous data as previous_data
            raw_record_size = raw_record_size + 1
            raw_record(raw_record_size) = control
            
            !===All sets of characters to be cleaned from data file
            IF (TRIM(ADJUSTL(control(1:8)))=="%%%%%%%%") THEN
               CYCLE
            ELSE IF (TRIM(ADJUSTL(control(1:8)))=="||||||||") THEN
               CYCLE
            ELSE IF (TRIM(ADJUSTL(control(1:8)))=="========") THEN
               CYCLE
            ELSE
               record_size_clean = record_size_clean + 1
               record(record_size_clean) = control
            END IF
         END DO
         100 CONTINUE
         CLOSE(in_unit)
         !===Rewriting data as previous_data
         OPEN(UNIT = in_unit, FILE = 'previous_data', FORM = 'formatted', STATUS = 'unknown')
         DO j = 1, raw_record_size
            IF (TRIM(ADJUSTL(raw_record(j)))=='') CYCLE
            WRITE(in_unit,'(A)') TRIM(ADJUSTL(raw_record(j)))
         END DO
         CLOSE(in_unit)
         
         !===Rewriting data after cleaning
         OPEN(UNIT = in_unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
         DO j = 1, record_size_clean
            IF (TRIM(ADJUSTL(record(j)))=='') CYCLE
            WRITE(in_unit,'(A)') TRIM(ADJUSTL(record(j)))
         END DO
         CLOSE(in_unit)

      END IF
      CALL MPI_BARRIER(PETSC_COMM_WORLD, code)
   END SUBROUTINE clean_data_once
   
   !===================================================================================
   !========== MANDATORY INITIALIZING DATA READING
   !===================================================================================
   
   SUBROUTINE read_data_init_list(raw_section_name)
      USE PETSC
      IMPLICIT NONE
      
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL       :: raw_section_name
      CHARACTER(LEN=rec_length)                    :: section_bounds, section_name
      
      CHARACTER(LEN=rec_length) :: control, string
      CHARACTER(LEN=5)   :: fmt
      INTEGER            :: length_section_name, line_section
      INTEGER            :: code
      
      !========== MANDATORY INITIALIZING record_info_from_data AND list_info_for_new_data
      
      IF (.NOT. data_cleaned) THEN
         CALL error_petsc('BUG in character_strings.F90 (read_data_init_list): you should have called "clean_data_once" before reading data for the first time')
      END IF
      
      IF (ALLOCATED(record_info_from_data) .OR. ALLOCATED(list_info_for_new_data)) THEN
         CALL error_petsc('BUG in character_strings.F90 (read_data_in_record): &
      record_info_from_data or list_info_for_new_data is allocated &
      , you might have forgotten to deallocate (for instance by calling "rewrite_data_from_list_record")')
      ELSE 
         ALLOCATE(record_info_from_data(list_length))
         record_info_from_data = ""
         ALLOCATE(list_info_for_new_data(list_length))
         list_info_for_new_data = ""
         index_list_info_data = 0
         record_size = 0
      END IF
      
      !=======================================================
      !========== TYPING IN SECTION NAME (if REPEAT('=', ...) is modified, make sure to include the modification in "clean_data_once") 
      !=======================================================
      
      IF (PRESENT(raw_section_name)) THEN
         
         section_name = '!' // repeat(' ', 11) // TRIM(ADJUSTL(raw_section_name))
         section_bounds = REPEAT('=', LEN(TRIM(ADJUSTL(section_name)))+10)
         
         index_list_info_data = index_list_info_data + 1
         list_info_for_new_data(index_list_info_data) = section_bounds
         index_list_info_data = index_list_info_data + 1
         list_info_for_new_data(index_list_info_data) = TRIM(ADJUSTL(section_name))
         index_list_info_data = index_list_info_data + 1
         list_info_for_new_data(index_list_info_data) = section_bounds
         
      END IF
      
      
      !========== READING CURRENT INFORMATION FROM DATA FILE
      
      OPEN(UNIT = in_unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')
      IF (PRESENT(raw_section_name)) THEN
         length_section_name = LEN(TRIM(ADJUSTL(section_name))) 
      END IF
      
      line_section = -1
      
      !===Read data file into record
      DO
         READ(in_unit,'(A)',END=100) control
         IF (TRIM(ADJUSTL(control))=='') CYCLE
         record_size = record_size+1
         record_info_from_data(record_size)=control
         IF (PRESENT(raw_section_name)) THEN
            WRITE(fmt, '("(A", I0, ")")') length_section_name
            WRITE(string,fmt) TRIM(ADJUSTL(control))
            IF (string==section_name) THEN
               record_info_from_data(record_size) = ""
               record_size = record_size - 1
            END IF
         END IF
      END DO
      100 CONTINUE
      CLOSE(in_unit)
      
      CALL MPI_BARRIER(PETSC_COMM_WORLD, code)

   END SUBROUTINE read_data_init_list
   
   SUBROUTINE finalize_rewrite_data
#include "petsc/finclude/petsc.h"
      USE petsc
      IMPLICIT NONE
      
      INTEGER                                    :: j, rank, code
      PetscErrorCode :: ierr
      
      CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
      
      !=== WRITING NEW DATA IN DATA FILE
      OPEN(unit=in_unit,file='data',FORM='FORMATTED',STATUS='UNKNOWN')
      IF (rank == 0) THEN 
         DO j = 1, record_size
            IF (TRIM(ADJUSTL(record_info_from_data(j)))=='') CYCLE
            WRITE(in_unit,'(A)') TRIM(ADJUSTL(record_info_from_data(j)))
         END DO
         DO j = 1, index_list_info_data
            WRITE(in_unit,'(A)') TRIM(ADJUSTL(list_info_for_new_data(j)))
         END DO
      END IF
      CLOSE(in_unit)
      
      !=== REINITIALIZING DATA TO ZERO AND FALSE BY DEFAULT
      DEALLOCATE(record_info_from_data)
      DEALLOCATE(list_info_for_new_data)
      index_list_info_data = 0
      record_size = 0
      
      !=== WAITING ALL PROCESSES TO FINISH WRITING
      CALL MPI_BARRIER(PETSC_COMM_WORLD, code)

   END SUBROUTINE finalize_rewrite_data
   
   !==========================================================================================================================
   !==========================================================================================================================
   !========          SUBROUTINES TO READ DATA IN A RECORD DEPENDING ON THE TYPE OF DATA          ============================
   !==========================================================================================================================
   !==========================================================================================================================

   SUBROUTINE read_real_data(argument, val_in_out)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: argument
      CHARACTER(LEN=rec_length)    :: string_default
      REAL(KIND=8), INTENT(INOUT)  :: val_in_out
      INTEGER                      :: end_idx_record
      LOGICAL                      :: okay

      WRITE(string_default,*) val_in_out
      CALL compare_string_to_record(argument, string_default, okay, end_idx_record)
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record)

   END SUBROUTINE read_real_data
   
   
   SUBROUTINE read_integer_data(argument, val_in_out)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: argument
      CHARACTER(LEN=rec_length)    :: string_default
      LOGICAL                      :: okay
      INTEGER, INTENT(INOUT)       :: val_in_out
      INTEGER                      :: end_idx_record

      WRITE(string_default,*) val_in_out
      CALL compare_string_to_record(argument, string_default, okay, end_idx_record)
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record)

   END SUBROUTINE read_integer_data
   
   SUBROUTINE read_integer_array_data(argument, val_in_out, opt_skip_data)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: argument
      CHARACTER(LEN=rec_length)    :: string_default
      LOGICAL                      :: okay
      INTEGER, INTENT(INOUT)       :: val_in_out(:)
      INTEGER                      :: end_idx_record
      LOGICAL                      :: str_added, raw_okay
      LOGICAL, OPTIONAL            :: opt_skip_data

      string_default = "0 0"
      IF (PRESENT(opt_skip_data)) THEN
         IF (.NOT. opt_skip_data) WRITE(string_default,*) val_in_out
      ELSE
         WRITE(string_default,*) val_in_out
      END IF

      CALL compare_string_to_record(argument, string_default, raw_okay, end_idx_record)
      IF (PRESENT(opt_skip_data)) THEN
         okay = raw_okay .AND. (.NOT. opt_skip_data)
         IF (.NOT. okay) end_idx_record = end_idx_record - 1
      END IF
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record, str_added=str_added)

   END SUBROUTINE read_integer_array_data
   
   SUBROUTINE read_character_data(argument, val_in_out)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)             :: argument
      CHARACTER(LEN=rec_length)                :: string_default
      LOGICAL                                  :: okay
      CHARACTER(LEN=rec_length), INTENT(INOUT) :: val_in_out
      INTEGER                                  :: end_idx_record

      WRITE(string_default,*) TRIM(ADJUSTL(val_in_out))
      CALL compare_string_to_record(argument, string_default, okay, end_idx_record)
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record)

   END SUBROUTINE read_character_data
   
   SUBROUTINE read_logical_data(argument, val_in_out)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)             :: argument
      CHARACTER(LEN=rec_length)                :: string_default
      LOGICAL                                  :: okay
      LOGICAL, INTENT(INOUT)                   :: val_in_out
      INTEGER                                  :: end_idx_record

      WRITE(string_default,*) val_in_out
      CALL compare_string_to_record(argument, string_default, okay, end_idx_record)
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record)

   END SUBROUTINE read_logical_data
   
   !==========================================================================================================================
   !==========================================================================================================================
   !========          VERY SPECIFIC SET OF SUBROUTINES TO READ DATA IN A RECORD          =====================================
   !==========================================================================================================================
   !==========================================================================================================================

   !=== This subroutine reads '=== Indices of periodic boundaries and corresponding vectors on ' // trim(adjustl(this%name)) // '? ==='
   SUBROUTINE read_periodic_data(argument_list_periodic, nb_bords, list_periodic, vect_e)
      IMPLICIT NONE
      CHARACTER(LEN=rec_length)              :: string_default, string
      LOGICAL                                :: okay
      INTEGER                                :: i, k, j
      INTEGER,          INTENT(IN)           :: nb_bords
      CHARACTER(LEN=*), INTENT(IN)           :: argument_list_periodic
      INTEGER,          INTENT(INOUT)        :: list_periodic(:, :)
      REAL(KIND=8),     INTENT(INOUT)        :: vect_e(:, :)
      INTEGER                                :: k_dim, end_idx_record

      k_dim = SIZE(vect_e,1)
      string = argument_list_periodic
      string_default = "0 0 0.d0"
      SELECT CASE(k_dim)
      CASE(1)
         string_default = "0 0 0.d0"
      CASE(2)
         string_default = "0 0 0.d0 0.d0"
      CASE DEFAULT
         ! CALL error_petsc('BUG in character_strings.F90 (read_periodic_data):&
         !  k_dim should be 1 or 2 not '//itoa(k_dim))
          write(*,*) "pb in read_periodic"
      END SELECT

      okay = .FALSE.
      index_list_info_data = index_list_info_data+1
      list_info_for_new_data(index_list_info_data) = string
      DO i = 1, SIZE(record_info_from_data)
         !=== detecting if there is Periodic BC
         IF (TRIM(ADJUSTL(record_info_from_data(i)))==list_info_for_new_data(index_list_info_data)) THEN
            j = i
            record_info_from_data(j) = ''
            okay = .TRUE.
            EXIT
         END IF
      END DO
      !=== reading all Periodic BC if detected
      DO k=1, MAX(nb_bords, 1)
         index_list_info_data = index_list_info_data + 1
         IF (okay) THEN
            list_info_for_new_data(index_list_info_data) = record_info_from_data(j+k)
            record_info_from_data(j+k) = ''
         ELSE
            !=== default value if no Periodic BC detected
            list_info_for_new_data(index_list_info_data) = string_default
         END IF
         IF (nb_bords /= 0) THEN
            READ(list_info_for_new_data(index_list_info_data), *) list_periodic(:, k), vect_e(:, k)
         END IF
      END DO
      IF (okay) THEN
         end_idx_record = j + nb_bords
         CALL add_dummy_string_to_record(end_idx_record)
      END IF

   END SUBROUTINE read_periodic_data
   
END MODULE character_strings
