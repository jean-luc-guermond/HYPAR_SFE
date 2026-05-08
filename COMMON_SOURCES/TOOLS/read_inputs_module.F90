MODULE read_inputs_module

    IMPLICIT NONE
    PUBLIC :: clean_data_once, read_data_init_list, finalize_rewrite_data, read_data, read_periodic_data

    INTERFACE read_data
        MODULE PROCEDURE read_real_data, read_integer_data, read_integer_array_data, read_character_data, read_logical_data
    END INTERFACE read_data

    INTEGER, PARAMETER, PUBLIC :: rec_length=200,list_length=200
    CHARACTER(LEN=rec_length), DIMENSION(:), ALLOCATABLE, PRIVATE :: record_info_from_data, list_info_for_new_data
    INTEGER, PARAMETER, PRIVATE :: in_unit=21
    INTEGER, PRIVATE :: index_list_info_data, record_size
    LOGICAL, PRIVATE :: data_cleaned = .FALSE., if_regression_test
    INTEGER, PRIVATE :: num_test, num_data_file
    REAL(KIND=8), PRIVATE :: idx_data_file
    CHARACTER(LEN=:), ALLOCATABLE, PRIVATE :: file_in, file_save, file_out
    CHARACTER(LEN=*), PARAMETER, PRIVATE :: file_in_par = 'data', file_save_par = 'previous_data', file_out_par = 'data'


CONTAINS
   
   SUBROUTINE compare_string_to_record(argument, string_default, okay, end_idx_record, opt_add)
      use my_util, ONLY: pack_opt
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)   :: argument
      CHARACTER(LEN=*), INTENT(IN)   :: string_default
      LOGICAL, INTENT(OUT)           :: okay
      INTEGER, INTENT(OUT)           :: end_idx_record
      INTEGER                        :: i
      LOGICAL, OPTIONAL, INTENT(IN)  :: opt_add
      LOGICAL                        :: to_add
      
      okay = .TRUE.
      DO i = 1, SIZE(record_info_from_data)
         IF (TRIM(ADJUSTL(record_info_from_data(i)))==TRIM(ADJUSTL(argument))) THEN
            record_info_from_data(i) = ''
            index_list_info_data = index_list_info_data+1
            list_info_for_new_data(index_list_info_data) = argument
            index_list_info_data = index_list_info_data + 1
            list_info_for_new_data(index_list_info_data) = record_info_from_data(i+1)
            record_info_from_data(i+1) = ''
            end_idx_record = i+1
            RETURN
         END IF
      END DO
!=== okay = false ===> decide whether argument is added to data or not
      okay = .FALSE.
      CALL pack_opt(to_add, .TRUE., opt_add)
!=== warning message if to_add = true
      IF (to_add) THEN
         WRITE(*,*) ' File reading error for list(index_list_info_data) =   '&
               , argument
         index_list_info_data = index_list_info_data+1
         list_info_for_new_data(index_list_info_data) = argument
         index_list_info_data = index_list_info_data+1
         list_info_for_new_data(index_list_info_data) = string_default
      END IF
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
      USE post_processing_debug_MODULE, ONLY: get_num_test
      USE my_util, ONLY: to_str
      IMPLICIT NONE
      INTEGER                                            :: rank, code, record_size_clean, raw_record_size, j
      CHARACTER(LEN=rec_length), DIMENSION(list_length)  :: record, raw_record
      CHARACTER(LEN=rec_length)                          :: control 
      
      PetscErrorCode :: ierr
      
      !===Make sure data is cleaned only once
      data_cleaned = .TRUE.
      
      CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    
!==== FIX VB 08/05/2026 ==> DATA REWRITING HANDLING if_regression_test
      CALL getarg(1, control)
      IF (trim(adjustl(control))=='regression') THEN
         if_regression_test = .TRUE.
         CALL get_num_test(num_test)
         CALL MPI_Comm_Size(PETSC_COMM_WORLD, num_data_file, code)
         file_in = 'data_'//to_str(num_test)
         file_save = 'previous_data_'//to_str(num_test)
         !=== comment/uncomment to replace current data file in regression
         file_out = 'data_regression_'//to_str(num_test)//'_NPROC_'//to_str(num_data_file)
         ! file_out = file_in
         !=== comment/uncomment to replace current data file in regression
      ELSE
         if_regression_test = .FALSE.
         file_in   = file_in_par
         file_save = file_save_par
         file_out  = file_out_par
      END IF
!==== FIX VB 08/05/2026 ==> DATA REWRITING HANDLING if_regression_test

      IF (rank == 0) THEN 
         raw_record_size = 0
         record_size_clean = 0
         
         !===Cleaning data
         OPEN(UNIT = in_unit, FILE = file_in, FORM = 'formatted', STATUS = 'unknown')
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
         OPEN(UNIT = in_unit, FILE = file_save, FORM = 'formatted', STATUS = 'unknown')
         DO j = 1, raw_record_size
            IF (TRIM(ADJUSTL(raw_record(j)))=='') CYCLE
            WRITE(in_unit,'(A)') TRIM(ADJUSTL(raw_record(j)))
         END DO
         CLOSE(in_unit)
         
         !===Rewriting data after cleaning
         OPEN(UNIT = in_unit, FILE = file_out, FORM = 'formatted', STATUS = 'unknown')
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
      USE my_util, ONLY: error_Petsc
      IMPLICIT NONE
      
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL       :: raw_section_name
      CHARACTER(LEN=rec_length)                    :: section_bounds, section_name
      
      CHARACTER(LEN=rec_length) :: control, string
      CHARACTER(LEN=5)   :: fmt
      INTEGER            :: length_section_name, line_section
      INTEGER            :: code
      
      !========== MANDATORY INITIALIZING record_info_from_data AND list_info_for_new_data
      
      IF (.NOT. data_cleaned) THEN
         CALL error_petsc('BUG in character_strings.F90 (read_data_init_list): &
                 you should have called "clean_data_once" before reading data for the first time')
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
         
         length_section_name = LEN(TRIM(ADJUSTL(section_name))) 
      END IF
      
      
      !========== READING CURRENT INFORMATION FROM DATA FILE
      
      OPEN(UNIT = in_unit, FILE = file_out, FORM = 'formatted', STATUS = 'unknown')
      
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
      OPEN(unit=in_unit,file=file_out,FORM='FORMATTED',STATUS='UNKNOWN')
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

   SUBROUTINE concatenate_argument_name(raw_argument, name, argument)
      USE my_util, ONLY: error_petsc
      USE character_strings, ONLY: last_of_string
      IMPLICIT NONE
      CHARACTER(LEN=*),          INTENT(IN)  :: raw_argument
      CHARACTER(LEN=*), OPTIONAL,INTENT(IN)  :: name
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: argument
      INTEGER :: size1, size2, size3

      IF (raw_argument(1:3) /= '===') THEN
         CALL error_petsc("wrong argument "//TRIM(ADJUSTL(raw_argument))//", &
         make sure it starts with '==='")
      END IF

      IF (.NOT. PRESENT(name)) THEN
         argument = raw_argument
      ELSE
         IF (TRIM(ADJUSTL(name))=="") THEN
            argument = raw_argument
            RETURN
         END IF

         size1 = last_of_string(TRIM(ADJUSTL(raw_argument)))
         size2 = last_of_string(TRIM(ADJUSTL(name)))
         size3 = LEN(" ===")
         IF (raw_argument(size1-size3+1:size1) /= " ===") THEN
            CALL error_petsc("BUG in concatenate_argument_name: wrong argument"//TRIM(ADJUSTL(raw_argument))//&
            ", it should end with ' ==='")
         END IF
         argument = raw_argument(:size1-size3) // " " // name(:size2) // " ==="

      END IF
   END SUBROUTINE concatenate_argument_name


   SUBROUTINE read_real_data(raw_argument, val_in_out, opt_name, opt_add)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)           :: raw_argument
      CHARACTER(LEN=:), ALLOCATABLE          :: argument
      CHARACTER(LEN=rec_length)              :: string_default
      REAL(KIND=8),     INTENT(INOUT)        :: val_in_out
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_name
      INTEGER                      :: end_idx_record
      LOGICAL                      :: okay
      LOGICAL, OPTIONAL, INTENT(IN)  :: opt_add

!=== define default value to set in data if inexistant
      WRITE(string_default,*) val_in_out
!=== opt_name => e.g for an argument supposed to appear several times but for different fields
      CALL concatenate_argument_name(raw_argument, opt_name, argument)
!=== opt_name => e.g for an argument supposed to appear several times but for different fields
      CALL compare_string_to_record(argument, string_default, okay, end_idx_record, opt_add=opt_add)
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record)

   END SUBROUTINE read_real_data
   
   
   SUBROUTINE read_integer_data(raw_argument, val_in_out, opt_name, opt_add)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)           :: raw_argument
      CHARACTER(LEN=:), ALLOCATABLE          :: argument
      CHARACTER(LEN=rec_length)    :: string_default
      LOGICAL                      :: okay
      INTEGER, INTENT(INOUT)       :: val_in_out
      INTEGER                      :: end_idx_record
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_name
      LOGICAL, OPTIONAL, INTENT(IN)  :: opt_add

      WRITE(string_default,*) val_in_out
      CALL concatenate_argument_name(raw_argument, opt_name, argument)
      CALL compare_string_to_record(argument, string_default, okay, end_idx_record, opt_add=opt_add)
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record)

   END SUBROUTINE read_integer_data
   
   SUBROUTINE read_integer_array_data(raw_argument, val_in_out, opt_skip_data, opt_name, opt_add)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)           :: raw_argument
      CHARACTER(LEN=:), ALLOCATABLE          :: argument
      CHARACTER(LEN=rec_length)    :: string_default
      LOGICAL                      :: okay
      INTEGER, INTENT(INOUT)       :: val_in_out(:)
      INTEGER                      :: end_idx_record
      LOGICAL                      :: str_added, raw_okay
      LOGICAL, OPTIONAL            :: opt_skip_data
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_name
      LOGICAL, OPTIONAL, INTENT(IN)  :: opt_add

      string_default = "0 0"
      IF (PRESENT(opt_skip_data)) THEN
         IF (.NOT. opt_skip_data) WRITE(string_default,*) val_in_out
      ELSE
         WRITE(string_default,*) val_in_out
      END IF

      CALL concatenate_argument_name(raw_argument, opt_name, argument)
      CALL compare_string_to_record(argument, string_default, raw_okay, end_idx_record, opt_add=opt_add)
      IF (PRESENT(opt_skip_data)) THEN
         okay = raw_okay .AND. (.NOT. opt_skip_data)
         IF (.NOT. okay) end_idx_record = end_idx_record - 1
      END IF
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record, str_added=str_added)

   END SUBROUTINE read_integer_array_data
   
   SUBROUTINE read_character_data(raw_argument, val_in_out, opt_name, opt_add)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)             :: raw_argument
      CHARACTER(LEN=:), ALLOCATABLE            :: argument
      CHARACTER(LEN=rec_length)                :: string_default
      LOGICAL                                  :: okay
      CHARACTER(LEN=rec_length), INTENT(INOUT) :: val_in_out
      INTEGER                                  :: end_idx_record
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: opt_name
      LOGICAL, OPTIONAL, INTENT(IN)  :: opt_add

      WRITE(string_default,*) TRIM(ADJUSTL(val_in_out))
      CALL concatenate_argument_name(raw_argument, opt_name, argument)

      CALL compare_string_to_record(argument, string_default, okay, end_idx_record, opt_add=opt_add)
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record)

   END SUBROUTINE read_character_data
   
   SUBROUTINE read_logical_data(raw_argument, val_in_out, opt_name, opt_add)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)             :: raw_argument
      CHARACTER(LEN=:), ALLOCATABLE            :: argument
      CHARACTER(LEN=rec_length)                :: string_default
      LOGICAL                                  :: okay
      LOGICAL, INTENT(INOUT)                   :: val_in_out
      INTEGER                                  :: end_idx_record
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: opt_name
      LOGICAL, OPTIONAL, INTENT(IN)  :: opt_add

      WRITE(string_default,*) val_in_out
      CALL concatenate_argument_name(raw_argument, opt_name, argument)

      CALL compare_string_to_record(argument, string_default, okay, end_idx_record, opt_add=opt_add)
      IF (okay) READ(list_info_for_new_data(index_list_info_data),*) val_in_out

      CALL add_dummy_string_to_record(end_idx_record)

   END SUBROUTINE read_logical_data
   
   !==========================================================================================================================
   !==========================================================================================================================
   !========          VERY SPECIFIC SET OF SUBROUTINES TO READ DATA IN A RECORD          =====================================
   !==========================================================================================================================
   !==========================================================================================================================

   !=== This subroutine reads '=== Indices of periodic boundaries and corresponding vectors on ' // trim(adjustl(this%name)) // '? ==='
   SUBROUTINE read_periodic_data(argument_list_periodic, nb_bords, list_periodic, vect_e, opt_name)
      USE space_dim
      USE my_util, ONLY:to_str, error_petsc
      IMPLICIT NONE
      CHARACTER(LEN=:), ALLOCATABLE          :: string_default, argument
      LOGICAL                                :: okay
      INTEGER                                :: i, k, j
      INTEGER,          INTENT(IN)           :: nb_bords
      CHARACTER(LEN=*), INTENT(IN)           :: argument_list_periodic
      INTEGER,          INTENT(INOUT)        :: list_periodic(:, :)
      REAL(KIND=8),     INTENT(INOUT)        :: vect_e(:, :)
      INTEGER                                :: end_idx_record
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_name

      CALL concatenate_argument_name(argument_list_periodic, opt_name, argument)

      string_default = "0 0 0.d0"
      SELECT CASE(k_dim)
      CASE(1)
         string_default = "0 0 0.d0"
      CASE(2)
         string_default = "0 0 0.d0 0.d0"
      CASE DEFAULT
         CALL error_petsc('BUG in read_inputs_module.F90 (read_periodic_data):&
          k_dim should be 1 or 2 not '//to_str(k_dim))
      END SELECT

      okay = .FALSE.
      ! index_list_info_data = index_list_info_data+1
      ! list_info_for_new_data(index_list_info_data) = argument
      DO i = 1, SIZE(record_info_from_data)
         !=== detecting if there is Periodic BC
         IF (TRIM(ADJUSTL(record_info_from_data(i)))==argument) THEN
            j = i
            record_info_from_data(j) = ''
            okay = .TRUE.
            EXIT
         END IF
      END DO
      !=== reading all Periodic BC if detected
      IF (okay .OR. nb_bords > 0) THEN !=== write pbc if nb_bords > 0 or was already in data file
         !=== add argument to data
         index_list_info_data = index_list_info_data+1
         list_info_for_new_data(index_list_info_data) = argument
         !=== warning message if wasn't there
         IF (.NOT. okay) WRITE(*,*) ' File reading error for list(index_list_info_data) =   '&
               , argument
         !=== attempting to read pbc data entered by user
         DO k=1, MAX(nb_bords, 1)
            index_list_info_data = index_list_info_data + 1
            !=== effectively reading data as CHARACTER
            IF (okay) THEN
               list_info_for_new_data(index_list_info_data) = record_info_from_data(j+k)
               record_info_from_data(j+k) = ''
            !=== add default value if no Periodic BC detected
            ELSE
               list_info_for_new_data(index_list_info_data) = string_default
            END IF
            !=== converting data to INTEGER/REAL
            IF (nb_bords /= 0) THEN
               READ(list_info_for_new_data(index_list_info_data), *) list_periodic(:, k), vect_e(:, k)
            END IF
         END DO
      END IF
      !=== copy/paste all following lines until the next argument
      IF (okay) THEN
         end_idx_record = j + nb_bords
         CALL add_dummy_string_to_record(end_idx_record)
      END IF

   END SUBROUTINE read_periodic_data
   
END MODULE read_inputs_module
