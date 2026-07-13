MODULE my_util
  IMPLICIT NONE

  INTERFACE to_str
      MODULE PROCEDURE to_str_int, to_str_real, to_str_int_array_1D, to_str_char_array_1D
  END INTERFACE to_str

  INTERFACE pack_opt
      MODULE PROCEDURE pack_opt_logical, pack_opt_integer, pack_opt_real, pack_opt_character
  END INTERFACE pack_opt

  PUBLIC :: user_time, error_Petsc, local_error_petsc, to_str, write_rank_0, pack_opt, get_tab_idx_char
  PRIVATE
CONTAINS
  !
  !Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
  !

    FUNCTION user_time() RESULT(time)
      IMPLICIT NONE
      REAL(KIND=8) :: time
      INTEGER :: count, count_rate, count_max
      CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
      time = (1.d0*count)/count_rate
    END FUNCTION user_time

    SUBROUTINE get_tab_idx_char(val_in, tab_in, idx_out)
      IMPLICIT NONE
      CHARACTER(LEN=*),               INTENT(IN) :: val_in
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: tab_in
      INTEGER,                       INTENT(OUT) :: idx_out
      
      DO idx_out=1, SIZE(tab_in)
        IF (TRIM(ADJUSTL(tab_in(idx_out)))==TRIM(ADJUSTL(val_in))) RETURN
      END DO
      CALL error_petsc("BUG in get_tab_idx => could not find "//val_in//" inside "//to_str(tab_in))
    END SUBROUTINE get_tab_idx_char


   !========================================================================
   !=========== write on rank 0 subs =======================================
   !========================================================================

    SUBROUTINE error_Petsc(string)
#include "petsc/finclude/petsc.h"
      USE petsc
      IMPLICIT NONE
      CHARACTER(LEN=*),  INTENT(IN) :: string
      PetscErrorCode :: ierr
      CALL write_rank_0(string)
      CALL PetscFinalize(ierr)
      STOP
    END SUBROUTINE error_Petsc

    SUBROUTINE local_error_petsc(string)
#include "petsc/finclude/petsc.h"
      USE petsc
      IMPLICIT NONE
      CHARACTER(LEN=*),  INTENT(IN) :: string
      PetscErrorCode :: ierr
      WRITE(*,*) string
      CALL MPI_Abort(PETSC_COMM_WORLD, 1, ierr)
    END SUBROUTINE local_error_petsc

    SUBROUTINE write_rank_0(string)
#include "petsc/finclude/petsc.h"
      USE petsc
      IMPLICIT NONE
      CHARACTER(LEN=*),  INTENT(IN) :: string
      INTEGER                            :: rank
      PetscErrorCode :: ierr
      CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      IF (rank==0) WRITE(*,*) string
    END SUBROUTINE write_rank_0

   !============================================================================
   !=========== pack_opt interfaces (handle optionals in a more compact way) ===
   !============================================================================

    SUBROUTINE pack_opt_logical(val_out, val_default, opt_val_in)
      IMPLICIT NONE
      LOGICAL, INTENT(IN)  :: val_default
      LOGICAL, INTENT(OUT) :: val_out
      LOGICAL, OPTIONAL, INTENT(IN) :: opt_val_in

      IF (PRESENT(opt_val_in)) THEN
        val_out = opt_val_in
      ELSE
        val_out = val_default
      END IF

    END SUBROUTINE pack_opt_logical

    SUBROUTINE pack_opt_integer(val_out, val_default, opt_val_in)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: val_default
      INTEGER, INTENT(OUT) :: val_out
      INTEGER, OPTIONAL, INTENT(IN) :: opt_val_in

      IF (PRESENT(opt_val_in)) THEN
        val_out = opt_val_in
      ELSE
        val_out = val_default
      END IF

    END SUBROUTINE pack_opt_integer

    SUBROUTINE pack_opt_real(val_out, val_default, opt_val_in)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN)  :: val_default
      REAL(KIND=8), INTENT(OUT) :: val_out
      REAL(KIND=8), OPTIONAL, INTENT(IN) :: opt_val_in

      IF (PRESENT(opt_val_in)) THEN
        val_out = opt_val_in
      ELSE
        val_out = val_default
      END IF

    END SUBROUTINE pack_opt_real

    SUBROUTINE pack_opt_character(val_out, val_default, opt_val_in)
      IMPLICIT NONE
      CHARACTER(LEN=*),              INTENT(IN)  :: val_default
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: val_out
      CHARACTER(LEN=*), OPTIONAL,     INTENT(IN) :: opt_val_in

      IF (PRESENT(opt_val_in)) THEN
        val_out = TRIM(ADJUSTL(opt_val_in))
      ELSE
        val_out = TRIM(ADJUSTL(val_default))
      END IF

    END SUBROUTINE pack_opt_character

   !========================================================================
   !=========== to_str interfaces ==========================================
   !========================================================================
   
    FUNCTION to_str_int(i) RESULT (str)
        INTEGER, INTENT(IN) :: i
        CHARACTER(LEN=:), ALLOCATABLE :: str
        CHARACTER(LEN=32) :: tmp

        WRITE(tmp, '(I0)') i
        str = trim(tmp)
    END FUNCTION to_str_int

    FUNCTION to_str_real(i, opt_precision) RESULT (str)
        REAL(KIND=8), INTENT(IN) :: i
        CHARACTER(LEN=:), ALLOCATABLE :: str
        CHARACTER(LEN=32) :: tmp
        INTEGER, INTENT(IN), OPTIONAL :: opt_precision
        
        IF (PRESENT(opt_precision)) THEN
          IF (ABS(i) < 1.d0) THEN
            WRITE(tmp, '(F'//to_str_int(opt_precision+5)//'.'//to_str_int(opt_precision)//')') i
          ELSE
            WRITE(tmp, '(F0.'//to_str_int(opt_precision)//')') i
          END IF
        ELSE
          WRITE(tmp, '(F0.6)') i
        END IF
        ! WRITE(tmp, '(F0.6)') i
        str = trim(tmp)
    END FUNCTION to_str_real

    FUNCTION to_str_int_array_1D(i) RESULT (str)
        INTEGER, DIMENSION(:), INTENT(IN) :: i
        INTEGER :: n, size_i
        CHARACTER(LEN=:), ALLOCATABLE :: str
        CHARACTER(LEN=32) :: tmp

        size_i = SIZE(i)
        DO n=1, size_i
          WRITE(tmp, '(I0)') i(n)
          IF (n==1) THEN
            str = trim(tmp)
          ELSE
            str = trim(adjustl(str)) // ', '//trim(adjustl(tmp))
          END IF
        END DO
    END FUNCTION to_str_int_array_1D

    FUNCTION to_str_char_array_1D(i) RESULT (str)
        CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: i
        INTEGER :: n, size_i
        CHARACTER(LEN=:), ALLOCATABLE :: str
        CHARACTER(LEN=32) :: tmp

        size_i = SIZE(i)
        DO n=1, size_i
          IF (n==1) THEN
            str = TRIM(ADJUSTL(tmp))
          ELSE
            str = TRIM(ADJUSTL(str)) // ', '//TRIM(ADJUSTL(tmp))
          END IF
        END DO
    END FUNCTION to_str_char_array_1D

END MODULE my_util
