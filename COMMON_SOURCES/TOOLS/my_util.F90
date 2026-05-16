MODULE my_util
  IMPLICIT NONE

  INTERFACE to_str
      MODULE PROCEDURE to_str_int, to_str_real
  END INTERFACE to_str

  INTERFACE pack_opt
      MODULE PROCEDURE pack_opt_logical, pack_opt_real
  END INTERFACE pack_opt

  PUBLIC :: user_time, error_Petsc, to_str, write_rank_0, pack_opt
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

   FUNCTION to_str_real(i) RESULT (str)
      REAL(KIND=8), INTENT(IN) :: i
      CHARACTER(LEN=:), ALLOCATABLE :: str
      CHARACTER(LEN=32) :: tmp

      WRITE(tmp, '(F0.6)') i
      str = trim(tmp)
   END FUNCTION to_str_real

END MODULE my_util
