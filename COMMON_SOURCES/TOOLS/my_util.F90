MODULE my_util
  IMPLICIT NONE
  PUBLIC :: user_time, error_Petsc, to_str, WRITE_rank_0
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

  SUBROUTINE error_Petsc(string)
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    CHARACTER(LEN=*),       INTENT(IN) :: string
    INTEGER                            :: rank
    PetscErrorCode :: ierr
    CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
    IF (rank==0) WRITE(*,*) string
    CALL PetscFinalize(ierr)
    STOP
  END SUBROUTINE error_Petsc

  SUBROUTINE WRITE_rank_0(string)
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    CHARACTER(LEN=*),       INTENT(IN) :: string
    INTEGER                            :: rank
    PetscErrorCode :: ierr
    CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
    IF (rank==0) WRITE(*,*) string
  END SUBROUTINE WRITE_rank_0

   !========================================================================
   
   FUNCTION to_str(i) RESULT (str)
      INTEGER, INTENT(IN) :: i
      CHARACTER(LEN=:), ALLOCATABLE :: str
      CHARACTER(LEN=32) :: tmp

      WRITE(tmp, '(I0)') i
      str = trim(tmp)
   END FUNCTION to_str

END MODULE my_util
