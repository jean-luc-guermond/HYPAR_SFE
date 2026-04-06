MODULE post_processing_debug_MODULE
  IMPLICIT NONE
CONTAINS
  
  SUBROUTINE regression(absolute_error)
     REAL(KIND=8), DIMENSION(:), INTENT(IN) :: absolute_error
     REAL(KIND=8)                           :: regression_ref, relative_error
     REAL(KIND=8)                           :: tol=1.d-7
     INTEGER, PARAMETER                     :: in_unit=21, out_unit=22
     INTEGER                                :: n
     LOGICAL                                :: test_passed=.TRUE.

     OPEN(in_unit, FILE='regression_reference', STATUS='UNKNOWN', FORM='FORMATTED')
     OPEN(out_unit, FILE='current_regression_reference', STATUS='UNKNOWN', FORM='FORMATTED')
     
     DO n=1, SIZE(absolute_error)
        READ(in_unit, *) regression_ref
        WRITE(out_unit, *) absolute_error(n)
        relative_error = ABS(absolute_error(n) - regression_ref)/ABS(regression_ref)
        IF (relative_error < tol) THEN
            WRITE(*,*) "Regression test number", n, "passed"
            WRITE(*,*) "Relative error density ", relative_error
        ELSE
            WRITE(*,*) "Regression test number", n, "failed"
            WRITE(*,*) "Relative error density ", relative_error
            test_passed = .FALSE.
        END IF
     END DO

     CLOSE(in_unit)
     CLOSE(out_unit)

     IF (test_passed) THEN
         WRITE(*,*) "Regression test passed", "(1234567891)"
     ELSE
         WRITE(*,*) "Regression test failed"
     END IF
     
  END SUBROUTINE regression

END MODULE post_processing_debug_MODULE