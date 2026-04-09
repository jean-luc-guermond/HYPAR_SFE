MODULE post_processing_debug_MODULE
  IMPLICIT NONE
CONTAINS
  
  SUBROUTINE regression(absolute_error, opt_num_test)
     INTEGER, OPTIONAL                      :: opt_num_test
     REAL(KIND=8), DIMENSION(:), INTENT(IN) :: absolute_error
     REAL(KIND=8)                           :: regression_ref, relative_error
     REAL(KIND=8)                           :: tol=1.d-7
     INTEGER, PARAMETER                     :: in_unit=21, out_unit=22
     INTEGER                                :: n, num_test
     CHARACTER(LEN=1)                       :: str_num_test
     CHARACTER(LEN=10)                      :: regex
     LOGICAL                                :: test_passed=.TRUE.

     IF (PRESENT(opt_num_test)) THEN
        num_test = opt_num_test
     ELSE
        num_test = 1
     END IF
     WRITE(str_num_test, '(I0)') num_test
     IF (num_test==1) THEN
        regex = '1234567891'
     ELSE IF (num_test==2) THEN
        regex = '2345678912'
     ELSE IF (num_test==3) THEN
        regex = '3456789123'
     ELSE
        WRITE(*,*) "Invalid test number ", num_test, ". Allowed: 1, 2, 3"
        RETURN
     END IF

     OPEN(in_unit, FILE='regression_reference_'//trim(adjustl(str_num_test)), STATUS='UNKNOWN', FORM='FORMATTED')
     OPEN(out_unit, FILE='current_regression_reference_'//trim(adjustl(str_num_test)), STATUS='UNKNOWN', FORM='FORMATTED')
     
     DO n=1, SIZE(absolute_error)
        READ(in_unit, *) regression_ref
        WRITE(out_unit, *) absolute_error(n)
        relative_error = ABS(absolute_error(n) - regression_ref)/ABS(regression_ref)
        IF (relative_error < tol) THEN
            WRITE(*,*) "Regression test component", n, "passed"
            WRITE(*,*) "Relative error ", relative_error
        ELSE
            WRITE(*,*) "Regression test component", n, "failed"
            WRITE(*,*) "Relative error ", relative_error
            test_passed = .FALSE.
        END IF
     END DO

     CLOSE(in_unit)
     CLOSE(out_unit)

     IF (test_passed) THEN
         WRITE(*,*) "Regression test number ", num_test, "passed ", regex
     ELSE
         WRITE(*,*) "Regression test number ", num_test, "failed"
     END IF
     
  END SUBROUTINE regression

END MODULE post_processing_debug_MODULE