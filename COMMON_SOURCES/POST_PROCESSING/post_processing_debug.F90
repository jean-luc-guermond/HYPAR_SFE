MODULE post_processing_debug_MODULE
  IMPLICIT NONE
CONTAINS
  
  SUBROUTINE get_num_test(num_test)
     IMPLICIT NONE
      INTEGER, INTENT(OUT) :: num_test
      CHARACTER(LEN=100) :: string

      CALL getarg(2, string)
      IF (trim(adjustl(string))=='1') THEN
         num_test = 1
      ELSE IF (trim(adjustl(string))=='2') THEN
         num_test = 2
      ELSE IF (trim(adjustl(string))=='3') THEN
         num_test = 3
      ELSE
         WRITE(*,*) "Invalid test number ", trim(adjustl(string)), ". Allowed: 1, 2, 3"
         STOP
      END IF
  END SUBROUTINE get_num_test

  SUBROUTINE regression(absolute_error, opt_num_test)
      INTEGER, OPTIONAL                      :: opt_num_test
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: absolute_error
      REAL(KIND=8), DIMENSION(SIZE(absolute_error)) :: regression_ref, relative_error
      REAL(KIND=8)                           :: tol=1.d-7
      INTEGER, PARAMETER                     :: in_unit=21, out_unit=22
      INTEGER                                :: n, num_test, size_regression_ref
      CHARACTER(LEN=1)                       :: str_num_test
      CHARACTER(LEN=10)                      :: regex
      LOGICAL                                :: test_passed=.TRUE.

!==== Defining regression test number
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

!==== Opening regression reference file + seeing if existent or not
      OPEN(in_unit, FILE='regression_reference_'//trim(adjustl(str_num_test)), STATUS='UNKNOWN', FORM='FORMATTED')
      OPEN(out_unit, FILE='current_regression_reference_'//trim(adjustl(str_num_test)), STATUS='UNKNOWN', FORM='FORMATTED')
      size_regression_ref = 0
      
      DO WHILE(size_regression_ref < size(absolute_error))
         READ(in_unit,*,END=100) regression_ref(size_regression_ref+1)
         size_regression_ref = size_regression_ref + 1
      END DO
      100 CONTINUE

!==== Write current_regression anyway
      DO n=1, SIZE(absolute_error)
         WRITE(out_unit, *) absolute_error(n)
      END DO

!==== Regression test
      test_passed = .TRUE.
      !==== If regression reference file too small, skip test and write error
      IF (size_regression_ref < SIZE(absolute_error)) THEN
         WRITE(*,*) "ERROR IN REGRESSION: size(reference)<size(absolute_error)", size_regression_ref, SIZE(absolute_error)
         test_passed = .FALSE.
      !==== If regression reference large enough, perform regression test
      ELSE
         relative_error = ABS(absolute_error(:) - regression_ref(:size_regression_ref))/ABS(regression_ref(:size_regression_ref))
         DO n=1, SIZE(absolute_error)
             IF (relative_error(n) < tol) THEN
               WRITE(*,*) "Regression test component", n, "passed"
               WRITE(*,*) "Relative error ", relative_error(n)
             ELSE
               WRITE(*,*) "Regression test component", n, "failed"
               WRITE(*,*) "Relative error ", relative_error(n)
               test_passed = .FALSE.
             END IF
         END DO
      END IF

      CLOSE(in_unit)
      CLOSE(out_unit)

      IF (test_passed) THEN
            WRITE(*,*) "Regression test number ", num_test, "passed ", regex
      ELSE
            WRITE(*,*) "Regression test number ", num_test, "failed"
      END IF
     
  END SUBROUTINE regression

END MODULE post_processing_debug_MODULE