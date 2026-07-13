MODULE post_processing_debug_MODULE
  IMPLICIT NONE
  PUBLIC :: regression
  PRIVATE
  CHARACTER(LEN=10), PARAMETER :: list_regex(9) = (/ &
   "1234567891",&
   "2345678912",&
   "3456789123",&
   "4567891234",&
   "5678912345",&
   "6789123456",&
   "7891234567",&
   "8912345678",&
   "9123456789"&
  /)

CONTAINS
  
!   SUBROUTINE get_num_test(num_test)
!      IMPLICIT NONE
!       INTEGER, INTENT(OUT) :: num_test
!       CHARACTER(LEN=100) :: string

!       CALL getarg(2, string)
!       IF (trim(adjustl(string))=='1') THEN
!          num_test = 1
!       ELSE IF (trim(adjustl(string))=='2') THEN
!          num_test = 2
!       ELSE IF (trim(adjustl(string))=='3') THEN
!          num_test = 3
!       ELSE
!          WRITE(*,*) "Invalid test number ", trim(adjustl(string)), ". Allowed: 1, 2, 3"
!          STOP
!       END IF
!   END SUBROUTINE get_num_test

  SUBROUTINE regression(absolute_error, opt_num_test, opt_tol)
      USE my_util, ONLY: pack_opt, error_petsc, to_str
      INTEGER, OPTIONAL                      :: opt_num_test
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: absolute_error
      REAL(KIND=8), DIMENSION(SIZE(absolute_error)) :: regression_ref, relative_error
      REAL(KIND=8)                           :: raw_tol=1.d-7, tol
      REAL(KIND=8), OPTIONAL,     INTENT(IN) :: opt_tol
      INTEGER, PARAMETER                     :: in_unit=21, out_unit=22
      INTEGER                                :: n, num_test, size_regression_ref
      CHARACTER(LEN=1)                       :: str_num_test
      CHARACTER(LEN=10)                      :: regex
      LOGICAL                                :: test_passed=.TRUE.

      CALL pack_opt(tol, raw_tol, opt_tol)

!==== Defining regression test number
      CALL pack_opt(num_test, 1, opt_num_test)

      WRITE(str_num_test, '(I0)') num_test
      IF (num_test > SIZE(list_regex)) THEN
         CALL error_petsc("BUG in regression: invalid num_test = "//to_str(num_test)//&
         &", only valid up to "//to_str(SIZE(list_regex)))
      ELSE
         regex = list_regex(num_test)
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
         WRITE(*,*) "BUG IN REGRESSION: size(reference)<size(absolute_error)", size_regression_ref, SIZE(absolute_error)
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
               WRITE(*,*) "Relative error ", relative_error(n), " for tol = ", tol
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