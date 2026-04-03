MODULE post_processing_debug_MODULE
  IMPLICIT NONE
CONTAINS
  
  SUBROUTINE regression(absolute_error_density)
     REAL(KIND=8), INTENT(IN) :: absolute_error_density
     REAL(KIND=8)             :: regression_ref, relative_error_density
     REAL(KIND=8)             :: tol=1.d-7
     INTEGER, PARAMETER       :: in_unit=21, out_unit=22

     OPEN(in_unit, FILE='regression_reference', STATUS='UNKNOWN', FORM='FORMATTED')
     OPEN(out_unit, FILE='current_regression_reference', STATUS='UNKNOWN', FORM='FORMATTED')
     READ(in_unit, *) regression_ref
     WRITE(out_unit, *) absolute_error_density
     CLOSE(in_unit)
     CLOSE(out_unit)

     relative_error_density = ABS(absolute_error_density - regression_ref)/ABS(regression_ref)
     IF (relative_error_density < tol) THEN
         WRITE(*,*) "Regression test passed", "1234567891"
         WRITE(*,*) "Relative error density ", relative_error_density
     ELSE
         WRITE(*,*) "Regression test failed"
         WRITE(*,*) "Relative error density ", relative_error_density
     END IF
     
  END SUBROUTINE regression

END MODULE post_processing_debug_MODULE