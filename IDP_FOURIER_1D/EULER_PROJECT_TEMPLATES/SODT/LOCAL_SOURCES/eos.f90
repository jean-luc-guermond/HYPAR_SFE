MODULE eos_module
  IMPLICIT NONE
  TYPE eos_type
   CONTAINS
     PROCEDURE, PUBLIC, NOPASS :: pressure => pressure_eos
     PROCEDURE, PUBLIC, NOPASS :: entropy => entropy_eos
  END type eos_type
     REAL(KIND = 8), PUBLIC, PARAMETER :: gamma=1.4d0 !FIX ME<===MUST DISAPPEAR
   CONTAINS
   FUNCTION pressure_eos(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      vv = rho * e * (gamma - 1)
    END FUNCTION pressure_eos
    FUNCTION entropy_eos(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      vv = e/rho**(gamma - 1)
    END FUNCTION entropy_eos
END MODULE eos_module
  
