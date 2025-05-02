MODULE eos
   IMPLICIT NONE
   REAL(KIND = 8) :: gamma !FIX ME<===MUST DISAPEAR
CONTAINS
   FUNCTION pressure(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      vv = rho * e * (gamma - 1)
   END FUNCTION pressure
END MODULE eos
  
