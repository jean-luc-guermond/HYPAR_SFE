MODULE eos
  IMPLICIT NONE
  REAL(KIND=8) :: gamma !FIX ME<===MUST DISAPEAR
CONTAINS
  FUNCTION pressure(un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(SIZE(un,1))            :: vv
    vv=(un(:,3)-0.5d0*(un(:,2)**2)/un(:,1))*(gamma-1)
  END FUNCTION pressure
END MODULE eos
  
