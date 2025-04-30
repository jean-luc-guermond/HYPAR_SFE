MODULE setup
  USE space_dim
  USE eos
  PUBLIC :: sol_anal, init, rho_anal, press_anal, mt_anal, E_anal
  PRIVATE
  REAL(KIND=8) :: x0, x1
CONTAINS
  
  SUBROUTINE init(un,time,rr)
    USE def_of_gamma
    USE lambda_module
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),                INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2),k_dim+2), INTENT(OUT):: un
    REAL(KIND=8),                                INTENT(IN) :: time
    gamma = 1.4d0
    x0=0.1
    x1=0.3
    un(:,1) = rho_anal(time,rr)
    un(:,2) = mt_anal(1,time,rr)
    un(:,3) = E_anal(time,rr)
    CALL set_gamma_for_riemann_solver(gamma)
  END SUBROUTINE init

  FUNCTION rho_anal(time,rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: rr
    REAL(KIND=8),                         INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(SIZE(rr,2))              :: vv
    INTEGER :: n
    IF (SIZE(vv)==0) RETURN
    DO n = 1, SIZE(vv)
       IF ((rr(1,n)-time)<x0 .OR. (rr(1,n)-time)>x1) THEN
          vv(n) = 1.d0
       ELSE
          vv(n) = 1 + (2/(x1-x0))**6*(rr(1,n)-time-x0)**3*(x1-rr(1,n)+time)**3
       END IF
    END DO
  END FUNCTION rho_anal

  FUNCTION press_anal(time,rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8),                        INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    IF (SIZE(vv)==0) RETURN
    vv = 1.d0
  END FUNCTION press_anal

  FUNCTION vit_anal(comp,time,rr) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8),                        INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    IF (SIZE(vv)==0) RETURN
    vv = 1.d0
  END FUNCTION vit_anal

  FUNCTION E_anal(time,rr) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8),                        INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    vv = press_anal(time,rr)/(gamma-1.d0) &
         + rho_anal(time,rr)*(vit_anal(1,time,rr)**2)/2
  END FUNCTION E_anal

  FUNCTION mt_anal(comp,time,rr) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8),                        INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    vv = rho_anal(time,rr)*vit_anal(comp,time,rr)
  END FUNCTION mt_anal



  FUNCTION sol_anal(comp,time,rr) RESULT(vv)
    IMPLICIT NONE 
    INTEGER,                             INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr
    REAL(KIND=8),                        INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv
    SELECT CASE(comp)
    CASE(1)
       vv = rho_anal(time,rr)
    CASE(2)
       vv = mt_anal(1,time,rr)
    CASE(3)
       vv = E_anal(time,rr)
    CASE DEFAULT
       WRITE(*,*) ' BUG in sol_anal'
       STOP
    END SELECT
  END FUNCTION sol_anal
END MODULE setup
