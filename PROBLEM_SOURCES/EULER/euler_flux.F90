MODULE euler_flux
  USE space_dim
  USE eos !<===FIX ME
  PUBLIC :: flux
CONTAINS
    FUNCTION flux(comp,un) RESULT(vv)
    IMPLICIT NONE 
    REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
    INTEGER,                            INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(SIZE(un,1),k_dim)      :: vv
    REAL(KIND=8), DIMENSION(SIZE(un,1))            :: H, u
    INTEGER :: k
    SELECT CASE(comp)
    CASE(1)
       DO k = 1, k_dim
          vv(:,k) = un(:,k+1)
       END DO
    CASE(2:k_dim+1)
       u = un(:,comp)/un(:,1)
       DO k = 1, k_dim
          vv(:,k) = un(:,k+1)*u
       END DO
       vv(:,comp) = vv(:,comp) + pressure(un)   
    CASE(k_dim+2) 
       H = pressure(un) + un(:,comp)
       DO k = 1, k_dim
          vv(:,k) = (un(:,k+1)/un(:,1))*H
       END DO
    CASE DEFAULT
       WRITE(*,*) ' BUG in flux'
       STOP
    END SELECT
  END FUNCTION flux
END MODULE euler_flux
