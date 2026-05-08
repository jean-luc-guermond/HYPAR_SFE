MODULE limiting_bounds_euler_module

CONTAINS

    FUNCTION psi_rho_min(x,psi_m) RESULT(v)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
        REAL(KIND=8), INTENT(IN) :: psi_m
        REAL(KIND=8) :: v
        v = x(1)-psi_m
    END FUNCTION psi_rho_min

    FUNCTION zero_of_psi_rho_min(psi_m,u0,P) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u0, P
       REAL(KIND=8), INTENT(IN) :: psi_m
       REAL(KIND=8) :: v
       v = (psi_m-u0(1))/P(1)
    END FUNCTION zero_of_psi_rho_min

    FUNCTION psi_rho_max(x,psi_m) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
       REAL(KIND=8), INTENT(IN) :: psi_m
       REAL(KIND=8) :: v
       v = psi_m-x(1)
    END FUNCTION psi_rho_max

    FUNCTION zero_of_psi_rho_max(psi_m,u0,P) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u0, P
       REAL(KIND=8), INTENT(IN) :: psi_m
       REAL(KIND=8) :: v
       v = (psi_m-u0(1))/P(1)
    END FUNCTION zero_of_psi_rho_max

END MODULE limiting_bounds_euler_module