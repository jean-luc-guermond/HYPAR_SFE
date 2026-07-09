     PURE FUNCTION psi_scalar_rho_max(x,psi_m) RESULT(v)
          IMPLICIT NONE
          REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: x
          REAL(KIND=8), DIMENSION(:), INTENT(IN) :: psi_m
          REAL(KIND=8), DIMENSION(SIZE(psi_m)) :: v
          INTEGER :: n
          DO n=1, SIZE(psi_m)
               v(n) = -psi_m(n)-x(n,1) !<== note the minus sign (max(x) = min(-x))
          END DO
     END FUNCTION psi_scalar_rho_max

     PURE FUNCTION scal_zero_of_psi_scalar_rho_max(psi_m,u0,P) RESULT(v)
          IMPLICIT NONE
          REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u0, P
          REAL(KIND=8), INTENT(IN)  :: psi_m
          REAL(KIND=8) :: v
          v = (-psi_m-u0(1))/P(1)!<== note the minus sign (max(x) = min(-x))
     END FUNCTION scal_zero_of_psi_scalar_rho_max

     PURE FUNCTION zero_of_psi_scalar_rho_max(psi_m,u0,P) RESULT(v)
          IMPLICIT NONE
          REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: u0, P
          REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: psi_m
          REAL(KIND=8), DIMENSION(SIZE(psi_m)) :: v
          INTEGER :: n
          DO n=1, SIZE(psi_m)
               v(n) = (-psi_m(n)-u0(n,1))/P(n,1)!<== note the minus sign (max(x) = min(-x))
          END DO
     END FUNCTION zero_of_psi_scalar_rho_max