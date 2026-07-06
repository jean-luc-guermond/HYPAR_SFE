     PURE FUNCTION psi_scalar_rho_min(x,psi_m) RESULT(v)
          IMPLICIT NONE
          REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: x
          REAL(KIND=8), DIMENSION(:), INTENT(IN) :: psi_m
          REAL(KIND=8), DIMENSION(SIZE(psi_m)) :: v
          INTEGER :: n
          DO n=1, SIZE(x, 1)
               v(n) = x(n,1)-psi_m(n)
          END DO
     END FUNCTION psi_scalar_rho_min
     
     PURE FUNCTION scal_zero_of_psi_scalar_rho_min(psi_m,u0,P) RESULT(v)
          IMPLICIT NONE
          REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u0, P
          REAL(KIND=8), INTENT(IN)  :: psi_m
          REAL(KIND=8) :: v
          v = (psi_m-u0(1))/P(1)
     END FUNCTION scal_zero_of_psi_scalar_rho_min

     PURE FUNCTION zero_of_psi_scalar_rho_min(psi_m,u0,P) RESULT(v)
          IMPLICIT NONE
          REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: u0, P
          REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: psi_m
          REAL(KIND=8), DIMENSION(SIZE(psi_m)) :: v
          INTEGER :: n
          DO n=1, SIZE(u0, 1)
               v(n) = (psi_m(n)-u0(n,1))/P(n,1)
          END DO
     END FUNCTION zero_of_psi_scalar_rho_min