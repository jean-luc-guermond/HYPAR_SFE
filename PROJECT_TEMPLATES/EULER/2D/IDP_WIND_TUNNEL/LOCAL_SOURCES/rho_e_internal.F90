MODULE limiting_rho_e_internal_module
    
CONTAINS

    FUNCTION zero_of_psi_euler_rho_e_internal(psi_m,u0,P) RESULT(v)
        USE space_dim
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: u0, P
        REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: psi_m
        REAL(KIND=8), DIMENSION(SIZE(psi_m)) :: a, b, c, discriminant, v
        REAL(KIND=8) :: r1, r2
        INTEGER :: k

        a = 0.5d0*(2.d0*P(:, 1)*P(:, k_dim+2) - SUM(P(:, 2:k_dim+1)**2,DIM=2))
        ! b = u0(:, k_dim+2)*P(:, 1) - SUM(u0(:, 2:k_dim+1)*P(:, 2:k_dim+1),DIM=2) + u0(:, 1)*P(:, k_dim+2)
        b = (u0(:, k_dim+2) - psi_m)*P(:, 1) - SUM(u0(:, 2:k_dim+1)*P(:, 2:k_dim+1),DIM=2) + u0(:, 1)*P(:, k_dim+2)
        c = psi_euler_rho_e_internal(u0, psi_m)

        discriminant = b**2-4.d0*a*c

        DO k=1, k_dim+1
            IF(discriminant(k) < 0) THEN
                v(k) = 1
            ELSE
                r1 = (2.d0*c(k))/(-b(k)-discriminant(k))
                IF (r1 < 0.d0) r1 = 1.d0
                r2 = (-b(k)-SQRT(discriminant(k)))/(2.d0*a(k))
                IF (r2 < 0.d0) r2 = 1.d0
                v(k) = MIN(r1, r2)
            ENDIF
        END DO

    END FUNCTION zero_of_psi_euler_rho_e_internal

    FUNCTION psi_euler_rho_e_internal(x,psi_m) RESULT(v)
        USE space_dim
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: x
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: psi_m
        REAL(KIND=8), DIMENSION(SIZE(psi_m)) :: v

        ! v = x(:, 1)*x(:, k_dim+2) - 0.5d0*SUM(x(:, 2:k_dim+1)**2,DIM=2) - psi_m
        v = x(:, k_dim+2) - 0.5d0*SUM(x(:, 2:k_dim+1)**2,DIM=2)/x(:, 1) - psi_m

    END FUNCTION psi_euler_rho_e_internal

END MODULE limiting_rho_e_internal_module