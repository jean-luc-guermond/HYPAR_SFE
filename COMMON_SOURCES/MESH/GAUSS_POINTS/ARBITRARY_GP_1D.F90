MODULE arbitrary_GP_1D_module
    
CONTAINS

    FUNCTION weight_gauss_points(root_legendre, k_legendre) RESULT(weight)
        USE elementary_polynomials_module, ONLY: D_legendre
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: root_legendre
        INTEGER,      INTENT(IN) :: k_legendre
        REAL(KIND=8) :: weight

        weight = 2 / (1-root_legendre**2)
        weight = weight / (D_legendre(root_legendre, k_legendre))**2 

    END FUNCTION weight_gauss_points


    SUBROUTINE  element_1d_arbitrary_pk (w, d, p, n_ws, l_Gs, k)
        !===One-dimensional element with linear interpolation
        !===and two Gauss integration points
        !===w(n_w, l_G)    : values of shape functions at Gauss points
        !===d(1, n_w, l_G) : derivatives values of shape functions at Gauss points
        !===p(l_G)         : weight for Gaussian quadrature at Gauss points
        USE elementary_polynomials_module, ONLY: zero_legendre, Lagrange_1D, D_Lagrange_1D
        IMPLICIT NONE
        INTEGER,                                INTENT(IN)  :: n_ws, l_Gs, k
        REAL(KIND=8), DIMENSION(   n_ws, l_Gs), INTENT(OUT) :: w
        REAL(KIND=8), DIMENSION(1, n_ws, l_Gs), INTENT(OUT) :: d
        REAL(KIND=8), DIMENSION(l_Gs),          INTENT(OUT) :: p
        REAL(KIND=8), DIMENSION(l_Gs) :: xx
        REAL(KIND=8), DIMENSION(n_ws) :: tab_nodes
        INTEGER :: i, j, n
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: legendre_pts

        tab_nodes(1) = -1
        tab_nodes(2) = 1
        IF (n_ws > 2) THEN
            tab_nodes(3:) = [(-1+2.d0*n/(n_ws-1.d0), n=1, n_ws-2)]
        END IF
        CALL zero_legendre(k+1, legendre_pts)

        DO j = 1, l_Gs
            DO i=1, n_ws
                w(i, j)    = Lagrange_1D(tab_nodes, i, legendre_pts(j))! Lagrange associated to tab_nodes(i), evaluated in legendre(j)
                d(1, i, j) = D_Lagrange_1D(tab_nodes, i, legendre_pts(j))! D_Lagrange associated to tab_nodes(i), evaluated in legendre(j)
            END DO
            p(j) = weight_gauss_points(legendre_pts(j), k+1)
        ENDDO
    END SUBROUTINE  element_1d_arbitrary_pk

END MODULE arbitrary_GP_1D_module