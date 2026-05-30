MODULE elementary_polynomials_module

CONTAINS

!========================================================
!=========== LAGRANGE POLYNOMIALS =======================
!========================================================

    FUNCTION Lagrange_1D(tab_xj, i, x) RESULT(fx)
        !> Computes the Lagrange polynomials associated to tab_xj
        !! equal to 1 if evaluated on tab_xj(i), equal to zero if on any other tab_xj(j)
        !! input x
        !! output fx
        IMPLICIT NONE
        INTEGER,                    INTENT(IN) :: i
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: tab_xj
        REAL(KIND=8),               INTENT(IN) :: x             
        REAL(KIND=8)                           :: fx

        INTEGER :: k, j

        k = SIZE(tab_xj)
        fx = 1.d0
        DO j=1, k
            IF (j==i) CYCLE
            fx = fx * (x-tab_xj(j))/(tab_xj(i)-tab_xj(j))
        END DO

    END FUNCTION Lagrange_1D

    FUNCTION D_Lagrange_1D(tab_xj, i, x) RESULT(fx)
        !> Computes the derivative of the Lagrange polynomials associated to tab_xj
        !! input x
        !! output fx
        IMPLICIT NONE
        INTEGER,                    INTENT(IN) :: i
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: tab_xj
        REAL(KIND=8),               INTENT(IN) :: x             
        REAL(KIND=8)                           :: fx
        REAL(KIND=8)                           :: lth_term
        INTEGER :: k, j, l

        k = SIZE(tab_xj)
        fx = 0.d0
        DO j=1,k
            IF (j==i) CYCLE
            lth_term = 1.d0/(tab_xj(i)-tab_xj(j))
            DO l=1, k
                IF (l==i) CYCLE
                IF (l==j) CYCLE
                lth_term = lth_term * (x-tab_xj(l))/(tab_xj(i)-tab_xj(l))
            END DO
            fx = fx + lth_term
        END DO

    END FUNCTION D_Lagrange_1D

!========================================================
!=========== LEGENDRE POLYNOMIALS =======================
!========================================================

    FUNCTION legendre(x, m) RESULT(fx)
        !> legendre polynomial of order m
        !! Computation uses the recursive formula
        USE my_util, ONLY: error_petsc, to_str
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: x
        INTEGER,      INTENT(IN) :: m

        INTEGER                  :: k
        REAL(KIND=8)             :: term_k1, term_k2, fx, f_k1, f_k2

        SELECT CASE(m)
        CASE(0)
            fx = 1.d0
        CASE(1)
            fx = x
        CASE(2:)
            f_k2 = 1.d0
            f_k1 = x
            DO k=2, m
                term_k1 = (2.d0*k - 1.d0) * x * f_k1
                term_k2 = (1.d0*k - 1.d0) * f_k2
                fx = 1.d0/k * (term_k1 - term_k2)

                f_k2 = f_k1
                f_k1 = fx
            END DO
        CASE DEFAULT
            CALL error_petsc('BUG in legendre => index should m >= 0, not '//to_str(m))
        END SELECT
    END FUNCTION legendre

    FUNCTION D_legendre(x, m) RESULT(fx)
        !> derivative of legendre polynomial of order m
        !!! WARNING: not accurate enough for m >= 100
        USE my_util, ONLY: error_petsc, to_str
        IMPLICIT NONE
        REAL(KIND=8), INTENT(IN) :: x
        INTEGER,      INTENT(IN) :: m

        REAL(KIND=8)             :: f_k1, f_k2, f_p_k1, f_p_k2, fx, term_k2, term_k1
        INTEGER :: k


        SELECT CASE(m)
        CASE(0)
            fx = 0.d0
        CASE(1)
            fx = 1.d0
        CASE(2:)
            f_p_k1 = 1.d0
            f_p_k2 = 0.d0
            f_k1 = x            
            f_k2 = 1.d0
            DO k=2, m
                term_k1 = (2.d0*k - 1.d0)*(f_k1 + x*f_p_k1)
                term_k2 = (k - 1.d0)*f_p_k2
                fx = 1.d0/k * (term_k1 - term_k2)
                f_p_k2 = f_p_k1
                f_p_k1 = fx

                term_k1 = (2.d0*k - 1.d0) * x * f_k1
                term_k2 = (1.d0*k - 1.d0) * f_k2
                f_k2 = f_k1
                f_k1 = 1.d0/k * (term_k1 - term_k2)
            END DO

        CASE DEFAULT
            CALL error_petsc('BUG in D_legendre => index should m >= 0, not '//to_str(m))
        END SELECT
    END FUNCTION D_legendre

    SUBROUTINE zero_legendre(k, zeros)
        USE my_util,       ONLY: error_petsc, to_str
        USE Newton_method, ONLY: newton
         !> Computes the zeros of the legendre polynomial of order k
         !! output zeros in increasing order in the array zeros
         !! input k, number of zeros to compute (i.e. order of legendre polynomial)
         !! output zeros, array of size k containing the zeros in increasing order
         !! uses a Newton method to find the zeros, with a good initial guess based on the asymptotic distribution of the zeros of legendre polynomials
         !! for k <= 5, the zeros are hard-coded, for k > 5, the zeros are computed using a Newton method with a good initial guess based on the asymptotic distribution of the zeros of legendre pol
        IMPLICIT NONE
        INTEGER,                                 INTENT(IN) :: k
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE,INTENT(OUT) :: zeros
        REAL(KIND=8) :: delta, x0, x1, xmid, val0, val1
        REAL(KIND=8) :: pi=ACOS(-1.d0), tol=1.d-13
        INTEGER      :: i, l


        ALLOCATE(zeros(k))
        SELECT CASE(k)
        CASE(1)
            zeros(1) = 0.d0
        CASE(2)
            zeros(1) = -1.d0/SQRT(3.d0)
            zeros(2) = +1.d0/SQRT(3.d0)
        CASE(3)
            zeros(1) = -SQRT(3.d0/5.d0)
            zeros(2) = 0.d0
            zeros(3) = +SQRT(3.d0/5.d0)
        CASE(4)
            zeros(1) =  -SQRT((15.d0+2*SQRT(30.d0))/35.d0)
            zeros(2) =  -SQRT((15.d0-2*SQRT(30.d0))/35.d0)
            zeros(3) =  +SQRT((15.d0-2*SQRT(30.d0))/35.d0)
            zeros(4) =  +SQRT((15.d0+2*SQRT(30.d0))/35.d0)
        CASE(5)
            delta = 70.d0**2 - 4.d0*15.d0*63.d0
            zeros(1) = -SQRT((70.d0 + SQRT(delta)) / (2*63.d0))
            zeros(2) = -SQRT((70.d0 - SQRT(delta)) / (2*63.d0))
            zeros(3) = 0.d0
            zeros(4) = +SQRT((70.d0 - SQRT(delta)) / (2*63.d0))
            zeros(5) = +SQRT((70.d0 + SQRT(delta)) / (2*63.d0))
        CASE(6:)
            !=== Init array, assign by the way the value zero, necessary if k is odd
            zeros = 0.d0 
            DO i=1, k/2 !compute roots in [1,0)
                !=== Init knowing the intervals of all zeros of legendre polynomials
                l = (k - 1) - i + 1
                x0 = COS(pi * (l+0.5d0) / (k+0.5d0))
                x1 = COS(pi * (l+1.d0) / (k+0.5d0))
                val0 = legendre(x0, k)
                val1 = legendre(x1, k)
                IF (val0 * val1 > 0) THEN
                    CALL error_petsc("BUG in zero_legendre for m = "//to_str(k)//": wrong init for finding zeros")
                END IF
                !=== Newton loop
                xmid = (x0+x1)/2.d0
                CALL newton(xmid, legendre_k, Dlegendre_k, tol)
                zeros(i) = xmid
                !=== Check the roots are in increasing order
                IF (i > 1) THEN
                    IF (zeros(i) < zeros(i-1)) THEN
                        WRITE(*,*) zeros
                        CALL error_petsc("BUG in zero_legendre => zeros are not in increasing order")
                    END IF
                END IF
            !=== Set roots in the other half-space (0,1]
                zeros(k-i+1) = -zeros(i)
            END DO
        CASE DEFAULT
            CALL error_petsc('BUG in zero_legendre => case k = '//to_str(k)//' not coded yet, only 1/2/3')
        END SELECT

    CONTAINS
        FUNCTION legendre_k(x) RESULT(fx)
            REAL(KIND=8), INTENT(IN) :: x
            REAL(KIND=8)             :: fx
            fx = legendre(x, k)
        END FUNCTION legendre_k
        FUNCTION Dlegendre_k(x) RESULT(fx)
            REAL(KIND=8), INTENT(IN) :: x
            REAL(KIND=8)             :: fx
            fx = D_legendre(x, k)
        END FUNCTION Dlegendre_k

    END SUBROUTINE zero_legendre

END MODULE elementary_polynomials_module