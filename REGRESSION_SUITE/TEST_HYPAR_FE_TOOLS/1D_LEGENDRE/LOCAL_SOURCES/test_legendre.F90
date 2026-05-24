PROGRAM test_lagrange
    USE elementary_polynomials_module, ONLY: legendre, zero_legendre, D_legendre
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tab_xj
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fx
    REAL(KIND=8) :: dx, epsilon=1.d-9, epsilon_loose=1.d-5, numerator, denominator, result, err_relat, ap, am, a
    INTEGER :: i, j, n, l, m, N_max
    LOGICAL :: test

    n = 50

    !====== TEST OF Legendre vs zero_legendre 1d ======!
    test = .TRUE.
    out_do_1: DO m=1, n
        CALL zero_legendre(m, tab_xj)
        ALLOCATE(fx(size(tab_xj)))
        DO j=1, m
            fx(j) = legendre(tab_xj(j), m)
            IF (ABS(fx(j)) > epsilon) THEN
                test = .FALSE.
                EXIT out_do_1
            END IF
        END DO
        DEALLOCATE(fx, tab_xj)
        write(*,*) 'Success for legendre m = ', m, ' vs zero_legendre '
    END DO out_do_1

    IF (.NOT. test) THEN
        write(*,*) 'BUG in legendre => no match between zeros and legendre for m = ', m, ', zero, ', j, fx(j)
        stop
    END IF

    !===== TEST OF Legendre vs D_legendre

    N_max = 10
    dx = 1.d0/(N_max-1)

    ALLOCATE(fx(N), tab_xj(N_max))
    tab_xj = [(dx*i, i=0,N_max-1)]

    test = .TRUE.
    out_do_2: DO m=1, n
        fx = 0.d0
        DO j=1, N_max
            fx(j) = D_Legendre(tab_xj(j), m)
            !=== estimate derivative with basic function
            result = (Legendre(tab_xj(j)+epsilon_loose, m) - Legendre(tab_xj(j)-epsilon_loose, m)) / (2*epsilon_loose)
            err_relat = ABS(fx(j)-result)/ABS(result)
            IF (err_relat>epsilon_loose) THEN
                test = .FALSE.
                EXIT out_do_2
            END IF
            
        END DO
    END DO out_do_2

    IF (.NOT. test) THEN
        WRITE(*,*) 'BUG in D_legendre test 2 for m = ', m, ' at point x = ', tab_xj(j) 
        write(*,*) 'values = ', fx(j), result, err_relat
        stop
    ELSE
        write(*,*) 'Success for D_legendre vs legendre'
    END IF

    WRITE(*,*) "Overall success for Legendre", "1234567891"

END PROGRAM test_lagrange