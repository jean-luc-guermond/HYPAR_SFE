PROGRAM test_lagrange
    USE elementary_polynomials_module, ONLY: Lagrange_1D, D_Lagrange_1D
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tab_xj
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fx
    REAL(KIND=8) :: dx, epsilon=1.d-11, epsilon_loose=1.d-5, numerator, denominator, result, err_relat
    INTEGER :: i, j, n, l
    LOGICAL :: test

    n = 10
    dx = 1.d0/(n-1)

    ALLOCATE(fx(n), tab_xj(n))
    tab_xj = [(dx*i, i=0,n-1)]

    !====== TEST OF LAGRANGE 1d ======!
    test = .TRUE.
    out_do_1: DO i=1, n
        fx = 0.d0
        DO j=1, n
            fx(j) = Lagrange_1D(tab_xj, i, tab_xj(j))
            IF ((i==j) .AND. (ABS(fx(j)-1) > epsilon)) THEN
                test = .FALSE.
                EXIT out_do_1
            ELSE IF ((i/=j) .AND. (ABS(fx(j)) > epsilon)) THEN
                test = .FALSE.
                EXIT out_do_1
            END IF
        END DO
    END DO out_do_1

    IF (.NOT. test) THEN
        write(*,*) 'BUG in Lagrange_1D => wrong value at special point (i,j)=', i, j, fx(j)
        stop
    ELSE
        write(*,*) 'Success for Lagrange_1D'
    END IF

    !====== TEST 1 OF D_LAGRANGE 1d ======!
    test = .TRUE.
    out_do_2: DO i=1, n
        fx = 0.d0
        DO j=1, n
            fx(j) = D_Lagrange_1D(tab_xj, i, tab_xj(j))
            !=== compute analytical special value
            IF ((i/=j)) THEN
                numerator = 1.d0
                denominator = 1.d0
                DO l=1, n
                    IF (l==i) CYCLE
                    denominator = denominator * (tab_xj(i) - tab_xj(l))
                    IF (l==j) CYCLE
                    numerator = numerator * (tab_xj(j) - tab_xj(l)) 
                END DO
                result = numerator / denominator
                err_relat = ABS(fx(j)-result)/ABS(result)
                IF (err_relat>epsilon) THEN
                    test = .FALSE.
                    EXIT out_do_2
                END IF
            
            ELSE IF ((i==j)) THEN
                result = 0.d0
                DO l=1, n
                    IF (l==i) CYCLE
                    result = result + 1/(tab_xj(i) - tab_xj(l))
                END DO
                err_relat = ABS(fx(j)-result)/ABS(result)
                IF (err_relat>epsilon) THEN
                    test = .FALSE.
                    EXIT out_do_2
                END IF
            END IF
        END DO
    END DO out_do_2

    IF (.NOT. test) THEN
        write(*,*) 'BUG in D_Lagrange_1D test 1 => wrong value at special point (i,j)=', i, j, fx(j), result, err_relat
        stop
    ELSE
        write(*,*) 'Success for D_Lagrange_1D test 1'
    END IF

    !====== TEST 2 OF D_LAGRANGE 1d ======!
    test = .TRUE.
    out_do_3: DO i=1, n
        fx = 0.d0
        DO j=1, n
            fx(j) = D_Lagrange_1D(tab_xj, i, tab_xj(j))
            !=== compute analytical special value
            result = (Lagrange_1D(tab_xj, i, tab_xj(j)+epsilon) - Lagrange_1D(tab_xj, i, tab_xj(j)-epsilon)) / (2*epsilon)
            err_relat = ABS(fx(j)-result)/ABS(result)
            IF (err_relat>epsilon_loose) THEN
                test = .FALSE.
                EXIT out_do_3
            END IF
            
        END DO
    END DO out_do_3

    IF (.NOT. test) THEN
        write(*,*) 'BUG in D_Lagrange_1D test 2 => wrong value at special point (i,j)=', i, j, fx(j), result, err_relat
        stop
    ELSE
        write(*,*) 'Success for D_Lagrange_1D test 2'
    END IF

    WRITE(*,*) "Overall success for Lagrange 1D", "1234567891"

END PROGRAM test_lagrange