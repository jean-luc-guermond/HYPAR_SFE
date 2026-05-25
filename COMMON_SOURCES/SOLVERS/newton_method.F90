MODULE Newton_method
    PUBLIC :: newton
    PRIVATE
CONTAINS
    SUBROUTINE newton(T,psi_func,dpsi_func,norm_tol, opt_iter_max, opt_stop_on_neg)
        !> Newton method to solve psi_func(T) = 0 using dpsi_func its derivative.
        !! T is the initial guess and the output
        !! opt_stop_on_neg: set to True to break the code if negative T is reached
        USE my_util, ONLY: pack_opt, error_petsc, to_str
        IMPLICIT NONE
        INTERFACE
            FUNCTION psi_func(T) RESULT(vout)
            IMPLICIT NONE
            REAL(KIND=8), INTENT(IN) :: T
            REAL(KIND=8) :: vout
            END FUNCTION psi_func
            FUNCTION dpsi_func(T) RESULT(vout)
            IMPLICIT NONE
            REAL(KIND=8), INTENT(IN) :: T
            REAL(KIND=8) :: vout
            END FUNCTION dpsi_func
        END INTERFACE
        REAL(KIND=8), INTENT(INOUT) :: T
        REAL(KIND=8)    :: norm_tol
        REAL(KIND=8) :: v_tol, psi, dp, dv, Tin, psio, To
        INTEGER :: iter, itermax, j
        LOGICAL :: stop_on_neg
        INTEGER, INTENT(IN), OPTIONAL :: opt_iter_max
        LOGICAL, INTENT(IN), OPTIONAL :: opt_stop_on_neg

        !=== 100 it by default
        CALL pack_opt(itermax, 100, opt_iter_max)
        !=== By default, do not stop if T becomes negative (if T is a temperature or a density, we do want to stop)
        CALL pack_opt(stop_on_neg, .FALSE., opt_stop_on_neg)

        iter = 0
        psi =  psi_func(T)
        Tin = T
        !=== DEBUGGING ===!
        !write(*,*) 'before Newton loop', T, abs(psi/norm_tol), 'relative tolerance', norm_tol
        !=== DEBUGGING ===!
        DO WHILE (ABS(psi).GE.norm_tol .AND. iter<itermax)
            To = T
            dp = dpsi_func(T)
            dv = psi/dp
            T = T-dv
            iter = iter+1
            !=== DEBUGGING ===!
            !write(*,*) 'in Newton iter, iter, To, T', iter,  To,  T
            !write(*,*) 'in Newton iter, iter, psio, dpo, dvo, psio,dp' , psi, dp, dv, psi/dp
            !=== DEBUGGING ===!
            psi = psi_func(T)
            !=== DEBUGGING ===!
            !write(*,*) 'in Newton iter, iter, psi, psi/tol', psi, ABS(psi)/norm_tol
            !=== DEBUGGING ===!
            IF (iter==itermax) THEN
                CALL error_petsc('Max iter in Newton reached without convergence, &
                T = '//to_str(T)//', Tin = '//to_str(Tin)//', abs(psi) = '//to_str(ABS(psi))//&
                ' norm_tol = '//to_str(norm_tol))
            END IF
            IF ((T<0) .AND. stop_on_neg) THEN
                CALL error_petsc('Negative T in Newton, T = '//to_str(T)//', Tin = '//to_str(Tin)//&
                ', abs(psi) = '//to_str(ABS(psi))//' norm_tol = '//to_str(norm_tol))
            END IF
        END DO
    END SUBROUTINE newton

END MODULE Newton_method