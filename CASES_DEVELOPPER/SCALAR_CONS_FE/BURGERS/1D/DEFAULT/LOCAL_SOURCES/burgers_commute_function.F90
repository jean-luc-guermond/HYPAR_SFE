MODULE burgers_eta_commute

CONTAINS

    FUNCTION default_eta_commute(un) RESULT(eta)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:, :), INTENT(IN)    :: un
        REAL(KIND=8), DIMENSION(SIZE(un,1))         :: eta
        eta = un(:,1) + 3.d0
    END FUNCTION default_eta_commute

END MODULE burgers_eta_commute
