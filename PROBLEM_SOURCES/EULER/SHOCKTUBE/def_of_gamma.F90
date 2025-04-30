MODULE def_of_gamma 
  REAL(KIND=8), PUBLIC :: gamma_ 
CONTAINS
  SUBROUTINE set_gamma_for_riemann_solver(gamma)
    IMPLICIT NONE
    REAL(KIND=8) :: gamma
    gamma_ = gamma
  END SUBROUTINE set_gamma_for_riemann_solver
END MODULE def_of_gamma 
