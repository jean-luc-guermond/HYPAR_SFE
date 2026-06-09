MODULE setup

    USE create_laplace_solver_ksp_module, ONLY: abstract_laplace_solver_type
    TYPE, EXTENDS(abstract_laplace_solver_type) :: laplace_solver_type
    CONTAINS
        PROCEDURE :: dir_bc => ex_sol
    END TYPE

CONTAINS

  FUNCTION source(rr) RESULT(uu)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:, :) :: rr
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: uu
    REAL(KIND = 8) :: kmax=4*ACOS(-1.d0)
    uu = (1+kmax**2)*COS(kmax * rr(1, :) +  .7d0)
  END FUNCTION source

  FUNCTION ex_sol(this, rr) RESULT(uu)
    IMPLICIT NONE
    CLASS(laplace_solver_type),    INTENT(INOUT) :: this
    REAL(KIND = 8), DIMENSION(:, :),  INTENT(IN) :: rr
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: uu
    REAL(KIND = 8) :: kmax=4*ACOS(-1.d0)
    uu = COS(kmax * rr(1, :) +  .7d0 )
    !=== dummy to avoid warning
    RETURN
    uu = this%viscosity
    !=== dummy
  END FUNCTION ex_sol

END MODULE setup
