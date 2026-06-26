MODULE setup

    USE create_laplace_solver_ksp_module, ONLY: abstract_laplace_solver_type
    TYPE, EXTENDS(abstract_laplace_solver_type) :: laplace_solver_type
    CONTAINS
        PROCEDURE :: dir_bc => ex_sol
    END TYPE

    PUBLIC :: rho, velocity, momentum

    ! INTEGER, PARAMETER :: kx = -1, ky = 1
    INTEGER, PARAMETER :: kx = 1, ky = 2
    REAL(KIND = 8), PARAMETER :: pi=ACOS(-1.d0)
CONTAINS

  FUNCTION ex_sol(this, rr) RESULT(uu)
    IMPLICIT NONE
    CLASS(laplace_solver_type),    INTENT(INOUT) :: this
    REAL(KIND = 8), DIMENSION(:, :),  INTENT(IN) :: rr
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: uu
    REAL(KIND = 8) :: kxmax=ACOS(-1.d0), kymax=3*ACOS(-1.d0)
    uu = SIN(kxmax * rr(1, :))*SIN(kymax * rr(2, :) + .1d0)
    !=== dummy to avoid warning
    RETURN
    uu = this%viscosity
    !=== dummy
  END FUNCTION ex_sol

  FUNCTION rho(rr, time) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(size(rr,2)) :: vv
    vv = COS(2*pi*kx*(rr(1,:)-time))*COS(2*pi*ky*(rr(2,:)+time)) + 1.1d0
  END FUNCTION rho

  FUNCTION velocity(comp, rr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    INTEGER, INTENT(IN) :: comp
    REAL(KIND=8), DIMENSION(size(rr,2)) :: vv
    IF (comp==1) THEN
      vv(:) = COS(2*pi*kx*rr(1,:))*SIN(2*pi*ky*rr(2,:))
    ELSEIF (comp==2) THEN
      vv(:) = SIN(2*pi*kx*rr(1,:))*COS(2*pi*ky*rr(2,:))
    END IF
  END FUNCTION velocity

  FUNCTION momentum(comp, rr, time) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    INTEGER, INTENT(IN) :: comp
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(size(rr,2)) :: vv
    vv = velocity(comp, rr)*rho(rr, time)
  END FUNCTION momentum

  FUNCTION source(rr,mu,lambda_in) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(size(rr,2),size(rr,1)) :: vv
    REAL(KIND=8) :: mu, lambda_in, lambda
    INTEGER :: k
    lambda = lambda_in - 2.d0/3.d0*mu
    vv(:,1) = (2*pi)**2*((2*kx**2+ky*(kx+ky))*mu + kx*(kx+ky)*lambda)*COS(2*pi*kx*(rr(1,:)))*SIN(2*pi*ky*(rr(2,:)))
    vv(:,2) = (2*pi)**2*((2*ky**2+kx*(kx+ky))*mu + ky*(kx+ky)*lambda)*SIN(2*pi*kx*(rr(1,:)))*COS(2*pi*ky*(rr(2,:)))
  END FUNCTION
END MODULE setup
