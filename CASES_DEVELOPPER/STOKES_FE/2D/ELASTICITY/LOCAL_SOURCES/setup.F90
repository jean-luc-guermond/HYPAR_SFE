MODULE setup

    INTEGER :: test_case

    INTEGER,        PARAMETER :: kx = 1, ky = 2
    REAL(KIND = 8), PARAMETER :: pi=ACOS(-1.d0), cv=2.5d0
    REAL(KIND = 8) :: mu_viscosity, lambda_viscosity, exp_decay
    PUBLIC :: init_state_functions, Etot !<=== FIXME Etot
    PRIVATE

CONTAINS


!========================================
!==== BUILDING POINTERS =================
!========================================

  SUBROUTINE init_state_functions(stokes, test_case_in)
    USE stokes_parabolic_module
    IMPLICIT NONE
    TYPE(stokes_parabolic_type), INTENT(INOUT) :: stokes 
    INTEGER, INTENT(IN) :: test_case_in

    mu_viscosity     = stokes%mu_viscosity
    lambda_viscosity = stokes%lambda_viscosity

    stokes%bc%vit_anal  => vit_anal_stokes_2D
    stokes%bc%temp_anal => temp_anal_stokes
    stokes%forcing      => stokes_forcing

    !=== Specific to Stokes pb solved without coupling (e.g without Euler)
    stokes%bc%rho_imposed => rho
    !=== Specific to Stokes pb solved without coupling (e.g without Euler)

    test_case = test_case_in

    SELECT CASE(test_case)
    CASE(5,6)
      exp_decay = 0.1d0
    CASE DEFAULT
      exp_decay = 0.d0
    END SELECT

  END SUBROUTINE init_state_functions

  FUNCTION temp_anal_stokes(time, rr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
    REAL(KIND = 8),                  INTENT(IN) :: time
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
    vv = 1.d0
    !=== dummy to avoid warnings
    RETURN
    vv = time
    !=== dummy to avoid warnings
  END FUNCTION temp_anal_stokes

  FUNCTION stokes_forcing(comp, rr, time) RESULT(vv)
    USE stokes_bc_arrays
    USE space_dim 
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), INTENT(IN) :: time
    INTEGER,      INTENT(IN) :: comp
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2))  :: vv
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2), k_dim)  :: vv_complete

    vv_complete = source(rr, mu_viscosity, lambda_viscosity, time)
    vv = vv_complete(:, comp)
  
  END FUNCTION stokes_forcing

!========================================
!==== ANALYTICAL TESTING CASES ==========
!========================================

  FUNCTION rho(time, rr) RESULT(vv)
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(size(rr,2)) :: vv

    SELECT CASE(test_case)
    CASE(1,5)
      vv = 1.d0
    CASE(2)
      vv = 2.d0
    CASE(3,4,6)
      vv = COS(2*pi*kx*(rr(1,:)-time))*COS(2*pi*ky*(rr(2,:)+time)) + 1.1d0
    CASE DEFAULT
      CALL error_petsc("BUG in rho: test_case "//to_str(test_case)//" unknown")
    END SELECT
  END FUNCTION rho

  FUNCTION vit_anal_stokes_2D(comp, time, rr) RESULT(vv)
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    INTEGER, INTENT(IN) :: comp
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(size(rr,2)) :: vv

    SELECT CASE(test_case)
    CASE(5)
          vv = 1.d0*exp(-exp_decay*time)
    CASE DEFAULT
          SELECT CASE(comp)
          CASE(1)
            vv(:) = COS(2*pi*kx*rr(1,:))*SIN(2*pi*ky*rr(2,:))*exp(-exp_decay*time)
          CASE(2)
            vv(:) = SIN(2*pi*kx*rr(1,:))*COS(2*pi*ky*rr(2,:))*exp(-exp_decay*time)
          CASE DEFAULT
            CALL error_petsc("BUG in vit_anal_stokes_2D: wrong comp "//to_str(comp))
          END SELECT  
    END SELECT
  END FUNCTION vit_anal_stokes_2D

  FUNCTION momentum(comp, rr, time) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    INTEGER, INTENT(IN) :: comp
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(size(rr,2)) :: vv
    vv = vit_anal_stokes_2D(comp, time, rr)*rho(time, rr)
  END FUNCTION momentum

  FUNCTION source_elast_2D(rr, mu, lambda_in, time) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(size(rr,2),size(rr,1)) :: vv
    REAL(KIND=8) :: mu, lambda_in, lambda
    lambda = lambda_in - 2.d0/3.d0*mu
    vv(:,1) = (2*pi)**2*((2*kx**2+ky*(kx+ky))*mu + kx*(kx+ky)*lambda)*COS(2*pi*kx*(rr(1,:)))*SIN(2*pi*ky*(rr(2,:)))*exp(-exp_decay*time)
    vv(:,2) = (2*pi)**2*((2*ky**2+kx*(kx+ky))*mu + ky*(kx+ky)*lambda)*SIN(2*pi*kx*(rr(1,:)))*COS(2*pi*ky*(rr(2,:)))*exp(-exp_decay*time)
  END FUNCTION source_elast_2D

  FUNCTION source(rr, mu, lambda_in, time) RESULT(vv)
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(size(rr,2),size(rr,1)) :: vv
    REAL(KIND=8) :: mu, lambda_in
    INTEGER :: k

    vv = source_elast_2D(rr, mu, lambda_in, time)
    DO k=1, k_dim
      vv(:, k) = vv(:, k) + (-exp_decay)*momentum(k, rr, time)
    END DO
  END FUNCTION source

  FUNCTION temperature(rr, time) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(size(rr,2)) :: vv
    vv(:) = 1.d0
    !=== dummy to avoid warnings
    RETURN
    vv = time
    !=== dummy to avoid warnings
  END FUNCTION temperature

  FUNCTION Etot(rr, time) RESULT(vv)
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(size(rr,2)) :: vv
    INTEGER :: k

    vv = 0.d0
    DO k=1, k_dim
      vv = vv + vit_anal_stokes_2D(k, time, rr)**2
    END DO
    vv = 0.5d0*rho(time, rr)*vv**2 + rho(time, rr)*cv*temperature(rr, time)
  END FUNCTION Etot

END MODULE setup
