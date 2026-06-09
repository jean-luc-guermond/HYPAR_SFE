MODULE setup
   USE space_dim, ONLY : k_dim
   USE linear_transport_bc_arrays, ONLY: linear_transport_bc_type
   USE linear_transport_type_module, ONLY: linear_transport_type
   USE linear_transport_eta_commute, ONLY: default_eta_commute

   PUBLIC :: init_state_functions, sbr_linear_transport!cst_linear_transport

   TYPE, EXTENDS(linear_transport_type) :: sbr_linear_transport
   CONTAINS
      PROCEDURE :: transport => sbr_transport 
   END TYPE sbr_linear_transport

   PRIVATE
   REAL(KIND=8), PARAMETER :: x0=0.0d0, y0=0.d0, a = 0.3d0
CONTAINS

!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================
   
   SUBROUTINE init_state_functions(linear_transport)
      IMPLICIT NONE
      CLASS(linear_transport_type), INTENT(INOUT) :: linear_transport

      linear_transport%bc%rho_anal => rho_paper_circle
      ! linear_transport%bc%rho_anal   => sine_rho_anal
      ! linear_transport%bc%rho_anal   => step_rho_anal
      linear_transport%eta_commute   => default_eta_commute
      
   END SUBROUTINE init_state_functions

   FUNCTION rho_paper_circle(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8),                  INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      INTEGER :: n, k_x, k_y
      REAL(KIND=8) :: xt, yt, pi=ACOS(-1.d0), length = 2.d0

         ! k = floor((rr(1, n) - time)/length)
         ! x = rr(1, n) - time -k*length  

      DO n=1, SIZE(rr,2)

         ! xt = x0 + 1.d0*time + 1.d0
         ! yt = y0 + 2.d0*time + 1.d0

         k_x = floor((rr(1, n) + 1.d0 - 1.d0*time)/length)
         k_y = floor((rr(2, n) + 1.d0 - 2.d0*time)/length)
         ! k_x = FLOOR(xt/length)
         ! k_y = FLOOR(yt/length)

         xt = rr(1, n) - 1.d0*time -k_x*length
         yt = rr(2, n) - 2.d0*time -k_y*length

         vv(n) = 0.5d0*(1.d0-tanh(((xt)**2+(yt)**2)/a**2 - 1.d0))
      END DO

   END FUNCTION rho_paper_circle

   ! FUNCTION step_rho_anal(this, time, rr) RESULT(vv)
   !    IMPLICIT NONE
   !    CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
   !    REAL(KIND = 8),                  INTENT(IN) :: time
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
   !    INTEGER :: n

   !    DO n=1, SIZE(rr,2)
   !       IF (rr(1, n)-cx*time < x0) THEN
   !          vv(n) = rhol
   !       ELSE
   !          vv(n) = rhor
   !       END IF
   !    END DO

   ! END FUNCTION step_rho_anal

   ! FUNCTION sine_rho_anal(this, time, rr) RESULT(vv)
   !    IMPLICIT NONE
   !    CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
   !    REAL(KIND = 8),                  INTENT(IN) :: time
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
   !    REAL(KIND = 8),                   PARAMETER :: pi=ACOS(-1.d0)
   !    INTEGER                                     :: n

   !    DO n=1, SIZE(rr,2)
   !       vv(n) = SIN((rr(1,n)-cx*time)*2*pi) + 2.d0
   !    END DO

   ! END FUNCTION sine_rho_anal


   ! FUNCTION cst_transport(this, rr, comp) RESULT(vv)
   !    USE my_util, ONLY: error_petsc, to_str
   !    IMPLICIT NONE
   !    CLASS(cst_linear_transport)                  :: this
   !    REAL(KIND = 8), DIMENSION(k_dim), INTENT(IN) :: rr
   !    INTEGER,                          INTENT(IN) :: comp
   !    REAL(KIND = 8)                               :: vv

   !    SELECT CASE(comp)
   !    CASE(1)
   !       vv = cx
   !    CASE(2)
   !       vv = cy
   !    CASE DEFAULT
   !       CALL error_petsc("BUG in cst_transport => wrong comp = "//to_str(comp))
   !    END SELECT
   ! END FUNCTION cst_transport

   FUNCTION sbr_transport(this, rr, comp) RESULT(vv)
      USE my_util, ONLY: error_petsc, to_str
      USE space_dim
      IMPLICIT NONE
      CLASS(sbr_linear_transport)                  :: this
      REAL(KIND = 8), DIMENSION(k_dim), INTENT(IN) :: rr
      INTEGER,                          INTENT(IN) :: comp
      REAL(KIND = 8)                               :: vv
      REAL(KIND = 8), PARAMETER :: pi=ACOS(-1.d0)

      SELECT CASE(comp)
      CASE(1)
         vv = 1.d0
      CASE(2)
         vv = 2.d0
      CASE DEFAULT
         CALL error_petsc("BUG in cst_transport => wrong comp = "//to_str(comp))
      END SELECT
   END FUNCTION sbr_transport

END MODULE setup
