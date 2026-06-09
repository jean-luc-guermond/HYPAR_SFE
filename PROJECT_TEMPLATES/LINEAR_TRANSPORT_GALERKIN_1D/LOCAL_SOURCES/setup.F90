MODULE setup
   USE space_dim, ONLY : k_dim
   USE linear_transport_bc_arrays, ONLY: linear_transport_bc_type
   USE linear_transport_type_module, ONLY: linear_transport_type
   USE linear_transport_eta_commute, ONLY: default_eta_commute

   PUBLIC :: init_state_functions, cst_linear_transport

   TYPE, EXTENDS(linear_transport_type) :: cst_linear_transport
   CONTAINS
      PROCEDURE :: transport => cst_transport 
   END TYPE cst_linear_transport

   PRIVATE
   REAL(KIND=8), PARAMETER :: x0=0.5d0, rhol = 1.d0, rhor = 2.d0, cx = 1.d0, cy = 0.d0
CONTAINS

!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================
   
   SUBROUTINE init_state_functions(linear_transport)
      IMPLICIT NONE
      CLASS(linear_transport_type), INTENT(INOUT) :: linear_transport

      linear_transport%bc%rho_anal   => sine_rho_anal
      ! linear_transport%bc%rho_anal   => step_rho_anal
      linear_transport%eta_commute   => default_eta_commute
      
   END SUBROUTINE init_state_functions

   FUNCTION step_rho_anal(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8),                  INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      INTEGER :: n

      DO n=1, SIZE(rr,2)
         IF (rr(1, n)-cx*time < x0) THEN
            vv(n) = rhol
         ELSE
            vv(n) = rhor
         END IF
      END DO

   END FUNCTION step_rho_anal

   FUNCTION sine_rho_anal(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8),                  INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      REAL(KIND = 8),                   PARAMETER :: pi=ACOS(-1.d0)
      INTEGER                                     :: n

      DO n=1, SIZE(rr,2)
         vv(n) = SIN((rr(1,n)-cx*time)*2*pi) + 2.d0
      END DO

   END FUNCTION sine_rho_anal


   FUNCTION cst_transport(this, rr, comp) RESULT(vv)
      USE my_util, ONLY: error_petsc, to_str
      IMPLICIT NONE
      CLASS(cst_linear_transport)                  :: this
      REAL(KIND = 8), DIMENSION(k_dim), INTENT(IN) :: rr
      INTEGER,                          INTENT(IN) :: comp
      REAL(KIND = 8)                               :: vv

      SELECT CASE(comp)
      CASE(1)
         vv = cx
      CASE(2)
         vv = cy
      CASE DEFAULT
         CALL error_petsc("BUG in cst_transport => wrong comp = "//to_str(comp))
      END SELECT
   END FUNCTION cst_transport

END MODULE setup
