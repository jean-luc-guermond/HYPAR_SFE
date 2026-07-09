MODULE setup
   USE space_dim, ONLY : k_dim
   USE burgers_bc_arrays, ONLY: burgers_bc_type
   USE burgers_type_module, ONLY: burgers_type
   USE burgers_eta_commute, ONLY: default_eta_commute

   PUBLIC :: init_state_functions

   PRIVATE
   REAL(KIND=8), PARAMETER :: x0=0.5d0, length = 1.d0

   INTEGER, PRIVATE, PARAMETER :: INIT_SINE_RHO=1, INIT_STEP_RHO=2
   CHARACTER(LEN=20), DIMENSION(2), PRIVATE, PARAMETER  :: list_init_rho = &
               [CHARACTER(LEN=20) :: 'sine', 'step']
CONTAINS

!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================
   
   SUBROUTINE init_state_functions(burgers, init_case)
      USE my_util, ONLY: error_petsc, to_str
      IMPLICIT NONE
      CLASS(burgers_type), INTENT(INOUT) :: burgers
      INTEGER, INTENT(IN) :: init_case
      
      SELECT CASE(init_case)
      CASE(INIT_SINE_RHO)
         burgers%bc%rho_anal   => sine_rho_anal
      CASE(INIT_STEP_RHO)
         burgers%bc%rho_anal   => step_rho_anal
      CASE DEFAULT
         CALL error_petsc("BUG in init_state_functions => wrong init case "//to_str(init_case))
      END SELECT
      burgers%eta_commute   => default_eta_commute
      
   END SUBROUTINE init_state_functions

   FUNCTION step_rho_anal(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(burgers_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8),                  INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      INTEGER :: n

      DO n = 1, SIZE(rr,2)
         IF (rr(1, n)< Length/2.d0) THEN
            vv(n) = rr(1, n)/(1.d0+time)
         ELSE
            vv(n) = (rr(1, n)-Length)/(1.d0+time)
         END IF
      END DO
      !=== Dummy to avoid warning
      RETURN
      vv = SUM(this%rho_bc%jsd)*1.d0
      !=== Dummy to avoid warning
   END FUNCTION step_rho_anal

   FUNCTION sine_rho_anal(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(burgers_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8),                  INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      REAL(KIND = 8),                   PARAMETER :: pi=ACOS(-1.d0)
      INTEGER                                     :: n

      vv = SIN(2*pi*rr(1,:)/Length)

      !=== Dummy to avoid warning
      RETURN
      vv = SUM(this%rho_bc%jsd)*1.d0
      !=== Dummy to avoid warning
   END FUNCTION sine_rho_anal

END MODULE setup
