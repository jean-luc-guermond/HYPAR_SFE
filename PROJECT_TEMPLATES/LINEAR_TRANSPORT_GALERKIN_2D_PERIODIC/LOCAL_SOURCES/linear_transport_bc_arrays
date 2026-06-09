MODULE linear_transport_bc_arrays
   USE dirichlet_type_module, ONLY : dirichlet_bc
   USE def_type_mesh,         ONLY : mesh_type, petsc_csr_LA

   TYPE :: linear_transport_bc_type
      TYPE(dirichlet_bc) :: rho_bc
      PROCEDURE(template_scal_anal), POINTER :: rho_anal   => NULL()
   CONTAINS
      PROCEDURE :: sol_anal          => sol_anal_linear_transport
      PROCEDURE :: initial_condition => init_anal
   END TYPE linear_transport_bc_type

   ABSTRACT INTERFACE
      FUNCTION template_scal_anal(this, time, rr) RESULT(vv)
         IMPORT :: linear_transport_bc_type
         IMPLICIT NONE
         CLASS(linear_transport_bc_type),            INTENT(INOUT) :: this
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
         REAL(KIND = 8),                  INTENT(IN) :: time
         REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      END FUNCTION template_scal_anal

      FUNCTION template_vect_anal(this, comp, time, rr) RESULT(vv)
         IMPORT :: linear_transport_bc_type
         IMPLICIT NONE
         CLASS(linear_transport_bc_type),            INTENT(INOUT) :: this
         INTEGER,                         INTENT(IN) :: comp
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
         REAL(KIND = 8),                  INTENT(IN) :: time
         REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      END FUNCTION template_vect_anal
   END INTERFACE

   REAL(KIND=8), PRIVATE :: r_dummy
   INTEGER, PRIVATE :: int_dummy

CONTAINS

!======== GENERIC SOL_ANAL AND INIT BASED ON STATE FUNCTIONS

   FUNCTION sol_anal_linear_transport(this, comp, time, rr) RESULT(vv)
      USE space_dim, ONLY: k_dim
      USE my_util, ONLY: to_str, error_petsc
      IMPLICIT NONE
      CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv

      IF (comp /= 1) CALL error_petsc("BUG in sol_anal => linear transport is a scalar pb, so comp="//to_str(comp)//" is wrong")
      vv = this%rho_anal(time, rr)

   END FUNCTION sol_anal_linear_transport

   SUBROUTINE init_anal(this, un, time, rr)
      USE my_util,   ONLY: error_petsc, to_str
      USE space_dim, ONLY: k_dim
      IMPLICIT NONE
      CLASS(linear_transport_bc_type) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), 1), INTENT(OUT) :: un
      REAL(KIND = 8), INTENT(IN) :: time
      INTEGER :: comp

      DO comp=1, SIZE(un, 2)
         un(:, comp) = this%sol_anal(comp, time, rr)
      END DO
   END SUBROUTINE init_anal

!======== EXAMPLES OF STATE FUNCTIONS FOR linear_transport OBJECT

   ! FUNCTION sine_initialize(this, rr) RESULT(vv)
   !    IMPLICIT NONE
   !    CLASS(linear_transport_bc_type),            INTENT(INOUT) :: this
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN)    :: rr
   !    REAL(KIND = 8),                  INTENT(IN)    :: time
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2))         :: vv
   !    REAL(KIND = 8) :: x0, x1
   !    INTEGER :: n
   !    x0 = 0.1d0
   !    x1 = 0.3d0
   !    DO n=1, SIZE(rr,2)
   !       IF (x0 < rr(n) .AND. rr(n) < x1) THEN

   !       END IF
   !    END DO
   !    vv = 1.d0
   !    RETURN
   !    !===dummy to avoid warning in compilation===!
   !    r_dummy = time; r_dummy = this%gamma 
   !    !===dummy to avoid warning in compilation===!
   ! END FUNCTION sine_initialize

   FUNCTION scal_one(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(linear_transport_bc_type),            INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN)    :: rr
      REAL(KIND = 8),                  INTENT(IN)    :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2))         :: vv
      vv = 1.d0
      RETURN
      !===dummy to avoid warning in compilation===!
      r_dummy = time
      !===dummy to avoid warning in compilation===!
   END FUNCTION scal_one

   FUNCTION vect_one(this, comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
      INTEGER,              INTENT(IN)    :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8),                  INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = 1.d0
      RETURN
      !===dummy to avoid warning in compilation===!
      r_dummy = time; int_dummy = comp
      !===dummy to avoid warning in compilation===!
   END FUNCTION vect_one

END MODULE linear_transport_bc_arrays
