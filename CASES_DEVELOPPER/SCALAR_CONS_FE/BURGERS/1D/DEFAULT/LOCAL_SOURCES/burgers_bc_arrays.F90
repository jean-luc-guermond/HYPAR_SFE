MODULE burgers_bc_arrays
   USE dirichlet_type_module, ONLY : dirichlet_bc
   USE def_type_mesh,         ONLY : mesh_type

   TYPE :: burgers_bc_type
      TYPE(dirichlet_bc) :: rho_bc
      PROCEDURE(template_scal_anal), POINTER :: rho_anal   => NULL()
   CONTAINS
      PROCEDURE :: sol_anal          => sol_anal_burgers
      PROCEDURE :: initial_condition => init_anal
   END TYPE burgers_bc_type

   ABSTRACT INTERFACE
      FUNCTION template_scal_anal(this, time, rr) RESULT(vv)
         IMPORT :: burgers_bc_type
         IMPLICIT NONE
         CLASS(burgers_bc_type),            INTENT(INOUT) :: this
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
         REAL(KIND = 8),                  INTENT(IN) :: time
         REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      END FUNCTION template_scal_anal

      FUNCTION template_vect_anal(this, comp, time, rr) RESULT(vv)
         IMPORT :: burgers_bc_type
         IMPLICIT NONE
         CLASS(burgers_bc_type),            INTENT(INOUT) :: this
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

   FUNCTION sol_anal_burgers(this, comp, time, rr) RESULT(vv)
      USE space_dim, ONLY: k_dim
      USE my_util, ONLY: to_str, error_petsc
      IMPLICIT NONE
      CLASS(burgers_bc_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv

      IF (comp /= 1) CALL error_petsc("BUG in sol_anal => linear transport is a scalar pb, so comp="//to_str(comp)//" is wrong")
      vv = this%rho_anal(time, rr)

   END FUNCTION sol_anal_burgers

   SUBROUTINE init_anal(this, un, time, rr)
      USE my_util,   ONLY: error_petsc, to_str
      USE space_dim, ONLY: k_dim
      IMPLICIT NONE
      CLASS(burgers_bc_type) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), 1), INTENT(OUT) :: un
      REAL(KIND = 8), INTENT(IN) :: time
      INTEGER :: comp

      DO comp=1, SIZE(un, 2)
         un(:, comp) = this%sol_anal(comp, time, rr)
      END DO
   END SUBROUTINE init_anal

!======== EXAMPLES OF STATE FUNCTIONS FOR burgers OBJECT

   FUNCTION scal_one(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(burgers_bc_type),            INTENT(INOUT) :: this
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
      CLASS(burgers_bc_type), INTENT(INOUT) :: this
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

END MODULE burgers_bc_arrays
