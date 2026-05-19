MODULE euler_bc_arrays
   USE dirichlet_type_module, ONLY : dirichlet_bc
   USE def_type_mesh,         ONLY : mesh_type, petsc_csr_LA

   TYPE :: euler_bc_type
      TYPE(dirichlet_bc) :: rho_bc, u_bc(2), whole_bdy_bc, udotn_bc
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: udotn_normal_vtx
      REAL(KIND=8) :: gamma
      PROCEDURE(template_vect_anal), POINTER :: mt_anal    => NULL()
      PROCEDURE(template_vect_anal), POINTER :: vit_anal   => NULL()
      PROCEDURE(template_scal_anal), POINTER :: rho_anal   => NULL()
      PROCEDURE(template_scal_anal), POINTER :: press_anal => NULL()
      PROCEDURE(template_scal_anal), POINTER :: E_anal     => NULL()
   CONTAINS
      PROCEDURE :: sol_anal          => sol_anal_euler
      PROCEDURE :: initial_condition => init_anal
   END TYPE euler_bc_type

   ABSTRACT INTERFACE
      FUNCTION template_scal_anal(this, time, rr) RESULT(vv)
         IMPORT :: euler_bc_type
         IMPLICIT NONE
         CLASS(euler_bc_type),            INTENT(INOUT) :: this
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
         REAL(KIND = 8),                  INTENT(IN) :: time
         REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      END FUNCTION template_scal_anal

      FUNCTION template_vect_anal(this, comp, time, rr) RESULT(vv)
         IMPORT :: euler_bc_type
         IMPLICIT NONE
         CLASS(euler_bc_type),            INTENT(INOUT) :: this
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

   FUNCTION sol_anal_euler(this, comp, time, rr) RESULT(vv)
      USE space_dim, ONLY: k_dim
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      SELECT CASE(comp)
      CASE(1)
         vv = this%rho_anal(time, rr)
      CASE(2:k_dim+1)
         vv = this%mt_anal(comp-1, time, rr)
      CASE(k_dim+2)
         vv = this%E_anal(time, rr)
      CASE DEFAULT
         WRITE(*, *) ' BUG in sol_anal, comp=', comp, 'should be <=', k_dim+2
         STOP
      END SELECT
   END FUNCTION sol_anal_euler

   SUBROUTINE init_anal(this, un, time, rr)
      USE my_util,   ONLY: error_petsc, to_str
      USE space_dim, ONLY: k_dim
      IMPLICIT NONE
      CLASS(euler_bc_type) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), k_dim + 2), INTENT(OUT) :: un
      REAL(KIND = 8), INTENT(IN) :: time
      INTEGER :: comp

      DO comp=1, SIZE(un, 2)
         un(:, comp) = this%sol_anal(comp, time, rr)
      END DO
   END SUBROUTINE init_anal

!======== EXAMPLES OF STATE FUNCTIONS FOR EULER OBJECT

   FUNCTION E_anal_ideal_gas(this, time, rr) RESULT(vv)
      USE space_dim
      IMPLICIT NONE
      CLASS(euler_bc_type),            INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN)    :: rr
      REAL(KIND = 8),                  INTENT(IN)    :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2))         :: vv, v_sqr
      INTEGER :: k

      v_sqr = 0.d0
      DO k=1, k_dim
         v_sqr = v_sqr + (this%vit_anal(k, time,  rr)**2)
      END DO

      vv = this%press_anal(time, rr) / (this%gamma - 1.d0) &
           + this%rho_anal(time, rr) * v_sqr / 2
   END FUNCTION E_anal_ideal_gas

   FUNCTION mt_anal_rho_times_vit(this, comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = this%rho_anal(time, rr) * this%vit_anal(comp, time, rr)
   END FUNCTION mt_anal_rho_times_vit


   FUNCTION scal_one(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type),            INTENT(INOUT) :: this
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
      CLASS(euler_bc_type), INTENT(INOUT) :: this
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

END MODULE euler_bc_arrays
