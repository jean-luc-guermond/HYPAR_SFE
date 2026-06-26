MODULE stokes_bc_arrays
   USE dirichlet_type_module, ONLY : dirichlet_bc
   USE def_type_mesh,         ONLY : mesh_type, petsc_csr_LA
   USE space_dim

   TYPE :: stokes_bc_type
      TYPE(dirichlet_bc) :: temp, vel(k_dim)
      PROCEDURE(template_vect_anal), POINTER :: vit_anal   => NULL()
      PROCEDURE(template_scal_anal), POINTER :: temp_anal  => NULL()
   CONTAINS
      ! PROCEDURE :: sol_anal          => sol_anal_stokes
      ! PROCEDURE :: initial_condition => init_anal
   END TYPE stokes_bc_type

   ABSTRACT INTERFACE
      FUNCTION template_scal_anal(this, time, rr) RESULT(vv)
         IMPORT :: stokes_bc_type
         IMPLICIT NONE
         CLASS(stokes_bc_type),            INTENT(INOUT) :: this
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
         REAL(KIND = 8),                  INTENT(IN) :: time
         REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      END FUNCTION template_scal_anal

      FUNCTION template_vect_anal(this, comp, time, rr) RESULT(vv)
         IMPORT :: stokes_bc_type
         IMPLICIT NONE
         CLASS(stokes_bc_type),            INTENT(INOUT) :: this
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

   ! FUNCTION sol_anal_stokes(this, comp, time, rr) RESULT(vv)
   !    USE space_dim, ONLY: k_dim
   !    IMPLICIT NONE
   !    CLASS(stokes_bc_type), INTENT(INOUT) :: this
   !    INTEGER, INTENT(IN) :: comp
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
   !    REAL(KIND = 8), INTENT(IN) :: time
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
   !    SELECT CASE(comp)
   !    CASE(1)
   !       vv = this%rho_anal(time, rr)
   !    CASE(2:k_dim+1)
   !       vv = this%mt_anal(comp-1, time, rr)
   !    CASE(k_dim+2)
   !       vv = this%E_anal(time, rr)
   !    CASE DEFAULT
   !       WRITE(*, *) ' BUG in sol_anal, comp=', comp, 'should be <=', k_dim+2
   !       STOP
   !    END SELECT
   ! END FUNCTION sol_anal_stokes

   ! SUBROUTINE init_anal(this, un, time, rr)
   !    USE my_util,   ONLY: error_petsc, to_str
   !    USE space_dim, ONLY: k_dim
   !    IMPLICIT NONE
   !    CLASS(stokes_bc_type) :: this
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2), k_dim + 2), INTENT(OUT) :: un
   !    REAL(KIND = 8), INTENT(IN) :: time
   !    INTEGER :: comp

   !    DO comp=1, SIZE(un, 2)
   !       un(:, comp) = this%sol_anal(comp, time, rr)
   !    END DO
   ! END SUBROUTINE init_anal

END MODULE stokes_bc_arrays
