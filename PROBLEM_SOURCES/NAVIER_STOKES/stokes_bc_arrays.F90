MODULE stokes_bc_arrays
   USE dirichlet_type_module, ONLY : dirichlet_bc
   USE def_type_mesh,         ONLY : mesh_type
   USE space_dim

   TYPE :: stokes_bc_type
      !> Type containing:
      !! Dirichlet BCs
      !! vit_anal & temp_anal => impose Dirichlet BC
      !! rho_imposed => Only used inside stokes%update, i.e only relevant if
      !!                imposing a given rho in a problem where rho does not come from another pb (say e.g only solve Stokes without Euler update)
      TYPE(dirichlet_bc) :: temp, vel(k_dim)
      PROCEDURE(template_vect_anal), NOPASS, POINTER :: vit_anal    => NULL()
      PROCEDURE(template_scal_anal), NOPASS, POINTER :: temp_anal   => NULL()
      PROCEDURE(template_scal_anal), NOPASS, POINTER :: rho_imposed => NULL()
   CONTAINS
   !=== FIXME ===! (to remove?)
      ! PROCEDURE :: sol_anal          => sol_anal_stokes
      ! PROCEDURE :: initial_condition => init_anal
   !=== FIXME ===!
   END TYPE stokes_bc_type

   ABSTRACT INTERFACE
      FUNCTION template_scal_anal(time, rr) RESULT(vv)
         IMPLICIT NONE
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
         REAL(KIND = 8),                  INTENT(IN) :: time
         REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      END FUNCTION template_scal_anal

      FUNCTION template_vect_anal(comp, time, rr) RESULT(vv)
         IMPLICIT NONE
         INTEGER,                         INTENT(IN) :: comp
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
         REAL(KIND = 8),                  INTENT(IN) :: time
         REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      END FUNCTION template_vect_anal

   END INTERFACE

   REAL(KIND=8), PRIVATE :: r_dummy
   INTEGER, PRIVATE :: int_dummy

CONTAINS

END MODULE stokes_bc_arrays
