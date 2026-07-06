MODULE linear_transport_type_module
!>> limited global uses to avoid unexpected behaviors
#include "petsc/finclude/petsc.h"
   USE petsc
   USE linear_transport_abstract_hyperbolic_module,           ONLY: hyperbolic_type
   USE linear_transport_bc_arrays,           ONLY: linear_transport_bc_type
   USE Butcher_tableau
   ! USE def_type_mesh,                        ONLY: mesh_type
   ! USE petsc_csr_LA_module,                  ONLY: petsc_csr_LA

   USE read_inputs_module,                   ONLY: rec_length
   USE space_dim,                            ONLY: k_dim
!>> limited global uses to avoid unexpected behaviors

   IMPLICIT NONE

   TYPE argument_linear_transport_type
      CHARACTER(LEN=rec_length) :: c = "=== Linear transport velocity (1D) ==="
   END TYPE argument_linear_transport_type

   TYPE, EXTENDS(hyperbolic_type), ABSTRACT  :: linear_transport_type
      !===Parameters read from data
      REAL(KIND = 8)                          :: c = 1.d0
      TYPE(linear_transport_bc_type)          :: bc
      CHARACTER(len=5), DIMENSION(:), ALLOCATABLE :: name_comp
   CONTAINS
      PROCEDURE, PUBLIC  :: init_linear_transport
      PROCEDURE, PRIVATE :: read_linear_transport_data
      PROCEDURE :: flux           => flux_linear_transport
      PROCEDURE :: compute_lambda => lambda_linear_transport
      PROCEDURE :: construct_bc   => construct_linear_transport_bc
      PROCEDURE :: impose_bc      => impose_bc_linear_transport
      PROCEDURE(function_template_transport), DEFERRED :: transport
   END TYPE linear_transport_type

   ABSTRACT INTERFACE
      FUNCTION function_template_transport(this, rr, comp) RESULT(vv)
         USE space_dim
         IMPORT :: linear_transport_type
         IMPLICIT NONE
         CLASS(linear_transport_type)                 :: this
         REAL(KIND = 8), DIMENSION(k_dim), INTENT(IN) :: rr
         INTEGER,                          INTENT(IN) :: comp
         REAL(KIND = 8)                               :: vv
      END FUNCTION function_template_transport
   END INTERFACE

CONTAINS
   SUBROUTINE init_linear_transport(this, name)
      USE space_dim
      USE my_util, ONLY: error_petsc, to_str
      IMPLICIT NONE
      CLASS(linear_transport_type), INTENT(INOUT) :: this
      CHARACTER(LEN=*),  INTENT(IN)    :: name

      this%name = name
      this%syst_dim = 1
      ALLOCATE(this%name_comp(this%syst_dim))

      this%name_comp(1) = 'rho'

      CALL this%read_linear_transport_data(trim(adjustl(name))//" PARAMETERS")

   END SUBROUTINE init_linear_transport

   SUBROUTINE read_linear_transport_data(this, section_name)
     USE read_inputs_module
     USE space_dim
     IMPLICIT NONE
     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name

     CLASS(linear_transport_type), INTENT(INOUT) :: this
     TYPE(argument_linear_transport_type)        :: argument_data


     !================
     !=== MANDATORY Reading all data file
     !================
     IF (PRESENT(section_name)) THEN
        CALL read_data_init_list(section_name)
     ELSE
        CALL read_data_init_list()
     END IF

     !================
     !=== We now find the relevant information for this specific linear_transport data
     !================

     !===c_transport
     CALL read_data(argument_data%c, this%c, opt_name=this%name, opt_add=k_dim==1)

     !================
     !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
     !================

     CALL finalize_rewrite_data
   END SUBROUTINE read_linear_transport_data

!====================================================================
!====================================================================
!====== MANDATORY PROCEDURES FOR DEFINING HYPERBOLIC OBJECT =========
!====================================================================
!====================================================================

   FUNCTION flux_linear_transport(this, comp, un) RESULT(vv)  
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      IMPLICIT NONE
      CLASS(linear_transport_type),               INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
      INTEGER,                         INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(SIZE(un, 1), k_dim) :: vv

      INTEGER :: k, n

      DO k=1, k_dim
         DO n=1, SIZE(this%mesh%rr, 2)
            vv(n, k) = this%transport(this%mesh%rr(:,n), k) * un(n, comp)
         END DO
      END DO

   END FUNCTION flux_linear_transport

   SUBROUTINE lambda_linear_transport(this, un, i, j, nij, lambda_max)
      USE space_dim
      IMPLICIT NONE
      CLASS(linear_transport_type),                      INTENT(INOUT) :: this
      REAL(KIND=8), DIMENSION(this%mesh%np, 1), INTENT(IN) :: un
      INTEGER,                                              INTENT(IN) :: i, j
      REAL(KIND=8), DIMENSION(k_dim),                       INTENT(IN) :: nij
      REAL(KIND=8), DIMENSION(2),                          INTENT(OUT) :: lambda_max

      INTEGER, DIMENSION(1) :: i_t, j_t
      INTEGER :: k, ierr
      REAL(KIND = 8), DIMENSION(1, k_dim) :: nij_c
      REAL(KIND = 8), DIMENSION(2)        :: u
      REAL(KIND = 8), DIMENSION(2, k_dim) :: u_transport

      i_t = i
      j_t = j

      DO k = 1, k_dim
         u_transport(1, k) = this%transport(this%mesh%rr(:,i), k)
         u_transport(2, k) = this%transport(this%mesh%rr(:,j), k)
      END DO

      u(1) = SUM(u_transport(1, :) * nij(:))
      u(2) = SUM(u_transport(2, :) * nij(:))

      lambda_max = ABS(u)
      !=== Dummy to avoid warnings
      RETURN
      lambda_max = SUM(un)
      !=== Dummy to avoid warnings
   END SUBROUTINE lambda_linear_transport

   SUBROUTINE construct_linear_transport_bc(this, mesh, LA)
      USE petsc
#include "petsc/finclude/petsc.h"
      USE def_type_mesh,                        ONLY: mesh_type
      USE petsc_csr_LA_module,                  ONLY: petsc_csr_LA
      IMPLICIT NONE
      CLASS(linear_transport_type), INTENT(INOUT) :: this
      TYPE(mesh_type)                            :: mesh
      TYPE(petsc_csr_LA)                         :: LA

      CALL this%bc%rho_bc%set(mesh, "rho "//TRIM(ADJUSTL(this%name)), "DIRICHLET BC PARAMETERS FOR "//TRIM(ADJUSTL(this%name)))
      !===dummy to avoid warnings 
      RETURN
      WRITE(*,*) LA%loc_to_glob
      !===dummy to avoid warnings 
   END SUBROUTINE construct_linear_transport_bc


   SUBROUTINE impose_bc_linear_transport(this, un, mesh, time)
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(linear_transport_type), INTENT(INOUT) :: this
      TYPE(mesh_type) :: mesh
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
      INTEGER :: comp_sys

      !=== Simple Dirichlet boundary conditions
      comp_sys = 1
      un(this%bc%rho_bc%jsd, comp_sys) = this%bc%rho_anal(time, mesh%rr(:, this%bc%rho_bc%jsd))

   END SUBROUTINE impose_bc_linear_transport

 END MODULE linear_transport_type_module
