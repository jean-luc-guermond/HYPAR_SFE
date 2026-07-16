MODULE burgers_type_module
!>> limited global uses to avoid unexpected behaviors

   USE burgers_abstract_hyperbolic_module,   ONLY: hyperbolic_type
   USE burgers_bc_arrays,                    ONLY: burgers_bc_type

   USE read_inputs_module,                   ONLY: rec_length
   USE space_dim,                            ONLY: k_dim
!>> limited global uses to avoid unexpected behaviors

   IMPLICIT NONE

   TYPE argument_burgers_type
   END TYPE argument_burgers_type

   TYPE, EXTENDS(hyperbolic_type)  :: burgers_type
      !===Parameters read from data
      TYPE(burgers_bc_type)          :: bc
      CHARACTER(len=5), DIMENSION(:), ALLOCATABLE :: name_comp
   CONTAINS
      PROCEDURE, PUBLIC  :: init_burgers
      PROCEDURE, PRIVATE :: read_burgers_data
      PROCEDURE :: flux           => flux_burgers
      PROCEDURE :: compute_lambda => lambda_burgers
      PROCEDURE :: construct_bc   => construct_burgers_bc
      PROCEDURE :: impose_bc      => impose_bc_burgers
   END TYPE burgers_type

   ABSTRACT INTERFACE
      FUNCTION function_template_transport(this, rr, comp) RESULT(vv)
         USE space_dim
         IMPORT :: burgers_type
         IMPLICIT NONE
         CLASS(burgers_type)                          :: this
         REAL(KIND = 8), DIMENSION(k_dim), INTENT(IN) :: rr
         INTEGER,                          INTENT(IN) :: comp
         REAL(KIND = 8)                               :: vv
      END FUNCTION function_template_transport
   END INTERFACE

CONTAINS
   SUBROUTINE init_burgers(this, name)
      USE space_dim
      USE my_util, ONLY: error_petsc, to_str
      IMPLICIT NONE
      CLASS(burgers_type), INTENT(INOUT) :: this
      CHARACTER(LEN=*),  INTENT(IN)    :: name

      IF (k_dim == 2) THEN
         CALL error_petsc("BUG in init_burgers => problem should be dimension 1, not "//to_str(k_dim))
      END IF

      this%name = name
      this%syst_dim = 1
      ALLOCATE(this%name_comp(this%syst_dim))

      this%name_comp(1) = 'rho'

      CALL this%read_burgers_data(trim(adjustl(name))//" PARAMETERS")

   END SUBROUTINE init_burgers

   SUBROUTINE read_burgers_data(this, section_name)
     USE read_inputs_module
     USE space_dim
     IMPLICIT NONE
     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name

     CLASS(burgers_type), INTENT(INOUT) :: this
     TYPE(argument_burgers_type)        :: argument_data


     !================
     !=== MANDATORY Reading all data file
     !================
     IF (PRESENT(section_name)) THEN
        CALL read_data_init_list(section_name)
     ELSE
        CALL read_data_init_list()
     END IF

     !================
     !=== We now find the relevant information for this specific burgers data
     !================

     !================
     !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
     !================

     CALL finalize_rewrite_data
   END SUBROUTINE read_burgers_data

!====================================================================
!====================================================================
!====== MANDATORY PROCEDURES FOR DEFINING HYPERBOLIC OBJECT =========
!====================================================================
!====================================================================

   FUNCTION flux_burgers(this, comp, un) RESULT(vv)  
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      IMPLICIT NONE
      CLASS(burgers_type),               INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
      INTEGER,                         INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(SIZE(un, 1), k_dim) :: vv
      INTEGER :: n

      DO n=1, SIZE(this%mesh%rr, 2)
         vv(n, 1) = un(n, comp)**2/2.d0
      END DO

   END FUNCTION flux_burgers

   SUBROUTINE lambda_burgers(this, un, i, j, nij, lambda_max)
      USE space_dim
      IMPLICIT NONE
      CLASS(burgers_type),                      INTENT(INOUT) :: this
      REAL(KIND=8), DIMENSION(this%mesh%np, 1), INTENT(IN) :: un
      INTEGER,                                              INTENT(IN) :: i, j
      REAL(KIND=8), DIMENSION(k_dim),                       INTENT(IN) :: nij
      REAL(KIND=8), DIMENSION(2),                          INTENT(OUT) :: lambda_max

      INTEGER, DIMENSION(1) :: i_t, j_t
      REAL(KIND = 8), DIMENSION(1, k_dim) :: nij_c
      REAL(KIND = 8), DIMENSION(2)        :: u
      REAL(KIND = 8), DIMENSION(2, k_dim) :: u_transport

      i_t = i
      j_t = j

      u(1) = un(i, 1)
      u(2) = un(j, 1)

      !=== works but naive
      ! lambda_max = ABS(u)
      
      !=== More accurate
      IF (u(1) < u(2)) THEN
         lambda_max = ABS(u(1) + u(2))/2
      ELSE
         lambda_max = MAXVAL(ABS(u))
      END IF

      !=== Dummy to avoid warnings
      RETURN
      lambda_max = SUM(nij_c)
      !=== Dummy to avoid warnings
   END SUBROUTINE lambda_burgers

   SUBROUTINE construct_burgers_bc(this, mesh, LA)

      USE def_type_mesh,                        ONLY: mesh_type
      USE petsc_csr_LA_module,                  ONLY: petsc_csr_LA
      IMPLICIT NONE
      CLASS(burgers_type), INTENT(INOUT) :: this
      TYPE(mesh_type)                            :: mesh
      TYPE(petsc_csr_LA)                         :: LA

      CALL this%bc%rho_bc%set(mesh, "rho "//TRIM(ADJUSTL(this%name)), "DIRICHLET BC PARAMETERS FOR "//TRIM(ADJUSTL(this%name)))
      !===dummy to avoid warnings 
      RETURN
      WRITE(*,*) LA%loc_to_glob
      !===dummy to avoid warnings 
   END SUBROUTINE construct_burgers_bc


   SUBROUTINE impose_bc_burgers(this, un, mesh, time)
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(burgers_type), INTENT(INOUT) :: this
      TYPE(mesh_type) :: mesh
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
      INTEGER :: comp_sys

      !=== Simple Dirichlet boundary conditions
      comp_sys = 1
      un(this%bc%rho_bc%jsd, comp_sys) = this%bc%rho_anal(time, mesh%rr(:, this%bc%rho_bc%jsd))

   END SUBROUTINE impose_bc_burgers

 END MODULE burgers_type_module
