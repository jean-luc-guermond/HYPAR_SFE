MODULE abstract_hyperbolic_module

!>> limited global uses to avoid unexpected behaviors
#include "petsc/finclude/petsc.h"
   USE petsc
   USE petsc_tools,                          ONLY: array_to_petsc_vec
   USE Butcher_tableau
   USE hyperbolic_matrices_module,           ONLY: hyperbolic_matrices_type
   USE hyperbolic_bc_tools,                  ONLY: construct_udotn
   USE cell_limiting_engine_parallel_module, ONLY: limiting_type, limiting_functional_type, limiting_all_functional_type
   USE def_type_mesh,                        ONLY: mesh_type, petsc_csr_LA
   USE read_inputs_module,                   ONLY: rec_length
   USE space_dim,                            ONLY: k_dim
!>> limited global uses to avoid unexpected behaviors

   IMPLICIT NONE
   INTEGER, PRIVATE, PARAMETER :: LUMPED_MASS=1, QUASI_CONSISTENT_MASS=2, CONSISTENT_MASS=3
   CHARACTER(LEN=20), DIMENSION(3), PRIVATE, PARAMETER :: list_mass = &
               [CHARACTER(LEN=20) :: 'lumped', 'quasi_consistent', 'consistent']
   INTEGER, PRIVATE, PARAMETER :: METHOD_VISCOUS=1, METHOD_HIGH=2, METHOD_GALERKIN=3
   CHARACTER(LEN=20), DIMENSION(3), PRIVATE, PARAMETER  :: list_method = &
               [CHARACTER(LEN=20) :: 'viscous', 'high', 'galerkin']

   TYPE argument_hyperbolic_type
      CHARACTER(LEN=rec_length) :: CFL                     = '=== CFL ? ==='
      CHARACTER(LEN=rec_length) :: char_method             = '=== Which method to solve conservation equation (viscous,high,galerkin)? ==='
      CHARACTER(LEN=rec_length) :: which_high_method       = '=== Which high order method (1/2)? ==='
      CHARACTER(LEN=rec_length) :: if_hybrid_mesh_limiting = '=== Do we use hybrid Pk/P1 meshes for limiting? ==='
      CHARACTER(LEN=rec_length) :: erk_sv                  = '=== ERK ? ==='
      CHARACTER(LEN=rec_length) :: char_which_mass         = '=== Which mass matrix: lumped, quasi_consistent, consistent? ==='
      CHARACTER(LEN=rec_length) :: nb_correction_mass      = '=== For quasi_consistent mass, how many corrections? (0=lumped_mass) ==='
   END TYPE argument_hyperbolic_type

   TYPE, ABSTRACT :: hyperbolic_type
      !===Parameters read from data
      REAL(KIND=8)                 :: CFL                     = 0.5d0
      CHARACTER(LEN=rec_length)    :: char_method             = 'viscous'
      INTEGER                      :: method                  = METHOD_VISCOUS
      INTEGER                      :: which_high_method       = 1
      INTEGER                      :: erk_sv                  = -21
      LOGICAL                      :: if_hybrid_mesh_limiting = .TRUE.
      CHARACTER(LEN=rec_length)    :: char_which_mass         = 'lumped'
      INTEGER                      :: which_mass              = LUMPED_MASS
      INTEGER                      :: nb_correction_mass      = 1
      !===Parameters built along way
      MPI_Comm :: communicator
      Vec          :: x1vec, x2vec, x2_ghost, vec_loc, x3vec
      Vec          :: x4vec, x5vec !!!!!Conveniance vectors to be used only inside procedures!!!!
      CHARACTER(LEN=:), ALLOCATABLE :: name
      INTEGER                       :: syst_dim
      REAL(KIND = 8) :: dt, time, final_time
      INTEGER, DIMENSION(:), ALLOCATABLE :: tab
      TYPE(mesh_type),     POINTER :: mesh, mesh_L
      TYPE(petsc_csr_LA),  POINTER :: LA,   LA_L
      TYPE(BT),             PUBLIC :: ERK
      TYPE(hyperbolic_matrices_type), POINTER :: matrices, matrices_L
      TYPE(limiting_type)            :: limiting
      ! CLASS(limiting_functional_type), DIMENSION(:), POINTER :: limiting_functionals
      TYPE(limiting_all_functional_type) :: limiting_all_functionals
      PROCEDURE(template_eta_commute), NOPASS,        POINTER :: eta_commute
   CONTAINS
      PROCEDURE, PUBLIC   :: init_hyperbolic
      PROCEDURE, PRIVATE  :: read_hyperbolic_data, init_vectors
      PROCEDURE, PUBLIC   :: update
      PROCEDURE, PRIVATE  :: compute_dij, compute_dt, invert_mass, commutator, apply_limiting
      PROCEDURE(template_flux),         DEFERRED :: flux
      PROCEDURE(template_lambda),       DEFERRED :: compute_lambda
      PROCEDURE, NOPASS                          :: construct_udotn => construct_udotn
      PROCEDURE(template_construct_bc), DEFERRED :: construct_bc
      PROCEDURE(template_impose_bc),    DEFERRED :: impose_bc
   END TYPE hyperbolic_type


   ABSTRACT INTERFACE

      FUNCTION template_eta_commute(un) RESULT(eta)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(:, :), INTENT(IN)    :: un
         REAL(KIND=8), DIMENSION(SIZE(un,1))         :: eta
      END FUNCTION template_eta_commute

      SUBROUTINE template_construct_bc(this, mesh, LA)
         USE def_type_mesh
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type), INTENT(INOUT) :: this
         TYPE(mesh_type)           :: mesh
         TYPE(petsc_csr_LA)        :: LA
      END SUBROUTINE template_construct_bc

      SUBROUTINE template_impose_bc(this, un, mesh, time)
         USE def_type_mesh
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type),       INTENT(INOUT) :: this
         TYPE(mesh_type)                                :: mesh
         REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
         REAL(KIND = 8), INTENT(IN)                     :: time
      END SUBROUTINE template_impose_bc

      SUBROUTINE template_lambda(this, un, i, j, nij, lambda_max)
         USE space_dim
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type),                               INTENT(INOUT) :: this
         REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(IN) :: un
         INTEGER,                                              INTENT(IN) :: i, j
         REAL(KIND=8), DIMENSION(k_dim),                       INTENT(IN) :: nij
         REAL(KIND=8), DIMENSION(2),                           INTENT(OUT) :: lambda_max
      END SUBROUTINE template_lambda

      FUNCTION template_flux(this, comp, un) RESULT(vv)
         USE space_dim
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type),          INTENT(INOUT) :: this
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
         INTEGER,                         INTENT(IN) :: comp
         REAL(KIND = 8), DIMENSION(SIZE(un, 1), k_dim) :: vv
      END FUNCTION template_flux
   END INTERFACE

CONTAINS

   SUBROUTINE init_hyperbolic(this, communicator, name, mesh, LA, times, limiting_functionals)
      USE my_util,            ONLY: error_petsc, to_str
      USE space_dim
      USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
      USE construct_mesh, ONLY: generate_boundary_structure, refine_mesh_1D2D_Pk2P1, create_gauss_points_1D_2D


      IMPLICIT NONE
      CLASS(hyperbolic_type), INTENT(INOUT) :: this
      MPI_Comm,                   INTENT(IN) :: communicator
      CHARACTER(100),             INTENT(IN) :: name
      TYPE(mesh_type), TARGET,    INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
      REAL(KIND = 8), DIMENSION(2) :: times
      TYPE(limiting_functional_type), DIMENSION(:), TARGET :: limiting_functionals

      this%name = name
      this%mesh => mesh
      this%communicator = communicator
      this%LA   => LA

      this%time = times(1) !<==initial_time
      this%final_time = times(2) !<==final_time

      CALL this%read_hyperbolic_data("HYPERBOLIC PARAMETERS FOR "//trim(adjustl(this%name)))

      !=== Build ERK structure
      this%ERK%sv = this%erk_sv
      CALL this%ERK%init()

      !===Matrices
      ALLOCATE(this%matrices)
      this%matrices%method     = this%method
      this%matrices%which_mass = this%which_mass
      CALL this%matrices%construct(this%communicator, this%mesh, this%LA)

      !===Hybrid meshes
      IF (mesh%info%type_fe>1) THEN
         IF (mesh%rank==0) WRITE(*,*) "Building hybrid mesh for ", this%name
         !=== creating low order stencil
         ALLOCATE(this%mesh_L)
         ALLOCATE(this%LA_L)
         ALLOCATE(this%matrices_L)
         CALL refine_mesh_1D2D_Pk2P1(this%mesh, this%mesh_L)
         !=== generating associated P1 gauss points
         CALL create_gauss_points_1D_2D(this%mesh_L, 1)
         !=== adding surface elements (for Dirichlet and periodicity)
         CALL generate_boundary_structure(this%mesh_L)
         !=== Generating sparse structures
         CALL st_aij_csr_glob_block_with_extra_layer(this%communicator, 1, this%mesh_L, this%LA_L)
         this%matrices_L%method = METHOD_VISCOUS
         CALL this%matrices_L%construct(this%communicator, this%mesh_L, this%LA_L)
      ELSE
         this%mesh_L     => this%mesh
         this%LA_L       => this%LA
         this%matrices_L => this%matrices
      END IF

      !===Goshting structures
      CALL this%init_vectors

      !===Build boundary conditions
      CALL this%construct_bc(this%mesh, this%LA)

      !===Limiting
      IF (this%if_hybrid_mesh_limiting) THEN
         CALL this%limiting%init(this%name, this%mesh_L, this%LA_L)
         CALL this%limiting_all_functionals%init(limiting_functionals, name, this%LA_L, this%mesh_L)
      ELSE
         CALL this%limiting%init(this%name, this%mesh, this%LA)
         CALL this%limiting_all_functionals%init(limiting_functionals, name, this%LA, this%mesh)
      END IF
   END SUBROUTINE init_hyperbolic

   SUBROUTINE read_hyperbolic_data(this, section_name)
      USE read_inputs_module
      USE my_util, ONLY: get_tab_idx_char
      IMPLICIT NONE
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name

      CLASS(hyperbolic_type), INTENT(INOUT) :: this
      TYPE(argument_hyperbolic_type)        :: argument_data


      !================
      !=== MANDATORY Reading all data file
      !================
      IF (PRESENT(section_name)) THEN
         CALL read_data_init_list(section_name)
      ELSE
         CALL read_data_init_list()
      END IF

      !================
      !=== We now find the relevant information for this specific Euler data
      !================
      !===CFL
      CALL read_data(argument_data%CFL, this%CFL, opt_name=this%name)

      !===ERK
      CALL read_data(argument_data%erk_sv, this%erk_sv, opt_name=this%name)

      !===Method order
      CALL read_data(argument_data%char_method, this%char_method, opt_name=this%name)
      CALL get_tab_idx_char(this%char_method, list_method, this%method)

      !===which high order method
      CALL read_data(argument_data%which_high_method, this%which_high_method, opt_name=this%name, &
                     opt_add=this%method==METHOD_HIGH)

      !===mass matrix
      CALL read_data(argument_data%char_which_mass, this%char_which_mass, opt_name=this%name)
      CALL get_tab_idx_char(this%char_which_mass, list_mass, this%which_mass)

      !===nb_correction_mass
      CALL read_data(argument_data%nb_correction_mass, this%nb_correction_mass, opt_name=this%name)

      !===if_hybrid_mesh_limiting
      CALL read_data(argument_data%if_hybrid_mesh_limiting, this%if_hybrid_mesh_limiting, &
                     opt_name=this%name, opt_add=this%method==METHOD_HIGH)

      !================
      !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
      !================

      CALL finalize_rewrite_data
   END SUBROUTINE read_hyperbolic_data

   SUBROUTINE update(this,un_in)
     IMPLICIT NONE
     CLASS(hyperbolic_type)                                             :: this
     REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim)               :: un_in
     REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s+1) :: urk
     REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s)   :: flux_rk_at_dof
     INTEGER  :: stage
     urk(:,:,1) = un_in
     DO stage = 2, this%ERK%s+1
        CALL one_step_ERK(this,stage,urk,flux_rk_at_dof)
     END DO
     un_in = urk(:,:,this%ERK%s+1)
     this%time = this%time + this%dt
   END SUBROUTINE update

   SUBROUTINE one_step_ERK(this,stage,urk,flux_rk_at_dof)
      USE my_util, ONLY : error_petsc, to_str
      IMPLICIT NONE
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s+1)  :: urk
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s)    :: flux_rk_at_dof
      INTEGER, INTENT(IN) :: stage

      SELECT CASE(this%method)
      CASE(METHOD_GALERKIN)
         CALL one_step_ERK_galerkin(this,stage,urk,flux_rk_at_dof)
      CASE(METHOD_VISCOUS)
         CALL one_step_ERK_viscous(this,stage,urk)
      CASE(METHOD_HIGH)
         SELECT CASE(this%which_high_method)
         CASE(1)
            CALL one_step_ERK_high(this,stage,urk,flux_rk_at_dof)
         CASE(2)
            CALL one_step_ERK_high_2(this,stage,urk,flux_rk_at_dof)
         CASE DEFAULT
            CALL error_petsc("Bug in one_step_ERK: wrong high order method "//to_str(this%which_high_method))
         END SELECT
      CASE DEFAULT
         CALL error_petsc("wrong method in one_step_ERK")
      END SELECT
   END SUBROUTINE one_step_ERK

   SUBROUTINE one_step_ERK_galerkin(this,stage,urk,flux_rk_at_dof)
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      USE cell_limiting_engine_module
      USE sub_plot
      USE compute_periodic, ONLY : periodic_rhs_petsc, periodic_vector_petsc
      USE st_matrix, ONLY: extract_through_ghost
      USE fem_tn
      IMPLICIT NONE
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s+1)  :: urk
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s)    :: flux_rk_at_dof
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim)                :: un_temp
      REAL(KIND = 8), DIMENSION(this%mesh%np)                               :: rk
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim, this%syst_dim)         :: flux_array
      REAL(KIND = 8), DIMENSION(this%mesh%np,this%limiting_all_functionals%nl) :: bounds
      INTEGER, INTENT(IN) :: stage
      INTEGER :: comp, k, l, ierr, stage_prime, stage0
      REAL(KIND = 8) :: time_stage

      !===flux_array: flux at l=stage
      DO comp=1, this%syst_dim
         flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage-1)) !<== Notice: Flux at l=stage
      END DO

      DO comp = 1, this%syst_dim
         !=== flux_rk_at_dof(:,comp,stage-1): -sum_j(f(uj(stage)) cj) and store in flux_rk_at_dof(stage)
         CALL VecSet(this%x3vec, 0.d0, ierr)
         DO k = 1, k_dim
            !=== set flux_k in x1vec
            CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%LA, 'insert')
            !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
            CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
            !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
            CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
         END DO
         CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, flux_rk_at_dof(:,comp,stage-1), opt_assemble=.FALSE.) !<== Store sum_j(f(uj)cij) at l=stage
      END DO

      !=== dt
      IF (stage==2) THEN !< ==Compute time step only once per ERK step
         CALL this%compute_dij(flux_array,urk(:,:,stage_prime), bounds) !<== Notice: State at l=stage'
         CALL this%compute_dt
         ! IF (this%time+this%dt.GE.this%final_time) THEN
         !   this%dt = this%final_time-this%time
         ! END IF
      END IF
      time_stage = this%time+this%ERK%C(stage)*this%dt !<== Wait for dt to be computed

      !=== stage0
      IF (this%ERK%sv<0) THEN
         stage0 = stage_prime
      ELSE
         stage0 = 1
      END IF

      !=== ERK update for each component
      DO comp = 1, this%syst_dim

         !=== rk: sum_l MatRK_(s,l) f_l
         rk =0.d0
         DO l = 1, stage-1
            rk =  rk + this%ERK%MatRK(stage,l)*flux_rk_at_dof(:,comp,l)
         END DO

         !=== rk in x3vec
         CALL array_to_petsc_vec(rk, this%x2vec, this%LA, 'insert')
         CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x2vec, this%LA)

         !=== Inverting mass matrix and updating un with dt
         CALL this%invert_mass(this%x2vec, this%x3vec, this%which_mass)

         !=== set un(comp) at l=stage0 in x1vec
         CALL array_to_petsc_vec(urk(:,comp,stage0), this%x1vec, this%LA, 'insert') !<== Notice un at l=stage0
         !=== x3 <-- dt*x3 + un (x1 <-- un a few lines above)
         CALL VecAYPX(this%x3vec, this%dt, this%x1vec, ierr)
         !=== Manually make un periodic and extract the result
         CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x3vec, this%LA)
         CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
      END DO

      !===Boundary conditions
      CALL this%impose_bc(urk(:,:,stage), this%mesh, this%time)
   END SUBROUTINE one_step_ERK_galerkin

   SUBROUTINE one_step_ERK_viscous(this,stage,urk)
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      USE cell_limiting_engine_module
      USE sub_plot
      USE compute_periodic, ONLY : periodic_rhs_petsc, periodic_vector_petsc
      USE st_matrix, ONLY: extract_through_ghost
      USE fem_tn
      IMPLICIT NONE
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, this%syst_dim, this%ERK%s+1)  :: urk
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, this%syst_dim)                :: un_temp
      REAL(KIND = 8), DIMENSION(this%mesh_L%np)                               :: rk
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, k_dim, this%syst_dim)         :: flux_array
      REAL(KIND = 8), DIMENSION(this%mesh_L%np,this%limiting_all_functionals%nl) :: bounds
      INTEGER, INTENT(IN) :: stage
      INTEGER :: comp, k, l, ierr, stage_prime, stage0
      REAL(KIND = 8) :: time_stage

      stage_prime = this%ERK%lp_of_l(stage) !<== Method should work with incremental ERK method
      un_temp = urk(:,:,stage_prime)      !<== Method should work with incremental ERK method
      DO comp=1, this%syst_dim
         flux_array(:, :, comp) = this%flux(comp, un_temp)
      END DO

      !===compute dijL and dt
      CALL this%compute_dij(flux_array,un_temp, bounds)
      IF (stage==2) THEN !< ==Compute time step
         CALL this%compute_dt
         ! IF (this%time+this%dt.GE.this%final_time) THEN
         !   this%dt = this%final_time-this%time
         ! END IF
      END IF
      time_stage = this%time+this%ERK%C(stage)*this%dt !<== Wait for dt to be computed

      DO comp = 1, this%syst_dim
         CALL VecZeroEntries(this%x3vec, ierr)
         DO k = 1, k_dim
            !=== set flux_k in x1vec
            CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%LA_L, 'insert')
            !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
            CALL MatMult(this%matrices_L%cij(k), this%x1vec, this%x2vec, ierr)
            !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
            CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
         END DO

         !=== set un(comp) in x1vec
         CALL array_to_petsc_vec(un_temp(:, comp), this%x1vec, this%LA_L, 'insert')
         !=== add dij un(comp)to x3vec in x2vec
         CALL MatMultAdd(this%matrices_L%dijL, this%x1vec, this%x3vec, this%x2vec, ierr)
         CALL periodic_rhs_petsc(this%mesh_L%per%nb_bords, this%mesh_L%per%list, this%mesh_L%per%perlist, this%x2vec, this%LA_L)

         ! !=== x3 <-- x2 / lumped_mass by default in viscous method
         CALL this%invert_mass(this%x2vec, this%x3vec, LUMPED_MASS)

         !=== x3 <-- un + x3*dt   (x1 <-- un few lines above)
         CALL VecAYPX(this%x3vec, this%ERK%inc_C(stage)*this%dt, this%x1vec, ierr) !<==time step is ERK%inc_C(stage)*this%dt
         CALL periodic_vector_petsc(this%mesh_L%per%nb_bords, this%mesh_L%per%list, this%mesh_L%per%perlist, this%x3vec, this%LA_L)
         !=== un+1 <-- x3
         CALL extract_through_ghost(this%x3vec, 1, 1, this%LA_L, urk(:, comp, stage), opt_assemble=.FALSE.)
      END DO

      CALL this%impose_bc(urk(:,:,stage), this%mesh_L, time_stage)

   END SUBROUTINE one_step_ERK_viscous

   SUBROUTINE one_step_ERK_high(this,stage,urk,flux_rk_at_dof)
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      USE cell_limiting_engine_module
      USE sub_plot
      USE compute_periodic, ONLY : periodic_rhs_petsc, periodic_vector_petsc
      USE st_matrix, ONLY: extract_through_ghost
      USE fem_tn
      IMPLICIT NONE
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s+1)  :: urk
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s)    :: flux_rk_at_dof
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim)                :: un_temp
      REAL(KIND = 8), DIMENSION(this%mesh%np)                               :: rk
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim, this%syst_dim)         :: flux_array
      REAL(KIND = 8), DIMENSION(this%mesh%np,this%limiting_all_functionals%nl) :: bounds
      INTEGER, INTENT(IN) :: stage
      INTEGER :: comp, k, l, ierr, stage_prime, stage0
      REAL(KIND = 8) :: time_stage

      !===flux_array: flux at l=stage
      DO comp=1, this%syst_dim
         flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage-1)) !<== Notice: Flux at l=stage
      END DO

      DO comp = 1, this%syst_dim
         !=== flux_rk_at_dof(:,comp,stage-1): -sum_j(f(uj(stage)) cj) and store in flux_rk_at_dof(stage)
         CALL VecSet(this%x3vec, 0.d0, ierr)
         DO k = 1, k_dim
            !=== set flux_k in x1vec
            CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%LA, 'insert')
            !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
            CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
            !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
            CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
         END DO
         CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, flux_rk_at_dof(:,comp,stage-1), opt_assemble=.FALSE.) !<== Store sum_j(f(uj)cij) at l=stage
      END DO
      !=== dij
      stage_prime = this%ERK%lp_of_l(stage) !<== Method should work with incremental ERK method
      IF (stage-1 .NE. stage_prime) THEN
         DO comp=1, this%syst_dim
            flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage_prime)) !<== Notice: Flux at l=stage'
         END DO
      END IF
      CALL this%compute_dij(flux_array,urk(:,:,stage_prime), bounds) !<== Notice: State at l=stage'

      !=== dt
      IF (stage==2) THEN !< ==Compute time step only once per ERK step
         CALL this%compute_dt
         ! IF (this%time+this%dt.GE.this%final_time) THEN
         !   this%dt = this%final_time-this%time
         ! END IF
      END IF
      time_stage = this%time+this%ERK%C(stage)*this%dt !<== Wait for dt to be computed

      !=== stage0
      IF (this%ERK%sv<0) THEN
         stage0 = stage_prime
      ELSE
         stage0 = 1
      END IF

      !=== ERK update for each component
      DO comp = 1, this%syst_dim
         !=== rk: sum_l MatRK_(s,l) f_l
         rk =0.d0
         DO l = 1, stage-1
            rk =  rk + this%ERK%MatRK(stage,l)*flux_rk_at_dof(:,comp,l)
         END DO

         !=== rk in x3vec
         CALL array_to_petsc_vec(rk, this%x3vec, this%LA, 'insert')
         !=== set un(comp) at l=stage' in x1vec to compute viscous contribution
         CALL array_to_petsc_vec(urk(:,comp,stage_prime), this%x1vec, this%LA, 'insert') !<== l=stage'

         !=== add dij un(comp) to x3vec in x2vec
         CALL MatMultAdd(this%matrices%dijH, this%x1vec, this%x3vec, this%x2vec, ierr)
         CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x2vec, this%LA)
         !=== Inverting mass matrix and updating un with dt
         CALL this%invert_mass(this%x2vec, this%x3vec, this%which_mass)

         !=== set un(comp) at l=stage0 in x1vec
         CALL array_to_petsc_vec(urk(:,comp,stage0), this%x1vec, this%LA, 'insert') !<== Notice un at l=stage0
         !=== x3 <-- dt*x3 + un (x1 <-- un a few lines above)
         CALL VecAYPX(this%x3vec, this%dt, this%x1vec, ierr)
         !=== Manually make un periodic and extract the result
         CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x3vec, this%LA)
         CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
      END DO
      !===Limiting
      !CALL ns_mass_PAR(this%mesh, urk(:,1,stage), mass_before, this%communicator)
      IF (this%limiting%if_limiting) THEN
         CALL this%apply_limiting(urk(:,:,stage), bounds)
      END IF
      !CALL ns_mass_PAR(this%mesh, urk(:,1,stage), mass_after, this%communicator)
      !IF (this%mesh%rank==0) THEN
      !   WRITE(*,*) 'mass before limiting = ', mass_before
      !   WRITE(*,*) 'mass after limiting = ', mass_after
      !END IF

      !===Periodicity
      DO comp = 1, this%syst_dim
         CALL array_to_petsc_vec(urk(:,comp,stage), this%x1vec, this%LA, 'insert')
         CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x1vec, this%LA)
         CALL extract_through_ghost(this%x1vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
      END DO

      !===Boundary conditions
      CALL this%impose_bc(urk(:,:,stage), this%mesh, this%time)
   END SUBROUTINE one_step_ERK_high


   SUBROUTINE one_step_ERK_high_2(this,stage,urk,flux_rk_at_dof)
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      USE cell_limiting_engine_module
      USE sub_plot
      USE compute_periodic, ONLY : periodic_rhs_petsc, periodic_vector_petsc
      USE st_matrix, ONLY: extract_through_ghost
      USE fem_tn
      IMPLICIT NONE
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s+1)  :: urk
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s)    :: flux_rk_at_dof
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim)                :: un_temp, u_L
      REAL(KIND = 8), DIMENSION(this%mesh%np)                               :: rk, zero_vec
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim, this%syst_dim)         :: flux_array
      REAL(KIND = 8), DIMENSION(this%mesh%np,this%limiting_all_functionals%nl) :: bounds
      REAL(KIND = 8), DIMENSION(this%mesh%np) :: temp_bounds, min_bounds
      INTEGER, INTENT(IN) :: stage
      INTEGER :: comp, k, l, ierr, stage_prime, stage0, nl, ni, nj, i, m, me, nw
      REAL(KIND = 8) :: time_stage

      !=== Bounds computations using u_L
      CALL one_step_ERK_viscous(this,stage,urk)
      u_L = urk(:, :, stage)
      IF (this%limiting%if_limiting) THEN
         zero_vec = 0.d0
         DO nl=1, this%limiting_all_functionals%nl
            temp_bounds = this%limiting_all_functionals%limiting_functionals(nl)%psi(u_L(:, :), zero_vec)
            
            min_bounds = temp_bounds
            me = this%mesh_L%me
            nw = this%mesh_L%gauss%n_w
            DO m = 1, me
               DO ni = 1, nw
                  i = this%mesh_L%jj(ni,m)
                  min_bounds(i) = MIN(min_bounds(i),MINVAL(temp_bounds(this%mesh_L%jj(:,m))))
               END DO
            END DO
            CALL array_to_petsc_vec(min_bounds, this%x3vec, this%LA, "min")
            CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, bounds(:, nl), opt_assemble=.FALSE.)
         END DO
      END IF
      !=== Bounds computations using u_L
         
      !===flux_array: flux at l=stage
      DO comp=1, this%syst_dim
         flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage-1)) !<== Notice: Flux at l=stage
      END DO

      DO comp = 1, this%syst_dim
         !=== flux_rk_at_dof(:,comp,stage-1): -sum_j(f(uj(stage)) cj) and store in flux_rk_at_dof(stage)
         CALL VecSet(this%x3vec, 0.d0, ierr)
         DO k = 1, k_dim
            !=== set flux_k in x1vec
            CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%LA, 'insert')
            !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
            CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
            !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
            CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
         END DO
         CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, flux_rk_at_dof(:,comp,stage-1), opt_assemble=.FALSE.) !<== Store sum_j(f(uj)cij) at l=stage
      END DO
      !===flux_array: flux at l=stage

      ! !=== dij
      ! IF (stage-1 .NE. stage_prime) THEN
      !    DO comp=1, this%syst_dim
      !       flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage_prime)) !<== Notice: Flux at l=stage'
      !    END DO
      ! END IF
      ! CALL this%compute_dij(flux_array,urk(:,:,stage_prime), bounds) !<== Notice: State at l=stage'
      ! !=== dij
      
      !=== dt
      stage_prime = this%ERK%lp_of_l(stage) !<== Method should work with incremental ERK method
      IF (stage==2) THEN !< ==Compute time step only once per ERK step
         CALL this%compute_dt
         ! IF (this%time+this%dt.GE.this%final_time) THEN
         !   this%dt = this%final_time-this%time
         ! END IF
      END IF
      time_stage = this%time+this%ERK%C(stage)*this%dt !<== Wait for dt to be computed

      !=== stage0
      IF (this%ERK%sv<0) THEN
         stage0 = stage_prime
      ELSE
         stage0 = 1
      END IF

      !=== ERK update for each component
      DO comp = 1, this%syst_dim
         !=== rk: sum_l MatRK_(s,l) f_l
         rk =0.d0
         DO l = 1, stage-1
            rk =  rk + this%ERK%MatRK(stage,l)*flux_rk_at_dof(:,comp,l)
         END DO

         !=== rk in x3vec
         CALL array_to_petsc_vec(rk, this%x3vec, this%LA, 'insert')
         !=== set un(comp) at l=stage' in x1vec to compute viscous contribution
         CALL array_to_petsc_vec(urk(:,comp,stage_prime), this%x1vec, this%LA, 'insert') !<== l=stage'
         !=== add dij un(comp) to x3vec in x2vec
         CALL MatMultAdd(this%matrices%dijH, this%x1vec, this%x3vec, this%x2vec, ierr)
         CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x2vec, this%LA)
         !=== Inverting mass matrix and updating un with dt
         CALL this%invert_mass(this%x2vec, this%x3vec, this%which_mass)

         !=== set un(comp) at l=stage0 in x1vec
         CALL array_to_petsc_vec(urk(:,comp,stage0), this%x1vec, this%LA, 'insert') !<== Notice un at l=stage0
         !=== x3 <-- dt*x3 + un (x1 <-- un a few lines above)
         CALL VecAYPX(this%x3vec, this%dt, this%x1vec, ierr)
         !=== Manually make un periodic and extract the result
         CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x3vec, this%LA)
         CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
      END DO
      !===Limiting
      !CALL ns_mass_PAR(this%mesh, urk(:,1,stage), mass_before, this%communicator)
      IF (this%limiting%if_limiting) THEN
         CALL this%apply_limiting(urk(:,:,stage), bounds)
      END IF
      !CALL ns_mass_PAR(this%mesh, urk(:,1,stage), mass_after, this%communicator)
      !IF (this%mesh%rank==0) THEN
      !   WRITE(*,*) 'mass before limiting = ', mass_before
      !   WRITE(*,*) 'mass after limiting = ', mass_after
      !END IF

      !===Periodicity
      DO comp = 1, this%syst_dim
         CALL array_to_petsc_vec(urk(:,comp,stage), this%x1vec, this%LA, 'insert')
         CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x1vec, this%LA)
         CALL extract_through_ghost(this%x1vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
      END DO

      !===Boundary conditions
      CALL this%impose_bc(urk(:,:,stage), this%mesh, this%time)
   END SUBROUTINE one_step_ERK_high_2

!========================================================
!========== PRIVATE PROCEDURES ==========================
!========================================================


   SUBROUTINE compute_dt(this)
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8)         :: dt_min_glob
      INTEGER                :: ierr

      CALL MatGetDiagonal(this%matrices_L%dijL, this%x1vec, ierr)
      CALL VecAbs(this%x1vec, ierr)
      CALL VecPointWiseDivide(this%x2vec, this%matrices_L%lump_mass_vec, this%x1vec, ierr)
      CALL VecMin(this%x2vec, PETSC_NULL_INTEGER, dt_min_glob, ierr)

      this%dt = this%ERK%s * this%CFL * dt_min_glob / 2.d0
      !===Notice rescale of time step with this%ERK%s

   END SUBROUTINE compute_dt

   SUBROUTINE compute_dij(this, flux_array, un, bounds)
      USE space_dim
      USE petsc
      USE def_type_mesh
      USE arbitrary_eos_lambda_module
      USE compute_periodic
      USE st_matrix, ONLY: extract_through_ghost
      USE my_util
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      TYPE(mesh_type), POINTER :: mesh
      TYPE(petsc_csr_LA), POINTER :: LA
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, k_dim, this%syst_dim) :: flux_array
      REAL(KIND = 8), DIMENSION(:, :) :: un, bounds
      INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge, nl, comp
      INTEGER, DIMENSION(1) :: i_t, j_t, idx, jdx
      REAL(KIND = 8), DIMENSION(1, k_dim) :: nij_c
      REAL(KIND = 8), DIMENSION(1) :: norm_c, dijL_c, dijH_c
      REAL(KIND = 8), DIMENSION(2) :: lambda_max
      REAL(KIND = 8), DIMENSION(this%mesh_L%gauss%n_w,this%mesh_L%gauss%n_w) :: mat_loc, mat_h_loc, norm_c_loc
      REAL(KIND = 8), DIMENSION(this%mesh_L%gauss%n_w,this%mesh_L%gauss%n_w, k_dim) :: nij_c_loc
      REAL(KIND = 8), DIMENSION(this%mesh_L%gauss%n_w*this%mesh_L%gauss%n_w, k_dim) :: nij_c_loc_vec
      INTEGER,        DIMENSION(this%mesh_L%gauss%n_w) :: idx_loc

      REAL(KIND = 8) :: max_lambda
      REAL(KIND = 8), DIMENSION(1, this%syst_dim) :: uijbar  !<=== FIXME
      REAL(KIND = 8), DIMENSION(1), PARAMETER   :: zero=0.d0
      REAL(KIND = 8), DIMENSION(this%mesh_L%np) :: zero_vec
      LOGICAL, DIMENSION(this%mesh_L%medge) :: virgin_edge
      REAL(KIND = 8), DIMENSION(this%mesh_L%np)  :: alpha !<==commutator in (0,1)

      mesh => this%mesh_L
      LA   => this%LA_L

      !===Compute dijL
      CALL MatZeroEntries(this%matrices_L%dijL, ierr)
      IF (this%method==METHOD_HIGH) THEN
         CALL this%commutator(un, alpha)
         CALL MatZeroEntries(this%matrices%dijH, ierr)
         IF (this%limiting%if_limiting .AND. this%which_high_method==1) THEN
            zero_vec = 0.d0
            DO nl=1, this%limiting_all_functionals%nl
               bounds(:, nl) = this%limiting_all_functionals%limiting_functionals(nl)%psi(un(:, :), zero_vec)
            END DO
         END IF
      END IF


      virgin_edge = .TRUE.
      nw = mesh%gauss%n_w

      DO m = 1, mesh%me
         mat_loc = 0.d0
         mat_h_loc = 0.d0
         DO ni = 1, nw
            i = mesh%jj(ni, m)
            idx_loc(ni) = LA%loc_to_glob(1, i) - 1
         END DO

         norm_c_loc(:,:) = this%matrices_L%cij_norm_loc_array(:,:,m)
         nij_c_loc(:, :, :) = this%matrices_L%nij_loc_array(:,:,:,m)

         DO n = 1, mesh%gauss%n_e
            IF (mesh%attr_e(mesh%jce(n, m))) THEN
               edge = mesh%jce_loc(n, m)
               IF (.NOT. virgin_edge(edge)) CYCLE
               virgin_edge(edge) = .FALSE.

               ni = MOD(n, nw) + 1
               nj = MOD(n + 1, nw) + 1
               i = mesh%jj(ni, m)
               j = mesh%jj(nj, m)
               i_t = i
               j_t = j

               CALL this%compute_lambda(un, i, j, nij_c_loc(ni, nj, :), lambda_max)

               max_lambda = MAXVAL(lambda_max)
               dijL_c = max_lambda * norm_c_loc(ni, nj) 

               IF (mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j
                     CALL this%compute_lambda(un, j, i, nij_c_loc(nj, ni, :), lambda_max)
                     norm_c(1) = norm_c_loc(nj, ni)
                     dijL_c = MAX(dijL_c, MAXVAL(lambda_max) * norm_c)
                     max_lambda = MAX(max_lambda,MAXVAL(lambda_max))
               END IF

               mat_loc(ni, nj) = dijL_c(1)
               mat_loc(nj, ni) = dijL_c(1)

               IF (this%method==METHOD_HIGH) THEN
                  dijH_c = dijL_c*(alpha(i)+alpha(j))/2
                  mat_h_loc(ni, nj) = dijH_c(1)
                  mat_h_loc(nj, ni) = dijH_c(1)

                  !===Compute low-order update to estimate bounds
                  IF (this%limiting%if_limiting .AND. this%which_high_method==1) THEN

                     nij_c(1, :) = nij_c_loc(ni, nj, :)

                     DO comp=1, this%syst_dim
                        uijbar(1,comp) =  (un(i, comp)+un(j, comp))/2 - &
                        SUM((flux_array(j, :, comp)-flux_array(i, :, comp))*nij_c(1, :))/(2*max_lambda)
                     END DO

                     DO nl=1, this%limiting_all_functionals%nl
                        bounds(i:i, nl) = MIN(bounds(i, nl), this%limiting_all_functionals%limiting_functionals(nl)%psi(uijbar, zero))
                        bounds(j:j, nl) = MIN(bounds(j, nl), this%limiting_all_functionals%limiting_functionals(nl)%psi(uijbar, zero))
                     END DO

                     !===End compute low-order update to estimate bounds
                  END IF
               END IF ! END METHOD_HIGH
            END IF ! END mesh_attr
         END DO ! END mesh%gauss%n_e
         ! <=== WARNING filling by blocks here works because those blocks are symmetric
         !=== Matrices in petsc are row-oriented
         CALL MatSetValues(this%matrices_L%dijL, nw, idx_loc, nw, idx_loc, mat_loc, ADD_VALUES, ierr) 
         IF (this%method==METHOD_HIGH) THEN
            ! <=== WARNING filling by blocks here works because those blocks are symmetric
            !=== Matrices in petsc are row-oriented
            CALL MatSetValues(this%matrices%dijH, nw, idx_loc, nw, idx_loc, mat_h_loc, ADD_VALUES, ierr)
         END IF
      END DO

      CALL MatAssemblyBegin(this%matrices_L%dijL, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd  (this%matrices_L%dijL, MAT_FINAL_ASSEMBLY, ierr)

      !===add value on diagonal
      CALL MatGetRowSum(this%matrices_L%dijL, this%x4vec, ierr)
      CALL VecScale(this%x4vec, -1.d0, ierr)
      CALL MatDiagonalSet(this%matrices_L%dijL, this%x4vec, INSERT_VALUES, ierr)
      !===add value on diagonal

      IF (this%method==METHOD_HIGH) THEN
         CALL MatAssemblyBegin(this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)
         CALL MatAssemblyEnd  (this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)

         !===add value on diagonal
         CALL MatGetRowSum(this%matrices%dijH, this%x4vec, ierr)
         CALL VecScale(this%x4vec, -1.d0, ierr)
         CALL MatDiagonalSet(this%matrices%dijH, this%x4vec, INSERT_VALUES, ierr)
         !===add value on diagonal

         IF (this%limiting%if_limiting .AND. this%which_high_method==1) THEN
            DO nl=1, this%limiting_all_functionals%nl
               CALL array_to_petsc_vec(bounds(:,nl), this%x1vec, this%LA, 'min')
               CALL extract_through_ghost(this%x1vec, 1, 1, this%LA, bounds(:, nl), opt_assemble=.FALSE.)
            END DO
         END IF
      END IF
   END SUBROUTINE compute_dij

   SUBROUTINE invert_mass(this, vec_rhs, vec_lhs, which_mass)
      USE my_util,          ONLY: error_petsc, to_str
      USE solver_petsc,     ONLY: solver
      USE compute_periodic, ONLY: periodic_vector_petsc
      IMPLICIT NONE
      CLASS(hyperbolic_type), INTENT(INOUT) :: this
      INTEGER, INTENT(IN) :: which_mass
      INTEGER :: it, ierr
      Vec     :: vec_rhs, vec_lhs

      SELECT CASE(which_mass)
      CASE(LUMPED_MASS)
         !=== vec_lhs <-- vec_rhs / lumped_mass
         CALL VecPointWiseDivide(vec_lhs, vec_rhs, this%matrices%lump_mass_vec, ierr)
      CASE(QUASI_CONSISTENT_MASS)
         CALL VecCopy(vec_rhs, vec_lhs, ierr)
         CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, vec_lhs, this%LA)
         !=== First division by lumped mass
         CALL VecCopy(vec_lhs, this%x5vec, ierr)
         !=== Loop on number of corrections
         DO it=1, this%nb_correction_mass
            !=== x5 <-- Ar@x5 
            CALL MatMult(this%matrices%Al_mass, this%x5vec, this%x4vec, ierr)
            CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x4vec, this%LA)
            CALL VecCopy(this%x4vec, this%x5vec, ierr)
            !=== vec_lhs <-- vec_lhs + x5
            CALL VecAYPX(vec_lhs, 1.d0, this%x5vec, ierr)
         END DO
         CALL VecPointWiseDivide(vec_lhs, vec_lhs, this%matrices%lump_mass_vec, ierr)
      CASE(CONSISTENT_MASS)
         CALL solver(this%matrices%ksp_consistent_mass, vec_rhs, vec_lhs, &
                     reinit = .FALSE., verbose = .FALSE.)
      CASE DEFAULT
         CALL error_petsc("BUG in hyperbolic%update, invert_mass => wrong value "&
         &//to_str(which_mass)//", should be 1, 2, 3")
      END SELECT
   END SUBROUTINE invert_mass

   SUBROUTINE commutator(this, un, alpha)
      USE space_dim
      USE sub_plot
      USE st_matrix, ONLY: extract_through_ghost
      USE fem_tn
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: un
      REAL(KIND = 8), DIMENSION(:), INTENT(OUT):: alpha
      REAL(KIND = 8), DIMENSION(this%mesh%np)  :: rk, rk_norm, eta, logeta
      INTEGER :: k, ierr, np_tot
      REAL(KIND = 8) :: norm_diff, norm_log
      CHARACTER(5) :: char
      LOGICAL :: if_commutator_H=.FALSE.
      PetscReal :: norm
      CALL VecGetSize(this%x5vec, np_tot, ierr)
      !===
      CALL VecSet(this%x4vec, 0.d0, ierr)
      CALL VecSet(this%x5vec, 0.d0, ierr)
      eta = this%eta_commute(un)
      logeta = log(abs(eta))
      norm_diff = 0.d0
      norm_log = 0.d0

      IF (if_commutator_H) THEN
         DO k = 1, k_dim
            CALL array_to_petsc_vec(logeta, this%x1vec, this%LA, 'insert') !<==v1 = log(eta)
            CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr) !<==v2 = dk(log(eta))

            CALL array_to_petsc_vec(eta,    this%x1vec, this%LA, 'insert') !<==v1 = eta
            CALL VecPointwiseMult(this%x3vec,  this%x1vec, this%x2vec, ierr) !<==v3 = eta*dk(log(eta))

            CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)          !<==v2 = dk(eta))
            CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr) !<==v3 = eta*dk(log(eta)) - dk(eta)

            CALL VecNorm(this%x3vec, Norm_1, norm, ierr)
            norm_diff = norm_diff + norm

            CALL VecNorm(this%x2vec, Norm_1, norm, ierr)
            norm_log = norm_log + norm

            CALL VecAbs(this%x3vec,ierr)
            CALL VecAXPY(this%x4vec, 1.d0, this%x3vec, ierr) !<==v4 = sum_k |dk(eta)-eta*dk(log(eta))|

            CALL VecAbs(this%x2vec,ierr)
            CALL VecAXPY(this%x5vec, 1.d0, this%x2vec, ierr) !<==v5 = sum_k |dk(eta)
         END DO
         CALL extract_through_ghost(this%x4vec, 1, 1, this%LA, rk, opt_assemble=.FALSE.)

         CALL extract_through_ghost(this%x5vec, 1, 1, this%LA, rk_norm, opt_assemble=.FALSE.)
      ELSE
         DO k = 1, k_dim
            CALL array_to_petsc_vec(logeta, this%x1vec, this%LA_L, 'insert') !<==v1 = log(eta)
            CALL MatMult(this%matrices_L%cij(k), this%x1vec, this%x2vec, ierr) !<==v2 = dk(log(eta))

            CALL array_to_petsc_vec(eta,    this%x1vec, this%LA_L, 'insert') !<==v1 = eta
            CALL VecPointwiseMult(this%x3vec,  this%x1vec, this%x2vec, ierr) !<==v3 = eta*dk(log(eta))

            CALL MatMult(this%matrices_L%cij(k), this%x1vec, this%x2vec, ierr)          !<==v2 = dk(eta))
            CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr) !<==v3 = eta*dk(log(eta)) - dk(eta)

            CALL VecNorm(this%x3vec, Norm_1, norm, ierr)
            norm_diff = norm_diff + norm

            CALL VecNorm(this%x2vec, Norm_1, norm, ierr)
            norm_log = norm_log + norm

            CALL VecAbs(this%x3vec,ierr)
            CALL VecAXPY(this%x4vec, 1.d0, this%x3vec, ierr) !<==v4 = sum_k |dk(eta)-eta*dk(log(eta))|

            CALL VecAbs(this%x2vec,ierr)
            CALL VecAXPY(this%x5vec, 1.d0, this%x2vec, ierr) !<==v5 = sum_k |dk(eta)
         END DO
      END IF

      CALL extract_through_ghost(this%x4vec, 1, 1, this%LA_L, rk, opt_assemble=.FALSE.)

      CALL extract_through_ghost(this%x5vec, 1, 1, this%LA_L, rk_norm, opt_assemble=.FALSE.)

      norm_log = norm_log/np_tot

      rk = abs(rk)/max(abs(rk_norm),1.d-1*norm_log,1d-30)
      alpha = MIN(1*rk,1.d0)
      alpha = threshold(alpha)

      ! CALL ns_l1_PAR(this%mesh, alpha, norm, this%mesh%comm)
      ! WRITE(*,*) 'norm commutator = ', norm

      !IF (this%time+1.1*this%dt>this%final_time .AND. stage==this%ERK%s+1) THEN
      IF (this%time+1.5*this%dt>this%final_time) THEN
         WRITE(char, '(I5)') this%mesh%rank
         CALL plot_scalar_field(this%mesh%jj, this%mesh%rr, alpha, 'a'//trim(adjustl(char))//'.plt')
         CALL plot_scalar_field(this%mesh%jj, this%mesh%rr, eta, 'eta'//trim(adjustl(char))//'.plt')
      END IF
   END SUBROUTINE commutator

   FUNCTION threshold(x) RESULT(g)
      IMPLICIT NONE
      INTEGER, PARAMETER :: exp=3
      REAL(KIND=8), DIMENSION(:)  :: x
      REAL(KIND=8), DIMENSION(SIZE(x))  :: z, t, zp, relu, f, g
      REAL(KIND=8), PARAMETER :: x0 = .5d0
      REAL(KIND=8), PARAMETER :: x1=SQRT(3.d0)*x0
      SELECT CASE(exp)
      CASE(2)
         !===Quadratic threshold
         z = x-x0
         zp = x-2*x0
         relu = (zp+ABS(zp))/2
         f = -z*(z**2-x1**2)  + relu*(z-x0)*(z+2*x0)
         g = (f + 2*x0**3)/(4*x0**3)
      CASE(3)
         !===Cubic threshold
         relu = ((x-2*x0)+ABS(x-2*x0))/2
         t = x/(2*x0)
         g = t**3*(10-15*t+6*t**2) - relu*(t-1)**2*(6*t**2+3*t+1)/(2*x0)
      END SELECT
      RETURN
   END FUNCTION threshold

   SUBROUTINE apply_limiting(this, un, bounds)
      USE space_dim
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8), DIMENSION(:,:),       INTENT(INOUT) :: un
      REAL(KIND = 8), DIMENSION(:,:),       INTENT(INOUT) :: bounds
real(kind = 8), DIMENSION(SIZE(bounds, 1)) :: loc_min_bis
      REAL(KIND = 8), DIMENSION(SIZE(un,1),SIZE(un,2))    :: un_temp
      INTEGER :: it, i

      DO i=1, size(bounds, 2)
!TESTTTTTT
!          loc_min_bis = bounds(:,i)
!          CALL this%limiting_all_functionals%RELAX_BOUNDS(un, loc_min_bis, i)
! ! WRITE(*,*) "after", loc_min
! ! WRITE(*,*) "MAXVAL DIFF", SUM(ABS(loc_min_bis-bounds(:,i)))/SUM(ABS(bounds(:,i)))
! bounds(:,i) = loc_min_bis 
!TESTTTTTT
         CALL this%limiting_all_functionals%RELAX_BOUNDS(un, bounds(:,i), i)
         DO it = 1, this%limiting%limit_max
            CALL this%limiting%iterative_cell_limiting_procedure(un,bounds(:,i),&
               this%limiting_all_functionals%limiting_functionals(i), un_temp)
            un(:,:) = un_temp(:,:)
         END DO
      END DO

   END SUBROUTINE apply_limiting

   SUBROUTINE init_vectors(this)
      USE st_matrix, ONLY : create_my_ghost
      USE petsc
#include "petsc/finclude/petsc.h"

      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      INTEGER, POINTER, DIMENSION(:) :: ifrom
      INTEGER :: n, ierr

      CALL create_my_ghost(this%mesh, this%LA, ifrom)
      CALL VecCreateGhost(this%communicator, this%mesh%dom_np, &
           PETSC_DETERMINE, SIZE(ifrom), ifrom, this%x1vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x2vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x3vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x4vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x5vec, ierr)

      CALL VecCreateSeq(PETSC_COMM_SELF, this%mesh%dom_np, this%vec_loc, ierr)

      ALLOCATE(this%tab(this%mesh%dom_np))
      DO n = 1, this%mesh%dom_np
         this%tab(n) = n - 1
      END DO

   END SUBROUTINE init_vectors

END MODULE abstract_hyperbolic_module
