MODULE abstract_hyperbolic_module

!>> limited global uses to avoid unexpected behaviors
#include "petsc/finclude/petsc.h"
   USE petsc
   USE petsc_tools,                          ONLY: array_to_petsc_vec
   USE Butcher_tableau
   USE hyperbolic_matrices_module,           ONLY: hyperbolic_matrices_type
   USE hyperbolic_bc_tools,                  ONLY: construct_udotn
   USE cell_limiting_engine_parallel_module, ONLY: limiting_type, limiting_functionals_type
   USE def_type_mesh,                        ONLY: mesh_type, petsc_csr_LA
   USE read_inputs_module,                   ONLY: rec_length
   USE space_dim,                            ONLY: k_dim
!>> limited global uses to avoid unexpected behaviors

   INTEGER, PRIVATE, PARAMETER :: LUMPED_MASS=1, QUASI_CONSISTENT_MASS=2, CONSISTENT_MASS=3
   CHARACTER(LEN=20), DIMENSION(3), PRIVATE, PARAMETER :: list_mass = &
               [CHARACTER(LEN=20) :: 'lumped', 'quasi_consistent', 'consistent']
   INTEGER, PRIVATE, PARAMETER :: METHOD_VISCOUS=1, METHOD_HIGH=2, METHOD_GALERKIN=3
   CHARACTER(LEN=20), DIMENSION(3), PRIVATE, PARAMETER  :: list_method = &
               [CHARACTER(LEN=20) :: 'viscous', 'high', 'galerkin']

   TYPE argument_hyperbolic_type
      CHARACTER(LEN=rec_length) :: CFL                     = '=== CFL ? ==='
      CHARACTER(LEN=rec_length) :: char_method             = '=== Which method to solve conservation equation (viscous,high,galerkin)? ==='
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
      CLASS(limiting_functionals_type), DIMENSION(:), POINTER :: limiting_functionals
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

      SUBROUTINE template_lambda(this, un, i, j, lambda_max)
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type),                               INTENT(INOUT) :: this
         REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(IN) :: un
         INTEGER,                                              INTENT(IN) :: i, j
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
      TYPE(limiting_functionals_type), DIMENSION(:), TARGET :: limiting_functionals

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
         CALL this%limiting%init(this%communicator, this%name, this%mesh_L, this%LA_L)
      ELSE
         CALL this%limiting%init(this%communicator, this%name, this%mesh, this%LA)
      END IF
      IF (this%limiting%if_limiting) THEN
         this%limiting_functionals => limiting_functionals
      ELSE
         ALLOCATE(this%limiting_functionals(0))
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
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      USE cell_limiting_engine_module
      USE sub_plot
      USE compute_periodic, ONLY : periodic_rhs_petsc, periodic_vector_petsc
      USE st_matrix, ONLY: extract_through_ghost
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s+1)  :: urk
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s)    :: flux_rk_at_dof
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim)                :: un_temp
      REAL(KIND = 8), DIMENSION(this%mesh%np)                               :: rk
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim, this%syst_dim)         :: flux_array
      REAL(KIND = 8), DIMENSION(this%mesh%np,SIZE(this%limiting_functionals)) :: bounds
      INTEGER, INTENT(IN) :: stage
      INTEGER :: comp, k, l, ierr, stage_prime, stage0
      REAL(KIND = 8) :: time_stage

      SELECT CASE(this%method)
      CASE(METHOD_GALERKIN)
         !===flux_array: flux at l=stage
         DO comp=1, this%syst_dim
            flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage-1)) !<== Notice: Flux at l=stage
         END DO

         !=== dij
         stage_prime = this%ERK%lp_of_l(stage) !<== Method should work with incremental ERK method
         IF (stage-1 .NE. stage_prime) THEN
            DO comp=1, this%syst_dim
               flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage_prime)) !<== Notice: Flux at l=stage'
            END DO
         END IF
         
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

            !=== flux_rk_at_dof(:,comp,stage-1): -sum_j(f(uj(stage)) cj) and store in flux_rk_at_dof(stage)
            CALL VecSet(this%x3vec, 0.d0, ierr)
            DO k = 1, k_dim
               !=== set flux_k in x1vec
               CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%mesh, this%LA, 'insert')
               !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
               CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
               !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
               CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
            END DO
            CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, flux_rk_at_dof(:,comp,stage-1), opt_assemble=.FALSE.) !<== Store sum_j(f(uj)cij) at l=stage
            !=== rk: sum_l MatRK_(s,l) f_l
            rk =0.d0
            DO l = 1, stage-1
               rk =  rk + this%ERK%MatRK(stage,l)*flux_rk_at_dof(:,comp,l)
            END DO

            !=== rk in x3vec
            CALL array_to_petsc_vec(rk, this%x2vec, this%mesh, this%LA, 'insert')
            CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x2vec, this%LA)

            !=== Inverting mass matrix and updating un with dt
            CALL this%invert_mass(this%x2vec, this%x3vec)

            !=== set un(comp) at l=stage0 in x1vec
            CALL array_to_petsc_vec(urk(:,comp,stage0), this%x1vec, this%mesh, this%LA, 'insert') !<== Notice un at l=stage0
            !=== x3 <-- dt*x3 + un (x1 <-- un a few lines above)
            CALL VecAYPX(this%x3vec, this%dt, this%x1vec, ierr)
            !=== Manually make un periodic and extract the result
            CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x3vec, this%LA)
            CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
         END DO

         !===Periodicity
         DO comp = 1, this%syst_dim
            CALL array_to_petsc_vec(urk(:,comp,stage), this%x1vec, this%mesh, this%LA, 'insert')
            CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x1vec, this%LA)
            CALL extract_through_ghost(this%x1vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
         END DO

         !===Boundary conditions
         CALL this%impose_bc(urk(:,:,stage), this%mesh, this%time)

      CASE(METHOD_VISCOUS)
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
               CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%mesh, this%LA, 'insert')
               !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
               CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
               !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
               CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
            END DO

            !=== set un(comp) in x1vec
            CALL array_to_petsc_vec(un_temp(:, comp), this%x1vec, this%mesh, this%LA, 'insert')
            !=== add dij un(comp)to x3vec in x2vec
            CALL MatMultAdd(this%matrices%dijL, this%x1vec, this%x3vec, this%x2vec, ierr)
            CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x2vec, this%LA)

            ! !=== x3 <-- x2 / lumped_mass by default in viscous method
            CALL this%invert_mass(this%x2vec, this%x3vec)

            !=== x3 <-- un + x3*dt   (x1 <-- un few lines above)
            CALL VecAYPX(this%x3vec, this%ERK%inc_C(stage)*this%dt, this%x1vec, ierr) !<==time step is ERK%inc_C(stage)*this%dt
            CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x3vec, this%LA)
            !=== un+1 <-- x3
            CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
         END DO

         CALL this%impose_bc(urk(:,:,stage), this%mesh, time_stage)

      CASE(METHOD_HIGH)
         !===flux_array: flux at l=stage
         DO comp=1, this%syst_dim
            flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage-1)) !<== Notice: Flux at l=stage
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

            !=== flux_rk_at_dof(:,comp,stage-1): -sum_j(f(uj(stage)) cj) and store in flux_rk_at_dof(stage)
            CALL VecSet(this%x3vec, 0.d0, ierr)
            DO k = 1, k_dim
               !=== set flux_k in x1vec
               CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%mesh, this%LA, 'insert')
               !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
               CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
               !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
               CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
            END DO
            CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, flux_rk_at_dof(:,comp,stage-1), opt_assemble=.FALSE.) !<== Store sum_j(f(uj)cij) at l=stage
            !=== rk: sum_l MatRK_(s,l) f_l
            rk =0.d0
            DO l = 1, stage-1
               rk =  rk + this%ERK%MatRK(stage,l)*flux_rk_at_dof(:,comp,l)
            END DO

            !=== rk in x3vec
            CALL array_to_petsc_vec(rk, this%x3vec, this%mesh, this%LA, 'insert')
            !=== set un(comp) at l=stage' in x1vec to compute viscous contribution
            CALL array_to_petsc_vec(urk(:,comp,stage_prime), this%x1vec, this%mesh, this%LA, 'insert') !<== l=stage'

            !=== add dij un(comp) to x3vec in x2vec
            CALL MatMultAdd(this%matrices%dijH, this%x1vec, this%x3vec, this%x2vec, ierr)
            CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x2vec, this%LA)
            !=== Inverting mass matrix and updating un with dt
            CALL this%invert_mass(this%x2vec, this%x3vec)

            !=== set un(comp) at l=stage0 in x1vec
            CALL array_to_petsc_vec(urk(:,comp,stage0), this%x1vec, this%mesh, this%LA, 'insert') !<== Notice un at l=stage0
            !=== x3 <-- dt*x3 + un (x1 <-- un a few lines above)
            CALL VecAYPX(this%x3vec, this%dt, this%x1vec, ierr)
            !=== Manually make un periodic and extract the result
            CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x3vec, this%LA)
            CALL extract_through_ghost(this%x3vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
         END DO
         !===Limiting
         IF (this%limiting%if_limiting) THEN
            CALL this%apply_limiting(urk(:,:,stage), bounds)
         END IF
         !===Periodicity
         DO comp = 1, this%syst_dim
            CALL array_to_petsc_vec(urk(:,comp,stage), this%x1vec, this%mesh, this%LA, 'insert')
            CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x1vec, this%LA)
            CALL extract_through_ghost(this%x1vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
            ! DO k = 1, this%mesh%per%nb_bords
            !    urk(this%mesh%per%list(k)%DIL, comp, stage) = urk(this%mesh%per%perlist(k)%DIL, comp, stage)
            ! END DO
         END DO

         !===Boundary conditions
         CALL this%impose_bc(urk(:,:,stage), this%mesh, this%time)
      CASE DEFAULT
         CALL error_petsc('Wrong method '//to_str(this%method)//' in '//TRIM(ADJUSTL(this%name))//&
         ' update, should be "viscous(1)" or "high(2)"')
      END SELECT

   END SUBROUTINE one_step_ERK

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
      REAL(KIND = 8) :: max_lambda
      REAL(KIND = 8), DIMENSION(this%syst_dim) :: uijbar
      LOGICAL, DIMENSION(this%mesh_L%medge) :: virgin_edge
      REAL(KIND = 8), DIMENSION(this%mesh_L%np)  :: alpha !<==commutator in (0,1)

      mesh => this%mesh_L
      LA   => this%LA_L

      !===Compute dijL
      CALL MatZeroEntries(this%matrices_L%dijL, ierr)
      IF (this%method==METHOD_HIGH) THEN
         CALL this%commutator(un, alpha)
         CALL MatZeroEntries(this%matrices%dijH, ierr)
         IF (this%limiting%if_limiting) THEN
            DO nl=1, Size(this%limiting_functionals)
               DO i=1, mesh%np
                  bounds(i, nl) = this%limiting_functionals(nl)%psi(un(i, :), 0.d0)
               END DO
            END DO
         END IF
      END IF


      virgin_edge = .TRUE.
      nw = mesh%gauss%n_w

      DO m = 1, mesh%me
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

               CALL this%compute_lambda(un, i, j, lambda_max)
               CALL MatGetValues(this%matrices_L%cij_norm_loc, 1, i_t - 1, 1, j_t - 1, norm_c, ierr)

               max_lambda = MAXVAL(lambda_max)
               dijL_c = max_lambda * norm_c

               IF (mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j
                     CALL this%compute_lambda(un, j, i, lambda_max)
                     CALL MatGetValues(this%matrices_L%cij_norm_loc, 1, j_t - 1, 1, i_t - 1, norm_c, ierr)

                     dijL_c = MAX(dijL_c, MAXVAL(lambda_max) * norm_c)
                     max_lambda = MAX(max_lambda,MAXVAL(lambda_max))
               END IF

               idx = LA%loc_to_glob(1, i) - 1
               jdx = LA%loc_to_glob(1, j) - 1

               CALL MatSetValues(this%matrices_L%dijL, 1, idx, 1, jdx, dijL_c, ADD_VALUES, ierr)
               CALL MatSetValues(this%matrices_L%dijL, 1, jdx, 1, idx, dijL_c, ADD_VALUES, ierr)

               CALL MatSetValues(this%matrices_L%dijL, 1, idx, 1, idx, -dijL_c, ADD_VALUES, ierr) !===add value on diagonal
               CALL MatSetValues(this%matrices_L%dijL, 1, jdx, 1, jdx, -dijL_c, ADD_VALUES, ierr) !===add value on diagonal
               IF (this%method==METHOD_HIGH) THEN
                  dijH_c = dijL_c*(alpha(i)+alpha(j))/2
                  CALL MatSetValues(this%matrices%dijH, 1, idx, 1, jdx, dijH_c, ADD_VALUES, ierr)
                  CALL MatSetValues(this%matrices%dijH, 1, jdx, 1, idx, dijH_c, ADD_VALUES, ierr)
                  CALL MatSetValues(this%matrices%dijH, 1, idx, 1, idx, -dijH_c, ADD_VALUES, ierr) !===add value on diagonal
                  CALL MatSetValues(this%matrices%dijH, 1, jdx, 1, jdx, -dijH_c, ADD_VALUES, ierr) !===add value on diagonal
                  !===Compute low-order update to estimate bounds
                  IF (this%limiting%if_limiting) THEN
                     DO k = 1, k_dim
                        CALL MatGetValues(this%matrices_L%nij_loc(k), 1, i_t - 1, 1, j_t - 1, nij_c(:, k), ierr)
                     END DO

                     DO comp=1, this%syst_dim
                        uijbar(comp) =  (un(i, comp)+un(j, comp))/2 - &
                        SUM((flux_array(j, :, comp)-flux_array(i, :, comp))*nij_c(1, :))/(2*max_lambda)
                     END DO

                     DO nl=1, SIZE(this%limiting_functionals)
                        bounds(i, nl) = MIN(bounds(i, nl), this%limiting_functionals(nl)%psi(uijbar, 0.d0))
                        bounds(j, nl) = MIN(bounds(j, nl), this%limiting_functionals(nl)%psi(uijbar, 0.d0))
                     END DO

                     !===End compute low-order update to estimate bounds
                  END IF
               END IF
            END IF
         END DO
      END DO
      CALL MatAssemblyBegin(this%matrices_L%dijL, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd  (this%matrices_L%dijL, MAT_FINAL_ASSEMBLY, ierr)

      IF (this%method==METHOD_HIGH) THEN
         CALL MatAssemblyBegin(this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)
         CALL MatAssemblyEnd  (this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)

         IF (this%limiting%if_limiting) THEN
            DO nl=1, SIZE(this%limiting_functionals)
               CALL array_to_petsc_vec(bounds(:,nl), this%x1vec, this%mesh, this%LA, 'min')
               CALL extract_through_ghost(this%x1vec, 1, 1, this%LA, bounds(:, nl), opt_assemble=.FALSE.)
            END DO
         END IF
      END IF
   END SUBROUTINE compute_dij

   SUBROUTINE invert_mass(this, vec_rhs, vec_lhs)
      USE my_util,          ONLY: error_petsc, to_str
      USE solver_petsc,     ONLY: solver
      USE compute_periodic, ONLY: periodic_vector_petsc
      IMPLICIT NONE
      CLASS(hyperbolic_type), INTENT(INOUT) :: this
      INTEGER :: it, ierr
      Vec     :: vec_rhs, vec_lhs

      SELECT CASE(this%which_mass)
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
         &//to_str(this%which_mass)//", should be 1, 2, 3")
      END SELECT
   END SUBROUTINE invert_mass

   SUBROUTINE commutator(this, un, alpha)
      USE space_dim
      USE sub_plot
      USE st_matrix, ONLY: extract_through_ghost
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: un
      REAL(KIND = 8), DIMENSION(:), INTENT(OUT):: alpha
      REAL(KIND = 8), DIMENSION(this%mesh%np)  :: rk, rk_norm, eta, logeta
      INTEGER :: k, ierr, np_tot
      REAL(KIND = 8) :: norm_diff, norm_log
      CHARACTER(5) :: char
      PetscReal :: norm
      CALL VecGetSize(this%x5vec, np_tot, ierr)
      !===
      CALL VecSet(this%x4vec, 0.d0, ierr)
      CALL VecSet(this%x5vec, 0.d0, ierr)
      !eta = pressure_from_state(this, un)/un(:,1)**1.4
      !eta = pressure_from_state(this, un)
      eta = this%eta_commute(un)
      ! eta = un(:,1)
      logeta = log(abs(eta))
      norm_diff = 0.d0
      norm_log = 0.d0

      ! DO k = 1, k_dim
      !    CALL array_to_petsc_vec(logeta, this%x1vec, this%mesh, this%LA, 'insert') !<==v1 = log(eta)
      !    CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr) !<==v2 = dk(log(eta))

      !    CALL array_to_petsc_vec(eta,    this%x1vec, this%mesh, this%LA, 'insert') !<==v1 = eta
      !    CALL VecPointwiseMult(this%x3vec,  this%x1vec, this%x2vec, ierr) !<==v3 = eta*dk(log(eta))

      !    CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)          !<==v2 = dk(eta))
      !    CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr) !<==v3 = eta*dk(log(eta)) - dk(eta)

      !    CALL VecNorm(this%x3vec, Norm_1, norm, ierr)
      !    norm_diff = norm_diff + norm

      !    CALL VecNorm(this%x2vec, Norm_1, norm, ierr)
      !    norm_log = norm_log + norm

      !    CALL VecAbs(this%x3vec,ierr)
      !    CALL VecAXPY(this%x4vec, 1.d0, this%x3vec, ierr) !<==v4 = sum_k |dk(eta)-eta*dk(log(eta))|

      !    CALL VecAbs(this%x2vec,ierr)
      !    CALL VecAXPY(this%x5vec, 1.d0, this%x2vec, ierr) !<==v5 = sum_k |dk(eta)
      ! END DO
      ! CALL extract_through_ghost(this%x4vec, 1, 1, this%LA, rk, opt_assemble=.FALSE.)

      ! CALL extract_through_ghost(this%x5vec, 1, 1, this%LA, rk_norm, opt_assemble=.FALSE.)
      DO k = 1, k_dim
         CALL array_to_petsc_vec(logeta, this%x1vec, this%mesh_L, this%LA_L, 'insert') !<==v1 = log(eta)
         CALL MatMult(this%matrices_L%cij(k), this%x1vec, this%x2vec, ierr) !<==v2 = dk(log(eta))

         CALL array_to_petsc_vec(eta,    this%x1vec, this%mesh_L, this%LA_L, 'insert') !<==v1 = eta
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

      CALL extract_through_ghost(this%x4vec, 1, 1, this%LA_L, rk, opt_assemble=.FALSE.)

      CALL extract_through_ghost(this%x5vec, 1, 1, this%LA_L, rk_norm, opt_assemble=.FALSE.)



      norm_log = norm_log/np_tot

      rk = abs(rk)/max(abs(rk_norm),1.d-1*norm_log)
      alpha = MIN(10*rk,1.d0)
      alpha = threshold(alpha)

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
      REAL(KIND = 8), DIMENSION(:,:),       INTENT(IN)    :: bounds
      REAL(KIND = 8), DIMENSION(SIZE(un,1),SIZE(un,2))    :: un_temp
      INTEGER :: it, i

      DO i=1, size(bounds, 2)
         DO it = 1, this%limiting%limit_max
            CALL this%limiting%iterative_cell_limiting_procedure(un,bounds(:,i),&
               this%limiting_functionals(i), un_temp)
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
