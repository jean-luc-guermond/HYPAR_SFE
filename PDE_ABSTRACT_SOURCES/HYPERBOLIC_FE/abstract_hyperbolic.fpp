#:set n_lim = n_lim_input
#:set syst_dim = syst_dim_input
MODULE template_abstract_hyperbolic_module

!>> limited global uses to avoid unexpected behaviors
#include "petsc/finclude/petsc.h"
   USE petsc
   USE petsc_tools,                          ONLY: array_to_petsc_vec
   USE Butcher_tableau
   USE hyperbolic_matrices_module,           ONLY: hyperbolic_matrices_type
   USE hyperbolic_bc_tools,                  ONLY: construct_udotn
   USE template_cell_limiting_engine_parallel_module, ONLY: limiting_type, limiting_functional_type, limiting_all_functional_type
   USE def_type_mesh,                        ONLY: mesh_type
   USE petsc_csr_LA_module,                  ONLY: petsc_csr_LA
   USE read_inputs_module,                   ONLY: rec_length
   USE space_dim,                            ONLY: k_dim
   ! use profiler_module, ONLY: profiler_type
!>> limited global uses to avoid unexpected behaviors

!==============================================!
!========= DEFINITIONS AND INTERFACES =========!
!==============================================!

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
      CHARACTER(LEN=rec_length) :: char_which_mass         = '=== Which mass matrix: lumped, quasi_consistent, consistent? ==='
      CHARACTER(LEN=rec_length) :: nb_correction_mass      = '=== For quasi_consistent mass, how many corrections? (0=lumped_mass) ==='
   END TYPE argument_hyperbolic_type

   TYPE, ABSTRACT, PUBLIC :: hyperbolic_type
      !===Parameters read from data
      REAL(KIND=8)                 :: CFL                     = 0.5d0
      CHARACTER(LEN=rec_length)    :: char_method             = 'viscous'
      INTEGER                      :: method                  = METHOD_VISCOUS
      INTEGER                      :: which_high_method       = 1
      INTEGER                      :: erk_sv                  = -31
      LOGICAL                      :: if_hybrid_mesh_limiting = .TRUE.
      CHARACTER(LEN=rec_length)    :: char_which_mass         = 'lumped'
      INTEGER                      :: which_mass              = LUMPED_MASS
      INTEGER                      :: nb_correction_mass      = 1

      !===Parameters built along way
      MPI_Comm :: communicator
      Vec          :: x1vec, x2vec, x2_ghost, x3vec
      Vec          :: x4vec, x5vec, x6vec !!!!!Conveniance vectors to be used only inside procedures!!!!
      Vec, DIMENSION(:,:), ALLOCATABLE :: flux_rk_at_dof
      CHARACTER(LEN=:), ALLOCATABLE :: name
      INTEGER :: syst_dim = ${syst_dim}$
      REAL(KIND = 8) :: dt, time, final_time
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: dij_diag
      INTEGER                              :: n_edge_in, n_edge_out
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: edge_in, edge_out
      TYPE(mesh_type),     POINTER :: mesh, mesh_L
      TYPE(petsc_csr_LA),  POINTER :: LA,   LA_L
      TYPE(BT),             PUBLIC :: ERK
      TYPE(hyperbolic_matrices_type) :: matrices
      TYPE(limiting_type)            :: limiting
      TYPE(limiting_all_functional_type) :: limiting_all_functionals
      PROCEDURE(interface_eta_commute),           NOPASS, POINTER :: eta_commute
      PROCEDURE(interface_compute_bounds_uijbar),         POINTER :: compute_bounds_uijbar => default_compute_bounds_uijbar
      ! type(profiler_type) :: profiler
!=== WORKSPACE ===!
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: x1, x2, x3, x4, x5, x6
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xx1, xx2, xx3, xx4, xx5, xx6
      REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: xxx1, xxx2, xxx3, xxx4, xxx5, xxx6

!=== WORKSPACE ===!

   CONTAINS
      PROCEDURE, PUBLIC   :: init_hyperbolic
      PROCEDURE, PUBLIC   :: set_times => set_times_hyperbolic
      PROCEDURE, PRIVATE  :: read_hyperbolic_data, init_vectors, init_edges_dij, init_F90_arrays
      PROCEDURE, PUBLIC   :: update, one_step_ERK
      PROCEDURE, PRIVATE  :: compute_dt, compute_dij
      PROCEDURE, PRIVATE  :: invert_mass, commutator, apply_limiting
      PROCEDURE(interface_flux),         DEFERRED :: flux
      PROCEDURE(interface_lambda),       DEFERRED :: compute_lambda
      PROCEDURE, NOPASS                           :: construct_udotn => construct_udotn
      PROCEDURE(interface_construct_bc), DEFERRED :: construct_bc
      PROCEDURE(interface_impose_bc),    DEFERRED :: impose_bc
   END TYPE hyperbolic_type

   ABSTRACT INTERFACE

      FUNCTION interface_eta_commute(un) RESULT(eta)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(:, :), INTENT(IN)    :: un
         REAL(KIND=8), DIMENSION(SIZE(un,1))         :: eta
      END FUNCTION interface_eta_commute

      SUBROUTINE interface_construct_bc(this, mesh, LA)
         USE def_type_mesh
         USE petsc_csr_LA_module
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type), INTENT(INOUT) :: this
         TYPE(mesh_type)           :: mesh
         TYPE(petsc_csr_LA)        :: LA
      END SUBROUTINE interface_construct_bc

      SUBROUTINE interface_impose_bc(this, un, mesh, time)
         USE def_type_mesh
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type),       INTENT(INOUT) :: this
         TYPE(mesh_type)                                :: mesh
         REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
         REAL(KIND = 8), INTENT(IN)                     :: time
      END SUBROUTINE interface_impose_bc

      SUBROUTINE interface_lambda(this, un, i, j, nij, lambda_max)
         USE space_dim
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type),                               INTENT(INOUT) :: this
         REAL(KIND=8), DIMENSION(this%mesh%np, ${syst_dim}$), INTENT(IN) :: un
         INTEGER,                                              INTENT(IN) :: i, j
         REAL(KIND=8), DIMENSION(k_dim),                       INTENT(IN) :: nij
         REAL(KIND=8), DIMENSION(2),                           INTENT(OUT) :: lambda_max
      END SUBROUTINE interface_lambda

      FUNCTION interface_flux(this, comp, un) RESULT(vv)
         USE space_dim
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type),          INTENT(INOUT) :: this
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
         INTEGER,                         INTENT(IN) :: comp
         REAL(KIND = 8), DIMENSION(SIZE(un, 1), k_dim) :: vv
      END FUNCTION interface_flux

      SUBROUTINE interface_compute_bounds_uijbar(this, flux_array, un, bounds)
         USE space_dim
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type) :: this
         REAL(KIND = 8), DIMENSION(this%mesh_L%np, k_dim, ${syst_dim}$), INTENT(IN) :: flux_array
         REAL(KIND = 8), DIMENSION(this%mesh_L%np, ${syst_dim}$),        INTENT(IN)  :: un
         REAL(KIND = 8), DIMENSION(this%mesh_L%np, ${n_lim}$),        INTENT(OUT) :: bounds
      END SUBROUTINE interface_compute_bounds_uijbar

   END INTERFACE

CONTAINS

!===================================!
!========= INITIALIZATIONS =========!
!===================================!

   SUBROUTINE init_hyperbolic(this, communicator, name, mesh, limiting_functionals)  
      USE my_util,            ONLY: error_petsc, to_str
      USE space_dim
      USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
      USE construct_mesh, ONLY: generate_boundary_structure, refine_mesh_1D2D_Pk2P1, create_gauss_points_1D_2D
      USE omp_lib

      IMPLICIT NONE
      CLASS(hyperbolic_type), INTENT(INOUT) :: this
      MPI_Comm,                   INTENT(IN) :: communicator
      CHARACTER(LEN=*),             INTENT(IN) :: name
      TYPE(mesh_type), TARGET,    INTENT(IN) :: mesh
      TYPE(limiting_functional_type), DIMENSION(:), TARGET :: limiting_functionals
      ! character(len=10), dimension(:), allocatable :: list_profiling

      this%name = name
      this%mesh => mesh
      this%communicator = communicator

      !===Construct LA
      ALLOCATE(this%LA)
      CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, this%LA)

      CALL this%read_hyperbolic_data("HYPERBOLIC PARAMETERS FOR "//trim(adjustl(this%name)))

      !=== Build ERK structure
      CALL this%ERK%init(this%erk_sv)

      !===Hybrid meshes
      IF (mesh%info%type_fe>1) THEN
         IF (mesh%rank==0) WRITE(*,*) "Building hybrid mesh for ", this%name
         !=== creating low order stencil
         ALLOCATE(this%mesh_L)
         ALLOCATE(this%LA_L)

         CALL refine_mesh_1D2D_Pk2P1(this%mesh, this%mesh_L)
         !=== generating associated P1 gauss points
         CALL create_gauss_points_1D_2D(this%mesh_L, 1)
         !=== adding surface elements (for Dirichlet and periodicity)
         CALL generate_boundary_structure(this%mesh_L)
         !=== generating petscSF structure
         CALL this%mesh_L%build_petscSF
         !=== Generating sparse structures
         CALL st_aij_csr_glob_block_with_extra_layer(this%communicator, 1, this%mesh_L, this%LA_L)
      ELSE
         this%mesh_L     => this%mesh
         this%LA_L       => this%LA
      END IF

      CALL this%LA_L%init_mat_loc_to_glob(this%mesh_L)
      ! CALL this%LA%init_mat_loc_to_glob(this%mesh)

      !===Matrices
      this%matrices%method     = this%method
      this%matrices%which_mass = this%which_mass
      CALL this%matrices%construct(this%communicator, this%mesh, this%LA, this%mesh_L, this%LA_L)

      !===Goshting structures
      CALL this%init_vectors
      CALL this%init_F90_arrays
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

      !===Edges for dij computations
      CALL this%init_edges_dij
      ALLOCATE(this%dij_diag(this%mesh_L%np))
      
      !===Profiler
      ! ALLOCATE(list_profiling(6))
      ! list_profiling(1) = 'tot'
      ! list_profiling(2) = 'other'
      ! list_profiling(3) = 'dij_riemann'
      ! list_profiling(4) = 'dij_reorder'
      ! list_profiling(5) = 'dij_matset'
      ! list_profiling(6) = 'limiting'
      ! CALL this%profiler%init(list_profiling)

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

      !===Method
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

   SUBROUTINE set_times_hyperbolic(this, times)
      IMPLICIT NONE
      CLASS(hyperbolic_type)                 :: this
      REAL(KIND=8), DIMENSION(2), INTENT(IN) :: times
      this%time = times(1) !<==initial_time
      this%final_time = times(2) !<==final_time
   END SUBROUTINE set_times_hyperbolic

   SUBROUTINE init_edges_dij(this)
      IMPLICIT NONE
      CLASS(hyperbolic_type)                :: this
      LOGICAL, DIMENSION(this%mesh_L%medge) :: virgin_edge
      INTEGER :: m, n

      ASSOCIATE(mesh => this%mesh_L)

      virgin_edge = .TRUE.
      this%n_edge_in = 0
      this%n_edge_out = 0
      
      !=== First loop to count the number of edges in and out
      DO m=1, mesh%me
         DO n=1, mesh%gauss%n_e
            IF (mesh%attr_e(mesh%jce(n, m))) THEN
               IF (.NOT. virgin_edge(mesh%jce_loc(n, m))) CYCLE
               virgin_edge(mesh%jce_loc(n, m)) = .FALSE.
               IF (mesh%side_edge(n, m)) THEN
                  this%n_edge_out = this%n_edge_out + 1
               ELSE
                  this%n_edge_in = this%n_edge_in + 1
               END IF
            END IF
         END DO
      END DO
      ALLOCATE(this%edge_in(2,this%n_edge_in))
      ALLOCATE(this%edge_out(2,this%n_edge_out))

      !=== Second loop to fill the edge_in and edge_out arrays
      virgin_edge = .TRUE.
      this%n_edge_in = 0
      this%n_edge_out = 0

      DO m=1, mesh%me
         DO n=1, mesh%gauss%n_e
            IF (mesh%attr_e(mesh%jce(n, m))) THEN
               IF (.NOT. virgin_edge(mesh%jce_loc(n, m))) CYCLE
               virgin_edge(mesh%jce_loc(n, m)) = .FALSE.
               IF (mesh%side_edge(n, m)) THEN
                  this%n_edge_out = this%n_edge_out + 1
                  this%edge_out(1,this%n_edge_out) = m
                  this%edge_out(2,this%n_edge_out) = n
               ELSE
                  this%n_edge_in = this%n_edge_in + 1
                  this%edge_in(1,this%n_edge_in) = m
                  this%edge_in(2,this%n_edge_in) = n
               END IF
            END IF
         END DO
      END DO

      END ASSOCIATE

   END SUBROUTINE init_edges_dij

   SUBROUTINE init_vectors(this)
      USE st_matrix, ONLY : create_my_ghost
      USE petsc
#include "petsc/finclude/petsc.h"

      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      INTEGER, POINTER, DIMENSION(:) :: ifrom
      INTEGER :: comp, stage, ierr

      CALL create_my_ghost(this%mesh, this%LA, ifrom)
      CALL VecCreateGhost(this%communicator, this%mesh%dom_np, &
           PETSC_DETERMINE, SIZE(ifrom), ifrom, this%x1vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x2vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x3vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x4vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x5vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x6vec, ierr)

      ALLOCATE(this%flux_rk_at_dof(${syst_dim}$, this%ERK%s+1))
      DO comp = 1, ${syst_dim}$
         DO stage=1, this%ERK%s+1
            CALL VecDuplicate(this%x1vec, this%flux_rk_at_dof(comp, stage), ierr)
         END DO
      END DO
   END SUBROUTINE init_vectors

   SUBROUTINE init_F90_arrays(this)

      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      INTEGER :: nw, me, np

      nw = this%mesh_L%gauss%n_w
      me = this%mesh_L%me
      np = this%mesh_L%np

      ALLOCATE(this%xxx1(nw, nw, me), this%xxx2(nw, nw, me))
      ALLOCATE(this%x1(np), this%x2(np))

   END SUBROUTINE init_F90_arrays

!===================================!
!========= TIMESTEPPING SUB ========!
!===================================!

   SUBROUTINE update(this,un_in)
     IMPLICIT NONE
     CLASS(hyperbolic_type)                                            :: this
     REAL(KIND=8), DIMENSION(this%mesh%np, ${syst_dim}$)               :: un_in
     REAL(KIND=8), DIMENSION(this%mesh%np, ${syst_dim}$, this%ERK%s+1) :: urk
     INTEGER  :: stage
     urk(:,:,1) = un_in
     DO stage = 2, this%ERK%s+1
        CALL this%one_step_ERK(stage,urk)
     END DO
     un_in = urk(:,:,this%ERK%s+1)
     this%time = this%time + this%dt
   END SUBROUTINE update

   SUBROUTINE one_step_ERK(this,stage,urk)
      USE my_util, ONLY : error_petsc, to_str
      IMPLICIT NONE
      CLASS(hyperbolic_type)                                               :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, ${syst_dim}$, this%ERK%s+1)  :: urk
      INTEGER, INTENT(IN) :: stage

      SELECT CASE(this%method)
      CASE(METHOD_GALERKIN)
         CALL one_step_ERK_galerkin(this,stage,urk)
      CASE(METHOD_VISCOUS)
         CALL one_step_ERK_viscous(this,stage,urk)
      CASE(METHOD_HIGH)
         CALL one_step_ERK_high(this,stage,urk)
      CASE DEFAULT
         CALL error_petsc("wrong method in one_step_ERK "//to_str(this%method))
      END SELECT
   END SUBROUTINE one_step_ERK

   SUBROUTINE one_step_ERK_galerkin(this,stage,urk)
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      USE cell_limiting_engine_module
      USE sub_plot
      USE compute_periodic, ONLY : periodic_rhs_petsc, periodic_vector_petsc
      USE st_matrix, ONLY: extract_through_ghost
      USE fem_tn
      IMPLICIT NONE
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, ${syst_dim}$, this%ERK%s+1)  :: urk
      REAL(KIND = 8), DIMENSION(this%mesh%np, ${syst_dim}$)                :: un_temp
      REAL(KIND = 8), DIMENSION(this%mesh%np)                               :: rk
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim, ${syst_dim}$)         :: flux_array
      REAL(KIND = 8), DIMENSION(this%mesh%np,this%limiting_all_functionals%nl) :: bounds
      INTEGER, INTENT(IN) :: stage
      INTEGER :: comp, k, l, ierr, stage_prime, stage0
      REAL(KIND = 8) :: time_stage

      !===flux_array: flux at l=stage
      DO comp=1, ${syst_dim}$
         flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage-1)) !<== Notice: Flux at l=stage
      END DO

      DO comp = 1, ${syst_dim}$
         !=== flux_rk_at_dof(:,comp,stage-1): -sum_j(f(uj(stage)) cj) and store in flux_rk_at_dof(stage)
         CALL VecZeroEntries(this%flux_rk_at_dof(comp, stage-1), ierr)
         DO k = 1, k_dim
            !=== set flux_k in x1vec
            CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%LA, 'insert')
            !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
            CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
            !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
            CALL VecAXPY(this%flux_rk_at_dof(comp, stage-1), -1.d0, this%x2vec, ierr)
         END DO
      END DO

      !=== dt
      stage_prime = this%ERK%lp_of_l(stage) !<== Method should work with incremental ERK method
      IF (stage==2) THEN !< ==Compute time step only once per ERK step
         CALL this%compute_dij(urk(:,:,stage_prime)) !<== Notice: State at l=stage'
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
      DO comp = 1, ${syst_dim}$

         !=== rk: sum_l MatRK_(s,l) f_l
         CALL VecZeroEntries(this%x2vec, ierr)
         CALL VecMAXPY(this%x2vec, stage-1, this%ERK%MatRK(stage,1:stage-1), this%flux_rk_at_dof(comp, 1:stage-1), ierr)

         !=== rk in x3vec
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
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, ${syst_dim}$, this%ERK%s+1)  :: urk
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, ${syst_dim}$)                :: un_temp
      REAL(KIND = 8), DIMENSION(this%mesh_L%np)                               :: rk
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, k_dim, ${syst_dim}$)         :: flux_array
      INTEGER, INTENT(IN) :: stage
      INTEGER :: comp, k, l, ierr, stage_prime, stage0
      REAL(KIND = 8) :: time_stage

      stage_prime = this%ERK%lp_of_l(stage) !<== Method should work with incremental ERK method
      un_temp = urk(:,:,stage_prime)      !<== Method should work with incremental ERK method
      DO comp=1, ${syst_dim}$
         flux_array(:, :, comp) = this%flux(comp, un_temp)
      END DO

      !===compute dijL (ONLY IF METHOD==VISCOUS) and dt
      IF (this%method==METHOD_VISCOUS) THEN
         CALL this%compute_dij(un_temp)
      END IF
      IF (stage==2) THEN !< ==Compute time step
         CALL this%compute_dt
         ! IF (this%time+this%dt.GE.this%final_time) THEN
         !   this%dt = this%final_time-this%time
         ! END IF
      END IF
      time_stage = this%time+this%ERK%C(stage)*this%dt !<== Wait for dt to be computed

      DO comp = 1, ${syst_dim}$
         CALL VecZeroEntries(this%x3vec, ierr)
         DO k = 1, k_dim
            !=== set flux_k in x1vec
            CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%LA_L, 'insert')
            !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
            CALL MatMult(this%matrices%cijL(k), this%x1vec, this%x2vec, ierr)
            !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
            CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
         END DO

         !=== set un(comp) in x1vec
         CALL array_to_petsc_vec(un_temp(:, comp), this%x1vec, this%LA_L, 'insert')
         !=== add dij un(comp)to x3vec in x2vec
         CALL MatMultAdd(this%matrices%dijL, this%x1vec, this%x3vec, this%x2vec, ierr)
         ! CALL MatMultAdd(this%matrices_L%dijL, this%x1vec, this%x3vec, this%x2vec, ierr)
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

   SUBROUTINE one_step_ERK_high(this,stage,urk)
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      USE cell_limiting_engine_module
      USE sub_plot
      USE compute_periodic, ONLY : periodic_rhs_petsc, periodic_vector_petsc
      USE st_matrix, ONLY: extract_through_ghost
      USE fem_tn
      IMPLICIT NONE
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, ${syst_dim}$, this%ERK%s+1)  :: urk
      REAL(KIND = 8), DIMENSION(this%mesh%np, ${syst_dim}$)                :: un_temp, u_L
      REAL(KIND = 8), DIMENSION(this%mesh%np)                               :: rk, zero_vec, temp_bounds, min_bounds
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim, ${syst_dim}$)         :: flux_array
      REAL(KIND = 8), DIMENSION(this%mesh%np,this%limiting_all_functionals%nl) :: bounds
      INTEGER, INTENT(IN) :: stage
      INTEGER :: comp, k, l, ierr, stage_prime, stage0, i, m, me, ni, nl, nw
      REAL(KIND = 8) :: time_stage, min_temp
      !===flux_array: flux at l=stage
      DO comp=1, ${syst_dim}$
         flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage-1)) !<== Notice: Flux at l=stage
      END DO

      DO comp = 1, ${syst_dim}$
         !=== flux_rk_at_dof(:,comp,stage-1): -sum_j(f(uj(stage)) cj) and store in flux_rk_at_dof(stage)
         CALL VecZeroEntries(this%flux_rk_at_dof(comp, stage-1), ierr)
         DO k = 1, k_dim
            !=== set flux_k in x1vec
            CALL array_to_petsc_vec(flux_array(:, k, comp), this%x1vec, this%LA, 'insert')
            !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
            CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
            !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
            CALL VecAXPY(this%flux_rk_at_dof(comp, stage-1), -1.d0, this%x2vec, ierr)
         END DO
      END DO
      !=== dij
      stage_prime = this%ERK%lp_of_l(stage) !<== Method should work with incremental ERK method

      CALL this%compute_dij(urk(:,:,stage_prime)) !<== Notice: State at l=stage'
      IF (this%limiting%if_limiting) THEN
         SELECT CASE(this%which_high_method)
         CASE(1)
            IF (stage-1 .NE. stage_prime) THEN
               DO comp=1, ${syst_dim}$
                  flux_array(:, :, comp) = this%flux(comp, urk(:,:,stage_prime)) !<== Notice: Flux at l=stage'
               END DO
            END IF
            CALL this%compute_bounds_uijbar(flux_array, urk(:,:,stage_prime), bounds)
         CASE(2)
            CALL one_step_ERK_viscous(this,stage,urk)
            u_L = urk(:, :, stage)
            zero_vec = 0.d0
            DO nl=1, this%limiting_all_functionals%nl
               temp_bounds = this%limiting_all_functionals%limiting_functionals(nl)%psi(u_L(:, :), zero_vec)
               
               min_bounds = temp_bounds
               me = this%mesh_L%me
               nw = this%mesh_L%gauss%n_w
               DO m = 1, me
                  min_temp = MINVAL(temp_bounds(this%mesh_L%jj(:,m)))
                  DO ni = 1, nw
                     i = this%mesh_L%jj(ni,m)
                     min_bounds(i) = MIN(min_bounds(i), min_temp)
                  END DO
               END DO
               CALL this%mesh_L%reduce_through_ghost(min_bounds, MPI_MIN)
               bounds(:, nl) = min_bounds
            END DO
         CASE DEFAULT
            CALL error_petsc('BUG in one_step_ERK_high: wrong method for bounds '//to_str(this%which_high_method))
         END SELECT
      END IF

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
      DO comp = 1, ${syst_dim}$
         !=== rk: sum_l MatRK_(s,l) f_l
         CALL VecZeroEntries(this%x3vec, ierr)
         CALL VecMAXPY(this%x3vec, stage-1, this%ERK%MatRK(stage,1:stage-1), this%flux_rk_at_dof(comp, 1:stage-1), ierr)

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
      DO comp = 1, ${syst_dim}$
         CALL array_to_petsc_vec(urk(:,comp,stage), this%x1vec, this%LA, 'insert')
         CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x1vec, this%LA)
         CALL extract_through_ghost(this%x1vec, 1, 1, this%LA, urk(:, comp, stage), opt_assemble=.FALSE.)
      END DO

      !===Boundary conditions
      CALL this%impose_bc(urk(:,:,stage), this%mesh, this%time)
   END SUBROUTINE one_step_ERK_high

!=======================================!
!========= CONVENIANCE SUBROUTINES =====!
!=======================================!

   SUBROUTINE compute_dt(this)
      USE my_util
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8)         :: dt_min_glob
      INTEGER                :: ierr

      CALL VecZeroEntries(this%x1vec, ierr)
      CALL array_to_petsc_vec(this%dij_diag, this%x1vec, this%LA_L, 'add')
      CALL VecAbs(this%x1vec, ierr)
      CALL VecPointWiseDivide(this%x2vec, this%matrices%lump_mass_vec_L, this%x1vec, ierr)
      CALL VecMin(this%x2vec, PETSC_NULL_INTEGER, dt_min_glob, ierr)

      this%dt = this%ERK%s * this%CFL * dt_min_glob / 2.d0
      !===Notice rescale of time step with this%ERK%s

   END SUBROUTINE compute_dt

   SUBROUTINE compute_dij(this, un)
      USE space_dim
      USE petsc
      USE def_type_mesh
      USE compute_periodic
      USE st_matrix, ONLY: extract_through_ghost
      USE my_util
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, ${syst_dim}$),       INTENT(IN)  :: un
      INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge, nl, comp, edge_idx, idx
      REAL(KIND = 8), DIMENSION(k_dim)    :: nij_c_scal
      REAL(KIND = 8), DIMENSION(2) :: lambda_max
      INTEGER,        DIMENSION(this%mesh_L%gauss%n_w) :: idx_loc
      REAL(KIND = 8) :: max_lambda, norm_c_scal, dijL_c_scal, dijH_c_scal

      
      ASSOCIATE(mesh => this%mesh_L, LA => this%LA_L, &
         lim_func => this%limiting_all_functionals%limiting_functionals,&
         mat_L_glob => this%xxx1, max_lambda_array => this%xxx2,&  ! DIMENSION(this%mesh_L%gauss%n_w,this%mesh_L%gauss%n_w, this%mesh_L%me)
         alpha => this%x2)                    ! DIMENSION(this%mesh_L%np)  !<==commutator in (0,1)
      nw = mesh%gauss%n_w
      IF (this%method==METHOD_HIGH) THEN
         CALL this%commutator(un, alpha, 1)
      END IF

      mat_L_glob = 0.d0
      this%dij_diag = 0.d0
      !================================================!
      !=== Loop on inner edges (cij symmetric here) ===!
      !================================================!
      DO edge_idx=1, this%n_edge_in
         m = this%edge_in(1, edge_idx)
         n = this%edge_in(2, edge_idx)

         ni = MOD(n, nw) + 1
         nj = MOD(n + 1, nw) + 1
         i = mesh%jj(ni, m)
         j = mesh%jj(nj, m)

         nij_c_scal = this%matrices%nijL_loc_array(ni,nj,:,m)
         norm_c_scal = this%matrices%cijL_norm_loc_array(ni,nj,m)

         CALL this%compute_lambda(un, i, j, nij_c_scal, lambda_max)
         max_lambda = MAXVAL(lambda_max)
         dijL_c_scal = max_lambda * norm_c_scal

         mat_L_glob(ni, nj, m) = dijL_c_scal
         mat_L_glob(nj, ni, m) = dijL_c_scal
         max_lambda_array(ni, nj, m) = max_lambda
         max_lambda_array(nj, ni, m) = max_lambda

         this%dij_diag(i) = this%dij_diag(i) - dijL_c_scal
         this%dij_diag(j) = this%dij_diag(j) - dijL_c_scal

      END DO


      !=======================================================!
      !=== Loop on boundary edges (cij not symmetric here) ===!
      !=======================================================!
      DO edge_idx=1, this%n_edge_out
         m = this%edge_out(1, edge_idx)
         n = this%edge_out(2, edge_idx)

         ni = MOD(n, nw) + 1
         nj = MOD(n + 1, nw) + 1
         i = mesh%jj(ni, m)
         j = mesh%jj(nj, m)

         !=== (i,j) term
         nij_c_scal = this%matrices%nijL_loc_array(ni,nj,:,m)
         norm_c_scal = this%matrices%cijL_norm_loc_array(ni,nj,m)
         CALL this%compute_lambda(un, i, j, nij_c_scal, lambda_max)
         max_lambda = MAXVAL(lambda_max)
         dijL_c_scal = MAXVAL(lambda_max) * norm_c_scal

         !=== (j,i) term
         nij_c_scal = this%matrices%nijL_loc_array(nj,ni,:,m)
         norm_c_scal = this%matrices%cijL_norm_loc_array(nj,ni,m)
         CALL this%compute_lambda(un, j, i, nij_c_scal, lambda_max)

         !=== Cumulate dijL_c_scal
         dijL_c_scal = MAX(dijL_c_scal, MAXVAL(lambda_max) * norm_c_scal)
         mat_L_glob(ni, nj, m) = dijL_c_scal
         mat_L_glob(nj, ni, m) = dijL_c_scal

         max_lambda = MAX(MAXVAL(lambda_max), max_lambda)

         max_lambda_array(ni, nj, m) = max_lambda
         max_lambda_array(nj, ni, m) = max_lambda

         this%dij_diag(i) = this%dij_diag(i) - dijL_c_scal
         this%dij_diag(j) = this%dij_diag(j) - dijL_c_scal

      END DO

      !=================================================!
      !=== Computation of dijL on low  order stencil ===!
      !=================================================!
!=== avoid constructing matrices when possible
      IF ((this%method == METHOD_VISCOUS .OR. (this%method == METHOD_HIGH .AND. this%which_high_method==2))) THEN
         ASSOCIATE(arr_L => LA%zz_contig_1, mat_loc_to_glob => LA%mat_loc_to_glob)
            arr_L(:) = 0.d0
            DO m=1, mesh%me
               DO ni=1, nw
                     DO nj=1, nw
                        idx = mat_loc_to_glob(ni, nj, m)
                        arr_L(idx) = arr_L(idx) + mat_L_glob(nj, ni, m) !=== WARNING ROW ORIENTATION (not a pb here because dij is symmetric)
                     END DO
               END DO
            END DO
            CALL MatZeroEntries(this%matrices%dijL, ierr)
            CALL LA%fill_mat   (this%matrices%dijL, arr_L)
         END ASSOCIATE

         !===add value on diagonal
         CALL MatGetRowSum(this%matrices%dijL, this%x4vec, ierr)
         CALL VecScale(this%x4vec, -1.d0, ierr)
         CALL MatDiagonalSet(this%matrices%dijL, this%x4vec, INSERT_VALUES, ierr)
         !===add value on diagonal
      END IF

      !=================================================!
      !=== Computation of dijH on high order stencil ===!
      !=================================================!
      IF (this%method==METHOD_HIGH) THEN
         ASSOCIATE(arr_H => LA%zz_contig_2, mat_loc_to_glob => LA%mat_loc_to_glob)
         arr_H(:) = 0.d0
         CALL this%commutator(un, alpha, 2)
         DO m=1, mesh%me
            DO ni=1, nw
               i = mesh%jj(ni, m)
               DO nj=1, nw
                  j = mesh%jj(nj, m)
                  idx = mat_loc_to_glob(ni, nj, m)
                  arr_H(idx) = arr_H(idx) + mat_L_glob(nj, ni, m)*(alpha(i)+alpha(j))/2.d0 !=== WARNING ROW ORIENTATION (not a pb here because dij is symmetric)
               END DO
            END DO
         END DO
         CALL MatZeroEntries(this%matrices%dijH, ierr)
         CALL LA%fill_mat   (this%matrices%dijH, arr_H)
         END ASSOCIATE
         
         !===add value on diagonal
         CALL MatGetRowSum(this%matrices%dijH, this%x4vec, ierr)
         CALL VecScale(this%x4vec, -1.d0, ierr)
         CALL MatDiagonalSet(this%matrices%dijH, this%x4vec, INSERT_VALUES, ierr)
         !===add value on diagonal

      END IF

      END ASSOCIATE

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

   SUBROUTINE commutator(this, un, alpha, int_assembly)
      USE space_dim
      USE sub_plot
      USE st_matrix, ONLY: extract_through_ghost
      USE fem_tn
      USE my_util
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: un
      REAL(KIND = 8), DIMENSION(:), INTENT(OUT):: alpha
      INTEGER,                      INTENT(IN) :: int_assembly
      REAL(KIND = 8), DIMENSION(this%mesh%np)  :: rk, rk_norm, eta, logeta
      REAL(KIND = 8), DIMENSION(:), POINTER :: dummy
      INTEGER :: k, ierr, np_tot
      REAL(KIND = 8) :: norm_diff, norm_log
      CHARACTER(5) :: char
      LOGICAL :: if_commutator_H=.FALSE.
      PetscReal :: norm

      SELECT CASE(int_assembly)
      CASE(1)
         
         eta = this%eta_commute(un)
         CALL VecSetValues(this%x1vec, SIZE(eta), this%LA_L%loc_to_glob(1,:)-1, eta, INSERT_VALUES, ierr)
         CALL VecAssemblyBegin(this%x1vec, ierr)

         logeta = log(abs(eta))
         CALL VecSetValues(this%x2vec, SIZE(logeta), this%LA_L%loc_to_glob(1,:)-1, logeta, INSERT_VALUES, ierr)
         CALL VecAssemblyBegin(this%x2vec, ierr)
      CASE(2)
         CALL VecAssemblyEnd(this%x1vec, ierr)
         CALL VecAssemblyEnd(this%x2vec, ierr)
         
         CALL VecGetSize(this%x5vec, np_tot, ierr)
         !===
         CALL VecZeroEntries(this%x5vec, ierr)
         CALL VecZeroEntries(this%x6vec, ierr)
         norm_diff = 0.d0
         norm_log = 0.d0

         DO k = 1, k_dim
            CALL MatMult(this%matrices%cijL(k), this%x2vec, this%x3vec, ierr) !<==v3 = dk(log(eta))

            CALL VecPointwiseMult(this%x3vec,  this%x1vec, this%x3vec, ierr) !<==v3 = eta*dk(log(eta))

            CALL MatMult(this%matrices%cijL(k), this%x1vec, this%x4vec, ierr)          !<==v4 = dk(eta))
            CALL VecAXPY(this%x3vec, -1.d0, this%x4vec, ierr) !<==v3 = eta*dk(log(eta)) - dk(eta)

            CALL VecNorm(this%x3vec, Norm_1, norm, ierr)
            norm_diff = norm_diff + norm

            CALL VecNorm(this%x4vec, Norm_1, norm, ierr)
            norm_log = norm_log + norm

            CALL VecAbs(this%x3vec,ierr)
            CALL VecAXPY(this%x5vec, 1.d0, this%x3vec, ierr) !<==v5 = sum_k |dk(eta)-eta*dk(log(eta))|

            CALL VecAbs(this%x4vec,ierr)
            CALL VecAXPY(this%x6vec, 1.d0, this%x4vec, ierr) !<==v6 = sum_k |dk(eta)|
         END DO
      
      
         norm_log = norm_log/np_tot

         CALL VecGhostGetLocalForm(this%x6vec, this%x2_ghost, ierr)
         CALL VecGetArrayF90(this%x2_ghost, dummy, ierr)
         dummy = max(dummy,1.d-1*norm_log,1d-30)
         CALL VecRestoreArrayF90(this%x2_ghost, dummy, ierr)
         CALL VecGhostUpdateBegin(this%x6vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(this%x6vec, INSERT_VALUES, SCATTER_FORWARD, ierr)

         CALL VecPointWiseDivide(this%x5vec, this%x5vec, this%x6vec, ierr)
         CALL extract_through_ghost(this%x5vec, 1, 1, this%LA_L, rk, opt_assemble=.FALSE.)

         alpha = MIN(1*rk,1.d0)
         alpha = threshold(alpha)

         !IF (this%time+1.1*this%dt>this%final_time .AND. stage==this%ERK%s+1) THEN
         IF (this%time+1.5*this%dt>this%final_time) THEN
            WRITE(char, '(I5)') this%mesh%rank
            CALL plot_scalar_field(this%mesh%jj, this%mesh%rr, alpha, 'a'//trim(adjustl(char))//'.plt')
            CALL plot_scalar_field(this%mesh%jj, this%mesh%rr, eta, 'eta'//trim(adjustl(char))//'.plt')
         END IF

      CASE DEFAULT
         CALL error_petsc('BUG in commutator => wrong case '//to_str(int_assembly))
      END SELECT

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


!===================================================!
!=== Computation of bounds for limiting method 1 ===!
!===================================================!

   SUBROUTINE default_compute_bounds_uijbar(this, flux_array, un, bounds)
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, k_dim, ${syst_dim}$), INTENT(IN) :: flux_array
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, ${syst_dim}$),        INTENT(IN)  :: un
      REAL(KIND = 8), DIMENSION(this%mesh_L%np, ${n_lim}$), INTENT(OUT) :: bounds
      REAL(KIND = 8), DIMENSION(1, ${syst_dim}$) :: uijbar  !<=== FIXME
      REAL(KIND = 8), DIMENSION(1), PARAMETER   :: zero=0.d0
      REAL(KIND = 8), DIMENSION(k_dim)    :: nij_c_scal
      REAL(KIND = 8)                      :: norm_c_scal, max_lambda
      INTEGER :: nl, m, n, ni, nj, i, j, edge_idx, nw, comp

      ASSOCIATE(mesh => this%mesh_L,&
         lim_func => this%limiting_all_functionals%limiting_functionals, &
         max_lambda_array => this%xxx2, zero_vec => this%x1) ! DIMENSION(this%mesh_L%gauss%n_w,this%mesh_L%gauss%n_w, this%mesh_L%me)
         
      nw = mesh%gauss%n_w
      zero_vec = 0.d0
      DO nl=1, this%limiting_all_functionals%nl
         bounds(:, nl) = lim_func(nl)%psi(un(:, :), zero_vec)
      END DO

      DO edge_idx=1, this%n_edge_in
         m = this%edge_in(1, edge_idx)
         n = this%edge_in(2, edge_idx)

         ni = MOD(n, nw) + 1
         nj = MOD(n + 1, nw) + 1
         i = mesh%jj(ni, m)
         j = mesh%jj(nj, m)
         nij_c_scal(:) = this%matrices%nijL_loc_array(ni,nj,:,m)
         norm_c_scal = this%matrices%cijL_norm_loc_array(ni,nj,m)

         max_lambda = max_lambda_array(ni, nj, m)

         DO comp=1, ${syst_dim}$
            uijbar(1,comp) =  (un(i, comp)+un(j, comp))/2.d0 - &
            SUM((flux_array(j, :, comp)-flux_array(i, :, comp))*nij_c_scal(:))/(2.d0*max_lambda)
         END DO

         DO nl=1, this%limiting_all_functionals%nl
            bounds(i:i, nl) = MIN(bounds(i, nl), lim_func(nl)%psi(uijbar, zero))
            bounds(j:j, nl) = MIN(bounds(j, nl), lim_func(nl)%psi(uijbar, zero))
         END DO
      END DO

      DO edge_idx=1, this%n_edge_out
         m = this%edge_out(1, edge_idx)
         n = this%edge_out(2, edge_idx)

         ni = MOD(n, nw) + 1
         nj = MOD(n + 1, nw) + 1
         i = mesh%jj(ni, m)
         j = mesh%jj(nj, m)
         nij_c_scal(:) = this%matrices%nijL_loc_array(ni,nj,:,m)
         norm_c_scal = this%matrices%cijL_norm_loc_array(ni,nj,m)

         max_lambda = max_lambda_array(ni, nj, m)

         DO comp=1, ${syst_dim}$
            uijbar(1,comp) =  (un(i, comp)+un(j, comp))/2.d0 - &
            SUM((flux_array(j, :, comp)-flux_array(i, :, comp))*nij_c_scal(:))/(2.d0*max_lambda)
         END DO

         DO nl=1, this%limiting_all_functionals%nl
            bounds(i:i, nl) = MIN(bounds(i, nl), lim_func(nl)%psi(uijbar, zero))
            bounds(j:j, nl) = MIN(bounds(j, nl), lim_func(nl)%psi(uijbar, zero))
         END DO
      END DO

      DO nl=1, this%limiting_all_functionals%nl
         CALL mesh%reduce_through_ghost(bounds(:,nl), MPI_MIN)
      END DO

      END ASSOCIATE
   END SUBROUTINE default_compute_bounds_uijbar


   SUBROUTINE apply_limiting(this, un, bounds)
      USE space_dim
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8), DIMENSION(:,:),       INTENT(INOUT) :: un
      REAL(KIND = 8), DIMENSION(:,:),       INTENT(INOUT) :: bounds
      REAL(KIND = 8), DIMENSION(SIZE(un,1),SIZE(un,2))    :: un_temp
      INTEGER :: it, i

      DO i=1, size(bounds, 2)
         !TESTTTTTT
         !loc_min_bis = bounds(:,i)
         !CALL this%limiting_all_functionals%RELAX_BOUNDS(un, loc_min_bis, i)
         !WRITE(*,*) "after", loc_min
         !WRITE(*,*) "MAXVAL DIFF", SUM(ABS(loc_min_bis-bounds(:,i)))/SUM(ABS(bounds(:,i)))
         !bounds(:,i) = loc_min_bis 
         !TESTTTTTT

         CALL this%limiting_all_functionals%RELAX_BOUNDS(un, bounds(:,i), i)
         DO it = 1, this%limiting%limit_max
            ASSOCIATE(lim_func => this%limiting_all_functionals%limiting_functionals(i))
               IF (ASSOCIATED(lim_func%spe_iterative_cell_limiting_procedure)) THEN
                  CALL lim_func%spe_iterative_cell_limiting_procedure(this%limiting, un, bounds(:,i))
               ELSE
                  CALL this%limiting%iterative_cell_limiting_procedure(un, bounds(:,i), lim_func, un_temp)
                  un(:,:) = un_temp(:,:)
               END IF
            END ASSOCIATE
         END DO
      END DO

   END SUBROUTINE apply_limiting

END MODULE template_abstract_hyperbolic_module
