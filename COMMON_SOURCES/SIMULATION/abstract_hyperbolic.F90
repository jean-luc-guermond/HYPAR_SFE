MODULE abstract_hyperbolic_module

!>> limited global uses to avoid unexpected behaviors
#include "petsc/finclude/petsc.h"
   USE petsc
   USE petsc_tools,     ONLY : array_to_petsc_vec
   USE Butcher_tableau
   USE hyperbolic_matrices_module, ONLY: hyperbolic_matrices_type, init_my_vectors
   USE hyperbolic_bc_tools, ONLY: construct_udotn
   USE cell_limiting_engine_parallel_module, ONLY : limiting_type, limiting_bounds_type
   USE def_type_mesh, ONLY : mesh_type, petsc_csr_LA
   USE read_inputs_module,    ONLY : rec_length
   USE space_dim, ONLY: k_dim
!>> limited global uses to avoid unexpected behaviors


   TYPE argument_hyperbolic_type
      CHARACTER(LEN=rec_length) :: CFL    = '=== CFL ? ==='
      CHARACTER(LEN=rec_length) :: method = '=== Which method to solve Euler (viscous, high) ? ==='
      CHARACTER(LEN=rec_length) :: erk_sv = '=== ERK ? ==='
   END TYPE argument_hyperbolic_type

   TYPE, ABSTRACT :: hyperbolic_type
      !===Parameters read from data
      REAL(KIND=8)                 :: CFL       = 0.5d0
      CHARACTER(LEN=rec_length)    :: method = 'viscous'
      INTEGER                      :: erk_sv    = -21
      !===Parameters built along way
      MPI_Comm :: communicator
      Vec, POINTER :: x1vec, x2vec, x2_ghost, vec_loc
      Vec          :: x3vec, x4vec, x5vec
      CHARACTER(LEN=:), ALLOCATABLE :: name
      INTEGER                       :: syst_dim
      REAL(KIND = 8) :: dt, time, final_time
      INTEGER, DIMENSION(:), ALLOCATABLE :: tab
      TYPE(mesh_type),     POINTER :: mesh
      TYPE(petsc_csr_LA),  POINTER :: LA
      TYPE(BT),             PUBLIC :: ERK
      TYPE(hyperbolic_matrices_type) :: matrices
      TYPE(limiting_type)            :: limiting
      CLASS(limiting_bounds_type), POINTER :: limiting_bounds => NULL()
   CONTAINS
      PROCEDURE, PUBLIC   :: init_hyperbolic
      PROCEDURE, PRIVATE  :: read_hyperbolic_data
      PROCEDURE, PUBLIC   :: update
      PROCEDURE, PRIVATE  :: compute_dij, compute_dt, commutator
      PROCEDURE(template_flux),         DEFERRED :: flux
      PROCEDURE(template_lambda),       DEFERRED :: compute_lambda
      PROCEDURE, NOPASS                          :: construct_udotn => construct_udotn
      PROCEDURE(template_construct_bc), DEFERRED :: construct_bc
      PROCEDURE(template_impose_bc),    DEFERRED :: impose_bc
      ! PROCEDURE, PRIVATE :: compute_dK, compute_dt_from_dK
   END TYPE hyperbolic_type


   ABSTRACT INTERFACE
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

      SUBROUTINE template_lambda(this, un, i, j, lambda_max, on_edge)
         IMPORT :: hyperbolic_type
         IMPLICIT NONE
         CLASS(hyperbolic_type),                               INTENT(INOUT) :: this
         REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(IN) :: un
         INTEGER,                                              INTENT(IN) :: i, j
         REAL(KIND=8), DIMENSION(2),                           INTENT(OUT) :: lambda_max
         LOGICAL,                                              INTENT(IN) :: on_edge
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

   SUBROUTINE init_hyperbolic(this, communicator, name, mesh, LA, times, opt_limiting_bounds)
      USE space_dim, ONLY: k_dim
      USE hyperbolic_matrices_module, ONLY : init_my_vectors,&
      x1vec, x2vec, x2_ghost, vec_loc
      USE st_matrix, ONLY : create_my_ghost
      USE my_util, ONLY : error_petsc
      IMPLICIT NONE
      CLASS(hyperbolic_type), INTENT(INOUT) :: this
      MPI_Comm, INTENT(IN) :: communicator
      CHARACTER(100), INTENT(IN) :: name
      TYPE(mesh_type), TARGET, INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
      INTEGER :: ierr, n
      REAL(KIND = 8), DIMENSION(2) :: times
      TYPE(limiting_bounds_type), OPTIONAL, TARGET :: opt_limiting_bounds

      this%syst_dim = k_dim + 2

      this%name = name
      this%mesh => mesh
      this%communicator = communicator
      this%LA => LA

      this%time = times(1) !<==initial_time
      this%final_time = times(2) !<==final_time

      CALL this%read_hyperbolic_data("HYPERBOLIC PARAMETERS FOR "//trim(adjustl(this%name)))

      !=== new Butcher module
      this%ERK%sv = this%erk_sv
      CALL this%ERK%init()
      !=== end new Butcher module

      this%matrices%method = this%method

      !===Periodic boundary if any
      CALL this%matrices%construct(this%communicator, this%mesh, this%LA)

      !===Goshting structures
      this%x1vec => x1vec
      this%x2vec => x2vec
      this%x2_ghost => x2_ghost
      this%vec_loc => vec_loc
      CALL VecDuplicate(this%x1vec, this%x3vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x4vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x5vec, ierr)

      ALLOCATE(this%tab(this%mesh%dom_np))
      DO n = 1, this%mesh%dom_np
         this%tab(n) = n - 1
      END DO

      !===Build boundary conditions
      CALL this%construct_bc(this%mesh, this%LA)

      !===Limiting
      CALL this%limiting%init(this%communicator, this%name, this%mesh, this%LA)
      IF (this%limiting%if_limiting) THEN
         IF (PRESENT(opt_limiting_bounds)) THEN
            this%limiting_bounds => opt_limiting_bounds
         ELSE
            CALL error_petsc("BUG in init_hyperbolic: if_limiting set to TRUE in data &
                             but you forgot to define add optional opt_limiting_bounds &
                             in init_hyperbolic.")
         END IF
      END IF
   END SUBROUTINE init_hyperbolic

   SUBROUTINE read_hyperbolic_data(this, section_name)
     USE read_inputs_module
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
     CALL read_data(argument_data%method, this%method, opt_name=this%name)

     !================
     !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
     !================

     CALL finalize_rewrite_data
   END SUBROUTINE read_hyperbolic_data


   SUBROUTINE update(this, un)
      USE space_dim
      USE my_util, ONLY : error_petsc
      USE cell_limiting_engine_module
      USE sub_plot
      USE compute_periodic, ONLY : periodic_rhs_petsc, periodic_vector_petsc
      USE st_matrix, ONLY: extract_through_ghost
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim)                :: un_temp
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim) :: ff
      REAL(KIND = 8), DIMENSION(this%mesh%np)                       :: rk
      REAL(KIND = 8), DIMENSION(this%mesh%np,2)                     :: bounds
      INTEGER :: comp, k, ierr, it

      un_temp = un
      SELECT CASE(this%method)
      CASE('viscous')
         !===compute dijL and dt
         CALL this%compute_dij(un_temp, bounds)
         CALL this%compute_dt
         this%time = this%time + this%dt
         DO comp = 1, this%syst_dim
            ff = 0.d0
            ff = this%flux(comp, un_temp)
            CALL VecZeroEntries(this%x3vec, ierr)
            DO k = 1, k_dim
               !=== set flux_k in x1vec
               CALL array_to_petsc_vec(ff(:, k), this%x1vec, this%mesh, this%LA, 'insert')
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

            !=== x3 <-- x2 / lumped_mass
            CALL VecPointWiseDivide(this%x3vec, this%x2vec, this%matrices%lump_mass_vec, ierr)
            !=== x3 <-- un + x3*dt   (x1 <-- un few lines above)
            CALL VecAYPX(this%x3vec, this%dt, this%x1vec, ierr)
            CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x3vec, this%LA)
            !=== un+1 <-- x3
            CALL extract_through_ghost(this%x3vec, this%x2_ghost, 1, 1, this%LA, un(:, comp), &
                                    'insert', opt_assemble=.FALSE.)


            CALL this%impose_bc(un, this%mesh, this%time)
         END DO
      CASE('high')
         !===compute dijL and dt
         CALL this%compute_dij(un_temp, bounds)
         CALL this%compute_dt

         this%time = this%time + this%dt

         DO comp = 1, this%syst_dim
            ff = this%flux(comp, un_temp)

            CALL VecSet(this%x3vec, 0.d0, ierr)
            DO k = 1, k_dim
               !=== set flux_k in x1vec
               CALL array_to_petsc_vec(ff(:, k), this%x1vec, this%mesh, this%LA, 'insert')
               !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
               CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
               !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
               CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
            END DO

            !=== set un(comp) in x1vec
            CALL array_to_petsc_vec(un_temp(:, comp), this%x1vec, this%mesh, this%LA, 'insert')
            !TEST Low order
!!$           IF (comp==1) THEN
!!$              CALL MatMultAdd(this%matrices%dijL, this%x1vec, this%x3vec, this%x4vec, ierr)
!!$              CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x4vec, this%LA)
!!$              CALL VecGhostGetLocalForm(this%x4vec, this%x2_ghost, ierr)
!!$              CALL VecGhostUpdateBegin(this%x4vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
!!$              CALL VecGhostUpdateEnd(this%x4vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
!!$              CALL extract(this%x2_ghost, 1, 1, this%LA, rk)
!!$              rk = rk * this%dt / this%matrices%lumped_mass
!!$              un(:, comp) = un_temp(:, comp) + rk
!!$              bounds(:,2) = un(:, 1)
!!$           END IF
           !END TEST
            !=== add dij un(comp)to x3vec in x2vec
            CALL MatMultAdd(this%matrices%dijH, this%x1vec, this%x3vec, this%x2vec, ierr)

            CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x2vec, this%LA)
            !=== Inverting mass matrix and updating un with dt
!======================== USING LUMPED MASS =========================!
            ! !=== x3 <-- x2 / lumped_mass
            ! CALL VecPointWiseDivide(this%x3vec, this%x2vec, this%matrices%lump_mass_vec, ierr)
            ! !=== x3 <-- un + x3*dt   (x1 <-- un few lines above)
            ! CALL VecAYPX(this%x3vec, this%dt, this%x1vec, ierr)
            ! CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x3vec, this%LA)
            ! !=== un+1 <-- x3
            ! CALL extract_through_ghost(this%x3vec, this%x2_ghost, 1, 1, this%LA, un(:, comp), &
            !                         'insert', opt_assemble=.FALSE.)
!======================== USING FULL MASS =========================!
            !=== x2 = rk
            !=== x3 <-- lump_inv @ rk
            CALL VecPointWiseDivide(this%x3vec, this%x2vec, this%matrices%lump_mass_vec, ierr)
            !=== x4 <-- Mass @ x3 (x4 <-- Mass@lump_inv@rk)
            CALL MatMult(this%matrices%mass, this%x3vec, this%x4vec, ierr)
            !=== x2 <-- lump_inv @ x4 (x2 <-- lump_inv @ Mass @ lump_inv @ rk)
            CALL VecPointWiseDivide(this%x2vec, this%x4vec, this%matrices%lump_mass_vec, ierr)
            !=== x3 <-- 2*x3 - x2 = (2I - lump_inv @ Mass) @ lump_inv @ rk
            !=== in petsc, VecAXPBY(y_vec, alpha, beta, x_vec) ==> y <-- y*beta + x*alpha ... ...
            CALL VecAXPBY(this%x3vec, -1.d0, 2.d0, this%x2vec, ierr)
            !=== x3 <-- dt*x3 + un (x1 <-- un a few lines above)
            CALL VecAYPX(this%x3vec, this%dt, this%x1vec, ierr)
            !=== Manually make un periodic and extract the result
            CALL periodic_vector_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%x3vec, this%LA)
            CALL extract_through_ghost(this%x3vec, this%x2_ghost, 1, 1, this%LA, un(:, comp), &
                                    'insert', opt_assemble=.FALSE.)
         END DO
         !===Limiting
         IF (this%limiting%if_limiting) THEN
            DO it = 1, this%limiting%limit_max
               CALL this%limiting%iterative_cell_limiting_procedure(un,bounds(:,1),&
                  this%limiting_bounds, 'MIN', un_temp)
               un(:,:) = un_temp(:,:)
            END DO
            DO it = 1, this%limiting%limit_max
               CALL this%limiting%iterative_cell_limiting_procedure(un,bounds(:,2),&
                  this%limiting_bounds, 'MAX', un_temp)
               un(:,:) = un_temp(:,:)
            END DO
         END IF
         !===Periodicity
         DO comp = 1, this%syst_dim 
            DO k = 1, this%mesh%per%nb_bords
               un(this%mesh%per%list(k)%DIL, comp) = un(this%mesh%per%perlist(k)%DIL, comp)
            END DO
         END DO

         !===Boundary conditions
         CALL this%impose_bc(un, this%mesh, this%time)
      CASE DEFAULT
         CALL error_petsc('Wrong method '//this%method//' in '//TRIM(ADJUSTL(this%name))//&
         ' update, should be "viscous" or "high"')
      END SELECT

   END SUBROUTINE update

!========================================================
!========== PRIVATE PROCEDURES ==========================
!========================================================


   SUBROUTINE compute_dt(this)
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      REAL(KIND = 8), DIMENSION(this%mesh%dom_np) :: dijL_diag
      REAL(KIND = 8) :: dt_min_loc, dt_min_glob
      INTEGER :: ierr

      CALL MatGetDiagonal(this%matrices%dijL, this%x1vec, ierr)
      CALL VecAbs(this%x1vec, ierr)
      CALL VecPointWiseDivide(this%x2vec, this%matrices%lump_mass_vec, this%x1vec, ierr)
      CALL VecMin(this%x2vec, PETSC_NULL_INTEGER, dt_min_glob, ierr)

      this%dt = this%CFL * dt_min_glob / 2
   END SUBROUTINE compute_dt

   SUBROUTINE compute_dij(this, un, bounds)
      USE space_dim
      USE petsc
      ! USE my_util, ONLY : error_petsc
      USE def_type_mesh
      USE arbitrary_eos_lambda_module
      USE compute_periodic
      USE st_matrix, ONLY: extract_through_ghost
      IMPLICIT NONE
      CLASS(hyperbolic_type) :: this
      TYPE(mesh_type), POINTER :: mesh
      TYPE(petsc_csr_LA), POINTER :: LA
      REAL(KIND = 8), DIMENSION(:, :) :: un, bounds
      REAL(KIND = 8), DIMENSION(this%mesh%np) :: arr
      INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge, n_size
      INTEGER, DIMENSION(1) :: i_t, j_t, idx, jdx
      REAL(KIND = 8), DIMENSION(1, k_dim) :: nij_c
      REAL(KIND = 8), DIMENSION(1) :: norm_c, dijL_c
      REAL(KIND = 8), DIMENSION(1) :: dijH_c
      REAL(KIND = 8), DIMENSION(2) :: u, rho, ie, p, lambda_max
      REAL(KIND = 8) :: pstar, max_lambda, uijbar
      LOGICAL, DIMENSION(this%mesh%medge) :: virgin_edge
      REAL(KIND = 8), DIMENSION(this%mesh%np)  :: alpha !<==commutator in (0,1)
      real(kind=8) :: norm
      !===Compute commutator if needed
      IF (this%method=='high') THEN
         CALL this%commutator(un, alpha)
      END IF

      !===Compute dijL
      CALL MatZeroEntries(this%matrices%dijL, ierr)
      IF (this%method=='high') THEN
         CALL MatZeroEntries(this%matrices%dijH, ierr)
      END IF
      mesh => this%mesh
      LA => this%LA

      virgin_edge = .TRUE.
      nw = mesh%gauss%n_w
      bounds(:,1) = un(:,1)
      bounds(:,2) = un(:,1)
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

               CALL this%compute_lambda(un, i, j, lambda_max, on_edge = .FALSE.)
               CALL MatGetValues(this%matrices%cij_norm_loc, 1, i_t - 1, 1, j_t - 1, norm_c, ierr)

               max_lambda = MAXVAL(lambda_max)
               dijL_c = max_lambda * norm_c

               IF (mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j
                     CALL this%compute_lambda(un, i, j, lambda_max, on_edge = .TRUE.)

                     dijL_c = MAX(dijL_c, MAXVAL(lambda_max) * norm_c)
                     max_lambda = MAX(max_lambda,MAXVAL(lambda_max))
               END IF

               idx = LA%loc_to_glob(1, i) - 1
               jdx = LA%loc_to_glob(1, j) - 1

               CALL MatSetValues(this%matrices%dijL, 1, idx, 1, jdx, dijL_c, ADD_VALUES, ierr)
               CALL MatSetValues(this%matrices%dijL, 1, jdx, 1, idx, dijL_c, ADD_VALUES, ierr)

               CALL MatSetValues(this%matrices%dijL, 1, idx, 1, idx, -dijL_c, ADD_VALUES, ierr) !===add value on diagonal
               CALL MatSetValues(this%matrices%dijL, 1, jdx, 1, jdx, -dijL_c, ADD_VALUES, ierr) !===add value on diagonal
               IF (this%method=='high') THEN
                  dijH_c = dijL_c*(alpha(i)+alpha(j))/2
                  CALL MatSetValues(this%matrices%dijH, 1, idx, 1, jdx, dijH_c, ADD_VALUES, ierr)
                  CALL MatSetValues(this%matrices%dijH, 1, jdx, 1, idx, dijH_c, ADD_VALUES, ierr)
                  CALL MatSetValues(this%matrices%dijH, 1, idx, 1, idx, -dijH_c, ADD_VALUES, ierr) !===add value on diagonal
                  CALL MatSetValues(this%matrices%dijH, 1, jdx, 1, jdx, -dijH_c, ADD_VALUES, ierr) !===add value on diagonal
                  !===Compute low-order update to estimate bounds
                  DO k = 1, k_dim
                     CALL MatGetValues(this%matrices%nij_loc(k), 1, i_t - 1, 1, j_t - 1, nij_c(:, k), ierr)
                  END DO
                  uijbar = (un(i, 1)+un(j, 1))/2 &
                        - SUM((un(j, 2:k_dim+1) - un(i, 2:k_dim+1))*nij_c(1, :))/(2*max_lambda)
                  bounds(i,1) = MIN(bounds(i,1),uijbar)
                  bounds(i,2) = MAX(bounds(i,2),uijbar)
                  bounds(j,1) = MIN(bounds(j,1),uijbar)
                  bounds(j,2) = MAX(bounds(j,2),uijbar)
                  !===End compute low-order update to estimate bounds
               END IF
            END IF
         END DO
      END DO
      CALL MatAssemblyBegin(this%matrices%dijL, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd  (this%matrices%dijL, MAT_FINAL_ASSEMBLY, ierr)

      IF (this%method=='high') THEN
         CALL MatAssemblyBegin(this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)
         CALL MatAssemblyEnd  (this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)

         CALL array_to_petsc_vec(bounds(:,1), this%x1vec, this%mesh, this%LA, 'insert')
         CALL extract_through_ghost(this%x1vec, this%x2_ghost, 1, 1, this%LA, bounds(:, 1), &
                                 'min', opt_assemble=.FALSE.)

         CALL array_to_petsc_vec(bounds(:,2), this%x1vec, this%mesh, this%LA, 'insert')
         CALL extract_through_ghost(this%x1vec, this%x2_ghost, 1, 1, this%LA, bounds(:, 2), &
                                 'max', opt_assemble=.FALSE.)
      END IF
   END SUBROUTINE compute_dij


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
      eta = un(:,1)
      logeta = log(abs(eta))
      norm_diff = 0.d0
      norm_log = 0.d0

      DO k = 1, k_dim
         CALL array_to_petsc_vec(logeta, this%x1vec, this%mesh, this%LA, 'insert') !<==v1 = log(eta)
         CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr) !<==v2 = dk(log(eta))
        
         CALL array_to_petsc_vec(eta,    this%x1vec, this%mesh, this%LA, 'insert') !<==v1 = eta
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
     
      CALL extract_through_ghost(this%x4vec, this%x2_ghost, 1, 1, this%LA, rk, &
                              'insert', opt_assemble=.FALSE.)

      CALL extract_through_ghost(this%x5vec, this%x2_ghost, 1, 1, this%LA, rk_norm, &
                              'insert', opt_assemble=.FALSE.)
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

END MODULE abstract_hyperbolic_module