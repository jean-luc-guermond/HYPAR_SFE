MODULE abstract_hyperbolic_module

    USE hyperbolic_matrices_module, ONLY: hyperbolic_matrices_type, init_my_vectors

    TYPE, ABSTRACT :: hyperbolic_type
        !===Parameters read from data
        REAL(KIND=8)                 :: CFL       = 0.5d0
        CHARACTER(LEN=rec_length)    :: method = 'viscous'
        REAL(KIND = 8), DIMENSION(1) :: eos_param = 0.d0
        INTEGER                      :: erk_sv    = -21
        !===Parameters built along way
        MPI_Comm :: communicator
        Vec, POINTER :: x1vec, x2vec, x2_ghost, vec_loc
        Vec          :: x3vec, x4vec, x5vec
        CHARACTER(100) :: name
        LOGICAL :: no_iter
        INTEGER :: syst_dim
        REAL(KIND = 8) :: dt, time, final_time, in_tol
        INTEGER, DIMENSION(:), POINTER :: tab
        TYPE(mesh_type),     POINTER :: mesh
        TYPE(petsc_csr_LA),  POINTER :: LA
        TYPE(BT),             PUBLIC :: ERK
        TYPE(hyperbolic_bc_type)     :: bc
        ! TYPE(euler_bc_type)          :: euler_bc
        TYPE(hyperbolic_bc_type)       :: bc
        TYPE(hyperbolic_matrices_type) :: matrices
        TYPE(limiting_type)            :: limiting
        CLASS(limiting_bounds_type), ALLOCATABLE :: limiting_bounds(:)
        ! PROCEDURE(function_template_pressure),  NOPASS, POINTER :: pressure
        ! PROCEDURE(function_template_impose_bc), NOPASS, POINTER :: impose_bc
   CONTAINS
        PROCEDURE, PUBLIC  :: init => init_hyperbolic
        ! PROCEDURE, PUBLIC  :: read_euler_data
        PROCEDURE, PUBLIC   :: update
        PROCEDURE, PRIVATE  :: compute_dij, compute_dt
        PROCEDURE, DEFERRED :: flux
        PROCEDURE, DEFERRED :: impose_bc
        ! PROCEDURE, PRIVATE :: compute_dK, compute_dt_from_dK
    END TYPE hyperbolic_type

CONTAINS

   SUBROUTINE init_hyperbolic(this, communicator, name, mesh, LA, times)! pressure, impose_bc, times)
      USE space_dim, ONLY: k_dim
      USE euler_matrices_module, ONLY : init_my_vectors,&
      x1vec, x2vec, x2_ghost, vec_loc
      USE st_matrix, ONLY : create_my_ghost
      IMPLICIT NONE
      CLASS(hyperbolic_type), INTENT(INOUT) :: this
      MPI_Comm, INTENT(IN) :: communicator
      CHARACTER(100) :: name
      TYPE(mesh_type), TARGET, INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
      INTEGER :: ierr, n
      REAL(KIND = 8), DIMENSION(2) :: times
      PROCEDURE(function_template_pressure) :: pressure
      PROCEDURE(function_template_impose_bc) :: impose_bc

      this%syst_dim = k_dim + 2

      this%name = name
      this%mesh => mesh
      this%communicator = communicator
      this%LA => LA
    !   this%pressure => pressure
    !   this%impose_bc => impose_bc
      this%bc%syst_dim = this%syst_dim
      ALLOCATE(this%limiting_bounds(this%syst_dim))
    !   this%euler_bc%syst_dim = this%syst_dim
      this%time = times(1) !<==initial_time
      this%final_time = times(2) !<==final_time

      !===Parameters for lambda_arbitrary_eos
      this%in_tol = 1.d-2
      this%no_iter = .TRUE.

      CALL this%read_euler_data("EULER PARAMETERS")

      !=== new Butcher module
      this%ERK%sv = this%erk_sv
      CALL this%ERK%init()
      !=== end new Butcher module

      this%matrices%method = this%method !<==transfer this%method to this%matrices

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

      !===Boundary conditions
      CALL this%bc%construct_bc(this%mesh, this%LA)
    !   CALL this%euler_bc%construct_euler_bc(this%mesh, this%LA)

      !===Limiting
      CALL this%limiting%init(this%communicator, this%name, this%mesh, this%LA)
   END SUBROUTINE init_hyperbolic


   SUBROUTINE update(this, un)
      USE space_dim
      USE my_util, ONLY : error_petsc
      USE cell_limiting_engine_module
      USE sub_plot
      USE compute_periodic, ONLY : periodic_rhs_petsc, periodic_vector_petsc
      CLASS(hyperbolic_type)                                                :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim)                :: un_temp
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim) :: ff
      REAL(KIND = 8), DIMENSION(this%mesh%np)                       :: rk
      REAL(KIND = 8), DIMENSION(this%mesh%np,2)                     :: bounds
      INTEGER :: comp, k, ierr, it

      INTEGER, PARAMETER :: limit_max = 2 !<<FIXME
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

            CALL this%impose_bc(un, this%bc, this%mesh, this%time)
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
            DO it = 1, limit_max
                  CALL this%limiting%iterative_cell_limiting_procedure(un,bounds(:,1),&
                     psi_rho_min,zero_of_psi_rho_min,un_temp)
                  un(:,:) = un_temp(:,:)
            END DO
            DO it = 1, limit_max
                  CALL this%limiting%iterative_cell_limiting_procedure(un,bounds(:,2),&
                     psi_rho_max,zero_of_psi_rho_max,un_temp)
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
         CALL this%impose_bc(un, this%bc, this%mesh, this%time)
      CASE DEFAULT
         CALL error_petsc('Wrong method '//this%method//' in '//TRIM(ADJUSTL(this%name))//&
         ' update, should be "viscous" or "high"')
      END SELECT

   CONTAINS
     FUNCTION psi_rho_min(x,psi_m) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:) :: x
       REAL(KIND=8) :: psi_m, v
       v = x(1)-psi_m
     END FUNCTION psi_rho_min

     FUNCTION zero_of_psi_rho_min(psi_m,u0,P) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:) :: u0, P
       REAL(KIND=8) :: psi_m, v
       v = (psi_m-u0(1))/P(1)
     END FUNCTION zero_of_psi_rho_min

     FUNCTION psi_rho_max(x,psi_m) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:) :: x
       REAL(KIND=8) :: psi_m, v
       v = psi_m-x(1)
     END FUNCTION psi_rho_max

     FUNCTION zero_of_psi_rho_max(psi_m,u0,P) RESULT(v)
       IMPLICIT NONE
       REAL(KIND=8), DIMENSION(:) :: u0, P
       REAL(KIND=8) :: psi_m, v
       v = (psi_m-u0(1))/P(1)
     END FUNCTION zero_of_psi_rho_max

   END SUBROUTINE update


END MODULE abstract_hyperbolic_module