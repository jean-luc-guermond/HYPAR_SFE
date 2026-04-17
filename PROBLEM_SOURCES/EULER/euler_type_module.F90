MODULE euler_type_MODULE
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE euler_bc_arrays
   USE Butcher_tableau
   USE euler_matrices_module
   USE mesh_parameters
   USE periodic_data_module
   IMPLICIT NONE

   INTEGER, PARAMETER, PRIVATE :: rec_length = 200

   ABSTRACT INTERFACE
      FUNCTION function_template_pressure(rho, ie) RESULT(vv)
         REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, ie
         REAL(KIND = 8), DIMENSION(SIZE(rho, 1)) :: vv
      END FUNCTION function_template_pressure
   END INTERFACE

   ABSTRACT INTERFACE
      SUBROUTINE function_template_impose_bc(un, euler_bc, mesh, time)
         USE def_type_mesh
         USE euler_bc_arrays
         REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
         TYPE(mesh_type)     :: mesh
         TYPE(euler_bc_type) :: euler_bc
         REAL(KIND = 8)      :: time
      END SUBROUTINE function_template_impose_bc
   END INTERFACE

   TYPE argument_euler_type
      CHARACTER(LEN=rec_length) :: CFL       = '=== CFL ? ==='
      CHARACTER(LEN=rec_length) :: method = '=== Which method to solve Euler (viscous, high) ? ==='
      CHARACTER(LEN=rec_length) :: erk_sv
      CHARACTER(LEN=rec_length) :: eos_param 
   END TYPE argument_euler_type
   
   TYPE euler_type
      !===Parameters read from data
      REAL(KIND=8)                 :: CFL       = 0.5d0
      CHARACTER(LEN=rec_length)    :: method = 'viscous'
      REAL(KIND = 8), DIMENSION(1) :: eos_param = 0.d0
      INTEGER                      :: erk_sv    = -21
      !===Parameters built along way
      MPI_Comm :: communicator
      CHARACTER(100) :: name
      TYPE(mesh_type),     POINTER :: mesh
      TYPE(petsc_csr_LA),  POINTER :: LA
      TYPE(periodic_type), POINTER :: per
      PROCEDURE(function_template_pressure),  NOPASS, POINTER :: pressure
      PROCEDURE(function_template_impose_bc), NOPASS, POINTER :: impose_bc
      TYPE(BT), PUBLIC :: ERK
      TYPE(euler_bc_type) :: euler_bc
      TYPE(euler_matrices_type) :: matrices
      REAL(KIND = 8) :: dt, time, in_tol
      LOGICAL :: no_iter
      INTEGER :: syst_dim
      Vec, PRIVATE :: x1vec, x2vec, x3vec, x2_ghost, vec_loc
      INTEGER, DIMENSION(:), POINTER :: tab
   CONTAINS
      PROCEDURE, PUBLIC  :: init => init_euler
      PROCEDURE, PUBLIC  :: read_euler_data
      PROCEDURE, PUBLIC  :: update
      PROCEDURE, PRIVATE :: compute_dij, compute_dt
      PROCEDURE, PRIVATE :: flux
      PROCEDURE, PRIVATE :: compute_dK, compute_dt_from_dK
   END TYPE euler_type

CONTAINS
   SUBROUTINE init_euler(this, communicator, name, mesh, LA, per, pressure, impose_bc, time_init)
      USE st_matrix
      USE periodic_data_module
      CLASS(euler_type), INTENT(INOUT) :: this
      MPI_Comm, INTENT(IN) :: communicator
      CHARACTER(100) :: name
      TYPE(mesh_type), TARGET, INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
      TYPE(periodic_type), TARGET, INTENT(IN) :: per
      INTEGER :: ierr, n
      REAL(KIND = 8) :: time_init
      PROCEDURE(function_template_pressure) :: pressure
      PROCEDURE(function_template_impose_bc) :: impose_bc
      INTEGER, POINTER, DIMENSION(:) :: ifrom

      this%syst_dim = mesh_data_info%k_dim + 2

      this%name = name
      this%mesh => mesh
      this%communicator = communicator
      this%LA => LA
      this%per => per
      this%pressure => pressure
      this%impose_bc => impose_bc
      this%euler_bc%syst_dim = this%syst_dim
      this%time = time_init

      !===Parameters for lambda_arbitrary_eos
      this%in_tol = 1.d-2
      this%no_iter = .TRUE.
      !this%eos_param(1) = 0.d0 !===b_covolume
      !===CFL number
      !this%CFL = 0.5d0

      CALL this%read_euler_data("EULER PARAMETERS")
      !=== new Butcher module
      this%ERK%sv = this%erk_sv
      CALL this%ERK%init()
      !CALL this%ERK%init(this%erk_sv)
      !=== new Butcher module
      this%matrices%method = this%method !<==transfer this%method to this%matrices
      CALL this%euler_bc%construct_euler_bc(this%mesh)
      CALL this%matrices%construct(this%communicator, this%mesh, this%LA, this%per)
    
      CALL create_my_ghost(this%mesh, this%LA, ifrom)
      CALL VecCreateGhost(this%communicator, this%mesh%dom_np, &
           PETSC_DETERMINE, SIZE(ifrom), ifrom, this%x1vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x2vec, ierr)
      CALL VecDuplicate(this%x1vec, this%x3vec, ierr)
      CALL VecGhostGetLocalForm(this%x2vec, this%x2_ghost, ierr)

      CALL VecCreateSeq(PETSC_COMM_SELF, this%mesh%dom_np, this%vec_loc, ierr)
      ALLOCATE(this%tab(this%mesh%dom_np))
      DO n = 1, this%mesh%dom_np
         this%tab(n) = n - 1
      END DO
      
   END SUBROUTINE init_euler

   SUBROUTINE read_euler_data(this, section_name)
     USE character_strings
     IMPLICIT NONE
     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name

     CLASS(euler_type), INTENT(INOUT) :: this
     TYPE(argument_euler_type)        :: argument_data


     !===Initialize data arguments (depends on the name)

     argument_data%eos_param = '===' // TRIM(ADJUSTL(this%name)) // ': b_covolume ?==='
     argument_data%erk_sv = '===' // TRIM(ADJUSTL(this%name)) // ': ERK ?==='

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
     CALL read_data(argument_data%CFL, this%CFL)

     !===b_covolume
     CALL read_data(argument_data%eos_param, this%eos_param(1))

     !===ERK
     CALL read_data(argument_data%erk_sv, this%erk_sv)

     !===Method order
     CALL read_data(argument_data%method, this%method)

     !================
     !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
     !================

     CALL finalize_rewrite_data

   END SUBROUTINE read_euler_data

   SUBROUTINE update(this, un)
     USE petsc_tools
     USE st_matrix
     USE my_util
     CLASS(euler_type) :: this
     REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
     REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim) :: un_temp
     REAL(KIND = 8), DIMENSION(this%mesh%np, mesh_data_info%k_dim) :: ff
     REAL(KIND = 8), DIMENSION(this%mesh%np) :: rk
     INTEGER :: k, comp, ierr
     un_temp = un

     SELECT CASE(this%method)
     CASE('viscous')
        !===compute dijL and dt
        CALL this%compute_dij(un_temp)
        CALL this%compute_dt(un_temp)

        this%time = this%time + this%dt

        DO comp = 1, this%syst_dim
           ff = 0.d0
           ff = this%flux(comp, un_temp)

           CALL VecSet(this%x3vec, 0.d0, ierr)
           DO k = 1, mesh_data_info%k_dim
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
           CALL periodic_rhs_petsc(this%per%nb_bords, this%per%list, this%per%perlist, this%x2vec, this%LA)
           CALL VecGhostGetLocalForm(this%x2vec, this%x2_ghost, ierr)
           CALL VecGhostUpdateBegin(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
           CALL VecGhostUpdateEnd(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
           CALL extract(this%x2_ghost, 1, 1, this%LA, rk)

           rk = rk * this%dt / this%matrices%lumped_mass

           un(:, comp) = un_temp(:, comp) + rk

           DO k = 1, this%per%nb_bords
              un(this%per%list(k)%DIL, comp) = un(this%per%perlist(k)%DIL, comp)
           END DO

           CALL this%impose_bc(un, this%euler_bc, this%mesh, this%time)
        END DO
     CASE('high')
         !===compute dijL and dt
        CALL this%compute_dij(un_temp)
        CALL this%compute_dt(un_temp)

        this%time = this%time + this%dt

        DO comp = 1, this%syst_dim
           ff = 0.d0
           ff = this%flux(comp, un_temp)

!!$           CALL VecSet(this%x3vec, 0.d0, ierr)
!!$           DO k = 1, mesh_data_info%k_dim
!!$              !=== set flux_k in x1vec
!!$              CALL array_to_petsc_vec(ff(:, k), this%x1vec, this%mesh, this%LA, 'insert')
!!$              !=== compute sum_j (cij_k * fluxj_k) and store into x2vec
!!$              CALL MatMult(this%matrices%cij(k), this%x1vec, this%x2vec, ierr)
!!$              !=== compute sum_k (sum_j (cij_k * flux_k)) and store into x3vec
!!$              CALL VecAXPY(this%x3vec, -1.d0, this%x2vec, ierr)
!!$           END DO

           CALL compute_flux(this, ff, this%x3vec)

           !=== set un(comp) in x1vec
           CALL array_to_petsc_vec(un_temp(:, comp), this%x1vec, this%mesh, this%LA, 'insert')
           !=== add dij un(comp)to x3vec in x2vec
           CALL MatMultAdd(this%matrices%dijL, this%x1vec, this%x3vec, this%x2vec, ierr)
           CALL periodic_rhs_petsc(this%per%nb_bords, this%per%list, this%per%perlist, this%x2vec, this%LA)
           CALL VecGhostGetLocalForm(this%x2vec, this%x2_ghost, ierr)
           CALL VecGhostUpdateBegin(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
           CALL VecGhostUpdateEnd(this%x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
           CALL extract(this%x2_ghost, 1, 1, this%LA, rk)

           rk = rk * this%dt / this%matrices%lumped_mass

           un(:, comp) = un_temp(:, comp) + rk

           DO k = 1, this%per%nb_bords
              un(this%per%list(k)%DIL, comp) = un(this%per%perlist(k)%DIL, comp)
           END DO

           CALL this%impose_bc(un, this%euler_bc, this%mesh, this%time)
        END DO
     CASE DEFAULT
        CALL error_petsc('Wrong method in euler update')
     END SELECT

   END SUBROUTINE update

   FUNCTION flux(this, comp, un) RESULT(vv)  
      IMPLICIT NONE
      CLASS(euler_type)                           :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
      INTEGER,                         INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(SIZE(un, 1), mesh_data_info%k_dim) :: vv
      REAL(KIND = 8), DIMENSION(SIZE(un, 1))                       :: H, u, ie
      INTEGER :: k

      ! SELECT CASE(comp)
      ! CASE(1)
      IF (comp == 1) THEN
         DO k = 1, mesh_data_info%k_dim
            vv(:, k) = un(:, k + 1)
         END DO
      ! CASE(2:k_dim + 1)
      ELSE IF ((comp>=2) .AND. (comp<=mesh_data_info%k_dim+1)) THEN
         u = un(:, comp) / un(:, 1)
         DO k = 1, mesh_data_info%k_dim
            vv(:, k) = un(:, k + 1) * u
         END DO
         ie = un(:, mesh_data_info%k_dim + 2) / un(:, 1)
         DO k = 1, mesh_data_info%k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO
         vv(:, comp - 1) = vv(:, comp - 1) + this%pressure(un(:, 1), ie)
      ! CASE(mesh_data_info%k_dim + 2)
      ELSE IF (comp == mesh_data_info%k_dim + 2) THEN
         ie = un(:, mesh_data_info%k_dim + 2) / un(:, 1)
         DO k = 1, mesh_data_info%k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO

         H = un(:, comp) + this%pressure(un(:, 1), ie)
         DO k = 1, mesh_data_info%k_dim
            vv(:, k) = (un(:, k + 1) / un(:, 1)) * H
         END DO
      ! CASE DEFAULT
      ELSE
         WRITE(*, *) ' BUG in flux, wrong comp = ', comp
         WRITE(*,*) 'pb dimension is ', mesh_data_info%k_dim
         STOP
      ! END SELECT
      END IF
   END FUNCTION flux
   

!========================================================
!========== PRIVATE PROCEDURES ==========================
!========================================================


   SUBROUTINE compute_dt(this, un)
      IMPLICIT NONE
      CLASS(euler_type) :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
      REAL(KIND = 8), DIMENSION(this%mesh%dom_np) :: dijL_diag
      REAL(KIND = 8) :: dt_min_loc, dt_min_glob
      INTEGER :: ierr
   
      CALL MatGetDiagonal(this%matrices%dijL, this%vec_loc, ierr)
      CALL VecGetValues(this%vec_loc, this%mesh%dom_np, this%tab, dijL_diag, ierr)
      dijL_diag = this%matrices%lumped_mass(1:this%mesh%dom_np) / ABS(dijL_diag)

      dt_min_loc = MINVAL(dijL_diag) / 2.d0

      CALL MPI_ALLREDUCE(dt_min_loc, dt_min_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PETSC_COMM_WORLD, ierr)
      this%dt = this%CFL * dt_min_glob
   END SUBROUTINE compute_dt


   SUBROUTINE compute_dij(this, un)
     USE mesh_parameters
     ! USE space_dim
     USE petsc
     USE my_util
     USE def_type_mesh
     USE arbitrary_eos_lambda_module
     USE compute_periodic
     IMPLICIT NONE
     CLASS(euler_type) :: this
     TYPE(mesh_type), POINTER :: mesh
     TYPE(petsc_csr_LA), POINTER :: LA
     REAL(KIND = 8), DIMENSION(:, :) :: un
     INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge
     INTEGER, DIMENSION(1) :: i_t, j_t, idx, jdx
     REAL(KIND = 8), DIMENSION(1, mesh_data_info%k_dim) :: nij_c
     REAL(KIND = 8), DIMENSION(1) :: norm_c, dijL_c
!!$     REAL(KIND = 8), DIMENSION(1) :: dijH_c
     REAL(KIND = 8), DIMENSION(2) :: u, rho, ie, p, lambda_max
     REAL(KIND = 8) :: pstar
     LOGICAL, DIMENSION(this%mesh%medge) :: virgin_edge
!!$     REAL(KIND = 8), DIMENSION(this%mesh%np)  :: e, alpha !<==commutator in (0,1)
     
!!$     !===Compute commutator if needed
!!$     IF (this%method=='high') THEN
!!$        e = 0.d0
!!$        DO k = 1, mesh_data_info%k_dim
!!$           e = e + un(:,k+1)**2
!!$        END DO
!!$        e = un(:,mesh_data_info%k_dim+2)/un(:,1) - 0.5d0*e/un(:,1)**2
!!$        write(*,*) 'before commutator'
!!$        CALL commutator(this%pressure(un(:,1),e), this%matrices%stiffL, this%LA, alpha)
!!$     END IF
!!$     
     !===Compute dijL
     CALL MatZeroEntries(this%matrices%dijL, ierr)

     mesh => this%mesh
     LA => this%LA

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

              DO k = 1, mesh_data_info%k_dim
                 CALL MatGetValues(this%matrices%nij_loc(k), 1, i_t - 1, 1, j_t - 1, nij_c(:, k), ierr)
              END DO

              rho(1) = un(i, 1)
              rho(2) = un(j, 1)

              u(1) = SUM(un(i, 2:1 + mesh_data_info%k_dim) * nij_c(1, :)) / rho(1)
              u(2) = SUM(un(j, 2:1 + mesh_data_info%k_dim) * nij_c(1, :)) / rho(2)

              ie(1) = un(i, mesh_data_info%k_dim + 2) / rho(1) - 0.5d0 * u(1) * u(1)
              ie(2) = un(j, mesh_data_info%k_dim + 2) / rho(2) - 0.5d0 * u(2) * u(2)

              p = this%pressure(rho, ie)

              CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, lambda_max, pstar)
              CALL MatGetValues(this%matrices%cij_norm_loc, 1, i_t - 1, 1, j_t - 1, norm_c, ierr)

              dijL_c = MAXVAL(lambda_max) * norm_c

              IF (mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j

                 DO k = 1, mesh_data_info%k_dim
                    CALL MatGetValues(this%matrices%nij_loc(k), 1, j_t - 1, 1, i_t - 1, nij_c(:, k), ierr)
                 END DO

                 u(1) = SUM(un(i, 2:1 + mesh_data_info%k_dim) * nij_c(1, :)) / rho(1)
                 u(2) = SUM(un(j, 2:1 + mesh_data_info%k_dim) * nij_c(1, :)) / rho(2)

                 rho = (/rho(2), rho(1)/)
                 ie = (/ie(2), ie(1)/)
                 p = (/p(2), p(1)/)

                 CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, lambda_max, pstar)

                 dijL_c = MAX(dijL_c, MAXVAL(lambda_max) * norm_c)

              END IF
              
              idx = LA%loc_to_glob(1, i) - 1
              jdx = LA%loc_to_glob(1, j) - 1

              CALL MatSetValues(this%matrices%dijL, 1, idx, 1, jdx, dijL_c, ADD_VALUES, ierr)
              CALL MatSetValues(this%matrices%dijL, 1, jdx, 1, idx, dijL_c, ADD_VALUES, ierr)
              CALL MatSetValues(this%matrices%dijL, 1, idx, 1, idx, -dijL_c, ADD_VALUES, ierr) !===add value on diagonal
              CALL MatSetValues(this%matrices%dijL, 1, jdx, 1, jdx, -dijL_c, ADD_VALUES, ierr) !===add value on diagonal
!!$              IF (this%method=='high') THEN
!!$                 dijH_c = dijL_c*(alpha(i)+alpha(j))/2
!!$                 CALL MatSetValues(this%matrices%dijH, 1, idx, 1, jdx, dijH_c, ADD_VALUES, ierr)
!!$                 CALL MatSetValues(this%matrices%dijH, 1, jdx, 1, idx, dijH_c, ADD_VALUES, ierr)
!!$              END IF
           END IF

        END DO

     END DO

     CALL MatAssemblyBegin(this%matrices%dijL, MAT_FINAL_ASSEMBLY, ierr)
     CALL MatAssemblyEnd  (this%matrices%dijL, MAT_FINAL_ASSEMBLY, ierr)
     
!!$     IF (this%method=='high') THEN
!!$        CALL MatAssemblyBegin(this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)
!!$        CALL MatAssemblyEnd  (this%matrices%dijH, MAT_FINAL_ASSEMBLY, ierr)
!!$     END IF

   END SUBROUTINE compute_dij

   SUBROUTINE compute_dk (this, un)
     USE arbitrary_eos_lambda_module
     IMPLICIT NONE
     CLASS(euler_type) :: this
     REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
     INTEGER, DIMENSION(1) :: i_t, j_t
     REAL(KIND = 8), DIMENSION(1, this%mesh%gauss%k_d) :: nij_c
     REAL(KIND = 8), DIMENSION(1) :: norm_c, dijL_c
     REAL(KIND = 8), DIMENSION(2) :: u, rho, ie, p, lambda_max
     LOGICAL, DIMENSION(this%mesh%medge) :: virgin_edge
     REAL(KIND = 8) :: pstar
     LOGICAL :: bug
     INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge, divider, nb_shared_cell
     nw = this%mesh%gauss%n_w

     bug = .FALSE.
     SELECT CASE(this%mesh%gauss%k_d)
     CASE(1)
        nb_shared_cell = 1
        IF (this%mesh%gauss%n_w/=2) bug=.TRUE.
     CASE(2)
        nb_shared_cell = 2
        IF (this%mesh%gauss%n_w/=3) bug=.TRUE.
     END SELECT
     IF (bug) THEN
        CALL error_petsc('Wrong polynomial degree for low-order viscosity')
     END IF

     DO m = 1, this%mesh%dom_me
        DO n = 1, this%mesh%gauss%n_e
           IF (this%mesh%attr_e(this%mesh%jce(n, m))) THEN
              edge = this%mesh%jce_loc(n, m)
              IF (.NOT. virgin_edge(edge)) CYCLE
              virgin_edge(edge) = .FALSE.
              ni = MOD(n, nw) + 1
              nj = MOD(n + 1, nw) + 1
              i = this%mesh%jj(ni, m)
              j = this%mesh%jj(nj, m)
              i_t = i
              j_t = j
              DO k = 1, this%mesh%gauss%k_d
                 CALL MatGetValues(this%matrices%nij_loc(k), 1, i_t - 1, 1, j_t - 1, &
                      nij_c(:, k), ierr)
              END DO
              rho(1) = un(i, 1)
              rho(2) = un(j, 1)
              u(1) = SUM(un(i, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(1)
              u(2) = SUM(un(j, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(2)
              ie(1) = un(i, this%mesh%gauss%k_d + 2) / rho(1) - 0.5d0 * u(1) * u(1)
              ie(2) = un(j, this%mesh%gauss%k_d + 2) / rho(2) - 0.5d0 * u(2) * u(2)
              p = this%pressure(rho, ie)
              CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, &
                   lambda_max, pstar)
              dijL_c = MAXVAL(lambda_max) * norm_c
              divider = nb_shared_cell
              
              IF (this%mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j
                 DO k = 1, this%mesh%gauss%k_d
                    CALL MatGetValues(this%matrices%nij_loc(k), 1, j_t - 1, 1, i_t - 1, &
                         nij_c(:, k), ierr)
                 END DO
                 u(1) = SUM(un(i, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(1)
                 u(2) = SUM(un(j, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(2)
                 rho = (/rho(2), rho(1)/)
                 ie = (/ie(2), ie(1)/)
                 p = (/p(2), p(1)/)
                 CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, &
                      lambda_max, pstar)
                 dijL_c = MAX(dijL_c, MAXVAL(lambda_max) * norm_c)
                 divider = 1
              END IF
              
              this%matrices%dK(m) = MAX(this%matrices%dK(m),dijL_c(1)/divider)
           END IF
        END DO
     END DO
   END SUBROUTINE compute_dk

   SUBROUTINE compute_dt_from_dK(this)
     IMPLICIT NONE
     CLASS(euler_type) :: this
     REAL(KIND = 8), DIMENSION(this%mesh%dom_np) :: dijL_diag
     REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w) :: v_loc
     INTEGER, DIMENSION(this%mesh%gauss%n_w) :: idxm
     INTEGER :: i, m, ni, iglob
     REAL(KIND = 8) :: dt_min_loc, dt_min_glob
     Vec                                         :: vect
     PetscErrorCode                              :: ierr
     CALL VecSet(vect, 0.d0, ierr)

     DO m = 1, this%mesh%me
        v_loc = 0.d0
        DO ni = 1, this%mesh%gauss%n_w
           i = this%mesh%jj(ni, m)
           iglob = this%LA%loc_to_glob(1, i)
           idxm(ni) = iglob - 1
           v_loc(ni) = v_loc(ni) + this%matrices%dK(m)
        ENDDO
        CALL VecSetValues(vect, this%mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
     ENDDO
     CALL VecAssemblyBegin(vect, ierr)
     CALL VecAssemblyEnd(vect, ierr)

     CALL VecGetValues(this%vec_loc, this%mesh%dom_np, this%tab, dijL_diag, ierr)
     dijL_diag = this%matrices%lumped_mass(1:this%mesh%dom_np) / ABS(dijL_diag)

     dt_min_loc = MINVAL(dijL_diag) / 2.d0

     CALL MPI_ALLREDUCE(dt_min_loc, dt_min_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PETSC_COMM_WORLD, ierr)
     this%dt = this%CFL * dt_min_glob
   END SUBROUTINE compute_dt_from_dK

   SUBROUTINE compute_flux(this, ff, Vect)
     IMPLICIT NONE
     CLASS(euler_type) :: this
     REAL(KIND = 8), DIMENSION(this%mesh%np, mesh_data_info%k_dim) :: ff 
     REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w) :: v_loc
     REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w, mesh_data_info%k_dim) :: f_loc
     REAL(KIND = 8), DIMENSION(this%mesh%np) :: v_glb
     REAL(KIND = 8), DIMENSION(this%mesh%me) :: volK
     INTEGER, DIMENSION(this%mesh%gauss%n_w) :: idxm, jj_loc
     INTEGER :: i, k, m, ni, nj, iglob
     Vec                                         :: vect
     PetscErrorCode                              :: ierr

     CALL VecSet(vect, 0.d0, ierr)
     v_glb = 0.d0
     DO m = 1, this%mesh%dom_me
        jj_loc = this%mesh%jj(:, m)
        f_loc = ff(jj_loc,:)
        !<==recompute cij on the fly
        v_loc = 0.d0
        DO ni = 1, this%mesh%gauss%n_w
           DO k = 1, this%mesh%gauss%k_d
              DO nj = 1, this%mesh%gauss%n_w
                 v_loc(ni) = v_loc(ni) + f_loc(nj,k)* &
                 SUM(this%mesh%gauss%dw(k,nj,:,m)*this%mesh%gauss%ww(ni,:)*this%mesh%gauss%rj(:,m))
              ENDDO
           ENDDO
        ENDDO
        idxm = this%LA%loc_to_glob(1, jj_loc) -1
        v_loc = -v_loc
        CALL VecSetValues(vect, this%mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
!!$
!!$        v_glb(jj_loc) = v_glb(jj_loc) - v_loc

     ENDDO

!!$     CALL VecSetValues(vect, this%mesh%np, this%LA%loc_to_glob(1,1:this%mesh%np)-1, v_glb, INSERT_VALUES, ierr)

     CALL VecAssemblyBegin(vect, ierr)
     CALL VecAssemblyEnd(vect, ierr)
   END SUBROUTINE compute_flux


!!$   SUBROUTINE commutator(un, stiff, LA,  alpha)
!!$     IMPLICIT NONE
!!$     REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: un
!!$     TYPE(petsc_csr_LA) :: LA
!!$     REAL(KIND = 8), DIMENSION(:), INTENT(OUT) :: alpha
!!$     INTEGER :: ierr, nrows, i, nz_cols
!!$     INTEGER, DIMENSION(1) :: iglob
!!$     INTEGER, DIMENSION(1000) :: nz_cols_idx
!!$     REAL(KIND = 8), DIMENSION(1000) :: nz_vals
!!$     REAL(KIND = 8) :: num, den
!!$     Mat :: stiff
!!$     CALL MatGetLocalSize(stiff, nrows, PETSC_NULL_INTEGER, ierr)
!!$     write(*,*) 'nb of rows', nrows
!!$     DO i = 1, nrows
!!$
!!$       !write(*,*) 'TEST'
!!$       !CALL MatGetLocalSize(this%stiffL, k, PETSC_NULL_INTEGER, ierr)
!!$       !write(*,*) 'k', k, mesh%dom_np, mesh%np
!!$       !DO  i = 1, 10  
!!$       !   k = LA%loc_to_glob(1, i) - 1
!!$       !   write(*,*) 'k', k
!!$       !   CALL MatGetRow(this%stiffL, k, nz_cols, nz_cols_idx, nz_vals, ierr)
!!$       !   CALL MatRestoreRow(this%stiffL, k, nz_cols, nz_cols_idx, nz_vals, ierr)
!!$       !END DO
!!$       !stop
!!$big mess
!!$        write(*,*) ' i', i, size(LA%loc_to_glob,1), size(LA%loc_to_glob,2)
!!$        iglob = LA%loc_to_glob(1, i) - 1
!!$        write(*,*) ' verif', iglob(1)
!!$        CALL MatGetRow(stiff, iglob(1), nz_cols, nz_cols_idx, nz_vals, ierr)
!!$        write(*,*) ' verif', nz_cols, LA%glob_to_loc(1, nz_cols_idx+1)
!!$        num = SUM(nz_vals*(un(nz_cols_idx+1)))
!!$        den = SUM(abs(nz_vals*un(nz_cols_idx+1)))
!!$        alpha(i) =num/(max(den,1.d-20))
!!$        CALL MatRestoreRow(stiff, iglob(1), nz_cols, nz_cols_idx, nz_vals, ierr)
!!$     END DO
!!$     stop
!!$   END SUBROUTINE commutator

 END MODULE euler_type_MODULE
