MODULE euler_type_module
!>> limited global uses to avoid unexpected behaviors
#include "petsc/finclude/petsc.h"
   USE petsc
   USE abstract_hyperbolic_module, ONLY: hyperbolic_type
   USE euler_bc_arrays, ONLY : euler_bc_type
   USE petsc_tools,     ONLY : array_to_petsc_vec
   USE Butcher_tableau
   USE cell_limiting_engine_parallel_module, ONLY : limiting_type
   USE def_type_mesh, ONLY : mesh_type, petsc_csr_LA
   USE read_inputs_module,    ONLY : rec_length
!>> limited global uses to avoid unexpected behaviors

   IMPLICIT NONE

   ABSTRACT INTERFACE
      FUNCTION function_template_pressure(rho, ie) RESULT(vv)
         REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, ie
         REAL(KIND = 8), DIMENSION(SIZE(rho, 1))  :: vv
      END FUNCTION function_template_pressure
   END INTERFACE

   TYPE argument_euler_type

      CHARACTER(LEN=rec_length) :: no_iter   = '=== No iteration for lambda solver? (t/f) ==='
      CHARACTER(LEN=rec_length) :: in_tol    = '=== Tolerance for lambda solver ==='
      CHARACTER(LEN=rec_length) :: eos_param = '=== b_covolume? ==='
   END TYPE argument_euler_type

   TYPE, EXTENDS(hyperbolic_type)  :: euler_type
      !===Parameters read from data
      REAL(KIND = 8), DIMENSION(1) :: eos_param = 0.d0
      REAL(KIND = 8)               :: in_tol = 1.d-2
      LOGICAL                      :: no_iter = .TRUE.
      TYPE(euler_bc_type)          :: bc
      CHARACTER(len=5), DIMENSION(:), ALLOCATABLE :: name_comp
      PROCEDURE(function_template_pressure),  NOPASS, POINTER :: pressure => NULL()
   CONTAINS
      PROCEDURE, PUBLIC  :: init_euler
      PROCEDURE, PRIVATE :: read_euler_data
      PROCEDURE :: flux           => flux_euler
      PROCEDURE :: compute_lambda => lambda_euler
      PROCEDURE :: construct_bc => construct_euler_bc
      PROCEDURE :: impose_bc => impose_bc_euler
   END TYPE euler_type

CONTAINS
   SUBROUTINE init_euler(this, name)
      USE space_dim
      USE my_util, ONLY: error_petsc, to_str
      IMPLICIT NONE
      CLASS(euler_type), INTENT(INOUT) :: this
      CHARACTER(LEN=*),  INTENT(IN)    :: name

      this%name = name
      this%syst_dim = k_dim + 2
      ALLOCATE(this%name_comp(this%syst_dim))

      this%name_comp(1) = 'rho'
      this%name_comp(k_dim + 2) = 'E'

      SELECT CASE(k_dim)
      CASE(1)
         this%name_comp(2) = 'ux'
      CASE(2)
         this%name_comp(2) = 'ux'
         this%name_comp(3) = 'uy'
      CASE DEFAULT
         CALL error_petsc("BUG in init_euler: Wrong k_dim "//to_str(k_dim))
      END SELECT
      CALL this%read_euler_data(trim(adjustl(name))//" PARAMETERS")

   END SUBROUTINE init_euler

   SUBROUTINE read_euler_data(this, section_name)
     USE read_inputs_module
     IMPLICIT NONE
     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name

     CLASS(euler_type), INTENT(INOUT) :: this
     TYPE(argument_euler_type)        :: argument_data


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

     !===b_covolume
     CALL read_data(argument_data%eos_param, this%eos_param(1), opt_name=this%name)

     !===no_iter for lambda
     CALL read_data(argument_data%no_iter, this%no_iter, opt_name=this%name)

     !===tol for lambda
     CALL read_data(argument_data%in_tol, this%in_tol, opt_name=this%name)

     !================
     !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
     !================

     CALL finalize_rewrite_data
   END SUBROUTINE read_euler_data

!====================================================================
!====================================================================
!====== MANDATORY PROCEDURES FOR DEFINING HYPERBOLIC OBJECT =========
!====================================================================
!====================================================================

   FUNCTION flux_euler(this, comp, un) RESULT(vv)  
      USE space_dim
      USE my_util, ONLY : error_petsc, to_str
      IMPLICIT NONE
      CLASS(euler_type),               INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
      INTEGER,                         INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(SIZE(un, 1), k_dim) :: vv

      REAL(KIND = 8), DIMENSION(SIZE(un, 1))      :: H, u, ie
      INTEGER :: k

      SELECT CASE(comp)
      CASE(1)
         DO k = 1, k_dim
            vv(:, k) = un(:, k + 1)
         END DO
      CASE(2:k_dim + 1)
         u = un(:, comp) / un(:, 1)
         DO k = 1, k_dim
            vv(:, k) = un(:, k + 1) * u
         END DO
         ie = un(:, k_dim + 2) / un(:, 1)
         DO k = 1, k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO
         vv(:, comp - 1) = vv(:, comp - 1) + this%pressure(un(:, 1), ie)
      CASE(k_dim + 2)
         ie = un(:, k_dim + 2) / un(:, 1)
         DO k = 1, k_dim
            ie = ie - 0.5d0 * (un(:, k + 1) / un(:, 1))**2
         END DO

         H = un(:, comp) + this%pressure(un(:, 1), ie)
         DO k = 1, k_dim
            vv(:, k) = (un(:, k + 1) / un(:, 1)) * H
         END DO
      CASE DEFAULT
         CALL error_petsc(' BUG in flux, wrong comp = '//to_str(comp)//" with k_dim="//to_str(k_dim))
      END SELECT
   END FUNCTION flux_euler


   SUBROUTINE lambda_euler(this, un, i, j, lambda_max)
      USE arbitrary_eos_lambda_module
      USE space_dim
      IMPLICIT NONE
      CLASS(euler_type),                                 INTENT(INOUT) :: this
      REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(IN) :: un
      INTEGER,                                              INTENT(IN) :: i, j
      REAL(KIND=8), DIMENSION(2),                          INTENT(OUT) :: lambda_max

      INTEGER, DIMENSION(1) :: i_t, j_t
      INTEGER :: k, ierr
      REAL(KIND = 8), DIMENSION(1, k_dim) :: nij_c
      REAL(KIND = 8), DIMENSION(2) :: u, rho, ie, p
      REAL(KIND = 8) :: pstar


      i_t = i
      j_t = j

      DO k = 1, k_dim
         CALL MatGetValues(this%matrices%nij_loc(k), 1, i_t - 1, 1, j_t - 1, nij_c(:, k), ierr)
      END DO

      rho(1) = un(i, 1)
      rho(2) = un(j, 1)

      u(1) = SUM(un(i, 2:1 + k_dim) * nij_c(1, :)) / rho(1)
      u(2) = SUM(un(j, 2:1 + k_dim) * nij_c(1, :)) / rho(2)

      ie(1) = un(i, k_dim + 2) / rho(1) - 0.5d0 * SUM(un(i, 2:1 + k_dim)**2) / rho(1)**2
      ie(2) = un(j, k_dim + 2) / rho(2) - 0.5d0 * SUM(un(j, 2:1 + k_dim)**2) / rho(2)**2

      p = this%pressure(rho, ie)
      CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, lambda_max, pstar)

   END SUBROUTINE lambda_euler

   SUBROUTINE construct_euler_bc(this, mesh, LA)
      USE petsc
#include "petsc/finclude/petsc.h"
      !  USE sub_plot
   
      USE space_dim,           ONLY: k_dim
      IMPLICIT NONE
      CLASS(euler_type), INTENT(INOUT)        :: this
      TYPE(mesh_type)                            :: mesh
      TYPE(petsc_csr_LA)                         :: LA

      CALL this%bc%rho_bc%set(mesh, "density", "DIRICHLET BC PARAMETERS FOR "//TRIM(ADJUSTL(this%name)))

      CALL this%bc%u_bc(1)%set(mesh, "ux")
      
      IF (k_dim>1) THEN
         CALL this%bc%u_bc(2)%set(mesh, "uy")
         CALL this%bc%whole_bdy_bc%set(mesh, "whole boundary")
         CALL this%bc%udotn_bc%set(mesh, "u.n=0")
      END IF
      CALL this%construct_udotn(mesh, LA, this%bc%udotn_bc, this%bc%udotn_normal_vtx)

   END SUBROUTINE construct_euler_bc


   SUBROUTINE impose_bc_euler(this, un, mesh, time)
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(euler_type), INTENT(INOUT) :: this
      TYPE(mesh_type) :: mesh
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
      INTEGER :: comp_sys, comp_u
      REAL(KIND=8), DIMENSION(SIZE(this%bc%udotn_bc%jsd)) :: mdotn

      !=== Simple Dirichlet boundary conditions
      comp_sys = 1
      un(this%bc%rho_bc%jsd, comp_sys) = this%bc%rho_anal(time, mesh%rr(:, this%bc%rho_bc%jsd))

      DO comp_u=1, k_dim
         comp_sys = comp_sys + 1
         un(this%bc%u_bc(comp_sys-1)%jsd, comp_sys) = this%bc%mt_anal(comp_sys - 1, time, mesh%rr(:, this%bc%u_bc(comp_sys-1)%jsd))
      END DO

      comp_sys = comp_sys + 1
      un(this%bc%rho_bc%jsd, comp_sys) = this%bc%E_anal(time, mesh%rr(:, this%bc%rho_bc%jsd))

      !=== u.n boundary conditions
      IF (size(this%bc%udotn_bc%jsd).NE.0) THEN
         mdotn = this%bc%udotn_normal_vtx(:,1)*un(this%bc%udotn_bc%jsd,2) &
         +  this%bc%udotn_normal_vtx(:,2)*un(this%bc%udotn_bc%jsd,3)
         un(this%bc%udotn_bc%jsd,2) = un(this%bc%udotn_bc%jsd,2) - mdotn*this%bc%udotn_normal_vtx(:,1)
         un(this%bc%udotn_bc%jsd,3) = un(this%bc%udotn_bc%jsd,3) - mdotn*this%bc%udotn_normal_vtx(:,2)
         
         mdotn = this%bc%udotn_normal_vtx(:,1)*un(this%bc%udotn_bc%jsd,2) &
         +  this%bc%udotn_normal_vtx(:,2)*un(this%bc%udotn_bc%jsd,3)
      END IF
   END SUBROUTINE impose_bc_euler

!   !GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE
!   !GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE
!   !GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE
!   SUBROUTINE compute_dk (this, un)
!     USE arbitrary_eos_lambda_module
!     USE my_util, ONLY : error_petsc
!     IMPLICIT NONE
!     CLASS(euler_type) :: this
!     REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
!     INTEGER, DIMENSION(1) :: i_t, j_t
!     REAL(KIND = 8), DIMENSION(1, this%mesh%gauss%k_d) :: nij_c
!     REAL(KIND = 8), DIMENSION(1) :: norm_c, dijL_c
!     REAL(KIND = 8), DIMENSION(2) :: u, rho, ie, p, lambda_max
!     LOGICAL, DIMENSION(this%mesh%medge) :: virgin_edge
!     REAL(KIND = 8) :: pstar
!     LOGICAL :: bug
!     INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge, divider, nb_shared_cell
!     nw = this%mesh%gauss%n_w

!     bug = .FALSE.
!     SELECT CASE(this%mesh%gauss%k_d)
!     CASE(1)
!        nb_shared_cell = 1
!        IF (this%mesh%gauss%n_w/=2) bug=.TRUE.
!     CASE(2)
!        nb_shared_cell = 2
!        IF (this%mesh%gauss%n_w/=3) bug=.TRUE.
!     END SELECT
!     IF (bug) THEN
!        CALL error_petsc('Wrong polynomial degree for low-order viscosity')
!     END IF

!     DO m = 1, this%mesh%me
!        DO n = 1, this%mesh%gauss%n_e
!           IF (this%mesh%attr_e(this%mesh%jce(n, m))) THEN
!              edge = this%mesh%jce_loc(n, m)
!              IF (.NOT. virgin_edge(edge)) CYCLE
!              virgin_edge(edge) = .FALSE.
!              ni = MOD(n, nw) + 1
!              nj = MOD(n + 1, nw) + 1
!              i = this%mesh%jj(ni, m)
!              j = this%mesh%jj(nj, m)
!              i_t = i
!              j_t = j
!              DO k = 1, this%mesh%gauss%k_d
!                 CALL MatGetValues(this%matrices%nij_loc(k), 1, i_t - 1, 1, j_t - 1, &
!                      nij_c(:, k), ierr)
!              END DO
!              rho(1) = un(i, 1)
!              rho(2) = un(j, 1)
!              u(1) = SUM(un(i, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(1)
!              u(2) = SUM(un(j, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(2)
!              ie(1) = un(i, this%mesh%gauss%k_d + 2) / rho(1) - 0.5d0 * u(1) * u(1)
!              ie(2) = un(j, this%mesh%gauss%k_d + 2) / rho(2) - 0.5d0 * u(2) * u(2)
!              p = this%pressure(rho, ie)
!              CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, &
!                   lambda_max, pstar)
!              dijL_c = MAXVAL(lambda_max) * norm_c
!              divider = nb_shared_cell

!              IF (this%mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j
!                 DO k = 1, this%mesh%gauss%k_d
!                    CALL MatGetValues(this%matrices%nij_loc(k), 1, j_t - 1, 1, i_t - 1, &
!                         nij_c(:, k), ierr)
!                 END DO
!                 u(1) = SUM(un(i, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(1)
!                 u(2) = SUM(un(j, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(2)
!                 rho = (/rho(2), rho(1)/)
!                 ie = (/ie(2), ie(1)/)
!                 p = (/p(2), p(1)/)
!                 CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, &
!                      lambda_max, pstar)
!                 dijL_c = MAX(dijL_c, MAXVAL(lambda_max) * norm_c)
!                 divider = 1
!              END IF

!              this%matrices%dK(m) = MAX(this%matrices%dK(m),dijL_c(1)/divider)
!           END IF
!        END DO
!     END DO
!   END SUBROUTINE compute_dk

!    SUBROUTINE compute_dt_from_dK(this)
!      IMPLICIT NONE
!      CLASS(euler_type) :: this
!      REAL(KIND = 8), DIMENSION(this%mesh%dom_np) :: dijL_diag
!      REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w) :: v_loc
!      INTEGER, DIMENSION(this%mesh%gauss%n_w) :: idxm
!      INTEGER :: i, m, ni, iglob
!      REAL(KIND = 8) :: dt_min_loc, dt_min_glob
!      Vec                                         :: vect
!      PetscErrorCode                              :: ierr
!      CALL VecSet(vect, 0.d0, ierr)

!      WRITE(*,*) "VB: WARNING (20/04/2026) ==> this subroutine does not call any ghost points???"
!      STOP

!      DO m = 1, this%mesh%me
!         v_loc = 0.d0
!         DO ni = 1, this%mesh%gauss%n_w
!            i = this%mesh%jj(ni, m)
!            iglob = this%LA%loc_to_glob(1, i)
!            idxm(ni) = iglob - 1
!            v_loc(ni) = v_loc(ni) + this%matrices%dK(m)
!         ENDDO
!         CALL VecSetValues(vect, this%mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
!      ENDDO
!      CALL VecAssemblyBegin(vect, ierr)
!      CALL VecAssemblyEnd(vect, ierr)

!      CALL VecGetValues(this%vec_loc, this%mesh%dom_np, this%tab, dijL_diag, ierr)

!      WRITE(*,*) "VB: WARNING (01/05/2026) ==> this subroutine (not used right now) uses lumped_mass.",&
!      " Must be rewritten using lump_mass_vec instead"
!      STOP
!    !   dijL_diag = this%matrices%lumped_mass(1:this%mesh%dom_np) / ABS(dijL_diag)

!      dt_min_loc = MINVAL(dijL_diag) / 2.d0

!      CALL MPI_ALLREDUCE(dt_min_loc, dt_min_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PETSC_COMM_WORLD, ierr)
!      this%dt = this%CFL * dt_min_glob
!    END SUBROUTINE compute_dt_from_dK

!    SUBROUTINE compute_flux(this, ff, Vect)
!       USE space_dim
!       IMPLICIT NONE
!       CLASS(euler_type) :: this
!       REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim) :: ff
!       REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w) :: v_loc
!       REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w, k_dim) :: f_loc
!       REAL(KIND = 8), DIMENSION(this%mesh%np) :: v_glb
!       INTEGER, DIMENSION(this%mesh%gauss%n_w) :: idxm, jj_loc
!       REAL(KIND = 8) :: x
!       INTEGER :: k, m, ni, nj
!       Vec                                         :: vect
!       PetscErrorCode                              :: ierr
!       CALL VecSet(vect, 0.d0, ierr)
!       v_glb = 0.d0
!       DO m = 1, this%mesh%me
!          jj_loc = this%mesh%jj(:, m)
!          f_loc = ff(jj_loc,:)
!          !<==recompute cij on the fly
!          DO ni = 1, this%mesh%gauss%n_w
!             !wwrj = this%mesh%gauss%ww(ni,:)*this%mesh%gauss%rj(:,m)
!             x = 0.d0
!             DO k = 1, this%mesh%gauss%k_d
!                DO nj = 1, this%mesh%gauss%n_w
!                   x = x + f_loc(nj,k)* &
!                         !SUM(this%mesh%gauss%dw(k,nj,:,m)*wwrj)
!                      SUM(this%mesh%gauss%dw(k,nj,:,m)*this%mesh%gauss%ww(ni,:)*this%mesh%gauss%rj(:,m))
!                ENDDO
!             ENDDO
!             v_loc(ni) = x
!          ENDDO
!          idxm = this%LA%loc_to_glob(1, jj_loc) -1
!          v_loc = -v_loc
!          CALL VecSetValues(vect, this%mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
!    !!$        v_glb(jj_loc) = v_glb(jj_loc) - v_loc
!       ENDDO
!    !!$     CALL VecSetValues(vect, this%mesh%np, this%LA%loc_to_glob(1,:)-1, v_glb, INSERT_VALUES, ierr)
!       CALL VecAssemblyBegin(vect, ierr)
!       CALL VecAssemblyEnd(vect, ierr)
!    END SUBROUTINE compute_flux

 END MODULE euler_type_module
