MODULE cell_limiting_engine_parallel_module
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE st_matrix, ONLY : extract_through_ghost, create_my_ghost
   USE read_inputs_module
   USE periodic_data_module

   TYPE argument_limiting_type
      CHARACTER(LEN=rec_length) :: if_limiting       = '=== Apply cell-limiting (T/F) ? ==='
      CHARACTER(LEN=rec_length) :: limit_max         = '=== How many limiting iterations? ==='
   END TYPE argument_limiting_type

   TYPE limiting_type
      CHARACTER(100) :: name
      LOGICAL                   :: if_limiting       = .False.
      INTEGER                   :: limit_max         = 2
      INTEGER, DIMENSION(:,:), POINTER :: jj
      REAL(KIND=8) :: mass_eps
      REAL(KIND=8) :: epsilon = 1.d-8
      REAL(KIND = 8), DIMENSION(:,:), POINTER :: localized_mass
      REAL(KIND = 8), DIMENSION(:),   POINTER :: lumped_mass
      TYPE(petsc_csr_LA),  POINTER :: LA
      TYPE(periodic_type), POINTER :: per
      Vec, PRIVATE :: xvect1, xvect2, x_ghost

   CONTAINS
      PROCEDURE, PUBLIC  :: init => init_limiting
      PROCEDURE, PUBLIC  :: read => read_limiting_data
      PROCEDURE, PUBLIC  :: iterative_cell_limiting_procedure
   END TYPE limiting_type

   ABSTRACT INTERFACE
      FUNCTION template_zero_of_psi_vec(psi_m,u0,P) RESULT(v)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: u0, P
         REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: psi_m
         REAL(KIND=8), DIMENSION(SIZE(psi_m)) :: v
      END FUNCTION template_zero_of_psi_vec

      FUNCTION template_psi_vec(x,psi_m) RESULT(v)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: x
         REAL(KIND=8), DIMENSION(:), INTENT(IN) :: psi_m
         REAL(KIND=8), DIMENSION(SIZE(psi_m)) :: v
      END FUNCTION template_psi_vec
   END INTERFACE

   TYPE argument_limiting_functional_type
      CHARACTER(LEN=rec_length) :: char_relaxation_method   = '=== Relaxation method (none/avg/minmod) ==='
      CHARACTER(LEN=rec_length) :: relaxation_global_bound = '=== Global bounds for relaxation ==='
   END TYPE argument_limiting_functional_type

   TYPE ::  limiting_functional_type
      PROCEDURE(template_psi_vec),         POINTER, NOPASS :: psi => NULL()
      PROCEDURE(template_zero_of_psi_vec), POINTER, NOPASS :: zero_of_psi => NULL()
      CHARACTER(LEN=rec_length) :: name
      REAL(KIND=8)              :: relaxation_global_bound = -1.d10
      CHARACTER(len=rec_length) :: char_relaxation_method   = 'none'
      INTEGER                   :: relaxation_method
   CONTAINS
      PROCEDURE, PRIVATE :: read => read_limiting_functional_data 
   END TYPE limiting_functional_type

   TYPE :: limiting_all_functional_type
      TYPE(limiting_functional_type), DIMENSION(:), POINTER :: limiting_functionals
      TYPE(petsc_csr_LA),      POINTER :: LA
      TYPE(mesh_type), POINTER      :: mesh
      INTEGER                       :: nl
      CHARACTER(LEN=:), ALLOCATABLE :: name
      Vec :: xvect1, xvect2, xvect3, stiffness_RowSumAbs, avg_row_sum, ones_row_sum
      Mat :: stiffness, ones, avg
   CONTAINS
      PROCEDURE, PUBLIC :: init => init_limiting_functionals
      PROCEDURE :: RELAX_BOUNDS
      PROCEDURE :: init_mat_relax_avg, init_mat_relax_minmod
   END TYPE limiting_all_functional_type
   
   INTEGER, PRIVATE, PARAMETER :: RELAX_NONE=1, RELAX_AVG=2, RELAX_MINMOD=3
   CHARACTER(LEN=20), DIMENSION(3), PRIVATE, PARAMETER  :: list_relax_method = &
               [CHARACTER(LEN=20) :: 'none', 'avg', 'minmod']

CONTAINS

!=====================================
!======= Limiting initialization =====
!=====================================

   SUBROUTINE init_limiting(this, name, mesh, LA)
#include "petsc/finclude/petsc.h"
      USE petsc 
      USE solver_petsc
      USE def_type_mesh
      USE compute_periodic
      USE fem_M
      USE fem_petsc_matrix_factory_module
      IMPLICIT NONE
      CLASS(limiting_type),    INTENT(INOUT) :: this
      CHARACTER(LEN=*),             INTENT(IN) :: name
      TYPE(mesh_type),    TARGET, INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
      REAL(KIND=8), DIMENSION(mesh%np)         :: vol_of_Ti
      REAL(KIND=8), DIMENSION(mesh%gauss%n_w)  :: vol_of_Ti_loc
      INTEGER, DIMENSION(SIZE(mesh%jj,1))      :: idxm

      INTEGER, POINTER, DIMENSION(:) :: ifrom
      INTEGER :: m, n, ierr, k
      REAL(KIND=8) :: volK
      Mat :: mass
    
      !===Start reading limiting data
      this%name = name
      CALL this%read("LIMITING PARAMETERS FOR "//TRIM(ADJUSTL(name)))
      this%jj => mesh%jj
      this%LA => LA
      this%per => mesh%per

    !===Petsc ghosting for cell-averaging
      CALL create_my_ghost(mesh, LA, ifrom)
      CALL VecCreateGhost(mesh%comm, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%xvect1, ierr)
      CALL VecGhostGetLocalForm(this%xvect1, this%x_ghost, ierr)
      CALL VecDuplicate(this%xvect1, this%xvect2, ierr)
      !===Compute lumped masse
      CALL create_local_petsc_matrix(mesh%comm, LA, mass, clean = .FALSE.)
      ALLOCATE(this%lumped_mass(mesh%np))
      CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, mass)
      CALL construct_lumped_mass_vector(mesh, LA, mass, this%xvect1)
      CALL periodic_add_vector_petsc(mesh%per%nb_bords, mesh%per%list, mesh%per%perlist, this%xvect1, LA)
      CALL extract_through_ghost(this%xvect1, 1, 1, this%LA, this%lumped_mass, opt_assemble=.TRUE.)
   
      !===Localized mass construction
      ALLOCATE(this%localized_mass(mesh%gauss%n_w,mesh%me))
      this%localized_mass = 0.d0

!VB CORRECTED VERSION WHEN SEVERAL PROCESSES
      vol_of_Ti_loc = 0.d0
      CALL VecZeroEntries(this%xvect1, ierr)
      DO m = 1, mesh%me
         volK = SUM(mesh%gauss%rj(:,m))
         DO n = 1, mesh%gauss%n_w
            vol_of_Ti_loc(n) = volK
         END DO
         idxm = this%LA%loc_to_glob(1, this%jj(:,m)) -1
         CALL VecSetValues(this%xvect1, mesh%gauss%n_w, idxm, vol_of_Ti_loc, ADD_VALUES, ierr)
      END DO
      CALL VecAssemblyBegin(this%xvect1, ierr)
      CALL VecAssemblyEnd(this%xvect1, ierr)

      CALL periodic_add_vector_petsc(mesh%per%nb_bords, mesh%per%list, mesh%per%perlist, this%xvect1, LA)
      CALL extract_through_ghost(this%xvect1, 1, 1, this%LA, vol_of_Ti, opt_assemble=.FALSE.)
!VB CORRECTED VERSION WHEN SEVERAL PROCESSES

      DO m = 1, mesh%me
         volK = SUM(mesh%gauss%rj(:,m))
         DO n = 1, mesh%gauss%n_w
            this%localized_mass(n,m) = this%lumped_mass(mesh%jj(n,m))*volK/vol_of_Ti(mesh%jj(n,m))
         END DO
      END DO

      this%mass_eps = this%epsilon*SUM(this%lumped_mass)/mesh%np
   
      CALL MatDestroy(mass, ierr)

   END SUBROUTINE init_limiting

   SUBROUTINE init_mat_relax_avg(this)
      USE solver_petsc, ONLY: create_local_petsc_matrix
      IMPLICIT NONE
      CLASS(limiting_all_functional_type), INTENT(INOUT) :: this
      INTEGER :: i, j, p, dom_np, ierr, nb
      INTEGER, DIMENSION(1) :: idx, jdx
      REAL(KIND=8), DIMENSION(1,1) :: mat_loc

      CALL VecDuplicate(this%xvect1, this%avg_row_sum, ierr)
      CALL create_local_petsc_matrix(this%mesh%comm, this%LA, this%avg, clean = .FALSE.)

      dom_np = SIZE(this%LA%ia)-1
      DO i=0, dom_np-1
         idx = this%LA%loc_to_glob(1, i+1)-1
         nb = this%LA%ia(i+1) - this%LA%ia(i) - 1
         DO p=this%LA%ia(i), this%LA%ia(i+1)-1
            j = this%LA%ja(p)
            jdx = j
            IF (idx(1) /= jdx(1)) THEN
               mat_loc = 1.d0
               CALL MatSetValues(this%avg, 1, idx, 1, jdx, mat_loc, INSERT_VALUES, ierr)
            ELSE
               mat_loc = nb
               CALL MatSetValues(this%avg, 1, idx, 1, jdx, mat_loc, INSERT_VALUES, ierr)
            END IF
         END DO
      END DO

      CALL MatAssemblyBegin(this%avg, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(this%avg, MAT_FINAL_ASSEMBLY, ierr)
      CALL VecSet(this%xvect1, 1.d0, ierr)
      CALL MatMult(this%avg, this%xvect1, this%avg_row_sum, ierr)

   END SUBROUTINE init_mat_relax_avg

   SUBROUTINE init_mat_relax_minmod(this)
      USE solver_petsc, ONLY: create_local_petsc_matrix
      IMPLICIT NONE
      CLASS(limiting_all_functional_type), INTENT(INOUT) :: this
      INTEGER :: i, j, p, dom_np, ierr, nb
      INTEGER, DIMENSION(1) :: idx, jdx
      REAL(KIND=8), DIMENSION(1,1) :: mat_loc

      CALL VecDuplicate(this%xvect1, this%ones_row_sum, ierr)
      CALL create_local_petsc_matrix(this%mesh%comm, this%LA, this%ones, clean = .FALSE.)

      dom_np = SIZE(this%LA%ia)-1
      DO i=0, dom_np-1
         idx = this%LA%loc_to_glob(1, i+1)-1
         nb = this%LA%ia(i+1) - this%LA%ia(i) - 1
         DO p=this%LA%ia(i), this%LA%ia(i+1)-1
            j = this%LA%ja(p)
            jdx = j
            mat_loc = 1.d0
            CALL MatSetValues(this%ones, 1, idx, 1, jdx, mat_loc, INSERT_VALUES, ierr)
         END DO
      END DO

      CALL MatAssemblyBegin(this%ones, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(this%ones, MAT_FINAL_ASSEMBLY, ierr)
      CALL VecSet(this%xvect1, 1.d0, ierr)
      CALL MatMult(this%ones, this%xvect1, this%ones_row_sum, ierr)

   END SUBROUTINE init_mat_relax_minmod

   SUBROUTINE init_limiting_functionals(this, limiting_functionals, name, LA, mesh)
#include "petsc/finclude/petsc.h"
      USE petsc 
      USE solver_petsc
      USE def_type_mesh
      USE compute_periodic
      USE fem_M
      USE fem_petsc_matrix_factory_module
      IMPLICIT NONE
      CLASS(limiting_all_functional_type) :: this
      TYPE(limiting_functional_type), DIMENSION(:), TARGET, INTENT(IN) :: limiting_functionals
      TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
      TYPE(mesh_type),    TARGET, INTENT(IN) :: mesh
      CHARACTER(LEN=*), INTENT(IN) :: name
      INTEGER, POINTER, DIMENSION(:) :: ifrom, list_relax
      INTEGER :: i, ierr

      !=== BUILD POINTERS
      this%name = name
      this%mesh => mesh
      this%LA   => LA
      this%limiting_functionals => limiting_functionals
      this%nl = SIZE(this%limiting_functionals)
      IF (this%nl == 0) RETURN


      !=== Read data
      ALLOCATE(list_relax(this%nl))
      DO i=1, this%nl
         CALL this%limiting_functionals(i)%read
         list_relax(i) = this%limiting_functionals(i)%relaxation_method
      END DO

      !=== BUILD MATRICES/VECTORS
      CALL create_my_ghost(this%mesh, LA, ifrom)
      CALL VecCreateGhost(this%mesh%comm, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%xvect1, ierr)
      CALL VecDuplicate(this%xvect1, this%xvect2, ierr)
      CALL VecDuplicate(this%xvect1, this%xvect3, ierr)

      IF (ANY(RELAX_AVG==list_relax) .OR. ANY(RELAX_MINMOD==list_relax)) THEN
         CALL create_local_petsc_matrix(this%mesh%comm, LA, this%stiffness, clean = .FALSE.)
         CALL qs_mass_diff_M (this%mesh, 0.d0, 1.d0, LA, this%stiffness)
         ! CALL periodic_matrix_petsc(mesh%per, LA, this%stiffness) !<-- FIXME PERIODICITY
         CALL VecDuplicate(this%xvect1, this%stiffness_RowSumAbs, ierr)
         CALL MatGetRowSumAbs(this%stiffness, this%stiffness_RowSumAbs, ierr)
         CALL MatGetDiagonal(this%stiffness, this%xvect1, ierr)
         CALL VecAXPY(this%stiffness_RowSumAbs, -1.d0, this%xvect1, ierr)
      END IF

      IF (ANY(RELAX_AVG==list_relax)) THEN
         CALL this%init_mat_relax_avg
      END IF 

      IF (ANY(RELAX_MINMOD==list_relax)) THEN
         CALL this%init_mat_relax_minmod
      END IF 

   END SUBROUTINE init_limiting_functionals

   SUBROUTINE read_limiting_data(this, section_name)
      USE read_inputs_module
      IMPLICIT NONE
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name
      CLASS(limiting_type),    INTENT(INOUT) :: this
      TYPE(argument_limiting_type)           :: argument_data
      !================
      !=== MANDATORY Reading all data file
      !================
      IF (PRESENT(section_name)) THEN
         CALL read_data_init_list(section_name)
      ELSE
         CALL read_data_init_list()
      END IF

      !================
      !=== We now find the relevant information for this specific limiting data
      !================
      !===if_limiting
      CALL read_data(argument_data%if_limiting , this%if_limiting, &
      opt_name=this%name)

      !===Number of limiting iterations
      CALL read_data(argument_data%limit_max , this%limit_max, &
      opt_name=this%name, opt_add=this%if_limiting)

      !================
      !=== MANDATORY to close data for the current section and
      !=== rewrite it with new information for the next sections
      !================
      CALL finalize_rewrite_data
   END SUBROUTINE read_limiting_data

   SUBROUTINE read_limiting_functional_data(this, section_name)
      USE read_inputs_module
      USE my_util
      IMPLICIT NONE
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name
      CLASS(limiting_functional_type),    INTENT(INOUT) :: this
      TYPE(argument_limiting_functional_type)           :: argument_data
      !================
      !=== MANDATORY Reading all data file
      !================
      ! IF (PRESENT(section_name)) THEN   !<== FIXME
      !    CALL read_data_init_list(section_name)
      ! ELSE
         CALL read_data_init_list()
      ! END IF

      !================
      !=== We now find the relevant information for this specific limiting data
      !================


      !===relaxation_method
      CALL read_data(argument_data%char_relaxation_method, this%char_relaxation_method, &
      opt_name=this%name)
      CALL get_tab_idx_char(this%char_relaxation_method, list_relax_method, this%relaxation_method)

      !===relaxation_global_bound
      CALL read_data(argument_data%relaxation_global_bound, this%relaxation_global_bound, &
      opt_name=this%name) 

      !================
      !=== MANDATORY to close data for the current section and
      !=== rewrite it with new information for the next sections
      !================
      CALL finalize_rewrite_data
   END SUBROUTINE read_limiting_functional_data

!=====================================
!======= Limiting procedures =========
!=====================================

   SUBROUTINE iterative_cell_limiting_procedure(this, xx_in, loc_min, lim_bounds, xx_out)  
#include "petsc/finclude/petsc.h"
      USE petsc 
      USE compute_periodic
      USE my_util, ONLY: error_petsc
      IMPLICIT NONE
      CLASS(limiting_type),             INTENT(IN) :: this
      CLASS(limiting_functional_type),  INTENT(IN) :: lim_bounds
      PROCEDURE(template_zero_of_psi_vec), POINTER :: zero_of_psi
      PROCEDURE(template_psi_vec)        , POINTER :: psi
      REAL(KIND=8), DIMENSION(:,:),                         INTENT(IN) :: xx_in
      REAL(KIND=8), DIMENSION(SIZE(xx_in,1),SIZE(xx_in,2)), INTENT(OUT):: xx_out
      REAL(KIND=8), DIMENSION(:)             :: loc_min
      ! REAL(KIND=8), DIMENSION(1)               :: uk_minus, uk_plus
      REAL(KIND=8), DIMENSION(SIZE(xx_in,2))               :: uk_minus, uk_plus
      ! REAL(KIND=8), DIMENSION(k_dim+1,1)    :: xx_loc, UU, PP
      REAL(KIND=8), DIMENSION(k_dim+1,SIZE(xx_in,2))    :: xx_loc, UU, PP
      REAL(KIND=8), DIMENSION(k_dim+1) :: lambda_minus, lambda_plus, dummy
      REAL(KIND=8), DIMENSION(k_dim+1) :: loc_min_loc, loc_min_up_vec, loc_min_down_vec
      INTEGER,      DIMENSION(k_dim+1) :: jloc
      INTEGER,      DIMENSION(k_dim+1) :: limit_zero, limit_plus, limit_minus
      INTEGER :: k, m, n, me, nw, syst_size, iminus, iplus, comp, np, i, ierr
      LOGICAL, DIMENSION(k_dim+1) :: mask_up, mask_down
      REAL(KIND=8) :: mass_plus, mass_minus, &
            lambda_K_minus, lambda_K_plus, &
            lambda_star_minus, lambda_star_plus


      zero_of_psi => lim_bounds%zero_of_psi
      psi         => lim_bounds%psi

      np = SIZE(xx_in, 1)
      me = SIZE(this%jj,2)
      nw = SIZE(this%jj,1)
      syst_size = SIZE(xx_in,2)
      xx_out = 0.d0

      DO m = 1, me
         limit_plus = 0
         limit_minus = 0
         jloc = this%jj(:,m)
         DO k = 1, syst_size
            xx_loc(:,k) = xx_in(jloc,k)
         END DO
         loc_min_loc = loc_min(jloc)

         loc_min_down_vec = loc_min_loc(:) - this%epsilon*ABS(loc_min_loc(:))
         loc_min_up_vec   = loc_min_loc(:) + this%epsilon*ABS(loc_min_loc(:))
         
         mask_up = psi(xx_loc,loc_min_up_vec)>0
         mask_down = psi(xx_loc,loc_min_down_vec)<0

         WHERE(mask_down)
            limit_minus(:) = 1
         ELSEWHERE(mask_up)
            limit_plus(:) = 1
         ENDWHERE
         limit_zero = 1 - limit_minus - limit_plus

         DO k = 1, syst_size
            uk_minus(k)=SUM(this%localized_mass(:,m)*xx_loc(:,k)*limit_minus)
            uk_plus(k)=SUM(this%localized_mass(:,m)*xx_loc(:,k)*limit_plus)
         END DO
         mass_minus = SUM(this%localized_mass(:,m)*limit_minus)
         mass_plus  = SUM(this%localized_mass(:,m)*limit_plus)
         
         uk_minus = uk_minus/max(mass_minus,this%mass_eps)
         uk_plus  = uk_plus/max(mass_plus,this%mass_eps)

         DO k = 1, syst_size
            UU(:,k) = uk_plus(k)
            PP(:,k) = xx_loc(:,k)-uk_plus(k)
         END DO
         lambda_minus = zero_of_psi(loc_min_loc,UU,PP)
         WHERE(limit_minus==0)
            lambda_minus = 1.d0
         ENDWHERE
         
         DO k = 1, syst_size
            UU(:,k) = xx_loc(:,k)
            PP(:,k) = uk_minus(k)-xx_loc(:,k)
         END DO
         lambda_plus = zero_of_psi(loc_min_loc,UU,PP)
         WHERE(limit_plus==0)
            lambda_plus = 1.d0
         ENDWHERE

         lambda_minus = MAX(MIN(lambda_minus,1.d0),0.d0)
         lambda_plus  = MAX(MIN(lambda_plus,1.d0),0.d0)
         Lambda_star_minus = MINVAL(lambda_minus)
         Lambda_star_plus  = MINVAL(lambda_plus)
         Lambda_K_minus = MAX(Lambda_star_minus, 1.d0-Lambda_star_plus*mass_plus/mass_minus)
         Lambda_K_plus  = MIN(Lambda_star_plus, (1.d0-Lambda_star_minus)*mass_minus/mass_plus)
         !=== DEBUGGING ===!
         !write(*,*)  'm possible limiting', m
         !write(*,*) lambda_minus,  lambda_plus
         !write(*,*) lambda_star_minus,  lambda_star_plus
         !write(*,*)  Lambda_K_minus, Lambda_K_plus
         !=== DEBUGGING ===!
         !    ! !!$ ===P2 fix
         !    ! IF (ABS(this%lumped_mass(jloc(n))).LE.this%mass_eps) THEN
         !    !    xx(n,m,:) = uk_plus(:)
         !    ! ELSE
         !    ! !!$ ===END fix
         DO k=1, syst_size
               dummy = xx_loc(:,k)  &
                     +limit_minus(:)*(1-Lambda_K_minus)*(uk_plus(k)-xx_loc(:,k))&
                     +limit_plus(:) *     Lambda_K_plus*(uK_minus(k)-xx_loc(:,k)) 
               ! dummy  = limit_zero(:)*xx_loc(:,k) &
               ! +limit_minus(:)*(xx_loc(:,k)+(1-Lambda_K_minus)*(uk_plus(k)-xx_loc(:,k)))&
               ! +limit_plus(:) *(xx_loc(:,k)+     Lambda_K_plus*(uK_minus(k)-xx_loc(:,k)))
               xx_out(jloc,k) = xx_out(jloc,k) + dummy*this%localized_mass(:,m)
         END DO
      END DO

   !===Now we average over the nodes=========

      DO k = 1, syst_size
         CALL VecZeroEntries(this%xvect1, ierr)
         CALL VecSetValues(this%xvect1, np, this%LA%loc_to_glob(1, :) - 1, xx_out(:,k), ADD_VALUES, ierr)
         CALL VecAssemblyBegin(this%xvect1, ierr)
         CALL VecAssemblyEnd(this%xvect1, ierr)

         CALL periodic_add_vector_petsc(this%per%nb_bords, this%per%list, this%per%perlist, this%xvect1, this%LA)
         CALL extract_through_ghost(this%xvect1, 1, 1, this%LA, xx_out(:,k), opt_assemble=.TRUE.)
         xx_out(:,k) = xx_out(:,k)/this%lumped_mass
      END DO    
      
   END SUBROUTINE iterative_cell_limiting_procedure

   SUBROUTINE RELAX_BOUNDS(this, state, psi_min, nth_fonctional)
      USE petsc_tools, ONLY: array_to_petsc_vec
      IMPLICIT NONE
      CLASS(limiting_all_functional_type),  INTENT(IN) :: this
      REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: state
      INTEGER,                              INTENT(IN) :: nth_fonctional
      REAL(KIND=8), DIMENSION(:),       INTENT(INOUT)  :: psi_min
      REAL(KIND=8) :: norm
      ! REAL(KIND=8), INTENT(IN)               :: glob_min, glob_max
      REAL(KIND=8), DIMENSION(SIZE(state, 1)) :: alpha, denom, denom_ref, mask
      INTEGER,      DIMENSION(SIZE(state, 1)) ::   beta
      REAL(KIND=8), DIMENSION(SIZE(state, 1)) ::   un 
      REAL(kind=8), DIMENSION(SIZE(state, 1)) :: dummy
      INTEGER      :: i, j, j_loc, p, ni, nj, m, me, nw, n, np, dom_np, ierr

      IF(this%limiting_functionals(nth_fonctional)%relaxation_method == RELAX_NONE) THEN
         psi_min = MAX(this%limiting_functionals(nth_fonctional)%relaxation_global_bound,psi_min)
         RETURN
      END IF

      np = SIZE(state, 1)
      dummy = 0.d0
      un(:) = this%limiting_functionals(nth_fonctional)%psi(state(:, :), dummy)

      me = SIZE(this%mesh%jj,2)
      nw = SIZE(this%mesh%jj,1)
      alpha = 0.d0
      beta = 0
      CALL array_to_petsc_vec(un, this%xvect1, this%LA, "insert")
      CALL MatMult(this%stiffness, this%xvect1, this%xvect2, ierr)
      CALL VecPointWiseDivide(this%xvect2, this%xvect2, this%stiffness_RowSumAbs, ierr)
      

      SELECT CASE(this%limiting_functionals(nth_fonctional)%relaxation_method)
      CASE(RELAX_AVG) !==Average
         CALL MatMult(this%avg, this%xvect2, this%xvect1, ierr)
         CALL VecPointWiseDivide(this%xvect1, this%xvect1, this%avg_row_sum, ierr)
         CALL extract_through_ghost(this%xvect1, 1, 1, this%LA, denom, opt_assemble=.FALSE.)
         
         !TESTTTT
         ! CALL extract_through_ghost(this%avg_row_sum, 1, 1, this%LA, alpha, opt_assemble=.FALSE.)
         ! dom_np = SIZE(this%LA%ia) - 1
         ! norm = MAXVAL((alpha(1:dom_np) - (this%LA%ia(1:) - this%LA%ia(0:dom_np-1) - 1)))
         ! WRITE(*,*) "max norm alpha = ", norm 
         ! norm = MINVAL((alpha(1:dom_np) - (this%LA%ia(1:) - this%LA%ia(0:dom_np-1) - 1)))
         ! WRITE(*,*) "min norm alpha = ", norm 
         ! WRITE(*,*) "mean alpha = ", SUM(alpha)/SIZE(alpha)
         ! WRITE(*,*) "dif ia", (this%LA%ia(1:) - this%LA%ia(0:dom_np-1))
         ! stop
         !TESTTTT
      CASE(RELAX_MINMOD) !===Minmod
         dom_np = SIZE(this%LA%ia)-1
         CALL extract_through_ghost(this%xvect2, 1, 1, this%LA, alpha, opt_assemble=.FALSE.)
         CALL VecNorm(this%xvect2, NORM_1, norm, ierr) 
         norm = norm * 1.d-30 ! <=== epsilon
         CALL VecSet(this%xvect3, norm, ierr)
         CALL VecPointwiseMaxAbs(this%xvect1, this%xvect2, this%xvect3, ierr) ! <===MAX(ABS(alpha), epsilon)(:)
         CALL VecPointWiseDivide(this%xvect1, this%xvect2, this%xvect1, ierr) ! <=== sign(alpha) = alpha / MAX(ABS(alpha), epsilon)(:)
         CALL MatMult(this%ones, this%xvect1, this%xvect3, ierr)
         CALL VecAbs(this%xvect3, ierr)
         CALL VecAYPX(this%xvect3, -1.d0, this%ones_row_sum, ierr)

         CALL extract_through_ghost(this%xvect3, 1, 1, this%LA, mask, opt_assemble=.FALSE.)

         !TEST
         WHERE(ABS(mask)>1.d-1)
            mask = 0.d0
         ELSEWHERE
            mask = 1.d0
         ENDWHERE
         denom = ABS(alpha)
         DO m = 1, me
            DO ni = 1, nw
               i = this%mesh%jj(ni,m)
               denom(i) = MIN(denom(i),MINVAL(ABS(alpha(this%mesh%jj(:,m)))))
            END DO
         END DO
         denom = denom*mask
         CALL array_to_petsc_vec(denom, this%xvect3, this%LA, "min")
         CALL extract_through_ghost(this%xvect3, 1, 1, this%LA, denom, opt_assemble=.FALSE.)
      CASE DEFAULT
         WRITE(*,*) ' BUG in relax, method not implemented: ', TRIM(ADJUSTL(this%limiting_functionals(nth_fonctional)%char_relaxation_method))
         STOP
      END SELECT
      psi_min = psi_min - 4.*ABS(denom)
      psi_min = MAX(this%limiting_functionals(nth_fonctional)%relaxation_global_bound,psi_min)
   END SUBROUTINE RELAX_BOUNDS

END MODULE cell_limiting_engine_parallel_module
