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
      ! CHARACTER(LEN=rec_length) :: if_relax_bounds   = '=== Apply bound relaxation for limiting (T/F) ? ==='
      ! CHARACTER(LEN=rec_length) :: relaxation_method = '=== Relaxation method (avg/minmod) ==='
   END TYPE argument_limiting_type

   TYPE limiting_type
      CHARACTER(100) :: name
      LOGICAL                   :: if_limiting       = .False.
      INTEGER                   :: limit_max         = 2
!       LOGICAL                   :: if_relax_bounds   = .False.
      ! CHARACTER(len=rec_length) :: relaxation_method ='minmod'
      INTEGER, DIMENSION(:,:), POINTER :: jj
      REAL(KIND=8) :: mass_eps
      REAL(KIND=8) :: epsilon = 1.d-8
      REAL(KIND = 8), DIMENSION(:,:), POINTER :: localized_mass
      REAL(KIND = 8), DIMENSION(:),   POINTER :: lumped_mass
      TYPE(petsc_csr_LA),  POINTER :: LA
      TYPE(periodic_type), POINTER :: per
      Vec, PRIVATE :: xvect1, xvect2, x_ghost!, stiffness_RowSumAbs, ones_row_sum
      ! Mat, PRIVATE                 :: stiffness, ones

   CONTAINS
      PROCEDURE, PUBLIC  :: init => init_limiting
      PROCEDURE, PUBLIC  :: read => read_limiting_data
      PROCEDURE, PUBLIC  :: iterative_cell_limiting_procedure
      PROCEDURE, PRIVATE :: cell_averaging!, RELAX_BOUNDS
   END TYPE limiting_type

   ABSTRACT INTERFACE
      FUNCTION template_zero_of_psi(psi_m,u0,P) RESULT(v)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u0, P
         REAL(KIND=8), INTENT(IN)  :: psi_m
         REAL(KIND=8) :: v
      END FUNCTION template_zero_of_psi

      FUNCTION template_psi(x,psi_m) RESULT(v)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
         REAL(KIND=8), INTENT(IN) :: psi_m
         REAL(KIND=8) :: v
      END FUNCTION template_psi
   END INTERFACE

   TYPE argument_limiting_functional_type
      CHARACTER(LEN=rec_length) :: char_relaxation_method   = '=== Relaxation method (none/avg/minmod) ==='
      CHARACTER(LEN=rec_length) :: relaxation_global_bound = '=== Global bounds for relaxation ==='
   END TYPE argument_limiting_functional_type

   TYPE ::  limiting_functional_type
      PROCEDURE(template_psi),         POINTER, NOPASS :: psi => NULL()
      PROCEDURE(template_zero_of_psi), POINTER, NOPASS :: zero_of_psi => NULL()
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
      INTEGER, DIMENSION(:,:), POINTER :: jj
      INTEGER                       :: nl
      CHARACTER(LEN=:), ALLOCATABLE :: name
      MPI_Comm, POINTER         :: comm
      Vec :: xvect1, xvect2, stiffness_RowSumAbs, ones_row_sum
      Mat :: stiffness, ones
   CONTAINS
      PROCEDURE, PUBLIC :: init => init_limiting_functionals
      PROCEDURE :: RELAX_BOUNDS
      PROCEDURE :: init_mat_relax_avg
      ! PROCEUDRE :: init_mat_relax_avg
   END TYPE limiting_all_functional_type
   
   INTEGER, PRIVATE, PARAMETER :: RELAX_NONE=1, RELAX_AVG=2, RELAX_MINMOD=3
   CHARACTER(LEN=20), DIMENSION(3), PRIVATE, PARAMETER  :: list_relax_method = &
               [CHARACTER(LEN=20) :: 'none', 'avg', 'minmod']

CONTAINS

!=====================================
!======= Limiting initialization =====
!=====================================

   SUBROUTINE init_limiting(this, communicator, name, mesh, LA)
#include "petsc/finclude/petsc.h"
      USE petsc 
      USE solver_petsc
      USE def_type_mesh
      USE compute_periodic
      USE fem_M
      USE fem_petsc_matrix_factory_module
      IMPLICIT NONE
      CLASS(limiting_type),    INTENT(INOUT) :: this
      MPI_Comm,                   INTENT(IN) :: communicator
      CHARACTER(100),             INTENT(IN) :: name
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
      CALL VecCreateGhost(communicator, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%xvect1, ierr)
      CALL VecGhostGetLocalForm(this%xvect1, this%x_ghost, ierr)
      CALL VecDuplicate(this%xvect1, this%xvect2, ierr)
      !===Compute lumped masse
      CALL create_local_petsc_matrix(communicator, LA, mass, clean = .FALSE.)
      ALLOCATE(this%lumped_mass(mesh%np))
      CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, mass)
      ! CALL periodic_matrix_petsc(mesh%per, LA, mass)
      CALL construct_lumped_mass_vector(mesh, LA, mass, this%xvect1)
      CALL periodic_add_vector_petsc(mesh%per%nb_bords, mesh%per%list, mesh%per%perlist, this%xvect1, LA)
      ! CALL periodic_vector_petsc(mesh%per%nb_bords, mesh%per%list, mesh%per%perlist, this%xvect1, LA)
      CALL extract_through_ghost(this%xvect1, 1, 1, this%LA, this%lumped_mass, opt_assemble=.TRUE.)
      ! CALL construct_lumped_mass(mesh, LA, mass, this%lumped_mass)
      ! DO k = 1, mesh%per%nb_bords
      !    this%lumped_mass(mesh%per%list(k)%DIL) = this%lumped_mass(mesh%per%perlist(k)%DIL)
      ! END DO
   
      !===Localized mass construction
      ALLOCATE(this%localized_mass(mesh%gauss%n_w,mesh%me))
      this%localized_mass = 0.d0

      ! IF (this%if_relax_bounds) THEN
      !    CALL create_local_petsc_matrix(communicator, LA, this%stiffness, clean = .FALSE.)
      !    CALL qs_mass_diff_M (mesh, 0.d0, 1.d0, LA, this%stiffness)
      !    ! CALL periodic_matrix_petsc(mesh%per, LA, this%stiffness) !<-- FIXME PERIODICITY
      !    CALL VecDuplicate(this%xvect1, this%stiffness_RowSumAbs, ierr)
      !    CALL MatGetRowSumAbs(this%stiffness, this%stiffness_RowSumAbs, ierr)
      !    CALL MatGetDiagonal(this%stiffness, this%xvect1, ierr)
      !    CALL VecAXPY(this%stiffness_RowSumAbs, -1.d0, this%xvect1, ierr)

      !    IF (TRIM(ADJUSTL(this%relaxation_method)) == 'avg') THEN
      !       CALL VecDuplicate(this%xvect1, this%ones_row_sum, ierr)
      !       CALL MatDuplicate(mass, MAT_DO_NOT_COPY_VALUES, this%ones, ierr)
      !       CALL init_mat_relax_avg(this)
      !    END IF   
      ! END IF

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
! write(*,*) 'vol of Ti', vol_of_Ti
! write(*,*) 'mass', this%lumped_mass
! stop
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
      INTEGER :: i, j, p, dom_np, ierr
      INTEGER, DIMENSION(1) :: idx, jdx
      REAL(KIND=8), DIMENSION(1,1) :: mat_loc

      CALL VecDuplicate(this%xvect1, this%ones_row_sum, ierr)
      CALL create_local_petsc_matrix(this%comm, this%LA, this%ones, clean = .FALSE.)

      dom_np = SIZE(this%LA%ia)-1
      DO i=0, dom_np-1
         idx = this%LA%loc_to_glob(1, i+1)-1
         DO p=this%LA%ia(i), this%LA%ia(i+1)-1
            j = this%LA%ja(p)
            jdx = j
            IF (idx(1) /= jdx(1)) THEN
               mat_loc = 1.d0
               CALL MatSetValues(this%ones, 1, idx, 1, jdx, mat_loc, INSERT_VALUES, ierr)
            ELSE
               mat_loc = 0.d0
               CALL MatSetValues(this%ones, 1, idx, 1, jdx, mat_loc, INSERT_VALUES, ierr)
            END IF
         END DO
      END DO

      CALL MatAssemblyBegin(this%ones, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(this%ones, MAT_FINAL_ASSEMBLY, ierr)
      CALL VecSet(this%xvect1, 1.d0, ierr)
      CALL MatMult(this%ones, this%xvect1, this%ones_row_sum, ierr)

   END SUBROUTINE init_mat_relax_avg

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
      this%comm => mesh%comm
      this%jj   => mesh%jj
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
      CALL create_my_ghost(mesh, LA, ifrom)
      CALL VecCreateGhost(mesh%comm, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%xvect1, ierr)
      CALL VecDuplicate(this%xvect1, this%xvect2, ierr)

      IF (ANY(RELAX_AVG==list_relax) .OR. ANY(RELAX_MINMOD==list_relax)) THEN
         CALL create_local_petsc_matrix(mesh%comm, LA, this%stiffness, clean = .FALSE.)
         CALL qs_mass_diff_M (mesh, 0.d0, 1.d0, LA, this%stiffness)
         ! CALL periodic_matrix_petsc(mesh%per, LA, this%stiffness) !<-- FIXME PERIODICITY
         CALL VecDuplicate(this%xvect1, this%stiffness_RowSumAbs, ierr)
         CALL MatGetRowSumAbs(this%stiffness, this%stiffness_RowSumAbs, ierr)
         CALL MatGetDiagonal(this%stiffness, this%xvect1, ierr)
         CALL VecAXPY(this%stiffness_RowSumAbs, -1.d0, this%xvect1, ierr)
      END IF

      IF (ANY(RELAX_AVG==list_relax)) THEN
         CALL this%init_mat_relax_avg
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

      ! !===if_relax_bounds
      ! CALL read_data(argument_data%if_relax_bounds, this%if_relax_bounds, &
      ! opt_name=this%name, opt_add=this%if_limiting)

      ! !===relaxation_method
      ! CALL read_data(argument_data%relaxation_method, this%relaxation_method, &
      ! opt_name=this%name, opt_add=(this%if_limiting .AND. this%if_relax_bounds))

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
      ! IF (PRESENT(section_name)) THEN
      !    CALL read_data_init_list(section_name)
      ! ELSE
         CALL read_data_init_list()
      ! END IF

      !================
      !=== We now find the relevant information for this specific limiting data
      !================

      ! !===if_relax_bounds
      ! CALL read_data(argument_data%if_relax_bounds, this%if_relax_bounds, &
      ! opt_name=this%name, opt_add=this%if_limiting)

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
      USE my_util, ONLY: error_petsc
      IMPLICIT NONE
      CLASS(limiting_type),             INTENT(IN) :: this
      CLASS(limiting_functional_type),  INTENT(IN) :: lim_bounds
      PROCEDURE(template_zero_of_psi), POINTER :: zero_of_psi
      PROCEDURE(template_psi)        , POINTER :: psi
      REAL(KIND=8), DIMENSION(:,:),                         INTENT(IN) :: xx_in
      REAL(KIND=8), DIMENSION(SIZE(xx_in,1),SIZE(xx_in,2)), INTENT(OUT):: xx_out

      REAL(KIND=8), DIMENSION(:)   :: loc_min
      REAL(KIND=8), DIMENSION(SIZE(loc_min)) :: loc_min_bis
      REAL(KIND=8), DIMENSION(SIZE(xx_in,2))               :: uk_minus, uk_plus
      REAL(KIND=8), DIMENSION(SIZE(xx_in,1))               :: psi_x
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1),SIZE(this%jj,2),SIZE(xx_in,2))    :: xx
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1),SIZE(xx_in,2))    :: xx_loc
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1)) :: lambda_minus, lambda_plus
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1)) :: loc_min_loc
      INTEGER,      DIMENSION(SIZE(this%jj,1)) :: jloc
      INTEGER,      DIMENSION(SIZE(this%jj,1)) :: limit_zero, limit_plus, limit_minus
      INTEGER :: k, m, n, me, nw, syst_size, iminus, iplus, comp, np, i
      REAL(KIND=8) :: loc_min_down, loc_min_up
      REAL(KIND=8) :: mass_plus, mass_minus, &
            lambda_K_minus, lambda_K_plus, &
            lambda_star_minus, lambda_star_plus

      np = SIZE(xx_in, 1)

      zero_of_psi => lim_bounds%zero_of_psi
      psi         => lim_bounds%psi

!       IF (this%if_relax_bounds) THEN
!          DO i=1, np
!             psi_x(i) = psi(xx_in(i, :), 0.d0)
!          END DO
! ! WRITE(*,*) "before", loc_min
!          ! loc_min_bis = loc_min
!          ! CALL this%RELAX_BOUNDS(psi_x, loc_min_bis)
!          CALL this%RELAX_BOUNDS(psi_x, loc_min)
! ! WRITE(*,*) "after", loc_min
! ! WRITE(*,*) "MAXVAL DIFF", SUM(ABS(loc_min_bis-loc_min))/SUM(ABS(loc_min))
!       END IF


      me = SIZE(this%jj,2)
      nw = SIZE(this%jj,1)
      syst_size = SIZE(xx_in,2)
      DO m = 1, me
         lambda_minus = 1.d0
         lambda_plus = 1.d0
         limit_zero = 0
         limit_minus = 0
         limit_plus = 0
         uk_minus = 0.d0
         uk_plus  = 0.d0
         mass_plus = 0.d0
         mass_minus = 0.d0
         jloc = this%jj(:,m)
         DO k = 1, syst_size
            xx_loc(:,k) = xx_in(jloc,k)
         END DO
         loc_min_loc = loc_min(jloc)
         iminus = 0
         iplus  = 0
         DO n = 1, nw
            !===P2 fix
            IF (ABS(this%lumped_mass(jloc(n))).LE.this%mass_eps) THEN
               limit_zero(n) = 1
               CYCLE
            END IF
            !===END fix

            loc_min_down = loc_min_loc(n) - this%epsilon*ABS(loc_min_loc(n))
            loc_min_up   = loc_min_loc(n) + this%epsilon*ABS(loc_min_loc(n))
            IF (psi(xx_loc(n,:),loc_min_down)<0) THEN
               iplus = iplus + 1
               uk_minus = uk_minus + this%localized_mass(n,m)*xx_loc(n,:)
               mass_minus = mass_minus + this%localized_mass(n,m)
               limit_minus(n) = 1
            ELSE IF (psi(xx_loc(n,:),loc_min_up)>0) THEN   
               iminus = iminus + 1
               uk_plus = uk_plus + this%localized_mass(n,m)*xx_loc(n,:)
               mass_plus = mass_plus + this%localized_mass(n,m)
               limit_plus(n) = 1
            ELSE
               limit_zero(n) = 1
            END IF
         END DO
         IF (SUM(limit_zero+limit_plus+limit_minus).NE.nw) THEN
            WRITE(*,*) ' BUG in iterative_cell_limiting_procedure:',&
                limit_zero,'+',limit_plus,'+',limit_minus,'.ne.', nw
            STOP
         END IF
         IF (iplus*iminus==0) THEN
            xx(:,m,:) = xx_loc
            CYCLE !===No limiting is possible/or no limiting necessary
         END IF
         mass_minus = mass_minus + 1.d-15 
         mass_plus = mass_plus + 1.d-15 
         uk_minus = uk_minus/mass_minus
         uk_plus  = uk_plus/mass_plus
         DO n = 1, nw
            !===Lambda_minus
            IF (limit_minus(n)==1) THEN
               lambda_minus(n) = zero_of_psi(loc_min_loc(n),uk_plus,xx_loc(n,:)-uk_plus)
            END IF
            !===Lambda_plus
            IF (limit_plus(n)==1) THEN
               lambda_plus(n) = zero_of_psi(loc_min_loc(n),xx_loc(n,:),uk_minus-xx_loc(n,:))
            END IF
         END DO
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
         DO n = 1, nw
            !!$ ===P2 fix
            IF (ABS(this%lumped_mass(jloc(n))).LE.this%mass_eps) THEN
               xx(n,m,:) = uk_plus(:)
            ELSE
            !!$ ===END fix
               xx(n,m,:) = xx_loc(n,:) &
                     +limit_minus(n)*(1-Lambda_K_minus)*(uk_plus(:)-xx_loc(n,:))&
                     +limit_plus(n) *     Lambda_K_plus*(uK_minus(:)-xx_loc(n,:))
            END IF
         END DO
      END DO

   !===Now we average over the nodes=========
      DO comp = 1, syst_size
         CALL this%cell_averaging(xx(:,:,comp), xx_out(:,comp))
      END DO     

   END SUBROUTINE iterative_cell_limiting_procedure

   SUBROUTINE cell_averaging(this,xx,xx_out)
#include "petsc/finclude/petsc.h"
      USE petsc 
      USE compute_periodic
      IMPLICIT NONE
      CLASS(limiting_type), INTENT(IN) :: this
      REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: xx
      ! REAL(KIND=8), DIMENSION(SIZE(this%jj,1),SIZE(this%jj,2)), INTENT(INOUT) :: xx
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT)  :: xx_out
      REAL(KIND=8), DIMENSION(SIZE(xx_out))    :: xx_inter
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1)) :: v_loc
      INTEGER, DIMENSION(SIZE(this%jj,1))      :: idxm
      INTEGER :: m, nw, me, ierr
      nw = SIZE(this%jj,1)
      me = SIZE(this%jj,2)
      xx_inter = 0.d0
      CALL VecZeroEntries(this%xvect1, ierr)
      DO m = 1, me
         WHERE(ABS(this%lumped_mass(this%jj(:,m))).GE.this%mass_eps)
            v_loc =  xx(:,m)*this%localized_mass(:,m)
         ELSEWHERE
            xx_out(this%jj(:,m)) = xx(:,m)
            v_loc = 0.d0
         END WHERE
         idxm = this%LA%loc_to_glob(1, this%jj(:,m)) -1
         CALL VecSetValues(this%xvect1, nw, idxm, v_loc, ADD_VALUES, ierr)
      END DO
      CALL VecAssemblyBegin(this%xvect1, ierr)
      CALL VecAssemblyEnd(this%xvect1, ierr)

      CALL periodic_add_vector_petsc(this%per%nb_bords, this%per%list, this%per%perlist, this%xvect1, this%LA)
      CALL extract_through_ghost(this%xvect1, 1, 1, this%LA, xx_inter, opt_assemble=.TRUE.)

      !===Rescaling
      WHERE (this%lumped_mass .GT.this%mass_eps)
         xx_out= xx_inter/this%lumped_mass
      END WHERE

   END SUBROUTINE cell_averaging

!    SUBROUTINE RELAX_BOUNDS(this, un, minn)
!       USE petsc_tools, ONLY: array_to_petsc_vec_bis
!       IMPLICIT NONE
!       CLASS(limiting_type),       INTENT(IN) :: this
!       REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: un
!       REAL(KIND=8), DIMENSION(:)             :: minn
! real(kind=8) :: norm
!       ! REAL(KIND=8), INTENT(IN)               :: glob_min, glob_max
!       REAL(KIND=8), DIMENSION(SIZE(un, 1)) :: alpha, denom
!       INTEGER,      DIMENSION(SIZE(un, 1)) ::   beta 
!       REAL(KIND=8), DIMENSION(SIZE(un, 1)) ::   psi_x 
!       INTEGER      :: i, j, ni, nj, m, me, nw, n, np, dom_np, ierr

!       np = SIZE(un)
!       DO i=1, np
!          psi_x(i) = psi(xx_in(i, :), 0.d0)
!       END DO

!       me = SIZE(this%jj,2)
!       nw = SIZE(this%jj,1)
!       alpha = 0.d0
!       beta = 0
!       ! DO m = 1, me
!       !    DO n = 1, nw
!       !       i = jj(n,m)
!       !       DO np = 1, nw
!       !          IF (n==np) CYCLE
!       !          j = jj(np,m)
!       !          alpha(i) = alpha(i) + (un(i) - un(j))
!       !          beta(i) = beta(i) + 1
!       !       END DO
!       !    END DO
!       ! END DO
!       ! alpha = alpha/beta
!       CALL array_to_petsc_vec_bis(un, this%xvect1, this%LA, "insert")
!       CALL MatMult(this%stiffness, this%xvect1, this%xvect2, ierr)
!       CALL VecPointWiseDivide(this%xvect2, this%xvect2, this%stiffness_RowSumAbs, ierr)
      

!       SELECT CASE(TRIM(ADJUSTL(this%relaxation_method)))
!       CASE('avg') !==Average
!          !denom = 0.d0
!          ! denom = alpha
!          ! beta = 0
!          ! DO m = 1, me
!          !    DO n = 1, nw
!          !       i = this%jj(n,m)
!          !       DO np = 1, nw
!          !          IF (n==np) CYCLE
!          !          j = this%jj(np,m) 
!          !          !denom(i) = denom(i) + alpha(i) + alpha(j)
!          !          denom(i) = denom(i) + alpha(j)
!          !          beta(i) = beta(i) + 1
!          !       END DO
!          !    END DO
!          ! END DO
!          ! !denom = denom/(2*beta)
!          ! denom = denom/(beta)
!          CALL MatMult(this%ones, this%xvect2, this%xvect1, ierr)
!          CALL VecPointWiseDivide(this%xvect1, this%xvect1, this%ones_row_sum, ierr)
!          CALL extract_through_ghost(this%xvect1, 1, 1, this%LA, denom, opt_assemble=.FALSE.)
         
!          !TESTTTT
!          ! CALL extract_through_ghost(this%ones_row_sum, 1, 1, this%LA, alpha, opt_assemble=.FALSE.)
!          ! dom_np = SIZE(this%LA%ia) - 1
!          ! norm = MAXVAL((alpha(1:dom_np) - (this%LA%ia(1:) - this%LA%ia(0:dom_np-1) - 1)))
!          ! WRITE(*,*) "max norm alpha = ", norm 
!          ! norm = MINVAL((alpha(1:dom_np) - (this%LA%ia(1:) - this%LA%ia(0:dom_np-1) - 1)))
!          ! WRITE(*,*) "min norm alpha = ", norm 
!          ! WRITE(*,*) "mean alpha = ", SUM(alpha)/SIZE(alpha)
!          ! WRITE(*,*) "dif ia", (this%LA%ia(1:) - this%LA%ia(0:dom_np-1))
!          ! stop
!          !TESTTTT
!       CASE('minmod') !===Minmod
!          denom = alpha 
!          DO m = 1, me
!             DO ni = 1, nw
!                i = this%jj(ni,m)
!                DO nj = 1, nw
!                   j = this%jj(nj,m)
!                   IF (denom(i)*alpha(j).LE.0.d0) THEN
!                      denom(i) = 0.d0
!                   ELSE IF (ABS(denom(i)) > ABS(alpha(j))) THEN
!                      denom(i) = alpha(j)
!                   END IF
!                END DO
!             END DO
!          END DO
!       CASE DEFAULT
!          WRITE(*,*) ' BUG in relax', TRIM(ADJUSTL(this%relaxation_method))
!          STOP
!       END SELECT
!       minn = minn - 4.*ABS(denom)
!       ! minn = MAX(glob_min,minn)
!    END SUBROUTINE RELAX_BOUNDS

   SUBROUTINE RELAX_BOUNDS(this, state, psi_min, nth_fonctional)
      USE petsc_tools, ONLY: array_to_petsc_vec_bis
      IMPLICIT NONE
      CLASS(limiting_all_functional_type),  INTENT(IN) :: this
      REAL(KIND=8), DIMENSION(:,:),         INTENT(IN) :: state
      INTEGER,                              INTENT(IN) :: nth_fonctional
      REAL(KIND=8), DIMENSION(:),       INTENT(INOUT)  :: psi_min
real(kind=8) :: norm
      ! REAL(KIND=8), INTENT(IN)               :: glob_min, glob_max
      REAL(KIND=8), DIMENSION(SIZE(state, 1)) :: alpha, denom
      INTEGER,      DIMENSION(SIZE(state, 1)) ::   beta 
      REAL(KIND=8), DIMENSION(SIZE(state, 1)) ::   un 
      INTEGER      :: i, j, ni, nj, m, me, nw, n, np, dom_np, ierr

      IF(this%limiting_functionals(nth_fonctional)%relaxation_method == RELAX_NONE) THEN
         psi_min = MAX(this%limiting_functionals(nth_fonctional)%relaxation_global_bound,psi_min)
         RETURN
      END IF

      np = SIZE(state, 1)
      DO i=1, np
         un(i) = this%limiting_functionals(nth_fonctional)%psi(state(i, :), 0.d0)
      END DO

      me = SIZE(this%jj,2)
      nw = SIZE(this%jj,1)
      alpha = 0.d0
      beta = 0
      ! DO m = 1, me
      !    DO n = 1, nw
      !       i = jj(n,m)
      !       DO np = 1, nw
      !          IF (n==np) CYCLE
      !          j = jj(np,m)
      !          alpha(i) = alpha(i) + (un(i) - un(j))
      !          beta(i) = beta(i) + 1
      !       END DO
      !    END DO
      ! END DO
      ! alpha = alpha/beta
      CALL array_to_petsc_vec_bis(un, this%xvect1, this%LA, "insert")
      CALL MatMult(this%stiffness, this%xvect1, this%xvect2, ierr)
      CALL VecPointWiseDivide(this%xvect2, this%xvect2, this%stiffness_RowSumAbs, ierr)
      

      SELECT CASE(this%limiting_functionals(nth_fonctional)%relaxation_method)
      CASE(RELAX_AVG) !==Average
         !denom = 0.d0
         ! denom = alpha
         ! beta = 0
         ! DO m = 1, me
         !    DO n = 1, nw
         !       i = this%jj(n,m)
         !       DO np = 1, nw
         !          IF (n==np) CYCLE
         !          j = this%jj(np,m) 
         !          !denom(i) = denom(i) + alpha(i) + alpha(j)
         !          denom(i) = denom(i) + alpha(j)
         !          beta(i) = beta(i) + 1
         !       END DO
         !    END DO
         ! END DO
         ! !denom = denom/(2*beta)
         ! denom = denom/(beta)
         CALL MatMult(this%ones, this%xvect2, this%xvect1, ierr)
         CALL VecPointWiseDivide(this%xvect1, this%xvect1, this%ones_row_sum, ierr)
         CALL extract_through_ghost(this%xvect1, 1, 1, this%LA, denom, opt_assemble=.FALSE.)
         
         !TESTTTT
         ! CALL extract_through_ghost(this%ones_row_sum, 1, 1, this%LA, alpha, opt_assemble=.FALSE.)
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
         denom = alpha 
         DO m = 1, me
            DO ni = 1, nw
               i = this%jj(ni,m)
               DO nj = 1, nw
                  j = this%jj(nj,m)
                  IF (denom(i)*alpha(j).LE.0.d0) THEN
                     denom(i) = 0.d0
                  ELSE IF (ABS(denom(i)) > ABS(alpha(j))) THEN
                     denom(i) = alpha(j)
                  END IF
               END DO
            END DO
         END DO
      CASE DEFAULT
         WRITE(*,*) ' BUG in relax, method not implemented: ', TRIM(ADJUSTL(this%limiting_functionals(nth_fonctional)%char_relaxation_method))
         STOP
      END SELECT
      psi_min = psi_min - 4.*ABS(denom)
      psi_min = MAX(this%limiting_functionals(nth_fonctional)%relaxation_global_bound,psi_min)
   END SUBROUTINE RELAX_BOUNDS


   ! SUBROUTINE relax_min_and_max(bound_relaxing,glob_min,glob_max,jj,un,maxn,minn)
   !    IMPLICIT NONE
   !    CHARACTER(*),               INTENT(IN) :: bound_relaxing
   !    INTEGER, DIMENSION(:,:),    INTENT(IN) :: jj
   !    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: un
   !    REAL(KIND=8), DIMENSION(:)             :: minn
   !    REAL(KIND=8), DIMENSION(:)             :: maxn
   !    REAL(KIND=8), INTENT(IN)               :: glob_min, glob_max
   !    REAL(KIND=8), DIMENSION(SIZE(un))      :: alpha, denom
   !    INTEGER, DIMENSION(SIZE(un)) ::   beta 
   !    INTEGER      :: i, j, m, me, nw, n, np

   !    me = SIZE(jj,2)
   !    nw = SIZE(jj,1)
   !    alpha = 0.d0
   !    beta = 0
   !    DO m = 1, me
   !       DO n = 1, nw
   !          i = jj(n,m)
   !          DO np = 1, nw
   !             IF (n==np) CYCLE
   !             j = jj(np,m)
   !             alpha(i) = alpha(i) + (un(i) - un(j))
   !             beta(i) = beta(i) + 1
   !          END DO
   !       END DO
   !    END DO
   !    alpha = alpha/beta
   !    SELECT CASE(TRIM(ADJUSTL(bound_relaxing)))
   !    CASE('avg') !==Average
   !       !denom = 0.d0
   !       denom = alpha
   !       beta = 0
   !       DO m = 1, me
   !          DO n = 1, nw
   !             i = jj(n,m)
   !             DO np = 1, nw
   !                IF (n==np) CYCLE
   !                j = jj(np,m) 
   !                !denom(i) = denom(i) + alpha(i) + alpha(j)
   !                denom(i) = denom(i) + alpha(j)
   !                beta(i) = beta(i) + 1
   !             END DO
   !          END DO
   !       END DO
   !       !denom = denom/(2*beta)
   !       denom = denom/(beta)
   !    CASE('minmod') !===Minmod
   !       denom = alpha 
   !       DO m = 1, me
   !          DO n = 1, nw
   !             i = jj(n,m)
   !             DO np = 1, nw
   !                j = jj(np,m)
   !                IF (denom(i)*alpha(j).LE.0.d0) THEN
   !                   denom(i) = 0.d0
   !                ELSE IF (ABS(denom(i)) > ABS(alpha(j))) THEN
   !                   denom(i) = alpha(j)
   !                END IF
   !             END DO
   !          END DO
   !       END DO
   !    CASE DEFAULT
   !       WRITE(*,*) ' BUG in relax', TRIM(ADJUSTL(bound_relaxing))
   !       STOP
   !    END SELECT
   !    maxn = maxn + 4.*ABS(denom)
   !    minn = minn - 4.*ABS(denom)
   !    maxn = MIN(glob_max,maxn)
   !    minn = MAX(glob_min,minn)
   ! END SUBROUTINE RELAX_MIN_AND_MAX

END MODULE cell_limiting_engine_parallel_module
