MODULE cell_limiting_engine_parallel_module
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE st_matrix, ONLY : extract_through_ghost, create_my_ghost
   USE read_inputs_module

   TYPE argument_limiting_type
      CHARACTER(LEN=rec_length) :: if_limiting       = '=== Apply cell-limiting (T/F) ? ==='
      CHARACTER(LEN=rec_length) :: if_relax_bounds   = '=== Apply bound relaxation for limiting (T/F) ? ==='
      CHARACTER(LEN=rec_length) :: relaxation_method = '=== Relaxation method (avg/minmod) ==='
   END TYPE argument_limiting_type

   TYPE limiting_type
      CHARACTER(100) :: name
      LOGICAL                   :: if_limiting       = .False.
      LOGICAL                   :: if_relax_bounds   = .False.
      CHARACTER(len=rec_length) :: relaxation_method ='minmod'
      INTEGER, DIMENSION(:,:), POINTER :: jj
      REAL(KIND=8) :: mass_eps
      REAL(KIND=8) :: epsilon = 1.d-8
      Vec, PRIVATE :: xvect, x_ghost
      REAL(KIND = 8), DIMENSION(:,:), POINTER :: localized_mass
      REAL(KIND = 8), DIMENSION(:),   POINTER :: lumped_mass
      TYPE(petsc_csr_LA),  POINTER :: LA
   CONTAINS
      PROCEDURE, PUBLIC  :: init => init_limiting
      PROCEDURE, PUBLIC  :: read => read_limiting_data
      PROCEDURE, PUBLIC  :: iterative_cell_limiting_procedure
   END TYPE limiting_type

   ABSTRACT INTERFACE
      FUNCTION template_zero_of_psi(psi_m,u0,P) RESULT(v)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u0, P
         REAL(KIND=8), INTENT(IN)  :: psi_m
         REAL(KIND=8), INTENT(OUT) :: v
      END FUNCTION template_zero_of_psi

      FUNCTION template_psi(x,psi_m) RESULT(v)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
         REAL(KIND=8), INTENT(IN) :: psi_m
         REAL(KIND=8), INTENT(OUT):: v
      END FUNCTION template_psi
   END INTERFACE

   TYPE :: limiting_bounds_type
   CONTAINS
      PROCEDURE(template_psi),         DEFERRED :: psi_min, psi_max
      PROCEDURE(template_zero_of_psi), DEFERRED :: zero_of_psi_min, zero_of_psi_max
   END TYPE limiting_bounds_type
   
CONTAINS

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
      CALL this%read("LIMITING PARAMETERS FOR "//TRIM(ADJUSTL(name)))
      this%name = 'limiting_for_'//TRIM(ADJUSTL(name))
      this%jj => mesh%jj
      this%LA => LA
    
    !===Petsc ghosting for cell-averaging
      CALL create_my_ghost(mesh, LA, ifrom)
      CALL VecCreateGhost(communicator, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%xvect, ierr)
      CALL VecGhostGetLocalForm(this%xvect, this%x_ghost, ierr)

      !===Compute lumped masse
      CALL create_local_petsc_matrix(communicator, LA, mass, clean = .FALSE.)
      ALLOCATE(this%lumped_mass(mesh%np))
      CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, mass)
      CALL periodic_matrix_petsc(mesh%per, LA, mass)
      CALL construct_lumped_mass(mesh, LA, mass, this%lumped_mass)
      DO k = 1, mesh%per%nb_bords
         this%lumped_mass(mesh%per%list(k)%DIL) = this%lumped_mass(mesh%per%perlist(k)%DIL)
      END DO
   
      !===Localized mass construction
      ALLOCATE(this%localized_mass(mesh%gauss%n_w,mesh%me))
      this%localized_mass = 0.d0

!VB CORRECTED VERSION WHEN SEVERAL PROCESSES
      vol_of_Ti_loc = 0.d0
      DO m = 1, mesh%me
         volK = SUM(mesh%gauss%rj(:,m))
         DO n = 1, mesh%gauss%n_w
            vol_of_Ti_loc(n) = volK
         END DO
         idxm = this%LA%loc_to_glob(1, this%jj(:,m)) -1
         CALL VecSetValues(this%xvect, mesh%gauss%n_w, idxm, vol_of_Ti_loc, ADD_VALUES, ierr)
      END DO

      CALL extract_through_ghost(this%xvect, this%x_ghost, 1, 1, this%LA, vol_of_Ti, &
                                'insert', opt_assemble=.TRUE.)
!VB CORRECTED VERSION WHEN SEVERAL PROCESSES

      DO m = 1, mesh%me
         volK = SUM(mesh%gauss%rj(:,m))
         DO n = 1, mesh%gauss%n_w
            this%localized_mass(n,m) = this%lumped_mass(mesh%jj(n,m))*volK/vol_of_Ti(mesh%jj(n,m))
         END DO
      END DO

      this%mass_eps = this%epsilon*SUM(this%lumped_mass)/mesh%np
   
   END SUBROUTINE init_limiting

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
      CALL read_data(argument_data%if_limiting , this%if_limiting)

      !===if_relax_bounds
      CALL read_data(argument_data%if_relax_bounds, this%if_relax_bounds)

      !===relaxation_method
      CALL read_data(argument_data%relaxation_method, this%relaxation_method)

      !================
      !=== MANDATORY to close data for the current section and
      !=== rewrite it with new information for the next sections
      !================
      CALL finalize_rewrite_data
   END SUBROUTINE read_limiting_data

   SUBROUTINE iterative_cell_limiting_procedure(this, xx_in, loc_min, lim_bounds, min_max, xx_out)  
   ! SUBROUTINE iterative_cell_limiting_procedure(this,xx_in,loc_min,psi,zero_of_psi,xx_out)  
      IMPLICIT NONE
      ! INTERFACE
      !    FUNCTION psi(x,psi_min) RESULT(v)
      !       REAL(KIND=8), DIMENSION(:) :: x
      !       REAL(KIND=8) :: psi_min, v
      !    END FUNCTION psi
      !    FUNCTION zero_of_psi(psi_min,u0,P) RESULT(v)
      !       REAL(KIND=8), DIMENSION(:) :: u0, P
      !       REAL(KIND=8) :: psi_min, v
      !    END FUNCTION zero_of_psi
      ! END INTERFACE
      CLASS(limiting_type),         INTENT(IN) :: this
      CLASS(limiting_bounds_type),  INTENT(IN) :: lim_bounds
      PROCEDURE(template_zero_of_psi), POINTER :: zero_of_psi
      PROCEDURE(template_psi)        , POINTER :: psi
      CHARACTER(LEN=*),             INTENT(IN) :: minmax
      REAL(KIND=8), DIMENSION(:,:),                         INTENT(IN) :: xx_in
      REAL(KIND=8), DIMENSION(SIZE(xx_in,1),SIZE(xx_in,2)), INTENT(OUT):: xx_out

      REAL(KIND=8), DIMENSION(:)   :: loc_min
      REAL(KIND=8), DIMENSION(SIZE(xx_in,2))               :: uk_minus, uk_plus
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1),SIZE(this%jj,2),SIZE(xx_in,2))    :: xx
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1),SIZE(xx_in,2))    :: xx_loc
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1)) :: lambda_minus, lambda_plus
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1)) :: loc_min_loc
      INTEGER,      DIMENSION(SIZE(this%jj,1)) :: jloc
      INTEGER,      DIMENSION(SIZE(this%jj,1)) :: limit_zero, limit_plus, limit_minus
      INTEGER :: k, m, n, me, nw, syst_size, iminus, iplus, comp
      REAL(KIND=8) :: loc_min_down, loc_min_up
      REAL(KIND=8) :: mass_plus, mass_minus, &
            lambda_K_minus, lambda_K_plus, &
            lambda_star_minus, lambda_star_plus

      SELECT CASE(minmax)
      CASE('MAX')
         zero_of_psi => lim_bounds%zero_of_psi_max
         psi         => lim_bounds%psi_max
      CASE('MIN')
         zero_of_psi => lim_bounds%zero_of_psi_min
         psi         => lim_bounds%psi_min
      CASE DEFAULT
         CALL error_petsc("BUG in iterative_cell_limiting_procedure: you selected "//minmax//&
         ", please select either MIN or MAX.")
      END SELECT

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
         CALL cell_averaging(this,xx(:,:,comp), xx_out(:,comp))
      END DO     

   END SUBROUTINE iterative_cell_limiting_procedure

   SUBROUTINE cell_averaging(this,xx,xx_out)
#include "petsc/finclude/petsc.h"
      USE petsc 
      IMPLICIT NONE
      CLASS(limiting_type), INTENT(INOUT) :: this
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1),SIZE(this%jj,2))    :: xx
      REAL(KIND=8), DIMENSION(:)               :: xx_out
      REAL(KIND=8), DIMENSION(SIZE(xx_out))    :: xx_inter
      REAL(KIND=8), DIMENSION(SIZE(this%jj,1)) :: v_loc
      INTEGER, DIMENSION(SIZE(this%jj,1))      :: idxm
      INTEGER :: m, n, i, nw, me, ierr
      nw = SIZE(this%jj,1)
      me = SIZE(this%jj,2)
      xx_inter = 0.d0
      CALL VecZeroEntries(this%xvect, ierr)
      DO m = 1, me
         WHERE(ABS(this%lumped_mass(this%jj(:,m))).GE.this%mass_eps)
            v_loc =  xx(:,m)*this%localized_mass(:,m)
         ELSEWHERE
            xx_out(this%jj(:,m)) = xx(:,m)
            v_loc = 0.d0
         END WHERE
         idxm = this%LA%loc_to_glob(1, this%jj(:,m)) -1
         CALL VecSetValues(this%xvect, nw, idxm, v_loc, ADD_VALUES, ierr)
      END DO

      CALL extract_through_ghost(this%xvect, this%x_ghost, 1, 1, this%LA, xx_inter, &
                                'insert', opt_assemble=.TRUE.)

      ! CALL VecAssemblyBegin(this%xvect, ierr)
      ! CALL VecAssemblyEnd(this%xvect, ierr)

      ! CALL VecGhostGetLocalForm(this%xvect, this%x_ghost, ierr)
      ! CALL VecGhostUpdateBegin(this%xvect, INSERT_VALUES, SCATTER_FORWARD, ierr)
      ! CALL VecGhostUpdateEnd  (this%xvect, INSERT_VALUES, SCATTER_FORWARD, ierr)
      ! CALL extract(this%x_ghost, 1, 1, this%LA, xx_inter)

      !===Rescaling
      WHERE (this%lumped_mass .GT.this%mass_eps)
         xx_out= xx_inter/this%lumped_mass
      END WHERE

   END SUBROUTINE cell_averaging

   SUBROUTINE relax_min_and_max(bound_relaxing,glob_min,glob_max,jj,un,maxn,minn)
      IMPLICIT NONE
      CHARACTER(*),               INTENT(IN) :: bound_relaxing
      INTEGER, DIMENSION(:,:),    INTENT(IN) :: jj
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: un
      REAL(KIND=8), DIMENSION(:)             :: minn
      REAL(KIND=8), DIMENSION(:)             :: maxn
      REAL(KIND=8), INTENT(IN)               :: glob_min, glob_max
      REAL(KIND=8), DIMENSION(SIZE(un))      :: alpha, denom
      INTEGER, DIMENSION(SIZE(un)) ::   beta 
      INTEGER      :: i, j, m, me, nw, n, np
      REAL(KIND=8) :: norm

      me = SIZE(jj,2)
      nw = SIZE(jj,1)
      alpha = 0.d0
      beta = 0
      DO m = 1, me
         DO n = 1, nw
            i = jj(n,m)
            DO np = 1, nw
               IF (n==np) CYCLE
               j = jj(np,m)
               alpha(i) = alpha(i) + (un(i) - un(j))
               beta(i) = beta(i) + 1
            END DO
         END DO
      END DO
      alpha = alpha/beta
      SELECT CASE(TRIM(ADJUSTL(bound_relaxing)))
      CASE('avg') !==Average
         !denom = 0.d0
         denom = alpha
         beta = 0
         DO m = 1, me
            DO n = 1, nw
               i = jj(n,m)
               DO np = 1, nw
                  IF (n==np) CYCLE
                  j = jj(np,m) 
                  !denom(i) = denom(i) + alpha(i) + alpha(j)
                  denom(i) = denom(i) + alpha(j)
                  beta(i) = beta(i) + 1
               END DO
            END DO
         END DO
         !denom = denom/(2*beta)
         denom = denom/(beta)
      CASE('minmod') !===Minmod
         denom = alpha 
         DO m = 1, me
            DO n = 1, nw
               i = jj(n,m)
               DO np = 1, nw
                  j = jj(np,m)
                  IF (denom(i)*alpha(j).LE.0.d0) THEN
                     denom(i) = 0.d0
                  ELSE IF (ABS(denom(i)) > ABS(alpha(j))) THEN
                     denom(i) = alpha(j)
                  END IF
               END DO
            END DO
         END DO
      CASE DEFAULT
         WRITE(*,*) ' BUG in relax', TRIM(ADJUSTL(bound_relaxing))
         STOP
      END SELECT
      maxn = maxn + 4.*ABS(denom)
      minn = minn - 4.*ABS(denom)
      maxn = MIN(glob_max,maxn)
      minn = MAX(glob_min,minn)
   END SUBROUTINE RELAX_MIN_AND_MAX

END MODULE cell_limiting_engine_parallel_module
