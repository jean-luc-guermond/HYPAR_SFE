MODULE euler_module
  USE Butcher_tableau
  USE fourier_param_module
  USE cell_limiting_engine_module
  USE eos_module
  INTEGER, PARAMETER :: rec_length = 200, list_length=200
  
  ABSTRACT INTERFACE
     FUNCTION pressure_template(rho,e) RESULT(vv)
       REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
       REAL(KIND = 8), DIMENSION(SIZE(rho,1)) :: vv
     END FUNCTION pressure_template
     FUNCTION flux_template(k,u) RESULT(vv)
       INTEGER, INTENT(IN)  :: k
       REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: u
       REAL(KIND = 8), DIMENSION(SIZE(u,1)) :: vv
     END FUNCTION  flux_template
     SUBROUTINE lambda_max_template(eos_param, in_rho, in_u, in_e, in_p, in_tol, no_iter, &
          lambda_max, pstar)                                                       
       IMPLICIT NONE                                                               
       REAL(KIND = 8), DIMENSION(:) :: eos_param !===b_covolume         
       REAL(KIND = 8), DIMENSION(2), INTENT(IN) :: in_rho, in_e, in_u, in_p        
       REAL(KIND = 8) :: in_tol                                                    
       LOGICAL, INTENT(IN) :: no_iter                                              
       REAL(KIND = 8), INTENT(OUT) :: pstar                                        
       REAL(KIND = 8), DIMENSION(2), INTENT(OUT) :: lambda_max
     END SUBROUTINE lambda_max_template
  END INTERFACE
  
  TYPE argument_euler
     CHARACTER(LEN=rec_length) :: CFL            = '=== CFL ? ==='
     CHARACTER(LEN=rec_length) :: erk_sv         = '=== ERK ? ==='
     CHARACTER(LEN=rec_length) :: method         = '=== Which method (viscous,high) ? ==='
     CHARACTER(LEN=rec_length) :: if_limiting    = '=== Limiting ? ==='
     CHARACTER(LEN=rec_length) :: glob_rho_min       = '=== Global rho min ? ==='
     CHARACTER(LEN=rec_length) :: glob_rho_max       = '=== Global rho max ? ==='
     CHARACTER(LEN=rec_length) :: bound_relaxing = '=== Bound relaxing method (avg,minmod,none) ? ==='
  END TYPE argument_euler

  TYPE euler_type
     INTEGER :: syst_size=3, nb_bounds=3
     REAL(KIND = 8), DIMENSION(1) :: eos_param = 0.d0
     REAL(KIND = 8) :: CFL = 0.5d0, glob_rho_min = 0.d0, glob_rho_max = 1.d20
     CHARACTER(LEN=40) :: bound_relaxing = 'none'
     CHARACTER(LEN=40) :: method = 'viscous'
     INTEGER :: erk_sv    = -21
     INTEGER :: Nmax, Nmax_real
     REAL(KIND = 8) :: time, dt, lumped, dx, in_tol, final_time
     REAL(KIND = 8), DIMENSION(:,:,:), POINTER :: un
     REAL(KIND = 8), DIMENSION(:,:),   POINTER :: cij
     !===limiting
     LOGICAL :: if_limiting = .FALSE. !<===To be read
     LOGICAL :: if_no_iter = .TRUE.   !<===To be read
     INTEGER :: it_limiting_max = 2
     INTEGER, DIMENSION(:,:), POINTER :: jj
     TYPE(mass_for_limiting) :: mass
     !===end limiting
     TYPE(BT),                  PUBLIC :: ERK
     TYPE(fourier_param_type), POINTER :: FP
     TYPE(eos_type),           POINTER :: eos
     PROCEDURE(pressure_template),   NOPASS, POINTER :: pressure
     PROCEDURE(pressure_template),   NOPASS, POINTER :: entropy
     PROCEDURE(lambda_max_template), NOPASS, POINTER :: lambda_max
   CONTAINS
     PROCEDURE, PUBLIC :: init => init_euler
     PROCEDURE, PUBLIC :: READ => read_euler
     PROCEDURE, PUBLIC :: update
     PROCEDURE, PRIVATE :: flux
  END TYPE euler_type

CONTAINS
  SUBROUTINE init_euler(this, eos, lambda_max, fourier_param, init_time, final_time)
    IMPLICIT NONE
    CLASS(euler_type), INTENT(INOUT) :: this
    PROCEDURE(flux_template) :: flux
    PROCEDURE(lambda_max_template) :: lambda_max
    TYPE(fourier_param_type), TARGET :: fourier_param
    TYPE(eos_type), TARGET :: eos
    REAL(KIND=8) :: init_time, final_time
    INTEGER :: Nmax, m, n, i, j

    CALL this%READ()
    this%ERK%sv = this%erk_sv
    CALL this%ERK%init
    this%lambda_max => lambda_max
    this%time = init_time
    this%final_time = final_time
    this%FP => fourier_param !<===FIXE ME: clear redundance with fourier_param

    !===Parameters for lambda_arbitrary_eos
    this%in_tol = 1.d-2
    this%if_no_iter = .true.
    
    Nmax = fourier_param%Nmax
    this%Nmax  = Nmax
    this%Nmax_real = 2*Nmax-1
    this%dx = fourier_param%dx
    this%lumped = fourier_param%dx
    ALLOCATE (this%un(Nmax,2,this%syst_size))

    !===limiting objects
    ALLOCATE(this%mass%localized_mass(2,this%Nmax_real))
    ALLOCATE(this%mass%lumped_mass(this%Nmax_real+1))
    this%mass%lumped_mass = this%dx
    this%mass%localized_mass=this%dx/2
    ALLOCATE(this%jj(2,this%Nmax_real))
    DO m = 1, this%Nmax_real !===loop over cells
       DO n = 1, 2
          this%jj(1,m) = m
          this%jj(2,m) = m+1
       END DO
    END DO
    this%jj(2,this%Nmax_real) = 1

    !===cij
    ALLOCATE(this%cij(this%Nmax_real,2))
    this%cij = 0.d0
    DO m = 1, this%Nmax_real !===loop over cells
       i = this%jj(1,m)
       j = this%jj(2,m)
       this%cij(i,2) = 0.5d0
       this%cij(j,1) = -0.5d0
    END DO
  END SUBROUTINE init_euler

  SUBROUTINE read_euler(this)
    USE character_strings
    IMPLICIT NONE
    CHARACTER(1), PARAMETER :: begin_section ='~'
    CHARACTER(1), PARAMETER :: end_section   ='~'
    INTEGER, PARAMETER :: in_unit = 21
    CLASS(euler_type), INTENT(INOUT):: this
    TYPE(argument_euler)  :: argument_data
    CHARACTER(LEN=rec_length), DIMENSION(list_length) :: list, record
    CHARACTER(LEN=rec_length)                         :: string_default
    LOGICAL :: okay
    INTEGER :: rank, record_size, i_list, j

    !===Initialize data to zero and false by default
    list = ''
    record = ''
    i_list = 1
    !===Initializing record
    CALL read_data_in_record(record_size, record, begin_section, end_section)

    !===CFL
    WRITE(string_default,*) this%CFL
    CALL compare_string(record, list, argument_data%CFL, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%CFL
    END IF

    !===ERK
    WRITE(string_default,*) this%erk_sv
    CALL compare_string(record, list, argument_data%erk_sv, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%erk_sv
    END IF

    !===Higher-order vs. low-order
    WRITE(string_default,*) this%method
    CALL compare_string(record, list, argument_data%method, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%method
    END IF

    !=========================
    !===Limiting parameters===
    !=========================
    !===if_limiting
    WRITE(string_default,*) this%if_limiting
    CALL compare_string(record, list, argument_data%if_limiting, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%if_limiting
    END IF

    !===Bound_relaxing
    WRITE(string_default,*) this%bound_relaxing
    CALL compare_string(record, list, argument_data%bound_relaxing, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%bound_relaxing
    END IF

    !===Closing unit
    rank = 0
    CALL rewrite_data_from_list_record(rank, list, record, i_list, record_size)
  END SUBROUTINE read_euler

  SUBROUTINE update(this, fourier_param)
    IMPLICIT NONE
    CLASS(euler_type), INTENT(INOUT):: this
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(This%Nmax,2,this%syst_size,this%ERK%s+1) :: urk
    REAL(KIND=8), DIMENSION(This%Nmax,2,this%syst_size,this%ERK%s)   :: flux_rk
    INTEGER  :: stage
    urk(:,:,:,1) = this%un
    DO stage = 2, this%ERK%s+1
       CALL one_step_ERK(this,stage,fourier_param,urk,flux_rk)
    END DO
    this%un = urk(:,:,:,this%ERK%s+1)
    this%time = this%time + this%dt
  END SUBROUTINE update

  SUBROUTINE one_step_ERK(this,stage,fourier_param,urk,flux_rk)
    USE fft_1D
    IMPLICIT NONE
    CLASS(euler_type), INTENT(INOUT):: this
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(this%Nmax,2,this%syst_size,this%ERK%s+1) :: urk
    REAL(KIND=8), DIMENSION(This%Nmax,2,this%syst_size,this%ERK%s)   :: flux_rk
    REAL(KIND=8), DIMENSION(this%Nmax,2,this%syst_size) :: cs_diff, cs_dflux, cs_zz
    REAL(KIND=8), DIMENSION(this%Nmax_real,this%nb_bounds) :: bounds
    REAL(KIND=8), DIMENSION(this%Nmax_real,this%syst_size) :: r_in, r_out
    REAL(KIND = 8) :: test
    INTEGER  :: stage, k, l, it
    !===Compute viscous flux, actual flux, and u_min, u_max (for limiting)
    CALL compute_dt_viscous_flux_min_max(this,stage,&
         urk(:,:,:,stage-1),urk(:,:,:,this%ERK%lp_of_l(stage)),&
         flux_rk(:,:,:,stage-1),cs_diff,bounds)

    !===Low-order
    IF (TRIM(ADJUSTL(this%method))=='viscous') THEN
       cs_dflux= this%ERK%inc_C(stage)*cs_diff
       urk(:,:,:,stage) = urk(:,:,:,this%ERK%lp_of_l(stage))+this%dt*cs_dflux

       !TESTTTT
       DO k = 1, this%syst_size
          CALL Fourier_to_real(urk(:,:,k,stage),r_in(:,k))
       END DO
       test=1.d20
       Do k =1, this%Nmax_real
          test = min(test,entropy_min(r_in(k,:),bounds(k,3)))
       END DO
       write(*,*) ' NEW TIME STEP', test
!!$       if (test<-1.d-7) THEN
!!$          write(*,*) 'complex vel', MAXVAL(ABS(urk(:,:,2,stage)))
!!$          Do k =1, this%Nmax_real
!!$             write(*,*) r_in(k,:)
!!$             write(*,*) entropy_min(r_in(k,:),bounds(k,3))
!!$          END DO
!!$          stop
!!$       END if
       !TESTTTT
       RETURN
    END IF

    !===Higher-order
    IF (TRIM(ADJUSTL(this%method)).NE.'high') THEN
       WRITE(*,*) 'Bug in one_step_ERK, wrong method', this%method
       STOP
    END IF

    cs_zz =0.d0
    DO l = 1, stage-1
       cs_zz =  cs_zz + this%ERK%MatRK(stage,l)*flux_rk(:,:,:,l)
    END DO
    DO k = 1, this%syst_size
       CALL fourier_derivative(cs_zz(:,:,k),cs_dflux(:,:,k),fourier_param%Length)
    END DO


    IF (this%erk_sv>0) THEN
       cs_dflux= -cs_dflux + this%ERK%C(stage)*cs_diff
       urk(:,:,:,stage) = urk(:,:,:,1)+this%dt*cs_dflux
    ELSE
       cs_dflux= -cs_dflux + this%ERK%inc_C(stage)*cs_diff
       urk(:,:,:,stage) = urk(:,:,:,this%ERK%lp_of_l(stage))+this%dt*cs_dflux
    END IF


    !TESTTTT
    DO k = 1, this%syst_size
       CALL Fourier_to_real(urk(:,:,k,stage),r_in(:,k))
    END DO
    test=1.d20
    Do k =1, this%Nmax_real
       test = min(test,entropy_min(r_in(k,:),bounds(k,3)))
    END DO
    write(*,*) ' NEW TIME STEP', test
!!$    if (test<-1.d-7) THEN
!!$       write(*,*) 'complex vel', MAXVAL(ABS(urk(:,:,2,stage)))
!!$       Do k =1, this%Nmax_real
!!$          write(*,*) r_in(k,:)
!!$          write(*,*) entropy_min(r_in(k,:),bounds(k,3))
!!$       END DO
!!$       stop
!!$    END if
    !TESTTTT

    
    !===Limiting on density
    IF (this%if_limiting) THEN
       DO k = 1, this%syst_size
          CALL Fourier_to_real(urk(:,:,k,stage),r_in(:,k))
       END DO
       DO it = 1, this%it_limiting_max      
          !CALL iterative_cell_limiting_procedure(this%mass,this%jj,r_in(:,1),bounds(:,1),&
          !     psi_min,zero_of_psi_min,r_in)
       END DO
       DO it = 1, this%it_limiting_max
          !CALL iterative_cell_limiting_procedure(this%mass,this%jj,r_in(:1,),bounds(:,2),&
          !     psi_max,zero_of_psi_max,r_in)
          !write(*,*) 'mass after', sum(this%lumped*r_in)

          r_out = r_in
          CALL iterative_cell_limiting_procedure(this%mass,this%jj,r_in,bounds(:,3),&
               entropy_min,zero_of_specific_entropy,r_in)
          write(*,*) 'AFTER'
          Do k =1, this%Nmax_real
             write(*,*) entropy_min(r_in(k,:),bounds(k,3)), entropy_min(r_out(k,:),bounds(k,3))
          END DO
       END DO
       CALL real_to_fourier(r_in(:,1),urk(:,:,1,stage))
    END IF
  END SUBROUTINE one_step_ERK
  
  SUBROUTINE compute_dt_viscous_flux_min_max(this,stage,urk_in,u_visc,&
       cs_flux,cs_diff,bounds)
    USE fft_1D
    IMPLICIT NONE
    CLASS(euler_type), INTENT(INOUT):: this
    INTEGER, INTENT(IN) :: stage
    REAL(KIND = 8), DIMENSION(this%Nmax,2,this%syst_size), INTENT(IN)  :: urk_in, u_visc
    REAL(KIND = 8), DIMENSION(this%Nmax,2,this%syst_size), INTENT(OUT) :: cs_diff, cs_flux
    REAL(KIND = 8), DIMENSION(this%Nmax_real,this%syst_size)  :: r_out
    REAL(KIND = 8), DIMENSION(this%Nmax_real,2) :: dijL
    REAL(KIND = 8), DIMENSION(this%Nmax_real)   :: diag_dijL
    REAL(KIND = 8), DIMENSION(this%Nmax_real) :: r_diff, r_flux, rhomax, rhomin, alpha, eta, rho, e
    REAL(KIND = 8), DIMENSION(this%Nmax_real,this%nb_bounds) :: bounds
    REAL(KIND = 8), DIMENSION(2) :: in_rho, in_u, in_e, in_p, lambda
    REAL(KIND = 8) :: x, y, ul, ur, cij, nij, uijbar, length, pstar
    INTEGER :: i, j, k, m, n, np

    bounds(:,3) = 1.d20
    DO k = 1, this%syst_size
       CALL Fourier_to_real(u_visc(:,:,k),r_out(:,k))
    END DO
    rho = r_out(:,1)
    e = r_out(:,3)/rho - 0.5d0*(r_out(:,2)/rho)**2
    
    cij =0.5d0
    diag_dijL=0.d0
    DO m = 1, this%Nmax_real !===loop over cells
       i = this%jj(1,m)
       j = this%jj(2,m)
       in_rho(1) = r_out(i,1)
       in_rho(2) = r_out(j,1)
       nij = this%cij(1,2)/abs(this%cij(1,2))
       in_u(1) = r_out(i,2)*nij/in_rho(1)
       in_u(2) = r_out(j,2)*nij/in_rho(2)
       in_e(1) = r_out(i,this%syst_size)/in_rho(1) - 0.5d0*in_u(1)*in_u(1)
       in_e(2) = r_out(j,this%syst_size)/in_rho(2) - 0.5d0*in_u(2)*in_u(2)
       in_p = this%eos%pressure(in_rho, in_e)       
       bounds(i,3) = MIN(bounds(i,3),MINVAL(this%eos%entropy(in_rho,in_e))) !<===entropy lower bound
       CALL this%lambda_max(this%eos_param, in_rho, in_u, in_e, in_p, &
            this%in_tol, this%if_no_iter, lambda, pstar)
       dijL(i,2) = cij*MAXVAL(lambda)
       dijL(j,1) = cij*MAXVAL(lambda)
       diag_dijL(i) = diag_dijL(i) - dijL(i,2)
       diag_dijL(j) = diag_dijL(j) - dijL(j,1)
    END DO

    IF (stage==2) THEN !<==Compute time step
       this%dt = this%ERK%s*0.5d0*this%CFL*this%lumped/MAXVAL(ABS(diag_dijL))
       IF (this%time+this%dt>this%final_time) THEN
          this%dt = this%final_time-this%time
       END IF
    END IF

    !===Only low-order
    IF (TRIM(ADJUSTL(this%method))=='viscous') THEN
       DO k = 1, this%syst_size
          r_diff = 0.d0
          r_flux = this%flux(k,r_out)
          DO m = 1, this%Nmax_real !===loop over cells
             DO n = 1, 2
                i = this%jj(n,m)
                np = MOD(n,2)+1
                j = this%jj(np,m)
                r_diff(i) = r_diff(i) + dijL(i,np)*(r_out(j,k)-r_out(i,k))
                r_diff(i) = r_diff(i) - this%cij(n,np)*(r_flux(j) - r_flux(i))
             END DO
          END DO
          r_diff = r_diff/this%lumped
          CALL real_to_fourier(r_diff,cs_diff(:,:,k))
       END DO
       RETURN
    END IF

    !===High-order
    !eta=this%eos%entropy(rho,e)
    eta=this%eos%pressure(rho,e)
    !eta = rho
    CALL entropy_commutator(this,stage,eta,alpha)
    !CALL entropy_commutator_pad(this,stage,2*this%Nmax,u_visc,eta,alpha)
    
    IF (this%ERK%lp_of_l(stage).NE.stage-1) THEN
       DO k = 1, this%syst_size
          CALL Fourier_to_real(urk_in(:,:,k),r_out(:,k))
       END DO
    END IF
    bounds(:,1) = r_out(:,1)
    bounds(:,2) = r_out(:,1)
    DO k = 1, this%syst_size
       r_diff = 0.d0
       r_flux = this%flux(k,r_out)
       DO m = 1, this%Nmax_real !===loop over cells
          DO n = 1, 2
             i = this%jj(n,m)
             np = MOD(n,2)+1
             j = this%jj(np,m)
             r_diff(i) = r_diff(i) + 0.5*(alpha(i)+alpha(j))*dijL(i,np)*(r_out(j,k)-r_out(i,k))
             IF (k==1) THEN !<== min and max on density
                uijbar = 0.5d0*((r_out(j,k)+r_out(i,k)) &
                     - (r_flux(j) - r_flux(i))*this%cij(n,np)/dijL(i,np))
                bounds(i,1) = MIN(bounds(i,1),uijbar)
                bounds(i,2) = MAX(bounds(i,2),uijbar)
             END IF
          END DO
       END DO
       r_diff = r_diff/this%lumped
       CALL real_to_fourier(r_diff,cs_diff(:,:,k))
       CALL real_to_fourier(r_flux,cs_flux(:,:,k))
    END DO

    IF (this%if_limiting .AND. TRIM(ADJUSTL(this%bound_relaxing)).NE.'none') THEN
       CALL relax_min_and_max(this%bound_relaxing,this%glob_rho_min,this%glob_rho_max,&
            this%jj,r_out(:,1),bounds(:,2),bounds(:,1))
       CALL relax_cmin(this,r_out,bounds(:,3))
    END IF
  END SUBROUTINE compute_dt_viscous_flux_min_max

  FUNCTION flux(this, comp, un) RESULT(vv)
    IMPLICIT NONE
    CLASS(euler_type), INTENT(INOUT):: this
    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
    INTEGER, INTENT(IN) :: comp
    REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: vv
    REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: u, ie
    SELECT CASE(comp)
    CASE(1)
       vv(:) = un(:,2)
    CASE(2)
       u = un(:,2)/un(:,1)
       vv(:) = un(:,2)*u
       ie = un(:,3)/un(:,1)
       ie = ie - 0.5d0*u**2
       vv(:) = vv(:) + this%eos%pressure(un(:,1), ie)
    CASE(3)
       ie = un(:,3)/un(:,1)
       ie = ie - 0.5d0 *(un(:,2)/un(:,1))**2
       vv(:) = (un(:,2)/un(:,1)) * (un(:,3) + this%eos%pressure(un(:,1), ie))
    CASE DEFAULT
       WRITE(*, *) ' BUG in flux'
       STOP
    END SELECT
  END FUNCTION flux
  
  FUNCTION psi_min(x,psi_m) RESULT(v)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:) :: x
    REAL(KIND=8) :: psi_m, v
    v = x(1)-psi_m
  END FUNCTION psi_min

  FUNCTION zero_of_psi_min(psi_m,u0,P) RESULT(v)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:) :: u0, P
    REAL(KIND=8) :: psi_m, v
    v = (psi_m-u0(1))/P(1)
  END FUNCTION zero_of_psi_min

  FUNCTION psi_max(x,psi_m) RESULT(v)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:) :: x
    REAL(KIND=8) :: psi_m, v
    v = psi_m-x(1)
  END FUNCTION psi_max

  FUNCTION zero_of_psi_max(psi_m,u0,P) RESULT(v)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:) :: u0, P
    REAL(KIND=8) :: psi_m, v
    v = (psi_m-u0(1))/P(1)
  END FUNCTION zero_of_psi_max

 FUNCTION entropy_min(u,cmin) RESULT(psi)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:) :: u
    REAL(KIND=8) :: cmin, psi
    psi = u(3) - u(2)**2/(2.d0*u(1)) - cmin*u(1)**gamma
  END FUNCTION entropy_min
  
  FUNCTION zero_of_specific_entropy(cmin,u0,Pij) RESULT(v)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:) :: u0, Pij
    REAL(KIND=8)   :: cmin, v
    REAL(KIND=8), DIMENSION(SIZE(u0)) :: ul, ur
    REAL(KIND=8), PARAMETER :: small=1.d-8
    REAL(KIND=8) :: Esmall, psir, psil, ll, lr, llold, lrold

    Esmall= small*u0(3)
    ur = u0 + Pij 

    psir = entropy_min(ur,cmin)
    IF (psir.GE.-Esmall) THEN
       write(*,*) ' OK', psir, -Esmall
       v = 1.d0
       RETURN
    END IF
    ll = 0.d0
    ul = u0
    psil = entropy_min(ul,cmin)
    write(*,*) 'psil', psil, psir, ABS(psil-psir), Esmall
    DO WHILE (ABS(psil-psir) .GT. Esmall)
       llold = ll
       lrold = lr
       ll = ll - psil*(lr-ll)/(psir-psil)
       lr = lr - psir/psi_prime_func(Pij,ur,cmin)
       IF (ll.GE.lr) THEN
          ll = lr
          EXIT
       END IF
       IF (ll< llold) THEN
          ll = llold
          EXIT
       END IF
       IF (lr > lrold) THEN
          lr = lrold
          EXIT
       END IF
       ul = u0 + ll*Pij
       ur = u0 + lr*Pij
       psil = entropy_min(ul,cmin)
       psir = entropy_min(ur,cmin)
    END DO
    IF (psir.GE.-Esmall) THEN
       v = lr
    ELSE
       v = ll
    END IF
  END FUNCTION zero_of_specific_entropy
  
  FUNCTION psi_prime_func(Pij,u,cmin) RESULT(psi)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:) :: u, Pij
    REAL(KIND=8)               :: cmin, psi
    psi = Pij(3) - u(2)*Pij(2)/u(1) + Pij(1)*u(2)**2/(2*u(1)**2) &
         - cmin*gamma*Pij(1)*u(1)**(gamma-1)
  END FUNCTION psi_prime_func

  SUBROUTINE entropy_commutator(this,stage,r_eta,alpha)
    USE fft_1D
    IMPLICIT NONE
    CLASS(euler_type), INTENT(INOUT):: this
    INTEGER :: stage
    REAL(KIND = 8), DIMENSION(this%Nmax_real):: r_eta
    REAL(KIND = 8), DIMENSION(this%Nmax,2)   :: cs_diff, cs_flux
    REAL(KIND = 8), DIMENSION(this%Nmax_real):: alpha, beta, log_eta
    REAL(KIND = 8) :: Length, x, avg
    INTEGER :: i

    !DO i = 1, this%Nmax_real
    !   IF (this%FP%rr(i) <0.5) THEN
    !      r_eta(i) = 1.d0
    !   ELSE
    !      r_eta(i) = .1d0
    !   END IF
    !END DO
    !===Padding does not improve the computation of the commutator
    !===Compute entropy commutator (eta*dx(log(eta(u)))-dx(eta(u)))
    length = this%dx*this%Nmax_real
    log_eta = LOG(ABS(r_eta))
    CALL real_to_fourier(log_eta,cs_diff)
    CALL fourier_derivative(cs_diff,cs_flux,length)
    CALL Fourier_to_real(cs_flux,alpha) !<==derivative of entropy
    alpha = r_eta*alpha !<==multiply by eta
    CALL real_to_fourier(r_eta,cs_flux)
    CALL fourier_derivative(cs_flux,cs_diff,length)
    CALL Fourier_to_real(cs_diff,beta) !<=derivative of flux
    !avg = SUM(ABS(alpha))/this%Nmax_real
    avg= MAXVAL(ABS(alpha))
    alpha = min(abs(alpha - beta)/avg,1.d0)
    !alpha = threshold(alpha) 
    IF (this%time+1.1*this%dt>this%final_time .AND. stage==this%ERK%s+1) THEN
       CALL this%FP%plot_1d(alpha, 'commutator.plt')
    END IF
    !TESTTTTTTTTTTTTT
    alpha =1.
  END SUBROUTINE entropy_commutator

  SUBROUTINE entropy_commutator_pad(this,stage,Nmax_pad,cs_u,r_eta,alpha)
    USE fft_1D
    IMPLICIT NONE
    CLASS(euler_type), INTENT(INOUT):: this
    INTEGER :: stage
    INTEGER :: Nmax_pad
    REAL(KIND = 8), DIMENSION(this%Nmax,2,3)   :: cs_u
    REAL(KIND = 8), DIMENSION(2*Nmax_pad-1,3) :: ru_pad
    REAL(KIND = 8), DIMENSION(this%Nmax_real):: r_eta
    REAL(KIND = 8), DIMENSION(this%Nmax,2)   :: cs1
    REAL(KIND = 8), DIMENSION(Nmax_pad,2)   :: cs1_pad, cs2_pad
    REAL(KIND = 8), DIMENSION(2*Nmax_pad-1) :: r_eta_pad, log_eta_pad, alpha_pad, beta_pad
    REAL(KIND = 8), DIMENSION(this%Nmax_real):: alpha
    REAL(KIND = 8) :: Length, x, avg
    INTEGER :: i, k
    length=this%dx*this%Nmax_real
    !r_eta = SIN(2*ACOS(-1.d0)*20*this%FP%rr)**2+1.01
    !CALL real_to_fourier(r_eta,cs1)
    !cs1_pad = 0.d0
    !cs1_pad(1:this%Nmax,:) = cs1
    !CALL fourier_to_real_padded(cs1_pad,r_eta_pad,Nmax_pad)
    cs1_pad = 0.d0
    DO k = 1, 3
       cs1_pad(1:this%Nmax,:) = cs_u(:,:,k)
       CALL fourier_to_real_padded(cs1_pad,ru_pad(:,k),Nmax_pad)
    END DO
    !r_eta_pad = ru_pad(:,1) 
    r_eta_pad = ru_pad(:,3) - 0.5d0*ru_pad(:,2)**2/ru_pad(:,1)
    !r_eta_pad = (ru_pad(:,3) - 0.5d0*ru_pad(:,2)**2/ru_pad(:,1))/ru_pad(:,1)**1.4
    !===Compute entropy commutator (eta*dx(log(eta(u)))-dx(eta(u)))
    log_eta_pad = LOG(ABS(r_eta_pad))
    CALL real_to_fourier(log_eta_pad,cs1_pad)
    CALL fourier_derivative(cs1_pad,cs2_pad,length)
    CALL Fourier_to_real(cs2_pad,alpha_pad) !<==derivative of entropy
    alpha_pad = r_eta_pad*alpha_pad !<==multiply by eta
    CALL real_to_fourier(r_eta_pad,cs1_pad)
    CALL fourier_derivative(cs1_pad,cs2_pad,length)
    CALL Fourier_to_real(cs2_pad,beta_pad) !<=derivative of flux
    avg = SUM(ABS(alpha_pad))/size(alpha_pad)
    alpha_pad = min(abs(alpha_pad - beta_pad)/avg,1.d0)
    !alpha = threshold(alpha)
    IF (this%time+1.1*this%dt>this%final_time .AND. stage==this%ERK%s+1) THEN
       DO i = 1, 2*Nmax_pad-1
          write(10,*) (i-1)*this%FP%length/(2*Nmax_pad-1), alpha_pad(i)
          write(11,*) (i-1)*this%FP%length/(2*Nmax_pad-1), r_eta_pad(i)
       END DO
    END IF
    !alpha = alpha_pad
    alpha = alpha_pad(1::2)
    !alpha=1
  END SUBROUTINE entropy_commutator_pad

  FUNCTION threshold(x) RESULT(g)
    IMPLICIT NONE
    INTEGER, PARAMETER :: exp=3
    REAL(KIND=8), DIMENSION(:)  :: x
    REAL(KIND=8), DIMENSION(SIZE(x))  :: z, t, zp, relu, f, g
    REAL(KIND=8), PARAMETER :: x0 = 0.25d0, x1=SQRT(3.d0)*x0
    SELECT CASE(exp)
    CASE(2)
       !===Quadratic threshold    
       z = x-x0
       zp = x-2*x0
       relu = (zp+abs(zp))/2
       f = -z*(z**2-x1**2)  + relu*(z-x0)*(z+2*x0)
       g = (f + 2*x0**3)/(4*x0**3)
    CASE(3)
       !===Cubic threshold
       relu = ((x-2*x0)+abs(x-2*x0))/2
       t = x/(2*x0)
       g = t**3*(10-15*t+6*t**2) - relu*(t-1)**2*(6*t**2+3*t+1)/(2*x0)
    END SELECT
    RETURN
  END FUNCTION threshold

  SUBROUTINE relax_cmin(this,un,cmin)
    IMPLICIT NONE
    CLASS(euler_type), INTENT(INOUT):: this
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(:)             :: cmin
    REAL(KIND=8), DIMENSION(SIZE(cmin))    :: dc
    REAL(KIND=8), DIMENSION(SIZE(cmin))    :: ul
    INTEGER      :: i, j, m, n, np, me, nw
    REAL(KIND=8), DIMENSION(1) :: e, rho, cl
    dc = 0.d0
    me = SIZE(this%jj,2)
    nw = SIZE(this%jj,1)
    DO m = 1, me
       DO n = 1, nw
          i = this%jj(n,m)
          DO np = 1, nw
             j = this%jj(np,m)
             IF (i==j) CYCLE
             ul(:) = (un(i,:)+un(j,:))/2
             rho(1) = ul(1)
             e(1) = ul(3)/ul(1) - 0.5d0*(ul(2)/ul(1))**2
             cl = this%eos%entropy(rho, e)
             dc(i) = MAX(dc(i),cl(1)-cmin(i))
          END DO
       END DO
    END DO
    cmin = cmin - dc
  END SUBROUTINE RELAX_cmin
END MODULE euler_module
