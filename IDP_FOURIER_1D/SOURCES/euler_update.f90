MODULE euler_module
  USE Butcher_tableau
  USE fourier_param_module
  USE cell_limiting_engine_module
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
     CHARACTER(LEN=rec_length) :: bound_relaxing = '=== Bound relaxing method (avg,minmod) ? ==='
  END TYPE argument_euler

  TYPE euler_type
     INTEGER :: syst_size=3, nb_bounds=2
     REAL(KIND = 8), DIMENSION(1) :: eos_param = 0.d0
     REAL(KIND=8) :: CFL = 0.5d0, glob_rho_min = 0.d0, glob_rho_max = 1.d20
     LOGICAL :: if_limiting = .FALSE.
     CHARACTER(LEN=40) :: bound_relaxing = 'minmod'
     CHARACTER(LEN=40) :: method = 'viscous'
     INTEGER      :: erk_sv    = -21
     TYPE(BT), PUBLIC :: ERK
     INTEGER :: Nmax, Nmax_real
     REAL(KIND=8) :: time, dt, lumped, dx, in_tol, final_time
     LOGICAL :: if_no_iter = .TRUE.
     REAL(KIND=8), DIMENSION(:,:,:), POINTER :: un
     REAL(KIND=8), DIMENSION(:,:), POINTER :: cij
     TYPE(fourier_param_type), POINTER :: FP
     !===limiting
     INTEGER, DIMENSION(:,:), POINTER :: jj
     TYPE(mass_for_limiting) :: mass
     !===end limiting
     PROCEDURE(pressure_template),   NOPASS, POINTER :: pressure
     PROCEDURE(flux_template),       NOPASS, POINTER :: flux
     PROCEDURE(lambda_max_template), NOPASS, POINTER :: lambda_max
   CONTAINS
     PROCEDURE, PUBLIC :: init => init_euler
     PROCEDURE, PUBLIC :: READ => read_euler
     PROCEDURE, PUBLIC :: update
  END TYPE euler_type

CONTAINS
  SUBROUTINE init_euler(this, flux, pressure, lambda_max, fourier_param, init_time, final_time)
    IMPLICIT NONE
    CLASS(euler_type), INTENT(INOUT) :: this
    PROCEDURE(flux_template) :: flux
    PROCEDURE(pressure_template) :: pressure
    PROCEDURE(lambda_max_template) :: lambda_max
    TYPE(fourier_param_type), TARGET :: fourier_param
    REAL(KIND=8) :: init_time, final_time
    INTEGER :: Nmax, m, n, i, j

    CALL this%READ()
    this%ERK%sv = this%erk_sv
    CALL this%ERK%init
    this%flux => flux
    this%pressure => pressure
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
    REAL(KIND=8), DIMENSION(this%Nmax_real,1) :: r_in, r_out
    INTEGER  :: stage, k, l, it
    !===Compute viscous flux, actual flux, and u_min, u_max (for limiting)
    CALL compute_dt_viscous_flux_min_max(this,stage,&
         urk(:,:,:,stage-1),urk(:,:,:,this%ERK%lp_of_l(stage)),&
         flux_rk(:,:,:,stage-1),cs_diff,bounds)

    !===Low-order
    IF (TRIM(ADJUSTL(this%method))=='viscous') THEN
       cs_dflux= this%ERK%inc_C(stage)*cs_diff
       urk(:,:,:,stage) = urk(:,:,:,this%ERK%lp_of_l(stage))+this%dt*cs_dflux
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
  
    !===Limiting on density
    IF (this%if_limiting) THEN
       CALL Fourier_to_real(urk(:,:,1,stage),r_in(:,1))
       DO it = 1, 1
          CALL iterative_cell_limiting_procedure(this%mass,this%jj,r_in,bounds(:,1),&
               psi_min,zero_of_psi_min,r_in)
          CALL iterative_cell_limiting_procedure(this%mass,this%jj,r_in,bounds(:,2),&
               psi_max,zero_of_psi_max,r_in)
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
    REAL(KIND = 8), DIMENSION(this%Nmax_real) :: r_diff, r_flux, rhomax, rhomin, alpha, eta
    REAL(KIND = 8), DIMENSION(this%Nmax_real,this%nb_bounds) :: bounds
    REAL(KIND = 8), DIMENSION(2) :: in_rho, in_u, in_e, in_p, lambda
    REAL(KIND = 8) :: x, y, ul, ur, cij, nij, uijbar, length, pstar
    INTEGER :: i, j, k, m, n, np

    DO k = 1, this%syst_size
       CALL Fourier_to_real(u_visc(:,:,k),r_out(:,k))
    END DO
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
       in_p = this%pressure(in_rho, in_e)
       eta(i) = in_p(1) !<==pressure for commutator
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
    CALL entropy_commutator(this,stage,eta,alpha)
    !CALL entropy_commutator(this,stage,r_out(:,1),alpha)
    
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
                !umax(i) = MAX(umax(i),r_out(j))
                !umin(i) = MIN(umin(i),r_out(j))
             END IF
          END DO
       END DO
       r_diff = r_diff/this%lumped
       CALL real_to_fourier(r_diff,cs_diff(:,:,k))
       CALL real_to_fourier(r_flux,cs_flux(:,:,k))
    END DO

    IF (this%if_limiting) THEN
       CALL relax_min_and_max(this%bound_relaxing,this%glob_rho_min,this%glob_rho_max,&
            this%jj,r_out(:,1),bounds(:,2),bounds(:,1))
    END IF
  END SUBROUTINE compute_dt_viscous_flux_min_max

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


    !TEST
    !r_eta = SIN(2*ACOS(-1.d0)*this%FP%rr)**2+1.01
    !TEST
    !FIX ME . Add padding to compute entropy commutator.
    
    !===Compute entropy commutator (eta*dx(log(eta(u)))-dx(eta(u)))
    length=this%dx*this%Nmax_real
    log_eta = LOG(ABS(r_eta))
    CALL real_to_fourier(log_eta,cs_diff)
    CALL fourier_derivative(cs_diff,cs_flux,length)
    CALL Fourier_to_real(cs_flux,alpha) !<==derivative of entropy
    alpha = r_eta*alpha !<==multiply by eta

    CALL real_to_fourier(r_eta,cs_flux)
    CALL fourier_derivative(cs_flux,cs_diff,length)
    CALL Fourier_to_real(cs_diff,beta) !<=derivative of flux

    avg = SUM(ABS(alpha))/this%Nmax_real
    alpha = min(abs(alpha - beta)/avg,1.d0)
    alpha = threshold(alpha)
    
    write(*,*) this%time+1.1*this%dt, this%final_time
    IF (this%time+1.1*this%dt>this%final_time .AND. stage==this%ERK%s+1) THEN
       write(*,*) 'commutator', stage, this%ERK%s
       CALL this%FP%plot_1d(alpha, 'commutator.plt')
    END IF

    
  END SUBROUTINE entropy_commutator

  FUNCTION threshold(x) RESULT(g)
    IMPLICIT NONE
    INTEGER, PARAMETER :: exp=3
    REAL(KIND=8), DIMENSION(:)  :: x
    REAL(KIND=8), DIMENSION(SIZE(x))  :: z, t, zp, relu, f, g
    REAL(KIND=8), PARAMETER :: x0 = 0.75d0, x1=SQRT(3.d0)*x0
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
END MODULE euler_module
