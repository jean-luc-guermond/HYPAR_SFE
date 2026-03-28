MODULE nl_scalar_cons_module
  USE Butcher_tableau
  USE fourier_param_module
  USE cell_limiting_engine_module
  INTEGER, PARAMETER :: rec_length = 200, list_length=200
  ABSTRACT INTERFACE
     FUNCTION flux_template(u) RESULT(vv)
       REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: u
       REAL(KIND = 8), DIMENSION(SIZE(u)) :: vv
     END FUNCTION  flux_template
     FUNCTION lambda_max_template(ul,ur) RESULT(vv)
       REAL(KIND = 8), INTENT(IN) :: ul, ur
       REAL(KIND = 8) :: vv
     END FUNCTION lambda_max_template
  END INTERFACE

  TYPE argument_nl_scalar_cons
     CHARACTER(LEN=rec_length) :: CFL            = '=== CFL ? ==='
     CHARACTER(LEN=rec_length) :: erk_sv         = '=== ERK ? ==='
     CHARACTER(LEN=rec_length) :: method         = '=== Which method (viscous,high) ? ==='
     CHARACTER(LEN=rec_length) :: if_limiting    = '=== Limiting ? ==='
     CHARACTER(LEN=rec_length) :: glob_min       = '=== Global min ? ==='
     CHARACTER(LEN=rec_length) :: glob_max       = '=== Global max ? ==='
     CHARACTER(LEN=rec_length) :: bound_relaxing = '=== Bound relaxing method (avg,minmod) ? ==='
  END TYPE argument_nl_scalar_cons

  TYPE nl_scalar_cons_type
     REAL(KIND=8) :: CFL = 0.5d0, glob_min = -1.d20, glob_max = 1.d20
     LOGICAL :: if_limiting = .FALSE.
     CHARACTER(LEN=40) :: bound_relaxing = 'minmod'
     CHARACTER(LEN=40) :: method = 'viscous'
     INTEGER      :: erk_sv    = -21
     TYPE(BT), PUBLIC :: ERK
     INTEGER :: Nmax, Nmax_real
     REAL(KIND=8) :: time, dt, lumped, dx
     REAL(KIND=8), DIMENSION(:,:), POINTER :: un
     REAL(KIND=8), DIMENSION(:,:), POINTER :: cij
     !===limiting
     INTEGER, DIMENSION(:,:), POINTER :: jj
     TYPE(mass_for_limiting) :: mass
     !===end limiting
     PROCEDURE(flux_template),       NOPASS, POINTER :: flux
     PROCEDURE(flux_template),       NOPASS, POINTER :: flux_prime
     PROCEDURE(lambda_max_template), NOPASS, POINTER :: lambda_max
   CONTAINS
     PROCEDURE, PUBLIC :: init => init_nl_scalar_cons
     PROCEDURE, PUBLIC :: READ => read_nl_scalar_cons
     PROCEDURE, PUBLIC :: update
  END TYPE nl_scalar_cons_type

CONTAINS
  SUBROUTINE init_nl_scalar_cons(this, flux, flux_prime, lambda_max, fourier_param, init_time)
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT) :: this
    PROCEDURE(flux_template) :: flux, flux_prime
    PROCEDURE(lambda_max_template) :: lambda_max
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8) :: init_time
    INTEGER :: Nmax, m, n, i, j

    CALL this%READ()
    this%ERK%sv = this%erk_sv
    CALL this%ERK%init
    this%flux => flux
    this%flux_prime => flux_prime
    this%lambda_max => lambda_max
    this%time = init_time
    Nmax = fourier_param%Nmax
    this%Nmax  = Nmax
    this%Nmax_real = 2*Nmax-1
    this%dx = fourier_param%dx
    this%lumped = fourier_param%dx
    ALLOCATE (this%un(Nmax,2))

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
  END SUBROUTINE init_nl_scalar_cons

  SUBROUTINE read_nl_scalar_cons(this)
    USE character_strings
    IMPLICIT NONE

    CHARACTER(1), PARAMETER :: begin_section ='~'
    CHARACTER(1), PARAMETER :: end_section   ='~'
    INTEGER, PARAMETER :: in_unit = 21
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: this
    TYPE(argument_nl_scalar_cons)  :: argument_data
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

    !===Glob min
    WRITE(string_default,*) this%glob_min
    CALL compare_string(record, list, argument_data%glob_min, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%glob_min
    END IF

    !===Glob max
    WRITE(string_default,*) this%glob_max
    CALL compare_string(record, list, argument_data%glob_max, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%glob_max
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
  END SUBROUTINE read_nl_scalar_cons

  SUBROUTINE update(this, fourier_param)
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: this
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(This%Nmax,2,this%ERK%s+1) :: urk
    REAL(KIND=8), DIMENSION(This%Nmax,2,this%ERK%s) :: flux_rk
    INTEGER  :: stage
    urk(:,:,1) = this%un
    DO stage = 2, this%ERK%s+1
       CALL one_step_ERK(this,stage,fourier_param,urk,flux_rk)
    END DO
    this%un = urk(:,:,this%ERK%s+1)
    this%time = this%time + this%dt
  END SUBROUTINE update

  SUBROUTINE one_step_ERK(this,stage,fourier_param,urk,flux_rk)
    USE fft_1D
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: this
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8), DIMENSION(this%Nmax,2,this%ERK%s+1) :: urk
    REAL(KIND=8), DIMENSION(This%Nmax,2,this%ERK%s) :: flux_rk
    REAL(KIND=8), DIMENSION(this%Nmax,2) :: cs_diff, cs_dflux, cs_zz
    REAL(KIND=8), DIMENSION(this%Nmax_real) :: u_max, u_min
    REAL(KIND=8), DIMENSION(this%Nmax_real,1) :: r_in, r_out
    INTEGER  :: stage, l, it
    !===Compute viscous flux, actual flux, and u_min, u_max (for limiting)
    CALL compute_dt_viscous_flux_min_max(this,stage,&
         urk(:,:,stage-1),urk(:,:,this%ERK%lp_of_l(stage)),&
         flux_rk(:,:,stage-1),cs_diff,u_max,u_min)

    !===Low-order
    IF (TRIM(ADJUSTL(this%method))=='viscous') THEN
       cs_dflux= this%ERK%inc_C(stage)*cs_diff
       urk(:,:,stage) = urk(:,:,this%ERK%lp_of_l(stage))+this%dt*cs_dflux
       RETURN
    END IF
    
    !===Higher-order
    IF (TRIM(ADJUSTL(this%method)).NE.'high') THEN
       WRITE(*,*) 'Bug in one_step_ERK, wrong method', this%method
       STOP
    END IF
    cs_zz =0.d0
    DO l = 1, stage-1
       cs_zz =  cs_zz + this%ERK%MatRK(stage,l)*flux_rk(:,:,l)
    END DO
    CALL fourier_derivative(cs_zz,cs_dflux,fourier_param%Length)

    IF (this%erk_sv>0) THEN
       cs_dflux= -cs_dflux + this%ERK%C(stage)*cs_diff
       urk(:,:,stage) = urk(:,:,1)+this%dt*cs_dflux
    ELSE
       cs_dflux= -cs_dflux + this%ERK%inc_C(stage)*cs_diff
       urk(:,:,stage) = urk(:,:,this%ERK%lp_of_l(stage))+this%dt*cs_dflux
    END IF

    !===Limiting
    IF (this%if_limiting) THEN
       !u_min = this%glob_min
       !u_max = this%glob_max
       CALL Fourier_to_real(urk(:,:,stage),r_in(:,1))
       DO it = 1, 1
          CALL iterative_cell_limiting_procedure(this%mass,this%jj,r_in,u_min,&
               psi_min,zero_of_psi_min,r_in)
          CALL iterative_cell_limiting_procedure(this%mass,this%jj,r_in,u_max,&
               psi_max,zero_of_psi_max,r_in)
       END DO
       CALL real_to_fourier(r_in(:,1),urk(:,:,stage))
    END IF
  END SUBROUTINE one_step_ERK

  SUBROUTINE compute_dt_viscous_flux_min_max(this,stage,urk_in,u_visc,cs_flux,cs_diff,umax,umin)
    USE fft_1D
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: this
    INTEGER, INTENT(IN) :: stage
    REAL(KIND = 8), DIMENSION(this%Nmax,2), INTENT(IN) :: urk_in, u_visc
    REAL(KIND = 8), DIMENSION(this%Nmax,2), INTENT(OUT) :: cs_diff, cs_flux
    REAL(KIND = 8), DIMENSION(this%Nmax_real)   :: r_out
    REAL(KIND = 8), DIMENSION(this%Nmax_real,2) :: dijL
    REAL(KIND = 8), DIMENSION(this%Nmax_real)   :: diag_dijL
    REAL(KIND = 8), DIMENSION(this%Nmax_real) :: r_diff, r_flux, umax, umin, alpha, beta, eta, etap
    REAL(KIND = 8) :: x, y, ul, ur, cij, lambda, uijbar, length
    INTEGER :: i, j, m, n, np

    CALL Fourier_to_real(u_visc,r_out)
    cij =0.5d0
    diag_dijL=0.d0
    DO m = 1, this%Nmax_real !===loop over cells
       i = this%jj(1,m)
       j = this%jj(2,m) 
       ul = r_out(i)
       ur = r_out(j)
       lambda = this%lambda_max(ul,ur)
       dijL(i,2) = cij*lambda
       dijL(j,1) = cij*lambda
       diag_dijL(i) = diag_dijL(i) - dijL(i,2)
       diag_dijL(j) = diag_dijL(j) - dijL(j,1)
    END DO

    IF (stage==2) THEN
       this%dt = this%ERK%s*0.5d0*this%CFL*this%lumped/MAXVAL(ABS(diag_dijL))
    END IF

    !===Only low-order 
    IF (TRIM(ADJUSTL(this%method))=='viscous') THEN
       r_diff = 0.d0
       r_flux = this%flux(r_out)
       DO m = 1, this%Nmax_real !===loop over cells
          DO n = 1, 2
             i = this%jj(n,m)
             np = MOD(n,2)+1
             j = this%jj(np,m)
             r_diff(i) = r_diff(i) + dijL(i,np)*(r_out(j)-r_out(i))
             r_diff(i) = r_diff(i) - this%cij(n,np)*(r_flux(j) - r_flux(i))
          END DO
       END DO
       r_diff = r_diff/this%lumped
       CALL real_to_fourier(r_diff,cs_diff)
       RETURN
    END IF
    
    !===High-order
    CALL entropy_commutator(this,r_out,alpha)

    IF (this%ERK%lp_of_l(stage).NE.stage-1) THEN
       CALL Fourier_to_real(urk_in,r_out)
    END IF
    umax = r_out
    umin = r_out   
    r_diff = 0.d0
    r_flux = this%flux(r_out)
    DO m = 1, this%Nmax_real !===loop over cells
       DO n = 1, 2
          i = this%jj(n,m)
          np = MOD(n,2)+1
          j = this%jj(np,m)
          r_diff(i) = r_diff(i) + 0.5*(alpha(i)+alpha(j))*dijL(i,np)*(r_out(j)-r_out(i))
          uijbar = 0.5d0*((r_out(j)+r_out(i)) &
               - (r_flux(j) - r_flux(i))*this%cij(n,np)/dijL(i,np))
          umax(i) = MAX(umax(i),uijbar)
          umin(i) = MIN(umin(i),uijbar)
          !umax(i) = MAX(umax(i),r_out(j))
          !umin(i) = MIN(umin(i),r_out(j))
       END DO
    END DO
    r_diff = r_diff/this%lumped
    CALL real_to_fourier(r_diff,cs_diff)
    r_flux = this%flux(r_out(1:this%Nmax_real)) !<-------FIX ME: remove this line
    CALL real_to_fourier(r_flux,cs_flux)

    IF (this%if_limiting) THEN
       CALL relax_min_and_max(this%bound_relaxing,this%glob_min,this%glob_max,this%jj,r_out,umax,umin)
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

  SUBROUTINE entropy_commutator(this,r_out,alpha)
    USE fft_1D
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: this
    REAL(KIND = 8), DIMENSION(this%Nmax_real):: r_out
    REAL(KIND = 8), DIMENSION(this%Nmax,2)   :: cs_diff, cs_flux
    REAL(KIND = 8), DIMENSION(this%Nmax_real):: alpha, beta, eta, etap, r_flux
    REAL(KIND = 8) :: Length, x, avg
    INTEGER :: i
    !===Compute entropy commutator (eta'*dx(f(u))-f'(u)*dx(eta(u)))
    length=this%dx*this%Nmax_real
    x = sum(r_out)/this%Nmax_real
    eta = (r_out-x)**2/2
    etap = (r_out-x)
    CALL real_to_fourier(eta,cs_diff)
    CALL fourier_derivative(cs_diff,cs_flux,length)
    CALL Fourier_to_real(cs_flux,alpha) !<==derivative of entropy
    alpha = alpha*this%flux_prime(r_out) !<==multiply by f'

    r_flux = this%flux(r_out)
    CALL real_to_fourier(r_flux,cs_flux)
    CALL fourier_derivative(cs_flux,cs_diff,length)
    CALL Fourier_to_real(cs_diff,beta) !<=derivative of flux
    beta= beta*etap !<==multiply by eta'

    avg = SUM(ABS(alpha))/this%Nmax_real
    alpha = min(abs(alpha - beta)/avg,1.d0)
    alpha = threshold(alpha)

  END SUBROUTINE entropy_commutator
  
    FUNCTION threshold(x) RESULT(g)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:)  :: x
    REAL(KIND=8), DIMENSION(SIZE(x))  :: z, t, zp, relu, f, g
    !===Quadratic threshold
!!$    REAL(KIND=8), PARAMETER :: x0 = 0.75d0, x1=SQRT(3.d0)*x0
!!$    z = x-x0
!!$    zp = x-2*x0
!!$    relu = (zp+abs(zp))/2
!!$    f = -z*(z**2-x1**2)  + relu*(z-x0)*(z+2*x0)
!!$    g = (f + 2*x0**3)/(4*x0**3)
!!$    !CALL plot_1d(x,g,'threshold1.plt') 
!!$
    !===Cubic threshold
    REAL(KIND=8), PARAMETER :: x0 = 0.75d0 !x0=0.1 (cubic threshold)
    relu = ((x-2*x0)+abs(x-2*x0))/2
    t = x/(2*x0)
    g = t**3*(10-15*t+6*t**2) - relu*(t-1)**2*(6*t**2+3*t+1)/(2*x0)
    RETURN
  END FUNCTION threshold
END MODULE nl_scalar_cons_module
