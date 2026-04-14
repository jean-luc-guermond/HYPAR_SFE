MODULE nl_scalar_cons_module
  USE Butcher_tableau
  USE fourier_param_module
  USE cell_limiting_engine_module
  PRIVATE
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

  TYPE, PUBLIC :: nl_scalar_cons_type
     REAL(KIND=8) :: CFL = 0.5d0, glob_min = -1.d20, glob_max = 1.d20
     LOGICAL :: if_limiting = .FALSE.
     CHARACTER(LEN=rec_length) :: bound_relaxing = 'minmod'
     CHARACTER(LEN=rec_length) :: method = 'viscous'
     INTEGER      :: erk_sv    = -21
     TYPE(BT), PUBLIC :: ERK
     TYPE(fourier_param_type), POINTER :: FP
     INTEGER :: Nmax, Nmax_real
     REAL(KIND=8) :: time, dt, lumped, dx, final_time
     REAL(KIND=8), DIMENSION(:,:), POINTER :: un
     REAL(KIND=8), DIMENSION(:,:), POINTER :: cij
     !===limiting
     INTEGER :: it_limiting_max = 2
     INTEGER, DIMENSION(:,:), POINTER :: jj
     TYPE(mass_for_limiting) :: mass
     !===end limiting
     PROCEDURE(flux_template),       NOPASS, POINTER :: flux
     PROCEDURE(flux_template),       NOPASS, POINTER :: flux_prime
     PROCEDURE(lambda_max_template), NOPASS, POINTER :: lambda_max
   CONTAINS
     PROCEDURE, PUBLIC :: init => init_nl_scalar_cons
     PROCEDURE, PUBLIC :: read => read_nl_scalar_cons
     PROCEDURE, PUBLIC :: update
  END TYPE nl_scalar_cons_type

CONTAINS
  SUBROUTINE init_nl_scalar_cons(this, flux, flux_prime, lambda_max, fourier_param, init_time, final_time)
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT) :: this
    PROCEDURE(flux_template) :: flux, flux_prime
    PROCEDURE(lambda_max_template) :: lambda_max
    TYPE(fourier_param_type), TARGET :: fourier_param
    REAL(KIND=8) :: init_time, final_time
    INTEGER :: Nmax, m, n, i, j

    CALL this%READ()
    this%ERK%sv = this%erk_sv
    CALL this%ERK%init
    this%flux => flux
    this%flux_prime => flux_prime
    this%lambda_max => lambda_max
    this%time = init_time
    this%final_time = final_time
    this%FP => fourier_param !<===FIXE ME: clear redundance with fourier_para
    Nmax = this%FP%Nmax
    this%Nmax  = Nmax
    this%Nmax_real = 2*Nmax-1
    this%dx = this%FP%dx
    this%lumped = this%FP%dx
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

    CHARACTER(LEN=rec_length) :: section_name='NL SCALAR CONS PARAMETERS'

    CLASS(nl_scalar_cons_type), INTENT(INOUT):: this
    TYPE(argument_nl_scalar_cons)  :: argument_data
    CHARACTER(LEN=rec_length)                         :: string

!================
!=== MANDATORY Reading all data file
!================
    CALL read_data_init_list(section_name)

!================
!=== We now find the relevant information for this setup
!================

    !===CFL
    CALL read_data(argument_data%CFL, this%CFL)

    !===ERK
    CALL read_data(argument_data%erk_sv, this%erk_sv)

    !===Higher-order vs. low-order
    CALL read_data(argument_data%method, this%method)

    !=========================
    !===Limiting parameters===
    !=========================
    !===if_limiting
    CALL read_data(argument_data%if_limiting, this%if_limiting)

    !===Glob min
    CALL read_data(argument_data%glob_min, this%glob_min)

    !===Glob max
    CALL read_data(argument_data%glob_max, this%glob_max)

    !===Bound_relaxing
    CALL read_data(argument_data%bound_relaxing, this%bound_relaxing)

!================
!=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
!================
     CALL finalize_rewrite_data

!================
!=== TEST DATA
!================

     IF (this%if_limiting .AND. this%method=='viscous') THEN
        WRITE(*,*) 'WARNING in input data: Limiting cannot be used with method=viscous'
     END IF

  END SUBROUTINE read_nl_scalar_cons

  SUBROUTINE update(this)
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: this
    REAL(KIND=8), DIMENSION(This%Nmax,2,this%ERK%s+1) :: urk
    REAL(KIND=8), DIMENSION(This%Nmax,2,this%ERK%s) :: flux_rk
    INTEGER  :: stage
    urk(:,:,1) = this%un
    DO stage = 2, this%ERK%s+1
       CALL one_step_ERK(this,stage,urk,flux_rk)
    END DO
    this%un = urk(:,:,this%ERK%s+1)
    this%time = this%time + this%dt
  END SUBROUTINE update

  SUBROUTINE one_step_ERK(nl_scalar_cons,stage,urk,flux_rk)
    USE fft_1D
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: nl_scalar_cons
    REAL(KIND=8), DIMENSION(nl_scalar_cons%Nmax,2,nl_scalar_cons%ERK%s+1) :: urk
    REAL(KIND=8), DIMENSION(nl_scalar_cons%Nmax,2,nl_scalar_cons%ERK%s) :: flux_rk
    REAL(KIND=8), DIMENSION(nl_scalar_cons%Nmax,2) :: cs_diff, cs_dflux, cs_zz
    REAL(KIND=8), DIMENSION(nl_scalar_cons%Nmax_real) :: u_max, u_min
    REAL(KIND=8), DIMENSION(nl_scalar_cons%Nmax_real,1) :: r_in, r_out
    INTEGER  :: stage, l, it
    !===Compute viscous flux, actual flux, and u_min, u_max (for limiting)
    CALL compute_dt_viscous_flux_min_max(nl_scalar_cons,stage,&
         urk(:,:,stage-1),urk(:,:,nl_scalar_cons%ERK%lp_of_l(stage)),&
         flux_rk(:,:,stage-1),cs_diff,u_max,u_min)

    !===Low-order
    IF (TRIM(ADJUSTL(nl_scalar_cons%method))=='viscous') THEN
         cs_dflux= nl_scalar_cons%ERK%inc_C(stage)*cs_diff
         urk(:,:,stage) = urk(:,:,nl_scalar_cons%ERK%lp_of_l(stage))+nl_scalar_cons%dt*cs_dflux
         RETURN
    !===Higher-order
    ELSE IF (TRIM(ADJUSTL(nl_scalar_cons%method))=='high') THEN
         cs_zz =0.d0
         DO l = 1, stage-1
            cs_zz =  cs_zz + nl_scalar_cons%ERK%MatRK(stage,l)*flux_rk(:,:,l)
         END DO
         CALL fourier_derivative(cs_zz,cs_dflux,nl_scalar_cons%FP%Length)

         IF (nl_scalar_cons%erk_sv>0) THEN
            cs_dflux= -cs_dflux + nl_scalar_cons%ERK%C(stage)*cs_diff
            urk(:,:,stage) = urk(:,:,1)+nl_scalar_cons%dt*cs_dflux
         ELSE
            cs_dflux= -cs_dflux + nl_scalar_cons%ERK%inc_C(stage)*cs_diff
            urk(:,:,stage) = urk(:,:,nl_scalar_cons%ERK%lp_of_l(stage))+nl_scalar_cons%dt*cs_dflux
         END IF

         !===Limiting
         IF (nl_scalar_cons%if_limiting) THEN
            CALL Fourier_to_real(urk(:,:,stage),r_in(:,1))
            DO it = 1, nl_scalar_cons%it_limiting_max   
               CALL iterative_cell_limiting_procedure(nl_scalar_cons%mass,nl_scalar_cons%jj,r_in,u_min,&
                     psi_min,zero_of_psi_min,r_in)
               CALL iterative_cell_limiting_procedure(nl_scalar_cons%mass,nl_scalar_cons%jj,r_in,u_max,&
                     psi_max,zero_of_psi_max,r_in)
            END DO
            CALL real_to_fourier(r_in(:,1),urk(:,:,stage))
         END IF
    ELSE
         WRITE(*,*) 'BUG in one_step_ERK, method should be "viscous" or "high", not ', nl_scalar_cons%method
         STOP
    END IF

  END SUBROUTINE one_step_ERK

  SUBROUTINE compute_dt_viscous_flux_min_max(nl_scalar_cons,stage,urk_in,u_visc,cs_flux,cs_diff,umax,umin)
    USE fft_1D
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: nl_scalar_cons
    INTEGER, INTENT(IN) :: stage
    REAL(KIND = 8), DIMENSION(nl_scalar_cons%Nmax,2), INTENT(IN) :: urk_in, u_visc
    REAL(KIND = 8), DIMENSION(nl_scalar_cons%Nmax,2), INTENT(OUT) :: cs_diff, cs_flux
    REAL(KIND = 8), DIMENSION(nl_scalar_cons%Nmax_real)   :: r_out
    REAL(KIND = 8), DIMENSION(nl_scalar_cons%Nmax_real,2) :: dijL
    REAL(KIND = 8), DIMENSION(nl_scalar_cons%Nmax_real)   :: diag_dijL
    REAL(KIND = 8), DIMENSION(nl_scalar_cons%Nmax_real) :: r_diff, r_flux, umax, umin, alpha, beta, eta, etap
    REAL(KIND = 8) :: x, y, ul, ur, cij, lambda, uijbar, length
    INTEGER :: i, j, m, n, np

    CALL Fourier_to_real(u_visc,r_out)
    cij =0.5d0
    diag_dijL=0.d0
    DO m = 1, nl_scalar_cons%Nmax_real !===loop over cells
       i = nl_scalar_cons%jj(1,m)
       j = nl_scalar_cons%jj(2,m) 
       ul = r_out(i)
       ur = r_out(j)
       lambda = nl_scalar_cons%lambda_max(ul,ur)
       dijL(i,2) = cij*lambda
       dijL(j,1) = cij*lambda
       diag_dijL(i) = diag_dijL(i) - dijL(i,2)
       diag_dijL(j) = diag_dijL(j) - dijL(j,1)
    END DO

    IF (stage==2) THEN
       nl_scalar_cons%dt = nl_scalar_cons%ERK%s*0.5d0*nl_scalar_cons%CFL*nl_scalar_cons%lumped/MAXVAL(ABS(diag_dijL))
    END IF

    !===Only low-order 
    IF (TRIM(ADJUSTL(nl_scalar_cons%method))=='viscous') THEN
       r_diff = 0.d0
       r_flux = nl_scalar_cons%flux(r_out)
       DO m = 1, nl_scalar_cons%Nmax_real !===loop over cells
          DO n = 1, 2
             i = nl_scalar_cons%jj(n,m)
             np = MOD(n,2)+1
             j = nl_scalar_cons%jj(np,m)
             r_diff(i) = r_diff(i) + dijL(i,np)*(r_out(j)-r_out(i))
             r_diff(i) = r_diff(i) - nl_scalar_cons%cij(n,np)*(r_flux(j) - r_flux(i))
          END DO
       END DO
       r_diff = r_diff/nl_scalar_cons%lumped
       CALL real_to_fourier(r_diff,cs_diff)
    !===High-order
    ELSE IF (TRIM(ADJUSTL(nl_scalar_cons%method))=='high') THEN
         CALL entropy_commutator(nl_scalar_cons,stage,r_out,alpha)

         IF (nl_scalar_cons%ERK%lp_of_l(stage).NE.stage-1) THEN
            CALL Fourier_to_real(urk_in,r_out)
         END IF
         umax = r_out
         umin = r_out   
         r_diff = 0.d0
         r_flux = nl_scalar_cons%flux(r_out)
         DO m = 1, nl_scalar_cons%Nmax_real !===loop over cells
            DO n = 1, 2
               i = nl_scalar_cons%jj(n,m)
               np = MOD(n,2)+1
               j = nl_scalar_cons%jj(np,m)
               r_diff(i) = r_diff(i) + 0.5*(alpha(i)+alpha(j))*dijL(i,np)*(r_out(j)-r_out(i))
               uijbar = 0.5d0*((r_out(j)+r_out(i)) &
                     - (r_flux(j) - r_flux(i))*nl_scalar_cons%cij(n,np)/dijL(i,np))
               umax(i) = MAX(umax(i),uijbar)
               umin(i) = MIN(umin(i),uijbar)
            END DO
         END DO
         r_diff = r_diff/nl_scalar_cons%lumped
         CALL real_to_fourier(r_diff,cs_diff)
         CALL real_to_fourier(r_flux,cs_flux)

         IF (nl_scalar_cons%if_limiting) THEN
            CALL relax_min_and_max(nl_scalar_cons%bound_relaxing,nl_scalar_cons%glob_min,nl_scalar_cons%glob_max,nl_scalar_cons%jj,r_out,umax,umin)
         END IF
      ELSE
         WRITE(*,*) 'BUG in compute_dt_viscous_flux_min_max, method should be &
         "viscous" or "high", not ', nl_scalar_cons%method
         STOP
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

  SUBROUTINE entropy_commutator(nl_scalar_cons,stage,r_out,alpha)
    USE fft_1D
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: nl_scalar_cons
    INTEGER :: stage
    REAL(KIND = 8), DIMENSION(nl_scalar_cons%Nmax_real):: r_out
    REAL(KIND = 8), DIMENSION(nl_scalar_cons%Nmax,2)   :: cs_diff, cs_flux
    REAL(KIND = 8), DIMENSION(nl_scalar_cons%Nmax_real):: alpha, beta, eta, etap, r_flux
    REAL(KIND = 8) :: Length, x, avg
    INTEGER :: i
    !===Compute entropy commutator (eta'*dx(f(u))-f'(u)*dx(eta(u)))
    length=nl_scalar_cons%dx*nl_scalar_cons%Nmax_real
    x = sum(r_out)/nl_scalar_cons%Nmax_real
    eta =  (r_out-x)**2/2
    etap = (r_out-x)
    CALL real_to_fourier(eta,cs_diff)
    CALL fourier_derivative(cs_diff,cs_flux,length)
    CALL Fourier_to_real(cs_flux,alpha) !<==derivative of entropy
    alpha = alpha*nl_scalar_cons%flux_prime(r_out) !<==multiply by f'

    r_flux = nl_scalar_cons%flux(r_out)
    CALL real_to_fourier(r_flux,cs_flux)
    CALL fourier_derivative(cs_flux,cs_diff,length)
    CALL Fourier_to_real(cs_diff,beta) !<=derivative of flux
    beta= beta*etap !<==multiply by eta'

    avg = SUM(ABS(alpha))/nl_scalar_cons%Nmax_real
    alpha = min(abs(alpha - beta)/avg,1.d0)
    alpha = threshold(alpha)
    IF (nl_scalar_cons%time+1.1*nl_scalar_cons%dt>nl_scalar_cons%final_time &
    .AND. stage==nl_scalar_cons%ERK%s+1) THEN
       CALL nl_scalar_cons%FP%plot_1d(alpha, 'commutator.plt')
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

END MODULE nl_scalar_cons_module
