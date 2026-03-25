MODULE nl_scalar_cons_module
  USE Butcher_tableau
  USE fourier_param_module
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
     CHARACTER(LEN=rec_length) :: CFL       = '=== CFL ? ==='
     CHARACTER(LEN=rec_length) :: erk_sv    = '=== ERK ? ==='
  END TYPE argument_nl_scalar_cons

  TYPE nl_scalar_cons_type
     REAL(KIND=8) :: CFL = 0.5d0
     INTEGER      :: erk_sv    = -21
     TYPE(BT), PUBLIC :: ERK
     INTEGER :: Nmax, Nmax_real
     REAL(KIND=8) :: time, dt, lumped, dx
     REAL(KIND=8) , DIMENSION(:,:), POINTER :: un
     PROCEDURE(flux_template),       NOPASS, POINTER :: flux
     PROCEDURE(lambda_max_template), NOPASS, POINTER :: lambda_max
   CONTAINS
     PROCEDURE, PUBLIC :: init => init_nl_scalar_cons
     PROCEDURE, PUBLIC :: READ => read_nl_scalar_cons
     PROCEDURE, PUBLIC :: update
  END TYPE nl_scalar_cons_type


CONTAINS
  SUBROUTINE init_nl_scalar_cons(this, flux, lambda_max, fourier_param, init_time)
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT) :: this
    PROCEDURE(flux_template) :: flux
    PROCEDURE(lambda_max_template) :: lambda_max
    TYPE(fourier_param_type) :: fourier_param
    REAL(KIND=8) :: init_time
    INTEGER :: Nmax
 
    
    CALL this%READ()
    CALL this%ERK%init(this%erk_sv)
    this%flux => flux
    this%lambda_max => lambda_max
    this%time = init_time
    Nmax = fourier_param%Nmax
    this%Nmax  = Nmax
    this%Nmax_real = 2*Nmax-1
    this%dx = fourier_param%dx
    this%lumped = fourier_param%dx
    ALLOCATE (this%un(Nmax,2))
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
    INTEGER  :: stage, l
    CALL compute_dt_viscous_flux_min_max(this,stage,urk(:,:,stage-1),flux_rk(:,:,stage-1),&
         cs_diff,u_max,u_min)
    !===u_max,u_min will be used to do the limiting (someday)

    cs_zz =0.d0
    DO l = 1, stage-1
       cs_zz =  cs_zz + this%ERK%MatRK(stage,l)*flux_rk(:,:,l)
    END DO
    CALL fourier_derivative(cs_zz,cs_dflux,fourier_param%Length)
    cs_dflux= -cs_dflux + this%ERK%inc_C(stage)*cs_diff
    urk(:,:,stage) = urk(:,:,1)+this%dt*cs_dflux
  END SUBROUTINE one_step_ERK

  SUBROUTINE compute_dt_viscous_flux_min_max(this,stage,urk_in,cs_flux,cs_diff,umax,umin)
    USE fft_1D
    IMPLICIT NONE
    CLASS(nl_scalar_cons_type), INTENT(INOUT):: this
    INTEGER, INTENT(IN) :: stage
    REAL(KIND=8), DIMENSION(this%Nmax,2), INTENT(IN) :: urk_in
    REAL(KIND=8), DIMENSION(this%Nmax,2), INTENT(OUT) :: cs_diff, cs_flux
    REAL(KIND = 8), DIMENSION(0:this%Nmax_real+1)   :: r_out
    REAL(KIND = 8), DIMENSION(0:this%Nmax_real+1,2) :: dijL
    REAL(KIND = 8), DIMENSION(0:this%Nmax_real+1)   :: diag_dijL
    REAL(KIND = 8), DIMENSION(this%Nmax_real) :: r_diff, r_flux, umax, umin
    REAL(KIND = 8) :: ul, ur, cij, lambda, dx
    INTEGER :: i, j, m

    CALL Fourier_to_real(urk_in,r_out(1:this%Nmax_real))
    r_out(0) = r_out(this%Nmax_real)
    r_out(this%Nmax_real+1)=r_out(1)

    cij =0.5d0
    diag_dijL=0.d0
    DO m = 0, this%Nmax_real
       i = m
       j = m+1
       ul = r_out(i)
       ur = r_out(j)
       lambda = this%lambda_max(ul,ur)
       dijL(i,2) = cij*lambda
       dijL(j,1) = cij*lambda
       diag_dijL(i) = diag_dijL(i) - dijL(i,2)
       diag_dijL(j) = diag_dijL(j) - dijL(j,1)
    END DO

    IF (stage==2) THEN
       this%dt = this%ERK%s*0.5d0*this%CFL*this%lumped/MAXVAL(ABS(diag_dijL(1:this%Nmax_real)))
    END IF

    umax = r_out(1:this%Nmax_real)
    umin = r_out(1:this%Nmax_real)
    DO i = 1, this%Nmax_real
       r_diff(i) = dijL(i,1)*(r_out(i-1)-r_out(i)) + dijL(i,2)*(r_out(i+1)-r_out(i))
       umax(i) = MAX(umax(i),r_out(i-1),r_out(i+1))
       umin(i) = MIN(umin(i),r_out(i-1),r_out(i+1))
    END DO
    r_diff = r_diff*this%dx**2/this%lumped
    CALL real_to_fourier(r_diff,cs_diff)

    r_flux = this%flux(r_out(1:this%Nmax_real))
    CALL real_to_fourier(r_flux,cs_flux)

  END SUBROUTINE compute_dt_viscous_flux_min_max

END MODULE nl_scalar_cons_module
