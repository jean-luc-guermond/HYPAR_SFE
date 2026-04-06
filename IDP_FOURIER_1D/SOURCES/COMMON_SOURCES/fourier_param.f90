MODULE fourier_param_module
  INTEGER, PARAMETER, PRIVATE :: rec_length = 200, list_length=200
  TYPE argument_fourier_param_type
     CHARACTER(LEN=rec_length) :: Nmax       = '=== Number of Fourier modes ? ==='
     CHARACTER(LEN=rec_length) :: Length     = '=== Fourier periodic domain Length ? ==='
  END type argument_fourier_param_type
  TYPE fourier_param_type
     INTEGER :: Nmax = 10
     REAL(KIND=8) :: Length=1.d0
     INTEGER :: Nmax_real
     REAL(KIND=8) :: dx
      REAL(KIND=8), DIMENSION(:), POINTER :: rr
   CONTAINS
     PROCEDURE, PUBLIC :: init => init_fourier_param
     PROCEDURE, PUBLIC :: read => read_fourier_param
     PROCEDURE, PUBLIC :: plot_1d => plot_1d_fourier_param
  END type fourier_param_type
CONTAINS
  SUBROUTINE init_fourier_param(this)
    IMPLICIT NONE
    CLASS(fourier_param_type), INTENT(INOUT) :: this
    INTEGER :: i
    CALL this%read()
    this%Nmax_real = 2*this%Nmax-1
    this%dx = this%length/this%Nmax_real
    ALLOCATE(this%rr(this%Nmax_real))
    DO i = 1, this%Nmax_real
       this%rr(i) = (i-1)*this%dx
    END DO
  END SUBROUTINE init_fourier_param

  SUBROUTINE read_fourier_param(this)
    USE character_strings
    IMPLICIT NONE

    CHARACTER(LEN=rec_length) :: section_name='FOURIER PARAMETERS'

    CLASS(fourier_param_type), INTENT(INOUT):: this
    TYPE(argument_fourier_param_type)  :: argument_data
    CHARACTER(LEN=rec_length)                         :: string

!================
!=== MANDATORY Reading all data file
!================
    CALL read_data_init_list(section_name)
    
!================
!=== We now find the relevant information for this setup
!================

    !===nb Fourier modes
    CALL read_data(argument_data%Nmax, this%Nmax)

    !===Length
    CALL read_data(argument_data%Length, this%Length)

!================
!=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
!================
     CALL finalize_rewrite_data
     
  END SUBROUTINE read_fourier_param

  SUBROUTINE plot_1d_fourier_param(this,un,file)
    IMPLICIT NONE
    CLASS(fourier_param_type), INTENT(INOUT):: this
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: un
    CHARACTER(*) :: file
    INTEGER :: n, unit=10
    OPEN(unit,FILE=TRIM(ADJUSTL(file)),FORM='formatted')
    !WRITE(unit,*) '%toplabel='' '''
    !WRITE(unit,*) '%xlabel='' '''
    !WRITE(unit,*) '%ylabel='' '''
    !WRITE(unit,*) '%ymax=', MAXVAL(un)
    !WRITE(unit,*) '%ymin=', MINVAL(un)
    !WRITE(unit,*) '%xyratio=1'
    !WRITE(unit,*) '%mt=4'
    !WRITE(unit,*) '%mc=2'
    !WRITE(unit,*) '%lc=2'
    DO n = 1, SIZE(this%rr)
       WRITE(unit,*) this%rr(n), un(n)
    END DO
    CLOSE(unit)
  END SUBROUTINE plot_1d_fourier_param
  
END MODULE fourier_param_module
