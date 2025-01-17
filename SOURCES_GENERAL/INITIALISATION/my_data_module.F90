MODULE my_data_module
  IMPLICIT NONE
  TYPE my_data
     LOGICAL                        :: if_regression_test
     CHARACTER(len=200)             :: directory
     CHARACTER(len=200)             :: file_name
     LOGICAL                        :: if_mesh_formatted
     CHARACTER(len=20)              :: mesh_structure
     LOGICAL                        :: if_restart
     REAL(KIND=8)                   :: checkpointing_freq
     INTEGER                        :: nb_dom
     INTEGER, DIMENSION(:), POINTER :: list_dom
     INTEGER                        :: type_fe
     INTEGER                        :: RKsv

     REAL(KIND=8)                   :: Tfinal
     REAL(KIND=8)                   :: CFL
     REAL(KIND=8)                   :: dt, time
     INTEGER                        :: syst_size
     LOGICAL                        :: if_lumped
     CHARACTER(LEN=30)              :: max_viscosity
     CHARACTER(LEN=30)              :: equation_of_state
     REAL(KIND=8)                   :: b_covolume
     CHARACTER(LEN=30)              :: method_type
     CHARACTER(LEN=30)              :: high_order_viscosity
     LOGICAL                        :: if_convex_limiting
     LOGICAL                        :: if_relax_bounds
     CHARACTER(len=30)              :: limiter_type
     REAL(KIND=8)                   :: ce
     CHARACTER(LEN=4)               :: type_test
     INTEGER                        :: rho_nb_Dir_bdy, ux_nb_Dir_bdy, uy_nb_Dir_bdy
     INTEGER, DIMENSION(:), POINTER :: rho_Dir_list, ux_Dir_list, uy_Dir_list
     INTEGER                        :: nb_udotn_zero
     INTEGER, DIMENSION(:), POINTER :: udotn_zero_list
     INTEGER                        :: nb_DIR
     INTEGER, DIMENSION(:), POINTER :: DIR_list
     !===Plot
     REAL(KIND=8)                   :: dt_plot
   CONTAINS
     PROCEDURE, PUBLIC              :: init
  END TYPE my_data
CONTAINS
  SUBROUTINE init(a)
    CLASS(my_data), INTENT(INOUT) :: a
    !===Logicals
    a%if_regression_test = .FALSE.
    a%if_mesh_formatted = .FALSE.
    a%if_restart = .FALSE.
    a%if_lumped = .TRUE.
    a%if_convex_limiting =.TRUE.
    a%if_relax_bounds = .TRUE.
    !===Reals
    a%ce=1.d0
    a%CFL=1.d0
    a%checkpointing_freq=1.d20
    a%dt_plot=1.d20
    a%b_covolume = 0.d0
    !===Characters
    a%directory='.'
    a%limiter_type='avg'
    a%file_name='gnu'
    a%type_test='gnu'
    a%equation_of_state='gamma-law'
    a%max_viscosity='none'
    a%mesh_structure=''
    !===Integers
    a%nb_dom=-1
    a%type_fe=-1
    a%RKsv=-1
  END SUBROUTINE init
END MODULE my_data_module

MODULE input_data
  USE my_data_module
  IMPLICIT NONE
  PUBLIC :: read_my_data
  TYPE(my_data), PUBLIC  :: inputs
  PRIVATE
CONTAINS
  SUBROUTINE read_my_data(data_fichier)
    USE character_strings
    USE space_dim
    IMPLICIT NONE
    INTEGER, PARAMETER           :: in_unit=21
    CHARACTER(len=*), INTENT(IN) :: data_fichier
    CHARACTER(LEN=100)           :: argument
    LOGICAL :: okay
    !===Initialize data to zero and false by default
    CALL inputs%init

    OPEN(UNIT = in_unit, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(in_unit, "===Name of directory for mesh file===")
    READ (in_unit,*) inputs%directory
    CALL read_until(in_unit, "===Name of mesh file===")
    READ (in_unit,*) inputs%file_name
    CALL read_until(in_unit, "===Is the mesh formatted? (True/False)===")
    READ (in_unit,*) inputs%if_mesh_formatted
    CALL read_until(in_unit, '===Number of subdomains in the mesh===')
    READ(21,*) inputs%nb_dom
    ALLOCATE(inputs%list_dom(inputs%nb_dom))
    CALL read_until(21, '===List of subdomain in the mesh===')
    READ(21,*) inputs%list_dom

    !===Regression test
    CALL getarg(1,argument)
    IF (trim(adjustl(argument))=='regression') THEN
       inputs%if_regression_test = .true.
    ELSE
       inputs%if_regression_test = .false.
    END IF

    CLOSE(in_unit)
  END SUBROUTINE read_my_data
END MODULE input_data

