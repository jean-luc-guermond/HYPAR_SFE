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
    CALL find_string(in_unit, "===Mesh structure ('','POWELL_SABIN','HCT')===",okay)
    IF (okay) THEN
       READ (in_unit,*) inputs%mesh_structure
    END IF
    
    CALL read_until(21, '===Type of finite element===')
    READ(21,*) inputs%type_fe
    CALL find_string(in_unit, "===Equation of state===",okay)
    IF (okay) THEN
       READ (in_unit,*) inputs%equation_of_state
    END IF

    IF (inputs%equation_of_state.NE.'gamma-law') THEN
       CALL read_until(21, '===Covolume coefficient b===')
       READ(in_unit,*) inputs%b_covolume
    END IF

    CALL read_until(in_unit, "===RK method type===")
    READ (in_unit,*) inputs%RKsv
    CALL read_until(in_unit, "===Final time===") 
    READ (in_unit,*) inputs%Tfinal
    CALL read_until(in_unit, "===CFL number===") 
    READ (in_unit,*) inputs%CFL

    CALL find_string(in_unit, "===Maximum wave speed? 'none', 'local_max', 'global_max'===",okay)
    IF (okay) THEN
       READ (in_unit,*) inputs%max_viscosity
    ELSE
       write(*,*) 'inputs%max_viscosity set to none'
       inputs%max_viscosity = 'none'
    END IF

    !===Restart
    CALL find_string(in_unit, "===Restart (true/false)===",okay)
    IF (okay) THEN
       READ (in_unit,*) inputs%if_restart
    END IF
    CALL find_string(in_unit, "===Checkpointing frequency===",okay)
    IF (okay) THEN
       READ (in_unit,*) inputs%checkpointing_freq
    END IF

    CALL read_until(in_unit, "===Method type (galerkin, viscous, high)===")
    READ (in_unit,*) inputs%method_type
    SELECT CASE(inputs%method_type)
    CASE('high')
       CALL read_until(in_unit, "===High-order method===")
       READ (in_unit,*) inputs%high_order_viscosity
       CALL read_until(in_unit, "===Convex limiting (true/false)===")
       READ (in_unit,*) inputs%if_convex_limiting
       IF (inputs%if_convex_limiting) THEN
          CALL read_until(in_unit, "===Limiter type (avg, minmod)===")
          READ (in_unit,*) inputs%limiter_type
          CALL read_until(in_unit, "===Relax bounds?===")
          READ (in_unit,*)  inputs%if_relax_bounds
       END IF
    END SELECT

    SELECT CASE(inputs%method_type)
    CASE('galerkin','high')
       CALL read_until(in_unit, "===Mass matrix lumped (true/false)===")
       READ (in_unit,*) inputs%if_lumped
    CASE DEFAULT
       inputs%if_lumped = .TRUE.
    END SELECT

    SELECT CASE(inputs%high_order_viscosity)
    CASE('EV(p)','EV(s)')
       CALL read_until(in_unit, "===ce coefficient===") 
       READ (in_unit,*) inputs%ce
    CASE DEFAULT
       inputs%ce=0.d0
    END SELECT

    CALL find_string(in_unit, "===Test case name===",okay)
    IF (okay) THEN
       READ (in_unit,*) inputs%type_test
    ELSE
       inputs%type_test='NONE'
    END IF

    !===Boundary conditions
    CALL read_until(in_unit, "===How many Dirichlet boundaries for rho?===")
    READ (in_unit,*)  inputs%rho_nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries for rho?===")
    ALLOCATE(inputs%rho_Dir_list(inputs%rho_nb_Dir_bdy))
    READ (in_unit,*) inputs%rho_Dir_list
    CALL read_until(in_unit, "===How many Dirichlet boundaries for ux?===")
    READ (in_unit,*)  inputs%ux_nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries for ux?===")
    ALLOCATE(inputs%ux_Dir_list(inputs%ux_nb_Dir_bdy))
    READ (in_unit,*) inputs%ux_Dir_list
    IF (k_dim==2) THEN
       CALL read_until(in_unit, "===How many Dirichlet boundaries for uy?===")
       READ (in_unit,*)  inputs%uy_nb_Dir_bdy
       CALL read_until(in_unit, "===List of Dirichlet boundaries for uy?===")
       ALLOCATE(inputs%uy_Dir_list(inputs%uy_nb_Dir_bdy))
       READ (in_unit,*) inputs%uy_Dir_list
    ELSE
       inputs%uy_nb_Dir_bdy=0
       ALLOCATE(inputs%uy_Dir_list(inputs%uy_nb_Dir_bdy))
    END IF
    
    CALL find_string(in_unit, "===How many Diriclet boundaries (subsonic/supersonic)?===",okay)
    IF (okay) THEN
       READ (in_unit,*) inputs%nb_DIR
       IF (inputs%nb_DIR>0) THEN
          CALL read_until(in_unit, "===List of Diriclet boundaries (subsonic/supersonic)?===")
          ALLOCATE(inputs%DIR_list(inputs%nb_DIR))
          READ (in_unit,*) inputs%DIR_list
       ELSE
          inputs%nb_DIR = 0
          ALLOCATE(inputs%DIR_list(inputs%nb_DIR))
       END IF
    ELSE
       inputs%nb_DIR = 0
       ALLOCATE(inputs%DIR_list(inputs%nb_DIR))
    END IF

    !===Normal boundary condition
    CALL find_string(in_unit, "===How many boundaries for u.n=0?===",okay)
    IF (okay) THEN
       READ (in_unit,*)  inputs%nb_udotn_zero
       IF (inputs%nb_udotn_zero>0) THEN
          CALL read_until(in_unit, "===List of boundarie for u.n=0?===")
          ALLOCATE(inputs%udotn_zero_list(inputs%nb_udotn_zero))
          READ (in_unit,*) inputs%udotn_zero_list
       ELSE
          inputs%nb_udotn_zero = 0
          ALLOCATE(inputs%udotn_zero_list(inputs%nb_udotn_zero))
       END IF
    ELSE
       inputs%nb_udotn_zero = 0
       ALLOCATE(inputs%udotn_zero_list(inputs%nb_udotn_zero))
    END IF

    !===Postprocessing
    CALL find_string(in_unit, "===dt_plot?===",okay)
    IF (okay) THEN
       READ (in_unit,*)  inputs%dt_plot
    END IF

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


