MODULE mesh_data_module
  IMPLICIT NONE
  !===chain of characters that should appear in data file
   INTEGER, PARAMETER, PRIVATE :: rec_length=200
   TYPE argument_mesh_data_type                                                  
      CHARACTER(len = rec_length) :: directory         = '=== Name of directory for mesh file ==='
      CHARACTER(len = rec_length) :: file_name         = '=== Name of mesh file ==='         
      CHARACTER(len = rec_length) :: if_mesh_formatted = '=== Is the mesh formatted? (True/False) ==='        
      CHARACTER(len = rec_length) :: if_read_partition = '=== Do we read metis partition? (true/false) ==='
      CHARACTER(len = rec_length) :: nb_dom            = '=== Number of subdomains in the mesh ==='
      CHARACTER(len = rec_length) :: list_dom          = '=== List of subdomains in the mesh ==='
      CHARACTER(len = rec_length) :: type_fe           = '=== Type of finite element ===' 
      CHARACTER(len = rec_length) :: nb_refinement     = '=== Number of refinement steps ==='
   END TYPE argument_mesh_data_type
   !===default value in simulation                       
   TYPE mesh_data_type                                                           
      CHARACTER(len = rec_length)    :: directory         = '.'                  
      CHARACTER(len = rec_length)    :: file_name         = 'mesh_name'          
      LOGICAL                        :: if_mesh_formatted = .true.                            
      LOGICAL                        :: if_read_partition = .false.              
      INTEGER                        :: nb_dom            = 1                    
      INTEGER, DIMENSION(:), POINTER :: list_dom                                 
      INTEGER                        :: type_fe           = 1                    
      INTEGER                        :: nb_refinement     = 0                    
   CONTAINS                                                                      
!      PROCEDURE, PUBLIC              :: init                                    
      PROCEDURE, PUBLIC              :: read => read_mesh_data                   
   END TYPE mesh_data_type   
!!$   TYPE mesh_data_type
!!$      CHARACTER(len = 200) :: directory
!!$      CHARACTER(len = 200) :: file_name
!!$      LOGICAL :: if_mesh_formatted, if_HCT, if_read_partition
!!$      CHARACTER(len = 20) :: mesh_structure
!!$      INTEGER :: nb_dom
!!$      INTEGER, DIMENSION(:), POINTER :: list_dom
!!$      INTEGER :: type_fe, nb_refinement
!!$   CONTAINS
!!$      PROCEDURE, PUBLIC :: init
!!$      PROCEDURE, PUBLIC :: read => read_mesh_data
!!$   END TYPE mesh_data_type
CONTAINS
!!$   SUBROUTINE init(a)
!!$      CLASS(mesh_data_type), INTENT(INOUT) :: a
!!$      !===Logicals
!!$      a%if_mesh_formatted = .FALSE.
!!$      a%if_HCT = .false.
!!$      !===Characters
!!$      a%directory = '.'
!!$      a%file_name = 'gnu'
!!$      !===Integers
!!$      a%nb_dom = -1
!!$      a%type_fe = -1
!!$      a%nb_refinement = 0
!!$      a%type_fe = 2
!!$   END SUBROUTINE init

  SUBROUTINE read_mesh_data(this)
    USE character_strings
    USE petsc
    IMPLICIT NONE
    CLASS(mesh_data_type), INTENT(INOUT) :: this
    TYPE(argument_mesh_data_type)        :: argument_data
    INTEGER, PARAMETER :: in_unit = 21
    INTEGER, PARAMETER :: list_length=200
    CHARACTER(LEN=rec_length), DIMENSION(list_length) :: list, record
    CHARACTER(LEN=rec_length) :: control, st, string_default
    CHARACTER(LEN=15) :: end_section='%%% END SECTION' 
    LOGICAL :: okay
    INTEGER :: rank, ierr, record_size=0, i_list=0, section_counter=0, last_section_line=0, j

    CALL MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    OPEN(UNIT = in_unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')

    !===Read data file into record
    DO
       READ(in_unit,'(A)',END=100) control
       IF (TRIM(ADJUSTL(control))=='') CYCLE
       record_size = record_size+1
       record(record_size)=control
       WRITE(st,'(A15)') TRIM(ADJUSTL(control))
       IF (st==end_section) THEN
          section_counter=section_counter+1
          last_section_line = record_size
       END IF
    END DO
100 CONTINUE
    CLOSE(in_unit)

    !===Copy record into list
    list(1:last_section_line)=record(1:last_section_line)

    !===Now we reorganize record
    i_list = last_section_line
 
    i_list = i_list+1
    list(i_list) = '%%% BEGIN SECTION: MESH %%%'

    WRITE(string_default,*) TRIM(ADJUSTL(this%directory))
    CALL compare_string(record(last_section_line+1:), list, argument_data%directory, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%directory  
    END IF
    
    WRITE(string_default,*) TRIM(ADJUSTL(this%file_name))
    CALL compare_string(record(last_section_line+1:), list, argument_data%file_name, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%file_name  
    END IF
    
    WRITE(string_default,*) this%if_mesh_formatted
    CALL compare_string(record(last_section_line+1:), list, argument_data%if_mesh_formatted, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%if_mesh_formatted 
    END IF
    
    WRITE(string_default,*) this%if_read_partition
    CALL compare_string(record(last_section_line+1:), list, argument_data%if_read_partition,string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%if_read_partition
    END IF
    
    WRITE(string_default,*) this%type_fe
    CALL compare_string(record(last_section_line+1:), list, argument_data%type_fe, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%type_fe
    END IF
    
    WRITE(string_default,*) this%nb_refinement
    CALL compare_string(record(last_section_line+1:), list, argument_data%nb_refinement, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%nb_refinement
    END IF
    
    WRITE(string_default,*) this%nb_dom
    CALL compare_string(record(last_section_line+1:), list, argument_data%nb_dom, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%nb_dom
    END IF

    ALLOCATE(this%list_dom(this%nb_dom))
    this%list_dom(1) = 1
    WRITE(string_default,*) this%list_dom
    CALL compare_string(record(last_section_line+1:), list, argument_data%list_dom, string_default, okay, i_list, j)
    IF (okay) THEN
       READ (list(i_list), *) this%list_dom
    END IF

    i_list = i_list+1
    list(i_list) = '%%% END SECTION: MESH %%%'
    
    i_list = i_list+1
    list(i_list) = ''

  !===Record reorganized data
  OPEN(unit=in_unit,file='data',FORM='FORMATTED',STATUS='UNKNOWN')
  DO j = 1, i_list
     WRITE(in_unit,'(A)') TRIM(ADJUSTL(list(j)))
  END DO
  DO j = last_section_line+1, record_size
     IF (TRIM(ADJUSTL(record(j)))=='') CYCLE
     WRITE(in_unit,'(A)') TRIM(ADJUSTL(record(j)))
  END DO
  CLOSE(in_unit)

  
!!$    argument = '===Name of directory for mesh file==='
!!$    CALL find_string(in_unit, argument, test)
!!$    IF (test) THEN
!!$       READ (in_unit, *) this%directory
!!$    ELSE
!!$       CALL default_data(rank, in_unit, argument, '.')
!!$       this%directory = '.'
!!$    END IF
!!$
!!$    argument = "===Name of mesh file==="
!!$    CALL find_string(in_unit, argument, test)
!!$    IF (test) THEN
!!$       READ (in_unit, *) this%file_name
!!$    ELSE
!!$       CALL default_data(rank, in_unit, argument, 'mesh_name')
!!$       this%file_name = 'mesh_name'
!!$       IF (rank == 0) WRITE(*, *) "No mesh_name specified." ; STOP
!!$    END IF
!!$
!!$    argument = "===Is the mesh formatted? (True/False)==="
!!$    CALL find_string(in_unit, argument, test)
!!$    IF (test) THEN
!!$       READ (in_unit, *) this%if_mesh_formatted
!!$    ELSE
!!$       CALL default_data(rank, in_unit, argument, '.t.')
!!$       this%if_mesh_formatted = .true.
!!$    END IF
!!$
!!$    argument = '===Do we read metis partition? (true/false)'
!!$    CALL find_string(in_unit, argument, test)
!!$    IF (test) THEN
!!$       READ (in_unit, *) this%if_read_partition
!!$    ELSE
!!$       CALL default_data(rank, in_unit, argument, '.f.')
!!$       this%if_read_partition = .false.
!!$    END IF
!!$
!!$    argument = '===Number of subdomains in the mesh==='
!!$    CALL find_string(in_unit, argument, test)
!!$    IF (test) THEN
!!$       READ (in_unit, *) this%nb_dom
!!$    ELSE
!!$       CALL default_data(rank, in_unit, argument, '1')
!!$       this%nb_dom = 1
!!$    END IF
!!$
!!$    argument = '===List of subdomain in the mesh==='
!!$    ALLOCATE(this%list_dom(this%nb_dom))
!!$    CALL find_string(in_unit, argument, test)
!!$    IF (test) THEN
!!$       READ (in_unit, *) this%list_dom
!!$    ELSE
!!$       CALL default_data(rank, in_unit, argument, '1')
!!$       this%list_dom(1) = 1
!!$    END IF
!!$
!!$    argument = '===Number of refinement steps==='
!!$    ALLOCATE(this%list_dom(this%nb_dom))
!!$    CALL find_string(in_unit, argument, test)
!!$    IF (test) THEN
!!$       READ (in_unit, *) this%nb_refinement
!!$    ELSE
!!$       CALL default_data(rank, in_unit, argument, '0')
!!$       this%nb_refinement = 0
!!$    END IF

!!$    argument = "===HCT mesh ?==="
!!$    CALL find_string(in_unit, argument, test)
!!$    IF (test) THEN
!!$       READ (in_unit, *) this%if_HCT
!!$       IF (this%if_HCT) write(*, *) "HCT mesh not inmplemented yet"
!!$    ELSE
!!$       this%if_HCT = .false.
!!$    END IF

!!$    argument = '===Type of finite element==='
!!$    CALL find_string(in_unit, argument, test)
!!$    IF (test) THEN
!!$       READ (in_unit, *) this%type_fe
!!$    ELSE
!!$       CALL default_data(rank, in_unit, argument, '1')
!!$       this%type_fe = 1
!!$    END IF

!!$    CLOSE(in_unit)
  END SUBROUTINE read_mesh_data

  SUBROUTINE compare_string(record, list, string, string_default, okay, i_list, j)
    IMPLICIT NONE
    CHARACTER(LEN=*), DIMENSION(:) :: record, list
    CHARACTER(LEN=*)               :: string, string_default
    LOGICAL                        :: okay
    INTEGER, INTENT(OUT)           :: j
    INTEGER                        :: i, i_list
    okay = .TRUE.
    i_list = i_list+1
    list(i_list) = string
    DO i = 1, SIZE(record)
       IF (TRIM(ADJUSTL(record(i)))==list(i_list)) THEN
          j = i
          record(j) = ''
          i_list = i_list + 1
          list(i_list) = record(j+1)
          record(j+1) = ''
          RETURN
       END IF
    END DO
    WRITE(*,*) ' File reading error '
    i_list = i_list+1
    list(i_list) = string_default
    okay = .FALSE.
    j = -1
    RETURN
  END SUBROUTINE compare_string
    
END MODULE mesh_data_module
