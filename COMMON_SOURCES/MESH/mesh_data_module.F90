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
    INTEGER, PARAMETER :: list_length=200, length_begin=27, length_end=25
    CHARACTER(LEN=rec_length), DIMENSION(list_length) :: list, record
    CHARACTER(LEN=rec_length) :: control, st, string_default
    CHARACTER(LEN=length_begin) :: begin_section='%%% BEGIN SECTION: MESH %%%' 
    CHARACTER(LEN=length_end) :: end_section='%%% END SECTION: MESH %%%'
    CHARACTER(LEN=5)         :: fmt 
    LOGICAL :: okay
    INTEGER :: rank, ierr, record_size=0, i_list=0, section_counter=0, last_section_line=0, j
    INTEGER :: line_begin_section=-1, line_end_section=-1
    CALL MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    OPEN(UNIT = in_unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')

    !===Read data file into record
    DO
       READ(in_unit,'(A)',END=100) control
       IF (TRIM(ADJUSTL(control))=='') CYCLE
       record_size = record_size+1
       record(record_size)=control
       write(fmt, '("(A", I0, ")")') length_begin
       WRITE(st,fmt) TRIM(ADJUSTL(control))
       IF (st==begin_section) THEN
          line_begin_section = record_size
       END IF
       write(fmt, '("(A", I0, ")")') length_end
       WRITE(st,fmt) TRIM(ADJUSTL(control))
       IF (st==end_section) THEN
          line_end_section = record_size
       END IF
    END DO
100 CONTINUE
    CLOSE(in_unit)


    IF (line_begin_section .NE. -1) THEN
        record(line_begin_section:record_size-1) = record(line_begin_section+1: record_size)
        record_size = record_size -1
        IF (line_end_section > line_begin_section) line_end_section = line_end_section -1
    ELSE
        line_begin_section = 1
    END IF
    IF (line_end_section .NE. -1) THEN
        record(line_end_section:record_size-1) = record(line_end_section+1: record_size)
        record_size = record_size -1
    ELSE
        line_end_section = record_size
    END IF

    !===Now we reorganize record

    i_list = 1
    list(i_list) = '%%% BEGIN SECTION: MESH %%%'

    WRITE(string_default,*) TRIM(ADJUSTL(this%directory))
    CALL compare_string(record, list, argument_data%directory, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%directory  
    END IF
    
    WRITE(string_default,*) TRIM(ADJUSTL(this%file_name))
    CALL compare_string(record, list, argument_data%file_name, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%file_name  
    END IF
    
    WRITE(string_default,*) this%if_mesh_formatted
    CALL compare_string(record, list, argument_data%if_mesh_formatted, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%if_mesh_formatted 
    END IF
    
    WRITE(string_default,*) this%if_read_partition
    CALL compare_string(record, list, argument_data%if_read_partition,string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%if_read_partition
    END IF
    
    WRITE(string_default,*) this%type_fe
    CALL compare_string(record, list, argument_data%type_fe, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%type_fe
    END IF
    
    WRITE(string_default,*) this%nb_refinement
    CALL compare_string(record, list, argument_data%nb_refinement, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%nb_refinement
    END IF
    
    WRITE(string_default,*) this%nb_dom
    CALL compare_string(record, list, argument_data%nb_dom, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%nb_dom
    END IF

    ALLOCATE(this%list_dom(this%nb_dom))
    this%list_dom(1) = 1
    WRITE(string_default,*) this%list_dom
    CALL compare_string(record, list, argument_data%list_dom, string_default, okay, i_list, j)
    IF (okay) THEN
       READ (list(i_list), *) this%list_dom
    END IF

    i_list = i_list+1
    list(i_list) = '%%% END SECTION: MESH %%%'
    
    i_list = i_list+1
    list(i_list) = ''

    OPEN(unit=in_unit,file='data',FORM='FORMATTED',STATUS='UNKNOWN')
    IF (rank == 0) THEN 
       DO j = 1, record_size
          IF (TRIM(ADJUSTL(record(j)))=='') CYCLE
          WRITE(in_unit,'(A)') TRIM(ADJUSTL(record(j)))
       END DO
       DO j = 1, i_list
          WRITE(in_unit,'(A)') TRIM(ADJUSTL(list(j)))
       END DO
    END IF
    CLOSE(in_unit)
  
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
