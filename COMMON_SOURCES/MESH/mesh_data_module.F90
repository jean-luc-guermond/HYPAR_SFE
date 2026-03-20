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
      PROCEDURE, PUBLIC              :: read => read_mesh_data                   
   END TYPE mesh_data_type   
CONTAINS

  SUBROUTINE read_mesh_data(this)
    USE character_strings
    USE petsc
    IMPLICIT NONE
    INTEGER, PARAMETER :: in_unit = 21
    INTEGER, PARAMETER :: list_length=200, length_begin=27, length_end=25
    CHARACTER(LEN=length_begin), PARAMETER :: begin_section ='%%% BEGIN SECTION: MESH %%%' 
    CHARACTER(LEN=length_end),   PARAMETER :: end_section   ='%%% END SECTION: MESH %%%'
    CHARACTER(LEN=length_begin), PARAMETER :: char_begin    ='%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    CHARACTER(LEN=length_end),   PARAMETER :: char_end      ='%%%%%%%%%%%%%%%%%%%%%%%%%'
    
    CLASS(mesh_data_type), INTENT(INOUT) :: this
    TYPE(argument_mesh_data_type)        :: argument_data
    
    CHARACTER(LEN=rec_length), DIMENSION(list_length) :: list, record
    CHARACTER(LEN=rec_length)   :: string_default
    LOGICAL :: okay
    INTEGER :: rank, ierr, record_size, i_list, j
  
    !===Initialize data to zero and false by default
    list = ""
    record = ""
    CALL MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    !===Initializing record
    CALL read_data_in_record(record_size, record, begin_section, end_section)

    !===Now we reorganize record
    i_list = 1
    WRITE(list(i_list), '(A)') REPEAT('|',70)
    i_list = i_list + 1
    list(i_list) = char_begin
    i_list = i_list + 1
    list(i_list) = begin_section
    i_list = i_list + 1
    list(i_list) = char_begin

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
    list(i_list) = char_end
    i_list = i_list+1
    list(i_list) = end_section
    i_list = i_list+1
    list(i_list) = char_end
    i_list = i_list+1
  
    !===Closing unit 
    CALL rewrite_data_from_list_record(rank, list, record, i_list, record_size)
  END SUBROUTINE read_mesh_data

END MODULE mesh_data_module
