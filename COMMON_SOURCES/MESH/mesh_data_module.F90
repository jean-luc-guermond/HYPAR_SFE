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
      CHARACTER(len = rec_length) :: k_dim             = '=== FE Space dimension ==='
      CHARACTER(len = rec_length) :: nb_refinement     = '=== Number of refinement steps ==='
      CHARACTER(len = rec_length) :: nb_bords = '=== How many pieces of periodic boundary? ==='
      CHARACTER(len = rec_length) :: list_periodic = '=== Indices of periodic boundaries and corresponding vectors==='
   END TYPE argument_mesh_data_type
   !===default value in simulation                       
   TYPE mesh_data_type                                                           
      CHARACTER(len = rec_length)    :: directory         = '.'                  
      CHARACTER(len = rec_length)    :: file_name         = 'mesh_name'          
      LOGICAL                        :: if_mesh_formatted = .TRUE.                            
      LOGICAL                        :: if_read_partition = .FALSE.              
      INTEGER                        :: nb_dom            = 1                    
      INTEGER, DIMENSION(:), POINTER :: list_dom                                 
      INTEGER                        :: type_fe           = 1    
      INTEGER                        :: k_dim             = -1
      INTEGER                        :: nb_refinement     = 0
      INTEGER                                  :: nb_bords      = 0
      INTEGER, DIMENSION(:, :), POINTER        :: list_periodic
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: vect_e
   CONTAINS                                                                      
      PROCEDURE, PUBLIC              :: READ => read_mesh_data
      PROCEDURE, PUBLIC              :: init => init_mesh_data                   
   END TYPE mesh_data_type   
CONTAINS

  SUBROUTINE init_mesh_data(this)
    CLASS(mesh_data_type), INTENT(INOUT) :: this
    CALL this%READ
  END SUBROUTINE init_mesh_data

  SUBROUTINE read_mesh_data(this)
    USE character_strings
    IMPLICIT NONE

    CHARACTER(LEN=rec_length)            :: section_name='MESH PARAMETERS'

    CLASS(mesh_data_type), INTENT(INOUT) :: this
    TYPE(argument_mesh_data_type)        :: argument_data

    !=== Reading all data file
    CALL read_data_init_list(section_name)

    !================
    !=== We now find the relevant information for the mesh
    !================

    !=== directory
    CALL read_data(argument_data%directory, this%directory)

    !=== mesh name
    CALL read_data(argument_data%file_name, this%file_name)

    !=== is mesh formatted
    CALL read_data(argument_data%if_mesh_formatted, this%if_mesh_formatted)

    !=== do we read metis partition
    CALL read_data(argument_data%if_read_partition, this%if_read_partition)

    !=== type of finite element
    CALL read_data(argument_data%type_fe, this%type_fe)

    !=== FE space dimension
    CALL read_data(argument_data%k_dim, this%k_dim)

    !=== number of refinement steps
    CALL read_data(argument_data%nb_refinement, this%nb_refinement)

    !=== number of subdomains in the mesh
    CALL read_data(argument_data%nb_dom, this%nb_dom)

    !=== list of subdomains in the mesh (special treatment since it is an array)
    ALLOCATE(this%list_dom(this%nb_dom))
    this%list_dom(1) = 1
    CALL read_data(argument_data%list_dom, this%list_dom)

    !=== Number of periodic boundaries 
    CALL read_data(argument_data%nb_bords, this%nb_bords)

    !=== List of periodic boundaries (has its specific subroutine, see character_strings.F90 module) 
    ALLOCATE(this%list_periodic(2, this%nb_bords))
    ALLOCATE(this%vect_e(this%k_dim, this%nb_bords))
    CALL read_periodic_data(argument_data%list_periodic, this%nb_bords, this%list_periodic, this%vect_e)
    !================
    !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
    !================
    CALL finalize_rewrite_data

  END SUBROUTINE read_mesh_data

END MODULE mesh_data_module
