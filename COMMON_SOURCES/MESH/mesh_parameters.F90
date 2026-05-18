MODULE mesh_data_module
   USE read_inputs_module, ONLY: rec_length
   IMPLICIT NONE
  !===chain of characters that should appear in data file

   TYPE argument_mesh_data_type
      CHARACTER(len = rec_length) :: directory         = '=== Name of directory for mesh file ==='
      CHARACTER(len = rec_length) :: file_name         = '=== Name of mesh file ==='
      CHARACTER(len = rec_length) :: if_mesh_formatted = '=== Is the mesh formatted? (True/False) ==='
      CHARACTER(len = rec_length) :: if_read_partition = '=== Do we read metis partition? (true/false) ==='
   END TYPE argument_mesh_data_type
   !===default value in simulation
   TYPE mesh_data_type
      CHARACTER(len = rec_length)    :: directory         = '.'
      CHARACTER(len = rec_length)    :: file_name         = 'mesh_name'
      LOGICAL                        :: if_mesh_formatted = .TRUE.
      LOGICAL                        :: if_read_partition = .FALSE.
      LOGICAL :: init_once = .FALSE.
   CONTAINS
      PROCEDURE, PUBLIC              :: read => read_mesh_data
      PROCEDURE, PUBLIC              :: init => init_mesh_data
   END TYPE mesh_data_type

   TYPE argument_mesh_info_type
      CHARACTER(len = rec_length) :: nb_dom            = '=== Number of subdomains in the mesh ==='
      CHARACTER(len = rec_length) :: list_dom          = '=== List of subdomains in the mesh ==='
      CHARACTER(len = rec_length) :: type_fe           = '=== Type of finite element ==='
      CHARACTER(len = rec_length) :: nb_refinement     = '=== Number of refinement steps ==='
      CHARACTER(len = rec_length) :: refinement_order  = '=== Order of refinement? (1d) ==='
      CHARACTER(len = rec_length) :: nb_bords          = '=== How many pieces of periodic boundary? ==='
      CHARACTER(len = rec_length) :: list_periodic     = '=== Indices of periodic boundaries and corresponding vectors ==='
   END TYPE argument_mesh_info_type
   !===default value in simulation
   TYPE mesh_info_type
      INTEGER                        :: nb_dom            = 1
      INTEGER, DIMENSION(:), ALLOCATABLE :: list_dom
      INTEGER                        :: type_fe           = 1
      INTEGER                        :: nb_refinement     = 0
      INTEGER                        :: refinement_order  = 0
      INTEGER                        :: nb_bords          = 0
      INTEGER,        DIMENSION(:, :), ALLOCATABLE :: list_periodic
      REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: vect_e
   CONTAINS
      PROCEDURE, PUBLIC              :: read => read_mesh_info
      PROCEDURE, PUBLIC              :: init => init_mesh_info
      PROCEDURE, PUBLIC              :: copy => copy_mesh_info
   END TYPE mesh_info_type

CONTAINS

   SUBROUTINE init_mesh_data(this)
      CLASS(mesh_data_type), INTENT(INOUT) :: this
      IF (.NOT. this%init_once) CALL this%READ
      this%init_once = .TRUE.
   END SUBROUTINE init_mesh_data

   SUBROUTINE read_mesh_data(this)
      USE space_dim
      USE read_inputs_module
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

      !================
      !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
      !================
      CALL finalize_rewrite_data

   END SUBROUTINE read_mesh_data

   SUBROUTINE init_mesh_info(this, opt_name)
      CLASS(mesh_info_type), INTENT(INOUT) :: this
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_name
      CALL this%READ(opt_name=opt_name)
   END SUBROUTINE init_mesh_info

   SUBROUTINE read_mesh_info(this, opt_name)
      USE space_dim
      USE read_inputs_module
      IMPLICIT NONE

      ! CHARACTER(LEN=rec_length), PARAMETER   :: section_name='MESH SETUP'
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_name

      CLASS(mesh_info_type), INTENT(INOUT) :: this
      TYPE(argument_mesh_info_type)        :: argument_data

      !=== Reading all data file
      IF (PRESENT(opt_name)) THEN
         ! write(*,*) section_name//" FOR "//opt_name
         CALL read_data_init_list("MESH SETUP FOR "//opt_name)
      ELSE
         ! write(*,*) section_name
         CALL read_data_init_list("MESH SETUP")
      END IF
      !================
      !=== We now find the relevant information for the mesh
      !================

      !=== type of finite element
      CALL read_data(argument_data%type_fe, this%type_fe, opt_name=opt_name)!, opt_add=(k_dim==2))

      !=== number of refinement steps
      CALL read_data(argument_data%nb_refinement, this%nb_refinement, opt_name=opt_name, opt_add=(k_dim==2))

      !=== order refinement
      CALL read_data(argument_data%refinement_order, this%refinement_order, opt_name=opt_name, opt_add=(k_dim==1))

      !=== number of subdomains in the mesh
      CALL read_data(argument_data%nb_dom, this%nb_dom, opt_name=opt_name)

      !=== list of subdomains in the mesh (special treatment since it is an array)
      ALLOCATE(this%list_dom(this%nb_dom))
      this%list_dom(1) = 1
      CALL read_data(argument_data%list_dom, this%list_dom, opt_name=opt_name)

      !=== Number of periodic boundaries
      CALL read_data(argument_data%nb_bords, this%nb_bords, opt_name=opt_name)

      !=== List of periodic boundaries (has its specific subroutine, see read_inputs_module.F90 module)
      ALLOCATE(this%list_periodic(2, this%nb_bords))
      ALLOCATE(this%vect_e(k_dim, this%nb_bords))
      CALL read_periodic_data(argument_data%list_periodic, this%nb_bords, &
                              this%list_periodic, this%vect_e, opt_name=opt_name)
      !================
      !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
      !================
      CALL finalize_rewrite_data

   END SUBROUTINE read_mesh_info

   SUBROUTINE copy_mesh_info(this, info_bis)
      IMPLICIT NONE
      CLASS(mesh_info_type), INTENT(OUT)  :: this
      TYPE(mesh_info_type),  INTENT(IN) :: info_bis

      this%nb_dom           = info_bis%nb_dom
      this%list_dom         = info_bis%list_dom
      this%type_fe          = info_bis%type_fe
      this%nb_refinement    = info_bis%nb_refinement
      this%refinement_order = info_bis%refinement_order
      this%nb_bords         = info_bis%nb_bords
      this%list_periodic    = info_bis%list_periodic
      this%vect_e           = info_bis%vect_e
   END SUBROUTINE copy_mesh_info

END MODULE mesh_data_module


MODULE mesh_parameters
   USE mesh_data_module
   TYPE(mesh_data_type) :: mesh_data_info
END MODULE mesh_parameters