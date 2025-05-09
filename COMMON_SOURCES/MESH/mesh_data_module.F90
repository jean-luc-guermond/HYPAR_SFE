MODULE mesh_data_module
   IMPLICIT NONE
   TYPE mesh_data_type
      CHARACTER(len = 200) :: directory
      CHARACTER(len = 200) :: file_name
      LOGICAL :: if_mesh_formatted, if_HCT, if_read_partition
      CHARACTER(len = 20) :: mesh_structure
      INTEGER :: nb_dom
      INTEGER, DIMENSION(:), POINTER :: list_dom
      INTEGER :: type_fe, nb_refinement
   CONTAINS
      PROCEDURE, PUBLIC :: init
      PROCEDURE, PUBLIC :: read => read_mesh_data
   END TYPE mesh_data_type
CONTAINS
   SUBROUTINE init(a)
      CLASS(mesh_data_type), INTENT(INOUT) :: a
      !===Logicals
      a%if_mesh_formatted = .FALSE.
      a%if_HCT = .false.
      !===Characters
      a%directory = '.'
      a%file_name = 'gnu'
      !===Integers
      a%nb_dom = -1
      a%type_fe = -1
      a%nb_refinement = 0
      a%type_fe = 2
   END SUBROUTINE init

   SUBROUTINE read_mesh_data(this)
      USE character_strings
      USE petsc
      IMPLICIT NONE
      CLASS(mesh_data_type), INTENT(INOUT) :: this
      INTEGER, PARAMETER :: in_unit = 21
      CHARACTER(LEN = 100) :: argument
      LOGICAL :: test
      INTEGER :: rank, ierr

      CALL MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

      !===Initialize data to zero and false by default
      CALL this%init
      OPEN(UNIT = in_unit, FILE = 'data', FORM = 'formatted', STATUS = 'unknown')

      argument = '===Name of directory for mesh file==='
      CALL find_string(in_unit, argument, test)
      IF (test) THEN
         READ (in_unit, *) this%directory
      ELSE
         CALL default_data(rank, in_unit, argument, '.')
         this%directory = '.'
      END IF

      argument = "===Name of mesh file==="
      CALL find_string(in_unit, argument, test)
      IF (test) THEN
         READ (in_unit, *) this%file_name
      ELSE
         CALL default_data(rank, in_unit, argument, 'mesh_name')
         this%file_name = 'mesh_name'
         IF (rank == 0) WRITE(*, *) "No mesh_name specified." ; STOP
      END IF

      argument = "===Is the mesh formatted? (True/False)==="
      CALL find_string(in_unit, argument, test)
      IF (test) THEN
         READ (in_unit, *) this%if_mesh_formatted
      ELSE
         CALL default_data(rank, in_unit, argument, '.t.')
         this%if_mesh_formatted = .true.
      END IF

      argument = '===Do we read metis partition? (true/false)'
      CALL find_string(in_unit, argument, test)
      IF (test) THEN
         READ (in_unit, *) this%if_read_partition
      ELSE
         CALL default_data(rank, in_unit, argument, '.f.')
         this%if_read_partition = .false.
      END IF

      argument = '===Number of subdomains in the mesh==='
      CALL find_string(in_unit, argument, test)
      IF (test) THEN
         READ (in_unit, *) this%nb_dom
      ELSE
         CALL default_data(rank, in_unit, argument, '1')
         this%nb_dom = 1
      END IF

      argument = '===List of subdomain in the mesh==='
      ALLOCATE(this%list_dom(this%nb_dom))
      CALL find_string(in_unit, argument, test)
      IF (test) THEN
         READ (in_unit, *) this%list_dom
      ELSE
         CALL default_data(rank, in_unit, argument, '1')
         this%list_dom(1) = 1
      END IF

      argument = '===Number of refinement steps==='
      ALLOCATE(this%list_dom(this%nb_dom))
      CALL find_string(in_unit, argument, test)
      IF (test) THEN
         READ (in_unit, *) this%nb_refinement
      ELSE
         CALL default_data(rank, in_unit, argument, '0')
         this%nb_refinement = 0
      END IF

      argument = "===HCT mesh ?==="
      CALL find_string(in_unit, argument, test)
      IF (test) THEN
         READ (in_unit, *) this%if_HCT
         IF (this%if_HCT) write(*, *) "HCT mesh not inmplemented yet"
      ELSE
         this%if_HCT = .false.
      END IF

      argument = '===Type of finite element==='
      CALL find_string(in_unit, argument, test)
      IF (test) THEN
         READ (in_unit, *) this%type_fe
      ELSE
         CALL default_data(rank, in_unit, argument, '1')
         this%type_fe = 1
      END IF

      CLOSE(in_unit)
   END SUBROUTINE read_mesh_data
END MODULE mesh_data_module
