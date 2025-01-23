MODULE mesh_data_module
  IMPLICIT NONE
  TYPE mesh_data_type
     CHARACTER(len=200)             :: directory
     CHARACTER(len=200)             :: file_name
     LOGICAL                        :: if_mesh_formatted
     CHARACTER(len=20)              :: mesh_structure
     INTEGER                        :: nb_dom
     INTEGER, DIMENSION(:), POINTER :: list_dom
     INTEGER                        :: type_fe
   CONTAINS
     PROCEDURE, PUBLIC              :: init
  END TYPE mesh_data_type
CONTAINS
  SUBROUTINE init(a)
    CLASS(mesh_data_type), INTENT(INOUT) :: a
    !===Logicals
    a%if_mesh_formatted = .FALSE.
    !===Characters
    a%directory='.'
    a%file_name='gnu'
    !===Integers
    a%nb_dom=-1
    a%type_fe=-1
  END SUBROUTINE init
END MODULE mesh_data_module

MODULE input_mesh_data
  USE mesh_data_module
  IMPLICIT NONE
  PUBLIC :: read_mesh_data
  TYPE(mesh_data_type), PUBLIC  :: mesh_data
  PRIVATE
CONTAINS
  SUBROUTINE read_mesh_data(data_fichier)
    USE character_strings
    IMPLICIT NONE
    INTEGER, PARAMETER           :: in_unit=21
    CHARACTER(len=*), INTENT(IN) :: data_fichier
    CHARACTER(LEN=100)           :: argument
    LOGICAL :: okay
    !===Initialize data to zero and false by default
    CALL mesh_data%init

    OPEN(UNIT = in_unit, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(in_unit, "===Name of directory for mesh file===")
    READ (in_unit,*) mesh_data%directory
    CALL read_until(in_unit, "===Name of mesh file===")
    READ (in_unit,*) mesh_data%file_name
    CALL read_until(in_unit, "===Is the mesh formatted? (True/False)===")
    READ (in_unit,*) mesh_data%if_mesh_formatted
    CALL read_until(in_unit, '===Number of subdomains in the mesh===')
    READ(21,*) mesh_data%nb_dom
    ALLOCATE(mesh_data%list_dom(mesh_data%nb_dom))
    CALL read_until(21, '===List of subdomain in the mesh===')
    READ(21,*) mesh_data%list_dom


    CLOSE(in_unit)
  END SUBROUTINE read_mesh_data
END MODULE input_mesh_data
