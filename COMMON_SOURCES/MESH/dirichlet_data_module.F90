MODULE dirichlet_data_module
   IMPLICIT NONE
   TYPE dirichlet_data_type
      INTEGER :: nb_dirichlet
      INTEGER, DIMENSION(:), POINTER :: list_dirichlet
   CONTAINS
      PROCEDURE, PUBLIC :: init
   END TYPE dirichlet_data_type
CONTAINS
   SUBROUTINE init(a)
      CLASS(dirichlet_data_type), INTENT(INOUT) :: a
      a%nb_dirichlet = 0
   END SUBROUTINE init
END MODULE dirichlet_data_module

MODULE input_dirichlet_data
   USE dirichlet_data_module
   IMPLICIT NONE
   PUBLIC :: read_dirichlet_data
   TYPE(dirichlet_data_type), PUBLIC :: dirichlet_data
   PRIVATE
CONTAINS
   SUBROUTINE read_dirichlet_data(data_fichier)
      USE character_strings
      USE space_dim
      IMPLICIT NONE
      INTEGER, PARAMETER :: in_unit = 21
      INTEGER :: k
      CHARACTER(len = *), INTENT(IN) :: data_fichier
      CHARACTER(LEN = 100) :: argument
      LOGICAL :: test
      !===Initialize data to zero and false by default
      CALL dirichlet_data%init
      OPEN(UNIT = in_unit, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')

      CALL find_string(21, '===How many pieces of dirichlet boundary?===', test)
      IF (test) THEN
         READ (21, *) dirichlet_data%nb_dirichlet
      ELSE
         dirichlet_data%nb_dirichlet_pairs = 0
      END IF
      ALLOCATE(dirichlet_data%list_dirichlet(dirichlet_data%nb_dirichlet))

      IF (dirichlet_data%nb_dirichlet_pairs > 0) THEN
         CALL read_until(21, '===List of dirichlet boundaries===')
         READ(21, *) dirichlet_data%list_dirichlet
      END IF

      CLOSE(in_unit)
   END SUBROUTINE read_dirichlet_data
END MODULE input_dirichlet_data
