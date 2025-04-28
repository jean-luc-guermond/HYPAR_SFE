MODULE periodic_data_module
   IMPLICIT NONE
   TYPE periodic_data_type
      INTEGER :: nb_periodic_pairs
      INTEGER, DIMENSION(:, :), POINTER :: list_periodic
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: vect_e
   CONTAINS
      PROCEDURE, PUBLIC :: init
   END TYPE periodic_data_type
CONTAINS
   SUBROUTINE init(a)
      CLASS(periodic_data_type), INTENT(INOUT) :: a
      a%nb_periodic_pairs = 0
   END SUBROUTINE init
END MODULE periodic_data_module

MODULE input_periodic_data
   USE periodic_data_module
   IMPLICIT NONE
   PUBLIC :: read_periodic_data
   TYPE(periodic_data_type), PUBLIC :: periodic_data
   PRIVATE
CONTAINS
   SUBROUTINE read_periodic_data(data_fichier)
      USE character_strings
      IMPLICIT NONE
      INTEGER, PARAMETER :: in_unit = 21
      INTEGER :: k
      CHARACTER(len = *), INTENT(IN) :: data_fichier
      CHARACTER(LEN = 100) :: argument
      LOGICAL :: test
      !===Initialize data to zero and false by default
      CALL periodic_data%init

      CALL find_string(21, '===How many pieces of periodic boundary?===', test)
      IF (test) THEN
         READ (21, *) periodic_data%nb_periodic_pairs
      ELSE
         periodic_data%nb_periodic_pairs = 0
      END IF

      ALLOCATE(periodic_data%list_periodic(2, periodic_data%nb_periodic_pairs))
      ALLOCATE(periodic_data%vect_e(2, periodic_data%nb_periodic_pairs))

      IF (periodic_data%nb_periodic_pairs > 0) THEN
         CALL read_until(21, '===Indices of periodic boundaries and corresponding vectors===')
         DO k = 1, periodic_data%nb_periodic_pairs
            READ(21, *) periodic_data%list_periodic(:, k), periodic_data%vect_e(:, k)
         END DO
      END IF

      CLOSE(in_unit)
   END SUBROUTINE read_periodic_data
END MODULE input_periodic_data
