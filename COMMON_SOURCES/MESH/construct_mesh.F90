MODULE construct_mesh
#include "petsc/finclude/petsc.h"
   USE def_type_mesh
   USE space_dim
   USE st_matrix
   USE petsc
   USE input_mesh_data
   USE mesh_data_module
   PUBLIC :: get_mesh
   PRIVATE
CONTAINS
   SUBROUTINE get_mesh(communicator, mesh, opt_fe, opt_edge_stab)
      USE mesh_1d
      USE mesh_distribution_1d
      IMPLICIT NONE
      LOGICAL, OPTIONAL :: opt_edge_stab
      INTEGER, OPTIONAL :: opt_fe
      LOGICAL :: edge_stab
      TYPE(mesh_type) :: mesh_glob, mesh
      MPI_Comm       :: communicator

      CALL read_mesh_data('data')

      IF (.NOT.PRESENT(opt_edge_stab)) THEN
         edge_stab = .FALSE.
      ELSE
         edge_stab = opt_edge_stab
      END IF

      SELECT CASE(k_dim)
      CASE(2)
         write(*, *) ' BUG in construct_mesh, k_dim = 2 not implemented'
         STOP
      CASE(1)
         CALL load_mesh_1d(mesh_glob)

         CALL extract_mesh_1d(communicator, mesh_glob, mesh)

         CALL GAUSS_POINT_1d(mesh)

      CASE DEFAULT
         write(*, *) ' BUG in construct_mesh, k_dim not correct'
         STOP
      END SELECT

   END SUBROUTINE get_mesh

END MODULE  construct_mesh