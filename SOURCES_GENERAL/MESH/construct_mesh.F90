MODULE construct_mesh
#include "petsc/finclude/petsc.h"
   USE def_type_mesh
   USE space_dim
   USE st_matrix
   USE petsc
   USE input_data
   PUBLIC :: get_mesh
   PRIVATE
CONTAINS
   SUBROUTINE get_mesh(communicator, mesh, LA, js_d_loc, fe_type, opt_edge_stab)
      USE mesh_1d
      USE mesh_distribution_1d
      IMPLICIT NONE
      LOGICAL, OPTIONAL :: opt_edge_stab
      LOGICAL :: edge_stab
      INTEGER :: fe_type
      INTEGER, POINTER, DIMENSION(:) :: js_d_loc
      TYPE(periodic_type) :: opt_per
      TYPE(mesh_type) :: mesh_glob, mesh
      TYPE(petsc_csr_LA) :: LA
      MPI_Comm       :: communicator

      IF (.NOT.PRESENT(opt_edge_stab)) THEN
         edge_stab = .FALSE.
      ELSE
         edge_stab = opt_edge_stab
      END IF
      SELECT CASE(k_dim)
      CASE(2)
         !CALL mesh_2d(mesh)
         write(*, *) ' BUG in construct_mesh, k_dim = 2 not implemented'
         STOP
      CASE(1)
         CALL load_mesh_1d(mesh_glob)
         CALL extract_mesh_1d(communicator, mesh_glob, mesh)
      CASE DEFAULT
         write(*, *) ' BUG in construct_mesh, k_dim not correct'
         STOP
      END SELECT

      CALL prep_periodic_scal(inputs%my_periodic, mesh, opt_per)
      CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA, opt_per = opt_per)
      CALL dirichlet_nodes_parallel(mesh, inputs%Dir_list, js_d_loc)
   END SUBROUTINE get_mesh
END MODULE  construct_mesh