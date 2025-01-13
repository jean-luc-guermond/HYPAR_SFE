MODULE mesh_2d
   USE def_type_mesh
   USE input_data
   PUBLIC :: mesh_2d
   PRIVATE

CONTAINS
   SUBROUTINE mesh_2d(mesh)
#include "petsc/finclude/petsc.h"
      USE petsc
      USE input_data
      USE prep_maill
      USE mod_gauss_points_2d
      USE metis_reorder_elements
      USE st_matrix
      USE dir_nodes_petsc
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_glob
      INTEGER :: dom_np
      REAL(KIND = 8) :: t1
      PetscErrorCode :: ierr
      PetscMPIInt    :: nb_proc, rank

      !===Number of procs
      CALL MPI_Comm_size(PETSC_COMM_WORLD, nb_proc, ierr)
      CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

      !===Prepare the grid
      t1 = user_time()
      CALL load_mesh_free_format_iso(inputs%directory, inputs%file_name, inputs%list_dom, &
           inputs%type_fe, mesh_glob, inputs%if_mesh_formatted, edge_stab = .FALSE.)
      IF (rank==0) WRITE(*, *) ' time load_mesh_free_format_iso', (user_time() - t1)
      !CALL incr_vrtx_indx_enumeration(mesh_glob,inputs%type_fe)

      !===Metis reorganizes the mesh
      t1 = user_time()
      CALL reorder_mesh(PETSC_COMM_WORLD, nb_proc, mesh_glob, uu_mesh)
      IF (rank==0) WRITE(*, *) ' reorder_mesh', (user_time() - t1)

      !===Deallocate global mesh
      CALL free_mesh_after(mesh_glob)

      !===Create Sparsity pattern for the matrix (structure)
      t1 = user_time()
      CALL st_aij_csr_glob_block_with_extra_layer(PETSC_COMM_WORLD, 1, uu_mesh, LA_u)
      IF (rank==0) WRITE(*, *) ' st_aij_csr_glob_block', (user_time() - t1)

      dom_np = SIZE(LA_u%ia) - 1
      max_nnz_per_row = MAXVAL(LA_u%ia(1:dom_np) - LA_u%ia(0:dom_np - 1))

      !===Start Gauss points generation
      uu_mesh%edge_stab = .FALSE.
      CALL gauss_points_2d(uu_mesh, inputs%type_fe)

      !===Boundary conditions on local domain
      CALL dirichlet_nodes_parallel(uu_mesh, inputs%Dir_list, js_d_loc_u)

   END SUBROUTINE mesh_2d
END MODULE mesh_2d


