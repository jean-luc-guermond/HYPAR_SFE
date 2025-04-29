MODULE construct_mesh
#include "petsc/finclude/petsc.h"
   USE def_type_mesh
   USE space_dim
   USE st_matrix
   USE petsc
   USE input_mesh_data
   USE input_periodic_data
   USE mesh_data_module
   USE mesh_tools
   USE periodic_data_module
   PUBLIC :: get_mesh
   PRIVATE
CONTAINS
   SUBROUTINE get_mesh(communicator, mesh, opt_fe, opt_edge_stab, opt_per)
      USE mesh_1d
      USE mesh_distribution_1d
      USE load_mesh_2d
      USE refine_mesh
      USE two_dim_metis_distribution
      USE gauss_points_2d
      IMPLICIT NONE
      LOGICAL, OPTIONAL :: opt_edge_stab, opt_per
      INTEGER, OPTIONAL :: opt_fe
      INTEGER, DIMENSION(1) :: list_dom = 1
      INTEGER, DIMENSION(0) :: list_inter
      INTEGER :: n, nb_proc, ierr
      LOGICAL :: edge_stab, per_bool
      TYPE(mesh_type) :: mesh_glob, mesh, mesh_r
      MPI_Comm       :: communicator

      CALL read_mesh_data('data')
      CALL MPI_Comm_SIZE(communicator, nb_proc, ierr)

      IF (PRESENT(opt_per)) THEN
         per_bool = opt_per
      ELSE
         per_bool = .false.
      END IF

      IF (per_bool) THEN
         CALL read_periodic_data('data')
      END IF

      IF (.NOT.PRESENT(opt_edge_stab)) THEN
         edge_stab = .FALSE.
      ELSE
         edge_stab = opt_edge_stab
      END IF

      SELECT CASE(k_dim)
      CASE(2)
         IF (per_bool) THEN
            !===load and re order mesh
            CALL load_dg_mesh_free_format(mesh_data%directory, mesh_data%file_name, &
                 list_dom, list_inter, mesh_glob, mesh_data%if_mesh_formatted)
            write(*, *) 'load done'
            CALL reorder_mesh(PETSC_COMM_WORLD, nb_proc, mesh_glob, mesh)
            write(*, *) 'reorder done'
            CALL free_mesh(mesh_glob)

            !===mesh refinements
            DO n = 1, mesh_data%nb_refinement
               !===Create refined mesh
               CALL refinement_iso_grid_distributed(mesh)
               write(*, *) 'refinement done', n
            END DO

            !===special meshes
            !      IF(mesh_data%if_HCT) THEN
            !         CALL HCT_iso_grid_distributed(mesh_p1, HCT_mesh_p1)
            !         CALL deallocate_mesh(mesh_p1)
            !         CALL copy_mesh(HCT_mesh_p1, mesh_p1)
            !         CALL deallocate_mesh(HCT_mesh_p1)
            !      END IF

            !===create finite elements polynome on mesh
            CALL create_iso_grid_distributed(mesh, mesh_r, mesh_data%type_fe)
            write(*,*) 'iso done'
            CALL free_mesh(mesh)
            CALL copy_mesh(mesh_r, mesh)
            CALL free_mesh(mesh_r)

            !===gauss points on mesh
            CALL create_gauss_points_2d(mesh, mesh_data%type_fe)
            write(*,*) 'gauss points done'

         END IF

      CASE(1)
         CALL load_mesh_1d(mesh_glob)

         CALL extract_mesh_1d(communicator, mesh_glob, mesh, opt_per)
         CALL free_mesh(mesh_glob)
         CALL GAUSS_POINT_1d(mesh)

      CASE DEFAULT
         write(*, *) ' BUG in construct_mesh, k_dim not correct'
         STOP
      END SELECT

   END SUBROUTINE get_mesh

END MODULE  construct_mesh