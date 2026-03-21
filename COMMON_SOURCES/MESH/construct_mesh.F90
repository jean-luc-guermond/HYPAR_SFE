MODULE construct_mesh
#include "petsc/finclude/petsc.h"
   USE def_type_mesh
   USE space_dim
   USE st_matrix
   USE petsc
   USE mesh_data_module
   USE mesh_tools
   USE periodic_data_module
   PUBLIC :: get_mesh
   PRIVATE
CONTAINS
   SUBROUTINE get_mesh(communicator, mesh, opt_fe, opt_edge_stab, opt_pers)
      USE mesh_1d
      USE mesh_distribution_1d
      USE load_mesh_2d
      USE refine_mesh
      USE two_dim_metis_distribution
      USE gauss_points_2d
      USE Dir_nodes
      IMPLICIT NONE
      LOGICAL, OPTIONAL :: opt_edge_stab
      INTEGER, OPTIONAL :: opt_fe
      TYPE(periodic_type), DIMENSION(:), OPTIONAL :: opt_pers
      TYPE(mesh_data_type) :: mesh_data
      INTEGER, DIMENSION(1) :: list_dom = 1
      INTEGER, DIMENSION(0) :: list_inter
      INTEGER, DIMENSION(:), ALLOCATABLE :: part
      INTEGER :: n, nb_proc, ierr, rank
      LOGICAL :: edge_stab
      CHARACTER(LEN = 100) :: mesh_part_name
      TYPE(mesh_type) :: mesh_glob, mesh, mesh_r
      MPI_Comm       :: communicator

      CALL mesh_data%read()
      CALL MPI_Comm_SIZE(communicator, nb_proc, ierr)
      CALL MPI_Comm_rank(communicator, rank, ierr)

      IF (.NOT.PRESENT(opt_edge_stab)) THEN
         edge_stab = .FALSE.
      ELSE
         edge_stab = opt_edge_stab
      END IF

      !=== FIXME mesh%rank to be transferred through refinement_iso_grid_distributed,
      !=== create_iso_grid_distributed, copy_mesh
      mesh_glob%rank = -1

      SELECT CASE(k_dim)
      CASE(2)
         !===load and re order mesh
         CALL load_dg_mesh_free_format(mesh_data%directory, mesh_data%file_name, &
              list_dom, list_inter, mesh_glob, mesh_data%if_mesh_formatted)
         ALLOCATE(part(mesh_glob%me))

         mesh_part_name = 'mesh_part.' // TRIM(ADJUSTL(mesh_data%file_name))
         IF (mesh_data%if_read_partition) THEN
            IF (rank == 0) WRITE(*, *) 'read partition'
            OPEN(UNIT = 51, FILE = mesh_part_name, STATUS = 'unknown', FORM = 'formatted')
            READ(51, *) part
            CLOSE(51)
         ELSE
            IF (rank == 0) WRITE(*, *) 'create partition'
            CALL part_mesh(nb_proc, mesh_glob, list_inter, part, opt_pers)
            IF (rank==0) THEN
               OPEN(UNIT = 51, FILE = mesh_part_name, STATUS = 'replace', FORM = 'formatted')
               WRITE(51, *) part
               CLOSE(51)
            END IF
         END IF

         CALL extract_mesh(communicator, nb_proc, mesh_glob, part, list_dom, mesh)
         CALL free_mesh(mesh_glob)
         DEALLOCATE(part)
         !===mesh refinements
         DO n = 1, mesh_data%nb_refinement
            !===Create refined mesh
            CALL refinement_iso_grid_distributed(mesh)
            IF(rank == 0) write(*, *) 'refinement done', n
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
         CALL free_mesh(mesh)
         CALL copy_mesh(mesh_r, mesh)
         CALL free_mesh(mesh_r)

         !===(JLG) Added March 21 2026 
         ALLOCATE(mesh%iis(SIZE(mesh%jjs,1),mesh%mes))
         CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
         CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
         mesh%nps = SIZE(mesh%j_s)
         !===END (JLG) Added March 21 2026 

         mesh%rank = rank  !=== petsc convention
         mesh%edge_stab = .false.
         !===gauss points on mesh
         CALL create_gauss_points_2d(mesh, mesh_data%type_fe)

      CASE(1)
         CALL load_mesh_1d(mesh_data%directory, mesh_data%file_name, mesh_glob, mesh_data%if_mesh_formatted)
         CALL extract_mesh_1d(communicator, mesh_glob, mesh, opt_pers)
         CALL free_mesh(mesh_glob)
         mesh%edge_stab = .false.
         CALL GAUSS_POINT_1d(mesh)
         mesh%rank = rank

      CASE DEFAULT
         IF(rank == 0) write(*, *) ' BUG in construct_mesh, k_dim not correct'
         STOP
      END SELECT

      mesh%edge_stab = .false. !TODO remove edge_stab

   END SUBROUTINE get_mesh

END MODULE  construct_mesh
