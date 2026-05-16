MODULE construct_mesh
#include "petsc/finclude/petsc.h"
  USE def_type_mesh
  USE st_matrix
  USE petsc
  USE mesh_data_module
  USE mesh_tools
  use gauss_points_1d, ONLY: create_gauss_points_1d
  PUBLIC :: get_mesh
  PRIVATE
CONTAINS
  SUBROUTINE get_mesh(communicator, mesh, opt_fe, opt_edge_stab)
    USE mesh_1d
    USE mesh_distribution_1d
    USE load_mesh_2d
    USE refine_mesh
    USE two_dim_metis_distribution
    USE gauss_points_2d
    USE Dir_nodes
    USE space_dim
    USE mesh_parameters
    USE my_util, ONLY : error_petsc, to_str

    IMPLICIT NONE
    LOGICAL, OPTIONAL :: opt_edge_stab
    INTEGER, OPTIONAL :: opt_fe
    INTEGER, DIMENSION(1) :: list_dom = 1
    INTEGER, DIMENSION(0) :: list_inter
    INTEGER, DIMENSION(:), ALLOCATABLE :: part
    INTEGER :: n, nb_proc, ierr, rank
    LOGICAL :: edge_stab
    CHARACTER(LEN = 100) :: mesh_part_name
    TYPE(mesh_type) :: mesh_glob, mesh, mesh_r
    MPI_Comm       :: communicator

    CALL mesh_data_info%init
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
       CALL load_dg_mesh_free_format(mesh_data_info%directory, mesh_data_info%file_name, &
            list_dom, list_inter, mesh_glob, mesh_data_info%if_mesh_formatted)
       ALLOCATE(part(mesh_glob%me))

       mesh_part_name = 'mesh_part.' // TRIM(ADJUSTL(mesh_data_info%file_name))
       IF (mesh_data_info%if_read_partition) THEN
          IF (rank == 0) WRITE(*, *) 'read partition'
          OPEN(UNIT = 51, FILE = mesh_part_name, STATUS = 'unknown', FORM = 'formatted')
          READ(51, *) part
          CLOSE(51)
       ELSE
          IF (rank == 0) WRITE(*, *) 'create partition'
          CALL part_mesh(nb_proc, mesh_glob, list_inter, part)
          IF (rank==0) THEN
             OPEN(UNIT = 51, FILE = mesh_part_name, STATUS = 'replace', FORM = 'formatted')
             WRITE(51, *) part
             CLOSE(51)
          END IF
       END IF

       CALL extract_mesh(communicator, mesh_glob, part, list_dom, mesh)
       CALL free_mesh(mesh_glob)
       DEALLOCATE(part)
       !===mesh refinements
       DO n = 1, mesh_data_info%nb_refinement
          !===Create refined mesh
          CALL refinement_iso_grid_distributed(mesh)
          IF(rank == 0) write(*, *) 'refinement done', n
       END DO

       !===special meshes
       !      IF(mesh_data_info%if_HCT) THEN
       !         CALL HCT_iso_grid_distributed(mesh_p1, HCT_mesh_p1)
       !         CALL deallocate_mesh(mesh_p1)
       !         CALL copy_mesh(HCT_mesh_p1, mesh_p1)
       !         CALL deallocate_mesh(HCT_mesh_p1)
       !      END IF

       !===create finite elements polynome on mesh
       CALL create_iso_grid_distributed(mesh, mesh_r, mesh_data_info%type_fe)
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
       CALL create_gauss_points_2d(mesh, mesh_data_info%type_fe)

    CASE(1)
       CALL load_mesh_1d(mesh_data_info%directory, mesh_data_info%file_name, mesh_glob, mesh_data_info%if_mesh_formatted)
       CALL extract_mesh_1d(communicator, mesh_glob, mesh)
       CALL free_mesh(mesh_glob)
       mesh%edge_stab = .false.
       mesh%rank = rank

       != new ==> create Pk mesh
       CALL create_Pk_mesh_1D(communicator, mesh, mesh_r, mesh_data_info%type_fe)
       CALL free_mesh(mesh)
       CALL copy_mesh(mesh_r, mesh)
       CALL free_mesh(mesh_r)

       != new ==> create Pk mesh
       CALL create_gauss_points_1d(mesh, mesh_data_info%type_fe)

       
    CASE DEFAULT
       CALL error_petsc("BUG in construct_mesh, incorrect k_dim="//to_str(k_dim)//&
                         "  should be 1 or 2")
    END SELECT

    mesh%edge_stab = .false. !FIXME remove edge_stab

    mesh%per%nb_bords = mesh_data_info%nb_bords
    ALLOCATE(mesh%per%list_periodic(SIZE(mesh_data_info%list_periodic,1),SIZE(mesh_data_info%list_periodic,2)))
    ALLOCATE(mesh%per%vect_e(SIZE(mesh_data_info%vect_e,1),SIZE(mesh_data_info%vect_e,2)))
    IF (mesh%per%nb_bords/=0) THEN 
       mesh%per%list_periodic = mesh_data_info%list_periodic
       mesh%per%vect_e = mesh_data_info%vect_e
       CALL prep_periodic_scal(mesh%per,mesh)
    END IF

  END SUBROUTINE get_mesh


  SUBROUTINE prep_periodic_scal(periodic, mesh)
      !=========================================
      USE character_strings
      IMPLICIT NONE
      TYPE(periodic_type) :: periodic
      TYPE(mesh_type) :: mesh
      INTEGER, DIMENSION(:), POINTER :: list_loc, perlist_loc, list_dom, perlist_dom
      INTEGER :: n, side1, side2, n_b, nx, i
      REAL(KIND = 8), DIMENSION(:), POINTER :: e

      IF (mesh%np == 0) THEN
         !WRITE(*, *) 'no mesh on this proc'
         RETURN
      END IF

      ALLOCATE(e(SIZE(periodic%vect_e, 1)))

      IF (periodic%nb_bords .GT. 20) THEN
         WRITE(*, *) 'PREP_MESH_PERIODIC: too many periodic pieces'
         STOP
      END IF

      DO n = 1, periodic%nb_bords

         side1 = periodic%list_periodic(1, n)
         side2 = periodic%list_periodic(2, n)
         e = periodic%vect_e(:, n)

         CALL list_periodic(mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
              list_loc, perlist_loc)

         n_b = SIZE(list_loc)
         ALLOCATE(list_dom(n_b), perlist_dom(n_b))
         nx = 0
         DO i = 1, n_b
            IF (MAX(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
               nx = nx + 1
               list_dom(nx) = list_loc(i)
               perlist_dom(nx) = perlist_loc(i)
            ELSE IF (MIN(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
               WRITE(*, *) 'BUG in prep_periodic_scal, one of the boundary point is not attributed the same processor.',&
                           'rank = ', mesh%rank, list_loc(i), perlist_loc(i), mesh%dom_np
               WRITE(*,*) "BUG IS OCURRING AT COORDINATES LIST = ", mesh%rr(:, list_loc(i)), "PERLIST = ", mesh%rr(:, perlist_loc(i))
               STOP
            END IF
         END DO
         IF (n_b /= nx) WRITE(*, *) 'WARNING on bord n=',n,', I have removed', n_b - nx, &
                                   ' periodic pairs in prep_periodic_scal for rank ', mesh%rank
         n_b = nx

         ALLOCATE (periodic%list(n)%DIL(n_b), periodic%perlist(n)%DIL(n_b))
         periodic%list(n)%DIL = list_dom(1:n_b)
         periodic%perlist(n)%DIL = perlist_dom(1:n_b)

         DEALLOCATE(list_loc, perlist_loc, list_dom, perlist_dom)
      END DO

      DEALLOCATE(e)

    END SUBROUTINE prep_periodic_scal

   SUBROUTINE prep_periodic_bloc(periodic, mesh, nb_bloc)
      !=========================================
      USE character_strings
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      TYPE(periodic_type) :: periodic
      INTEGER, INTENT(IN) :: nb_bloc
      INTEGER, DIMENSION(:), POINTER :: list_loc, perlist_loc, list_dom, perlist_dom
      INTEGER :: n, side1, side2, nsize, n_b
      INTEGER :: k, k_deb, k_fin, nx, i
      REAL(KIND = 8), DIMENSION(2) :: e

      WRITE (*, *) 'Loading periodic-data file ...'

      IF (mesh%np == 0) THEN
         WRITE(*, *) 'no mesh on this proc'
         RETURN
      END IF

      IF (periodic%nb_bords .GT. 20) THEN
         WRITE(*, *) 'PREP_MESH_PERIODIC: trop de bords periodiques'
         STOP
      END IF

      DO n = 1, periodic%nb_bords

         side1 = periodic%list_periodic(1, n)
         side2 = periodic%list_periodic(2, n)
         e = periodic%vect_e(:, n)

         CALL list_periodic(mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
              list_loc, perlist_loc)

         !n_b = SIZE(list_loc)
         n_b = SIZE(perlist_loc)
         ALLOCATE(list_dom(n_b), perlist_dom(n_b))
         nx = 0
         DO i = 1, n_b
            IF (MAX(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
               nx = nx + 1
               list_dom(nx) = list_loc(i)
               perlist_dom(nx) = perlist_loc(i)
            ELSE IF (MIN(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
               WRITE(*, *) 'BUG in prep_periodic_bloc'
               STOP
            END IF
         END DO
         IF (n_b /= nx) WRITE(*, *) 'WARNING, I have removed', n_b - nx, ' periodic pairs in prep_periodic_bloc'
         n_b = nx

         nsize = nb_bloc * n_b !SIZE(list_loc) !n_b

         ALLOCATE(periodic%list(n)%DIL(nsize), periodic%perlist(n)%DIL(nsize))

         DO k = 1, nb_bloc
            k_deb = (k - 1) * n_b + 1
            k_fin = k * n_b
            periodic%list(n)%DIL(k_deb:k_fin) = list_dom(1:n_b) + (k - 1) * mesh%dom_np ! First bloc
            periodic%perlist(n)%DIL(k_deb:k_fin) = perlist_dom(1:n_b) + (k - 1) * mesh%dom_np ! First bloc
         END DO

         DEALLOCATE(list_loc, perlist_loc, list_dom, perlist_dom)

      END DO

      WRITE (*, *) 'Treatment of periodic-data done'

   END SUBROUTINE prep_periodic_bloc
!!$    

!!$  SUBROUTINE get_mesh_(communicator, mesh, opt_fe, opt_edge_stab, opt_pers)
!!$    USE mesh_1d
!!$    USE mesh_distribution_1d
!!$    USE load_mesh_2d
!!$    USE refine_mesh
!!$    USE two_dim_metis_distribution
!!$    USE gauss_points_2d
!!$    USE Dir_nodes
!!$    ! VB 2/04/2026
!!$    ! USE space_dim
!!$    USE mesh_parameters
!!$    ! VB 2/04/2026
!!$
!!$    IMPLICIT NONE
!!$    LOGICAL, OPTIONAL :: opt_edge_stab
!!$    INTEGER, OPTIONAL :: opt_fe
!!$    TYPE(periodic_type), DIMENSION(:), OPTIONAL :: opt_pers
!!$    ! TYPE(mesh_data_type) :: mesh_data
!!$    INTEGER, DIMENSION(1) :: list_dom = 1
!!$    INTEGER, DIMENSION(0) :: list_inter
!!$    INTEGER, DIMENSION(:), ALLOCATABLE :: part
!!$    INTEGER :: n, nb_proc, ierr, rank
!!$    LOGICAL :: edge_stab
!!$    CHARACTER(LEN = 100) :: mesh_part_name
!!$    TYPE(mesh_type) :: mesh_glob, mesh, mesh_r
!!$    MPI_Comm       :: communicator
!!$
!!$    CALL mesh_data_info%init
!!$    CALL MPI_Comm_SIZE(communicator, nb_proc, ierr)
!!$    CALL MPI_Comm_rank(communicator, rank, ierr)
!!$
!!$    IF (.NOT.PRESENT(opt_edge_stab)) THEN
!!$       edge_stab = .FALSE.
!!$    ELSE
!!$       edge_stab = opt_edge_stab
!!$    END IF
!!$
!!$    !=== FIXME mesh%rank to be transferred through refinement_iso_grid_distributed,
!!$    !=== create_iso_grid_distributed, copy_mesh
!!$    mesh_glob%rank = -1
!!$
!!$    SELECT CASE(mesh_data_info%k_dim)
!!$       ! SELECT CASE(k_dim)
!!$    CASE(2)
!!$       !===load and re order mesh
!!$       CALL load_dg_mesh_free_format(mesh_data_info%directory, mesh_data_info%file_name, &
!!$            list_dom, list_inter, mesh_glob, mesh_data_info%if_mesh_formatted)
!!$       ALLOCATE(part(mesh_glob%me))
!!$
!!$       mesh_part_name = 'mesh_part.' // TRIM(ADJUSTL(mesh_data_info%file_name))
!!$       IF (mesh_data_info%if_read_partition) THEN
!!$          IF (rank == 0) WRITE(*, *) 'read partition'
!!$          OPEN(UNIT = 51, FILE = mesh_part_name, STATUS = 'unknown', FORM = 'formatted')
!!$          READ(51, *) part
!!$          CLOSE(51)
!!$       ELSE
!!$          IF (rank == 0) WRITE(*, *) 'create partition'
!!$          CALL part_mesh(nb_proc, mesh_glob, list_inter, part, opt_pers)
!!$          IF (rank==0) THEN
!!$             OPEN(UNIT = 51, FILE = mesh_part_name, STATUS = 'replace', FORM = 'formatted')
!!$             WRITE(51, *) part
!!$             CLOSE(51)
!!$          END IF
!!$       END IF
!!$
!!$       CALL extract_mesh(communicator, nb_proc, mesh_glob, part, list_dom, mesh)
!!$       CALL free_mesh(mesh_glob)
!!$       DEALLOCATE(part)
!!$       !===mesh refinements
!!$       DO n = 1, mesh_data_info%nb_refinement
!!$          !===Create refined mesh
!!$          CALL refinement_iso_grid_distributed(mesh)
!!$          IF(rank == 0) write(*, *) 'refinement done', n
!!$       END DO
!!$
!!$       !===special meshes
!!$       !      IF(mesh_data_info%if_HCT) THEN
!!$       !         CALL HCT_iso_grid_distributed(mesh_p1, HCT_mesh_p1)
!!$       !         CALL deallocate_mesh(mesh_p1)
!!$       !         CALL copy_mesh(HCT_mesh_p1, mesh_p1)
!!$       !         CALL deallocate_mesh(HCT_mesh_p1)
!!$       !      END IF
!!$
!!$       !===create finite elements polynome on mesh
!!$       CALL create_iso_grid_distributed(mesh, mesh_r, mesh_data_info%type_fe)
!!$       CALL free_mesh(mesh)
!!$       CALL copy_mesh(mesh_r, mesh)
!!$       CALL free_mesh(mesh_r)
!!$
!!$       !===(JLG) Added March 21 2026 
!!$       ALLOCATE(mesh%iis(SIZE(mesh%jjs,1),mesh%mes))
!!$       CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
!!$       CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
!!$       mesh%nps = SIZE(mesh%j_s)
!!$       !===END (JLG) Added March 21 2026 
!!$
!!$       mesh%rank = rank  !=== petsc convention
!!$       mesh%edge_stab = .false.
!!$       !===gauss points on mesh
!!$       CALL create_gauss_points_2d(mesh, mesh_data_info%type_fe)
!!$
!!$    CASE(1)
!!$       CALL load_mesh_1d(mesh_data_info%directory, mesh_data_info%file_name, mesh_glob, mesh_data_info%if_mesh_formatted)
!!$       CALL extract_mesh_1d(communicator, mesh_glob, mesh, opt_pers)
!!$       CALL free_mesh(mesh_glob)
!!$       mesh%edge_stab = .false.
!!$       CALL GAUSS_POINT_1d(mesh)
!!$       mesh%rank = rank
!!$
!!$    CASE(-1)
!!$       IF (rank == 0) WRITE(*, *) 'BUG in construct_mesh: k_dim = -1 => you did not set k_dim in the data file'
!!$       STOP
!!$
!!$    CASE DEFAULT
!!$       IF(rank == 0) write(*, *) ' BUG in construct_mesh, k_dim not correct'
!!$       STOP
!!$    END SELECT
!!$
!!$    mesh%edge_stab = .false. !TODO remove edge_stab
!!$
!!$  END SUBROUTINE get_mesh_

END MODULE  construct_mesh
