MODULE mesh_distribution_1d
#include "petsc/finclude/petsc.h"
   USE petsc
   USE mesh_tools
   PUBLIC :: extract_mesh_1d, create_Pk_mesh_1D
   PRIVATE
   LOGICAL, PRIVATE :: per_bool ! <== FIXME
CONTAINS
   SUBROUTINE extract_mesh_1d(communicator, mesh_glob, mesh_loc)
      USE def_type_mesh
      USE mesh_parameters
      USE my_util, ONLY: local_error_petsc, to_str
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_glob, mesh_loc
      INTEGER :: n, m, np_start, np_end, me_start, me_end, p, n_glob
      PetscErrorCode :: ierr
      PetscMPIInt    :: rank, nb_procs
      MPI_Comm       :: communicator
      CALL MPI_Comm_rank(communicator, rank, ierr)
      CALL MPI_COMM_SIZE(communicator, nb_procs, ierr)
      rank = rank + 1
      CALL mesh_glob%create_comm(communicator)
      CALL mesh_loc%info%copy(mesh_glob%info)

      per_bool = .FALSE.
      IF (mesh_glob%info%nb_bords/=0) THEN
         !DO n = 1, SIZE(opt_pers)
         !   per_bool = per_bool .OR. opt_pers(n)%nb_bords > 0
         !END DO
         per_bool = .TRUE.
      END IF
      IF  (nb_procs == 1) THEN
         IF (per_bool) THEN
            mesh_glob%nis = 0
            DEALLOCATE(mesh_glob%isolated_jjs, mesh_glob%isolated_interfaces)
            ALLOCATE(mesh_glob%isolated_jjs(mesh_glob%nis), mesh_glob%isolated_interfaces(mesh_glob%nis, 1))
         END IF
         CALL copy_mesh(mesh_glob, mesh_loc)
         RETURN
      END IF

      mesh_loc%mi = 0
      mesh_loc%medges = 0
      mesh_loc%mes_extra = 0
      mesh_loc%mes_int = 0
      mesh_loc%edge_stab = mesh_glob%edge_stab

      np_start = (rank - 1) * (mesh_glob%np / nb_procs) + 1
      IF (rank==1) THEN
         np_end = rank * (mesh_glob%np / nb_procs)
         me_start = 1
         mesh_loc%mextra = 1
         mesh_loc%mes = 1
      ELSE IF(rank==nb_procs) THEN
         np_end = mesh_glob%np
         me_start = np_start - 1
         mesh_loc%mextra = 0
         mesh_loc%mes = 1
      ELSE
         np_end = rank * (mesh_glob%np / nb_procs)
         me_start = np_start - 1
         mesh_loc%mextra = 1
         mesh_loc%mes = 0
      END IF
      me_end = np_end - 1
      mesh_loc%me = me_end - me_start + 1
      mesh_loc%dom_np = np_end - np_start + 1
      
      
      IF (rank == 1) THEN
         mesh_loc%np = mesh_loc%dom_np
      ELSE
         mesh_loc%np = mesh_loc%dom_np + 1
      END IF

      IF (per_bool) THEN
         IF (rank == 1) THEN
            mesh_loc%dom_np = mesh_loc%dom_np - 1
            mesh_loc%me = mesh_loc%me - 1
            mesh_loc%mextra = 2
            mesh_loc%np = mesh_loc%dom_np
            mesh_loc%mes = 0
            np_start = np_start + 1
            me_start = me_start + 1
         ELSE IF (rank == nb_procs) THEN
            mesh_loc%dom_np = mesh_loc%dom_np + 1
            mesh_loc%me = mesh_loc%me + 1
            mesh_loc%mextra = 0
            mesh_loc%np = mesh_loc%dom_np + 2
            mesh_loc%mes = 2
         END IF
      END IF

      mesh_loc%medge = mesh_loc%me

      mesh_loc%nis = 0

      ALLOCATE(mesh_loc%jj(2, mesh_loc%me), mesh_loc%jjs(1, mesh_loc%mes), mesh_loc%iis(0, 0))
      ALLOCATE(mesh_loc%jj_extra(2, mesh_loc%mextra), mesh_loc%jce_extra(0, mesh_loc%mextra), &
            mesh_loc%jjs_extra(0, mesh_loc%mes_extra))
      ALLOCATE(mesh_loc%jjs_int(0, 0), mesh_loc%jcc_extra(mesh_loc%mextra), mesh_loc%jce(1, mesh_loc%me))
      ALLOCATE(mesh_loc%jees(0), mesh_loc%jecs(0))
      ALLOCATE(mesh_loc%jji(0, 0, 0), mesh_loc%jjsi(0, 0), mesh_loc%j_s(0))
      ALLOCATE(mesh_loc%rr(1, mesh_loc%np), mesh_loc%rrs_extra(1, 2, 0))
      ALLOCATE(mesh_loc%neigh(2, mesh_loc%me), mesh_loc%neighi(0, 0))
      ALLOCATE(mesh_loc%sides(mesh_loc%mes), mesh_loc%neighs(mesh_loc%mes))
      ALLOCATE(mesh_loc%sides_extra(mesh_loc%mes_extra), mesh_loc%neighs_extra(mesh_loc%mes_extra))
      ALLOCATE(mesh_loc%sides_int(mesh_loc%mes_int), mesh_loc%neighs_int(2, mesh_loc%mes_int))
      ALLOCATE(mesh_loc%i_d(mesh_loc%me), mesh_loc%loc_to_glob(mesh_loc%np), mesh_loc%proc_np_loc(2, mesh_loc%np-mesh_loc%dom_np))
      ALLOCATE(mesh_loc%isolated_jjs(mesh_loc%nis), mesh_loc%isolated_interfaces(mesh_loc%nis, 1))

      mesh_loc%loc_to_glob(:) = np_start - 1 + [(n, n=1, mesh_loc%dom_np)]
      mesh_loc%jce(1, :)      = me_start - 1 + [(n, n=1, mesh_loc%me)]

   !FIX neigh (VB 06/05/2026)
      IF (mesh_loc%me /= 0) THEN
         mesh_loc%neigh(1, :) = [(m, m = 2, mesh_loc%me + 1)]
         mesh_loc%neigh(2, :) = [(m, m = 0, mesh_loc%me - 1)]
         !=== other procs before m=1 and m=mesh_loc%me
         mesh_loc%neigh(:, 1) = -1
         mesh_loc%neigh(:, mesh_loc%me) = -1
         !=== special treatment for boundaries in the case of extremal ranks
         IF (rank == 1) THEN
            mesh_loc%neigh(2, 1) = 0
         END IF
         IF (rank == nb_procs) THEN
            mesh_loc%neigh(1, mesh_loc%me) = 0
         END IF
      END IF
   !VB 06/05/2026

      IF (per_bool) THEN
         mesh_loc%loc_to_glob = mesh_loc%loc_to_glob - 1
         mesh_loc%jce = mesh_loc%jce - 1
      END IF


      mesh_loc%i_d(1:me_end - me_start + 1) = mesh_glob%i_d(me_start:me_end)
      mesh_loc%jj(:, 1:me_end - me_start + 1) = mesh_glob%jj(:, me_start:me_end) - np_start + 1
      !=== define 1:dom_np
      mesh_loc%rr(:, 1:np_end - np_start + 1) = mesh_glob%rr(:, np_start:np_end)
      !=== define dom_np+1:np
      IF (rank > 1) THEN 
         mesh_loc%rr(:, mesh_loc%dom_np + 1) = mesh_glob%rr(:, np_start - 1)
      END IF


      IF (per_bool) THEN
         IF (rank == nb_procs) THEN
            mesh_loc%i_d(mesh_loc%me) = mesh_glob%i_d(1)
            mesh_loc%jj(1, mesh_loc%me) = mesh_loc%dom_np
            mesh_loc%rr(:, mesh_loc%dom_np)     = mesh_glob%rr(:, 1)
            mesh_loc%rr(:, mesh_loc%dom_np + 2) = mesh_glob%rr(:, 2)
         END IF
      END IF

      IF (MINVAL(mesh_loc%jj(:, :)) < 1) THEN
         IF (per_bool .AND. rank /= nb_procs) THEN
            mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 2
         ELSE
            mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
         END IF
      END IF

      IF (per_bool .AND. rank==nb_procs) THEN
         mesh_loc%loc_to_glob(mesh_loc%np - 1) = np_start - 2
         mesh_loc%loc_to_glob(mesh_loc%np) = 1
         mesh_loc%jj(1, 1) = mesh_loc%np - 1
         mesh_loc%jj(2, mesh_loc%me) = mesh_loc%np
      END IF


      IF (rank == 1) THEN
         IF (per_bool) THEN
            mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1) - 1
            mesh_loc%jcc_extra(1) = me_end! + 1 VB FIX

            mesh_loc%jcc_extra(2) = mesh_glob%me
            mesh_loc%jj_extra(2, 2) = 1
            mesh_loc%jj_extra(1, 2) = mesh_glob%np
         ELSE
            mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1)
            mesh_loc%jcc_extra = me_end + 1
            mesh_loc%sides(1) = mesh_glob%sides(1)
            mesh_loc%neighs(1) = 1
            mesh_loc%jjs(1, 1) = 1
         END IF
      ELSE IF (rank == nb_procs) THEN

         IF (per_bool) THEN
            mesh_loc%sides(1) = mesh_glob%sides(2)
            mesh_loc%neighs(1) = mesh_loc%me - 1
            mesh_loc%jjs(1, 1) = mesh_loc%dom_np - 1

            mesh_loc%sides(2) = mesh_glob%sides(1)
            mesh_loc%neighs(2) = mesh_loc%me
            mesh_loc%jjs(1, 2) = mesh_loc%dom_np
         ELSE
            mesh_loc%sides(1) = mesh_glob%sides(2)
            mesh_loc%neighs(1) = mesh_loc%me
            mesh_loc%jjs(1, 1) = mesh_loc%dom_np
         END IF
      ELSE
         mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1)
         mesh_loc%jcc_extra = me_end + 1

         IF (per_bool) THEN
            mesh_loc%jj_extra = mesh_loc%jj_extra - 1
            mesh_loc%jcc_extra = mesh_loc%jcc_extra - 1
         END IF

      END IF

      WHERE (mesh_loc%jj(1,:) < 1)
         mesh_loc%jj(1,:) = mesh_loc%np
      ELSEWHERE (mesh_loc%jj(2,:) < 1)
         mesh_loc%jj(2,:) = mesh_loc%np
      END WHERE

      CALL mesh_loc%create_comm(communicator)
      CALL mesh_loc%gather_dom_np
      CALL mesh_loc%gather_me
      CALL mesh_loc%gather_medge

      DO n=1, mesh_loc%np-mesh_loc%dom_np
         n_glob = mesh_loc%loc_to_glob(mesh_loc%dom_np + n)
         p = mesh_loc%get_proc(n_glob, 'np')
         mesh_loc%proc_np_loc(1, n) = p
         mesh_loc%proc_np_loc(2, n) = n_glob - mesh_loc%disp(p) + 1
      END DO

      !=== DEBUGGING => not supposed to change loc_to_glob
      ! write(*,*) "I", mesh_loc%proc, mesh_loc%loc_to_glob
      ! CALL mesh_loc%build_loc_to_glob
      ! write(*,*) "II", mesh_loc%proc, mesh_loc%loc_to_glob
      !=== DEBUGGING => not supposed to change loc_to_glob

      IF (mesh_loc%dom_np /= 0 .AND. mesh_loc%me == 0) THEN
         CALL local_error_petsc(&
         &"BUG in extract_mesh_1D: some proc has a single node, but no element. &
         &Think about adding more points or using less parallelization... fyi, domnp per proc = "//to_str(mesh_loc%domnp))
      END IF

   END SUBROUTINE extract_mesh_1d

   SUBROUTINE create_Pk_mesh_1D(mesh_P1, mesh_Pk, type_fe)
      USE def_type_mesh
      USE my_util, ONLY: error_petsc, to_str, local_error_petsc
      IMPLICIT NONE
      TYPE(mesh_type)     :: mesh_P1, mesh_Pk
      INTEGER, INTENT(IN) :: type_fe
      REAL(KIND=8) :: dx, x0
      INTEGER      :: n, m, n_shift, i_shift, n_loc, p, cell_l, cell_g

      CALL copy_mesh(mesh_P1, mesh_Pk)
      ! IF (type_fe == 1) RETURN

      mesh_Pk%dom_np = mesh_P1%dom_np + (type_fe-1)*mesh_P1%me
      mesh_Pk%np = mesh_P1%np + (type_fe-1)*mesh_P1%me

      DEALLOCATE(mesh_Pk%jj, mesh_Pk%rr, mesh_Pk%loc_to_glob, mesh_Pk%proc_np_loc)
      DEALLOCATE(mesh_Pk%disp, mesh_Pk%domnp, mesh_Pk%disedge, mesh_Pk%domedge, mesh_Pk%discell, mesh_Pk%domcell)
      ALLOCATE(mesh_Pk%jj(type_fe+1, mesh_P1%me), mesh_Pk%rr(1, mesh_Pk%np))
      ALLOCATE(mesh_Pk%proc_np_loc(2, mesh_Pk%np-mesh_Pk%dom_np))

      CALL mesh_Pk%gather_dom_np
      CALL mesh_Pk%gather_me
      CALL mesh_Pk%gather_medge

      mesh_Pk%proc_np_loc(:,:) = mesh_P1%proc_np_loc(:, :) !should stay the same since there is no additional point between dom_np and np

!===== DEBUGGING
      ! write(*,*) "proc_np_loc on proc ", mesh_Pk%proc, " = ", mesh_Pk%proc_np_loc!, mesh_P1%np, mesh_P1%dom_np
      ! write(*,*) mesh_Pk%proc, ": nb points Pk = ", mesh_Pk%dom_np, mesh_Pk%np, "/ nb points P1 = ", mesh_P1%dom_np, mesh_P1%np
      ! write(*,*) "P1 coordinates for proc ", mesh_P1%proc, mesh_P1%rr
      ! write(*,*) "dom_np for P1 and Pk at proc = ", mesh_Pk%proc, mesh_P1%dom_np, mesh_Pk%dom_np
      ! write(*,*) "np for P1 and Pk at proc = ", mesh_Pk%proc, mesh_P1%np, mesh_Pk%np
      ! call mpi_barrier(communicator, ierr)
      ! IF (mesh_Pk%proc==1) write(*,*) "========================================================================================"
!===== DEBUGGING



      ! copy previous points
      mesh_Pk%rr = -1 ! VB dummy init
      mesh_Pk%rr(:, 1:mesh_P1%dom_np) = mesh_P1%rr(:, 1:mesh_P1%dom_np)
      mesh_Pk%rr(:, mesh_Pk%dom_np+1:mesh_Pk%dom_np + (mesh_P1%np-mesh_P1%dom_np)) = mesh_P1%rr(:,mesh_P1%dom_np+1:mesh_P1%np)
      
      n_shift = mesh_P1%dom_np + 1 ! start appending the new Pk nodes only after the P1 nodes
      ! rebuild P1 nodes in jj

      DO m=1, mesh_P1%me
         WHERE(mesh_P1%jj(:, m) <= mesh_P1%dom_np)
            mesh_Pk%jj(1:2, m) = mesh_P1%jj(:, m)
         ELSEWHERE
            mesh_Pk%jj(1:2, m) = mesh_P1%jj(:, m) + (mesh_Pk%dom_np-mesh_P1%dom_np)
         END WHERE
      END DO

      DO m=1, mesh_P1%me

         dx = MINVAL(ABS(mesh_P1%rr(1, mesh_P1%jj(2:, m)) - mesh_P1%rr(1, mesh_P1%jj(1, m))))
         x0 = MINVAL(mesh_P1%rr(1, mesh_P1%jj(:, m))) + dx

         i_shift = n_shift + (m-1)*(type_fe-1)
         mesh_Pk%rr(1,i_shift:i_shift+type_fe-2) = x0 - dx/(type_fe*1.d0)*[(n, n=1 ,type_fe-1)]
      END DO

      DO m=1, mesh_P1%me
         ! appending Pk nodes
         mesh_Pk%jj(3:type_fe+1, m) = n_shift + (m-1)*(type_fe - 1) + [(n, n=type_fe-2, 0, -1)]
      END DO

!=== DEBUGGING
      ! DO m=1, mesh_Pk%me
      !    write(*,*) "P1 element ", m, ' on proc ', mesh_P1%proc, mesh_P1%rr(1, mesh_P1%jj(:, m)), mesh_P1%jj(:, m)
      !    write(*,*) "Pk element ", m, ' on proc ', mesh_Pk%proc, mesh_Pk%rr(1, mesh_Pk%jj(:, m)), mesh_Pk%jj(:, m)
      ! END DO
      ! write(*,*) "P1 coords dom_np on proc ", mesh_P1%proc, mesh_P1%rr(1,:mesh_P1%dom_np)
      ! write(*,*) "P1 coords np on proc ", mesh_P1%proc, mesh_P1%rr(1,mesh_P1%dom_np+1:)
      ! write(*,*) "Pk coords dom_np on proc ", mesh_Pk%proc, mesh_Pk%rr(1,:mesh_Pk%dom_np)
      ! write(*,*) "Pk coords np on proc ", mesh_Pk%proc, mesh_Pk%rr(1,mesh_Pk%dom_np+1:)
!=== DEBUGGING
      
      CALL mesh_Pk%build_loc_to_glob

      !=== Construct extra layer stuff
      DEALLOCATE(mesh_Pk%rrs_extra)
      DEALLOCATE(mesh_Pk%jj_extra)
      ALLOCATE(mesh_Pk%rrs_extra(1, type_fe+1, 0))
      ALLOCATE(mesh_Pk%jj_extra(type_fe+1, mesh_Pk%mextra))

   
      DO m=1, mesh_P1%mextra
         !=== start with endpoints
         DO n=1, 2
            p = mesh_P1%get_proc(mesh_P1%jj_extra(n, m), 'np')
            n_loc = mesh_P1%jj_extra(n, m) - mesh_P1%disp(p) + 1
            mesh_Pk%jj_extra(n, m) = mesh_Pk%global_numbering(p, n_loc)
         END DO
         !=== Then new inner points
         cell_g = mesh_P1%jcc_extra(m)
         p = mesh_P1%get_proc(cell_g, 'me')
         cell_l = cell_g - (mesh_P1%discell(p) - 1)
         DO n=3, type_fe+1
            n_loc = mesh_P1%domnp(p) + (cell_l-1)*(type_fe-1) + (type_fe-2 + (3-n)) + 1
            mesh_Pk%jj_extra(n, m) = mesh_Pk%global_numbering(p, n_loc)
         END DO
      END DO

!=== DEBUGGING
      ! WRITE(*,*) "on proc ", mesh_P1%proc, " there are ", mesh_Pk%mextra, " extra cells"
      ! DO m=1,mesh_Pk%mextra
      !    write(*,*) type_fe, "P1 jj_extra ", m, ' on proc ', mesh_P1%proc, mesh_P1%jj_extra(:, m)
      !    write(*,*) type_fe, "Pk jj_extra ", m, ' on proc ', mesh_Pk%proc, mesh_Pk%jj_extra(:, m)
      !    write(*,*) type_fe, "P1 jcc_extra ", m, ' on proc ', mesh_Pk%proc, mesh_Pk%jcc_extra(m)
      ! END DO
      ! write(*,*) 'P1 loc_to_glob on proc ', mesh_P1%proc, mesh_P1%loc_to_glob
      ! write(*,*) 'P1 rr on proc ', mesh_P1%proc, mesh_P1%rr
      ! write(*,*) 'Pk loc_to_glob on proc ', mesh_Pk%proc, mesh_Pk%loc_to_glob
      ! write(*,*) 'Pk rr on proc ', mesh_Pk%proc, mesh_Pk%rr
!=== DEBUGGING
   END SUBROUTINE create_Pk_mesh_1D

END MODULE mesh_distribution_1d
