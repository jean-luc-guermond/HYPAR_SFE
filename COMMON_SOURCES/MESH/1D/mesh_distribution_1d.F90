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
      per_bool = .FALSE.
      IF (mesh_data_info%nb_bords/=0) THEN
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
            np_end = np_end + 1
            me_end = me_end + 1
         END IF
      END IF

      mesh_loc%medge = mesh_loc%me

      mesh_loc%nis = 0

      mesh_loc%dom_me = mesh_loc%me
      mesh_loc%dom_mes = mesh_loc%mes

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

      DO n = 1, mesh_loc%dom_np
         mesh_loc%loc_to_glob(n) = np_start - 1 + n
      END DO
      DO n = 1, mesh_loc%me
         mesh_loc%jce(1, n) = me_start - 1 + n
      END DO
   !FIX neigh (VB 06/05/2026)
      mesh_loc%neigh(1, :) = [(m, m = 2, mesh_loc%me + 1)]
      mesh_loc%neigh(2, :) = [(m, m = 0, mesh_loc%me - 1)]
      !=== other procs before m=1 and m=mesh_loc%me
      mesh_loc%neigh(:, 1) = -1
      mesh_loc%neigh(:, mesh_loc%me) = -1
      !=== special treatment for boundaries in the case of extremal ranks
      IF (rank == 1) THEN
         mesh_loc%neigh(2, 1) = 0
      ELSE IF (rank == nb_procs) THEN
         mesh_loc%neigh(1, mesh_loc%me) = 0
      END IF
   !VB 06/05/2026

      IF (per_bool) THEN
         mesh_loc%loc_to_glob = mesh_loc%loc_to_glob - 1
         mesh_loc%jce = mesh_loc%jce - 1
      END IF

      IF (rank == 1) THEN
         mesh_loc%i_d = mesh_glob%i_d(me_start:me_end)
         mesh_loc%jj = mesh_glob%jj(:, me_start:me_end) - np_start + 1
         mesh_loc%rr(:, 1:mesh_loc%dom_np) = mesh_glob%rr(:, np_start:np_end)

         IF (per_bool) THEN
            mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1) - 1
            mesh_loc%jcc_extra(1) = me_end + 1

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
            mesh_loc%i_d(1:mesh_loc%me - 1) = mesh_glob%i_d(me_start:me_end - 1)
            mesh_loc%i_d(mesh_loc%me) = mesh_glob%i_d(1)
            mesh_loc%jj(:, 1:mesh_loc%me - 1) = mesh_glob%jj(:, me_start:me_end - 1) - np_start + 1
            mesh_loc%jj(1, mesh_loc%me) = mesh_loc%dom_np
            mesh_loc%rr(:, 1:mesh_loc%dom_np - 1) = mesh_glob%rr(:, np_start:np_end - 1)
            mesh_loc%rr(:, mesh_loc%dom_np) = mesh_glob%rr(:, 1)
         ELSE
            mesh_loc%i_d = mesh_glob%i_d(me_start:me_end)
            mesh_loc%jj = mesh_glob%jj(:, me_start:me_end) - np_start + 1
            mesh_loc%rr(:, 1:mesh_loc%dom_np) = mesh_glob%rr(:, np_start:np_end)
         END IF

         IF (per_bool) THEN
            mesh_loc%rr(:, mesh_loc%np - 1) = mesh_glob%rr(:, np_start - 1)
            mesh_loc%rr(:, mesh_loc%np) = mesh_glob%rr(:, 2)

            mesh_loc%loc_to_glob(mesh_loc%np - 1) = np_start - 2
            mesh_loc%jj(1, 1) = mesh_loc%np - 1

            mesh_loc%loc_to_glob(mesh_loc%np) = 1
            mesh_loc%jj(2, mesh_loc%me) = mesh_loc%np

            mesh_loc%sides(1) = mesh_glob%sides(2)
            mesh_loc%neighs(1) = mesh_loc%me - 1
            mesh_loc%jjs(1, 1) = mesh_loc%dom_np - 1

            mesh_loc%sides(2) = mesh_glob%sides(1)
            mesh_loc%neighs(2) = mesh_loc%me
            mesh_loc%jjs(1, 2) = mesh_loc%dom_np
         ELSE
            mesh_loc%rr(:, mesh_loc%np) = mesh_glob%rr(:, np_start - 1)
            DO m = 1, mesh_loc%me
               IF (mesh_loc%jj(1, m) < 1) THEN
                  mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
                  mesh_loc%jj(1, m) = mesh_loc%np
               END IF
               IF (mesh_loc%jj(2, m) < 1) THEN
                  mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
                  mesh_loc%jj(2, m) = mesh_loc%np
               END IF
            END DO
            mesh_loc%sides(1) = mesh_glob%sides(2)
            mesh_loc%neighs(1) = mesh_loc%me
            mesh_loc%jjs(1, 1) = mesh_loc%dom_np
         END IF
      ELSE
         mesh_loc%i_d = mesh_glob%i_d(me_start:me_end)
         mesh_loc%jj = mesh_glob%jj(:, me_start:me_end) - np_start + 1
         mesh_loc%rr(:, 1:mesh_loc%dom_np) = mesh_glob%rr(:, np_start:np_end)

         mesh_loc%rr(:, mesh_loc%np) = mesh_glob%rr(:, np_start - 1)
         mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1)
         mesh_loc%jcc_extra = me_end + 1
         DO m = 1, mesh_loc%me
            IF (mesh_loc%jj(1, m) < 1) THEN
               IF (per_bool) THEN
                  mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 2
               ELSE
                  mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
               END IF
               mesh_loc%jj(1, m) = mesh_loc%np
            END IF
            IF (mesh_loc%jj(2, m) < 1) THEN
               IF (per_bool) THEN
                  mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 2
               ELSE
                  mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
               END IF
               mesh_loc%jj(2, m) = mesh_loc%np
            END IF
         END DO

         IF (per_bool) THEN
            mesh_loc%jj_extra = mesh_loc%jj_extra - 1
            mesh_loc%jcc_extra = mesh_loc%jcc_extra - 1
         END IF

      END IF

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

      !=== Testing => not supposed to change loc_to_glob
      ! write(*,*) "I", mesh_loc%proc, mesh_loc%loc_to_glob
      ! CALL mesh_loc%build_loc_to_glob
      ! write(*,*) "II", mesh_loc%proc, mesh_loc%loc_to_glob
      !=== Testing => not supposed to change loc_to_glob

   END SUBROUTINE extract_mesh_1d

   SUBROUTINE create_Pk_mesh_1D(communicator, mesh_P1, mesh_Pk, type_fe)
      USE def_type_mesh
      USE my_util, ONLY: error_petsc, to_str
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_P1, mesh_Pk
      REAL(KIND=8) :: dx, x0
      INTEGER :: n, m, i, n_shift, i_shift, n_loc, other_proc, other_m_loc, p, k, jj_shift, dom_np_diff, np_diff
      INTEGER, INTENT(IN) :: type_fe
      integer :: ierr
      MPI_Comm :: communicator

      CALL copy_mesh(mesh_P1, mesh_Pk)
      ! IF (type_fe == 1) RETURN

      mesh_Pk%dom_np = mesh_P1%dom_np + (type_fe-1)*mesh_P1%me
      mesh_Pk%np = mesh_P1%np + (type_fe-1)*mesh_P1%me

      DEALLOCATE(mesh_Pk%jj, mesh_Pk%rr, mesh_Pk%loc_to_glob, mesh_Pk%proc_np_loc)!, mesh_Pk%proc_np_loc)
      DEALLOCATE(mesh_Pk%disp, mesh_Pk%domnp, mesh_Pk%disedge, mesh_Pk%domedge, mesh_Pk%discell, mesh_Pk%domcell)
      ALLOCATE(mesh_Pk%jj(type_fe+1, mesh_P1%me), mesh_Pk%rr(1, mesh_Pk%np))
      ALLOCATE(mesh_Pk%proc_np_loc(2, mesh_Pk%np-mesh_Pk%dom_np))

      ! CALL mesh_Pk%create_comm(communicator)
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
      mesh_Pk%rr(:, mesh_Pk%dom_np+1:mesh_Pk%np + (mesh_P1%np-mesh_P1%dom_np)) = mesh_P1%rr(:,mesh_P1%dom_np+1:mesh_P1%np)
      
      n_shift = mesh_P1%dom_np + 1 ! start appending the new Pk nodes only after the P1 nodes
      ! rebuild P1 nodes in jj

      DO m=1, mesh_P1%me
         WHERE(mesh_P1%jj(:, m) <= mesh_P1%dom_np)
            mesh_Pk%jj(:, m) = mesh_P1%jj(:, m)
         ELSEWHERE
            mesh_Pk%jj(:, m) = mesh_P1%jj(:, m) + (mesh_Pk%dom_np-mesh_P1%dom_np)
         END WHERE
      END DO

      dx = mesh_P1%rr(1, 2) - mesh_P1%rr(1, 1) ! FIXME WARNING: this is assuming uniform mesh
      DO m=1, mesh_P1%me
         ! dx = ABS(mesh_P1%jj(2, m) - mesh_P1%jj(1, m)) !<== FIXME

         IF ((mesh_Pk%proc==1) .OR. (per_bool .AND. m==mesh_P1%me .AND. mesh_Pk%proc==mesh_Pk%nb_proc)) THEN
            x0 = mesh_P1%rr(1, m) + dx
         ELSE
            x0 = mesh_P1%rr(1, m)
         END IF
         i_shift = n_shift + (m-1)*(type_fe-1)
         mesh_Pk%rr(1,i_shift:i_shift+type_fe-2) = x0 - dx/(type_fe*1.d0)*[(n, n=1 ,type_fe-1)]
      END DO
      
      DO m=1, mesh_P1%me
         ! appending Pk nodes
         mesh_Pk%jj(3:type_fe+1, m) = n_shift + (m-1)*(type_fe - 1) + [(n, n=type_fe-2, 0, -1)] ! may fail for fe = 3 (see Gauss points ordering)?
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
      ALLOCATE(mesh_Pk%rrs_extra(1, type_fe+1, 0))

      DEALLOCATE(mesh_Pk%jj_extra)
      ALLOCATE(mesh_Pk%jj_extra(type_fe+1, mesh_Pk%mextra))
      DO m=1, mesh_Pk%mextra
         other_proc = -1
         other_m_loc = -1
         !=== start with old P1 nodes which local numbering is unchanged
         DO n=1, 2
            p = mesh_P1%get_proc(mesh_P1%jj_extra(n, m), 'np')
            n_loc = mesh_P1%jj_extra(n, m) - mesh_P1%disp(p) + 1
            mesh_Pk%jj_extra(n, m) = mesh_PK%global_numbering(p, n_loc)
            !=== Initial P1 ordering => right node ID = element ID
            other_m_loc = MAX(other_m_loc, n_loc)
            !=== Initial P1 ordering => right node ID = element ID
            IF (p /= mesh_Pk%proc) THEN
               other_proc = p
            END IF
         END DO
         !=== new Pk nodes
         IF (other_proc == -1) THEN
            CALL error_petsc("BUG in create_Pk_mesh_1d: element from extra layer of proc " // to_str(mesh_Pk%proc)&
            //" did not contain nodes from another proc")
         END IF

         IF (per_bool .AND. m==2) THEN ! different treatment for the periodic extra-cell
            p = other_proc
            n_shift = mesh_P1%domnp(p) + (mesh_P1%domcell(p) - 1)*(type_fe-1)
         ELSE
            p = other_proc
            n_shift = mesh_P1%domnp(p)
         END IF
         
         DO n = 3, type_fe + 1
            n_loc = n_shift + (type_fe + 2 - n)
            mesh_Pk%jj_extra(n, m) = mesh_PK%global_numbering(p, n_loc)
         END DO
      END DO

!=== DEBUGGING
      ! DO m=1, SIZE(mesh_Pk%jj_extra,2)
      !    write(*,*) "P1 jj_extra ", m, ' on proc ', mesh_P1%proc, mesh_P1%jj_extra(:, m)
      !    write(*,*) "Pk jj_extra ", m, ' on proc ', mesh_Pk%proc, mesh_Pk%jj_extra(:, m)
      ! END DO
      ! write(*,*) 'P1 loc_to_glob on proc ', mesh_P1%proc, mesh_P1%loc_to_glob
      ! write(*,*) 'P1 rr on proc ', mesh_P1%proc, mesh_P1%rr
      ! write(*,*) 'Pk loc_to_glob on proc ', mesh_Pk%proc, mesh_Pk%loc_to_glob
      ! write(*,*) 'Pk rr on proc ', mesh_Pk%proc, mesh_Pk%rr
!=== DEBUGGING
      ! stop
   END SUBROUTINE create_Pk_mesh_1D


   !====== Old subroutine, taking opt_pers as argument
   !  SUBROUTINE extract_mesh_1d_(communicator, mesh_glob, mesh_loc, opt_pers)
   !    USE def_type_mesh
   !    USE periodic_data_module
   !    IMPLICIT NONE
   !    TYPE(periodic_type), DIMENSION(:), OPTIONAL :: opt_pers
   !    LOGICAL :: per_bool
   !    TYPE(mesh_type) :: mesh_glob, mesh_loc
   !    INTEGER :: n, m, np_start, np_end, me_start, me_end
   !    PetscErrorCode :: ierr
   !    PetscMPIInt    :: rank, nb_procs
   !    MPI_Comm       :: communicator
   !    CALL MPI_Comm_rank(communicator, rank, ierr)
   !    CALL MPI_COMM_SIZE(communicator, nb_procs, ierr)
   !    rank = rank + 1

   !    per_bool = .FALSE.
   !    IF (PRESENT(opt_pers)) THEN
   !       DO n = 1, SIZE(opt_pers)
   !          per_bool = per_bool .OR. opt_pers(n)%nb_bords > 0
   !       END DO
   !    END IF
   !    IF  (nb_procs == 1) THEN
   !       IF (per_bool) THEN
   !          mesh_glob%nis = 0
   !          DEALLOCATE(mesh_glob%isolated_jjs, mesh_glob%isolated_interfaces)
   !          ALLOCATE(mesh_glob%isolated_jjs(mesh_glob%nis), mesh_glob%isolated_interfaces(mesh_glob%nis, 1))
   !       END IF
   !       CALL copy_mesh(mesh_glob, mesh_loc)
   !       RETURN
   !    END IF

   !    mesh_loc%mi = 0
   !    mesh_loc%medges = 0
   !    mesh_loc%mes_extra = 0
   !    mesh_loc%mes_int = 0
   !    mesh_loc%edge_stab = mesh_glob%edge_stab
   !    IF (rank < nb_procs) THEN
   !       np_start = (rank - 1) * (mesh_glob%np / nb_procs) + 1
   !       np_end = rank * (mesh_glob%np / nb_procs)
   !       mesh_loc%dom_np = np_end - np_start + 1
   !       IF (rank == 1) THEN
   !          me_start = 1
   !       ELSE
   !          me_start = np_start - 1
   !       END IF
   !       me_end = np_end - 1
   !       mesh_loc%me = me_end - me_start + 1
   !       mesh_loc%mextra = 1
   !    ELSE
   !       np_start = (rank - 1) * (mesh_glob%np / nb_procs) + 1
   !       np_end = mesh_glob%np
   !       mesh_loc%dom_np = np_end - np_start + 1
   !       me_start = np_start - 1
   !       me_end = np_end - 1
   !       mesh_loc%me = me_end - me_start + 1
   !       mesh_loc%mextra = 0
   !    END IF

   !    IF (rank == 1) THEN
   !       mesh_loc%np = mesh_loc%dom_np
   !       mesh_loc%mes = 1
   !    ELSE IF (rank == nb_procs) THEN
   !       mesh_loc%np = mesh_loc%dom_np + 1
   !       mesh_loc%mes = 1
   !    ELSE
   !       mesh_loc%np = mesh_loc%dom_np + 1
   !       mesh_loc%mes = 0
   !    END IF

   !    IF (per_bool) THEN
   !       IF (rank == 1) THEN
   !          mesh_loc%dom_np = mesh_loc%dom_np - 1
   !          mesh_loc%me = mesh_loc%me - 1
   !          mesh_loc%mextra = 2
   !          mesh_loc%np = mesh_loc%dom_np
   !          mesh_loc%mes = 0
   !          np_start = np_start + 1
   !          me_start = me_start + 1
   !       ELSE IF (rank == nb_procs) THEN
   !          mesh_loc%dom_np = mesh_loc%dom_np + 1
   !          mesh_loc%me = mesh_loc%me + 1
   !          mesh_loc%mextra = 0
   !          mesh_loc%np = mesh_loc%dom_np + 2
   !          mesh_loc%mes = 2
   !          np_end = np_end + 1
   !          me_end = me_end + 1
   !       END IF
   !    END IF

   !    mesh_loc%medge = mesh_loc%me

   !    mesh_loc%nis = 0

   !    mesh_loc%dom_me = mesh_loc%me
   !    mesh_loc%dom_mes = mesh_loc%mes

   !    ALLOCATE(mesh_loc%jj(2, mesh_loc%me), mesh_loc%jjs(1, mesh_loc%mes), mesh_loc%iis(0, 0))
   !    ALLOCATE(mesh_loc%jj_extra(2, mesh_loc%mextra), mesh_loc%jce_extra(0, mesh_loc%mextra), &
   !         mesh_loc%jjs_extra(0, mesh_loc%mes_extra))
   !    ALLOCATE(mesh_loc%jjs_int(0, 0), mesh_loc%jcc_extra(mesh_loc%mextra), mesh_loc%jce(1, mesh_loc%me))
   !    ALLOCATE(mesh_loc%jees(0), mesh_loc%jecs(0))
   !    ALLOCATE(mesh_loc%jji(0, 0, 0), mesh_loc%jjsi(0, 0), mesh_loc%j_s(0))
   !    ALLOCATE(mesh_loc%rr(1, mesh_loc%np), mesh_loc%rrs_extra(1, 2, 0))
   !    ALLOCATE(mesh_loc%neigh(2, mesh_loc%me), mesh_loc%neighi(0, 0))
   !    ALLOCATE(mesh_loc%sides(mesh_loc%mes), mesh_loc%neighs(mesh_loc%mes))
   !    ALLOCATE(mesh_loc%sides_extra(mesh_loc%mes_extra), mesh_loc%neighs_extra(mesh_loc%mes_extra))
   !    ALLOCATE(mesh_loc%sides_int(mesh_loc%mes_int), mesh_loc%neighs_int(2, mesh_loc%mes_int))
   !    ALLOCATE(mesh_loc%i_d(mesh_loc%me), mesh_loc%loc_to_glob(mesh_loc%np))
   !    ALLOCATE(mesh_loc%disp(nb_procs + 1), mesh_loc%disedge(nb_procs + 1), mesh_loc%discell(nb_procs + 1))
   !    ALLOCATE(mesh_loc%domnp(nb_procs), mesh_loc%domedge(nb_procs), mesh_loc%domcell(nb_procs))
   !    ALLOCATE(mesh_loc%isolated_jjs(mesh_loc%nis), mesh_loc%isolated_interfaces(mesh_loc%nis, 1))

   !    DO n = 1, mesh_loc%dom_np
   !       mesh_loc%loc_to_glob(n) = np_start - 1 + n
   !    END DO
   !    DO n = 1, mesh_loc%me
   !       mesh_loc%jce(1, n) = me_start - 1 + n
   !    END DO

   !    IF (per_bool) THEN
   !       mesh_loc%loc_to_glob = mesh_loc%loc_to_glob - 1
   !       mesh_loc%jce = mesh_loc%jce - 1
   !    END IF

   !    IF (rank == 1) THEN
   !       mesh_loc%i_d = mesh_glob%i_d(me_start:me_end)
   !       mesh_loc%jj = mesh_glob%jj(:, me_start:me_end) - np_start + 1
   !       mesh_loc%rr(:, 1:mesh_loc%dom_np) = mesh_glob%rr(:, np_start:np_end)

   !       IF (per_bool) THEN
   !          mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1) - 1
   !          mesh_loc%jcc_extra(1) = me_end + 1

   !          mesh_loc%jcc_extra(2) = mesh_glob%me
   !          mesh_loc%jj_extra(2, 2) = 1
   !          mesh_loc%jj_extra(1, 2) = mesh_glob%np
   !       ELSE
   !          mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1)
   !          mesh_loc%jcc_extra = me_end + 1
   !          mesh_loc%sides(1) = mesh_glob%sides(1)
   !          mesh_loc%neighs(1) = 1
   !          mesh_loc%jjs(1, 1) = 1
   !       END IF
   !    ELSE IF (rank == nb_procs) THEN
   !       IF (per_bool) THEN
   !          mesh_loc%i_d(1:mesh_loc%me - 1) = mesh_glob%i_d(me_start:me_end - 1)
   !          mesh_loc%i_d(mesh_loc%me) = mesh_glob%i_d(1)
   !          mesh_loc%jj(:, 1:mesh_loc%me - 1) = mesh_glob%jj(:, me_start:me_end - 1) - np_start + 1
   !          mesh_loc%jj(1, mesh_loc%me) = mesh_loc%dom_np
   !          mesh_loc%rr(:, 1:mesh_loc%dom_np - 1) = mesh_glob%rr(:, np_start:np_end - 1)
   !          mesh_loc%rr(:, mesh_loc%dom_np) = mesh_glob%rr(:, 1)
   !       ELSE
   !          mesh_loc%i_d = mesh_glob%i_d(me_start:me_end)
   !          mesh_loc%jj = mesh_glob%jj(:, me_start:me_end) - np_start + 1
   !          mesh_loc%rr(:, 1:mesh_loc%dom_np) = mesh_glob%rr(:, np_start:np_end)
   !       END IF

   !       IF (per_bool) THEN
   !          mesh_loc%rr(:, mesh_loc%np - 1) = mesh_glob%rr(:, np_start - 1)
   !          mesh_loc%rr(:, mesh_loc%np) = mesh_glob%rr(:, 2)

   !          mesh_loc%loc_to_glob(mesh_loc%np - 1) = np_start - 2
   !          mesh_loc%jj(1, 1) = mesh_loc%np - 1

   !          mesh_loc%loc_to_glob(mesh_loc%np) = 1
   !          mesh_loc%jj(2, mesh_loc%me) = mesh_loc%np

   !          mesh_loc%sides(1) = mesh_glob%sides(2)
   !          mesh_loc%neighs(1) = mesh_loc%me - 1
   !          mesh_loc%jjs(1, 1) = mesh_loc%dom_np - 1

   !          mesh_loc%sides(2) = mesh_glob%sides(1)
   !          mesh_loc%neighs(2) = mesh_loc%me
   !          mesh_loc%jjs(1, 2) = mesh_loc%dom_np
   !       ELSE
   !          mesh_loc%rr(:, mesh_loc%np) = mesh_glob%rr(:, np_start - 1)
   !          DO m = 1, mesh_loc%me
   !             IF (mesh_loc%jj(1, m) < 1) THEN
   !                mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
   !                mesh_loc%jj(1, m) = mesh_loc%np
   !             END IF
   !             IF (mesh_loc%jj(2, m) < 1) THEN
   !                mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
   !                mesh_loc%jj(2, m) = mesh_loc%np
   !             END IF
   !          END DO
   !          mesh_loc%sides(1) = mesh_glob%sides(2)
   !          mesh_loc%neighs(1) = mesh_loc%me
   !          mesh_loc%jjs(1, 1) = mesh_loc%dom_np
   !       END IF
   !    ELSE
   !       mesh_loc%i_d = mesh_glob%i_d(me_start:me_end)
   !       mesh_loc%jj = mesh_glob%jj(:, me_start:me_end) - np_start + 1
   !       mesh_loc%rr(:, 1:mesh_loc%dom_np) = mesh_glob%rr(:, np_start:np_end)

   !       mesh_loc%rr(:, mesh_loc%np) = mesh_glob%rr(:, np_start - 1)
   !       mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1)
   !       mesh_loc%jcc_extra = me_end + 1
   !       DO m = 1, mesh_loc%me
   !          IF (mesh_loc%jj(1, m) < 1) THEN
   !             IF (per_bool) THEN
   !                mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 2
   !             ELSE
   !                mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
   !             END IF
   !             mesh_loc%jj(1, m) = mesh_loc%np
   !          END IF
   !          IF (mesh_loc%jj(2, m) < 1) THEN
   !             IF (per_bool) THEN
   !                mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 2
   !             ELSE
   !                mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
   !             END IF
   !             mesh_loc%jj(2, m) = mesh_loc%np
   !          END IF
   !       END DO

   !       IF (per_bool) THEN
   !          mesh_loc%jj_extra = mesh_loc%jj_extra - 1
   !          mesh_loc%jcc_extra = mesh_loc%jcc_extra - 1
   !       END IF

   !    END IF

   !    CALL MPI_ALLGATHER(mesh_loc%dom_np, 1, MPI_INTEGER, mesh_loc%domnp, 1, &
   !         MPI_INTEGER, communicator, ierr)
   !    mesh_loc%disp(1) = 1
   !    DO n = 1, nb_procs
   !       mesh_loc%disp(n + 1) = mesh_loc%disp(n) + mesh_loc%domnp(n)
   !    END DO

   !    CALL MPI_ALLGATHER(mesh_loc%me, 1, MPI_INTEGER, mesh_loc%domcell, 1, &
   !         MPI_INTEGER, communicator, ierr)
   !    mesh_loc%discell(1) = 1
   !    DO n = 1, nb_procs
   !       mesh_loc%discell(n + 1) = mesh_loc%discell(n) + mesh_loc%domcell(n)
   !    END DO

   !    CALL MPI_ALLGATHER(mesh_loc%medge, 1, MPI_INTEGER, mesh_loc%domedge, 1, &
   !         MPI_INTEGER, PETSC_COMM_WORLD, ierr)
   !    mesh_loc%disedge(1) = 1
   !    DO n = 1, nb_procs
   !       mesh_loc%disedge(n + 1) = mesh_loc%disedge(n) + mesh_loc%domedge(n)
   !    END DO

   !    !write(*,*) mesh_loc%domnp, mesh_loc%domcell, mesh_loc%domedge
   !  END SUBROUTINE extract_mesh_1d_


END MODULE mesh_distribution_1d
