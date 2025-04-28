MODULE mesh_distribution_1d
#include "petsc/finclude/petsc.h"
   USE petsc
   USE mesh_tools
   USE input_periodic_data
   PUBLIC :: extract_mesh_1d
   PRIVATE
CONTAINS
   SUBROUTINE extract_mesh_1d(communicator, mesh_glob, mesh_loc, opt_per)
      USE def_type_mesh
      IMPLICIT NONE
      LOGICAL, OPTIONAL :: opt_per
      LOGICAL :: per_bool
      TYPE(mesh_type) :: mesh_glob, mesh_loc
      INTEGER :: n, m, np_start, np_end, me_start, me_end
      PetscErrorCode :: ierr
      PetscMPIInt    :: rank, nb_procs
      MPI_Comm       :: communicator
      CALL MPI_Comm_rank(communicator, rank, ierr)
      CALL MPI_COMM_SIZE(communicator, nb_procs, ierr)
      rank = rank + 1

      IF (PRESENT(opt_per)) THEN
         per_bool = opt_per
      ELSE
         per_bool = .false.
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
      mesh_loc%medge = 0
      mesh_loc%medges = 0
      mesh_loc%mes_extra = 0
      mesh_loc%mes_int = 0
      mesh_loc%edge_stab = mesh_glob%edge_stab

      IF (rank < nb_procs) THEN
         np_start = (rank - 1) * (mesh_glob%np / nb_procs) + 1
         np_end = rank * (mesh_glob%np / nb_procs)
         mesh_loc%dom_np = np_end - np_start + 1
         IF (rank == 1) THEN
            me_start = 1
         ELSE
            me_start = np_start - 1
         END IF
         me_end = np_end - 1
         mesh_loc%me = me_end - me_start + 1

         mesh_loc%mextra = 1
      ELSE
         np_start = (rank - 1) * (mesh_glob%np / nb_procs) + 1
         np_end = mesh_glob%np
         mesh_loc%dom_np = np_end - np_start + 1

         me_start = np_start - 1
         me_end = np_end - 1
         mesh_loc%me = me_end - me_start + 1

         mesh_loc%mextra = 0
      END IF

      IF (rank == 1) THEN
         mesh_loc%np = mesh_loc%dom_np
         mesh_loc%nis = 1
         mesh_loc%mes = 1
      ELSE IF (rank == nb_procs) THEN
         mesh_loc%np = mesh_loc%dom_np + 1
         mesh_loc%nis = 1
         mesh_loc%mes = 1
      ELSE
         mesh_loc%np = mesh_loc%dom_np + 1
         mesh_loc%nis = 0
         mesh_loc%mes = 0
      END IF

      IF (per_bool) THEN
         IF (rank == 1) THEN
            mesh_loc%dom_np = mesh_loc%dom_np - 1
            mesh_loc%me = mesh_loc%me - 1
            mesh_loc%mextra = 2
            mesh_loc%np = mesh_loc%dom_np
            mesh_loc%nis = 0
            mesh_loc%mes = 0
            np_start = np_start + 1
            me_start = me_start + 1
         ELSE IF (rank == nb_procs) THEN
            mesh_loc%dom_np = mesh_loc%dom_np + 1
            mesh_loc%me = mesh_loc%me + 1
            mesh_loc%mextra = 0
            mesh_loc%np = mesh_loc%dom_np + 2
            mesh_loc%nis = 0
            mesh_loc%mes = 2
            np_end = np_end + 1
            me_end = me_end + 1
         END IF
      END IF

      mesh_loc%dom_me = mesh_loc%me
      mesh_loc%dom_mes = mesh_loc%mes

      ALLOCATE(mesh_loc%jj(2, mesh_loc%me), mesh_loc%jjs(1, mesh_loc%mes), mesh_loc%iis(0, 0))
      ALLOCATE(mesh_loc%jj_extra(2, mesh_loc%mextra), mesh_loc%jce_extra(0, mesh_loc%medge), &
           mesh_loc%jjs_extra(0, mesh_loc%mes_extra))
      ALLOCATE(mesh_loc%jjs_int(0, 0), mesh_loc%jcc_extra(mesh_loc%mextra), mesh_loc%jce(0, 0))
      ALLOCATE(mesh_loc%jees(0), mesh_loc%jecs(0))
      ALLOCATE(mesh_loc%jji(0, 0, 0), mesh_loc%jjsi(0, 0), mesh_loc%j_s(0))
      ALLOCATE(mesh_loc%rr(1, mesh_loc%np), mesh_loc%rrs_extra(1, 2, 0))
      ALLOCATE(mesh_loc%neigh(2, mesh_loc%me), mesh_loc%neighi(0, 0))
      ALLOCATE(mesh_loc%sides(mesh_loc%mes), mesh_loc%neighs(mesh_loc%mes))
      ALLOCATE(mesh_loc%sides_extra(mesh_loc%mes_extra), mesh_loc%neighs_extra(mesh_loc%mes_extra))
      ALLOCATE(mesh_loc%sides_int(mesh_loc%mes_int), mesh_loc%neighs_int(2, mesh_loc%mes_int))
      ALLOCATE(mesh_loc%i_d(mesh_loc%me), mesh_loc%loc_to_glob(mesh_loc%np))
      ALLOCATE(mesh_loc%disp(nb_procs + 1), mesh_loc%disedge(nb_procs + 1), mesh_loc%discell(nb_procs + 1))
      ALLOCATE(mesh_loc%domnp(nb_procs), mesh_loc%domedge(nb_procs), mesh_loc%domcell(nb_procs))
      ALLOCATE(mesh_loc%isolated_jjs(mesh_loc%nis), mesh_loc%isolated_interfaces(mesh_loc%nis, 1))

      DO n = 1, mesh_loc%dom_np
         mesh_loc%loc_to_glob(n) = np_start - 1 + n
      END DO
      IF (per_bool) THEN
         mesh_loc%loc_to_glob = mesh_loc%loc_to_glob - 1
      END IF

      IF (rank == 1) THEN
         mesh_loc%i_d = mesh_glob%i_d(me_start:me_end)
         mesh_loc%jj = mesh_glob%jj(:, me_start:me_end) - np_start + 1
         mesh_loc%rr(:, 1:mesh_loc%dom_np) = mesh_glob%rr(:, np_start:np_end)

         IF (per_bool) THEN
            mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1)
            mesh_loc%jcc_extra(1) = me_end
            mesh_loc%jj_extra(1, 2) = 1
            mesh_loc%jj_extra(2, 2) = mesh_glob%np
            mesh_loc%jcc_extra(2) = mesh_glob%me
         ELSE
            mesh_loc%isolated_jjs(1) = 1
            mesh_loc%isolated_interfaces(1, 1) = mesh_glob%sides(1)
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
            mesh_loc%jj(:, 1:mesh_loc%dom_np - 1) = mesh_glob%jj(:, me_start:me_end - 1) - np_start + 1
            mesh_loc%jj(2, mesh_loc%dom_np) = mesh_glob%np
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

            mesh_loc%loc_to_glob(mesh_loc%np - 1) = np_start - 1
            mesh_loc%jj(1, 1) = mesh_loc%np - 1

            mesh_loc%loc_to_glob(mesh_loc%np) = 1
            mesh_loc%jj(1, mesh_loc%me) = mesh_loc%np

            mesh_loc%sides(1) = mesh_glob%sides(2)
            mesh_loc%neighs(1) = mesh_glob%np - 2
            mesh_loc%jjs(1, 1) = mesh_glob%np - 1

            mesh_loc%sides(2) = mesh_glob%sides(1)
            mesh_loc%neighs(2) = mesh_glob%np - 1
            mesh_loc%jjs(1, 2) = mesh_glob%np
         ELSE
            mesh_loc%isolated_jjs(1) = mesh_glob%np
            mesh_loc%isolated_interfaces(1, 1) = mesh_glob%sides(2)
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
            mesh_loc%neighs(1) = mesh_glob%np - 1
            mesh_loc%jjs(1, 1) = mesh_glob%np
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
               mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
               mesh_loc%jj(1, m) = mesh_loc%np
            END IF
            IF (mesh_loc%jj(2, m) < 1) THEN
               mesh_loc%loc_to_glob(mesh_loc%np) = np_start - 1
               mesh_loc%jj(2, m) = mesh_loc%np
            END IF
         END DO

         IF (per_bool) THEN
            mesh_loc%jj_extra = mesh_loc%jj_extra - 1
            mesh_loc%jcc_extra = mesh_loc%jcc_extra - 1
         END IF

      END IF

      CALL MPI_ALLGATHER(mesh_loc%dom_np, 1, MPI_INTEGER, mesh_loc%domnp, 1, &
           MPI_INTEGER, communicator, ierr)
      mesh_loc%disp(1) = 1
      DO n = 1, nb_procs
         mesh_loc%disp(n + 1) = mesh_loc%disp(n) + mesh_loc%domnp(n)
      END DO

      CALL MPI_ALLGATHER(mesh_loc%me, 1, MPI_INTEGER, mesh_loc%domcell, 1, &
           MPI_INTEGER, communicator, ierr)
      mesh_loc%discell(1) = 1
      DO n = 1, nb_procs
         mesh_loc%discell(n + 1) = mesh_loc%discell(n) + mesh_loc%domcell(n)
      END DO
      write(*,*) rank, mesh_loc%jj

   END SUBROUTINE extract_mesh_1d

END MODULE mesh_distribution_1d