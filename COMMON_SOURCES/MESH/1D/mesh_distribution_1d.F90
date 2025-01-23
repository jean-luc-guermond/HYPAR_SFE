MODULE mesh_distribution_1d
#include "petsc/finclude/petsc.h"
   USE petsc
   USE mesh_tools
   PUBLIC :: extract_mesh_1d
   PRIVATE
CONTAINS
   SUBROUTINE extract_mesh_1d(communicator, mesh_glob, mesh_loc)
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_glob, mesh_loc
      INTEGER :: n, m, np_start, np_end, me_start, me_end
      PetscErrorCode :: ierr
      PetscMPIInt    :: rank, nb_procs
      MPI_Comm       :: communicator
      CALL MPI_Comm_rank(communicator, rank, ierr)
      CALL MPI_COMM_SIZE(communicator, nb_procs, ierr)
      rank = rank + 1

      IF  (nb_procs == 1) THEN
         CALL copy_mesh(mesh_glob, mesh_loc)
      END IF

      mesh_loc%mi = 0
      mesh_loc%medge = 0
      mesh_loc%medges = 0
      mesh_loc%mes_extra = 0
      mesh_loc%mes_int = 0
      mesh_loc%mes = 0
      mesh_loc%edge_stab = mesh_glob%edge_stab

      IF (rank < nb_procs) THEN
         np_start = (rank - 1) * (mesh_glob%np / nb_procs) + 1
         np_end = rank * mesh_glob%np / nb_procs
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
      ELSE IF (rank == nb_procs) THEN
         mesh_loc%np = mesh_loc%dom_np + 1
         mesh_loc%nis = 1
      ELSE
         mesh_loc%np = mesh_loc%dom_np + 1
         mesh_loc%nis = 0
      END IF

      mesh_loc%dom_me = mesh_loc%me
      mesh_loc%dom_mes = mesh_loc%mes
      ALLOCATE(mesh_loc%jj(2, mesh_loc%me), mesh_loc%jjs(1, mesh_loc%mes), mesh_loc%iis(0, 0))
      ALLOCATE(mesh_loc%jj_extra(2, mesh_loc%mextra), mesh_loc%jce_extra(0, mesh_loc%medge), &
           mesh_loc%jjs_extra(0, mesh_loc%mes_extra))
      ALLOCATE(mesh_loc%jjs_int(0, 0), mesh_loc%jcc_extra(mesh_loc%mextra), mesh_loc%jce(0, 0))
      ALLOCATE(mesh_loc%jees(0), mesh_loc%jecs(0))
      ALLOCATE(mesh_loc%jji(0, 0, 0), mesh_loc%jjsi(0, 0), mesh_loc%j_s(0))
      ALLOCATE(mesh_loc%rr(2, mesh_loc%np), mesh_loc%rrs_extra(2, 2, 0))
      ALLOCATE(mesh_loc%neigh(2, mesh_loc%me), mesh_loc%neighi(0, 0))
      ALLOCATE(mesh_loc%sides(mesh_loc%mes), mesh_loc%neighs(mesh_loc%mes))
      ALLOCATE(mesh_loc%sides_extra(mesh_loc%mes_extra), mesh_loc%neighs_extra(mesh_loc%mes_extra))
      ALLOCATE(mesh_loc%sides_int(mesh_loc%mes_int), mesh_loc%neighs_int(2, mesh_loc%mes_int))
      ALLOCATE(mesh_loc%i_d(mesh_loc%me), mesh_loc%loc_to_glob(mesh_loc%np))
      ALLOCATE(mesh_loc%disp(nb_procs + 1), mesh_loc%disedge(nb_procs + 1), mesh_loc%discell(nb_procs + 1))
      ALLOCATE(mesh_loc%domnp(nb_procs), mesh_loc%domedge(nb_procs), mesh_loc%domcell(nb_procs))
      ALLOCATE(mesh_loc%isolated_jjs(mesh_loc%nis), mesh_loc%isolated_interfaces(mesh_loc%nis, 1))

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

      DO n = 1, mesh_loc%dom_np
         mesh_loc%loc_to_glob(n) = np_start - 1 + n
      END DO

      mesh_loc%i_d = mesh_glob%i_d(me_start:me_end)
      mesh_loc%jj = mesh_glob%jj(:, me_start:me_end) - np_start + 1
      mesh_loc%rr(:, 1:mesh_loc%dom_np) = mesh_glob%rr(:, np_start:np_end)

      IF (rank == 1) THEN
         mesh_loc%isolated_jjs(1) = 1
         mesh_loc%isolated_interfaces(1, 1) = mesh_glob%sides(1)
         mesh_loc%jj_extra(:, 1) = mesh_glob%jj(:, me_end + 1)
         mesh_loc%jcc_extra = me_end + 1
      ELSE IF (rank == nb_procs) THEN
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
      ELSE
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
      END IF

   END SUBROUTINE extract_mesh_1d

END MODULE mesh_distribution_1d