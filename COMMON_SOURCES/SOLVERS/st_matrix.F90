!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE st_matrix
!!!>>> Old subroutines removed here, can still be found in SFEMaNS
      ! PUBLIC :: st_aij_csr_glob_block, extract, create_my_ghost, st_aij_csr, tri_jlg, st_aij_csr_loc_block, &
      !   st_csr, st_csr_block, st_csr_bloc, st_aij_csr_glob_block_with_extra_layer
!!!>>> Old subroutines removed here, can still be found in SFEMaNS
   PUBLIC :: extract, extract_through_ghost, create_my_ghost, tri_jlg, st_aij_csr_glob_block_with_extra_layer
   PRIVATE
#include "petsc/finclude/petsc.h"
CONTAINS

   !=========================================================================

   SUBROUTINE create_my_ghost(mesh, LA, ifrom)
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      type(petsc_csr_LA) :: LA
      INTEGER, DIMENSION(:), POINTER :: ifrom
      INTEGER :: kmax, nifrom, start, fin, k
      kmax = SIZE(LA%loc_to_glob, 1)
      nifrom = mesh%np - mesh%dom_np
      ALLOCATE(ifrom(kmax * nifrom))
      IF (nifrom/=0) THEN
         DO k = 1, kmax
            start = (k - 1) * nifrom + 1
            fin = start + nifrom - 1
            ifrom(start:fin) = LA%loc_to_glob(k, mesh%dom_np + 1:) - 1
         END DO
      END IF
   END SUBROUTINE create_my_ghost

   SUBROUTINE extract(xghost, ks, ke, LA, phi)
      !> Subroutine to extract values of xghost(vec_ghost PETSc) to phi(HYPAR)
      !! ks/ke ==> positions inside xghost: if xghost is the concatenation of n components
      !!           then phi will receive components between ks-ke (both included)
#include "petsc/finclude/petscvec.h"
      use petsc
      USE def_type_mesh
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(OUT) :: phi
      type(petsc_csr_LA) :: LA
      INTEGER :: ks, ke
      INTEGER :: k, start, fin, nbghost, s, f
      Vec            :: xghost
      PetscErrorCode :: ierr
      PetscScalar, POINTER :: x_loc(:)
      CALL VecGetArrayF90(xghost, x_loc, ierr)
      DO k = ks, ke
         start = SUM(LA%dom_np(1:k - 1)) + 1
         fin = start + LA%dom_np(k) - 1
         s = SUM(LA%np(ks:k - 1)) + 1
         f = s + LA%dom_np(k) - 1
         phi(s:f) = x_loc(start:fin)
         nbghost = LA%np(k) - LA%dom_np(k)
         start = SUM(LA%dom_np) + SUM(LA%np(1:k - 1) - LA%dom_np(1:k - 1)) + 1
         fin = start + nbghost - 1
         s = SUM(LA%np(ks:k - 1)) + LA%dom_np(k) + 1
         f = s + nbghost - 1
         phi(s:f) = x_loc(start:fin)
      END DO
      CALL VecRestoreArrayF90(xghost, x_loc, ierr)
   END SUBROUTINE extract

   SUBROUTINE extract_through_ghost(xvec, xghost, ks, ke, LA, phi, operation_ghost, opt_assemble)
      !> VB 30/04/2026 => subroutine to simplify readability of code
      !! instead of always having to copy/paste the same lines
      !! see SUBROUTINE extract() for more information

      !! opt_assemble = .TRUE. if need assembling xvec, .FALSE. by default

      USE def_type_mesh, ONLY : petsc_csr_LA
      USE my_util,       ONLY : error_Petsc
#include "petsc/finclude/petsc.h"
      USE petsc
      IMPLICIT NONE
      INTEGER,                      INTENT(IN)  :: ks, ke
      CHARACTER(len=*),             INTENT(IN)  :: operation_ghost
      REAL(KIND = 8), DIMENSION(:), INTENT(OUT) :: phi
      LOGICAL, OPTIONAL                         :: opt_assemble
      TYPE(petsc_csr_LA)                        :: LA    
      INTEGER :: ierr
      Vec :: xvec, xghost

      IF (PRESENT(opt_assemble)) THEN
         IF (opt_assemble) THEN
            CALL VecAssemblyBegin(xvec, ierr)
            CALL VecAssemblyEnd(xvec, ierr)
         END IF
      END IF

      CALL VecGhostGetLocalForm(xvec, xghost, ierr)
      SELECT CASE(operation_ghost)
      CASE('insert')
         CALL VecGhostUpdateBegin(xvec, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(xvec, INSERT_VALUES, SCATTER_FORWARD, ierr)
      CASE('min')
         CALL VecGhostUpdateBegin(xvec, MIN_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(xvec, MIN_VALUES, SCATTER_FORWARD, ierr)
      CASE('max')
         CALL VecGhostUpdateBegin(xvec, MAX_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(xvec, MAX_VALUES, SCATTER_FORWARD, ierr)
      CASE DEFAULT
         CALL error_petsc("unavailable operation_ghost in extract_through_ghost '"//operation_ghost//"'.&
          Available for now => 'insert/min/max'")
      END SELECT
      CALL extract(xghost, ks, ke, LA, phi)
   END SUBROUTINE extract_through_ghost

   SUBROUTINE block_index(communicator, kmax, mesh, loc_to_glob_LA)
      USE def_type_mesh
      USE my_util
#include "petsc/finclude/petsc.h"
      use petsc
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      INTEGER, INTENT(IN) :: kmax
      INTEGER, DIMENSION(:, :), POINTER :: loc_to_glob_LA
      INTEGER, DIMENSION(:), POINTER :: dom_np, disp
      INTEGER :: code, nb_procs, rank
      INTEGER :: i, p, n, k, i_loc, proc, iglob
      MPI_Comm       :: communicator

      CALL MPI_COMM_SIZE(communicator, nb_procs, code)
      CALL MPI_COMM_RANK(communicator, rank, code)
      ALLOCATE(dom_np(nb_procs), disp(nb_procs + 1))
      CALL MPI_ALLGATHER(mesh%dom_np, 1, MPI_INTEGER, dom_np, 1, &
           MPI_INTEGER, communicator, code)
      disp(1) = 1
      DO n = 1, nb_procs
         disp(n + 1) = disp(n) + dom_np(n)
      END DO
      IF (ASSOCIATED(mesh%disp)) THEN
         NULLIFY(mesh%disp)
      END IF
      IF (ASSOCIATED(mesh%domnp)) THEN
         NULLIFY(mesh%domnp)
      END IF
      ALLOCATE(mesh%disp(nb_procs + 1))
      ALLOCATE(mesh%domnp(nb_procs))
      mesh%disp = disp
      mesh%domnp = dom_np

      ALLOCATE(loc_to_glob_LA(kmax, mesh%np))
      proc = rank + 1

      DO i = 1, mesh%dom_np
         DO k = 1, kmax
            loc_to_glob_LA(k, i) = kmax * (disp(proc) - 1) + (k - 1) * dom_np(proc) + i
         END DO
      END DO

      !!$!TEST
      !!$    DO i = 1, mesh%dom_np
      !!$       iglob = mesh%loc_to_glob(i)
      !!$       DO p = 2, nb_procs+1
      !!$          IF (disp(p) > iglob) THEN
      !!$             proc = p - 1
      !!$             EXIT
      !!$          END IF
      !!$       END DO
      !!$       IF (rank+1/=proc) THEN
      !!$          write(*,*) 'BUG2', rank+1, proc
      !!$          STOP
      !!$       END IF
      !!$
      !!$       DO k = 2, kmax
      !!$          IF (loc_to_glob_LA(k,i) - dom_np(proc) /= loc_to_glob_LA(k-1,i)) THEN
      !!$             write(*,*) ' BUG1 ', rank
      !!$             stop
      !!$          END IF
      !!$       END DO
      !!$    END DO
      !!$!TEST

      DO i = mesh%dom_np + 1, mesh%np
         iglob = mesh%loc_to_glob(i)
         DO p = 2, nb_procs + 1
            IF (disp(p) > iglob) THEN
               proc = p - 1
               EXIT
            END IF
         END DO
         i_loc = iglob - disp(proc) + 1
         DO k = 1, kmax
            loc_to_glob_LA(k, i) = kmax * (disp(proc) - 1) + (k - 1) * dom_np(proc) + i_loc
         END DO
      END DO

      !!$!TEST
      !!$    DO i = 1, mesh%np
      !!$       iglob = mesh%loc_to_glob(i)
      !!$       DO p = 2, nb_procs+1
      !!$          IF (disp(p) > iglob) THEN
      !!$             proc = p - 1
      !!$             EXIT
      !!$          END IF
      !!$       END DO
      !!$       DO k = 2, kmax
      !!$          IF (loc_to_glob_LA(k,i) - dom_np(proc) /= loc_to_glob_LA(k-1,i)) THEN
      !!$             write(*,*) ' BUG ', rank
      !!$             stop
      !!$          END IF
      !!$       END DO
      !!$    END DO
      !!$!TEST

      DEALLOCATE(dom_np, disp)

   END SUBROUTINE block_index

   SUBROUTINE st_aij_csr_glob_block_with_extra_layer(communicator, kmax, mesh, LA)
      !> UPDATE 16/04/2026: opt_per is now in mesh%per
      !!  input coefficient structure of the matrix and
      !!  perform the ordering of the unknowns
      !!  jj(nodes_per_element, number_of_elements)
      !!                  --->  node number in the grid
      USE def_type_mesh
      USE my_util
      USE periodic_data_module
      IMPLICIT NONE
      TYPE(mesh_type), INTENT(IN) :: mesh
      INTEGER, INTENT(IN) :: kmax
      TYPE(petsc_csr_LA), INTENT(OUT) :: LA
      INTEGER :: nparm = 200
      INTEGER :: me, nw, nmax, np, knp, ki, kj, k, njt, i1, i2
      INTEGER :: m, ni, nj, i, j, n_a_d, iloc, jloc, iglob, jglob, nb_procs, p, proc = -1
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: ja_work
      INTEGER, DIMENSION(:), ALLOCATABLE :: nja, a_d
      INTEGER, DIMENSION(:), ALLOCATABLE :: per_loc
      LOGICAL :: out

      !#include "petsc/finclude/petsc.h"
      MPI_Comm       :: communicator

      CALL block_index(communicator, kmax, mesh, LA%loc_to_glob)
      nw = SIZE(mesh%jj, 1)
      me = mesh%dom_me
      np = mesh%dom_np
      knp = kmax * np
      nb_procs = SIZE(mesh%domnp)

      LA%kmax = kmax
      ALLOCATE(LA%dom_np(kmax), LA%np(kmax))
      LA%dom_np(:) = mesh%dom_np
      LA%np(:) = mesh%np

      IF (np==0) THEN
         ALLOCATE(LA%ia(0:0), LA%ja(0))
         LA%ia(0) = 0
         RETURN
      END IF
      ALLOCATE (ja_work(knp, kmax * nparm), a_d(kmax * nparm), nja(knp))
      ALLOCATE (per_loc(knp))
      ja_work = 0
      nja = 1
      DO ki = 1, kmax
         DO i = 1, np
            ja_work((ki - 1) * np + i, 1) = LA%loc_to_glob(ki, i)
         END DO
      END DO

      DO m = 1, mesh%me
         DO ni = 1, nw
            iloc = mesh%jj(ni, m)
            IF (iloc>np) CYCLE
            DO ki = 1, kmax
               i = iloc + (ki - 1) * np
               DO nj = 1, nw
                  jloc = mesh%jj(nj, m)
                  jglob = mesh%loc_to_glob(jloc)
                  IF (jloc>np) THEN
                     DO p = 2, nb_procs + 1
                        IF (mesh%disp(p) > jglob) THEN
                           proc = p - 1
                           EXIT
                        END IF
                     END DO
                     out = .TRUE.
                     jloc = jglob - mesh%disp(proc) + 1
                  ELSE
                     out = .FALSE.
                  END IF
                  DO kj = 1, kmax
                     IF (out) THEN
                        j = kmax * (mesh%disp(proc) - 1) + (kj - 1) * mesh%domnp(proc) + jloc
                     ELSE
                        j = LA%loc_to_glob(kj, jloc)
                     END IF
                     IF (MINVAL(ABS(ja_work(i, 1:nja(i)) - j)) /= 0) THEN
                        nja(i) = nja(i) + 1
                        ja_work(i, nja(i)) = j
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      !===Loop over the extra layer
      DO m = 1, mesh%mextra
         DO ni = 1, nw
            iglob = mesh%jj_extra(ni, m)
            IF (iglob<mesh%loc_to_glob(1) .OR. iglob>mesh%loc_to_glob(1) + mesh%dom_np - 1) CYCLE
            iloc = iglob - mesh%loc_to_glob(1) + 1
            DO ki = 1, kmax
               i = iloc + (ki - 1) * np
               DO nj = 1, nw
                  jglob = mesh%jj_extra(nj, m)
                  jloc = jglob - mesh%loc_to_glob(1) + 1
                  IF (jloc<1 .OR. jloc>np) THEN
                     DO p = 2, nb_procs + 1
                        IF (mesh%disp(p) > jglob) THEN
                           proc = p - 1
                           EXIT
                        END IF
                     END DO
                     out = .TRUE.
                     jloc = jglob - mesh%disp(proc) + 1
                  ELSE
                     out = .FALSE.
                  END IF
                  DO kj = 1, kmax
                     IF (out) THEN
                        j = kmax * (mesh%disp(proc) - 1) + (kj - 1) * mesh%domnp(proc) + jloc
                     ELSE
                        j = LA%loc_to_glob(kj, jloc)
                     END IF
                     IF (MINVAL(ABS(ja_work(i, 1:nja(i)) - j)) /= 0) THEN
                        nja(i) = nja(i) + 1
                        ja_work(i, nja(i)) = j
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO
      IF (mesh%per%nb_bords /= 0) THEN
         DO k = 1, mesh%per%nb_bords
            DO i = 1, SIZE(mesh%per%list(k)%DIL)
               per_loc = 0
               i1 = mesh%per%list(k)%DIL(i)
               i2 = mesh%per%perlist(k)%DIL(i)
               njt = nja(i1) + nja(i2)
               IF (njt > kmax * nparm) THEN
                  CALL error_Petsc('BUG in st_aij_glob_block, SIZE(ja) not large enough')
               END IF
               per_loc(1:nja(i1)) = ja_work(i1, 1:nja(i1))
               per_loc(nja(i1) + 1:nja(i1) + nja(i2)) = ja_work(i2, 1:nja(i2))
               nja(i1) = njt
               nja(i2) = njt
               ja_work(i1, 1:njt) = per_loc(1:njt)
               ja_work(i2, 1:njt) = per_loc(1:njt)
            END DO
         END DO
      END IF

      IF (MAXVAL(nja)>nparm) THEN
         WRITE(*, *) 'ST_SPARSEKIT: dimension de ja doit etre >= ', nparm
         STOP
      END IF
      nmax = 0
      DO i = 1, knp
         nmax = nmax + nja(i)
      END DO
      ALLOCATE(LA%ia(0:knp), LA%ja(0:nmax - 1))
      LA%ia(0) = 0
      DO i = 1, knp
         CALL tri_jlg (ja_work(i, 1:nja(i)), a_d, n_a_d)
         IF (n_a_d /= nja(i)) THEN
            WRITE(*, *) ' BUG : st_p1_CSR'
            WRITE(*, *) 'n_a_d ', n_a_d, 'nja(i)', nja(i)
            STOP
         END IF
         LA%ia(i) = LA%ia(i - 1) + nja(i)
         LA%ja(LA%ia(i - 1):LA%ia(i) - 1) = a_d(1:nja(i)) - 1
      END DO
      DEALLOCATE (ja_work, nja, a_d)
      DEALLOCATE (per_loc)
   END SUBROUTINE st_aij_csr_glob_block_with_extra_layer

   SUBROUTINE tri_jlg (a, a_d, n_a_d)
      !  sort in ascending order of the integer array  a  and generation
      !  of the integer array  a_d  whose first  n_a_d  leading entries
      !  contain different values in ascending order, while all the
      !  remaining entries are set to zero
      !  sorting by Shell's method.
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: a
      INTEGER, DIMENSION(:), INTENT(OUT) :: a_d
      INTEGER, INTENT(OUT) :: n_a_d
      INTEGER :: n, na, inc, i, j, k, ia

      na = SIZE(a)

      !  sort phase

      IF (na == 0) THEN
         n_a_d = 0
         RETURN
      ENDIF

      inc = 1
      DO WHILE (inc <= na)
         inc = inc * 3
         inc = inc + 1
      ENDDO

      DO WHILE (inc > 1)
         inc = inc / 3
         DO i = inc + 1, na
            ia = a(i)
            j = i
            DO WHILE (a(j - inc) > ia)
               a(j) = a(j - inc)
               j = j - inc
               IF (j <= inc) EXIT
            ENDDO
            a(j) = ia
         ENDDO
      ENDDO

      !  compression phase

      n = 1
      a_d(n) = a(1)
      DO k = 2, na
         IF (a(k) > a(k - 1)) THEN
            n = n + 1
            a_d(n) = a(k)
            !TEST JUIN 13 2008
         ELSE
            WRITE(*, *) 'We have a problem in the compression phase of tri_jlg', a(k), a(k - 1)
            !TEST JUIN 13 2008
         ENDIF
      ENDDO

      n_a_d = n

      a_d(n_a_d + 1:na) = 0

   END SUBROUTINE tri_jlg

END MODULE st_matrix
