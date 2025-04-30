MODULE two_dim_metis_distribution
   !#include "petsc/finclude/petsc.h"
   USE petsc
   PUBLIC :: reorder_mesh
   PRIVATE
   REAL(KIND = 8) :: epsilon = 1.d-10
   !!$ Dummy for metis...
   INTEGER :: METIS_NOPTIONS = 40, METIS_OPTION_NUMBERING = 18
   !!$ Dummy for metis...
CONTAINS
   SUBROUTINE reorder_mesh(communicator, nb_proc, mesh, mesh_loc, list_of_interfaces)
      USE def_type_mesh
      USE my_util
      USE sub_plot
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh, mesh_loc
      INTEGER, DIMENSION(2) :: me_loc, mes_loc, np_loc
      INTEGER, DIMENSION(:), POINTER, OPTIONAL :: list_of_interfaces
      INTEGER, DIMENSION(mesh%me + 1) :: xind
      INTEGER, DIMENSION(mesh%me) :: vwgt, vsize, old_m_to_new, cellnum
      REAL(KIND = 8), DIMENSION(mesh%me) :: plot
      INTEGER, DIMENSION(mesh%me) :: part
      LOGICAL, DIMENSION(mesh%np) :: virgin
      INTEGER, DIMENSION(mesh%np) :: old_j_to_new
      LOGICAL, DIMENSION(mesh%medge) :: virginss
      INTEGER, DIMENSION(mesh%medge) :: old_edge_to_new
      INTEGER, DIMENSION(mesh%mes) :: old_ms_to_new, parts
      INTEGER, DIMENSION(2, mesh%mes) :: inter_news
      INTEGER, DIMENSION(SIZE(mesh%jjs, 1)) :: i_loc
      LOGICAL, DIMENSION(mesh%mes) :: virgins
      REAL(KIND = 8), DIMENSION(mesh%np) :: r_old_j_to_new
      INTEGER, DIMENSION(5) :: opts
      INTEGER, DIMENSION(:), ALLOCATABLE :: xadj, adjwgt, nblmt_per_proc, start, displ
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: inter
      INTEGER :: dim, nb_neigh, edge, m, ms, n, nb, numflag, p, wgtflag, new_dof, j, &
           news, ns, nws, msop, nsop, proc, iop, ncon, k
      LOGICAL :: test = .FALSE.
      !REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: tpwgts !Up to petsc.3.7.7
      !REAL(KIND=8), DIMENSION(1)               :: ubvec !Up to petsc.3.7.7
      REAL(KIND = 4), DIMENSION(:), ALLOCATABLE :: tpwgts !(JLG, Feb 20, 2019) Up from petsc.3.8.4,
      REAL(KIND = 4), DIMENSION(1) :: ubvec  !petsc peoples changed metis types
      INTEGER, DIMENSION(METIS_NOPTIONS) :: metis_opt
      CHARACTER(len = 4) :: tit
      !LOGICAL, SAVE :: once = .true.

#include "petsc/finclude/petsc.h"
      PetscErrorCode :: ierr
      PetscMPIInt    :: rank, nb_proc
      MPI_Comm       :: communicator
      CALL MPI_Comm_rank(communicator, rank, ierr)
      !CALL MPI_Comm_size(PETSC_COMM_WORLD,nb_proc,ierr)

      ALLOCATE(mesh_loc%disp(nb_proc + 1), mesh_loc%domnp(nb_proc))
      ALLOCATE(mesh_loc%discell(nb_proc + 1), mesh_loc%domcell(nb_proc))
      ALLOCATE(mesh_loc%disedge(nb_proc + 1), mesh_loc%domedge(nb_proc))
      !IF (once) THEN
      !   ALLOCATE(part(mesh%me))
      !END IF

!      IF (nb_proc==1) THEN
!         me_loc = (/ 1, mesh%me /)
!         mes_loc = (/ 1, mesh%mes /)
!         np_loc = (/ 1, mesh%np /)
!         !CALL create_local_mesh(mesh, mesh_loc, me_loc, mes_loc, np_loc)
!         CALL create_local_mesh_with_extra_layer(mesh, mesh_loc, me_loc, mes_loc, np_loc)
!         RETURN
!      END IF
!
!      ! Create the connectivity array based on neigh
!      nb_neigh = SIZE(mesh%neigh, 1)
!      xind(1) = 1
!      DO m = 1, mesh%me
!         nb = 0
!         DO n = 1, nb_neigh
!            IF (mesh%neigh(n, m)==0) CYCLE
!            nb = nb + 1
!         END DO
!         xind(m + 1) = xind(m) + nb
!      END DO
!      ALLOCATE(xadj(xind(mesh%me + 1) - 1))
!      p = 0
!      DO m = 1, mesh%me
!         DO n = 1, nb_neigh
!            IF (mesh%neigh(n, m)==0) CYCLE
!            p = p + 1
!            xadj(p) = mesh%neigh(n, m)
!         END DO
!      END DO
!      IF (p/=xind(mesh%me + 1) - 1) THEN
!         CALL error_Petsc('BUG, p/=xind(mesh%me+1)-1')
!      END IF
!      ! End create the connectivity array based on neigh
!
!      ALLOCATE(adjwgt(SIZE(xadj)))
!      opts = 0
!      adjwgt = 1
!      numflag = 1 ! Fortran numbering of processors
!      wgtflag = 2
!      vwgt = 1
!      vsize = 1
!      ncon = 1
!      ALLOCATE(tpwgts(nb_proc))
!      tpwgts = 1.d0 / nb_proc
!      CALL METIS_SetDefaultOptions(metis_opt)
!      metis_opt(METIS_OPTION_NUMBERING) = 1
!      ubvec = 1.01
!      !IF (once) THEN
!      CALL METIS_PartGraphRecursive(mesh%me, ncon, xind(1:), xadj(1:), vwgt(1:), vsize(1:), adjwgt(1:), &
!           nb_proc, tpwgts(1:), ubvec(1:), metis_opt(1:), edge, part(1:))
!      !   once = .false.
!      !END IF
!
!      IF (rank==0) THEN
!         WRITE(tit, '(i4)') SIZE(mesh%jj, 1)
!         plot = 0.d0
!         do m = 1, mesh%me
!            if (part(m) == 150) then
!               plot(m) = 1.d0
!            else if (part(m) == 152)  then
!               plot(m) = 2.d0
!            else if (part(m) == 153)  then
!               plot(m) = 3.d0
!            end if
!         end do
!         CALL plot_const_p1_label(mesh%jj, mesh%rr, plot, 'dd' // tit // '.plt')
!      END IF
!
!      ! Count elements on processors
!      ALLOCATE(nblmt_per_proc(nb_proc), start(nb_proc), displ(nb_proc))
!      nblmt_per_proc = 0
!      DO m = 1, mesh%me
!         nblmt_per_proc(part(m)) = nblmt_per_proc(part(m)) + 1
!      END DO
!      start(1) = 0
!      DO n = 2, nb_proc
!         start(n) = start(n - 1) + nblmt_per_proc(n - 1)
!      END DO
!      me_loc(1) = start(rank + 1) + 1
!      me_loc(2) = start(rank + 1) + nblmt_per_proc(rank + 1)
!      displ = start
!      ! End count elements on processors
!
!      ! Re-order elements
!      DO m = 1, mesh%me
!         start(part(m)) = start(part(m)) + 1
!         old_m_to_new(m) = start(part(m))
!      END DO
!      ! Re-order elements
!
!
!      !==Search on the boundary whether ms is on a cut.
!      nws = SIZE(mesh%jjs, 1)
!      news = 0
!      inter_news = 0
!      parts = part(mesh%neighs)
!      IF (PRESENT(list_of_interfaces)) THEN !==There is an interface
!         IF (SIZE(list_of_interfaces)/=0) THEN
!            virgins = .TRUE.
!            news = 0
!            DO ms = 1, mesh%mes
!               IF (.NOT.virgins(ms)) CYCLE
!               IF (MINVAL(ABS(mesh%sides(ms) - list_of_interfaces))/=0) CYCLE !==ms not on a cut
!               i_loc = mesh%jjs(:, ms)
!               DO msop = 1, mesh%mes
!                  IF (msop==ms .OR. .NOT.virgins(msop)) CYCLE
!                  IF (MINVAL(ABS(mesh%sides(msop) - list_of_interfaces))/=0) CYCLE !==msop not on a cut
!                  DO ns = 1, nws
!                     test = .FALSE.
!                     DO nsop = 1, nws
!                        iop = mesh%jjs(nsop, msop)
!                        IF (MAXVAL(ABS(mesh%rr(:, i_loc(ns)) - mesh%rr(:, iop))).LT.epsilon) THEN
!                           test = .TRUE.
!                           EXIT
!                        END IF
!                     END DO
!                     IF (.NOT.test) THEN
!                        EXIT !==This msop does not coincide with ms
!                     END IF
!                  END DO
!                  IF (test) EXIT
!               END DO
!               IF (.NOT.test) THEN
!                  CALL error_Petsc('BUG in create_local_mesh, .NOT.test ')
!               END IF
!               IF (part(mesh%neighs(ms)) == part(mesh%neighs(msop))) CYCLE !==ms is an internal cut
!               proc = MIN(part(mesh%neighs(ms)), part(mesh%neighs(msop)))
!               parts(ms) = proc
!               parts(msop) = proc
!               virgins(ms) = .FALSE.
!               virgins(msop) = .FALSE.
!               IF (proc /= rank + 1) CYCLE !==ms and msop do not touch the current proc
!               news = news + 1
!               inter_news(1, news) = ms
!               inter_news(2, news) = msop
!            END DO
!         END IF
!      END IF
!
!
!      ! Re-order jj
!      ALLOCATE(inter(SIZE(mesh%jj, 1), mesh%me))
!      DO m = 1, mesh%me
!         inter(:, old_m_to_new(m)) = mesh%jj(:, m)
!      END DO
!      mesh%jj = inter
!
!      virgin = .TRUE.
!      new_dof = 0
!      DO k = 1, nb_proc
!         DO m = displ(k) + 1, displ(k) + nblmt_per_proc(k)
!            DO n = 1, 3
!               j = mesh%jj(n, m)
!               IF (.NOT.virgin(j)) CYCLE
!               new_dof = new_dof + 1
!               virgin(j) = .FALSE.
!               old_j_to_new(j) = new_dof
!            END DO
!         END DO
!
!         IF (SIZE(mesh%jj, 1) == 6) THEn
!            DO m = displ(k) + 1, displ(k) + nblmt_per_proc(k)
!               DO n = 4, 6
!                  j = mesh%jj(n, m)
!                  IF (.NOT.virgin(j)) CYCLE
!                  new_dof = new_dof + 1
!                  virgin(j) = .FALSE.
!                  old_j_to_new(j) = new_dof
!               END DO
!            END DO
!         END IF
!         IF (SIZE(mesh%jj, 1) == 10) THEn
!            DO m = displ(k) + 1, displ(k) + nblmt_per_proc(k)
!               DO n = 4, 9
!                  j = mesh%jj(n, m)
!                  IF (.NOT.virgin(j)) CYCLE
!                  new_dof = new_dof + 1
!                  virgin(j) = .FALSE.
!                  old_j_to_new(j) = new_dof
!               END DO
!            END DO
!            DO m = displ(k) + 1, displ(k) + nblmt_per_proc(k)
!               j = mesh%jj(10, m)
!               IF (.NOT.virgin(j)) CYCLE
!               new_dof = new_dof + 1
!               virgin(j) = .FALSE.
!               old_j_to_new(j) = new_dof
!            END DO
!         END IF
!      END DO
!
!      DO m = 1, mesh%me
!         inter(:, m) = old_j_to_new(mesh%jj(:, m))
!      END DO
!      mesh%jj = inter
!      DEALLOCATE(inter)
!
!      IF (rank == 0) THEN
!         np_loc(1) = 1
!      ELSE
!         np_loc(1) = MAXVAL(mesh%jj(:, displ(rank) + 1:displ(rank) + nblmt_per_proc(rank))) + 1
!      END IF
!      np_loc(2) = MAXVAL(mesh%jj(:, me_loc(1):me_loc(2)))
!      ! End re-order jj
!
!      ! Re-order edge
!      ALLOCATE(inter(SIZE(mesh%jce, 1), mesh%me))
!      DO m = 1, mesh%me
!         inter(:, old_m_to_new(m)) = mesh%jce(:, m)
!      END DO
!      mesh%jce = inter
!
!      virginss = .TRUE.
!      new_dof = 0
!      DO m = 1, mesh%me
!         DO n = 1, SIZE(mesh%jce, 1)
!            j = mesh%jce(n, m)
!            IF (.NOT.virginss(j)) CYCLE
!            new_dof = new_dof + 1
!            virginss(j) = .FALSE.
!            old_edge_to_new(j) = new_dof
!         END DO
!      END DO
!
!      DO m = 1, mesh%me
!         inter(:, m) = old_edge_to_new(mesh%jce(:, m))
!      END DO
!
!      mesh%jce = inter
!      DEALLOCATE(inter)
!      ! End re-order edge
!
!      ! Re-order rr
!      DO  n = 1, SIZE(mesh%rr, 1)
!         r_old_j_to_new(old_j_to_new) = mesh%rr(n, :)
!         mesh%rr(n, :) = r_old_j_to_new(:)
!      END DO
!      ! Re-order rr
!
!      ! Re-order neigh
!      ALLOCATE(inter(SIZE(mesh%neigh, 1), mesh%me))
!      dim = SIZE(mesh%rr, 1)
!      DO m = 1, mesh%me
!         DO n = 1, dim + 1
!            IF (mesh%neigh(n, m) /=0) THEN
!               inter(n, old_m_to_new(m)) = old_m_to_new(mesh%neigh(n, m))
!            ELSE
!               inter(n, old_m_to_new(m)) = 0
!            END IF
!         END DO
!      END DO
!      mesh%neigh = inter
!      DEALLOCATE(inter)
!      ! End re-order neigh
!
!      ! Re-order i_d
!      DEALLOCATE(xadj); ALLOCATE(xadj(mesh%me))
!      xadj(old_m_to_new) = mesh%i_d
!      mesh%i_d = xadj
!      ! End Re-order i_d
!
!      ! Re-order jev
!      !      ALLOCATE(inter(SIZE(mesh%jev, 1), mesh%medge))
!      !      dim = SIZE(mesh%jev, 1)
!      !      DO m = 1, mesh%medge
!      !         DO n = 1, dim
!      !            inter(n, old_edge_to_new(m)) = old_j_to_new(mesh%jev(n, m))
!      !         END DO
!      !      END DO
!      !      mesh%jev = inter
!      !      DEALLOCATE(inter)
!      ! End Re-order jev
!
!      ! Re-order neighs
!      DEALLOCATE(xadj); ALLOCATE(xadj(mesh%mes))
!
!      nblmt_per_proc = 0
!      DO ms = 1, mesh%mes
!         n = parts(ms)
!         nblmt_per_proc(n) = nblmt_per_proc(n) + 1
!      END DO
!      start(1) = 0
!      DO n = 2, nb_proc
!         start(n) = start(n - 1) + nblmt_per_proc(n - 1)
!      END DO
!      mes_loc(1) = start(rank + 1) + 1
!      mes_loc(2) = start(rank + 1) + nblmt_per_proc(rank + 1)
!
!      DO ms = 1, mesh%mes
!         n = parts(ms)
!         start(n) = start(n) + 1
!         old_ms_to_new(ms) = start(n)
!      END DO
!      xadj(old_ms_to_new) = mesh%neighs
!      mesh%neighs = xadj
!      xadj = old_m_to_new(mesh%neighs)
!      mesh%neighs = xadj
!      ! End re-order neighs
!
!      ! Re-order inter_news
!      xadj(1:news) = old_ms_to_new(inter_news(1, 1:news))
!      inter_news(1, 1:news) = xadj(1:news)
!      xadj(1:news) = old_ms_to_new(inter_news(2, 1:news))
!      inter_news(2, 1:news) = xadj(1:news)
!      ! End re-order inter_news
!
!      ! Re-order sides
!      xadj(old_ms_to_new) = mesh%sides
!      mesh%sides = xadj
!      ! End re-order sides
!
!      ! Re-order jjs
!      DO n = 1, SIZE(mesh%jjs, 1)
!         xadj(old_ms_to_new) = old_j_to_new(mesh%jjs(n, :))
!         mesh%jjs(n, :) = xadj
!      END DO
!      ! End re-order jjs
!
!      !==We create the local mesh now
!      !CALL create_local_mesh(mesh, mesh_loc, me_loc, mes_loc, np_loc, news, inter_news(:,1:news))
!
!      CALL create_local_mesh_with_extra_layer(mesh, mesh_loc, me_loc, mes_loc, np_loc, news, inter_news(:, 1:news))
!      CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

      DEALLOCATE(xadj, adjwgt, nblmt_per_proc, start, displ, tpwgts)
   END SUBROUTINE reorder_mesh


   SUBROUTINE create_local_mesh_with_extra_layer(mesh, mesh_loc, me_loc, mes_loc, np_loc, news, inter_news)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh, mesh_loc
      INTEGER, DIMENSION(2), INTENT(IN) :: me_loc, mes_loc, np_loc
      INTEGER, DIMENSION(:, :), INTENT(IN), OPTIONAL :: inter_news
      INTEGER, OPTIONAL :: news
      INTEGER, DIMENSION(mesh%me) :: m_glob_to_loc, m_loc_to_glob
      INTEGER, DIMENSION(mesh%np) :: glob_to_loc, loc_to_glob
      LOGICAL, DIMENSION(mesh%np) :: virgin
      LOGICAL, DIMENSION(mesh%medge) :: virgins
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: virginss
      LOGICAL, DIMENSION(mesh%me) :: not_my_cells
      INTEGER, DIMENSION(SIZE(mesh%jj, 1)) :: jglob, eglob
      LOGICAL :: test
      INTEGER :: dim, nws, nw, m, ms, mop, msop, ns, msup, minf, dof, proc, &
           dom_me, nwc, dom_mes, dom_np, n, i, rank, ierr, dom_np_glob, nb_extra, nb_proc, e_glob, medge, medges, j

      dim = SIZE(mesh%rr, 1)
      nws = SIZE(mesh%jjs, 1)
      nw = SIZE(mesh%jj, 1)
      nwc = SIZE(mesh%neigh, 1)
      nb_proc = SIZE(mesh_loc%domnp)

      !==Test if one proc only
      IF (me_loc(2) - me_loc(1) + 1==mesh%me) THEN
         mesh_loc%me = mesh%me
         mesh_loc%np = mesh%np
         mesh_loc%mes = mesh%mes
         mesh_loc%dom_me = mesh%me
         mesh_loc%dom_np = mesh%np
         mesh_loc%dom_mes = mesh%mes
         mesh_loc%mextra = 0
         mesh_loc%medge = mesh%medge
         mesh_loc%medges = 0
         mesh_loc%nis = 0

         ALLOCATE(mesh_loc%jees(mesh_loc%medges))
         ALLOCATE(mesh_loc%jecs(mesh_loc%medges))

         ALLOCATE(mesh_loc%jj(nw, mesh%me))
         mesh_loc%jj = mesh%jj
         ALLOCATE(mesh_loc%loc_to_glob(mesh%np))
         DO n = 1, mesh%np
            mesh_loc%loc_to_glob(n) = n
         END DO
         ALLOCATE(mesh_loc%rr(dim, mesh%np))
         mesh_loc%rr = mesh%rr
         ALLOCATE(mesh_loc%neigh(nwc, mesh%me))
         mesh_loc%neigh = mesh%neigh
         ALLOCATE(mesh_loc%i_d(mesh%me))
         mesh_loc%i_d = mesh%i_d
         ALLOCATE(mesh_loc%neighs(mesh_loc%mes))
         mesh_loc%neighs = mesh%neighs
         ALLOCATE(mesh_loc%sides(mesh_loc%mes))
         mesh_loc%sides = mesh%sides
         ALLOCATE(mesh_loc%jjs(nws, mesh_loc%mes))
         mesh_loc%jjs = mesh%jjs

         ALLOCATE(mesh_loc%extra_jj(nw, mesh_loc%mextra))
         ALLOCATE(mesh_loc%extra_jce(nw, mesh_loc%mextra))
         ALLOCATE(mesh_loc%extra_jcc(mesh_loc%mextra))

         ALLOCATE(mesh_loc%isolated_interfaces(mesh_loc%nis, 2))
         ALLOCATE(mesh_loc%isolated_jjs(mesh_loc%nis))
         ALLOCATE(mesh_loc%jce(SIZE(mesh%jce, 1), mesh%me))
         mesh_loc%jce = mesh%jce
         !ALLOCATE(mesh_loc%jev(SIZE(mesh%jev, 1), mesh%medge))
         !mesh_loc%jev = mesh%jev

         mesh_loc%disp = (/ 1, mesh%np + 1 /)
         mesh_loc%domnp = (/ mesh%np /)
         mesh_loc%discell = (/ 1, mesh%me + 1 /)
         mesh_loc%domcell = (/ mesh%me /)
         mesh_loc%disedge = (/ 1, mesh%medge + 1 /)
         mesh_loc%domedge = (/ mesh%medge /)
         RETURN
      END IF
      !==End test if one proc only

      IF (.NOT.PRESENT(news) .OR. .NOT.PRESENT(inter_news)) THEN
         CALL error_Petsc('BUG in create_local_mesh .NOT.present(news) .OR. .NOT.present( inter_news)')
      END IF

      !==Create the new mesh
      dom_me = me_loc(2) - me_loc(1) + 1
      dom_mes = mes_loc(2) - mes_loc(1) + 1
      dom_np = np_loc(2) - np_loc(1) + 1
      mesh_loc%me = dom_me + news
      mesh_loc%mes = dom_mes
      mesh_loc%dom_me = dom_me
      mesh_loc%dom_np = dom_np
      mesh_loc%dom_mes = dom_mes
      CALL MPI_ALLREDUCE(dom_np, dom_np_glob, 1, MPI_INTEGER, &
           MPI_MIN, PETSC_COMM_WORLD, ierr)
      IF (dom_np_glob.LE.0) THEN
         CALL error_petsc('Pb in create_local_mesh, not enough cells per processors')
      END IF

      CALL MPI_ALLGATHER(mesh_loc%dom_np, 1, MPI_INTEGER, mesh_loc%domnp, 1, &
           MPI_INTEGER, PETSC_COMM_WORLD, ierr)
      mesh_loc%disp(1) = 1
      DO n = 1, nb_proc
         mesh_loc%disp(n + 1) = mesh_loc%disp(n) + mesh_loc%domnp(n)
      END DO

      CALL MPI_ALLGATHER(mesh_loc%me, 1, MPI_INTEGER, mesh_loc%domcell, 1, &
           MPI_INTEGER, PETSC_COMM_WORLD, ierr)
      mesh_loc%discell(1) = 1
      DO n = 1, nb_proc
         mesh_loc%discell(n + 1) = mesh_loc%discell(n) + mesh_loc%domcell(n)
      END DO

      !==Re-order jj
      virgin = .TRUE.
      dof = 0
      DO m = me_loc(1), me_loc(2)
         DO n = 1, nw
            i = mesh%jj(n, m)
            IF(.NOT.virgin(i) .OR. i.GE.np_loc(1)) CYCLE
            virgin(i) = .FALSE.
            dof = dof + 1
         END DO
      END DO
      ALLOCATE(mesh_loc%jj(nw, mesh_loc%me))

      m_glob_to_loc = 0
      virgin = .TRUE.
      dof = dom_np
      DO m = me_loc(1), me_loc(2)
         DO n = 1, nw
            i = mesh%jj(n, m)
            IF(virgin(i)) THEN
               IF (i .LT. np_loc(1)) THEN
                  IF (n<=3) THEN
                     virgin(i) = .FALSE.
                     dof = dof + 1
                     glob_to_loc(i) = dof
                     loc_to_glob(dof) = i
                  END IF
               ELSE
                  virgin(i) = .FALSE.
                  glob_to_loc(i) = i - np_loc(1) + 1
                  loc_to_glob(i - np_loc(1) + 1) = i
               END IF
            END IF
         END DO
         m_loc_to_glob(m - me_loc(1) + 1) = m
         m_glob_to_loc(m) = m - me_loc(1) + 1
      END DO

      IF (SIZE(mesh%jj, 1) == 6) THEn
         DO m = me_loc(1), me_loc(2)
            DO n = 4, 6
               j = mesh%jj(n, m)
               IF (.NOT.virgin(j)) CYCLE
               IF (j .LT. np_loc(1)) THEN
                  dof = dof + 1
                  virgin(j) = .FALSE.
                  glob_to_loc(j) = dof
                  loc_to_glob(dof) = j
               END IF
            END DO
         END DO
      END IF
      IF (SIZE(mesh%jj, 1) == 10) THEn
         DO m = me_loc(1), me_loc(2)
            DO n = 4, 9
               j = mesh%jj(n, m)
               IF (.NOT.virgin(j)) CYCLE
               IF (j .LT. np_loc(1)) THEN
                  dof = dof + 1
                  virgin(j) = .FALSE.
                  glob_to_loc(j) = dof
                  loc_to_glob(dof) = j
               END IF
            END DO
         END DO
         DO m = me_loc(1), me_loc(2)
            j = mesh%jj(10, m)
            IF (.NOT.virgin(j)) CYCLE
            IF (j .GT. np_loc(1)) THEN
               dof = dof + 1
               virgin(j) = .FALSE.
               glob_to_loc(j) = dof
               loc_to_glob(dof) = j
            END IF
         END DO
      END IF

      DO n = 1, nw
         mesh_loc%jj(n, 1:dom_me) = glob_to_loc(mesh%jj(n, me_loc(1):me_loc(2)))
      END DO

      DO ns = 1, news
         ms = inter_news(1, ns)
         msop = inter_news(2, ns)
         IF (mesh%neighs(ms) < me_loc(1) .OR. mesh%neighs(ms) > me_loc(2)) THEN
            m = mesh%neighs(ms)
         ELSE
            m = mesh%neighs(msop)
         END IF

         DO n = 1, nw
            i = mesh%jj(n, m)
            IF(virgin(i)) THEN
               virgin(i) = .FALSE.
               IF (i.GE.np_loc(1) .AND. i.LE.np_loc(2)) THEN
                  CALL error_Petsc('BUG in create_local_mesh, i.GE.np_loc(1) .AND. i.LE.np_loc(2)')
               END IF
               dof = dof + 1
               glob_to_loc(i) = dof
               loc_to_glob(dof) = i
            END IF
         END DO
         mesh_loc%jj(:, dom_me + ns) = glob_to_loc(mesh%jj(:, m))
         m_loc_to_glob(dom_me + ns) = m
         m_glob_to_loc(m) = dom_me + ns
      END DO
      !==End re-order jj

      !==Create mesh%loc_to_glob
      IF (MAXVAL(mesh_loc%jj)/=dof) THEN
         CALL error_Petsc('BUG in create_local_mesh, mesh_loc%jj)/=dof')
      END IF
      mesh_loc%np = dof
      ALLOCATE(mesh_loc%loc_to_glob(mesh_loc%np))
      mesh_loc%loc_to_glob = loc_to_glob(1:mesh_loc%np)
      !==End create mesh%loc_to_glob

      !==Re-order rr
      ALLOCATE(mesh_loc%rr(dim, mesh_loc%np))
      DO  n = 1, mesh_loc%np
         mesh_loc%rr(:, n) = mesh%rr(:, mesh_loc%loc_to_glob(n))
      END DO
      !==End re-order rr

      !==Re-order neigh
      not_my_cells = .TRUE. !JLG Aug 18 2015
      not_my_cells(m_loc_to_glob(1:mesh_loc%me)) = .FALSE. !JLG Aug 18 2015
      ALLOCATE(mesh_loc%neigh(nwc, mesh_loc%me))
      msup = MAXVAL(m_loc_to_glob)
      minf = MINVAL(m_loc_to_glob)
      DO m = 1, mesh_loc%me
         DO n = 1, nwc
            mop = mesh%neigh(n, m_loc_to_glob(m))
            IF (mop==0) THEN
               mesh_loc%neigh(n, m) = 0
               !ELSE IF(mop<minf .OR. mop>msup) THEN
            ELSE IF (not_my_cells(mop)) THEN !JLG Aug 18 2015
               mesh_loc%neigh(n, m) = -1 !JLG Aug 13 2015
            ELSE
               mesh_loc%neigh(n, m) = m_glob_to_loc(mop)
            END IF
         END DO
      END DO
      !==End re-order neigh

      !==Re-order i_d
      ALLOCATE(mesh_loc%i_d(mesh_loc%me))
      mesh_loc%i_d = mesh%i_d(m_loc_to_glob(1:mesh_loc%me))
      !==End re-order i_d

      !==Re-order neighs
      ALLOCATE(mesh_loc%neighs(mesh_loc%mes))
      mesh_loc%neighs = m_glob_to_loc(mesh%neighs(mes_loc(1):mes_loc(2)))
      !==End re-order neighs


      !==Re-order sides
      ALLOCATE(mesh_loc%sides(mesh_loc%mes))
      mesh_loc%sides = mesh%sides(mes_loc(1):mes_loc(2))
      !==End re-order sides

      !==Re-order jjs
      ALLOCATE(mesh_loc%jjs(nws, mesh_loc%mes))
      DO ns = 1, nws
         mesh_loc%jjs(ns, :) = glob_to_loc(mesh%jjs(ns, mes_loc(1):mes_loc(2)))
      END DO
      !==End re-order jjs

      !==Re-order jce
      ALLOCATE(mesh_loc%jce(SIZE(mesh%jce, 1), mesh_loc%me))
      mesh_loc%jce = mesh%jce(:, me_loc(1):me_loc(2))
      !==End re-order jce

      mesh_loc%medge = 0
      mesh_loc%medges = 0
      virgins = .TRUE.
      DO m = 1, mesh_loc%me
         DO n = 1, SIZE(mesh%jce, 1)
            e_glob = mesh%jce(n, m_loc_to_glob(m))
            IF (virgins(e_glob)) THEN
               IF (mesh%neigh(n, m_loc_to_glob(m)) >= me_loc(1).or. mesh%neigh(n, m_loc_to_glob(m)) == 0) THEN
                  mesh_loc%medge = mesh_loc%medge + 1
                  virgins(e_glob) = .FALSE.
               ELSE
                  mesh_loc%medges = mesh_loc%medges + 1
               END IF
            END IF
         END DO
      END DO

      ALLOCATE(mesh_loc%jees(mesh_loc%medges))
      ALLOCATE(mesh_loc%jecs(mesh_loc%medges))
      !ALLOCATE(mesh_loc%jevs(SIZE(mesh%jev, 1), mesh_loc%medges))
      virgins = .TRUE.
      medge = 0
      medges = 0
      DO m = 1, mesh_loc%me
         DO n = 1, SIZE(mesh%jce, 1)
            e_glob = mesh%jce(n, m_loc_to_glob(m))
            IF (virgins(e_glob)) THEN
               IF (mesh%neigh(n, m_loc_to_glob(m)) >= me_loc(1) .or. mesh%neigh(n, m_loc_to_glob(m)) == 0) THEN
                  virgins(mesh%jce(n, m_loc_to_glob(m))) = .FALSE.
                  medge = medge + 1
               ELSE
                  medges = medges + 1
                  mesh_loc%jecs(medges) = m
                  mesh_loc%jees(medges) = e_glob
               END IF
            END IF
         END DO
      END DO

      CALL MPI_ALLGATHER(mesh_loc%medge, 1, MPI_INTEGER, mesh_loc%domedge, 1, &
           MPI_INTEGER, PETSC_COMM_WORLD, ierr)
      mesh_loc%disedge(1) = 1
      DO n = 1, nb_proc
         mesh_loc%disedge(n + 1) = mesh_loc%disedge(n) + mesh_loc%domedge(n)
      END DO

      !==Re-order jev
      !ALLOCATE(mesh_loc%jev(SIZE(mesh%jev, 1), mesh_loc%medge))
      !mesh_loc%jev = mesh%jev(:, mesh_loc%ltg_edge(1:mesh_loc%medge))
      !==End re-order jev
      DO proc = 1, nb_proc
         IF (mesh_loc%loc_to_glob(1) <= mesh_loc%disp(proc))    EXIT
      END DO
      nb_extra = 0
      DO m = 1, mesh%me
         jglob = mesh%jj(:, m)
         eglob = mesh%jce(:, m)
         DO n = 1, 3
            IF (jglob(n) < mesh_loc%loc_to_glob(1)) jglob(n) = -1
            IF (jglob(n) > mesh_loc%loc_to_glob(1) + mesh_loc%dom_np - 1) jglob(n) = -1
            IF (eglob(n) < mesh_loc%disedge(proc)) eglob(n) = -1
            IF (eglob(n) >= mesh_loc%disedge(proc + 1)) eglob(n) = -1
         END DO
         IF (MAXVAL(jglob) < 0 .and. MAXVAL(eglob) < 0) cycle
         IF (m<me_loc(1)) THEN
            CALL ERROR_PETSC('BUG  create_local_mesh_with_extra_layer')
         ELSE IF (me_loc(2)<m) THEN
            nb_extra = nb_extra + 1
         END IF
      END DO

      mesh_loc%mextra = nb_extra
      ALLOCATE(mesh_loc%extra_jj(nw, nb_extra), mesh_loc%extra_jce(SIZE(mesh%jce, 1), nb_extra), mesh_loc%extra_jcc(nb_extra))
      nb_extra = 0
      DO m = 1, mesh%me
         jglob = mesh%jj(:, m)
         eglob = mesh%jce(:, m)
         DO n = 1, 3
            IF (jglob(n) < mesh_loc%loc_to_glob(1)) jglob(n) = -1
            IF (jglob(n) > mesh_loc%loc_to_glob(1) + mesh_loc%dom_np - 1) jglob(n) = -1
            IF (eglob(n) < mesh_loc%disedge(proc)) eglob(n) = -1
            IF (eglob(n) >= mesh_loc%disedge(proc + 1)) eglob(n) = -1
         END DO
         IF (MAXVAL(jglob) < 0  .and. MAXVAL(eglob) < 0) cycle
         IF (me_loc(2)<m) THEN
            nb_extra = nb_extra + 1
            mesh_loc%extra_jj(:, nb_extra) = mesh%jj(:, m)
            mesh_loc%extra_jce(:, nb_extra) = mesh%jce(:, m)
            mesh_loc%extra_jcc(nb_extra) = m
         END IF
      END DO

      mesh_loc%edge_stab = .FALSE.
      mesh_loc%mi = 0


      !===Find the isolated points on the border
      nb_extra = 0
      virgin = .TRUE.
      DO m = 1, mesh%mes
         DO i = 1, 2
            j = mesh%jjs(i, m)
            IF (np_loc(1)<=j .AND. j<=np_loc(2)) THEN
               test = .TRUE.
               DO ms = mes_loc(1), mes_loc(2)
                  IF (MINVAL(ABS(j - mesh%jjs(:, ms))) == 0) test = .FALSE.
               END DO
               IF (test .AND. virgin(j)) THEN
                  nb_extra = nb_extra + 1
                  virgins(j) = .FALSE.
               END IF
            END IF
         END DO
      END DO

      mesh_loc%nis = nb_extra
      ALLOCATE(mesh_loc%isolated_jjs(mesh_loc%nis), mesh_loc%isolated_interfaces(mesh_loc%nis, 2))
      mesh_loc%isolated_interfaces = -1

      nb_extra = 0
      virgin = .TRUE.
      DO m = 1, mesh%mes
         DO i = 1, 2
            j = mesh%jjs(i, m)
            IF (np_loc(1)<=j .AND. j<=np_loc(2)) THEN
               test = .TRUE.
               DO ms = mes_loc(1), mes_loc(2)
                  IF (MINVAL(ABS(j - mesh%jjs(:, ms))) == 0) test = .FALSE.
               END DO
               IF (test .AND. virgin(j)) THEN
                  nb_extra = nb_extra + 1
                  virgins(j) = .FALSE.
                  mesh_loc%isolated_jjs(nb_extra) = j
                  DO ms = 1, mesh%mes
                     n = 0
                     IF (MINVAL(ABS(j - mesh%jjs(:, ms))) == 0) THEN
                        n = n + 1
                        mesh_loc%isolated_interfaces(nb_extra, n) = mesh%sides(m)
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END DO
      !===END Find the isolated points on the border

   END SUBROUTINE create_local_mesh_with_extra_layer

END MODULE two_dim_metis_distribution
