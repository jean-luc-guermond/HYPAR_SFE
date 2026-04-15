MODULE two_dim_metis_distribution
#include "petsc/finclude/petsc.h"
   USE petsc
   USE mesh_tools
   PUBLIC :: part_mesh, extract_mesh
   PRIVATE
   REAL(KIND = 8) :: epsilon = 1.d-10
CONTAINS
  SUBROUTINE part_mesh(nb_proc, mesh, list_of_interfaces, part)!, my_periodics)
    USE def_type_mesh
    USE my_util
    USE sub_plot
    ! USE periodic_data_module
    USE mesh_parameters
    IMPLICIT NONE
!!$ Dummy for metis...
    INTEGER, PARAMETER :: METIS_NOPTIONS = 40
!!$ Dummy for metis...
    TYPE(mesh_type) :: mesh
    INTEGER, DIMENSION(mesh%me) :: part
    INTEGER, DIMENSION(:) :: list_of_interfaces
    ! TYPE(periodic_type), DIMENSION(:), TARGET, OPTIONAL :: my_periodics
    ! TYPE(periodic_type), POINTER :: my_periodic
    LOGICAL, DIMENSION(mesh%mes) :: virgins
    INTEGER, DIMENSION(3, mesh%me) :: neigh_new
    INTEGER, DIMENSION(5) :: opts
    INTEGER, DIMENSION(SIZE(mesh%jjs, 1)) :: i_loc
    INTEGER, DIMENSION(:), ALLOCATABLE :: xind_dom, xadj_dom
    INTEGER, DIMENSION(:), ALLOCATABLE :: vwgt, adjwgt
    INTEGER, DIMENSION(1) :: jm_loc
    INTEGER, DIMENSION(mesh%np, 3) :: per_pts
    INTEGER, DIMENSION(mesh%np) :: indicator
    INTEGER, DIMENSION(3) :: j_loc
    INTEGER :: nb_neigh, edge, m, ms, n, nb, numflag, p, wgtflag, j, &
         ns, nws, msop, nsop, proc, iop, mop, s2, k, me
    REAL(KIND = 8) :: err
    LOGICAL :: test
    !===(JLG) Feb 20, 2019. Petsc developpers decide to use REAL(KIND=4) to interface with metis
    !REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: tpwgts
    !REAL(KIND=8), DIMENSION(1)               :: ubvec
    REAL(KIND = 4), DIMENSION(:), ALLOCATABLE :: tpwgts
    REAL(KIND = 4), DIMENSION(1) :: ubvec
    REAL(KIND = 4) :: one_K4 = 1.0
    !===(JLG)Feb 20, 2019.
    INTEGER, DIMENSION(METIS_NOPTIONS) :: metis_opt
    PetscMPIInt    :: nb_proc
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED
    PetscErrorCode :: ierr
    PetscMPIInt    :: rank
    CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED

    me = mesh%me
    IF (me == 0) THEN
       RETURN
    END IF

    IF (nb_proc==1) THEN
       part = 1
       RETURN
    END IF

    neigh_new = mesh%neigh

    !===Create neigh_new for interfaces
    nws = SIZE(mesh%jjs, 1)
    IF (SIZE(list_of_interfaces)/=0) THEN
       virgins = .TRUE.
       DO ms = 1, mesh%mes
          IF (.NOT.virgins(ms)) CYCLE
          IF (MINVAL(ABS(mesh%sides(ms) - list_of_interfaces))/=0) CYCLE !==ms not on a cut
          i_loc = mesh%jjs(:, ms)
          DO msop = 1, mesh%mes
             IF (msop==ms .OR. .NOT.virgins(msop)) CYCLE
             IF (MINVAL(ABS(mesh%sides(msop) - list_of_interfaces))/=0) CYCLE !==msop not on a cut
             DO ns = 1, nws
                test = .FALSE.
                DO nsop = 1, nws
                   iop = mesh%jjs(nsop, msop)
                   IF (MAXVAL(ABS(mesh%rr(:, i_loc(ns)) - mesh%rr(:, iop))).LT.epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   EXIT !==This msop does not coincide with ms
                END IF
             END DO
             IF (test) EXIT
          END DO
          IF (.NOT.test) THEN
             CALL error_Petsc('BUG in part_mesh_M_T_H_phi, .NOT.test ')
          END IF
          DO n = 1, 3
             IF (neigh_new(n, mesh%neighs(msop))==0) THEN
                neigh_new(n, mesh%neighs(msop)) = mesh%neighs(ms)
             END IF
             IF (neigh_new(n, mesh%neighs(ms))==0) THEN
                neigh_new(n, mesh%neighs(ms)) = mesh%neighs(msop)
             END IF
          END DO
          virgins(ms) = .FALSE.
          virgins(msop) = .FALSE.
       END DO
    END IF
    !===End Create neigh_new for interfaces

    !===Create neigh_new for periodic faces
    IF (mesh_data_info%nb_bords/=0) THEN
       ! IF (PRESENT(my_periodics)) THEN
       ! DO k = 1, SIZE(my_periodics)
       !  my_periodic => my_periodics(k)
       !  IF (mesh_%nb_bords/=0) THEN
       DO ms = 1, mesh%mes
          m = mesh%neighs(ms)
          IF (MINVAL(ABS(mesh%sides(ms) - mesh_data_info%list_periodic(1, :))) == 0) THEN
             jm_loc = MINLOC(ABS(mesh%sides(ms) - mesh_data_info%list_periodic(1, :)))
             s2 =  mesh_data_info%list_periodic(2, jm_loc(1))
             test = .FALSE.
             DO msop = 1, mesh%mes
                IF (mesh%sides(msop) /= s2) CYCLE

                err = 0.d0
                DO ns = 1, SIZE( mesh_data_info%vect_e, 1)
                   err = err + ABS(SUM(mesh%rr(ns, mesh%jjs(:, ms)) - mesh%rr(ns, mesh%jjs(:, msop)) &
                        +mesh_data_info%vect_e(ns, jm_loc(1))))
                END DO

                IF (err .LE. epsilon) THEN
                   test = .TRUE.
                   EXIT
                END IF
             END DO
             IF (.NOT.test) THEN
                CALL error_Petsc('BUG in part_mesh_M_T_H_phi, mop not found')
             END IF
             mop = mesh%neighs(msop)
             DO n = 1, 3
                IF (neigh_new(n, m) == 0) THEN
                   neigh_new(n, m) = mop
                END IF
                IF (neigh_new(n, mop) == 0) THEN
                   neigh_new(n, mop) = m
                END IF
             END DO
          END IF
       END DO
       ! END IF
       !END DO
    END IF
    !===End Create neigh_new for periodic faces


    !===Create the connectivity arrays Xind and Xadj based on neigh (for Metis)
    nb_neigh = SIZE(mesh%neigh, 1)
    ALLOCATE(xind_dom(me + 1))

    xind_dom(1) = 1
    DO k = 1, me
       nb = 0
       DO n = 1, nb_neigh
          mop = neigh_new(n, k)
          IF (mop==0) CYCLE
          nb = nb + 1
       END DO
       xind_dom(k + 1) = xind_dom(k) + nb
    END DO

    ALLOCATE(xadj_dom(xind_dom(me + 1) - 1))
    p = 0
    DO k = 1, me
       DO n = 1, nb_neigh
          mop = neigh_new(n, k)
          IF (mop==0) CYCLE
          p = p + 1
          xadj_dom(p) =  neigh_new(n,k) !=== (= k) Bug corrected March 20, 2026
       END DO
    END DO
    !TESTTT
    !DO k = 1, me
    !   WRITE(*,*) 'xind_dom, xadj_dom', xind_dom(k), xadj_dom(k)
    !   !WRITE(*,*) mesh%neigh(:,k)
    !END DO
    !TESTT
    IF (p/=xind_dom(me + 1) - 1) THEN
       CALL error_Petsc('BUG in  part_mesh, p/=xind_dom(me+1)-1')
    END IF
    !===End Create the connectivity arrays Xind and Xadj based on neigh (for Metis)

    !===Create partitions
    opts = 0
    numflag = 1
    wgtflag = 2
    ALLOCATE(tpwgts(nb_proc))
    tpwgts = one_K4 / nb_proc
    CALL METIS_SetDefaultOptions(metis_opt)
    metis_opt(18) = 1
    ubvec = 1.001

    ALLOCATE(vwgt(me), adjwgt(SIZE(xadj_dom)))
    vwgt = 1
    adjwgt = 1
    CALL METIS_PartGraphRecursive(me, 1, xind_dom, xadj_dom, vwgt, vwgt, adjwgt, nb_proc, tpwgts, &
         ubvec, metis_opt, edge, part)
    !===End Create partitions

    !TESTTTTT
    !IF (rank==0) THEN
    !   CALL plot_const_p1_label(mesh%jj, mesh%rr, 1.d0 * part, 'dd.plt')
    !END IF
    !STOP
    !TEST

    !===Create parts and modify part
    !===Search on the boundary whether ms is on a cut.
    IF (SIZE(mesh%jj, 1)/=3) THEN
       WRITE(*, *) 'SIZE(mesh%jj,1)', SIZE(mesh%jj, 1)
       CALL error_Petsc('BUG in part_mesh_M_T_H_phi, SIZE(mesh%jj,1)/=3')
    END IF
    indicator = -1
    nws = SIZE(mesh%jjs, 1)
    IF (SIZE(list_of_interfaces)/=0) THEN
       virgins = .TRUE.
       DO ms = 1, mesh%mes
          IF (.NOT.virgins(ms)) CYCLE
          IF (MINVAL(ABS(mesh%sides(ms) - list_of_interfaces))/=0) CYCLE !==ms not on a cut
          i_loc = mesh%jjs(:, ms)
          DO msop = 1, mesh%mes
             IF (msop==ms .OR. .NOT.virgins(msop)) CYCLE
             IF (MINVAL(ABS(mesh%sides(msop) - list_of_interfaces))/=0) CYCLE !==msop not on a cut
             DO ns = 1, nws
                test = .FALSE.
                DO nsop = 1, nws
                   iop = mesh%jjs(nsop, msop)
                   IF (MAXVAL(ABS(mesh%rr(:, i_loc(ns)) - mesh%rr(:, iop))).LT.epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   EXIT !==This msop does not coincide with ms
                END IF
             END DO
             IF (test) EXIT
          END DO
          IF (.NOT.test) THEN
             CALL error_Petsc('BUG in part_mesh_M_T_H_phi, .NOT.test ')
          END IF
          IF (part(mesh%neighs(ms)) == part(mesh%neighs(msop))) CYCLE !==ms is an internal cut
          proc = MIN(part(mesh%neighs(ms)), part(mesh%neighs(msop)))
          part(mesh%neighs(ms)) = proc !make sure interface are internal
          part(mesh%neighs(msop)) = proc !make sure interface are internal
          virgins(ms) = .FALSE.
          virgins(msop) = .FALSE.
          indicator(mesh%jjs(:, ms)) = proc
          indicator(mesh%jjs(:, msop)) = proc
       END DO
    END IF
    !===Fix the partition so that all the cells having one vertex on an
    !===interface belong to the same processor as those sharing this vertices and
    !===having two vertices on the interface (JLG + DCQ July 22 2015)
    DO m = 1, mesh%me
       j_loc = mesh%jj(:, m)
       n = MAXVAL(indicator(j_loc))
       IF (n == -1) CYCLE
       IF (indicator(j_loc(1)) * indicator(j_loc(2)) * indicator(j_loc(3))<0) CYCLE
       part(m) = n
    END DO
    !===End create parts and modify part

    !===Move the two elements with one periodic face on same processor
    !IF (PRESENT(my_periodics)) THEN
    !DO k = 1, SIZE(my_periodics)
    !my_periodic => my_periodics(k)
    IF (mesh_data_info%nb_bords/=0) THEN
       DO j = 1, mesh%np
          per_pts(j, 1) = j
       END DO
       per_pts(:, 2:3) = 0
       DO ms = 1, mesh%mes
          m = mesh%neighs(ms)
          IF ((MINVAL(ABS(mesh%sides(ms) - mesh_data_info%list_periodic(1, :))) /=0) .AND. &
               (MINVAL(ABS(mesh%sides(ms) - mesh_data_info%list_periodic(2, :))) /=0)) CYCLE
          DO ns = 1, SIZE(mesh%jjs, 1)
             j = mesh%jjs(ns, ms)
             per_pts(j, 2) = m
             DO msop = 1, mesh%mes
                IF (MINVAL(ABS(mesh%sides(msop) - mesh_data_info%list_periodic(:, :))) /=0) CYCLE
                IF (msop == ms) CYCLE
                DO nsop = 1, SIZE(mesh%jjs, 1)
                   IF (mesh%jjs(nsop, msop)==j) THEN
                      per_pts(j, 3) = mesh%neighs(msop)
                   END IF
                END DO
             END DO
          END DO
       END DO
       CALL reassign_per_pts(mesh, part, per_pts)
       DO ms = 1, mesh%mes
          m = mesh%neighs(ms)
          IF (MINVAL(ABS(mesh%sides(ms) - mesh_data_info%list_periodic(1, :))) /= 0) CYCLE
          jm_loc = MINLOC(ABS(mesh%sides(ms) - mesh_data_info%list_periodic(1, :)))
          s2 = mesh_data_info%list_periodic(2, jm_loc(1))
          test = .FALSE.
          DO msop = 1, mesh%mes
             IF (mesh%sides(msop) /= s2) CYCLE
             err = 0.d0
             DO ns = 1, SIZE(mesh_data_info%vect_e, 1)
                err = err + ABS(SUM(mesh%rr(ns, mesh%jjs(:, ms)) - mesh%rr(ns, mesh%jjs(:, msop)) &
                     + mesh_data_info%vect_e(ns, jm_loc(1))))
             END DO
             IF (err .LE. epsilon) THEN
                test = .TRUE.
                EXIT
             END IF
          END DO
          IF (.NOT.test) THEN
             CALL error_Petsc('BUG in part_mesh_M_T_H_phi, mop not found')
          END IF
          IF (part(mesh%neighs(ms)) /= part(mesh%neighs(msop))) THEN !==ms is an internal cut
             proc = MIN(part(mesh%neighs(ms)), part(mesh%neighs(msop)))
             part(mesh%neighs(ms)) = proc !make sure interface are internal
             part(mesh%neighs(msop)) = proc !make sure interface are internal
          END IF
       END DO
    END IF
    !END DO
    !END IF
    !===End Move the two elements with one periodic face on same processor

!!$ WARNING, FL 1/2/13 : TO BE ADDED IF NEEDED
    !================================================
    IF (rank==0) THEN
       CALL plot_const_p1_label(mesh%jj, mesh%rr, 1.d0 * part, 'dd.plt')
    END IF
    !================================================
!!$ WARNING, FL 1/2/13 : TO BE ADDED IF NEEDED

    DEALLOCATE(vwgt, adjwgt)
    DEALLOCATE(xadj_dom)
    DEALLOCATE(xind_dom)

    DEALLOCATE(tpwgts)

  END SUBROUTINE part_mesh

   SUBROUTINE extract_mesh(communicator, nb_proc, mesh_glob, part, list_dom, mesh_loc)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_glob, mesh, mesh_loc
      INTEGER, DIMENSION(:) :: part, list_dom
      INTEGER, DIMENSION(mesh_glob%me) :: bat
      INTEGER, DIMENSION(mesh_glob%np) :: i_old_to_new
      INTEGER, DIMENSION(mesh_glob%medge) :: old_edge_to_new
      INTEGER, DIMENSION(mesh_glob%mes) :: parts
      INTEGER, DIMENSION(nb_proc) :: nblmt_per_proc, start, displ
      INTEGER, DIMENSION(2) :: np_loc, me_loc, mes_loc
      INTEGER, DIMENSION(:), ALLOCATABLE :: list_m, tab, tabs
      INTEGER :: nb_proc, ms, i, index, m, mop, n, j
      PetscErrorCode :: ierr
      PetscMPIInt    :: rank
      MPI_Comm       :: communicator
      CALL MPI_Comm_rank(communicator, rank, ierr)

      ALLOCATE(mesh_loc%disp(nb_proc + 1), mesh_loc%domnp(nb_proc))
      ALLOCATE(mesh_loc%discell(nb_proc + 1), mesh_loc%domcell(nb_proc))
      ALLOCATE(mesh_loc%disedge(nb_proc + 1), mesh_loc%domedge(nb_proc))

      ! Create parts
      parts = part(mesh_glob%neighs)
      ! End create parts

      ! Create list_m
      i = 0
      DO m = 1, mesh_glob%me
         IF (MINVAL(ABS(list_dom - mesh_glob%i_d(m)))/=0) CYCLE
         i = i + 1
      END DO
      mesh%me = i
      ALLOCATE (list_m(mesh%me))
      i = 0
      DO m = 1, mesh_glob%me
         IF (MINVAL(ABS(list_dom - mesh_glob%i_d(m)))/=0) CYCLE
         i = i + 1
         list_m(i) = m
      END DO
      !End create list_m

      ! Count elements on processors
      nblmt_per_proc = 0
      DO i = 1, mesh%me
         m = list_m(i)
         nblmt_per_proc(part(m)) = nblmt_per_proc(part(m)) + 1
      END DO
      start(1) = 0
      DO n = 2, nb_proc
         start(n) = start(n - 1) + nblmt_per_proc(n - 1)
      END DO
      me_loc(1) = start(rank + 1) + 1
      me_loc(2) = start(rank + 1) + nblmt_per_proc(rank + 1)
      displ = start
      ! End count elements on processors

      ! Re-order elements
      ALLOCATE(tab(mesh%me))
      bat = 0
      DO i = 1, mesh%me
         m = list_m(i)
         start(part(m)) = start(part(m)) + 1
         tab(start(part(m))) = m
         bat(m) = start(part(m))
      END DO
      ! Re-order elements

      ! Create mesh%jj
      mesh%gauss%n_w = SIZE(mesh_glob%jj, 1)
      ALLOCATE(mesh%jj(SIZE(mesh_glob%jj, 1), mesh%me))
      i_old_to_new = 0
      index = 0
      DO m = 1, mesh%me
         DO n = 1, SIZE(mesh_glob%jj, 1)
            i = mesh_glob%jj(n, tab(m))
            IF (i_old_to_new(i)/=0) THEN
               mesh%jj(n, m) = i_old_to_new(i)
            ELSE
               index = index + 1
               i_old_to_new(i) = index
               mesh%jj(n, m) = i_old_to_new(i)
            END IF
         END DO
      END DO
      mesh%np = index
      ! End Create mesh%jj

      ! Create mesh%rr
      ALLOCATE(mesh%rr(2, mesh%np))
      DO i = 1, mesh_glob%np
         IF (i_old_to_new(i)==0) CYCLE
         mesh%rr(:, i_old_to_new(i)) = mesh_glob%rr(:, i)
      END DO
      !End Create mesh%rr

      ! Re-order edge
      ALLOCATE(mesh%jce(SIZE(mesh_glob%jce, 1), mesh%me))
      old_edge_to_new = 0
      index = 0
      DO m = 1, mesh%me
         DO n = 1, SIZE(mesh%jce, 1)
            j = mesh_glob%jce(n, tab(m))
            IF (old_edge_to_new(j)/=0) THEN
               mesh%jce(n, m) = old_edge_to_new(j)
            ELSE
               index = index + 1
               old_edge_to_new(j) = index
               mesh%jce(n, m) = old_edge_to_new(j)
            END IF
         END DO
      END DO
      mesh%medge = index
      ! End re-order edge


      ! Create mesh%neigh
      ALLOCATE(mesh%neigh(3, mesh%me))
      DO m = 1, mesh%me
         DO n = 1, 3
            mop = mesh_glob%neigh(n, tab(m))
            IF (mop==0) THEN
               mesh%neigh(n, m) = 0
            ELSE
               mesh%neigh(n, m) = bat(mop)
            END IF
         END DO
      END DO
      ! End  Create mesh%neigh

      ! Create mesh%i_d
      ALLOCATE(mesh%i_d(mesh%me))
      mesh%i_d = mesh_glob%i_d(tab)
      ! End mesh%i_d

      ! Create np_loc
      IF (displ(rank + 1)/=0) THEN
         np_loc(1) = MAXVAL(mesh%jj(:, 1:displ(rank + 1))) + 1
      ELSE
         np_loc(1) = 1
      END IF
      np_loc(2) = np_loc(1) - 1
      IF (me_loc(1).LE.me_loc(2)) THEN
         np_loc(2) = MAXVAL(mesh%jj(:, me_loc(1):me_loc(2)))
      END IF
      IF (np_loc(2) .LT. np_loc(1) - 1) THEN
         np_loc(2) = np_loc(1) - 1
      END IF
      ! End create np_loc

      ! Create mes_loc
      nblmt_per_proc = 0
      DO ms = 1, mesh_glob%mes
         IF (MINVAL(ABS(list_dom - mesh_glob%i_d(mesh_glob%neighs(ms))))/=0) CYCLE
         n = parts(ms)
         nblmt_per_proc(n) = nblmt_per_proc(n) + 1
      END DO
      start(1) = 0
      DO n = 2, nb_proc
         start(n) = start(n - 1) + nblmt_per_proc(n - 1)
      END DO
      mes_loc(1) = start(rank + 1) + 1
      mes_loc(2) = start(rank + 1) + nblmt_per_proc(rank + 1)
      mesh%mes = SUM(nblmt_per_proc)
      ! End create mes_loc

      ! Create tabs and sbat
      ALLOCATE(tabs(mesh%mes))
      DO ms = 1, mesh_glob%mes
         IF (MINVAL(ABS(list_dom - mesh_glob%i_d(mesh_glob%neighs(ms))))/=0) CYCLE
         start(parts(ms)) = start(parts(ms)) + 1
         tabs(start(parts(ms))) = ms
      END DO
      ! End create tabs and sbat

      ! Create neighs
      ALLOCATE(mesh%neighs(mesh%mes))
      mesh%neighs = bat(mesh_glob%neighs(tabs))
      ! End create neighs

      ! Re-order sides
      ALLOCATE(mesh%sides(mesh%mes))
      mesh%sides = mesh_glob%sides(tabs)
      ! End re-order sides

      ! Re-order jjs
      mesh%gauss%n_ws = SIZE(mesh_glob%jjs, 1)
      ALLOCATE(mesh%jjs(SIZE(mesh_glob%jjs, 1), mesh%mes))

      DO n = 1, SIZE(mesh%jjs, 1)
         mesh%jjs(n, :) = i_old_to_new(mesh_glob%jjs(n, tabs))
      END DO
      ! End re-order jjs

      mesh%mes_int = 0
      ! Create mes_int_loc
      DO ms = 1, mesh_glob%mes_int
         IF (MAXVAL(bat(mesh_glob%neighs_int(:, ms))) > 0) THEN
            mesh%mes_int = mesh%mes_int + 1
         END IF
      END DO

      ! Create neighs_int sides_int jjs_int
      ALLOCATE(mesh%neighs_int(2, mesh%mes_int))
      ALLOCATE(mesh%sides_int(mesh%mes_int))
      ALLOCATE(mesh%jjs_int(SIZE(mesh_glob%jjs, 1), mesh%mes_int))

      m = 0
      DO ms = 1, mesh_glob%mes_int
         IF (MAXVAL(bat(mesh_glob%neighs_int(:, ms))) > 0) THEN
            m = m + 1
            mesh%neighs_int(1, m) = bat(mesh_glob%neighs_int(1, ms))
            mesh%neighs_int(2, m) = bat(mesh_glob%neighs_int(2, ms))
            IF (mesh%neighs_int(2, m) > mesh%neighs_int(1, m)) THEN
               mesh%neighs_int(:, m) = (/mesh%neighs_int(2, m), mesh%neighs_int(1, m) /)
            END IF
            mesh%sides_int(m) = mesh_glob%sides_int(ms)
            DO n = 1, SIZE(mesh%jjs, 1)
               mesh%jjs_int(n, m) = i_old_to_new(mesh_glob%jjs_int(n, ms))
            END DO
         END IF
      END DO
      ! End create neighs_int


      !==We create the local mesh now
      mesh%edge_stab = .FALSE.
      ALLOCATE(mesh%jees(0), mesh%jecs(0))

      mesh%mextra = 0
      ALLOCATE(mesh%jj_extra(3, 0))
      ALLOCATE(mesh%jce_extra(3, 0))
      ALLOCATE(mesh%jcc_extra(0))

      mesh%mes_extra = 0
      ALLOCATE(mesh%neighs_extra(0))
      ALLOCATE(mesh%sides_extra(0))
      ALLOCATE(mesh%jjs_extra(2, 0))
      ALLOCATE(mesh%rrs_extra(2, 3, 0))
      ALLOCATE(mesh%loc_to_glob(mesh%np))
      DO n = 1, mesh%np
         mesh%loc_to_glob(n) = n
      END DO
      mesh%nis = 0
      ALLOCATE(mesh%isolated_jjs(0), mesh%isolated_interfaces(0, 2))

      ALLOCATE(mesh%disp(nb_proc + 1), mesh%domnp(nb_proc))
      ALLOCATE(mesh%discell(nb_proc + 1), mesh%domcell(nb_proc))
      ALLOCATE(mesh%disedge(nb_proc + 1), mesh%domedge(nb_proc))
      !mesh%disp = (/ 1, mesh%np + 1 /)
      !mesh%domnp = (/ mesh%np /)
      !mesh%discell = (/ 1, mesh%me + 1 /)
      !mesh%domcell = (/ mesh%me /)
      !mesh%disedge = (/ 1, mesh%medge + 1 /)
      !mesh%domedge = (/ mesh%medge /)
      mesh%disp(1) = 1
      mesh%disp(2) = mesh%np + 1
      mesh%domnp(1) = mesh%np
      mesh%discell(1) = 1
      mesh%discell(2) = mesh%me + 1
      mesh%domcell(1) = mesh%me
      mesh%disedge(1) = 1
      mesh%disedge(2) = mesh%medge + 1
      mesh%domedge(1) = mesh%medge
      CALL create_local_mesh_with_extra_layer(communicator, mesh, mesh_loc, me_loc, mes_loc, np_loc)
      CALL free_mesh(mesh)
      DEALLOCATE(list_m, tab, tabs)

   END SUBROUTINE extract_mesh


   SUBROUTINE create_local_mesh_with_extra_layer(communicator, mesh, mesh_loc, me_loc, mes_loc, np_loc)
      USE def_type_mesh
      USE my_util
      USE sub_plot
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh, mesh_loc
      INTEGER, DIMENSION(2), INTENT(IN) :: me_loc, mes_loc, np_loc
      INTEGER, DIMENSION(2) :: is1, is2
      INTEGER, DIMENSION(mesh%me) :: m_glob_to_loc, m_loc_to_glob
      INTEGER, DIMENSION(mesh%np) :: glob_to_loc, loc_to_glob
      LOGICAL, DIMENSION(mesh%np) :: virgin
      LOGICAL, DIMENSION(mesh%medge) :: virgins
      LOGICAL, DIMENSION(mesh%me) :: virginss
      LOGICAL, DIMENSION(mesh%me) :: not_my_cells
      INTEGER, DIMENSION(SIZE(mesh%jj, 1)) :: jglob, eglob
      LOGICAL :: test
      INTEGER :: dim, nws, nw, m, ms, mop, ns, msup, minf, dof, proc, m2, &
           dom_me, nwc, dom_mes, dom_np, n, i, ierr, dom_np_glob, nb_extra, nb_proc, e_glob, medge, medges, j
      MPI_Comm :: communicator

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
         mesh_loc%mes_int = mesh%mes_int
         mesh_loc%dom_me = mesh%me
         mesh_loc%dom_np = mesh%np
         mesh_loc%dom_mes = mesh%mes
         mesh_loc%mextra = 0
         mesh_loc%mes_extra = 0
         mesh_loc%medge = mesh%medge
         mesh_loc%medges = 0
         mesh_loc%nis = 0
         mesh_loc%nps = 0

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

         ALLOCATE(mesh_loc%neighs_int(2, mesh_loc%mes_int))
         mesh_loc%neighs_int = mesh%neighs_int
         ALLOCATE(mesh_loc%sides_int(mesh_loc%mes_int))
         mesh_loc%sides_int = mesh%sides_int
         ALLOCATE(mesh_loc%jjs_int(nws, mesh_loc%mes_int))
         mesh_loc%jjs_int = mesh%jjs_int

         ALLOCATE(mesh_loc%neighs_extra(mesh_loc%mes_extra))
         ALLOCATE(mesh_loc%sides_extra(mesh_loc%mes_extra))
         ALLOCATE(mesh_loc%jjs_extra(nws, mesh_loc%mes_extra))
         ALLOCATE(mesh_loc%rrs_extra(dim, nw, mesh_loc%mes_extra))

         ALLOCATE(mesh_loc%jj_extra(nw, mesh_loc%mextra))
         ALLOCATE(mesh_loc%jce_extra(nw, mesh_loc%mextra))
         ALLOCATE(mesh_loc%jcc_extra(mesh_loc%mextra))

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

      !==Create the new mesh
      dom_me = me_loc(2) - me_loc(1) + 1
      dom_mes = mes_loc(2) - mes_loc(1) + 1
      dom_np = np_loc(2) - np_loc(1) + 1
      mesh_loc%me = dom_me
      mesh_loc%mes = dom_mes
      mesh_loc%dom_me = dom_me
      mesh_loc%dom_np = dom_np
      mesh_loc%dom_mes = dom_mes
      CALL MPI_ALLREDUCE(dom_np, dom_np_glob, 1, MPI_INTEGER, &
           MPI_MIN, communicator, ierr)
      IF (dom_np_glob.LE.0) THEN
         CALL error_petsc('Pb in create_local_mesh, not enough cells per processors')
      END IF

      CALL MPI_ALLGATHER(mesh_loc%dom_np, 1, MPI_INTEGER, mesh_loc%domnp, 1, &
           MPI_INTEGER, communicator, ierr)
      mesh_loc%disp(1) = 1
      DO n = 1, nb_proc
         mesh_loc%disp(n + 1) = mesh_loc%disp(n) + mesh_loc%domnp(n)
      END DO

      CALL MPI_ALLGATHER(mesh_loc%me, 1, MPI_INTEGER, mesh_loc%domcell, 1, &
           MPI_INTEGER, communicator, ierr)
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

      IF (SIZE(mesh%jj, 1) == 6) THEN
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
      IF (SIZE(mesh%jj, 1) == 10) THEN
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

      mesh_loc%mes_int = 0
      !===Count number of internal edges
      DO ms = 1, mesh%mes_int
         test = m_glob_to_loc(mesh%neighs_int(1, ms)) > 0
         IF (mesh%neighs_int(2, ms) > 0) THEN
            test = test .OR. m_glob_to_loc(mesh%neighs_int(2, ms)) > 0
         END IF
         IF (test)  THEN
            mesh_loc%mes_int = mesh_loc%mes_int + 1
         END IF
      END DO

      !==Re-order neighs_int
      ALLOCATE(mesh_loc%neighs_int(2, mesh_loc%mes_int))
      ALLOCATE(mesh_loc%sides_int(mesh_loc%mes_int))
      ALLOCATE(mesh_loc%jjs_int(nws, mesh_loc%mes_int))

      ms = 0
      DO m = 1, mesh%mes_int
         test = m_glob_to_loc(mesh%neighs_int(1, m)) > 0
         IF (mesh%neighs_int(2, m) > 0) THEN
            test = test .OR. m_glob_to_loc(mesh%neighs_int(2, m)) > 0
         END IF
         IF (test)  THEN
            ms = ms + 1
            mesh_loc%neighs_int(:, ms) = m_glob_to_loc(mesh%neighs_int(:, m))
            IF (mesh_loc%neighs_int(2, ms) > mesh_loc%neighs_int(1, ms)) THEN
               mesh_loc%neighs_int(:, ms) = (/mesh_loc%neighs_int(2, ms), mesh_loc%neighs_int(1, ms) /)
            END IF

            mesh_loc%sides_int(ms) = mesh%sides_int(m)
            DO ns = 1, nws
               mesh_loc%jjs_int(ns, ms) = glob_to_loc(mesh%jjs_int(ns, m))
            END DO
         END IF
      END DO
      !==End re-order neighs

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
               IF (mesh%neigh(n, m_loc_to_glob(m)) >= me_loc(1).OR. mesh%neigh(n, m_loc_to_glob(m)) == 0) THEN
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
               IF (mesh%neigh(n, m_loc_to_glob(m)) >= me_loc(1) .OR. mesh%neigh(n, m_loc_to_glob(m)) == 0) THEN
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
           MPI_INTEGER, communicator, ierr)
      mesh_loc%disedge(1) = 1
      DO n = 1, nb_proc
         mesh_loc%disedge(n + 1) = mesh_loc%disedge(n) + mesh_loc%domedge(n)
      END DO

      !==Re-order jev
      !ALLOCATE(mesh_loc%jev(SIZE(mesh%jev, 1), mesh_loc%medge))
      !mesh_loc%jev = mesh%jev(:, mesh_loc%ltg_edge(1:mesh_loc%medge))
      !==End re-order jev


      !==Building extra cells
      virginss = .TRUE.
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
         IF (MAXVAL(jglob) < 0 .AND. MAXVAL(eglob) < 0) CYCLE
         IF (m<me_loc(1)) THEN
            CALL ERROR_PETSC('BUG  create_local_mesh_with_extra_layer')
         ELSE IF (me_loc(2)<m) THEN
            nb_extra = nb_extra + 1
            virginss(m) = .FALSE.
            IF (MINVAL(ABS(mesh%neighs - m)) == 0) THEN
               is1 = 0
               is2 = 0
               CALL find_cell_interface(mesh, m, m2, is1, is2)
               IF (m2 > me_loc(2)) THEN
                  IF (virginss(m2)) THEN
                     nb_extra = nb_extra + 1
                     virginss(m2) = .FALSE.
                  END IF
               END IF
               DO i = 1, 2
                  IF (is1(i) < mesh_loc%loc_to_glob(1) .OR. is1(i) > mesh_loc%loc_to_glob(1) + mesh_loc%dom_np - 1) THEN
                     CYCLE
                  END IF
                  DO m2 = 1, mesh%me
                     IF (MINVAL(ABS(mesh%jj(:, m2) - is2(i))) == 0) THEN
                        IF (m2 > me_loc(2)) THEN
                           IF (virginss(m2)) THEN
                              nb_extra = nb_extra + 1
                              virginss(m2) = .FALSE.
                           END IF
                        END IF
                     END IF
                  END DO
               END DO
            END IF
         END IF
      END DO

      mesh_loc%mextra = nb_extra
      ALLOCATE(mesh_loc%jj_extra(nw, nb_extra), mesh_loc%jce_extra(SIZE(mesh%jce, 1), nb_extra), &
           mesh_loc%jcc_extra(nb_extra))
      nb_extra = 0
      virginss = .TRUE.
      DO m = 1, mesh%me
         jglob = mesh%jj(:, m)
         eglob = mesh%jce(:, m)
         DO n = 1, 3
            IF (jglob(n) < mesh_loc%loc_to_glob(1)) jglob(n) = -1
            IF (jglob(n) > mesh_loc%loc_to_glob(1) + mesh_loc%dom_np - 1) jglob(n) = -1
            IF (eglob(n) < mesh_loc%disedge(proc)) eglob(n) = -1
            IF (eglob(n) >= mesh_loc%disedge(proc + 1)) eglob(n) = -1
         END DO
         IF (MAXVAL(jglob) < 0  .AND. MAXVAL(eglob) < 0) CYCLE
         IF (me_loc(2)<m) THEN
            nb_extra = nb_extra + 1
            mesh_loc%jj_extra(:, nb_extra) = mesh%jj(:, m)
            mesh_loc%jce_extra(:, nb_extra) = mesh%jce(:, m)
            mesh_loc%jcc_extra(nb_extra) = m
            virginss(m) = .FALSE.
            IF (MINVAL(ABS(mesh%neighs - m)) == 0) THEN
               is1 = 0
               is2 = 0
               CALL find_cell_interface(mesh, m, m2, is1, is2)
               IF (m2 > me_loc(2)) THEN
                  IF (virginss(m2)) THEN
                     nb_extra = nb_extra + 1
                     mesh_loc%jj_extra(:, nb_extra) = mesh%jj(:, m2)
                     mesh_loc%jce_extra(:, nb_extra) = mesh%jce(:, m2)
                     mesh_loc%jcc_extra(nb_extra) = m2
                     virginss(m2) = .FALSE.
                  END IF
               END IF
               DO i = 1, 2
                  IF (is1(i) < mesh_loc%loc_to_glob(1) .OR. is1(i) > mesh_loc%loc_to_glob(1) + mesh_loc%dom_np - 1) THEN
                     CYCLE
                  END IF
                  DO m2 = 1, mesh%me
                     IF (MINVAL(ABS(mesh%jj(:, m2) - is2(i))) == 0) THEN
                        IF (m2 > me_loc(2)) THEN
                           IF (virginss(m2)) THEN
                              nb_extra = nb_extra + 1
                              mesh_loc%jj_extra(:, nb_extra) = mesh%jj(:, m2)
                              mesh_loc%jce_extra(:, nb_extra) = mesh%jce(:, m2)
                              mesh_loc%jcc_extra(nb_extra) = m2
                              virginss(m2) = .FALSE.
                           END IF
                        END IF
                     END IF
                  END DO
               END DO
            END IF
         END IF
      END DO

      mesh_loc%edge_stab = .FALSE.
      mesh_loc%mi = 0


      !==Building extra edges along interfaces
      nb_extra = 0
      DO ms = 1, mesh%mes
         m = mesh%neighs(ms)
         IF (MINVAL(ABS(mesh_loc%jcc_extra - m)) == 0) THEN
            nb_extra = nb_extra + 1
         END IF
      END DO
      mesh_loc%mes_extra = nb_extra

      ALLOCATE(mesh_loc%jjs_extra(nws, nb_extra), mesh_loc%neighs_extra(nb_extra), &
           mesh_loc%sides_extra(nb_extra), mesh_loc%rrs_extra(2, nw, nb_extra))
      nb_extra = 0
      DO ms = 1, mesh%mes
         m = mesh%neighs(ms)
         IF (MINVAL(ABS(mesh_loc%jcc_extra - m)) == 0) THEN
            nb_extra = nb_extra + 1

            mesh_loc%jjs_extra(:, nb_extra) = mesh%jjs(:, ms)
            mesh_loc%sides_extra(nb_extra) = mesh%sides(ms)
            mesh_loc%neighs_extra(nb_extra) = m
            mesh_loc%rrs_extra(1, :, nb_extra) = mesh%rr(1, mesh%jj(:, m))
            mesh_loc%rrs_extra(2, :, nb_extra) = mesh%rr(2, mesh%jj(:, m))
         END IF
      END DO

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


   SUBROUTINE find_cell_interface(mesh, m1, m2, is1, is2)
      USE def_type_mesh
      TYPE(mesh_type), INTENT(IN) :: mesh
      INTEGER :: m1, m2, ms1, ms2, k, ns
      INTEGER, DIMENSION(2) :: is1, is2
      REAL(KIND = 8) :: eps_ref = 1.d-7, r_norm, epsilon
      LOGICAL :: okay
      INTEGER, DIMENSION(2) :: list

      DO ms1 = 1, mesh%mes
         IF (mesh%neighs(ms1) == m1) EXIT
      END DO
      r_norm = SUM(ABS(mesh%rr(:, mesh%jjs(1, ms1)) - mesh%rr(:, mesh%jjs(2, ms1))))
      epsilon = eps_ref * r_norm
      okay = .FALSE.

      lp2 : DO ms2 = 1, mesh%mes
         DO k = 0, 2
            DO ns = 1, 2
               list(ns) = MODULO(ns - 1 + k, 2) + 1
            END DO
            IF (MAXVAL(ABS(mesh%rr(:, mesh%jjs(list, ms1)) - mesh%rr(:, mesh%jjs(1:2, ms2)))) >= epsilon) CYCLE

            m2 = mesh%neighs(ms2)
            r_norm = SUM(ABS(mesh%rr(:, mesh%jj(1:3, m1)) - mesh%rr(:, mesh%jj(1:3, m2))))
            IF (r_norm <= 1d-9) THEN
               CYCLE
            END IF
            is1 = mesh%jjs(list, ms1)
            is2 = mesh%jjs(1:2, ms2)
            okay = .TRUE.
            EXIT lp2

         END DO
      END DO lp2

      IF (.NOT. okay) m2 = -1

   END SUBROUTINE find_cell_interface

   SUBROUTINE reassign_per_pts(mesh, partition, list_pts)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE

      TYPE(mesh_type), INTENT(IN) :: mesh
      INTEGER, DIMENSION(mesh%me), INTENT(INOUT) :: partition
      INTEGER, DIMENSION(mesh%np, 3), INTENT(IN) :: list_pts

      INTEGER :: i, j_loc, proc_min, index, i_loc, m, mop, n, proc1, proc2
      INTEGER, DIMENSION(50) :: list_elmts
      LOGICAL :: okay

      list_elmts = 0
      DO i = 1, mesh%np
         IF (list_pts(i, 2)==0) CYCLE
         j_loc = list_pts(i, 1)
         list_elmts = 0
         index = 1
         list_elmts(index) = list_pts(i, 2)
         okay = .TRUE.
         DO WHILE (okay)
            m = list_elmts(index)
            okay = .FALSE.
            i_loc = index
            DO n = 1, 3
               mop = mesh%neigh(n, m)
               IF (mop == 0) CYCLE
               IF (MINVAL(ABS(mesh%jj(:, mop) - j_loc)) /=0) CYCLE
               IF (MINVAL(ABS(mop - list_elmts))==0) CYCLE
               okay = .TRUE.
               i_loc = i_loc + 1
               IF (i_loc - index==2) THEN
                  CALL error_Petsc('BUG in reassign_per_pts, how is that possible?')
               END IF
               list_elmts(i_loc) = mop
            END DO
            index = i_loc
         END DO
         !!$       WRITE(*,*) i, list_elmts(1:index)
         IF (list_pts(i, 3) == 0) THEN  ! point au bord du bord periodique, ou sur une arete
            proc_min = MINVAL(partition(list_elmts(1:index)))
            partition(list_elmts(1)) = proc_min
         ELSE ! deux elements du bord periodique touchent le point
            IF (list_elmts(index) /= list_pts(i, 3)) THEN
               CALL error_Petsc('BUG in reassign_per_pts, wrong element')
            END IF
            proc1 = partition(list_elmts(1))
            proc2 = partition(list_elmts(2))
            partition(list_elmts(2:index - 1)) = MIN(proc1, proc2)
         END IF
      END DO

   END SUBROUTINE reassign_per_pts

END MODULE two_dim_metis_distribution
