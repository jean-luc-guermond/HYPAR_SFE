MODULE HCT_mesh
  USE def_type_mesh
CONTAINS

  SUBROUTINE convert_mesh_to_HCT(dir, fil, list_dom, mesh_formatted, mesh)
    USE sub_plot
    USE input_data
    USE dir_nodes
    IMPLICIT NONE
    CHARACTER(len=200),    INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    LOGICAL,               INTENT(IN) :: mesh_formatted
    TYPE(mesh_type)                   :: mesh
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect
    INTEGER :: k, kd, nw, np, me, nws, mes, n, m, ms, n1, n2, m_new, next, mop, nop, ns
    LOGICAL :: ok
    REAL(KIND=8) :: norm
    REAL(KIND=8):: x1(2), x2(2), xop1(2), xop2(2)
    LOGICAL, POINTER, DIMENSION(:) :: ddir
    REAL(KIND=8), POINTER, DIMENSION(:) :: uu
    INTEGER, POINTER, DIMENSION(:) :: js_D

    IF (mesh_formatted) THEN
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    ELSE
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    END IF

    !===Read P1 mesh
    IF (mesh_formatted) THEN
       READ  (30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
    END IF

    !===Only 2D
    kd = 2
    IF (nw.NE.3) THEN
       WRITE(*,*) ' BUG in convert_mesh_to_HCT, nw.NE.3'
       STOP
    END IF

    !===Read mesh
    ALLOCATE (jj_lect(nw,me),neigh_lect(nw,me),i_d_lect(me))
    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
       END DO
    ELSE
       READ(30) jj_lect, neigh_lect, i_d_lect
    END IF

    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))
    IF (mesh_formatted) THEN
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
       END DO
    ELSE
       READ(30) jjs_lect, neighs_lect, sides_lect
    END IF

    ALLOCATE(rr_lect(kd,np))
    IF (mesh_formatted) THEN
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) rr_lect
    END IF
    CLOSE(30)

    !===Create HCT mesh
    mesh%gauss%n_w = nw
    mesh%gauss%n_ws = 2
    mesh%mes = mes
    mesh%me = nw*me
    mesh%np = np + me
    ALLOCATE (mesh%jj(nw,mesh%me), mesh%neigh(nw,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, me
       m_new = (m-1)*nw
       mesh%i_d(m_new+1:m_new+nw) = i_d_lect(m)
       DO n = 1, nw
          m_new = m_new + 1 !===new element
          next = MODULO(n-1,nw)+1
          mesh%jj(1,m_new) = jj_lect(next,m)
          next = MODULO(n,nw)+1
          mesh%jj(2,m_new) = jj_lect(next,m)
          mesh%jj(3,m_new) = np + m

          next = MODULO(n,nw)+1
          mesh%neigh(1,m_new) = (m-1)*nw + next
          next = MODULO(n+1,nw)+1
          mesh%neigh(2,m_new) = (m-1)*nw + next
       END DO
    END DO
    !===Now construct neigh(3,.)
    DO m = 1, me
       DO n = 1, nw
          n1 = MODULO(n,nw)+1
          n2 = MODULO(n+1,nw)+1
          mop = neigh_lect(n,m)
          IF (mop.LE.0) THEN
             m_new = (m-1)*nw + MODULO(n,3)+1
             mesh%neigh(3,m_new) = 0
             CYCLE
          END IF
          DO nop = 1, nw
             IF (ABS(jj_lect(nop,mop)-jj_lect(n1,m)).NE.0 .AND. &
                  ABS(jj_lect(nop,mop)-jj_lect(n2,m)).NE.0) THEN
                EXIT !===nop is the opposite index
             END IF
          END DO
          m_new = (m-1)*nw + MODULO(n,3)+1
          mesh%neigh(3,m_new) = (mop-1)*nw + MODULO(nop,3) + 1
       END DO
    END DO

    !===Boundary arrays
    ALLOCATE(mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), mesh%sides(mesh%mes))
    mesh%sides = sides_lect
    mesh%jjs = jjs_lect
    DO ms = 1, mesh%mes
       m = neighs_lect(ms)
       ok = .false.
       DO n = 1, nw
          n1 = MODULO(n,nw)+1
          n2 = MODULO(n+1,nw)+1
          mop = neigh_lect(n,m)
          IF (mop.LE.0) THEN
             IF ((ABS(jjs_lect(1,ms)-jj_lect(n1,m))==0 .AND. &
                  ABS(jjs_lect(2,ms)-jj_lect(n2,m))==0) .OR. &
                  (ABS(jjs_lect(2,ms)-jj_lect(n1,m))==0 .AND. &
                  ABS(jjs_lect(1,ms)-jj_lect(n2,m))==0)) THEN
                ok = .true.
                k = n1+n2
                SELECT CASE(k)
                CASE(3)
                   ns = 3
                CASE(5)
                   ns = 1
                CASE(4)
                   ns =2
                END SELECT
                m_new = (m-1)*nw + MODULO(ns,3) + 1
                mesh%neighs(ms) = m_new
             END IF
          END IF
       END DO
       IF (.NOT.ok) THEN
          WRITE(*,*) ' BUG in convert_mesh_to_Powell_Sabin: flag-2'
       END IF
    END DO

!!$    !===Boundary arrays
!!$    ALLOCATE(mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), mesh%sides(mesh%mes))
!!$    mesh%sides = sides_lect
!!$    mesh%jjs = jjs_lect
!!$    DO ms = 1, mesh%mes
!!$       m = neighs_lect(ms)
!!$       DO n = 1, nw
!!$          mop = neigh_lect(n,m)
!!$          IF (mop.LE.0) THEN
!!$             m_new = (m-1)*nw + MODULO(n,3)+1
!!$             mesh%neighs(ms) = m_new 
!!$          END IF
!!$       END DO
!!$    END DO

    !====Grid points
    ALLOCATE(mesh%rr(kd,mesh%np))
    mesh%rr(:,1:np) = rr_lect
    DO m = 1, me
       n = np + m
       DO k = 1, kd
          mesh%rr(k,n) = SUM(rr_lect(k,jj_lect(1:nw,m)))/nw
       END DO
    END DO

    !===Test
    DO ms = 1, mesh%mes
       m = mesh%neighs(ms)
       !write(*,*) ' mesh%neighs(ms)', m
       !write(*,*) mesh%jj(:,m)
       ok=.false.
       DO n = 1, nw
          mop = mesh%neigh(n,m)
          !write(*,*) 'm, n, mop', m, n, mop
          !IF (mop.GT.0) write(*,*) ' mop >0, mesh%jj(:,mop)', mesh%jj(:,mop)
          IF (mop.GT.0) CYCLE
          n1 = MODULO(n,3)+1
          n2 = MODULO(n+1,3)+1
          !WRITE(*,*) mesh%jj(n1,m), mesh%jj(n2,m), mesh%jjs(1,ms), mesh%jjs(2,ms)
          xop1   = mesh%rr(:,mesh%jj(n1,m))
          xop2   = mesh%rr(:,mesh%jj(n2,m))
          x1 = mesh%rr(:,mesh%jjs(1,ms))
          x2 = mesh%rr(:,mesh%jjs(2,ms))
          norm = SUM(ABS(x1-x2))
          IF ((SUM(ABS(x1-xop1))/norm .LE. 1d-7 .AND. SUM(ABS(x2-xop2))/norm.LE. 1d-7)  .OR. &
               (SUM(ABS(x2-xop1))/norm .LE. 1d-7 .AND. SUM(ABS(x1-xop2))/norm.LE. 1d-7))  THEN
             ok =.true.
             EXIT
          END IF
       END DO
       IF (.NOT.ok) THEN
          WRITE(*,*) ' BUG in convert_mesh_to_Powell_Sabin: flag0', ms, m
          STOP
       END IF
    END DO



    !===Store mesh
    OPEN(30,FILE='HCT.FEM',FORM='unformatted')
    WRITE(30) mesh%np, mesh%gauss%n_w, mesh%me, mesh%gauss%n_ws, mesh%mes
    WRITE(30) mesh%jj, mesh%neigh, mesh%i_d
    WRITE(30) mesh%jjs, mesh%neighs, mesh%sides
    WRITE(30) mesh%rr
    CLOSE(30)

    !===Test whether mesh%neigh is well done (these lines can be removed)
    DO m = 1, mesh%me
       DO n = 1, nw
          mop = mesh%neigh(n,m)
          IF (mop==0) CYCLE
          n1 = MODULO(n,nw)+1
          n2 = MODULO(n+1,nw)+1
          ok = .FALSE.
          DO nop = 1, nw
             IF (ABS(mesh%jj(nop,mop)-mesh%jj(n1,m)).NE.0 .AND. &
                  ABS(mesh%jj(nop,mop)-mesh%jj(n2,m)).NE.0) THEN
                OK = .true.
                EXIT
             END IF
          END DO
          IF (.NOT.ok) THEN
             WRITE(*,*) ' BUG in convert_mesh_to_HCT'
             STOP
          END IF
       END DO
    END DO

  END SUBROUTINE convert_mesh_to_HCT
END MODULE HCT_mesh
