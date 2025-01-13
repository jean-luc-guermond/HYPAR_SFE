MODULE Powell_Sabin_mesh
  USE def_type_mesh
CONTAINS

  SUBROUTINE convert_mesh_to_Powell_Sabin(dir, fil, list_dom, mesh_formatted, powsab_mesh)
    USE sub_plot
    USE input_data
    USE dir_nodes
    IMPLICIT NONE
    CHARACTER(len=200),    INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    LOGICAL,               INTENT(IN) :: mesh_formatted
    TYPE(mesh_type)                   :: mesh, powsab_mesh
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect
    INTEGER :: k, kd, nw, np, me, nws, mes, n, m, ms, n1, n2, m_new, next, mop, nop
    LOGICAL :: ok
    LOGICAL, POINTER, DIMENSION(:) :: ddir, virgin
    REAL(KIND=8), POINTER, DIMENSION(:) :: uu
    INTEGER, POINTER, DIMENSION(:) :: js_D
    REAL(KIND=8) :: a11, a12, a21, a22, b1, b2, t, dist, norm
    REAL(KIND=8):: x1(2), x2(2), xop1(2), xop2(2), xc(2) ,xopc(2)
    INTEGER :: capM, inew, m1, m2, mm1, mm2, ms1, ms2, nm0, nm1, nm2, mop1, in, &
         n_start, n_end, n_k1, n_k2, ns, ns1
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

    !===TEST
    DO ms = 1, mes
       m = neighs_lect(ms)
       ok=.FALSE.
       DO n = 1, nw
          mop = neigh_lect(n,m)
          IF (mop.GT.0) CYCLE
          n1 = MODULO(n,3)+1
          n2 = MODULO(n+1,3)+1
          !WRITE(*,*) jj_lect(n1,m), jj_lect(n2,m), jjs_lect(1,ms), jjs_lect(2,ms)
          xop1 = rr_lect(:,jj_lect(n1,m))
          xop2 = rr_lect(:,jj_lect(n2,m))
          x1 = rr_lect(:,jjs_lect(1,ms))
          x2 = rr_lect(:,jjs_lect(2,ms))
          norm = SUM(ABS(x1-x2))
          IF (SUM(ABS(x1-xop1))/norm .LE. 1d-7 .AND. SUM(ABS(x2-xop2))/norm.LE. 1d-7 ) THEN
             ok=.true.
          ELSE IF (SUM(ABS(x2-xop1))/norm .LE. 1d-7 .AND. SUM(ABS(x1-xop2))/norm.LE. 1d-7 )  THEN
             ok =.true.
          ELSE
          END IF
       END DO
       IF (.NOT.ok) THEN
          WRITE(*,*) ' BUG in convert_mesh_to_Powell_Sabin: flag-1', ms, neighs_lect(ms)
       END IF
    END DO

    !==================
    !===Create HCT mesh
    !==================
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
          m_new = (m-1)*nw + MODULO(n,3) + 1
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

    !===========================
    !===Create Poweel Sabin mesh
    !===========================
    powsab_mesh%gauss%n_w = nw
    powsab_mesh%gauss%n_ws = 2
    powsab_mesh%mes = 2*mesh%mes
    powsab_mesh%me = 2*mesh%me
    powsab_mesh%np = mesh%np + (3*mesh%me+mesh%mes)/2 - mesh%me
    ALLOCATE (powsab_mesh%jj(nw,powsab_mesh%me), &
         powsab_mesh%neigh(nw,powsab_mesh%me), &
         powsab_mesh%i_d(powsab_mesh%me), &
         powsab_mesh%rr(kd,powsab_mesh%np), &
         virgin(powsab_mesh%me))

    powsab_mesh%rr(:,1:mesh%np) = mesh%rr(:,1:mesh%np)
    virgin = .true.
    inew = mesh%np
    DO m = 1 , mesh%me
       m1 = 2*(m-1)+1
       m2 = 2*(m-1)+2

       powsab_mesh%i_d(m1) = mesh%i_d(m)
       powsab_mesh%i_d(m2) = mesh%i_d(m)

       powsab_mesh%jj(1,m1) = mesh%jj(1,m)
       powsab_mesh%jj(3,m1) = mesh%jj(3,m)
       powsab_mesh%jj(1,m2) = mesh%jj(2,m)
       powsab_mesh%jj(3,m2) = mesh%jj(3,m)

       !===Create  powsab_mesh%neigh(3,m1), powsab_mesh%neigh(3,m2)
       capM = (m-1-MODULO(m-1,3))/3 + 1
       nm0 = m - 3*(capM-1)
       nm1 = MODULO(nm0,3) +1
       nm2 = MODULO(nm0+1,3) +1
       mm1 = 3*(capM-1)+ nm1
       mm2 = 3*(capM-1)+ nm2
       !write(*,*) 'neigh', 2*(mm2-1) + 2, 2*(mm1-1) + 1, capM, nm1, nm2

       powsab_mesh%neigh(1,m1) = m2
       powsab_mesh%neigh(2,m1) = 2*(mm2-1) + 2   
       powsab_mesh%neigh(1,m2) = m1
       powsab_mesh%neigh(2,m2) = 2*(mm1-1) + 1

       mop = mesh%neigh(3,m) !===index 3 is the center node
       IF (mop.EQ.0) THEN
          inew = inew + 1
          !write(*,*) 'boundary', inew, powsab_mesh%np
          powsab_mesh%jj(2,m2) = inew
          powsab_mesh%jj(2,m1) = inew

          !===Create middle point
          x1 = mesh%rr(:,mesh%jj(1,m))
          x2 = mesh%rr(:,mesh%jj(2,m))
          powsab_mesh%rr(:,inew) = (x1+x2)/2

          powsab_mesh%neigh(3,m1) = 0
          powsab_mesh%neigh(3,m2) = 0
       ELSE
          IF (virgin(mop)) THEN
             inew = inew + 1
             in = inew
             !write(*,*) 'internal', in, powsab_mesh%np
          ELSE
             mop1 = 2*(mop-1) + 1 
             in = powsab_mesh%jj(2,mop1)
          END IF
          virgin(m) = .FALSE.
          powsab_mesh%jj(2,m2) = in
          powsab_mesh%jj(2,m1) = in

          xc   = mesh%rr(:,mesh%jj(3,m))
          x1   = mesh%rr(:,mesh%jj(1,m))
          x2   = mesh%rr(:,mesh%jj(2,m))
          xopc = mesh%rr(:,mesh%jj(3,mop))
          xop1 = mesh%rr(:,mesh%jj(1,mop))
          xop2 = mesh%rr(:,mesh%jj(2,mop))

          !===Find middle point
          !===Solve  x2 - t(x1-x2) = xopc + top(xc-xopc)
          !===t(x1-x2) + top(xc-xopc) = x2 - xopc 
          a11 = x1(1)-x2(1)
          a12 = x1(2)-x2(2)
          a21 = xc(1)-xopc(1)
          a22 = xc(2)-xopc(2)
          b1  = x2(1) - xopc(1) 
          b2  = x2(2) - xopc(2) 
          t   = (a22*b1-a21*b2)/(a11*a22-a21*a12)
          !top = (a11*b2-a12*b1)/(a11*a22-a21*a12)
          powsab_mesh%rr(:,in) = x2 -t*(x1-x2)

          !===Create  powsab_mesh%neigh(3,m1), powsab_mesh%neigh(3,m2)
          IF (SUM(ABS(x1-xop1))/SUM(ABS(x1-x2)) .LE. 1d-7) THEN 
             powsab_mesh%neigh(3,m1) = 2*(mop-1) + 1 !===Indices 1 and 1 on both sides
             powsab_mesh%neigh(3,m2) = 2*(mop-1) + 2
          ELSE
             powsab_mesh%neigh(3,m1) = 2*(mop-1) + 2 !===Swap indices
             powsab_mesh%neigh(3,m2) = 2*(mop-1) + 1
          END IF
       END IF
    END DO

    !===Boundary arrays
    ALLOCATE(powsab_mesh%jjs(nws,2*mesh%mes), powsab_mesh%neighs(2*mesh%mes), powsab_mesh%sides(2*mesh%mes))
    DO ms = 1, mesh%mes
       m = mesh%neighs(ms)
       DO n = 1, nw
          mop = mesh%neigh(n,m)
          IF (mop.LE.0) EXIT
       END DO
       ms1 = 2*(ms-1) + 1
       ms2 = 2*(ms-1) + 2
       m1 = 2*(m-1) + 1
       m2 = 2*(m-1) + 2
       powsab_mesh%neighs(ms1) = m1
       powsab_mesh%neighs(ms2) = m2
       powsab_mesh%sides(ms1)  = mesh%sides(ms)
       powsab_mesh%sides(ms2)  = mesh%sides(ms)
       n1 = MODULO(n,3)+1
       n2 = MODULO(n+1,3)+1
       xop1   = mesh%rr(:,mesh%jj(n1,m))
       xop2   = mesh%rr(:,mesh%jj(n2,m))
       x1 = mesh%rr(:,mesh%jjs(1,ms))
       x2 = mesh%rr(:,mesh%jjs(2,ms))
       IF (SUM(ABS(x1-xop1))/SUM(ABS(x1-x2)) .LE. 1d-7) THEN
          powsab_mesh%jjs(1,ms1)  = mesh%jjs(1,ms)
          powsab_mesh%jjs(1,ms2)  = mesh%jjs(2,ms)
       ELSE IF (SUM(ABS(x1-xop2))/SUM(ABS(x1-x2)) .LE. 1d-7) THEN
          powsab_mesh%jjs(1,ms1)  = mesh%jjs(2,ms)
          powsab_mesh%jjs(1,ms2)  = mesh%jjs(1,ms)
       ELSE
          WRITE(*,*) ' BUG in convert_mesh_to_Powell_Sabin: flag1'
          write(*,*) m1, m2, powsab_mesh%me
          WRITE(*,*) mesh%rr(:,mesh%jj(n1,m))
          WRITE(*,*) mesh%rr(:,mesh%jj(n2,m))
          WRITE(*,*) mesh%rr(:,mesh%jjs(1,mes))
          WRITE(*,*) mesh%rr(:,mesh%jjs(2,mes))
          STOP
       END IF
       powsab_mesh%jjs(2,ms1)  = powsab_mesh%jj(2,m1)
       powsab_mesh%jjs(2,ms2)  = powsab_mesh%jj(2,m2)
    END DO

   

    !===Lines below are not useful 
    !===Test whether mesh%neigh is well done (these lines can be removed)
    DO m = 1, powsab_mesh%me
       DO n = 1, nw
          mop = powsab_mesh%neigh(n,m)
          IF (mop==0) CYCLE
          n1 = MODULO(n,nw)+1
          n2 = MODULO(n+1,nw)+1
          ok = .FALSE.
          DO nop = 1, nw
             IF (ABS(powsab_mesh%jj(nop,mop)-powsab_mesh%jj(n1,m)).NE.0 .AND. &
                  ABS(powsab_mesh%jj(nop,mop)-powsab_mesh%jj(n2,m)).NE.0) THEN
                OK = .true.
                EXIT
             END IF
          END DO
          IF (.NOT.ok) THEN
             WRITE(*,*) ' BUG in convert_mesh_to_Powell_Sabin'
             STOP
          END IF
       END DO
    END DO

    virgin = .TRUE.
    in = 0
    DO m = 1, powsab_mesh%me
       virgin(m) = .FALSE.
       DO n = 1, 3
          mop = powsab_mesh%neigh(n,m)
          IF (mop==0) CYCLE !Edge on boundary
          IF (.NOT.virgin(mop)) CYCLE !Edge already done
          in = in + 1 !New edge
       END DO
    END DO
    IF (SIZE(powsab_mesh%rr,1)==2) THEN
       IF (in/=(3*powsab_mesh%me - powsab_mesh%mes)/2) THEN
          WRITE(*,*) ' BUG in prep_interfaces, internal edge/=(3*mesh%me + mesh%mes)/2'
          WRITE(*,*) ' internal edges ', in, (3*powsab_mesh%me - powsab_mesh%mes)/2
          WRITE(*,*) ' mesh%mes ', powsab_mesh%mes, ' mesh%me ',powsab_mesh%me
          STOP
       END IF
    END IF
    DO m = 1, powsab_mesh%me !===loop on the elements
       IF (MINVAL(powsab_mesh%neigh(:,m)).GT.0) CYCLE
       DO k = 1, 3 !===loop on the nodes (sides) of the element
          mop = powsab_mesh%neigh(k,m)
          IF (mop .LE. 0) EXIT
       END DO
       in = powsab_mesh%jj(k,m)
       n_k1 = MODULO(k,3) + 1
       n_k2 = MODULO(k+1,3) + 1
       n1 = powsab_mesh%jj(n_k1,m)
       n2 = powsab_mesh%jj(n_k2,m)
       IF (n1<n2) THEN !===Go from lowest global index to highest global index
          n_start = n1
          n_end   = n2
       ELSE
          n_start = n2
          n_end   = n1
       END IF

       DO ms = 1, powsab_mesh%mes
          DO ns = 1, 2
             !WRITE(*,*) powsab_mesh%rr(:,powsab_mesh%jjs(ns,ms))
             dist = SUM(ABS(powsab_mesh%rr(:,n1)-powsab_mesh%rr(:,powsab_mesh%jjs(ns,ms))))/ &
                  SUM(ABS(powsab_mesh%rr(:,n1)-powsab_mesh%rr(:,in)))
             !write(*,*) 'dist', dist
             IF (dist.LE.1.d-7) THEN
                ns1 = MODULO(ns,   SIZE(powsab_mesh%jjs,1)) + 1 
                dist = SUM(ABS(powsab_mesh%rr(:,n2)-powsab_mesh%rr(:,powsab_mesh%jjs(ns1,ms))))/ &
                     SUM(ABS(powsab_mesh%rr(:,n1)-powsab_mesh%rr(:,in)))
                IF (dist.LE.1.d-7) THEN
                   GO TO 100
                END IF
             END IF
          END DO
       END DO
       WRITE(*,*) ' BUG in create_iso_grid'
       !===Algorithm not designed yet for internal interfaces
100    continue
    END DO

  END SUBROUTINE convert_mesh_to_Powell_Sabin
END MODULE POWELL_SABIN_MESH
