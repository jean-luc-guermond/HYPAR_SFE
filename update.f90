MODULE update
#include "petsc/finclude/petsc.h"
  USE petsc
  USE matrix_type
  USE def_type_mesh
  PUBLIC:: euler
  PRIVATE
  REAL(KIND=8), DIMENSION(:),   POINTER :: lumped, norm_cij
  REAL(KIND=8), DIMENSION(:,:), POINTER :: nij1, nij2
  REAL(KIND=8), DIMENSION(:,:), POINTER :: norm_cij_s, nij1_s, nij2_s
  REAL(KIND=8)                          :: eps=1.d-15
!#include "finclude/petsc.h"
  Mat, DIMENSION(2)                 :: cij
  Mat                               :: dij
  Mat                               :: mass
  Vec                               :: xx, x1vec, x2vec, x3vec, x4vec, xghost, vec_one
  PetscErrorCode                    :: ierr
CONTAINS 

  SUBROUTINE contruct_matrices
    USE st_matrix
    USE mesh_handling
    USE fem_M
    USE st_matrix
    USE solve_petsc
    USE my_util
    IMPLICIT NONE
    INTEGER, POINTER, DIMENSION(:)  :: ifrom  ! for ghost structure
 
    !===Create ghost structure
    CALL create_my_ghost(mesh,LA,ifrom)
    CALL VecCreateGhost(PETSC_COMM_WORLD, mesh%dom_np, & 
         PETSC_DETERMINE, SIZE(ifrom), ifrom, xx, ierr)
    CALL VecGhostGetLocalForm(xx, xghost, ierr)

    !===Duplicate
    CALL VecDuplicate(xx, x1vec, ierr)
    CALL VecDuplicate(xx, x2vec, ierr)
    CALL VecDuplicate(xx, x3vec, ierr)
    CALL VecDuplicate(xx, x4vec, ierr)
    CALL VecDuplicate(xx, vec_one, ierr)

    !===Create matrices
    CALL create_local_petsc_matrix(PETSC_COMM_WORLD, LA, mass, clean=.FALSE.)
    CALL create_local_petsc_matrix(PETSC_COMM_WORLD, LA, cij(1), clean=.FALSE.)
    CALL create_local_petsc_matrix(PETSC_COMM_WORLD, LA, cij(2), clean=.FALSE.)
    CALL create_local_petsc_matrix(PETSC_COMM_WORLD, LA, dij, clean=.FALSE.)

    !===mass
    CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, mass)

    !===lumped
    ALLOCATE(lumped(mesh%np))
    call VecSet(vec_one,1.d0,ierr)
    CALL MatMult(mass,vec_one,xx,ierr)
    CALL VecGhostUpdateBegin(xx,INSERT_VALUES,SCATTER_FORWARD,ierr) 
    CALL VecGhostUpdateEnd(xx,INSERT_VALUES,SCATTER_FORWARD,ierr)
    CALL extract(xghost,1,1,LA,lumped)

    !===cij(1), cij(2)
    CALL construct_cij

  END SUBROUTINE contruct_matrices

  SUBROUTINE construct_cij
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w*mesh%gauss%n_w) :: vv1_rowise, vv2_rowise
    !REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: test1, test2
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idx
    INTEGER      :: m, ni, nj, l

!TEST
    !CALL MatSetOption (cij(1), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
    !CALL MatSetOption (cij(2), MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
!TEST
    CALL MatZeroEntries (cij(1), ierr)
    CALL MatZeroEntries (cij(2), ierr)
    DO m = 1, mesh%dom_me
       idx =  LA%loc_to_glob(1,mesh%jj(:,m)) - 1
       l = 0
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w 
             l = l + 1
             vv1_rowise(l) = -SUM(mesh%gauss%dw(1,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
             vv2_rowise(l) = -SUM(mesh%gauss%dw(2,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
!TEST
             !test1(ni,nj) = -SUM(mesh%gauss%dw(1,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
             !test2(ni,nj) = -SUM(mesh%gauss%dw(2,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
!TEST
          ENDDO
       ENDDO
       CALL MatSetValues(cij(1), mesh%gauss%n_w, idx, mesh%gauss%n_w, idx, vv1_rowise, ADD_VALUES, ierr)
       CALL MatSetValues(cij(2), mesh%gauss%n_w, idx, mesh%gauss%n_w, idx, vv2_rowise, ADD_VALUES, ierr)
       !CALL MatSetValues(cij(1), mesh%gauss%n_w, idx, mesh%gauss%n_w, idx, test1, ADD_VALUES, ierr)
       !CALL MatSetValues(cij(2), mesh%gauss%n_w, idx, mesh%gauss%n_w, idx, test2, ADD_VALUES, ierr)
    ENDDO
    CALL MatAssemblyBegin(cij(1),MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd  (cij(1),MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyBegin(cij(2),MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd  (cij(2),MAT_FINAL_ASSEMBLY,ierr)

  END SUBROUTINE construct_cij

  SUBROUTINE euler(un,unext,tnext)
    USE mesh_handling
    USE boundary_conditions
    USE input_data
    USE my_util
    USE st_matrix
    USE dir_nodes_petsc
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np), INTENT(OUT) :: unext
    REAL(KIND=8),                     INTENT(IN)  :: tnext
    REAL(KIND=8), DIMENSION(2,mesh%np)            :: ff
    REAL(KIND=8), DIMENSION(mesh%np)              :: rk
    REAL(KIND=8), DIMENSION(SIZE(js_D_glob))      :: un_D
    LOGICAL, SAVE :: once=.true.

    IF (once) THEN
       CALL contruct_matrices
       once=.FALSE.
    END IF

    !===Viscosity
    CALL compute_dij(un)
    !CALL MatZeroEntries(dij, ierr)
    IF (inputs%if_viscous) THEN
       !===Compute rk
       ff=flux(un)
       CALL vector(ff(1,:),x1vec,'insert')
       CALL MatMult(cij(1),x1vec,x3vec,ierr) !===cij(1)*ff(1,:)
       CALL vector(ff(2,:),x2vec,'insert')
       CALL MatMultAdd(cij(2),x2vec,x3vec,x4vec,ierr) !==cij(2)*ff(2,:)+cij(1)*ff(1,:)

       CALL vector(un,x1vec,'insert')
       CALL MatMultAdd(dij,x1vec,x4vec,xx,ierr) !===xx = dij*un+cij(2)*ff(2,:)+cij(1)*ff(1,:)

       rk = inputs%dt/lumped
       CALL vector(rk,x2vec,'insert')
       CALL VecPointwiseMult(xx, x2vec, xx, ierr) !===xx = xx*inputs%dt/lumped
       CALL VecAXPY(xx, 1.d0, x1vec, ierr) !===xx = xx + un

       !===Boundary conditions
       un_D =  sol_anal(mesh%rr(:,js_D_loc),tnext)
       CALL dirichlet_rhs(js_D_glob,un_D,xx) 

       CALL VecGhostUpdateBegin(xx,INSERT_VALUES,SCATTER_FORWARD,ierr) 
       CALL VecGhostUpdateEnd(xx,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(xghost,1,1,LA,unext)

       !===Check maximum principle
       CALL check_max_principle(unext,un)
       RETURN
    ELSE
       CALL error_petsc('Error: EV not implemented yet')
    END IF
  END SUBROUTINE euler

  SUBROUTINE vector(uu,xx,what)
      USE mesh_handling
      USE my_util
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(mesh%np), INTENT(IN)  :: uu
      CHARACTER(LEN=6),                 INTENT(IN)  :: what           
      INTEGER,      DIMENSION(mesh%np)              :: idxm
      INTEGER                                       :: i
!#include "finclude/petsc.h"
      Vec                                           :: xx
      DO i = 1, mesh%np
          idxm(i) = LA%loc_to_glob(1,i) -1
      END DO
      SELECT CASE (what)
      CASE('insert')
         CALL VecSetValues(xx, mesh%np, idxm, uu, INSERT_VALUES, ierr)
      CASE('add')
         CALL VecSetValues(xx, mesh%np, idxm, uu, ADD_VALUES, ierr)
      CASE DEFAULT
         CALL error_petsc('what is wrong')
      END SELECT
 END SUBROUTINE vector

 SUBROUTINE compute_dij(un)
   USE mesh_handling
   USE input_data
   USE boundary_conditions
   USE my_util
   USE st_matrix
   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(mesh%np), INTENT(IN)  :: un
   INTEGER,      DIMENSION(0:0)                  :: idx, jdx
   REAL(KIND=8), DIMENSION(0:0)                  :: dij_loc
   REAL(KIND=8), DIMENSION(2,mesh%np)            :: fpu
   REAL(KIND=8), DIMENSION(mesh%np)              :: diag_dij
   REAL(KIND=8) :: lambda, fpur, fpul, ul, ur, umax
   INTEGER :: i, j, k, m, n, ni, nj, ms,ierr, rank
   LOGICAL :: Kuz=.FALSE., if_test_diagonal=.FALSE.
   LOGICAL, SAVE :: once=.true.

   CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
   CALL MatZeroEntries(dij, ierr)
!!$
!!$   IF (Kuz) THEN
!!$   END IF
!!$

   IF (once) THEN
      once = .FALSE.
      ALLOCATE(nij1(1,mesh%mi),nij2(1,mesh%mi),norm_cij(mesh%mi))
      ALLOCATE(nij1_s(2,mesh%mes),nij2_s(2,mesh%mes),norm_cij_s(2,mesh%mes))

      DO ms = 1, mesh%mi
         i = mesh%jjsi(1,ms)
         j = mesh%jjsi(2,ms)
         k = MIN(i,j)
         j = MAX(i,j)
         i = k ! Make sure row i is on the processor (1\le i \le mesh%dom_np)
         idx(0) = LA%loc_to_glob(1,i) - 1
         jdx(0) = LA%loc_to_glob(1,j) - 1
         IF (i.LE.mesh%dom_np) THEN
            CALL MatGetValues(cij(1),1,idx,1,jdx,nij1(1,ms:ms),ierr)
            CALL MatGetValues(cij(2),1,idx,1,jdx,nij2(1,ms:ms),ierr)
            norm_cij(ms) = SQRT(nij1(1,ms)**2+nij2(1,ms)**2)
            nij1(1,ms) = nij1(1,ms)/norm_cij(ms)
            nij2(1,ms) = nij2(1,ms)/norm_cij(ms)
         END IF
      END DO

      DO ms = 1, mesh%mes
         m = mesh%neighs(ms)
         DO n = 1, 3
            if (mesh%neigh(n,m)==0) EXIT
         END DO
         ni = MODULO(n,3)+1
         nj = MODULO(n+1,3)+1
         nij1_s(1,ms) = -SUM(mesh%gauss%dw(1,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
         nij2_s(1,ms) = -SUM(mesh%gauss%dw(2,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
         norm_cij_s(1,ms) = SQRT(nij1_s(1,ms)**2+nij2_s(1,ms)**2)
         nij1_s(1,ms) = nij1_s(1,ms)/norm_cij_s(1,ms)
         nij2_s(1,ms) = nij2_s(1,ms)/norm_cij_s(1,ms)

         !===Second time
         n = nj 
         nj = ni
         ni = n
         nij1_s(2,ms) = -SUM(mesh%gauss%dw(1,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
         nij2_s(2,ms) = -SUM(mesh%gauss%dw(2,nj,:,m)*mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
         norm_cij_s(2,ms) = SQRT(nij1_s(2,ms)**2+nij2_s(2,ms)**2)
         nij1_s(2,ms) = nij1_s(2,ms)/norm_cij_s(2,ms)
         nij2_s(2,ms) = nij2_s(2,ms)/norm_cij_s(2,ms)
      END DO
   END IF

   CALL MatZeroEntries(dij, ierr)

   fpu=flux_prime(un)
   DO ms = 1, mesh%mi
      i = mesh%jjsi(1,ms)
      j = mesh%jjsi(2,ms)
      k = MIN(i,j)
      j = MAX(i,j)
      i = k ! Make sure row i is on the processor (1\le i \le mesh%dom_np)
      IF (i.LE.mesh%dom_np) THEN 
         idx(0) = LA%loc_to_glob(1,i) - 1
         jdx(0) = LA%loc_to_glob(1,j) - 1
         fpur = fpu(1,j)*nij1(1,ms)+fpu(2,j)*nij2(1,ms)
         fpul = fpu(1,i)*nij1(1,ms)+fpu(2,i)*nij2(1,ms)
         ul=un(i)
         ur=un(j)
         CALL compute_lambda(ur,ul,fpur,fpul,nij1(1,ms),nij2(1,ms),lambda)
         dij_loc(0) =  norm_cij(ms)*lambda
         CALL MatSetValue(dij, idx, jdx, dij_loc, INSERT_VALUES, ierr)
         CALL MatSetValue(dij, jdx, idx, dij_loc, INSERT_VALUES, ierr)
      END IF
   ENDDO


   !===Boundary 
   DO ms = 1, mesh%mes
      m = mesh%neighs(ms)
      DO n = 1, 3
         if (mesh%neigh(n,m)==0) EXIT
      END DO
      !===First time
      ni = MODULO(n,3)+1
      nj = MODULO(n+1,3)+1
      i = mesh%jj(ni,m)
      j = mesh%jj(nj,m)
      idx(0) = LA%loc_to_glob(1,i) - 1
      jdx(0) = LA%loc_to_glob(1,j) - 1
      fpur = fpu(1,j)*nij1_s(1,ms)+fpu(2,j)*nij2_s(1,ms)
      fpul = fpu(1,i)*nij1_s(1,ms)+fpu(2,i)*nij2_s(1,ms)
      ul=un(i)
      ur=un(j)
      CALL compute_lambda(ur,ul,fpur,fpul,nij1_s(1,ms),nij2_s(1,ms),lambda)
      dij_loc(0) =  norm_cij_s(1,ms)*lambda

      !===Second time
      n = nj 
      nj = ni
      ni = n
      i = mesh%jj(ni,m)
      j = mesh%jj(nj,m)
      fpur = fpu(1,j)*nij1_s(2,ms)+fpu(2,j)*nij2_s(2,ms)
      fpul = fpu(1,i)*nij1_s(2,ms)+fpu(2,i)*nij2_s(2,ms)
      ul=un(i)
      ur=un(j)
      CALL compute_lambda(ur,ul,fpur,fpul,nij1_s(2,ms),nij2_s(2,ms),lambda)
      dij_loc(0) = MAX(dij_loc(0),norm_cij_s(2,ms)*lambda)

      CALL MatSetValue(dij, idx(0), jdx(0), dij_loc(0), INSERT_VALUES, ierr)
      CALL MatSetValue(dij, jdx(0), idx(0), dij_loc(0), INSERT_VALUES, ierr)
   END DO

   CALL MatAssemblyBegin(dij,MAT_FINAL_ASSEMBLY,ierr)
   CALL MatAssemblyEnd  (dij,MAT_FINAL_ASSEMBLY,ierr)

   !===Diagonal
   CALL MatMult(dij,vec_one,xx,ierr)
   CALL VecGhostUpdateBegin(xx,INSERT_VALUES,SCATTER_FORWARD,ierr) 
   CALL VecGhostUpdateEnd(xx,INSERT_VALUES,SCATTER_FORWARD,ierr)
   CALL extract(xghost,1,1,LA,diag_dij)
   DO i = 1, mesh%dom_np
      idx(0) = LA%loc_to_glob(1,i) - 1
      dij_loc(0) = -diag_dij(i)
      CALL MatSetValue(dij, LA%loc_to_glob(1,i)-1, LA%loc_to_glob(1,i)-1, &
           dij_loc(0),INSERT_VALUES, ierr)
   ENDDO
   CALL MatAssemblyBegin(dij,MAT_FINAL_ASSEMBLY,ierr)
   CALL MatAssemblyEnd  (dij,MAT_FINAL_ASSEMBLY,ierr)


   !===Test digonal terms and transpose
   IF (if_test_diagonal) THEN
      CALL MatMult(dij,vec_one,xx,ierr)
      IF (umax.GT.1.d-15) THEN
         WRITE(*,*) rank, umax
         CALL error_petsc('BUG: diagonal is not zero')
      END IF
      CALL MatTranspose(dij, MAT_REUSE_MATRIX, dij, ierr)
      CALL MatMult(dij,vec_one,xx,ierr)
      CALL VecNorm(xx,NORM_INFINITY,umax,ierr)
      IF (umax.GT.1.d-15) THEN
         WRITE(*,*) rank, umax
         CALL error_petsc('BUG: diagonal is not zero on transpose')
      END IF
      CALL MatTranspose(dij, MAT_REUSE_MATRIX, dij, ierr)
   END IF

 END SUBROUTINE compute_dij
  
  SUBROUTINE check_max_principle(uu,un)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu, un
!!$    INTEGER :: i
!!$    REAL(KIND=8) :: unmin, unmax, eps=1.d-15
!!$    DO i = 1, mesh%np
!!$       unmax = MAXVAL(un(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
!!$       unmin = MINVAL(un(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
!!$       IF (uu(i)<unmin-eps) THEN
!!$          WRITE(*,*) ' Minimum violation ', uu(i),unmin, ABS(uu(i)-unmin)/unmin
!!$       ELSE IF (uu(i)>unmax+eps) THEN
!!$          WRITE(*,*) ' Maximum violation ',  uu(i),unmax, ABS(uu(i)-unmax)/unmax
!!$       END IF
!!$    END DO
  END SUBROUTINE check_max_principle

  SUBROUTINE compute_lambda(ur,ul,fpur,fpul,nij1,nij2,lambda)
    USE input_data
    USE boundary_conditions
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT) :: lambda
    REAL(KIND=8)              :: ur, ul, fpur, fpul, nij1, nij2
    REAL(KIND=8)              :: umin, umax, theta, fur, ful
    IF (inputs%type_test==1) THEN  !Linear transport
       lambda = MAX(ABS(fpul),ABS(fpur))
    ELSE IF (inputs%type_test==2) THEN  !KPP
       umax=MAX(ur,ul)
       umin=MIN(ur,ul)
       theta = ATAN2(nij2,nij1)
       IF (floor((umin+theta)/pi).NE.floor((umax+theta)/pi+eps)) THEN
          lambda=1.d0
       ELSE IF (fpur.GE.fpul-eps) THEN !Expansion
          lambda = MAX(ABS(fpul),ABS(fpur))
       ELSE !Shock
          fur = SIN(ur+theta)
          ful = SIN(ul+theta)
          lambda = ABS((fur-ful)/(ul-ur))
       END IF
    ELSE IF (inputs%type_test==3) THEN !Burgers
       IF (fpur>fpul) THEN !Expansion
          lambda = MAX(ABS(fpur),ABS(fpul))
       ELSE !Shock
          lambda = ABS(ur+ul)/2 
       END IF
    ELSE
       CALL error_petsc('BUG: type_test not programmed yet')
    END IF
    RETURN
  END SUBROUTINE compute_lambda
END MODULE update
