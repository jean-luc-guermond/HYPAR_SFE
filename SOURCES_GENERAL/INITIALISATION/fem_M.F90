MODULE fem_M
   USE def_type_mesh
#include "petsc/finclude/petsc.h"
   USE petsc
   USE input_data
CONTAINS

  SUBROUTINE qs_mass_diff_M (mesh, mass, visco, LA, matrix)
    !=================================================
    IMPLICIT NONE
    TYPE(mesh_type), TARGET :: mesh
    REAL(KIND = 8), INTENT(IN) :: mass, visco
    type(petsc_csr_LA) :: LA
    REAL(KIND = 8), DIMENSION(mesh%gauss%n_w, mesh%gauss%n_w) :: mat_loc

    INTEGER, DIMENSION(mesh%gauss%n_w) :: idxn
    INTEGER :: m, ni, nj
    REAL(KIND = 8), DIMENSION(mesh%gauss%l_G) :: al, bl
    Mat            :: matrix
    PetscErrorCode :: ierr
    CALL MatZeroEntries (matrix, ierr)

    DO m = 1, mesh%dom_me
       idxn = LA%loc_to_glob(1, mesh%jj(:, m)) - 1

       al = visco * mesh%gauss%rj(:, m)
       bl = mass * mesh%gauss%rj(:, m)
       DO nj = 1, mesh%gauss%n_w;
          DO ni = 1, mesh%gauss%n_w;
             mat_loc(nj, ni) = SUM((mesh%gauss%dw(1, nj, :, m) * mesh%gauss%dw(1, ni, :, m) &
                  + mesh%gauss%dw(2, nj, :, m) * mesh%gauss%dw(2, ni, :, m)) * al &
                  + mesh%gauss%ww(ni, :) * mesh%gauss%ww(nj, :) * bl)
          ENDDO
       ENDDO

       CALL MatSetValues(matrix, mesh%gauss%n_w, idxn, mesh%gauss%n_w, idxn, mat_loc, ADD_VALUES, ierr)
    ENDDO

    CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
    CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
  END SUBROUTINE qs_mass_diff_M

   SUBROUTINE qs_mass_diff_M_AL (uu_mesh, pp_mesh, mass, visco, lambda, LA, matrix)
      !=================================================
      USE basis_change
      USE st_matrix
      IMPLICIT NONE

      REAL(KIND = 8), DIMENSION(:, :), POINTER, SAVE :: aij_p1p2, aij_p2p3
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: aij
      LOGICAL :: once_p1p2 = .TRUE., once_p2p3 = .TRUE.
      TYPE(mesh_type) :: uu_mesh, pp_mesh
      REAL(KIND = 8), INTENT(IN) :: mass, visco, lambda
      type(petsc_csr_LA) :: LA
      REAL(KIND = 8), DIMENSION(uu_mesh%gauss%n_w, uu_mesh%gauss%n_w) :: mat_loc
      INTEGER, DIMENSION(uu_mesh%gauss%n_w) :: idxn
      INTEGER :: m, ni, nj, nl, nk, kk, kl, l
      REAL(KIND = 8), DIMENSION(uu_mesh%gauss%l_G) :: al, bl
      REAL(KIND = 8), DIMENSION(uu_mesh%gauss%k_d, uu_mesh%gauss%n_w, uu_mesh%gauss%l_G) :: dwm
      REAL(KIND = 8), DIMENSION(pp_mesh%gauss%n_w, uu_mesh%gauss%l_G) :: ww_c
      REAL(KIND = 8), DIMENSION(pp_mesh%gauss%n_w, pp_mesh%gauss%n_w) :: loc_mass_inv, loc_mass
      REAL(KIND = 8) :: divu_pl, divu_pk, hatK_over_K
      Mat            :: matrix
      PetscErrorCode :: ierr
      CALL MatZeroEntries (matrix, ierr)

      IF (inputs%type_fe==2) THEN
         IF (once_p1p2) THEN
            CALL p1_p2(aij_p1p2)
            once_p1p2 = .FALSE.
         END IF
         aij => aij_p1p2
      ELSE IF (inputs%type_fe==3) THEN
         IF (once_p2p3) THEN
            CALL p2_p3(aij_p2p3)
            once_p2p3 = .FALSE.
         END IF
         aij => aij_p2p3
      ELSE
         WRITE(*, *) 'BUG in Bt_Mminus_B, inputs&type_fe_u not correct'
         STOP
      END IF

      !===Construct w_c
      DO l = 1, uu_mesh%gauss%l_G
         DO ni = 1, pp_mesh%gauss%n_w
            ww_c(ni, l) = SUM(aij(ni, :) * uu_mesh%gauss%ww(:, l))
         END DO
      END DO



      !===Inverse local mass matrix

      hatK_over_K = 0.5d0 / SUM(pp_mesh%gauss%rj(:, 1))
      loc_mass = 0.d0
      DO ni = 1, pp_mesh%gauss%n_w
         DO nj = 1, pp_mesh%gauss%n_w
            DO l = 1, pp_mesh%gauss%l_G
               loc_mass(ni, nj) = loc_mass(ni, nj) + ww_c(ni, l) * ww_c(nj, l) * pp_mesh%gauss%rj(l, 1)
            END DO
         END DO
      END DO
      loc_mass = loc_mass * hatK_over_K
      loc_mass_inv = loc_mass
      CALL inverse_mat(loc_mass_inv, pp_mesh%gauss%n_w)
      write(*, *) ww_c

      !===Construct w_c
      DO l = 1, uu_mesh%gauss%l_G
         DO ni = 1, pp_mesh%gauss%n_w
            ww_c(ni, l) = SUM(aij(ni, :) * uu_mesh%gauss%ww(:, l))
         END DO
      END DO


      !===Building matrix
      DO m = 1, uu_mesh%dom_me
         idxn = LA%loc_to_glob(1, uu_mesh%jj(:, m)) - 1
         mat_loc = 0.d0
         al = visco * uu_mesh%gauss%rj(:, m)
         bl = mass * uu_mesh%gauss%rj(:, m)
         DO nj = 1, uu_mesh%gauss%n_w;
            DO ni = 1, uu_mesh%gauss%n_w;
               mat_loc(nj, ni) = SUM((uu_mesh%gauss%dw(1, nj, :, m) * uu_mesh%gauss%dw(1, ni, :, m) &
                    + uu_mesh%gauss%dw(2, nj, :, m) * uu_mesh%gauss%dw(2, ni, :, m)) * al &
                    + uu_mesh%gauss%ww(ni, :) * uu_mesh%gauss%ww(nj, :) * bl)
               !===BtMB = SUM_nl_nk Bt(nj,nl) * M-1(nl,nk) * B(nk,ni)
               !===bcs DG don't get out of cell m
               DO nk = 1, pp_mesh%gauss%n_w
                  divu_pk = 0.d0
                  DO kk = 1, 2
                     DO l = 1, pp_mesh%gauss%l_G
                        divu_pk = divu_pk + dwm(kk, ni, l) * ww_c(nk, l) * uu_mesh%gauss%rj(l, m)
                     END DO
                     DO nl = 1, pp_mesh%gauss%n_w
                        divu_pl = 0.d0
                        DO kl = 1, 2
                           DO l = 1, pp_mesh%gauss%l_G
                              divu_pl = divu_pl + dwm(kl, nj, l) * ww_c(nl, l) * uu_mesh%gauss%rj(l, m)
                           END DO
                           mat_loc(nj, ni) = mat_loc(nj, ni) + lambda * divu_pk * divu_pl * &
                                loc_mass_inv(nk, nl) * (0.5d0 / SUM(pp_mesh%gauss%rj(:, m)))
                        END DO
                     END DO
                  END DO
               END DO
            ENDDO
         ENDDO
         CALL MatSetValues(matrix, uu_mesh%gauss%n_w, idxn, uu_mesh%gauss%n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO
      CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
   END SUBROUTINE qs_mass_diff_M_AL

!    SUBROUTINE qs_diff_mass_vect_M (LA, mesh, visc, mass, lambda, matrix)
!      !=================================================
!      USE my_util
!      USE basis_change
!      USE st_matrix
!      IMPLICIT NONE
!      TYPE(mesh_type), INTENT(IN) :: mesh
!      REAL(KIND = 8), INTENT(IN) :: visc, mass, lambda
!      TYPE(petsc_csr_LA) :: LA
!      INTEGER, DIMENSION(mesh%gauss%n_w) :: jj_loc
!      INTEGER, DIMENSION(:), ALLOCATABLE :: idxm, idxn
!      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w, mesh%gauss%n_w, mesh%gauss%l_G) :: wwprod
!      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w, mesh%gauss%n_w) :: a22_loc, a11_loc, b12_loc, b21_loc
!      REAL(KIND = 8), DIMENSION(pp_mesh%gauss%n_w, pp_mesh%gauss%n_w) :: wwinvww
!      REAL(KIND = 8), DIMENSION(pp_mesh%gauss%n_w, uu_mesh%gauss%l_G) :: ww_c
!      REAL(KIND = 8), DIMENSION(pp_mesh%gauss%n_w, pp_mesh%gauss%n_w) :: loc_mass_inv, loc_mass
!      REAL(KIND = 8), DIMENSION(:, :), POINTER, SAVE :: aij_p1p2, aij_p2p3
!      REAL(KIND = 8), DIMENSION(:, :), POINTER :: aij
!      REAL(KIND = 8), DIMENSION(2, 2) :: dwdw
!
!      LOGICAL :: once_p1p2 = .TRUE., once_p2p3 = .TRUE.
!      REAL(KIND = 8) :: hatK_over_K
!
!      REAL(KIND = 8) :: rj, rjrj, mrj, vrj, vvtrj, vtrj
!      INTEGER :: m, l, ni, nj, i, j, iglob, jglob, k_max, n_w, ix, jx, ki, kj, nl, nk
!      REAL(KIND = 8), DIMENSION(2*mesh%gauss%n_w,2*mesh%gauss%n_w) :: mat_loc
!      REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w,mesh%gauss%n_w,2)  :: bmat, bpmat
!      INTEGER :: k_dim, n_b, k1
!      REAL(KIND=8) :: x, y, lambda_grad_div, lambda_dg
!      !#include "petsc/finclude/petsc.h"
!      Mat                                         :: matrix
!      PetscErrorCode                              :: ierr
!      CALL MatZeroEntries (matrix, ierr)
!      CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
!
!      k_max = 2 ! 2x2 Structure
!
!      IF (k_max/=SIZE(LA%loc_to_glob, 1)) THEN
!         CALL error_Petsc('BUG in qs_diff_mass_vect_petsc_M, k_max/=SIZE(LA%loc_to_glob,1)')
!      END IF
!
!      n_w = mesh%gauss%n_w
!      DO l = 1, mesh%gauss%l_G
!         DO ni = 1, n_w
!            DO nj = 1, n_w
!               wwprod(ni, nj, l) = mesh%gauss%ww(ni, l) * mesh%gauss%ww(nj, l)
!            END DO
!         END DO
!      END DO
!
!      ALLOCATE(idxm(k_max * n_w), idxn(k_max * n_w))
!
!      IF (inputs%type_fe==2) THEN
!         IF (once_p1p2) THEN
!            CALL p1_p2(aij_p1p2)
!            once_p1p2 = .FALSE.
!         END IF
!         aij => aij_p1p2
!      ELSE IF (inputs%type_fe==3) THEN
!         IF (once_p2p3) THEN
!            CALL p2_p3(aij_p2p3)
!            once_p2p3 = .FALSE.
!         END IF
!         aij => aij_p2p3
!      ELSE
!         WRITE(*, *) 'BUG in Bt_Mminus_B, inputs&type_fe_u not correct'
!         STOP
!      END IF
!
!      wwinvww = 0.d0
!      IF(inputs%if_pp_dg) THEN
!         !===Construct w_c
!         DO l = 1, uu_mesh%gauss%l_G
!            DO ni = 1, pp_mesh%gauss%n_w
!               ww_c(ni, l) = SUM(aij(ni, :) * uu_mesh%gauss%ww(:, l))
!            END DO
!         END DO
!
!         !===Inverse local mass matrix
!         hatK_over_K = 0.5d0 / SUM(pp_mesh%gauss%rj(:, 1))
!         loc_mass = 0.d0
!         DO ni = 1, pp_mesh%gauss%n_w
!            DO nj = 1, pp_mesh%gauss%n_w
!               DO l = 1, pp_mesh%gauss%l_G
!                  loc_mass(ni, nj) = loc_mass(ni, nj) + ww_c(ni, l) * ww_c(nj, l) * pp_mesh%gauss%rj(l, 1)
!               END DO
!            END DO
!         END DO
!
!         loc_mass = loc_mass * hatK_over_K
!         loc_mass_inv = loc_mass
!         CALL inverse_mat(loc_mass_inv, pp_mesh%gauss%n_w)
!         lambda_grad_div = 0.d0
!         lambda_dg = lambda
!      ELSE
!         loc_mass_inv = 0.d0
!         lambda_grad_div = lambda
!         lambda_dg = 0.d0
!      END IF
!
!      k_dim = 2
!      n_b = k_dim
!      DO m = 1, mesh%dom_me
!         bmat = 0.d0
!         bpmat = 0.d0
!         IF(inputs%if_pp_dg) THEN
!            DO ni = 1, pp_mesh%gauss%n_w
!               DO nj = 1, mesh%gauss%n_w
!                  DO kj = 1, k_dim
!                     DO l = 1, mesh%gauss%l_G
!                        bmat(ni,nj,kj) = bmat(ni,nj,kj) + ww_c(ni,l)*mesh%gauss%dw(kj,nj,l,m)*mesh%gauss%rj(l,m)
!                     END DO
!                  END DO
!               END DO
!            END DO
!            DO ni = 1, pp_mesh%gauss%n_w
!               DO nj = 1, mesh%gauss%n_w
!                  DO kj = 1, k_dim
!                      bpmat(ni,nj,kj) = SUM(loc_mass_inv(ni,:)*bmat(:,nj,kj))
!                  END DO
!               END DO
!            END DO
!            bpmat = bpmat*(0.5d0/SUM(pp_mesh%gauss%rj(:, m)))
!         END IF
!
!         jj_loc = mesh%jj(:, m)
!         mat_loc = 0.d0
!         DO ni = 1, n_w
!            i = jj_loc(ni)
!            DO ki = 1, k_max
!               iglob = LA%loc_to_glob(ki, i)
!               ix = (ki - 1) * n_w + ni
!               idxm(ix) = iglob - 1
!               DO nj = 1, n_w
!                  j = jj_loc(nj)
!                  DO kj = 1, k_max
!                     jglob = LA%loc_to_glob(kj, j)
!                     jx = (kj - 1) * n_w + nj
!                     idxn(jx) = jglob - 1
!                     x = lambda_dg*SUM(bmat(:,ni,ki)*bpmat(:,nj,kj))
!                     DO l = 1, mesh%gauss%l_G
!                        y =  visc*mesh%gauss%dw(kj,ni,l,m)*mesh%gauss%dw(ki,nj,l,m) &
!                             + lambda_grad_div*mesh%gauss%dw(ki,ni,l,m)*mesh%gauss%dw(kj,nj,l,m)
!                        IF (kj.EQ.ki) THEN
!                           y = y + mass*mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
!                           DO k1 = 1, k_dim
!                              y = y + visc*mesh%gauss%dw(k1,ni,l,m)*mesh%gauss%dw(k1,nj,l,m)
!                           END DO
!                        END IF
!                        x = x + y * mesh%gauss%rj(l,m)
!                     END DO
!                     mat_loc(ix,jx) = x
!                  END DO
!               END DO
!            END DO
!         END DO
!
!         CALL MatSetValues(matrix, k_max * n_w, idxm, k_max * n_w, idxn, mat_loc, ADD_VALUES, ierr)
!      ENDDO
!
!      CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
!      CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
!
!      DEALLOCATE(idxm, idxn)
!
!    END SUBROUTINE qs_diff_mass_vect_M


    
    SUBROUTINE qs_LAP_mass_vect_M (LA, mesh, visc, mass, matrix)
      !=================================================
      USE my_util
      USE basis_change
      USE st_matrix
      IMPLICIT NONE
      TYPE(mesh_type), INTENT(IN) :: mesh
      REAL(KIND = 8), INTENT(IN) :: visc, mass
      TYPE(petsc_csr_LA) :: LA
      INTEGER, DIMENSION(mesh%gauss%n_w) :: jj_loc
      INTEGER, DIMENSION(:), ALLOCATABLE :: idxm, idxn
      REAL(KIND = 8), DIMENSION(:, :), POINTER, SAVE :: aij_p1p2, aij_p2p3
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: aij
      REAL(KIND = 8), DIMENSION(2, 2) :: dwdw

      LOGICAL :: once_p1p2 = .TRUE., once_p2p3 = .TRUE.
      REAL(KIND = 8) :: hatK_over_K
      
      REAL(KIND = 8) :: rj, rjrj, mrj, vrj, vvtrj, vtrj
      INTEGER :: m, l, ni, nj, i, j, iglob, jglob, k_max, n_w, ix, jx, ki, kj, nl, nk
      REAL(KIND = 8), DIMENSION(2*mesh%gauss%n_w,2*mesh%gauss%n_w) :: mat_loc
      INTEGER :: k_dim, n_b, k1
      REAL(KIND=8) :: x, y, lambda_grad_div, lambda_dg
      !#include "petsc/finclude/petsc.h"
      Mat                                         :: matrix
      PetscErrorCode                              :: ierr
      CALL MatZeroEntries (matrix, ierr)
      CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

      k_max = 2 ! 2x2 Structure

      IF (k_max/=SIZE(LA%loc_to_glob, 1)) THEN
         CALL error_Petsc('BUG in qs_diff_mass_vect_petsc_M, k_max/=SIZE(LA%loc_to_glob,1)')
      END IF

      n_w = mesh%gauss%n_w

      ALLOCATE(idxm(k_max * n_w), idxn(k_max * n_w))

      IF (inputs%type_fe==2) THEN
         IF (once_p1p2) THEN
            CALL p1_p2(aij_p1p2)
            once_p1p2 = .FALSE.
         END IF
         aij => aij_p1p2
      ELSE IF (inputs%type_fe==3) THEN
         IF (once_p2p3) THEN
            CALL p2_p3(aij_p2p3)
            once_p2p3 = .FALSE.
         END IF
         aij => aij_p2p3
      ELSE
         WRITE(*, *) 'BUG in Bt_Mminus_B, inputs&type_fe_u not correct'
         STOP
      END IF

      k_dim = 2
      n_b = k_dim
      DO m = 1, mesh%dom_me
         jj_loc = mesh%jj(:, m)
         mat_loc = 0.d0
         DO ni = 1, n_w
            i = jj_loc(ni)
            DO ki = 1, k_max
               iglob = LA%loc_to_glob(ki, i)
               ix = (ki - 1) * n_w + ni
               idxm(ix) = iglob - 1
               DO nj = 1, n_w
                  j = jj_loc(nj)
                  DO kj = 1, k_max
                     jglob = LA%loc_to_glob(kj, j)
                     jx = (kj - 1) * n_w + nj
                     idxn(jx) = jglob - 1
                     x = 0.d0
                     DO l = 1, mesh%gauss%l_G
                        y=0.d0
                        IF (kj.EQ.ki) THEN
                           y = mass*mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
                           DO k1 = 1, k_dim
                              y = y + visc*mesh%gauss%dw(k1,ni,l,m)*mesh%gauss%dw(k1,nj,l,m)
                           END DO
                        END IF
                        x = x + y * mesh%gauss%rj(l,m)
                     END DO
                     mat_loc(ix,jx) = x
                  END DO
               END DO
            END DO
         END DO

         CALL MatSetValues(matrix, k_max * n_w, idxm, k_max * n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO

      CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)

      DEALLOCATE(idxm, idxn)

    END SUBROUTINE qs_LAP_mass_vect_M

!   SUBROUTINE qs_BBt(LA, matrix, bc)
!      USE my_util
!      USE basis_change
!      USE st_matrix
!      IMPLICIT NONE
!      TYPE(petsc_csr_LA) :: LA
!      INTEGER, DIMENSION(:), POINTER :: bc
!      INTEGER, DIMENSION(:), ALLOCATABLE :: neighs
!      INTEGER, DIMENSION(pp_mesh%gauss%n_w) :: idxm, idxn
!      REAL(KIND = 8), DIMENSION(pp_mesh%gauss%n_w, uu_mesh%gauss%l_G) :: ww_c
!      REAL(KIND = 8), DIMENSION(uu_mesh%gauss%n_w, uu_mesh%gauss%l_G, pp_mesh%dom_me) :: div
!      REAL(KIND = 8), DIMENSION(:, :), POINTER, SAVE :: aij_p1p2, aij_p2p3
!      REAL(KIND = 8), DIMENSION(:, :), POINTER :: aij
!
!      LOGICAL :: once_p1p2 = .TRUE., once_p2p3 = .TRUE.
!      REAL(KIND = 8), DIMENSION(pp_mesh%gauss%n_w, pp_mesh%gauss%n_w) :: mat_loc
!      INTEGER :: m, l, ni, nk, nl, nw, ix, jx, mn, nb_neighs, m_neigh
!      REAL(KIND = 8) :: viscolm, xij
!      !#include "petsc/finclude/petsc.h
!
!      Mat                                         :: matrix
!      PetscErrorCode                              :: ierr
!      CALL MatZeroEntries (matrix, ierr)
!      CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
!
!      IF (inputs%type_fe==2) THEN
!         IF (once_p1p2) THEN
!            CALL p1_p2(aij_p1p2)
!            once_p1p2 = .FALSE.
!         END IF
!         aij => aij_p1p2
!      ELSE IF (inputs%type_fe==3) THEN
!         IF (once_p2p3) THEN
!            CALL p2_p3(aij_p2p3)
!            once_p2p3 = .FALSE.
!         END IF
!         aij => aij_p2p3
!      ELSE
!         WRITE(*, *) 'BUG in Bt_Mminus_B, inputs&type_fe_u not correct'
!         STOP
!      END IF
!
!      !===Construct w_c
!      DO l = 1, uu_mesh%gauss%l_G
!         DO ni = 1, pp_mesh%gauss%n_w
!            ww_c(ni, l) = SUM(aij(ni, :) * uu_mesh%gauss%ww(:, l))
!         END DO
!      END DO
!
!      div = uu_mesh%gauss%dw(1, :, :, :) + uu_mesh%gauss%dw(2, :, :, :)
!      nw = pp_mesh%gauss%n_w
!      ALLOCATE(neighs(0))
!      DO m = 1, pp_mesh%dom_me
!         idxn = LA%loc_to_glob(1, pp_mesh%jj(:, m)) - 1
!         DO ni = 1, uu_mesh%gauss%n_w
!            CALL neighbours(uu_mesh, m, ni, neighs, nb_neighs)
!            write(*, *) nb_neighs, neighs
!            DO nk = 1, pp_mesh%gauss%n_w
!               DO mn = 1, nb_neighs
!                  m_neigh = neighs(mn)
!                  idxm = LA%loc_to_glob(1, pp_mesh%jj(:, m_neigh)) - 1
!                  mat_loc = 0.d0
!                  DO nl = 1, pp_mesh%gauss%n_w
!                     mat_loc(nk, nl) = mat_loc(nk, nl) + SUM(ww_c(nk, :) * div(ni, :, m) * uu_mesh%gauss%rj(:, m)&
!                          * ww_c(nl, :) * div(ni, :, m_neigh) * uu_mesh%gauss%rj(:, m_neigh))
!                  END DO
!                  CALL MatSetValues(matrix, nw, idxn, nw, idxm, mat_loc, ADD_VALUES, ierr)
!               END DO
!            END DO
!         END DO
!      ENDDO
!
!      CALL MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
!      CALL MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
!
!   END SUBROUTINE qs_BBt

   SUBROUTINE inverse_mat(a, n)       ! Invert matrix by Gauss method
      ! --------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: n
      REAL(8) :: a(n, n)

      ! - - - Local Variables - - -
      REAL(8) :: b(n, n), c, d, temp(n)
      INTEGER :: i, j, k, m, imax(1), ipvt(n)
      ! - - - - - - - - - - - - - -

      b = a
      ipvt = (/ (i, i = 1, n) /)

      DO k = 1, n
         imax = MAXLOC(ABS(b(k:n, k)))
         m = k - 1 + imax(1)

         IF (m /= k) THEN
            ipvt((/m, k/)) = ipvt((/k, m/))
            b((/m, k/), :) = b((/k, m/), :)
         END IF
         d = 1 / b(k, k)

         temp = b(:, k)
         DO j = 1, n
            c = b(k, j) * d
            b(:, j) = b(:, j) - temp * c
            b(k, j) = c
         END DO
         b(:, k) = temp * (-d)
         b(k, k) = d
      END DO

      a(:, ipvt) = b

   END SUBROUTINE inverse_mat
END MODULE fem_M
