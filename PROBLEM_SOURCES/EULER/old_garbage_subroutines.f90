
  !GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE
  !GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE
  !GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE GARBADGE
  SUBROUTINE compute_dk (this, un)
    USE arbitrary_eos_lambda_module
    USE my_util, ONLY : error_petsc
    IMPLICIT NONE
    CLASS(euler_type) :: this
    REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
    INTEGER, DIMENSION(1) :: i_t, j_t
    REAL(KIND = 8), DIMENSION(1, this%mesh%gauss%k_d) :: nij_c
    REAL(KIND = 8), DIMENSION(1) :: norm_c, dijL_c
    REAL(KIND = 8), DIMENSION(2) :: u, rho, ie, p, lambda_max
    LOGICAL, DIMENSION(this%mesh%medge) :: virgin_edge
    REAL(KIND = 8) :: pstar
    LOGICAL :: bug
    INTEGER :: m, ni, nj, nw, n, i, j, k, ierr, edge, divider, nb_shared_cell
    nw = this%mesh%gauss%n_w

    bug = .FALSE.
    SELECT CASE(this%mesh%gauss%k_d)
    CASE(1)
       nb_shared_cell = 1
       IF (this%mesh%gauss%n_w/=2) bug=.TRUE.
    CASE(2)
       nb_shared_cell = 2
       IF (this%mesh%gauss%n_w/=3) bug=.TRUE.
    END SELECT
    IF (bug) THEN
       CALL error_petsc('Wrong polynomial degree for low-order viscosity')
    END IF

    DO m = 1, this%mesh%dom_me
       DO n = 1, this%mesh%gauss%n_e
          IF (this%mesh%attr_e(this%mesh%jce(n, m))) THEN
             edge = this%mesh%jce_loc(n, m)
             IF (.NOT. virgin_edge(edge)) CYCLE
             virgin_edge(edge) = .FALSE.
             ni = MOD(n, nw) + 1
             nj = MOD(n + 1, nw) + 1
             i = this%mesh%jj(ni, m)
             j = this%mesh%jj(nj, m)
             i_t = i
             j_t = j
             DO k = 1, this%mesh%gauss%k_d
                CALL MatGetValues(this%matrices%nij_loc(k), 1, i_t - 1, 1, j_t - 1, &
                     nij_c(:, k), ierr)
             END DO
             rho(1) = un(i, 1)
             rho(2) = un(j, 1)
             u(1) = SUM(un(i, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(1)
             u(2) = SUM(un(j, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(2)
             ie(1) = un(i, this%mesh%gauss%k_d + 2) / rho(1) - 0.5d0 * u(1) * u(1)
             ie(2) = un(j, this%mesh%gauss%k_d + 2) / rho(2) - 0.5d0 * u(2) * u(2)
             p = this%pressure(rho, ie)
             CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, &
                  lambda_max, pstar)
             dijL_c = MAXVAL(lambda_max) * norm_c
             divider = nb_shared_cell

             IF (this%mesh%side_edge(n, m)) THEN !=== if on the boundary, switch i for j
                DO k = 1, this%mesh%gauss%k_d
                   CALL MatGetValues(this%matrices%nij_loc(k), 1, j_t - 1, 1, i_t - 1, &
                        nij_c(:, k), ierr)
                END DO
                u(1) = SUM(un(i, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(1)
                u(2) = SUM(un(j, 2:1 + this%mesh%gauss%k_d) * nij_c(1, :)) / rho(2)
                rho = (/rho(2), rho(1)/)
                ie = (/ie(2), ie(1)/)
                p = (/p(2), p(1)/)
                CALL lambda_arbitrary_eos(this%eos_param, rho, u, ie, p, this%in_tol, this%no_iter, &
                     lambda_max, pstar)
                dijL_c = MAX(dijL_c, MAXVAL(lambda_max) * norm_c)
                divider = 1
             END IF

             this%matrices%dK(m) = MAX(this%matrices%dK(m),dijL_c(1)/divider)
          END IF
       END DO
    END DO
  END SUBROUTINE compute_dk

   SUBROUTINE compute_dt_from_dK(this)
     IMPLICIT NONE
     CLASS(euler_type) :: this
     REAL(KIND = 8), DIMENSION(this%mesh%dom_np) :: dijL_diag
     REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w) :: v_loc
     INTEGER, DIMENSION(this%mesh%gauss%n_w) :: idxm
     INTEGER :: i, m, ni, iglob
     REAL(KIND = 8) :: dt_min_loc, dt_min_glob
     Vec                                         :: vect
     PetscErrorCode                              :: ierr
     CALL VecSet(vect, 0.d0, ierr)

     WRITE(*,*) "VB: WARNING (20/04/2026) ==> this subroutine does not call any ghost points???"
     STOP

     DO m = 1, this%mesh%me
        v_loc = 0.d0
        DO ni = 1, this%mesh%gauss%n_w
           i = this%mesh%jj(ni, m)
           iglob = this%LA%loc_to_glob(1, i)
           idxm(ni) = iglob - 1
           v_loc(ni) = v_loc(ni) + this%matrices%dK(m)
        ENDDO
        CALL VecSetValues(vect, this%mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
     ENDDO
     CALL VecAssemblyBegin(vect, ierr)
     CALL VecAssemblyEnd(vect, ierr)

     CALL VecGetValues(this%vec_loc, this%mesh%dom_np, this%tab, dijL_diag, ierr)

     WRITE(*,*) "VB: WARNING (01/05/2026) ==> this subroutine (not used right now) uses lumped_mass.",&
     " Must be rewritten using lump_mass_vec instead"
     STOP
   !   dijL_diag = this%matrices%lumped_mass(1:this%mesh%dom_np) / ABS(dijL_diag)

     dt_min_loc = MINVAL(dijL_diag) / 2.d0

     CALL MPI_ALLREDUCE(dt_min_loc, dt_min_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, PETSC_COMM_WORLD, ierr)
     this%dt = this%CFL * dt_min_glob
   END SUBROUTINE compute_dt_from_dK

   SUBROUTINE compute_flux(this, ff, Vect)
      USE space_dim
      IMPLICIT NONE
      CLASS(euler_type) :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim) :: ff
      REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w) :: v_loc
      REAL(KIND = 8), DIMENSION(this%mesh%gauss%n_w, k_dim) :: f_loc
      REAL(KIND = 8), DIMENSION(this%mesh%np) :: v_glb
      INTEGER, DIMENSION(this%mesh%gauss%n_w) :: idxm, jj_loc
      REAL(KIND = 8) :: x
      INTEGER :: k, m, ni, nj
      Vec                                         :: vect
      PetscErrorCode                              :: ierr
      CALL VecSet(vect, 0.d0, ierr)
      v_glb = 0.d0
      DO m = 1, this%mesh%dom_me
         jj_loc = this%mesh%jj(:, m)
         f_loc = ff(jj_loc,:)
         !<==recompute cij on the fly
         DO ni = 1, this%mesh%gauss%n_w
            !wwrj = this%mesh%gauss%ww(ni,:)*this%mesh%gauss%rj(:,m)
            x = 0.d0
            DO k = 1, this%mesh%gauss%k_d
               DO nj = 1, this%mesh%gauss%n_w
                  x = x + f_loc(nj,k)* &
                        !SUM(this%mesh%gauss%dw(k,nj,:,m)*wwrj)
                     SUM(this%mesh%gauss%dw(k,nj,:,m)*this%mesh%gauss%ww(ni,:)*this%mesh%gauss%rj(:,m))
               ENDDO
            ENDDO
            v_loc(ni) = x
         ENDDO
         idxm = this%LA%loc_to_glob(1, jj_loc) -1
         v_loc = -v_loc
         CALL VecSetValues(vect, this%mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
   !!$        v_glb(jj_loc) = v_glb(jj_loc) - v_loc
      ENDDO
   !!$     CALL VecSetValues(vect, this%mesh%np, this%LA%loc_to_glob(1,:)-1, v_glb, INSERT_VALUES, ierr)
      CALL VecAssemblyBegin(vect, ierr)
      CALL VecAssemblyEnd(vect, ierr)
   END SUBROUTINE compute_flux
