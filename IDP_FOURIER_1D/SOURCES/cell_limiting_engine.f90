MODULE cell_limiting_engine_module
  PUBLIC :: iterative_cell_limiting_procedure, iterative_GS_cell_limiting_procedure
  TYPE mass_for_limiting
     REAL(KIND=8), DIMENSION(:,:), POINTER :: localized_mass
     REAL(KIND=8), DIMENSION(:),   POINTER :: lumped_mass
     REAL(KIND=8) :: mass_eps
  END TYPE mass_for_limiting
  REAL(KIND=8), PARAMETER, PRIVATE :: epsilon=1.d-14
CONTAINS

  SUBROUTINE iterative_cell_limiting_procedure(masses,jj,xx_in,loc_min,psi,zero_of_psi,xx_out)  
    IMPLICIT NONE
    INTERFACE
       FUNCTION psi(x,psi_min) RESULT(v)
         REAL(KIND=8), DIMENSION(:) :: x
         REAL(KIND=8) :: psi_min, v
       END FUNCTION psi
       FUNCTION zero_of_psi(psi_min,u0,P) RESULT(v)
         REAL(KIND=8), DIMENSION(:) :: u0, P
         REAL(KIND=8) :: psi_min, v
       END FUNCTION zero_of_psi
    END INTERFACE
    TYPE(mass_for_limiting)      :: masses
    REAL(KIND=8), DIMENSION(:,:) :: xx_in, xx_out
    REAL(KIND=8), DIMENSION(:)   :: loc_min
    INTEGER,      DIMENSION(:,:) :: jj
    REAL(KIND=8), DIMENSION(SIZE(xx_in,2))               :: uk_minus, uk_plus
    REAL(KIND=8), DIMENSION(SIZE(jj,1),SIZE(jj,2),SIZE(xx_in,2))    :: xx
    REAL(KIND=8), DIMENSION(SIZE(jj,1),SIZE(xx_in,2))    :: xx_loc
    REAL(KIND=8), DIMENSION(SIZE(xx_in,1),SIZE(xx_in,2)) :: xx_inter
    REAL(KIND=8), DIMENSION(SIZE(jj,1)) :: lambda_minus, lambda_plus
    REAL(KIND=8), DIMENSION(SIZE(jj,1)) :: loc_min_loc
    INTEGER,      DIMENSION(SIZE(jj,1)) :: jloc
    INTEGER,      DIMENSION(SIZE(jj,1)) :: limit_zero, limit_plus, limit_minus
    INTEGER :: i, k, m, n, me, nw, syst_size, iminus, iplus
    REAL(KIND=8) :: loc_min_down, loc_min_up
    REAL(KIND=8) :: mass_plus, mass_minus, &
         lambda_K_minus, lambda_K_plus, &
         lambda_star_minus, lambda_star_plus

    me = SIZE(jj,2)
    nw = SIZE(jj,1)
    syst_size = SIZE(xx_in,2)
    DO m = 1, me
       lambda_minus = 1.d0
       lambda_plus = 1.d0
       limit_zero = 0
       limit_minus = 0
       limit_plus = 0
       uk_minus = 0.d0
       uk_plus  = 0.d0
       mass_plus = 0.d0
       mass_minus = 0.d0
       jloc = jj(:,m)
       DO k = 1, syst_size
          xx_loc(:,k) = xx_in(jloc,k)
       END DO
       loc_min_loc = loc_min(jloc)
       iminus = 0
       iplus  = 0
       DO n = 1, nw
          !===P2 fix
          IF (ABS(masses%lumped_mass(jloc(n))).LE.masses%mass_eps) THEN
             limit_zero(n) = 1
             CYCLE
          END IF
          !===END fix

          loc_min_down = loc_min_loc(n) - epsilon*ABS(loc_min_loc(n))
          loc_min_up   = loc_min_loc(n) + epsilon*ABS(loc_min_loc(n))
          IF (psi(xx_loc(n,:),loc_min_down)<0) THEN
             iplus = iplus + 1
             uk_minus = uk_minus + masses%localized_mass(n,m)*xx_loc(n,:)
             mass_minus = mass_minus + masses%localized_mass(n,m)
             limit_minus(n) = 1
          ELSE IF (psi(xx_loc(n,:),loc_min_up)>0) THEN   
             iminus = iminus + 1
             uk_plus = uk_plus + masses%localized_mass(n,m)*xx_loc(n,:)
             mass_plus = mass_plus + masses%localized_mass(n,m)
             limit_plus(n) = 1
          ELSE
             limit_zero(n) = 1
          END IF
       END DO
       IF (SUM(limit_zero+limit_plus+limit_minus).NE.nw) THEN
          WRITE(*,*) ' BUG ', limit_zero+limit_plus+limit_minus
          STOP
       END IF
       IF (iplus*iminus==0) THEN
          xx(:,m,:) = xx_loc
          CYCLE !===No limiting is possible/or no limiting necessary
       END IF
       uk_minus = uk_minus/mass_minus
       uk_plus  = uk_plus/mass_plus
       DO n = 1, nw
          !===Lambda_minus
          IF (limit_minus(n)==1) THEN
             lambda_minus(n) = zero_of_psi(loc_min_loc(n),uk_plus,xx_loc(n,:)-uk_plus)
          END IF
          !===Lambda_plus
          IF (limit_plus(n)==1) THEN
             lambda_plus(n) = zero_of_psi(loc_min_loc(n),xx_loc(n,:),uk_minus-xx_loc(n,:))
          END IF
       END DO
       lambda_minus = MAX(MIN(lambda_minus,1.d0),0.d0)
       lambda_plus  = MAX(MIN(lambda_plus,1.d0),0.d0)
       Lambda_star_minus = MINVAL(lambda_minus)
       Lambda_star_plus  = MINVAL(lambda_plus)
       Lambda_K_minus = MAX(Lambda_star_minus, 1.d0-Lambda_star_plus*mass_plus/mass_minus)
       Lambda_K_plus  = MIN(Lambda_star_plus, (1.d0-Lambda_star_minus)*mass_minus/mass_plus)

       DO n = 1, nw
          !===P2 fix
          IF (ABS(masses%lumped_mass(jloc(n))).LE.masses%mass_eps) THEN
             xx(n,m,:) = uk_plus
             CYCLE
          END IF
          !===END fix
          xx(n,m,:) = xx_loc(n,:) &
               +limit_minus(n)*(1-Lambda_K_minus)*(uk_plus-xx_loc(n,:))&
               +limit_plus(n) *     Lambda_K_plus*(uK_minus-xx_loc(n,:))
       END DO
    END DO
    xx_inter = 0.d0
    DO m = 1, me
       DO n = 1, nw
          i = jj(n,m)
          IF (ABS(masses%lumped_mass(i)).LE.masses%mass_eps) THEN
             xx_out(i,:) = xx(n,m,:)
             CYCLE
          END IF
          xx_inter(i,:) = xx_inter(i,:)+xx(n,m,:)*masses%localized_mass(n,m)
       END DO
    END DO
    DO i = 1, SIZE(xx_out,1)
       IF (masses%lumped_mass(i).GT.masses%mass_eps) THEN
          xx_out(i,:)= xx_inter(i,:)/masses%lumped_mass(i)
       END IF
    END DO
  END SUBROUTINE iterative_cell_limiting_procedure

  SUBROUTINE iterative_GS_cell_limiting_procedure(masses,jj,xx_in,loc_min,psi,zero_of_psi,xx_out)  
    IMPLICIT NONE
    INTERFACE
       FUNCTION psi(x,psi_min) RESULT(v)
         REAL(KIND=8), DIMENSION(:) :: x
         REAL(KIND=8) :: psi_min, v
       END FUNCTION psi
       FUNCTION zero_of_psi(psi_min,u0,P) RESULT(v)
         REAL(KIND=8), DIMENSION(:) :: u0, P
         REAL(KIND=8) :: psi_min, v
       END FUNCTION zero_of_psi
    END INTERFACE
    TYPE(mass_for_limiting)      :: masses
    REAL(KIND=8), DIMENSION(:,:) :: xx_in, xx_out
    REAL(KIND=8), DIMENSION(:)   :: loc_min
    INTEGER,      DIMENSION(:,:) :: jj
    REAL(KIND=8), DIMENSION(SIZE(xx_in,2))            :: uk_minus, uk_plus, err
    REAL(KIND=8), DIMENSION(SIZE(jj,1),SIZE(jj,2))    :: xx
    REAL(KIND=8), DIMENSION(SIZE(jj,1),SIZE(xx_in,2)) :: xx_loc
    REAL(KIND=8), DIMENSION(SIZE(jj,1)) :: loc_min_loc
    REAL(KIND=8), DIMENSION(SIZE(jj,1)) :: lambda_minus, lambda_plus
    INTEGER,      DIMENSION(SIZE(jj,1)) :: jloc
    INTEGER,      DIMENSION(SIZE(jj,1)) :: limit_zero, limit_plus, limit_minus
    INTEGER :: i, k, m, n, me, nw, syst_size, iminus, iplus
    REAL(KIND=8) :: loc_min_down, loc_min_up
    REAL(KIND=8) :: mass_plus, mass_minus, &
         lambda_K_minus, lambda_K_plus, &
         lambda_star_minus, lambda_star_plus
    REAL(KIND=8) :: xx_inter(SIZE(xx_in))
    me = SIZE(jj,2)
    nw = SIZE(jj,1)
    syst_size = SIZE(xx_in,2)
    DO m = 1, me
       lambda_minus = 1.d0
       lambda_plus = 1.d0
       limit_zero = 0
       limit_minus = 0
       limit_plus = 0
       uk_minus = 0.d0
       uk_plus  = 0.d0
       mass_plus = 0.d0
       mass_minus = 0.d0
       jloc = jj(:,m)
       DO k = 1, syst_size
          xx_loc(:,k) = xx_in(jloc,k)
       END DO
       loc_min_loc = loc_min(jloc)
       iminus = 0
       iplus  = 0
       DO n = 1, nw
          loc_min_down = loc_min_loc(n) - epsilon*ABS(loc_min_loc(n))
          loc_min_up   = loc_min_loc(n) + epsilon*ABS(loc_min_loc(n))
          IF (psi(xx_loc(n,:),loc_min_down)<0) THEN
             iplus = iplus + 1
             uk_minus = uk_minus + masses%lumped_mass(jloc(n))*xx_loc(n,:)
             mass_minus = mass_minus + masses%lumped_mass(jloc(n))
             limit_minus(n) = 1
          ELSE IF (psi(xx_loc(n,:),loc_min_up)>0) THEN   
             iminus = iminus + 1
             uk_plus = uk_plus + masses%lumped_mass(jloc(n))*xx_loc(n,:)
             mass_plus = mass_plus + masses%lumped_mass(jloc(n))
             limit_plus(n) = 1
          ELSE
             limit_zero(n) = 1
          END IF
       END DO
       IF (SUM(limit_zero+limit_plus+limit_minus).NE.nw) THEN
          WRITE(*,*) ' BUG ', limit_zero+limit_plus+limit_minus
          STOP
       END IF
       IF (iplus*iminus==0) THEN
          CYCLE !===No limiting is possible
       END IF
       uk_minus = uk_minus/mass_minus
       uk_plus  = uk_plus/mass_plus
       DO n = 1, nw
          !===Lambda_minus
          IF (limit_minus(n)==1) THEN
             lambda_minus(n) = zero_of_psi(loc_min_loc(n),uk_plus,xx_loc(n,:)-uk_plus)
          END IF
          !===Lambda_plus
          IF (limit_plus(n)==1) THEN
             lambda_plus(n) = zero_of_psi(loc_min_loc(n),xx_loc(n,:),uk_minus-xx_loc(n,:))
          END IF
       END DO
       lambda_minus = MAX(MIN(lambda_minus,1.d0),0.d0)
       lambda_plus  = MAX(MIN(lambda_plus,1.d0),0.d0)
       Lambda_star_minus = MINVAL(lambda_minus)
       Lambda_star_plus  = MINVAL(lambda_plus)
       Lambda_K_minus = MAX(Lambda_star_minus, 1.d0-Lambda_star_plus*mass_plus/mass_minus)
       Lambda_K_plus  = MIN(Lambda_star_plus, (1.d0-Lambda_star_minus)*mass_minus/mass_plus)
       !err =0.d0
       DO n = 1, nw
          i = jj(n,m)
          !===P2 fix
          IF (ABS(masses%lumped_mass(jloc(n))).LE.masses%mass_eps) THEN
             xx_in(i,:) = uk_plus
             CYCLE
          END IF
          !===END fix
          xx_in(i,:) = xx_loc(n,:) &
               +limit_minus(n)*(1-Lambda_K_minus)*(uk_plus-xx_loc(n,:))&
               +limit_plus(n) *     Lambda_K_plus*(uK_minus-xx_loc(n,:))
          !err = err + masses%lumped_mass(jloc(n))*(xx_in(i,:) - xx_loc(n,:))
       END DO
       !write(*,*) m, ' TESTT err', SUM(ABS(err))/sum(masses%lumped_mass(jloc(:)))
    END DO
    xx_out=xx_in
  END SUBROUTINE iterative_GS_cell_limiting_procedure

  SUBROUTINE relax_min_and_max(bound_relaxing,glob_min,glob_max,jj,un,maxn,minn)
    IMPLICIT NONE
    CHARACTER(*),               INTENT(IN) :: bound_relaxing
    INTEGER, DIMENSION(:,:),    INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(:)             :: minn
    REAL(KIND=8), DIMENSION(:)             :: maxn
    REAL(KIND=8), INTENT(IN)               :: glob_min, glob_max
    REAL(KIND=8), DIMENSION(SIZE(un))      :: alpha, denom
    INTEGER, DIMENSION(SIZE(un)) ::   beta 
    INTEGER      :: i, j, m, me, nw, n, np
    REAL(KIND=8) :: norm
 
    me = SIZE(jj,2)
    nw = SIZE(jj,1)
    alpha = 0.d0
    beta = 0
    DO m = 1, me
       DO n = 1, nw
          i = jj(n,m)
          DO np = 1, nw
             j = jj(np,m)
             alpha(i) = alpha(i) + (un(i) - un(j))
             beta(i) = beta(i) + 1
          END DO
       END DO
    END DO
    alpha = alpha/beta
    SELECT CASE(TRIM(ADJUSTL(bound_relaxing)))
    CASE('avg') !==Average
       denom = 0.d0
       beta = 0
       DO m = 1, me
          DO n = 1, nw
             i = jj(n,m)
             DO np = 1, nw
                j = jj(np,m)
                denom(i) = denom(i) + alpha(i) + alpha(j)
                beta(i) = beta(i) + 1
             END DO
          END DO
       END DO
       denom = denom/(2*beta)
    CASE('minmod') !===Minmod
       denom = alpha 
       DO m = 1, me
          DO n = 1, nw
             i = jj(n,m)
             DO np = 1, nw
                j = jj(np,m)

                IF (denom(i)*alpha(j).LE.0.d0) THEN
                   denom(i) = 0.d0
                ELSE IF (ABS(denom(i)) > ABS(alpha(j))) THEN
                   denom(i) = alpha(j)
                END IF
             END DO
          END DO
       END DO
    CASE DEFAULT
       WRITE(*,*) ' BUG in relax', TRIM(ADJUSTL(bound_relaxing))
       STOP
    END SELECT
    maxn = maxn + 4.*ABS(denom)
    minn = minn - 4.*ABS(denom)
    maxn = MIN(glob_max,maxn)
    minn = MAX(glob_min,minn)
  END SUBROUTINE RELAX_MIN_AND_MAX
  

  
END MODULE cell_limiting_engine_module
