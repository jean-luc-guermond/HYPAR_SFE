MODULE IDP_update_euler
  USE matrix_type
  USE space_dim
  USE mesh_handling
  USE pardiso_solve
  USE Butcher_tableau
  USE input_data
  
  PUBLIC:: IDP_euler, IDP_construct_euler_matrices, IDP_compute_dt
  INTEGER,  PUBLIC :: isolve_euler_pardiso
  TYPE(BT), PUBLIC :: ERK
  PRIVATE
  LOGICAL, PUBLIC                              :: plot
  TYPE(matrice_bloc), DIMENSION(k_dim), PUBLIC :: cij, testcij
  TYPE(matrice_bloc)                           :: dijL
  TYPE(matrice_bloc), DIMENSION(:), POINTER    :: FluxijL, Fluxij_Gal
  TYPE(matrice_bloc)                           :: dijH, lij
  TYPE(matrice_bloc), DIMENSION(:,:), POINTER  :: FluxijH
  TYPE(matrice_bloc), PUBLIC                   :: mass, stiff
  TYPE(matrice_bloc)                           :: mc_minus_ml
  TYPE(matrice_bloc), DIMENSION(k_dim+2)       :: fctmat
  !TYPE(matrice_bloc)                           :: resij
  INTEGER, DIMENSION(:), POINTER, PUBLIC       :: diag
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC  :: lumped
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC  :: urelaxi, drelaxi
  REAL(KIND=8), PARAMETER                      :: small=1.d-8
  REAL(KIND=8), PARAMETER                      :: urelax=1.001, drelax=.999d0   
CONTAINS 

  SUBROUTINE IDP_construct_euler_matrices
    USE st_matrix
    USE fem_s_M
    USE lin_alg_tools
    USE CSR_transpose
        
    IMPLICIT NONE
    INTEGER :: i, k, p, l
    REAL(KIND=8) :: dd
    !===Mass
    CALL compute_mass(mesh,mass)
    
    !===lumped
    CALL lumped_mass(mesh,mass,lumped)

    !===diag
    CALL diag_mat(mass%ia,mass%ja,diag)

    !===Mass - lumped
    CALL duplicate(mass,mc_minus_ml)
    mc_minus_ml%aa = mass%aa
    mc_minus_ml%aa(diag) = mc_minus_ml%aa(diag) - lumped

    !===FluxijL (for limiting)
    ALLOCATE(fluxijL(inputs%syst_size))
    DO k = 1, inputs%syst_size
       CALL duplicate(mass,fluxijL(k))
    END DO
    
    !===fluxijH (for limiting)
    ALLOCATE(fluxijH(inputs%syst_size,ERK%s))
    DO l = 1, ERK%s
        DO k = 1, inputs%syst_size
           CALL duplicate(mass,fluxijH(k,l))
        END DO
    END DO

   !===stiff (for smoothness indicator)
    CALL compute_stiffness(mesh,stiff)
    
    !===dijL
    CALL duplicate(mass,dijL)
    
    !===dijH
    CALL duplicate(mass,dijH)
    
    !===fctmat
    DO k = 1, k_dim+2
       CALL st_csr(mesh%jj, fctmat(k)%ia, fctmat(k)%ja)
       ALLOCATE(fctmat(k)%aa(SIZE(fctmat(k)%ja)))
       fctmat(k)%aa = 0.d0
    END DO

    !===lij
    CALL duplicate(mass,lij)

    !===cij
    CALL compute_cij(mesh,cij)

    !===urelaxi, drelaxi
    ALLOCATE(urelaxi(mesh%np), drelaxi(mesh%np))
    DD = SUM(lumped)
    IF (k_dim==1) THEN
       urelaxi = MIN(1.d0 + 2*(lumped/DD)*SQRT(lumped/DD),urelax)
       drelaxi = MAX(1.d0 - 2*(lumped/DD)*SQRT(lumped/DD),drelax)
    ELSE IF (k_dim==2) THEN
       urelaxi = MIN(1.d0 + 2*SQRT(SQRT(lumped/DD))**3,urelax)
       drelaxi = MAX(1.d0 - 2*SQRT(SQRT(lumped/DD))**3,drelax)
    END IF
    
  END SUBROUTINE IDP_construct_euler_matrices

  SUBROUTINE divide_by_lumped(rk)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: rk
    INTEGER :: k
    DO k = 1, inputs%syst_size
       rk(:,k) = rk(:,k)/lumped
    END DO
  END SUBROUTINE divide_by_lumped

  SUBROUTINE full_step_ERK(un)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,ERK%s+1) :: un
    INTEGER  :: stage
    DO stage = 2, ERK%s+1
       !CALL one_step_ERK(stage,un)
    END DO
  END SUBROUTINE full_step_ERK

  SUBROUTINE one_step_ERK(stage,un)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size,ERK%s+1) :: un
    INTEGER, INTENT(IN)  :: stage
    REAL(KIND=8) :: time_stage
    
    CALL compute_fluxes(un(:,:,stage-1))
    
  END SUBROUTINE one_step_ERK

  SUBROUTINE compute_fluxes(un)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: rk

    CALL rhs(un,rk)
    CALL gal_flux(un,fluxij_Gal)
    IF (inputs%method_type=='Galerkin') THEN
    END IF
    
  END SUBROUTINE COMPUTE_FLUXES
  
  SUBROUTINE gal_flux(un,fluxij)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN)  :: un
    TYPE(matrice_bloc), DIMENSION(:), POINTER                      :: Fluxij
    REAL(KIND=8), DIMENSION(mesh%np,k_dim,inputs%syst_size)  :: vv
    REAL(KIND=8) :: rkloc
    INTEGER :: comp, i, j, p, k
    DO comp = 1, inputs%syst_size
       vv(:,:,comp)=flux(comp,un)
       fluxij(comp)%aa=0.d0
    END DO
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1)-1
          j = mass%ja(p)
          DO comp = 1, inputs%syst_size
             rkloc = 0.d0
             DO k = 1, k_dim
                rkloc  = rkloc - cij(k)%aa(p)*(vv(j,k,comp))
             END DO
             fluxij(comp)%aa(p) = rkloc
          END DO
       END DO
    END DO
  END SUBROUTINE gal_flux
  
  SUBROUTINE rhs(un,rk)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(mesh%np,k_dim,inputs%syst_size)  :: vv
    REAL(KIND=8) :: rkloc
    INTEGER :: comp, i, j, p, k
    DO comp = 1, inputs%syst_size
       vv(:,:,comp)=flux(comp,un)
    END DO
    rk=0.d0
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1)-1
          j = mass%ja(p)
          DO comp = 1, inputs%syst_size
             rkloc = 0.d0
             DO k = 1, k_dim
                rkloc  = rkloc - cij(k)%aa(p)*(vv(j,k,comp))
             END DO
             rk(i, comp) = rk(i, comp) + rkloc
          END DO
       END DO
    END DO
  END SUBROUTINE rhs

  SUBROUTINE galerkin(un,unext)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(OUT) :: unext
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: rk
    REAL(KIND=8), DIMENSION(mesh%np)                   :: ff
    INTEGER :: k
    CALL rhs(un,rk)
    IF (inputs%if_lumped) THEN
       CALL divide_by_lumped(rk)
       unext = un+inputs%dt*rk
    ELSE
       DO k = 1, inputs%syst_size
          CALL solve_pardiso(mass%aa,mass%ia,mass%ja,rk(:,k),ff,isolve_euler_pardiso)
          isolve_euler_pardiso=ABS(isolve_euler_pardiso)
          unext(:,k) = un(:,k)+inputs%dt*ff
       END DO
    END IF
  END SUBROUTINE galerkin

  SUBROUTINE viscous(un,unext)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(OUT) :: unext
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: rk
    INTEGER :: i, j, p, k
    CALL rhs(un,rk)
    CALL compute_dij(un)
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1)-1
          j = mass%ja(p)
          DO k = 1, inputs%syst_size
             rk(i, k) = rk(i, k) + dijL%aa(p)*(un(j,k)-un(i,k))
          END DO
       END DO
    END DO
    CALL divide_by_lumped(rk) !===We use lumped mass always
    unext = un+inputs%dt*rk
  END SUBROUTINE viscous

  SUBROUTINE high_order(un,ulow,unext)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(OUT) :: ulow, unext
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size)  :: rk, rkL, rkH
    REAL(KIND=8), DIMENSION(mesh%np)                   :: ff
    INTEGER :: i, j, k, p

    !===Galerkin RHS
    CALL rhs(un,rk)

    !===High-order viscosity
    CALL compute_dij(un)
    SELECT CASE(inputs%high_order_viscosity)
    CASE('EV(p)','EV(s)') 
       CALL entropy_residual(un,rk)
    CASE('SM(p)','SM(s)','SM(rho)')  
       CALL smoothness_viscosity(un)
    CASE DEFAULT
       WRITE(*,*) ' high_order_viscosity not defined'
       STOP
    END SELECT

    !===Add viscosity
    rkH=rk
    rkL=rk
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          rkH(i,:) = rkH(i,:) + dijH%aa(p)*(un(j,:)-un(i,:))
          rkL(i,:) = rkL(i,:) + dijL%aa(p)*(un(j,:)-un(i,:))
       END DO
    END DO
    
    !===Solve and update
    CALL divide_by_lumped(rkL)
    ulow = un+inputs%dt*rkL
    IF (inputs%if_lumped) THEN
       CALL divide_by_lumped(rkH)
       unext = un+inputs%dt*rkH
    ELSE
       DO k = 1, inputs%syst_size
          CALL solve_pardiso(mass%aa,mass%ia,mass%ja,rkH(:,k),ff,isolve_euler_pardiso)
          isolve_euler_pardiso=ABS(isolve_euler_pardiso)
          unext(:,k) = un(:,k)+inputs%dt*ff
       END DO
    END IF
  END SUBROUTINE high_order
  
  SUBROUTINE IDP_euler(un,unext)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2)  :: un
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2), INTENT(OUT) :: unext
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2)              :: ulow

    SELECT CASE(inputs%method_type) 
    CASE('galerkin') !===Galerkin
       CALL galerkin(un,unext)
       RETURN
    CASE('viscous') !===Viscous, firs-order
       CALL viscous(un,unext)
       RETURN
    CASE ('high') !===High-order
       CALL high_order(un,ulow,unext)
       IF (inputs%if_convex_limiting) THEN
          CALL convex_limiting_proc(un,ulow,unext)
       END IF
       RETURN
    CASE DEFAULT
       WRITE(*,*) ' BUG in euler: inputs%method_type not correct'
       STOP
    END SELECT
  END SUBROUTINE IDP_euler

  SUBROUTINE convex_limiting_proc(un,ulow,unext)
    USE mesh_handling
    USE boundary_conditions
    USE pardiso_solve
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2) :: un, ulow, unext, du
    REAL(KIND=8), DIMENSION(mesh%np,k_dim) :: vv
    REAL(KIND=8), DIMENSION(mesh%np) :: cc, ccmin, rhomax, rhomin
    INTEGER :: i, j, p, k, it
    REAL(KIND=8) :: rkloc, rhoij
    !===Entropy lower bound
    DO i = 1, mesh%np
       cc(i) = (un(i,k_dim+2)- SUM(un(i,2:k_dim+1)**2)/(2*un(i,1)))/(un(i,1)**gamma)
    END DO
    DO i = 1, mesh%np
       ccmin(i) = MINVAL(cc(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
    END DO

    !===Density bounds
    rhomax = un(:,1)
    rhomin = rhomin
    vv = flux(1,un) !===mass flux only
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1)-1
          j = mass%ja(p)
          rkloc = 0.d0
          DO k = 1, k_dim
             rkloc  = rkloc - cij(k)%aa(p)*(vv(j,k)-vv(i,k))
          END DO
          rhoij = (rkloc/dijL%aa(p) + un(j,1)+un(i,1))/2
          rhomin(i) = min(rhomin(i), rhoij)
          rhomax(i) = max(rhomax(i), rhoij)
       END DO
    END DO

    !===Relax bounds
    IF (inputs%if_relax_bounds) THEN
       !return
       CALL relax(un(:,1),rhomin,rhomax)
       CALL relax_cmin(un,ccmin)
    END IF

    !===Time increment if consistent matrix is used
    IF (inputs%if_lumped) THEN
       du = 0.d0
    ELSE
       du = unext-un
    END IF

    !===Limiting matrix, viscosity + mass matrix correction
    CALL compute_fct_matrix_full(un,du)

    !write(*,*) 'rhomin before limiting', minval(ulow(1,:))
    !===Convex limiting
    DO it = 1, 2 ! 2 and more is good for large CFL
       lij%aa = 1.d0
       !===Limit density
       !CALL FCT_generic(ulow,rhomax,rhomin,fctmat(1),lumped,1)
       CALL local_limit(ulow,rhomax,rhomin,fctmat(1),lumped,1) !===Works best

       !===Limit rho*e - e^(s_min) rho^gamma
       CALL limit_specific_entropy(ulow,un,ccmin) !===Works best

       !===Tranpose lij
       CALL transpose_op(lij,'min')
       CALL update_fct_full(ulow,ulow)
       DO k = 1, inputs%syst_size
          fctmat(k)%aa = (1-lij%aa)*fctmat(k)%aa
       END DO
    END DO
    unext = ulow
    !===End of computation
    RETURN

  END SUBROUTINE convex_limiting_proc

  SUBROUTINE relax(un,minn,maxn)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:)              :: un
    REAL(KIND=8), DIMENSION(:)              :: minn
    REAL(KIND=8), DIMENSION(:)              :: maxn
    REAL(KIND=8), DIMENSION(SIZE(un))       :: alpha, denom
    INTEGER      :: i, j, p
    REAL(KIND=8) :: norm
    alpha = 0.d0
    DO i = 1, mesh%np
       norm = 0.d0
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          IF (i==j) CYCLE
          alpha(i) = alpha(i) + stiff%aa(p)*(un(i) - un(j))
          !!!alpha(i) = alpha(i) + (un(i) - un(j))
          norm = norm + stiff%aa(p)
       END DO
       alpha(i) = alpha(i)/norm
    END DO
    SELECT CASE(inputs%limiter_type)
    CASE('avg') !==Average
       denom = 0.d0
       DO i = 1, SIZE(un)
          DO p = mass%ia(i), mass%ia(i+1) - 1
             j = mass%ja(p)
             IF (i==j) CYCLE
             denom(i) = denom(i) + alpha(j) + alpha(i)
          END DO
       END DO
       DO i = 1, SIZE(un)
          denom(i) = denom(i)/(2*(mass%ia(i+1)-mass%ia(i)-1))
       END DO
       maxn = MIN(urelaxi*maxn,maxn + ABS(denom)/2)
       minn = MAX(drelaxi*minn,minn - ABS(denom)/2)
    CASE ('minmod') !===Minmod
       denom = alpha    
       DO i = 1, SIZE(un)
          DO p = mass%ia(i), mass%ia(i+1) - 1
             j = mass%ja(p)
             IF (i==j) CYCLE
             IF (denom(i)*alpha(j).LE.0.d0) THEN
                denom(i) = 0.d0
             ELSE IF (ABS(denom(i)) > ABS(alpha(j))) THEN
                denom(i) = alpha(j)
             END IF
          END DO
       END DO
       maxn = MIN(urelaxi*maxn,maxn + ABS(denom)/2)
       minn = MAX(drelaxi*minn,minn - ABS(denom)/2)
    CASE DEFAULT
       WRITE(*,*) ' BUG in relax'
       STOP
    END SELECT
  END SUBROUTINE RELAX

  SUBROUTINE relax_cmin(un,cmin)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:)            :: un
    REAL(KIND=8), DIMENSION(:)              :: cmin
    REAL(KIND=8), DIMENSION(SIZE(cmin))     :: dc
    REAL(KIND=8), DIMENSION(k_dim+2)        :: ul
    INTEGER      :: i, j, p
    REAL(KIND=8) :: cl
    dc = 0.d0
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          IF (i==j) CYCLE
          ul(:) = (un(i,:)+un(j,:))/2
          cl = (ul(k_dim+2)-SUM(ul(2:k_dim+1)**2)/(2.d0*ul(1)))/(ul(1)**gamma)
          dc(i) = MAX(dc(i),cl-cmin(i))
       END DO
    END DO
    cmin = MAX(drelaxi*cmin, cmin - dc)
  END SUBROUTINE RELAX_cmin


  SUBROUTINE limit_specific_entropy(ulow,un,cmin)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2), INTENT(IN) :: ulow, un
    REAL(KIND=8), DIMENSION(mesh%np)                     :: cmin
    REAL(KIND=8), DIMENSION(mesh%np)  :: Esmall
    REAL(KIND=8), DIMENSION(k_dim+2)  :: ul, ur, Pij
    REAL(KIND=8) :: lambdai, coeff, psir, psil, ll, lr, llold, lrold, Budget
    INTEGER      :: i, j, p, k

    !===
    DO i = 1, mesh%np
       Esmall(i)= small*MINVAL(un(mass%ja(mass%ia(i):mass%ia(i+1)-1),k_dim+2))
       lambdai = 1.d0/(mass%ia(i+1) - 1.d0 - mass%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))

       !===Budget
!!$       Qplus = 0.d0
!!$       Card  = 0
!!$       DO p = mass%ia(i), mass%ia(i+1) - 1
!!$          j = mass%ja(p)
!!$          IF (i==j) CYCLE
!!$          ul = ulow(:,i)
!!$          DO k = 1 , k_dim+2
!!$             Pij(k) = fctmat(k)%aa(p)*coeff
!!$             ur(k) = ulow(k,i) + lij%aa(p)*Pij(k) !===Density must be positive
!!$          END DO
!!$          dQplus = MIN(psi_func(ul,cmin(i),0.d0),psi_func(ur,cmin(i),0.d0))
!!$          IF (dQplus>0.d0) THEN
!!$             Qplus = Qplus + dQplus
!!$          ELSE
!!$             Card  = Card + 1
!!$          END IF
!!$       END DO
!!$       IF (Card.NE.0) THEN
!!$          Budget = -Qplus/Card
!!$       ELSE
!!$          Budget = -1d15*Qplus
!!$       END IF
       Budget =0.d0
       !===End Budget

       DO p = mass%ia(i), mass%ia(i+1) - 1
          j =  mass%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p)
          DO k = 1, k_dim+2
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(i,k) + lr*Pij(k) !===Density must be positive
          END DO
          if (ur(1)<0.d0) then
             lij%aa(p) = 0.d0  !===CFL is too large
             CYCLE
          end if
          psir = psi_func(ur,cmin(i),Budget)
          IF (psir.GE.-Esmall(i)) THEN
             lij%aa(p) = lij%aa(p)
             CYCLE
          END IF
          ll = 0.d0
          ul = ulow(i,:)
          psil = psi_func(ul,cmin(i),Budget)
          DO WHILE (ABS(psil-psir) .GT. Esmall(i))
             llold = ll
             lrold = lr
             ll = ll - psil*(lr-ll)/(psir-psil)
             lr = lr - psir/psi_prime_func(Pij,ur,cmin(i))
             IF (ll.GE.lr) THEN
                ll = lr !lold
                EXIT
             END IF
             IF (ll< llold) THEN
                ll = llold
                EXIT
             END IF
             IF (lr > lrold) THEN
                lr = lrold
                EXIT
             END IF
             ul = ulow(i,:) + ll*Pij
             ur = ulow(i,:) + lr*Pij
             psil = psi_func(ul,cmin(i),Budget)
             psir = psi_func(ur,cmin(i),Budget)
          END DO
          IF (psir.GE.-Esmall(i)) THEN
             lij%aa(p) = lr
          ELSE
             lij%aa(p) = ll
          END IF
       END DO
    END DO

  CONTAINS
    FUNCTION psi_func(u,cmin,Budget) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u
      REAL(KIND=8),                     INTENT(IN) :: cmin
      REAL(KIND=8)                                 :: psi, Budget
      psi = u(k_dim+2) - SUM(u(2:k_dim+1)**2)/(2.d0*u(1)) - cmin*u(1)**gamma - Budget
    END FUNCTION psi_func
    FUNCTION psi_prime_func(Pij,u,cmin) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u, Pij
      REAL(KIND=8),                     INTENT(IN) :: cmin
      REAL(KIND=8)                                 :: psi
      psi = Pij(k_dim+2) - SUM(u(2:k_dim+1)*Pij(2:k_dim+1))/u(1) &
           + Pij(1)*SUM(u(2:k_dim+1)**2)/(2*u(1)**2) &
           - cmin*gamma*Pij(1)*u(1)**(gamma-1.d0)
    END FUNCTION psi_prime_func
  END SUBROUTINE limit_specific_entropy



  SUBROUTINE compute_fct_matrix_full(un,du)
    USE mesh_handling
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2), INTENT(IN) :: un, du
    INTEGER :: i, j, k, p
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          DO k = 1 , k_dim+2
             fctmat(k)%aa(p) = inputs%dt*(dijH%aa(p)-dijL%aa(p))*(un(j,k)-un(i,k)) &
                  -mc_minus_ml%aa(p)*(du(j,k)-du(i,k))
          END DO
       END DO
    END DO
  END SUBROUTINE compute_fct_matrix_full

  SUBROUTINE update_fct_full(ulow,unext)
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2) :: ulow
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2) :: unext
    REAL(KIND=8), DIMENSION(k_dim+2) :: x
    INTEGER :: i, k, p
    DO i = 1, mesh%np
       x = 0.d0
       DO p = mass%ia(i), mass%ia(i+1) - 1
          DO k = 1, k_dim+2
             x(k) = x(k) + lij%aa(p)*fctmat(k)%aa(p)
          END DO
       END DO
       DO k = 1, k_dim+2
          unext(i,k) = ulow(i,k) + x(k)/lumped(i)
       END DO
    END DO
  END SUBROUTINE update_fct_full

  SUBROUTINE compute_dij(un)
    USE mesh_handling
    USE boundary_conditions
    USE lambda_module
    USE lambda_module_full
    USE arbitrary_eos_lambda_module
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2), INTENT(IN)  :: un
    INTEGER                                       :: i, p, j, k
    REAL(KIND=8)                                  :: norm_cij, lambda_max, lambdal, lambdar, pl, pr
    REAL(KIND=8)                                  :: ml, mr, ul, ur, el, er, rhol, rhor, iel, ier
    REAL(KIND=8), DIMENSION(k_dim)                :: nij
    REAL(KIND=8) :: tol=1.d-3
    REAL(KIND=8) :: pstar
    LOGICAL      :: no_iter=.FALSE.
    INTEGER      :: nb_iterations
    REAL(KIND=8) :: uu(3,2), pp(2)
    DO i = 1, mesh%np
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = cij(1)%ja(p)
          IF (i.NE.j) THEN
             DO k = 1, k_dim
                nij(k) = cij(k)%aa(p)
             END DO
             norm_cij = SQRT(SUM(nij**2))
             nij=nij/norm_cij
             ml = SUM(un(i,2:k_dim+1)*nij)
             mr = SUM(un(j,2:k_dim+1)*nij)
             rhol = un(i,1)
             rhor = un(j,1)
             El = un(i,k_dim+2) - 0.5d0*(SUM(un(i,2:k_dim+1)**2) - ml**2)/rhol     
             Er = un(j,k_dim+2) - 0.5d0*(SUM(un(j,2:k_dim+1)**2) - mr**2)/rhor
             ul = ml/rhol
             ur = mr/rhor
             IF (inputs%equation_of_state=='gamma-law') THEN
                pr = ABS(Er-0.5d0*rhor*ur**2)*(gamma-1)
                pl = ABS(El-0.5d0*rhol*ul**2)*(gamma-1)
                CALL lambda(gamma,rhol,ul,pl,rhor,ur,pr,lambdal,lambdar)
             ELSE !===arbitrary eos
                uu(1,1) = rhol
                uu(2,1) = ml
                uu(3,1) = el
                uu(1,2) = rhor
                uu(2,2) = mr
                uu(3,2) = er
                pp = pressure(uu)
                pl = pp(1)
                pr = pp(2)
                iel = El/rhol-0.5d0*ul**2
                ier = Er/rhor-0.5d0*ur**2
                IF (min(iel,ier)<0.d0) THEN
                   WRITE(*,*) ' NEGATIVE iternal energy', min(iel,ier)
                   STOP
                END IF
                IF (min(pl,pr)<0.d0) THEN
                   WRITE(*,*) ' NEGATIVE pressure', min(pl,pr)
                   STOP
                END IF
                no_iter = .true.
                CALL lambda_arbitrary_eos(rhol,ul,iel,pl,rhor,ur,ier,pr,tol,no_iter, &
                     lambdal,lambdar,pstar,nb_iterations)
                !CALL lambda(1.02d0,rhol,ul,pl,rhor,ur,pr,lambdal,lambdar)
             END IF
             !CALL lambda_full(1.d-5,rhol,ul,pl,rhor,ur,pr,lambdal,lambdar,k)
             !CALL lambda_VDW(rhol,ul,el/rhol,rhor,ur,er/rhor,lambdal,lambdar)
             !write(*,*) 'lambdal, lambdar', lambdal, lambdar
             lambda_max = MAX(ABS(lambdal), ABS(lambdar))
             !write(*,*) rhol, ul, pl, iel
             !write(*,*) rhol, rhor, lambdal, lambdar
             dijL%aa(p) = norm_cij*lambda_max
          ELSE
             dijL%aa(p) = 0.d0
          END IF
       END DO
    END DO
    
    !===More viscosity for nonconvex pressure
    IF ((inputs%max_viscosity)=='global_max') THEN
       pstar = MAXVAL(dijL%aa)
       dijL%aa = pstar
       DO i = 1, mesh%np
          dijL%aa(diag(i)) = 0.d0
       END DO
    ELSE IF((inputs%max_viscosity)=='local_max') THEN
       DO i = 1, mesh%np
          pstar = MAXVAL(dijL%aa(dijL%ia(i):dijL%ia(i+1) - 1))
          dijL%aa(dijL%ia(i):dijL%ia(i+1) - 1) = pstar
          dijL%aa(diag(i)) =0.d0
       END DO
    END IF
    CALL transpose_op(dijL,'max')
    !===End More viscosity
    
    DO i = 1, mesh%np
       dijL%aa(diag(i)) = -SUM(dijL%aa(dijL%ia(i):dijL%ia(i+1)-1))
    END DO
    
  END SUBROUTINE compute_dij


  SUBROUTINE entropy_residual(un,rkgal)
    USE mesh_handling
    USE boundary_conditions
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2), INTENT(IN)  :: un, rkgal
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2)              :: DS
    REAL(KIND=8), DIMENSION(mesh%np)  :: s, e, u2, res, en, press, pg, &
         absres1, absres2
    INTEGER :: comp, k, i, j, p
    REAL(KIND=8) :: zz
    REAL(KIND=8) :: small_res
    u2 = 0.d0
    DO k = 1, k_dim 
       u2 = u2 + (un(:,k+1)/un(:,1))**2
    END DO

    SELECT CASE(inputs%high_order_viscosity)
    CASE ('EV(p)')
       press = (gamma-1.d0)*(un(:,k_dim+2) - un(:,1)*u2/2)
       pg = (1/gamma)*press**(1.d0/gamma-1.d0)
       en = gamma*press*pg  !S=Pressure*(1/gamma)
       s = en/un(:,1)
       DS(:,1) = pg*(gamma-1.d0)*u2/2
       DO k = 1, k_dim 
          DS(:,k+1) = -pg*(gamma-1.d0)*(un(:,k+1)/un(:,1))
       END DO
       DS(:,k_dim+2) = pg*(gamma-1.d0)
    CASE ('EV(s)')
       e = un(:,k_dim+2)/un(:,1) - u2/2
       s = (1.d0/(gamma-1d0))*LOG(e) - LOG(un(:,1))
       DS(:,1) = s + (1.d0/(gamma-1d0))*(u2/(2*e)-gamma)
       DO k = 1, k_dim 
          DS(:,k+1) = -un(:,k+1)/((gamma-1)*un(:,1)*e)
       END DO
       DS(:,k_dim+2) = (1.d0/(gamma-1d0))/e
       en = un(:,1)*s  !S = rho*s
    CASE DEFAULT
       WRITE(*,*) ' Bug: high_order_viscosity', inputs%high_order_viscosity
       STOP
    END SELECT

    !================TEST IN PAPER DONE WITH THIS SETTING
    !================DO NOT CHANGE
    res = rkgal(:,1)*(DS(:,1)- s)
    absres1 = ABS(res)
    DO comp = 2, k_dim+2
       res = res + rkgal(:,comp)*DS(:,comp)
       absres1 =  absres1 + ABS(rkgal(:,comp)*DS(:,comp))
    END DO
    !===It is essential to take the absolute value on each components
    !===to get correct convergence in 1D on the expansion wave.
    !===This also gives the correct scaling for VdW in 1D.
    DO i = 1, mesh%np
       zz = 0.d0
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          DO k = 1, k_dim
             zz = zz + cij(k)%aa(p)*un(j,k+1)*(s(j)-s(i))
          END DO
       END DO
       res(i) = res(i) + zz
       absres2(i) = ABS(zz) !==Essential to have this normalization in 2D
    END DO
    
    small_res = small*MAXVAL(abs(res))
    !================TEST IN PAPER DONE WITH THIS SETTING
    !================DO NOT CHANGE
    s = MIN(1.d0,inputs%ce*abs(res)/(absres1+absres2+small_res))
    
    !================TEST IN PAPER DONE WITH THIS SETTING
    !================DO NOT CHANGE
    DO i = 1, mesh%np
       DO p = dijH%ia(i), dijH%ia(i+1) - 1
          j = dijH%ja(p)
          dijH%aa(p) = dijL%aa(p)*((s(i)+s(j))/2)
       END DO
    END DO
    !================TEST IN PAPER DONE WITH THIS SETTING.
    !================DO NOT CHANGE
    IF (inputs%time+inputs%dt.GE.inputs%Tfinal) THEN
       SELECT CASE(k_dim)
       CASE(1)
          CALL plot_1d(mesh%rr(1,:),abs(s),'res.plt')
       CASE DEFAULT
          CALL plot_scalar_field(mesh%jj, mesh%rr, abs(s), 'res.plt')
       END SELECT
    END IF
  END SUBROUTINE entropy_residual

  SUBROUTINE smoothness_viscosity(un)
    USE mesh_handling
    USE boundary_conditions
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(mesh%np) :: alpha, bbeta, u2, e, s, ent, press, pg, num, ddenom, denom
    INTEGER      :: i, j, k, p

    u2 = 0.d0
    DO k = 1, k_dim 
       u2 = u2 + (un(:,k+1)/un(:,1))**2
    END DO
    SELECT CASE (inputs%high_order_viscosity)
    CASE('SM(p)')
       press = (gamma-1.d0)*(un(:,k_dim+2) - un(:,1)*u2/2)
       pg = (1/gamma)*press**(1.d0/gamma-1.d0)
       ent = gamma*press*pg  !=== S=Pressure*(1/gamma)
    CASE('SM(s)') 
       e = un(:,k_dim+2)/un(:,1) - u2/2
       s = (1.d0/(gamma-1d0))*LOG(e) - LOG(un(:,1))
       ent = un(:,1)*s  !=== S=rho*s
    CASE('SM(rho)')  !===Density
       ent = un(:,1)
    CASE DEFAULT
       WRITE(*,*) ' BUG in smoothness_viscosity'
    END SELECT

    DO i = 1, mesh%np
       denom(i) =0.d0
       ddenom(i) =0.d0
       num(i)  = 0.d0
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          IF (i==j) CYCLE
          num(i)   = num(i)         + stiff%aa(p)*(ent(j) - ent(i))
          denom(i)  = denom(i)  + ABS(stiff%aa(p)*(ent(j) - ent(i)))
          ddenom(i) = ddenom(i) + ABS(stiff%aa(p))*(ABS(ent(j)) + ABS(ent(i)))
       END DO
       num(i) = num(i)/lumped(i)
       denom(i) = denom(i)/lumped(i)
       ddenom(i) = ddenom(i)/lumped(i)
       IF (denom(i).GE.1.d-7*ddenom(i)) THEN
          alpha(i) = ABS(num(i))/(denom(i))
          bbeta(i) = ABS(num(i))/(ddenom(i))
       ELSE
          alpha(i) = 0.d0
          bbeta(i) = 0.d0
       END IF
    END DO

    DO i = 1, mesh%np
       DO p = dijL%ia(i), dijL%ia(i+1) - 1
          j = dijL%ja(p)
          dijH%aa(p)= dijL%aa(p)*MAX(MIN(alpha(i),bbeta(i)),MIN(alpha(j),bbeta(j)))
          !dijH%aa(p)= dijL%aa(p)*MAX(alpha(i),alpha(j))
          !dijH%aa(p)= dijL%aa(p)*MAX(bbeta(i),bbeta(j))
       END DO
    END DO

    !===Plot
    IF (inputs%time+inputs%dt.GE.inputs%Tfinal) THEN
       SELECT CASE(k_dim)
       CASE(1)
          CALL plot_1d(mesh%rr(1,:),alpha,'alpha.plt')
       CASE DEFAULT
          CALL plot_scalar_field(mesh%jj, mesh%rr, alpha, 'alpha.plt')
       END SELECT
    END IF
    !===Plot

  END SUBROUTINE smoothness_viscosity

  SUBROUTINE LOCAL_limit(unext,maxn,minn,mat,mass,comp)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: mass, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)   :: mat
    INTEGER,                    INTENT(IN)   :: comp
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2) :: unext
    REAL(KIND=8), PARAMETER :: smallplus = 1.d-15
    REAL(KIND=8) :: maxni, minni, ui, uij, xij, lambdai, umax, usmall
    INTEGER      :: i, j, p

    !===Compute lij
    umax = MAX(MAXVAL(maxn),-MINVAL(minn))
    usmall = umax*smallplus
    lij%aa = 1.d0
    DO i = 1, mesh%np
       lambdai = 1.d0/(mat%ia(i+1) - 1.d0 - mat%ia(i))
       maxni = maxn(i)
       minni = minn(i)
       ui = unext(i,comp)
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          xij = mat%aa(p)/(mass(i)*lambdai)
          uij = ui + xij
          IF (uij>maxni) THEN
             lij%aa(p) = MIN(ABS(maxni - ui)/(ABS(xij)+usmall),1.0)
          ELSE IF (uij<minni) THEN
             lij%aa(p) = MIN(ABS(minni - ui)/(ABS(xij)+usmall),1.d0)
          END IF
       END DO
    END DO
    lij%aa(diag) = 0.d0
  END SUBROUTINE LOCAL_limit


  SUBROUTINE transpose_op(mat,TYPE)
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(INOUT):: mat
    CHARACTER(LEN=3),  INTENT(IN)   :: TYPE
    INTEGER, DIMENSION(SIZE(mat%ia)) :: iao
    INTEGER:: i, j, p, next
    IF  (TYPE/='min' .AND. TYPE/='max') THEN
       WRITE(*,*) ' BUG in tanspose_op'
       STOP
    END IF
    iao = mat%ia
    DO i = 1, SIZE(mat%ia)-1
       DO p = mat%ia(i), mat%ia(i+1)-1 
          j = mat%ja(p)
          next = iao(j)
          iao(j) = next+1
          IF (j.LE.i) CYCLE
          IF (TYPE=='min') THEN
             mat%aa(next) = MIN(mat%aa(p),mat%aa(next))
             mat%aa(p) = mat%aa(next)
          ELSE IF (TYPE=='max') THEN
             mat%aa(next) = MAX(mat%aa(p),mat%aa(next))
             mat%aa(p) = mat%aa(next)
          END IF
       END DO
    END DO
  END SUBROUTINE transpose_op

  SUBROUTINE IDP_COMPUTE_DT(u0)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size), INTENT(IN) :: u0
    CALL compute_dij(u0)
    inputs%dt = inputs%CFL*1/MAXVAL(ABS(dijL%aa(diag))/lumped)
  END SUBROUTINE IDP_COMPUTE_DT

  SUBROUTINE FCT_generic(ulow,maxn,minn,mat,dg,comp)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: dg, maxn, minn
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    INTEGER,                    INTENT(IN)  :: comp
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2):: ulow
    REAL(KIND=8), DIMENSION(mesh%np)        :: Qplus, Qminus, Pplus, Pminus, Rplus, Rminus
    REAL(KIND=8), PARAMETER :: smallplus = small, smallminus = -small
    REAL(KIND=8) :: fij
    INTEGER      :: i, j, p, jp, jm
    Qplus  = dg*(maxn-ulow(:,comp))
    Qminus = dg*(minn-ulow(:,comp))
    Pplus  = smallplus
    Pminus = smallminus
    DO i = 1, mesh%np

       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          jp = 0
          jm =0
          IF (fij.GE.0.d0) THEN
             jp = jp + 1
             Pplus(i)  = Pplus(i) + fij
          ELSE
             jm = jm + 1
             Pminus(i) = Pminus(i) + fij
          END IF
       END DO
       IF (jp>0) THEN
          Rplus(i)  =  MIN(Qplus(i)/Pplus(i),1.d0)
       ELSE
          RPLUS(i) = 1.d0
       END IF
       IF (jm>0) THEN
          Rminus(i) =  MIN(Qminus(i)/Pminus(i),1.d0)
       ELSE
          Rminus(i) =1.d0
       END IF
    END DO

    DO i = 1, mesh%np
       DO p = mat%ia(i), mat%ia(i+1) - 1
          j = mat%ja(p)
          fij = mat%aa(p)
          IF (fij.GE.0.d0) THEN
             lij%aa(p) = MIN(Rplus(i),Rminus(j))
          ELSE
             lij%aa(p) = MIN(Rminus(i),Rplus(j))
          END IF
       END DO
    END DO
  END SUBROUTINE FCT_generic


  SUBROUTINE compute_greedy_dij(un)
    USE boundary_conditions
    USE lambda_module
    IMPLICIT NONE
    REAL(KIND=8):: cmin
    REAL(KIND=8), DIMENSION(mesh%np,k_dim+2), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(mesh%np,k_dim,k_dim+2)       :: vv
    REAL(KIND=8), DIMENSION(k_dim)                       :: nij
    REAL(KIND=8), DIMENSION(k_dim+2)                     :: Plr, usafe
    REAL(KIND=8)                                         :: p_min, rho_min, p_max, rho_max
    REAL(KIND=8)                                         :: romin, romax
    REAL(KIND=8) :: mr, ml, rhor, rhol, er, el, pr, pl, ul, ur, rhosmall, norm_cij, E_small, E_small_relax, ttt
    REAL(KIND=8) :: lambda_max, lambda_max_safe, lambda_max_small, lbd1, lbd2
    REAL(KIND=8) :: Sl, Sr, S_small, S_small_relax, aaa_entrop, bbb_entrop
    REAL(KIND=8) :: small_relax = 1.d-12, small_Newton = 1.d-5
    INTEGER :: i, j, p, k, comp

    DO comp = 1, k_dim+2
       vv(:,:,comp)=flux(comp,un)
    END DO

    DO i = 1, mesh%np
       DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
          j = cij(1)%ja(p)
          IF (i.EQ.j) THEN
             dijL%aa(p) = 0.d0
             CYCLE
          END IF
          DO k = 1, k_dim
             nij(k) = cij(k)%aa(p)
          END DO
          norm_cij = SQRT(SUM(nij**2))
          nij=nij/norm_cij
          mr = SUM(un(j,2:k_dim+1)*nij)
          ml = SUM(un(i,2:k_dim+1)*nij)
          rhor = un(1,j)
          rhol = un(1,i)
          er = un(j,k_dim+2) - 0.5d0*(SUM(un(j,2:k_dim+1)**2) - mr**2)/rhor
          el = un(i,k_dim+2) - 0.5d0*(SUM(un(i,2:k_dim+1)**2) - ml**2)/rhol
          ul = ml/rhol
          ur = mr/rhor
          pr = ABS(er-0.5d0*rhor*ur**2)*(gamma-1.d0)
          pl = ABS(el-0.5d0*rhol*ul**2)*(gamma-1.d0)
          !===Compute p_min,rho_min,p_max,rho_max,lambda_max_safe
          CALL ptilde(gamma,rhol,ul,pl,rhor,ur,pr,p_min,rho_min,p_max,rho_max,lambda_max_safe)
          rhosmall = rho_max*small_relax
          lambda_max_small = lambda_max_safe*small_relax
          !===Compute bar states
          DO comp = 1, k_dim+2
             Plr(comp)= - SUM((vv(j,:,comp)-vv(i,:,comp))*nij)/2.d0 
          END DO
          usafe = (un(i,:)+un(j,:))/2.d0 + Plr/lambda_max_safe
          romin = MIN(rho_min,usafe(1))
          romax = MAX(rho_max,usafe(1))
          !===First lambda based on density
          !IF (ABS(rhol+rhor-2.d0*rho_min).LE.rhosmall) THEN
          !    lbd1 = lambda_max_small
          IF (ABS(rhol+rhor-2.d0*romin).LE.rhosmall) THEN
             !lbd1 = lambda_max_small
             lbd1 = lambda_max_safe
          ELSE
             !lbd1 = (mr-ml)/(rhol+rhor-2.d0*rho_min)
             lbd1 = lambda_max_safe*(rhol+rhor-2.d0*usafe(1))/(rhol+rhor-2.d0*romin)
          END IF
          !IF (ABS(rhol+rhor-2.d0*rho_max).LE.rhosmall) THEN
          !   lbd2 = lambda_max_small
          IF (ABS(rhol+rhor-2.d0*romax).LE.rhosmall) THEN
             !lbd2 = lambda_max_small
             lbd2 = lambda_max_safe
          ELSE
             !lbd2 = (ml-mr)/(2.d0*rho_max-rhol-rhor)
             lbd2 = lambda_max_safe*(2.d0*usafe(1)-rhol-rhor)/(2.d0*romax-rhol-rhor)
          END IF
          lambda_max = MAX(lbd1, lbd2, lambda_max_small)
          IF (lambda_max>lambda_max_safe) THEN
             WRITE(*,*) ' THERE IS A BUG 1', lambda_max, lambda_max_safe
          END IF
          !===Second lambda based on internal energy rho e
          !===Third lambda based on minimum principle on entropy
          !===Compute cmin
          cmin = (1.d0/(gamma-1.d0))*MIN(pl/rhol**(gamma),pr/rhor**(gamma)) !===gamma is defined in lambda_module
          IF (cmin>0.d0) THEN
             cmin = cmin*(1-small_relax)
          ELSE
             cmin = cmin*(1+small_relax)
          END IF
          E_small = small_Newton*(un(i,k_dim+2)+un(j,k_dim+2))
          E_small_relax = small_relax*(un(k_dim+2,i)+un(k_dim+2,j))
          IF (lambda_max.GE.lambda_max_small) THEN
             ttt = 1.d0/lambda_max
          ELSE
             WRITE(*,*) ' lambda_max seems to be equal to zero', lambda_max
             ttt = 1.d0/lambda_max_small
             STOP
          END IF
          !if (ttt.LE.0.d0) THEN
          !   WRITE(*,*) 'lambda_max==0.d0', lambda_max, ttt
          !   STOP
          !END if
          CALL Newton_secant(un(i,:),un(j,:),Plr,ttt,psi_entrop,psi_entrop_prime,E_small)
          IF (ttt.LE.0.d0) THEN
             WRITE(*,*) ' BUG ttt<0', ttt
             STOP
          END IF
          lambda_max = MAX(lambda_max,1.d0/ttt)
          !IF (lambda_max>lambda_max_safe*(1.d0+small)) THEN
          !   WRITE(*,*) ' THERE IS A BUG 3', lambda_max, lambda_max_safe
          !END IF
          go to 10
          !===Fourth lambda based on entropy inequality
          Sl = entrop(un(i,:),0d0,0.d0,0.d0)
          Sr = entrop(un(j,:),0d0,0.d0,0.d0)
          aaa_entrop = -(Sl+Sr)/2.d0
          bbb_entrop = (ur*Sr-ul*Sl)/2.d0
          IF (aaa_entrop>0.d0) THEN
             aaa_entrop = aaa_entrop*(1-small_relax)
          ELSE
             aaa_entrop = aaa_entrop*(1+small_relax)
          END IF
          S_small = small_Newton*(ABS(Sl)+ABS(Sr))
          S_small_relax = small_relax*(ABS(Sl)+ABS(Sr))
          ttt = 1.d0/lambda_max
          CALL entrop_ineq(un(i,:),un(j,:),Plr,ttt,aaa_entrop,bbb_entrop,S_small)
          lambda_max = MAX(lambda_max,1.d0/ttt)
          !IF (lambda_max>lambda_max_safe*(1.d0+small)) THEN
          !   WRITE(*,*) ' THERE IS A BUG 4', lambda_max, lambda_max_safe
          !END IF
          !===Definition of dij
10        dijL%aa(p) = norm_cij*lambda_max
       END DO
    END DO
    CALL transpose_op(dijL,'max')
    DO i = 1, mesh%np
       dijL%aa(diag(i)) = -SUM(dijL%aa(dijL%ia(i):dijL%ia(i+1)-1))
    END DO

  CONTAINS
    FUNCTION psi_entrop(u) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u
      REAL(KIND=8)                                 :: psi
      psi = u(k_dim+2) - SUM(u(2:k_dim+1)**2)/(2.d0*u(1)) - cmin*u(1)**gamma
    END FUNCTION psi_entrop

    FUNCTION psi_entrop_prime(Pij,u) RESULT(psi)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u, Pij
      REAL(KIND=8)                                 :: psi
      psi = Pij(k_dim+2) - SUM(u(2:k_dim+1)*Pij(2:k_dim+1))/u(1) &
           + Pij(1)*SUM(u(2:k_dim+1)**2)/(2*u(1)**2) &
           - cmin*gamma*Pij(1)*u(1)**(gamma-1.d0)
    END FUNCTION psi_entrop_prime
  END SUBROUTINE compute_greedy_dij


  SUBROUTINE Newton_secant(ui, uj, Pij,limiter,psi_func,psi_prime_func,psi_small)
    IMPLICIT NONE
    INTERFACE
       FUNCTION psi_func(u) RESULT(psi)
         USE space_dim
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN)  :: u
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_func
       FUNCTION psi_prime_func(Pij,u) RESULT(psi)
         USE space_dim
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN)  :: u, Pij
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_prime_func
    END INTERFACE
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: Pij
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: ui, uj
    REAL(KIND=8), INTENT(IN)    :: psi_small
    REAL(KIND=8), INTENT(INOUT) :: limiter
    REAL(KIND=8), DIMENSION(k_dim+2) ::  ul, ur, ulr
    REAL(KIND=8) :: psil, psir, ll, lr, llold, lrold
    LOGICAL :: once

    ul = ui
    ur = uj
    ulr = (ul+ur)/2

    lr = limiter
    ur = ulr + lr*Pij
    psir = psi_func(ur)

    IF (psir.GE.0.d0) THEN 
       !===input limiter is okay
       RETURN
    END IF
    ll = 0.d0
    ul = ulr
    psil = MAX(psi_func(ulr),0.d0)

    once=.TRUE.
    DO WHILE (ABS(psil-psir) .GT. psi_small .OR. once)
       once =.FALSE.
       llold = ll
       lrold = lr
       ll = ll - psil*(lr-ll)/(psir-psil)
       lr = lr - psir/psi_prime_func(Pij,ur)
       IF (ll.GE.lr) THEN
          ll = lr
          EXIT
       END IF
       IF (ll< llold) THEN
          ll = llold
          EXIT
       END IF
       IF (lr > lrold) THEN
          lr = lrold
          EXIT
       END IF
       ul = ulr + ll*Pij
       ur = ulr + lr*Pij
       psil = psi_func(ul)
       psir = psi_func(ur)
    END DO

    IF (psir.GE.0.d0) THEN 
       limiter = lr
    ELSE
       limiter = ll
    END IF

  END SUBROUTINE Newton_secant

  SUBROUTINE entrop_ineq(ui,uj,Pij,limiter,aaa,bbb,psi_small)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: Pij
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: ui, uj
    REAL(KIND=8), INTENT(IN)    :: aaa, bbb, psi_small
    REAL(KIND=8), INTENT(INOUT) :: limiter
    REAL(KIND=8), DIMENSION(k_dim+2) ::  ul, ur, ulr
    REAL(KIND=8) :: psil, psir, ll, lr, llold, lrold
    LOGICAL      :: once
    ul = ui
    ur = uj
    ulr = (ul+ur)/2
    !IF (entrop(ulr,aaa,bbb,0.d0)>0.d0) THEN
    !   write(*,*) ' BUG , entropy ineq violated by (ul+ur)/2 at 0', entrop(ulr,aaa,bbb,0.d0)
    !   stop
    !END IF

    lr = limiter
    ur = ulr + lr*Pij
    psir = entrop(ur,aaa,bbb,lr)    
    IF (psir.LE. 0.d0) THEN 
       !===input limiter is okay
       RETURN
    END IF
    ll = 0.d0
    ul = ulr
    psil = entrop(ul,aaa,bbb,ll)

    once =.TRUE.
    DO WHILE (ABS(psil-psir) .GT. psi_small  .OR. once)
       once=.FALSE.
       llold = ll
       lrold = lr
       ll = ll - psil*(lr-ll)/(psir-psil)
       lr = lr - psir/entrop_prime(Pij,ur,bbb)
!!$       IF (ll.GE.lr) THEN
!!$          ll = lr
!!$          EXIT
!!$       END IF
!!$       IF (ll< llold) THEN
!!$          ll = llold
!!$          EXIT
!!$       END IF
!!$       IF (lr > lrold) THEN
!!$          lr = lrold
!!$          EXIT
!!$       END IF
       ul = ulr + ll*Pij
       ur = ulr + lr*Pij
       psil = entrop(ul,aaa,bbb,ll)
       psir = entrop(ur,aaa,bbb,lr)
    END DO
    IF (psir.LE.0.d0) THEN   
       limiter = lr
       ul = ur
    ELSE
       limiter = ll
    END IF

    IF (entrop(ulr+limiter*pij,aaa,bbb,limiter)>2*psi_small) THEN
       WRITE(*,*) ' BUG  entrop ineq. violated',  entrop(ulr+limiter*pij,aaa,bbb,limiter), ABS(psil-psir)
       STOP
    END IF
  END SUBROUTINE entrop_ineq

  FUNCTION entrop(u,aaa,bbb,t) RESULT(psi)
    USE boundary_conditions !===Access to gamma
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u
    REAL(KIND=8)                                 :: aaa, bbb, t, psi
    REAL(KIND=8) :: S
    S = (u(1)/(gamma-1.d0))*LOG(u(1)*u(k_dim+2) - SUM(u(2:k_dim+1)**2)/2.d0) &
         - ((gamma+1.d0)/(gamma-1.d0))*u(1)*LOG(u(1)) !===rho/(gamma-1)*log(rho*E-M**2/2) -(gamma+1)/(gamma-1)*rho*log(rho)
    psi = -S + aaa + t*bbb !===The negative of the entropy is convex
  END FUNCTION entrop

  FUNCTION entrop_prime(Pij,u,bbb) RESULT(psi)
    USE boundary_conditions !===Access to gamma
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(k_dim+2), INTENT(IN) :: u, Pij
    REAL(KIND=8)                                 :: bbb, psi
    REAL(KIND=8), DIMENSION(k_dim+2) :: DS
    REAL(KIND=8)                     :: u2, e, s
    INTEGER                          :: k 
    u2 = SUM(u(2:k_dim+1)**2)/(u(1)**2) 
    e = u(k_dim+2)/u(1) - u2/2
    s = (1.d0/(gamma-1d0))*LOG(e) - LOG(u(1))
    DS(1) = s + (1.d0/(gamma-1d0))*(u2/(2*e)-gamma)
    DO k = 1, k_dim 
       DS(k+1) = -u(k+1)/((gamma-1)*u(1)*e)
    END DO
    DS(k_dim+2) = (1.d0/(gamma-1d0))/e
    psi = -SUM(Pij*DS) + bbb !===The negative of the entropy is convex
  END FUNCTION entrop_prime

END MODULE IDP_update_euler
