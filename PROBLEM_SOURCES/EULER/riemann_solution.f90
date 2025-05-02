! Authors: Jean-Luc Guermond and Bojan Popov, Texas A&M, April 5, 2016
MODULE riemann_solution_module
  IMPLICIT NONE
  PUBLIC                           :: lambda, riemann_solution_at_zero
  !INTEGER,      PRIVATE, PARAMETER :: Mgas=5 ! 3 for leblanc and 5 for the others
  !REAL(KIND=8), PUBLIC, PARAMETER :: gamma=(Mgas+2.d0)/Mgas, exp=1.d0/(Mgas+2)
  PRIVATE
  REAL(KIND=8) :: gamma, exp
  REAL(KIND=8) :: b=0.0d0
  REAL(KIND=8) :: al, capAl, capBl, covl, ar, capAr, capBr, covr
CONTAINS

  SUBROUTINE init(rhol,pl,rhor,pr)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: rhol, pl, rhor, pr
    al = SQRT(gamma*pl/rhol)
    capAl = 2/((gamma+1)*rhol)
    capBl = pl*(gamma-1)/(gamma+1)
    covl = SQRT(1-b*rhol)
    ar = SQRT(gamma*pr/rhor)
    capAr = 2/((gamma+1)*rhor)
    capBr = pr*(gamma-1)/(gamma+1)
    covr = SQRT(1-b*rhor)
  END SUBROUTINE init

  SUBROUTINE lambda(gamma_in,tol,rhol,ul,pl,rhor,ur,pr,lambda_max,pstar,k,v11,v32)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: gamma_in, tol, rhol, ul, pl, rhor, ur, pr
    REAL(KIND=8), INTENT(OUT):: lambda_max, pstar
    INTEGER,      INTENT(OUT):: k
    REAL(KIND=8)             :: capAmin, capBmin, acovmin, acovmax, ratio
    REAL(KIND=8)             :: phimax, phimin, ptilde
    REAL(KIND=8)             :: phi1, phi11, phi12, phi112, phi2, phi22, phi221
    REAL(KIND=8)             :: p1, p2 , pmin, pmax, rhomin, rhomax, v11, v12, v31, v32
    REAL(KIND=8)             :: err1, err3, aquad, bquad, cquad, xl, xr

    gamma = gamma_in
    exp = (gamma-1)/(2*gamma)
    !===Initialization
    CALL init(rhol,pl,rhor,pr)
    k = 0
    IF (pl.LE.pr) THEN
       pmin   = pl
       rhomin = rhol
       pmax   = pr
       rhomax = rhor
    ELSE
       pmin   = pr
       rhomin = rhor
       pmax   = pl
       rhomax = rhol
    END IF
    capAmin = 2/((gamma+1)*rhomin)
    capBmin = pmin*(gamma-1)/(gamma+1)
    acovmin = SQRT(gamma*pmin*(1-b*rhomin)/rhomin)
    acovmax = SQRT(gamma*pmax*(1-b*rhomax)/rhomax)
    ratio = (pmin/pmax)**exp
    phimin = (2/(gamma-1.d0))*acovmax*(ratio-1.d0) + ur-ul
    phimax = (pmax-pmin)*SQRT(capAmin/(pmax+capBmin)) + ur-ul
    ptilde = pmin*((acovmin+acovmax - (ur-ul)*(gamma-1)/2)&
         /(acovmin + acovmax*ratio))**(1/exp)
    IF (phimax<0.d0) THEN
       xl = SQRT(capAl/(1+capBl/pmax))
       xr = SQRT(capAr/(1+capBr/pmax))
       aquad = xl+xr
       bquad = ur-ul
       cquad = -pl*xl-pr*xr
       p2 = ((-bquad+SQRT(bquad**2-4*aquad*cquad))/(2*aquad))**2
       ptilde = min(ptilde,p2)
    END IF
    IF (phimax < 0.d0) THEN
       p1 = pmax
       p2 = ptilde
    ELSE
       p1=pmin
       p2 = MIN(pmax,ptilde)
    END IF
    !===Check for accuracy after initialization
    v11 = lambdaz(ul,pl,al/covl,p2,-1)
    v12 = lambdaz(ul,pl,al/covl,p1,-1)
    v31 = lambdaz(ur,pr,ar/covr,p1,1)
    v32 = lambdaz(ur,pr,ar/covr,p2,1)
    !===Test on full solution
    lambda_max = MAX(MAX(v32,0.d0),MAX(-v11,0.d0))
    err3 =  abs(v32 - v31)/lambda_max
    err1 =  abs(v12 - v11)/lambda_max
    IF (MAX(err1,err3).LE.tol) THEN
       pstar = p2
       RETURN
    END IF
    !pstar =p2
    !return
    !===Full solution
    p1 = MAX(p1,p2-phi(p2,ul,pl,ur,pr)/phi_prime(p2,pl,pr))
    !===Iterations
    DO WHILE(.TRUE.)
       !===Test on full solution
       v11 = lambdaz(ul,pl,al/covl,p2,-1)
       v12 = lambdaz(ul,pl,al/covl,p1,-1)
       v31 = lambdaz(ur,pr,ar/covr,p1,1)
       v32 = lambdaz(ur,pr,ar/covr,p2,1)
       lambda_max = MAX(MAX(v32,0.d0),MAX(-v11,0.d0))
       err3 =  abs(v32 - v31)/lambda_max
       err1 =  abs(v12 - v11)/lambda_max
       IF (MAX(err1,err3).LE.tol) THEN
          pstar = p2
          RETURN
       END IF
       !===Full solution
       phi1 =  phi(p1,ul,pl,ur,pr)
       phi11 = phi_prime(p1,pl,pr)
       phi2 =  phi(p2,ul,pl,ur,pr)
       phi22 = phi_prime(p2,pl,pr)
       IF (phi1>0.d0) THEN
          !lambda_max = lambda_min
          RETURN
       END IF
       IF (phi2<0.d0) RETURN
       phi12 = (phi2-phi1)/(p2-p1) 
       phi112 = (phi12-phi11)/(p2-p1)
       phi221 = (phi22-phi12)/(p2-p1)
       p1 = p1 - 2*phi1/(phi11 + SQRT(phi11**2 - 4*phi1*phi112))
       p2 = p2 - 2*phi2/(phi22 + SQRT(phi22**2 - 4*phi2*phi221))
       k = k+1
    END DO
  END SUBROUTINE lambda

  FUNCTION lambdaz(uz,pz,az,pstar,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: uz,pz,az,pstar
    INTEGER,      INTENT(IN) :: z
    REAL(KIND=8)             :: vv
    vv = uz + z*az*SQRT(1+MAX((pstar-pz)/pz,0.d0)*(gamma+1)/(2*gamma))
  END FUNCTION lambdaz

  FUNCTION phi(p,ul,pl,ur,pr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: p, ul, pl, ur, pr
    REAL(KIND=8)             :: vv, fl, fr
    IF (p>pl) THEN
       fl = (p-pl)*SQRT(capAl/(p+capBl))
    ELSE
       fl = (2*al/(gamma-1))*((p/pl)**exp-1)
    END IF
    IF (p>pr) THEN
       fr = (p-pr)*SQRT(capAr/(p+capBr))
    ELSE
       fr = (2*ar/(gamma-1))*((p/pr)**exp-1)
    END IF
    vv = fl*covl + fr*covr + ur - ul
  END FUNCTION phi

  FUNCTION phi_prime(p,pl,pr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: p, pl, pr
    REAL(KIND=8)             :: vv, fl, fr
    IF (p>pl) THEN
       fl = SQRT(capAl/(p+capBl))*(1-(p-pl)/(2*(capBl+p)))
    ELSE
       fl = (al/(gamma*pl))*(p/pl)**(-(gamma+1)/(2*gamma))
    END IF
    IF (p>pr) THEN
       fr = SQRT(capAr/(p+capBr))*(1-(p-pr)/(2*(capBr+p)))
    ELSE
       fr = (ar/(gamma*pr))*(p/pr)**(-(gamma+1)/(2*gamma))
    END IF
    vv = fl*covl + fr*covr
  END FUNCTION phi_prime

  FUNCTION ustar(p,ul,pl,ur,pr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: p, ul, pl, ur, pr
    REAL(KIND=8)             :: vv, fl, fr
    IF (p>pl) THEN
       fl = (p-pl)*SQRT(capAl/(p+capBl))
    ELSE
       fl = (2*al/(gamma-1))*((p/pl)**exp-1)
    END IF
    IF (p>pr) THEN
       fr = (p-pr)*SQRT(capAr/(p+capBr))
    ELSE
       fr = (2*ar/(gamma-1))*((p/pr)**exp-1)
    END IF
    vv = (-fl*covl + fr*covr + ur + ul)/2.d0
  END FUNCTION ustar

  SUBROUTINE riemann_solution_at_zero(pstar,rhol,ul,pl,rhor,ur,pr,rho0,u0,p0)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: pstar,rhol,ul,pl,rhor,ur,pr
    REAL(KIND=8), INTENT(OUT):: rho0,u0,p0
    REAL(KIND=8) :: al, lambda1_minus, rhostarl, astarl, ustarl, lambda1_plus
    REAL(KIND=8) :: ar, lambda3_minus, rhostarr, astarr, ustarr, lambda3_plus
    REAL(KIND=8) :: ustar
    !===left wave
    al = SQRT(gamma*pl/rhol)
    lambda1_minus = lambdaz(ul,pl,al,pstar,-1)
    IF (pl<pstar) THEN !===shock
       rhostarl = rhol*((pstar/pl) + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*(pstar/pl)+1)
       ustarl = ul - (pstar-pl)*SQRT(capAl/(pstar+capBl))
       lambda1_plus  = lambda1_minus
    ELSE
       rhostarl = rhol*(pstar/pl)**(1/gamma)
       astarl = al*(pstar/pl)**((gamma-1)/(2*gamma))
       ustarl = ul - (2*al/(gamma-1))*((pstar/pl)**exp-1)
       lambda1_plus  = ustarl - astarl
    END IF
    
    !===right wave
    ar = SQRT(gamma*pr/rhor)
    lambda3_plus  = lambdaz(ur,pr,ar,pstar,1)
    IF (pr<pstar) THEN !===shock
       rhostarr = rhor*((pstar/pr) + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*(pstar/pr)+1)
       ustarr = ur + (pstar-pr)*SQRT(capAr/(pstar+capBr))
       lambda3_minus  = lambda3_plus
    ELSE
       rhostarr = rhor*(pstar/pr)**(1/gamma)
       astarr = ar*(pstar/pr)**((gamma-1)/(2*gamma))
       ustarr = ur + (2*ar/(gamma-1))*((pstar/pr)**exp-1)
       lambda3_minus  = ustarr + astarr
    END IF
    ustar = (ustarl+ustarr)/2

    IF (0.d0<lambda1_minus) THEN
       rho0=rhol
       u0=ul
       p0=pl
    ELSE IF (0.d0<lambda1_plus) THEN
       rho0 = rhol*(2/(gamma+1) + (ul/al)*(gamma-1)/(gamma+1))**(2/(gamma-1))
       u0 = (2/(gamma+1))*(al+ul*(gamma-1)/2)
       p0 = pl*(rho0/rhol)**gamma
    ELSE IF (0.d0<ustar) THEN
       rho0 = rhostarl
       u0 = ustar
       p0 = pstar
    ELSE IF (0.d0<lambda3_minus) THEN
       rho0 = rhostarr
       u0 = ustar
       p0 = pstar
    ELSE IF (0.d0<lambda3_plus) THEN
       rho0 = rhor*(2/(gamma+1) - (ur/ar)*(gamma-1)/(gamma+1))**(2/(gamma-1))
       u0 = (2/(gamma+1))*(-ar+ur*(gamma-1)/2)
       p0 = pr*(rho0/rhor)**gamma
    ELSE
       rho0=rhor
       u0=ur
       p0=pr
    END IF
    RETURN
  END SUBROUTINE riemann_solution_at_zero
  
END MODULE riemann_solution_module

