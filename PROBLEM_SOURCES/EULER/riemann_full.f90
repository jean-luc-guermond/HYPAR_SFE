! Authors: Jean-Luc Guermond and Bojan Popov, Texas A&M, April 5, 2016
MODULE lambda_module_full
  PUBLIC                  :: lambda_full
  PRIVATE
  INTEGER,      PARAMETER :: Mgas=3
  REAL(KIND=8), PARAMETER :: gamma=(Mgas+2.d0)/Mgas, exp=1.d0/(Mgas+2), b=0.0d0
  REAL(KIND=8)            :: al, capAl, capBl, covl, ar, capAr, capBr, covr
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

  SUBROUTINE lambda_full(tol,rhol,ul,pl,rhor,ur,pr,lambdal,lambdar,k)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: tol, rhol, ul, pl, rhor, ur, pr
    REAL(KIND=8), INTENT(OUT):: lambdal, lambdar
    INTEGER, INTENT(OUT)     :: k
    REAL(KIND=8)             :: lambda_max, pstar
    REAL(KIND=8)             :: capAmin, capBmin, acovmin, acovmax, ratio
    REAL(KIND=8)             :: lambda_min, phimax, phimin, ptilde
    REAL(KIND=8)             :: phi1, phi11, phi12, phi112, phi2, phi22, phi221
    REAL(KIND=8)             :: p1, p2 , pmin, pmax, rhomin, rhomax, v11, v12, v31, v32
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
    IF (phimin.GE.0) THEN
       pstar = 0.d0
       lambdal = ul-al
       lambdar = ur+ar
       RETURN
    END IF
    phimax = (pmax-pmin)*SQRT(capAmin/(pmax+capBmin)) + ur-ul
    ptilde = pmin*((acovmin+acovmax - (ur-ul)*(gamma-1)/2)&
         /(acovmin + acovmax*ratio))**(Mgas+2)
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
    lambda_max = MAX(MAX(v32,0.d0),MAX(-v11,0.d0))
    lambda_min = MAX(MAX(MAX(v31,0.d0),MAX(-v12,0.d0)),0.d0)
    IF (lambda_min>0.d0) THEN
       IF (lambda_max/lambda_min -1.d0 .LE. tol) THEN
          lambdal = v11
          lambdar = v32
          RETURN
       END IF
    END IF
    !lambdal = ul - al*SQRT(1.d0+MAX((p2-pl)/pl,0.d0)*(1-exp))
    !lambdar = ur + ar*SQRT(1.d0+MAX((p2-pr)/pr,0.d0)*(1-exp))
    !return
    p1 = MAX(p1,p2-phi(p2,ul,pl,ur,pr)/phi_prime(p2,pl,pr))
    !===Iterations
    DO WHILE(.TRUE.) 
       v11 = lambdaz(ul,pl,al/covl,p2,-1)
       v12 = lambdaz(ul,pl,al/covl,p1,-1)
       v31 = lambdaz(ur,pr,ar/covr,p1,1)
       v32 = lambdaz(ur,pr,ar/covr,p2,1)
       lambda_max = MAX(MAX(v32,0.d0),MAX(-v11,0.d0))
       lambda_min = MAX(MAX(MAX(v31,0.d0),MAX(-v12,0.d0)),0.d0)
       IF (lambda_min>0.d0) THEN
          IF (lambda_max/lambda_min -1.d0 .LE. tol) THEN
             lambdal = v11
             lambdar = v32
             RETURN
          END IF
       END IF
       phi1 =  phi(p1,ul,pl,ur,pr)
       phi11 = phi_prime(p1,pl,pr)
       phi2 =  phi(p2,ul,pl,ur,pr)
       phi22 = phi_prime(p2,pl,pr)
       IF (phi1>0.d0) THEN
          lambdal = v12
          lambdar = v31
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
  END SUBROUTINE lambda_full

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
END MODULE lambda_module_full

