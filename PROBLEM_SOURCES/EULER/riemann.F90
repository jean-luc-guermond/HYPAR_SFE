! Authors: Jean-Luc Guermond and Bojan Popov, Texas A&M, Jan 22, 2016
MODULE lambda_module
  PUBLIC :: lambda, ptilde
CONTAINS

  SUBROUTINE lambda(gamma,rhol,ul,pl,rhor,ur,pr,lambdal,lambdar)
    USE arbitrary_eos_lambda_module
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: rhol, ul, pl, rhor, ur, pr
    REAL(KIND=8), INTENT(OUT):: lambdal, lambdar 
    REAL(KIND=8)             :: pstar
    REAL(KIND=8)             :: pmin, pmax
    REAL(KIND=8)             :: amin, amax
    REAL(KIND=8)             :: al, ar
    REAL(KIND=8)             :: gamma, exponent, one_m_exponent
    REAL(KIND=8)             :: Mgas, Mgasptwo

    Mgas = 2/(gamma-1)
    Mgasptwo = Mgas+2
    exponent=(gamma-1)/(2*gamma)
    one_m_exponent = 1-exponent
    
    al=SQRT(gamma*pr/rhor)
    ar=SQRT(gamma*pl/rhol)
    pstar = al+ar-0.5d0*(gamma-1.0d0)*(ur-ul)
    if (pstar>0.0d0) THEN
       pmin=pr
       pmax=pl
       amin=ar
       amax=al
       if (pl < pr) THEN 
          pmin=pl
          pmax=pr
          amin=al
          amax=ar
       END if
       pstar = pstar/(amin + amax*(pmin/pmax)**exponent)
       pstar = pstar**(Mgasptwo)*pmin
    END if

    lambdal = ul - al*SQRT(1.d0+MAX((pstar-pl)/pl,0.d0)*one_m_exponent)
    lambdar = ur + ar*SQRT(1.d0+MAX((pstar-pr)/pr,0.d0)*one_m_exponent)
  END SUBROUTINE lambda

!!$  SUBROUTINE lambda_VDW(rhol,ul,el,rhor,ur,er,lambdal,lambdar)
!!$    USE boundary_conditions
!!$    IMPLICIT NONE
!!$    REAL(KIND=8) :: rhol,ul,el,rhor,ur,er,lambdal,lambdar
!!$    REAL(KIND=8) :: un(3,1), ppl(1), ppr(1)
!!$    REAL(KIND=8) :: pl, al, gammal, xl, capAl, capBl, pr, ar, gammar, xr, capAr, capBr
!!$    REAL(KIND=8) :: pmin, pmax, gammamin, gammamax, gammacapM, gammalwcm
!!$    REAL(KIND=8) :: num, denom, pstar, ratio, expo, alpha
!!$    REAL(KIND=8) :: amax, amin 
!!$    un(1,1) = rhol
!!$    un(2,1) = rhol*ul
!!$    un(3,1) = rhol*(el + 0.5*ul**2)
!!$    ppl = pressure(un)
!!$    pl = ppl(1)
!!$    gammal=pl/(el*rhol)+1
!!$    al = SQRT(gammal*pl/rhol)
!!$    capAL = 2/(rhol*(gammal+1))
!!$    capBL = pl*(gammal-1)/(gammal+1)
!!$    un(1,1) = rhor
!!$    un(2,1) = rhor*ur
!!$    un(3,1) = rhor*(er + 0.5*ur**2)
!!$    ppr = pressure(un)
!!$    pr =ppr(1)
!!$    gammar=pr/(er*rhor)+1
!!$    ar = SQRT(gammar*pr/rhor)
!!$    capAr = 2/(rhor*(gammar+1))
!!$    capBr = pl*(gammar-1)/(gammar+1)
!!$    if (pr.le.0 .or. pl.le.0) THEN
!!$       write(*,*) 'l ', rhol, ul, el, pl
!!$       write(*,*) 'r ', rhor, ur, er, pr
!!$       write(*,*) (gamma-1.d0)*(un(3,1)-0.5d0*(un(2,1)**2)/un(1,1) + un(1,1)**2)/(1-un(1,1))&
!!$            - un(1,:)**2
!!$       stop
!!$    end if
!!$    IF (pl<pr) THEN
!!$       pmin = pl
!!$       amin = al
!!$       gammamin = gammal
!!$       pmax = pr
!!$       amax = ar
!!$       gammamax = gammar
!!$    ELSE
!!$       pmin = pr
!!$       amin = ar
!!$       gammamin = gammar
!!$       pmax = pl
!!$       amax = al
!!$       gammamax = gammal
!!$    END IF
!!$    gammacapM = max(gammal,gammar)
!!$    gammalwcm = min(gammal,gammar)
!!$    IF (phi(pmin)>0) THEN
!!$       xl = 2*ar/(gammar-1)
!!$       xr = 2*al/(gammal-1)
!!$       num = xl+xr-(ur-ul)
!!$       denom = xl/pl**expg(1/gammacapM) + xr/pr**expg(1/gammacapM)
!!$       pstar = (num/denom)**expg(gammacapM)
!!$    ELSE IF(phi(pmax)>0) THEN
!!$       IF (gammamin==gammalwcm) THEN
!!$          expo = gammacapM
!!$          ratio = pmin/pmax
!!$       ELSE
!!$          expo = gammalwcm
!!$          ratio = 1
!!$       END IF
!!$       alpha = (gammacapM-gammalwcm)/(2*gammacapM*gammalwcm)
!!$       xl = 2*amin/(gammamin-1)
!!$       xr = 2*amax/(gammamax-1)
!!$       num = xl+xr-(ur-ul)
!!$       denom = xl*(1/pmin**expg(1/expo))*(ratio)**(alpha) + xr/pmax**expg(1/expo)
!!$       pstar = (num/denom)**expg(expo)
!!$    ELSE
!!$       xl = 2*ar/(gammar-1)
!!$       xr = 2*al/(gammal-1)
!!$       num = xl+xr-(ur-ul)
!!$       denom = xl/pl**expg(1/gammalwcm) + xr/pr**expg(1/gammalwcm)
!!$       pstar = (num/denom)**expg(gammalwcm)
!!$    END IF
!!$    lambdal = ul - al*SQRT(1.d0+MAX((pstar-pl)/pl,0.d0)*(gammal+1)/(2*gammal))
!!$    lambdar = ur + ar*SQRT(1.d0+MAX((pstar-pr)/pr,0.d0)*(gammar+1)/(2*gammar))
!!$  CONTAINS
!!$    FUNCTION expg(gamma) RESULT(vv)
!!$      IMPLICIT NONE
!!$      REAL(KIND=8) :: gamma, vv
!!$      vv = (gamma-1)/(2*gamma)
!!$    END FUNCTION expg
!!$    FUNCTION phi(p) RESULT(vv)
!!$      REAL(KIND=8) ::  p, vv
!!$      vv = ur-ul
!!$      IF (p<pl) THEN
!!$         vv = vv + ((p/pl)**expg(gammal)-1)*2*al/(gammal-1)
!!$      ELSE
!!$         vv = vv + (p-pl)*SQRT(capAl/(p+capBl))
!!$      END IF
!!$      IF (p<pr) THEN
!!$         vv = vv + ((p/pr)**expg(gammar)-1)*2*ar/(gammar-1)
!!$      ELSE
!!$         vv = vv + (p-pr)*SQRT(capAr/(p+capBr))
!!$      END IF
!!$    END FUNCTION phi
!!$  END SUBROUTINE lambda_VDW

  SUBROUTINE ptilde(gamma,rhol,ul,pl,rhor,ur,pr,p_min,rho_min,p_max,rho_max,lambda_max)
    !USE input_data
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: rhol, ul, pl, rhor, ur, pr
    REAL(KIND=8), INTENT(OUT):: p_min,rho_min,p_max,rho_max,lambda_max
    REAL(KIND=8)             :: pmin, pmax, rhomin, rhomax, amin, amax, phimin, phimax
    REAL(KIND=8)             :: al, capAl, capBl, ar, capAr, capBr, capAmin, capBmin, p1, p2
    REAL(KIND=8)             :: ratio, ptil, lambdal,  lambdar 
    REAL(KIND=8)             :: gamma, exponent, one_m_exponent, gm1_over_gp1
    REAL(KIND=8)             :: Mgas, Mgasptwo
 
    Mgas = 2/(gamma-1)
    Mgasptwo = Mgas+2
    exponent=(gamma-1)/(2*gamma)
    one_m_exponent = 1-exponent
    gm1_over_gp1 =(gamma-1)/(gamma+1)
    
    al = SQRT(gamma*pl/rhol)
    capAl = 2/((gamma+1)*rhol)
    capBl = pl*(gamma-1)/(gamma+1)

    ar = SQRT(gamma*pr/rhor)
    capAr = 2/((gamma+1)*rhor)
    capBr = pr*(gamma-1)/(gamma+1)

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
    amin = SQRT(gamma*pmin/rhomin)
    amax = SQRT(gamma*pmax/rhomax)
    ratio = (pmin/pmax)**exponent
    phimin = (2/(gamma-1.d0))*amax*(ratio-1.d0) + ur-ul
    !===Check phi(pmin)>0. Phi(p) defined by (4.1) JLG/BP JCP 321 (2016) 908-926
    IF (phimin.GE.0) THEN !===Two expansions
       p2 = pmin*((amin + amax - (ur-ul)*(gamma-1)/2)/(amin + amax*ratio))**(Mgas+2) !=== =p* (exact)
       rho_min = min(rhol*(p2/pl)*(1/gamma),rhor*(p2/pr)*(1/gamma))
       !TEST
       !rho_min = min((rhol*(p2/pl)*(1/gamma) + rhol)/2,(rhor*(p2/pr)*(1/gamma)+rhor)/2)
       !TEST
       rho_max = max(rhol,rhor)
       p_min = p2
       p_max = pmax
       lambdal = ul - al
       lambdar = ur + ar
       lambda_max = MAX(ABS(lambdal),ABS(lambdar))
       RETURN
    END IF
    !===We continue with phimin>0
    phimax = (pmax-pmin)*SQRT(capAmin/(pmax+capBmin)) + ur-ul !===(3.3) with (3.4) for p=pmax
    ptil = pmin*((amin + amax - (ur-ul)*(gamma-1)/2)/(amin + amax*ratio))**(Mgas+2)
    IF (phimax < 0.d0) THEN
       p1 = pmax
       p2 = ptil
    ELSE
       p1 = pmin
       p2 = MIN(pmax,ptil)
    END IF
    p1 = MAX(p1,p2-phi(p2,ul,pl,ur,pr)/phi_prime(p2,pl,pr))
    !===Assign min and max
    rho_min = min(rhol,rhor)
    rho_max = max(rhol,rhor)
    !===First wave
    !===Phi(pmin)<0: shock
    rho_max = max(rho_max,rhomin*(gm1_over_gp1 + p1/pmin)/(gm1_over_gp1*(p1/pmin) +1.d0))  !===(4.19) Toro, p.122
    !===Second wave
    IF (phimax < 0.d0) THEN
       !===Phi(pmax)<0: shock
       rho_max = max(rho_max,rhomax*(gm1_over_gp1 + p1/pmax)/(gm1_over_gp1*(p1/pmax) +1.d0))
       p_min = pmin
       p_max = p1
    ELSE
       !===Phi(pmax)>0: expansion
       rho_min = min(rho_min,rhomax*(p2/pmax)*(1/gamma))
       p_min = pmin
       p_max = pmax
    END IF

    lambdal = ul - al*SQRT(1.d0+MAX((p2-pl)/pl,0.d0)*one_m_exponent)
    lambdar = ur + ar*SQRT(1.d0+MAX((p2-pr)/pr,0.d0)*one_m_exponent)
    lambda_max = MAX(ABS(lambdal),ABS(lambdar))
  CONTAINS

    FUNCTION phi(p,ul,pl,ur,pr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: p, ul, pl, ur, pr
      REAL(KIND=8)             :: vv, fl, fr
      IF (p>pl) THEN
         fl = (p-pl)*SQRT(capAl/(p+capBl))
      ELSE
         fl = (2*al/(gamma-1))*((p/pl)**exponent-1)
      END IF
      IF (p>pr) THEN
         fr = (p-pr)*SQRT(capAr/(p+capBr))
      ELSE
         fr = (2*ar/(gamma-1))*((p/pr)**exponent-1)
      END IF
      vv = fl + fr + ur - ul
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
      vv = fl + fr
    END FUNCTION phi_prime
  END SUBROUTINE ptilde
END MODULE lambda_module

