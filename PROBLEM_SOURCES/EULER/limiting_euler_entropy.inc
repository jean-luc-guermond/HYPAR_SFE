
! FUNCTION zero_of_specific_entropy(cmin,u0,Pij) RESULT(v)
!      IMPLICIT NONE
!      REAL(KIND=8), DIMENSION(:) :: u0, Pij
!      REAL(KIND=8)   :: cmin, v
!      REAL(KIND=8), DIMENSION(SIZE(u0)) :: ul, ur
!      REAL(KIND=8) :: Esmall, psir, psil, ll, lr, llold, lrold

!      Esmall= small*u0(inputs%syst_size)
!      ur = u0 + Pij 

!      psir = entropy_min(ur,cmin)
!      IF (psir.GE.-Esmall) THEN
!           v = 1.d0
!           RETURN
!      END IF
!      ll = 0.d0
!      ul = u0
!      psil = entropy_min(ul,cmin)
!      DO WHILE (ABS(psil-psir) .GT. Esmall)
!           llold = ll
!           lrold = lr
!           ll = ll - psil*(lr-ll)/(psir-psil)
!           lr = lr - psir/psi_prime_func(Pij,ur,cmin)
!           IF (ll.GE.lr) THEN
!                ll = lr
!                EXIT
!           END IF
!           IF (ll< llold) THEN
!                ll = llold
!                EXIT
!           END IF
!           IF (lr > lrold) THEN
!                lr = lrold
!                EXIT
!           END IF
!           ul = u0 + ll*Pij
!           ur = u0 + lr*Pij
!           psil = entropy_min(ul,cmin)
!           psir = entropy_min(ur,cmin)
!      END DO
!      IF (psir.GE.-Esmall) THEN
!           v = lr
!      ELSE
!           v = ll
!      END IF
!      END FUNCTION zero_of_specific_entropy

!   FUNCTION entropy_min(u,psi_m) RESULT(psi)
!     USE boundary_conditions
!     IMPLICIT NONE
!     REAL(KIND=8), DIMENSION(:) :: u
!     REAL(KIND=8) :: psi_m, psi
!     INTEGER :: syst_size
!     syst_size = k_dim+2
!     !psi = u(:,syst_size) - SUM(u(:,2:syst_size-1)**2,dim=2)/(2.d0*u(:,1)) - psi_m*u(:,1)**gamma
!     psi = (u(:,syst_size) - SUM(u(:,2:syst_size-1)**2,dim=2)/(2.d0*u(:,1)))/u(:,1)**gamma - psi_m
!   END FUNCTION entropy_min
  
!   FUNCTION psi_prime_func(Pij,u,psi_m) RESULT(psi)
!     USE boundary_conditions
!     IMPLICIT NONE
!     REAL(KIND=8), DIMENSION(:) :: u, Pij
!     REAL(KIND=8)               :: psi_m, psi
!     INTEGER :: syst_size
!     syst_size = k_dim+2
!     psi = Pij(syst_size) - SUM(u(2:syst_size-1)*Pij(2:syst_size-1))/u(1) &
!          + Pij(1)*SUM(u(2:syst_size-1)**2)/(2*u(1)**2) &
!          - psi_m*gamma*Pij(1)*u(1)**(gamma-1.d0)
!   END FUNCTION psi_prime_func