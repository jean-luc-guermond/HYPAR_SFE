PROGRAM test
  USE fft_1d
  IMPLICIT NONE
  REAL(KIND=8), PARAMETER ::pi=ACOS(-1.d0)
  INTEGER, PARAMETER :: Nmaxc = 20
  REAL(KIND=8), DIMENSION(2*Nmaxc-1) :: r_u, r_out
  REAL(KIND=8),  DIMENSION(Nmaxc,2)  :: cs_u
  REAL(KIND=8) :: theta
  INTEGER :: i
  
  DO i = 1, SIZE(r_u)
     theta =  2*pi*(i-1)/(2*Nmaxc-1)
     r_u(i)  = 3 + COS(5*theta) + 0.5*SIN(19*theta) 
  END DO
  CALL real_to_fourier(r_u,cs_u)

  DO i = 1, Nmaxc 
     WRITE(*,*) cs_u(i,1), cs_u(i,2)
  END DO

  CALL fourier_to_real(cs_u,r_out)

  WRITE(*,*) 'ERROR', SUM(ABS(r_u-r_out))/SUM(ABS(r_u))
  
END PROGRAM test
