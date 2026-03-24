PROGRAM test
  USE fft_1d
  IMPLICIT NONE
  REAL(KIND=8), PARAMETER ::pi=ACOS(-1.d0), Length=2.d0
  INTEGER, PARAMETER :: Nmaxc = 20, k1=5, k2=19
  REAL(KIND=8), DIMENSION(Nmaxc) :: ar, ai
  REAL(KIND=8), DIMENSION(2*Nmaxc-1) :: r_u, r_out
  REAL(KIND=8),  DIMENSION(Nmaxc,2)  :: cs_u
  REAL(KIND=8) :: theta, h,  x
  INTEGER :: i, k

  CALL RANDOM_NUMBER(ar)
  CALL RANDOM_NUMBER(ai)
  ai(1) = 0
  DO i = 1, SIZE(r_u)
     theta =  2*pi*(i-1)/(2*Nmaxc-1)
     DO k = 1, Nmaxc 
        r_u(i)  = ar(k)*COS((k-1)*theta) + ai(k)*SIN((k-1)*theta)
     END DO
  END DO
  CALL real_to_fourier(r_u,cs_u)

  DO i = 1, Nmaxc 
     WRITE(*,*) cs_u(i,1), cs_u(i,2)
  END DO

  CALL fourier_to_real(cs_u,r_out)

  WRITE(*,*) 'ERROR FFT', SUM(ABS(r_u-r_out))/SUM(ABS(r_u))

  h = Length/(2*Nmaxc-1)
  DO i = 1, SIZE(r_u)
     x = (i-1)*h
     theta =  2*pi*x/Length
     DO k = 1, Nmaxc 
        r_u(i)  = ar(k)*COS((k-1)*theta) + ai(k)*SIN((k-1)*theta)
     END DO
  END DO

  CALL real_derivative(r_u,r_out,Length)

  DO i = 1, SIZE(r_u)
     x = (i-1)*h
     theta =  2*pi*x/Length
     DO k = 1, Nmaxc 
        r_u(i)  = -ar(k)*(k-1)*SIN((k-1)*theta) + ai(k)*(k-1)*COS((k-1)*theta)
     END DO
  END DO
  r_u = 2*pi*r_u/length

  WRITE(*,*) 'ERROR on Derivative', SUM(ABS(r_u-r_out))/SUM(ABS(r_u))
  DO i = 1, SIZE(r_u)
     write(*,*) r_u(i), r_out(i)
  END DO
END PROGRAM test
