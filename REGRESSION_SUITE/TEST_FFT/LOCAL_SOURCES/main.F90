PROGRAM test_fourier_1D
  USE fft_1d
  USE fourier_param_module
#include "petsc/finclude/petsc.h"
  USE petsc
  USE character_strings, ONLY : clean_data_once

  IMPLICIT NONE
  TYPE(fourier_param_type)  :: fourier_param
  REAL(KIND=8), PARAMETER ::pi=ACOS(-1.d0)!, Length=2.d0
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ar, ai
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: r_u, r_out, r1, r2
  REAL(KIND=8),  DIMENSION(:,:), ALLOCATABLE  :: cs_u
  REAL(KIND=8) :: theta, h,  x, r_tol, err
  INTEGER :: i, k, rank
  LOGICAL :: test_passed
  MPI_Comm       :: communicator
  PetscErrorCode :: ierr

!===Start PETSC and MPI (mandatory)
  CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  communicator = PETSC_COMM_WORLD
  CALL MPI_Comm_rank(communicator, rank, ierr)

!===Clean data once
  CALL clean_data_once

!=== Init Fourier setup
  CALL fourier_param%init
  ALLOCATE(ar(fourier_param%Nmax), ai(fourier_param%Nmax))
  ALLOCATE(r_u(fourier_param%Nmax_real), r_out(fourier_param%Nmax_real))
  ALLOCATE(cs_u(fourier_param%Nmax,2))

!=== TEST OF FFT AND INVERSE FFT
  test_passed = .TRUE.
  r_tol = 1.d-14
  
  
  CALL RANDOM_NUMBER(ar)
  CALL RANDOM_NUMBER(ai)
  ai(1) = 0
  DO i = 1, SIZE(r_u)
     theta =  2*pi*(i-1)/(fourier_param%Nmax_real)
     DO k = 1, fourier_param%Nmax 
        r_u(i)  = ar(k)*COS((k-1)*theta) + ai(k)*SIN((k-1)*theta)
     END DO
  END DO
  CALL real_to_fourier(r_u,cs_u)



  CALL fourier_to_real(cs_u,r_out)
  err = SUM(ABS(r_u-r_out))/SUM(ABS(r_u))
  IF (err > r_tol) THEN
     test_passed = .FALSE.
     WRITE(*,*) 'ERROR FFT-IFFT', err, r_tol
     DO i = 1, fourier_param%Nmax 
        WRITE(*,*) cs_u(i,1), cs_u(i,2)
     END DO
  END IF
!=== TEST OF FFT AND INVERSE FFT

!=== TEST OF DERIVATIVE
  h = fourier_param%Length/(fourier_param%Nmax_real)
  DO i = 1, SIZE(r_u)
     x = (i-1)*h
     theta =  2*pi*x/fourier_param%Length
     DO k = 1, fourier_param%Nmax 
        r_u(i)  = ar(k)*COS((k-1)*theta) + ai(k)*SIN((k-1)*theta)
     END DO
  END DO

  CALL real_derivative(r_u,r_out,fourier_param%Length)

  DO i = 1, SIZE(r_u)
     x = (i-1)*h
     theta =  2*pi*x/fourier_param%Length
     DO k = 1, fourier_param%Nmax 
        r_u(i)  = -ar(k)*(k-1)*SIN((k-1)*theta) + ai(k)*(k-1)*COS((k-1)*theta)
     END DO
  END DO
  r_u = 2*pi*r_u/fourier_param%Length
  err = SUM(ABS(r_u-r_out))/SUM(ABS(r_u))

  IF (err > r_tol) THEN
     test_passed = .FALSE.
     WRITE(*,*) 'ERROR on Derivative', err, r_tol
     DO i = 1, SIZE(r_u)
        write(*,*) r_u(i), r_out(i)
     END DO
  END IF

  IF (test_passed) THEN
     WRITE(*,*) 'test fft_1D SUCCESS'
     WRITE(*,*) '1234567891'
  ELSE
     WRITE(*,*) 'test fft_1D FAILURE'
  END IF


!=== TEST OF DERIVATIVE

!=== TEST OF PRODUCT (todo??)

!=== Finalize PETSC and MPI (mandatory)
  CALL PetscFinalize(ierr) 

END PROGRAM test_fourier_1D
