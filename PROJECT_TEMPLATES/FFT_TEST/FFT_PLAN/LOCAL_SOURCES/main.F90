PROGRAM test_fourier_1D_plan
   USE fft_1d
   USE fft_1d_plan
   USE fourier_param_module
   ! USE fourier_param_plan_module
#include "petsc/finclude/petsc.h"
   USE petsc
   USE read_inputs_module, ONLY : clean_data_once
   USE my_util, ONLY : error_petsc, itoa
   USE def_type_mesh
   USE st_matrix
   USE construct_mesh

   IMPLICIT NONE
   TYPE(fourier_param_type)  :: fourier_param
   ! TYPE(fourier_param_plan_type)  :: fourier_param_plan
   REAL(KIND=8), PARAMETER ::pi=ACOS(-1.d0)!, Length=2.d0
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: ar, ai
   REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: kvec
   REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: r_u, r_out, r1, r2
   REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: cs_u
   REAL(KIND=8) :: theta, h,  x, r_tol, err
   INTEGER :: i, k, rank, size_glob, dim_1
   LOGICAL :: test_passed

   TYPE(mesh_type)    :: mesh
   TYPE(petsc_csr_LA) :: LA

   MPI_Comm       :: communicator
   PetscErrorCode :: ierr

!===Start PETSC and MPI (mandatory)
   CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   communicator = PETSC_COMM_WORLD
   CALL MPI_Comm_rank(communicator, rank, ierr)
   CALL MPI_Comm_size(communicator, size_glob, ierr)

!===Clean data once
   CALL clean_data_once

!===Construct mesh & LA
   CALL get_mesh(communicator, mesh)
   CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA)

!=== Init Fourier setup
   CALL fourier_param%init
   ! CALL fourier_param_plan%init

   ! IF (MOD(fourier_param%Nmax, fourier_param_plan%size_F) /= 0) THEN
   !    CALL error_petsc("BUG in fourier param => incompatible MF="//itoa(fourier_param%Nmax)//&
   !    " and size_F="//itoa(fourier_param_plan%size_F))
   ! END IF
   ! IF (MOD(size_glob, fourier_param_plan%size_F) /= 0) THEN
   !    CALL error_petsc("BUG in fourier param => incompatible ntasks="//itoa(size_glob)//&
   !    " and size_F="//itoa(fourier_param_plan%size_F))
   ! END IF

   ! dim_1 = mesh%np
   dim_1 = 1000!fourier_param_plan%dim_1

   ALLOCATE(ar(dim_1,fourier_param%Nmax), ai(dim_1,fourier_param%Nmax), source=0.d0)
   ALLOCATE(kvec(fourier_param%Nmax), source=0.d0)
   ALLOCATE(r_u(dim_1,fourier_param%Nmax_real), r_out(dim_1,fourier_param%Nmax_real), source=0.d0)
   ALLOCATE(cs_u(dim_1,fourier_param%Nmax,2), source=0.d0)

!=== TEST OF FFT AND INVERSE FFT
   test_passed = .TRUE.
   r_tol = 1.d-14
  
   CALL RANDOM_NUMBER(ar)
   CALL RANDOM_NUMBER(ai)
   ai(:,1) = 0.d0
   kvec = [(k-1.d0, k=1, fourier_param%Nmax)]

   DO i = 1, SIZE(r_u,2)
      theta =  2*pi*(i-1)/(fourier_param%Nmax_real)
      DO k=1,fourier_param%Nmax
         r_u(:,i) = r_u(:,i) + ar(:,k)*COS(kvec(k)*theta) + ai(:,k)*SIN(kvec(k)*theta)
      END DO
   END DO

   CALL real_to_fourier_plan(r_u, cs_u)
   CALL fourier_to_real_plan(cs_u, r_out)


   err = SUM(ABS(r_u-r_out))/SUM(ABS(r_u))
   IF (err > r_tol) THEN
      test_passed = .FALSE.
      WRITE(*,*) 'ERROR FFT-IFFT', err, r_tol
   ELSE
      WRITE(*,*) "test passed with err=", err, r_tol
   END IF
!=== TEST OF DERIVATIVES

   h = fourier_param%Length/(fourier_param%Nmax_real)
   r_u = 0.d0
   DO i = 1, SIZE(r_u, 2)
      x = (i-1)*h
      theta =  2*pi*x/fourier_param%Length
      DO k = 1, fourier_param%Nmax 
         r_u(:,i) = r_u(:,i) + ar(:,k)*COS(kvec(k)*theta) + ai(:,k)*SIN(kvec(k)*theta)
      END DO
   END DO

   CALL real_derivative_plan(r_u,r_out,fourier_param%Length)

   r_u = 0.d0
   DO i = 1, SIZE(r_u, 2)
      x = (i-1)*h
      theta =  2*pi*x/fourier_param%Length
      DO k = 1, fourier_param%Nmax 
         r_u(:,i) = r_u(:,i) -ar(:,k)*(k-1)*SIN((k-1)*theta) + ai(:,k)*(k-1)*COS((k-1)*theta)
      END DO
   END DO

   r_u = 2*pi*r_u/fourier_param%Length
   err = SUM(ABS(r_u-r_out))/SUM(ABS(r_u))

   IF (err > r_tol) THEN
      test_passed = .FALSE.
      WRITE(*,*) 'ERROR on Derivative', err, r_tol
      ! DO i = 1, SIZE(r_u)
      !    write(*,*) r_u(i), r_out(i)
      ! END DO
   ELSE
      WRITE(*,*) 'SUCCESS on Derivative', err, r_tol
   END IF

   IF (test_passed) THEN
      WRITE(*,*) 'test fft_1D SUCCESS'
      WRITE(*,*) '1234567891'
   ELSE
      WRITE(*,*) 'test fft_1D FAILURE'
   END IF

   CALL error_petsc("end program")

END PROGRAM test_fourier_1D_plan
