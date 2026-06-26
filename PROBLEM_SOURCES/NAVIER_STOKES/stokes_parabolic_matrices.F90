MODULE stokes_parabolic_matrices_module
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE solver_petsc
   USE periodic_data_module
   USE compute_periodic
   USE fem_petsc_matrix_factory_module, &
              ONLY : construct_lumped_mass_vector, construct_cij
TYPE stokes_parabolic_matrices_type
      INTEGER                          :: method, which_mass
      REAL(KIND = 8)                   :: thermal_diffusivity
      REAL(KIND = 8)                   :: mu_viscosity
      REAL(KIND = 8)                   :: lambda_viscosity
      REAL(KIND = 8)                   :: cv
      REAL(KIND=8), DIMENSION(:), POINTER :: scal_lumped_mass
      TYPE(petsc_csr_LA), POINTER      :: LA_vel, LA_temp
      TYPE(solver_data_type)           :: elasticity_solver_param, temperature_solver_param
      Mat :: vel_mass_mat,  vel_diff_mat,  vel_mat,  precond_vel_mat
      Mat :: temp_mass_mat, temp_diff_mat, temp_mat, precond_temp_mat
      Vec :: temp_lumped_mass_vec
      KSP :: vel_ksp, temp_ksp 
   CONTAINS
      PROCEDURE, PUBLIC :: construct => construct_stokes_parabolic_matrices
      PROCEDURE, PUBLIC :: var_mass_M
      PROCEDURE, PUBLIC :: elasticity_M, mass_vel_M

   END TYPE stokes_parabolic_matrices_type

CONTAINS

   SUBROUTINE construct_stokes_parabolic_matrices(this, communicator, mesh, LA_vel, LA_temp)
      USE space_dim
      USE fem_M
      USE st_matrix, ONLY: create_my_ghost, extract_through_ghost
      IMPLICIT NONE
      CLASS(stokes_parabolic_matrices_type) :: this
      TYPE(mesh_type),    INTENT(IN) :: mesh
      type(petsc_csr_LA), TARGET, INTENT(IN) :: LA_vel, LA_temp
      INTEGER                        :: k, ierr
      INTEGER, DIMENSION(:), POINTER :: ifrom
      MPI_Comm                       :: communicator
      
      this%LA_temp => LA_temp
      this%LA_vel => LA_vel

   !===Mat allocations temperature (temp_diff_mat)
      CALL create_local_petsc_matrix(communicator, LA_temp, this%temp_mass_mat, clean = .FALSE.)
      CALL MatDuplicate(this%temp_mass_mat, MAT_DO_NOT_COPY_VALUES, this%temp_diff_mat, ierr)
      CALL MatDuplicate(this%temp_mass_mat, MAT_DO_NOT_COPY_VALUES, this%temp_mat, ierr)

      !===Temperature diffusion Mat construction (temp_diff_mat)
      CALL qs_mass_diff_M (mesh, 0.d0, this%thermal_diffusivity, LA_temp, this%temp_diff_mat)! <=== construct this%temp_diff_mat
      ! CALL periodic_matrix_petsc(mesh%per, LA_temp, this%temp_diff_mat)
      !===Temperature diffusion Mat construction

      !===Temperature mass construction (temp_mass_mat)
      CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA_temp, this%temp_mass_mat) !<=== construct this%temp_mass_mat
      ! CALL periodic_matrix_petsc(mesh%per, LA_temp, this%temp_mass_mat)
      !===Temperature mass construction

      !===Temperature preconditioner
      CALL MatDuplicate(this%temp_mass_mat, MAT_SHARE_NONZERO_PATTERN, this%precond_temp_mat, ierr)
      CALL MatSetOption (this%precond_temp_mat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      !===Temperature preconditioner

      !===Temperature lumped mass construction (temp_lumped_mass_vec)
      CALL create_my_ghost(mesh, LA_temp, ifrom)
      CALL VecCreateGhost(communicator, mesh%dom_np, &
           PETSC_DETERMINE, SIZE(ifrom), ifrom, this%temp_lumped_mass_vec, ierr)
      CALL construct_lumped_mass_vector(mesh, LA_temp, this%temp_mass_mat, this%temp_lumped_mass_vec, opt_per=.FALSE.)
      !===Temperature lumped mass construction
   !===Mat allocations temperature (temp_diff_mat)

   !===Mat allocations & construction Velocity (vel_diff_mat)
      !===vel_dif_mat
      CALL create_local_petsc_matrix(communicator, LA_vel, this%vel_diff_mat, clean = .FALSE.)
      CALL MatSetOption (this%vel_diff_mat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      CALL this%elasticity_M(mesh) !<=== construct this%vel_diff_mat

      !===vel_mass_mat
      CALL MatDuplicate(this%vel_diff_mat, MAT_SHARE_NONZERO_PATTERN, this%vel_mass_mat, ierr)
      !CALL MatSetOption (this%vel_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
      CALL MatSetOption (this%vel_mass_mat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      CALL this%mass_vel_M (mesh) !<=== construct this%vel_mass_mat

      !===vel_mass_mat
      CALL MatDuplicate(this%vel_diff_mat, MAT_SHARE_NONZERO_PATTERN, this%vel_mat, ierr)
      !CALL MatSetOption (this%vel_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
      CALL MatSetOption (this%vel_mat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

      !===precond_vel_mat
      CALL MatDuplicate(this%vel_diff_mat, MAT_SHARE_NONZERO_PATTERN, this%precond_vel_mat, ierr)
      !CALL MatSetOption (this%precond_vel_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
      CALL MatSetOption (this%precond_vel_mat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)


      !===Init KSP vel
      CALL this%elasticity_solver_param%init('elasticity')
      CALL MatCopy(this%vel_diff_mat, this%vel_mat, SAME_NONZERO_PATTERN, ierr)
      CALL MatCopy(this%vel_diff_mat, this%precond_vel_mat, SAME_NONZERO_PATTERN, ierr)
      CALL init_solver(communicator, this%elasticity_solver_param, this%vel_ksp, this%vel_mat, opt_mat_pre=this%precond_vel_mat) 
      !===Init KSP vel
   

      !scal_lumped_mass
      ALLOCATE(this%scal_lumped_mass(mesh%np))
      CALL extract_through_ghost(this%temp_lumped_mass_vec, 1, 1, LA_temp, this%scal_lumped_mass, opt_assemble=.FALSE.)
   !===End Mat Velocity (vel_diff_mat)

   END SUBROUTINE construct_stokes_parabolic_matrices

   SUBROUTINE var_mass_M (this, mesh, mass, rho)
      !=================================================
      IMPLICIT NONE
      CLASS(stokes_parabolic_matrices_type) :: this
      TYPE(mesh_type)    :: mesh
      REAL(KIND = 8),               INTENT(IN) :: mass
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w, mesh%gauss%n_w) :: mat_loc

      INTEGER, DIMENSION(mesh%gauss%n_w) :: idxn
      INTEGER :: m, ni, nj, k, l
      REAL(KIND = 8), DIMENSION(mesh%gauss%l_G) :: bl
      PetscErrorCode :: ierr

      CALL MatZeroEntries (this%vel_mass_mat, ierr)
      DO m = 1, mesh%me
         idxn = this%LA_temp%loc_to_glob(1, mesh%jj(:, m)) - 1
         DO l = 1, mesh%gauss%l_G
            bl(l) =  SUM(rho(mesh%jj(:,m))*mesh%gauss%ww(:,l)) * mass * mesh%gauss%rj(l, m)
         END DO
         DO nj = 1, mesh%gauss%n_w;
            DO ni = 1, mesh%gauss%n_w;
               mat_loc(nj, ni) = SUM(mesh%gauss%ww(ni, :) * mesh%gauss%ww(nj, :) * bl)
            ENDDO
         ENDDO
         CALL MatSetValues(this%vel_mass_mat, mesh%gauss%n_w, idxn, mesh%gauss%n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO

      CALL MatAssemblyBegin(this%vel_mass_mat, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(this%vel_mass_mat, MAT_FINAL_ASSEMBLY, ierr)
   END SUBROUTINE var_mass_M

   SUBROUTINE elasticity_M (this, mesh)
      !=================================================
      USE space_dim
      IMPLICIT NONE
      CLASS(stokes_parabolic_matrices_type) :: this
      TYPE(mesh_type)                       :: mesh
      INTEGER, DIMENSION(mesh%gauss%n_w) :: jj_loc
      INTEGER, DIMENSION(k_dim*mesh%gauss%n_w) :: idxm, idxn
      INTEGER :: m, mi, i, ki, iglob, ix, nj, j, kj, jx, l, n_w, ni, k1, jglob
      REAL(KIND=8) :: x, y, lambda
      REAL(KIND = 8), DIMENSION(k_dim*mesh%gauss%n_w, k_dim*mesh%gauss%n_w) :: mat_loc
      PetscErrorCode                     :: ierr


      lambda = this%lambda_viscosity-2.d0/3.d0*this%mu_viscosity

      CALL MatZeroEntries (this%vel_diff_mat, ierr)
      CALL MatSetOption (this%vel_diff_mat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      CALL MatSetOption (this%vel_diff_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
      n_w = mesh%gauss%n_w
      DO m=1, mesh%me
         jj_loc = mesh%jj(:, m)
         mat_loc = 0.d0
         DO ni = 1, n_w
            i = jj_loc(ni)
            DO ki = 1, k_dim
               iglob = this%LA_vel%loc_to_glob(ki, i)
               ix = (ki - 1) * n_w + ni
               idxm(ix) = iglob - 1
               DO nj = 1, n_w
                  j = jj_loc(nj)
                  DO kj = 1, k_dim
                     jglob = this%LA_vel%loc_to_glob(kj, j)
                     jx = (kj - 1) * n_w + nj
                     idxn(jx) = jglob - 1
                     x = 0
                     DO l = 1, mesh%gauss%l_G
                        y =  this%mu_viscosity*mesh%gauss%dw(kj,ni,l,m)*mesh%gauss%dw(ki,nj,l,m) &
                       + lambda*mesh%gauss%dw(ki,ni,l,m)*mesh%gauss%dw(kj,nj,l,m)
                        IF (kj.EQ.ki) THEN
                           DO k1 = 1, k_dim
                              y = y + this%mu_viscosity*mesh%gauss%dw(k1,ni,l,m)*mesh%gauss%dw(k1,nj,l,m)
                           END DO
                        END IF
                        x = x + y * mesh%gauss%rj(l,m)
                     END DO
                     mat_loc(ix,jx) = x
                  END DO
               END DO
            END DO
         END DO
         CALL MatSetValues(this%vel_diff_mat, k_dim * n_w, idxm, k_dim * n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO

      CALL MatAssemblyBegin(this%vel_diff_mat, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(this%vel_diff_mat, MAT_FINAL_ASSEMBLY, ierr)

   END SUBROUTINE elasticity_M

   SUBROUTINE mass_vel_M (this, mesh)
      !=================================================
      USE space_dim
      IMPLICIT NONE
      CLASS(stokes_parabolic_matrices_type) :: this
      TYPE(mesh_type)                       :: mesh
      INTEGER, DIMENSION(mesh%gauss%n_w) :: jj_loc
      INTEGER, DIMENSION(k_dim*mesh%gauss%n_w) :: idxm, idxn
      INTEGER :: m, mi, i, ki, iglob, ix, nj, j, kj, jx, l, n_w, ni, k1, jglob
      REAL(KIND=8) :: x, y
      REAL(KIND = 8), DIMENSION(k_dim*mesh%gauss%n_w, k_dim*mesh%gauss%n_w) :: mat_loc
      PetscErrorCode                     :: ierr

      CALL MatZeroEntries (this%vel_mass_mat, ierr)
      CALL MatSetOption (this%vel_mass_mat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      CALL MatSetOption (this%vel_mass_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
      n_w = mesh%gauss%n_w
      DO m=1, mesh%me
         jj_loc = mesh%jj(:, m)
         mat_loc = 0.d0
         DO ni = 1, n_w
            i = jj_loc(ni)
            DO ki = 1, k_dim
               iglob = this%LA_vel%loc_to_glob(ki, i)
               ix = (ki - 1) * n_w + ni
               idxm(ix) = iglob - 1
               DO nj = 1, n_w
                  j = jj_loc(nj)
                  kj = ki
                  jglob = this%LA_vel%loc_to_glob(kj, j)
                  jx = (kj - 1) * n_w + nj
                  idxn(jx) = jglob - 1
                  x = SUM(mesh%gauss%ww(ni,:)*mesh%gauss%ww(nj,:)*mesh%gauss%rj(:,m))
                  mat_loc(ix,jx) = x
               END DO
            END DO
         END DO
         CALL MatSetValues(this%vel_mass_mat, k_dim * n_w, idxm, k_dim * n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO

      CALL MatAssemblyBegin(this%vel_mass_mat, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(this%vel_mass_mat, MAT_FINAL_ASSEMBLY, ierr)

   END SUBROUTINE mass_vel_M



END MODULE stokes_parabolic_matrices_module