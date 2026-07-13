#include "petsc/finclude/petsc.h"
MODULE stokes_parabolic_matrices_module
   !> Module to build matrices used for solving 
   !! the velocity and the temperature successively
   USE petsc
   USE def_type_mesh
   USE petsc_csr_LA_module

   USE solver_petsc
   USE periodic_data_module
   USE compute_periodic
   USE fem_petsc_matrix_factory_module, &
              ONLY : construct_lumped_mass_vector, construct_elasticity_M
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
   END TYPE stokes_parabolic_matrices_type

CONTAINS

   SUBROUTINE construct_stokes_parabolic_matrices(this, communicator, mesh, LA_vel, LA_temp)
      !> Construction of matrices for Stokes
      !! TODO: periodic conditions; consistent matrices; implicit lumped diffusion?
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

      !===Init KSP temp
      CALL this%temperature_solver_param%init('temperature')
      CALL init_solver(communicator, this%temperature_solver_param, this%temp_ksp, this%temp_mat, opt_mat_pre=this%precond_temp_mat) 
      !===Init KSP temp

   !===Mat allocations temperature (temp_diff_mat)

   !===Mat allocations & construction Velocity (vel_diff_mat)
      !===vel_dif_mat
      CALL create_local_petsc_matrix(communicator, LA_vel, this%vel_diff_mat, clean = .FALSE.)
      CALL MatSetOption (this%vel_diff_mat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      CALL construct_elasticity_M(mesh, this%LA_vel, this%vel_diff_mat, &
                                 this%lambda_viscosity, this%mu_viscosity) !<=== construct this%vel_diff_mat
      ! CALL this%elasticity_M(mesh) !<=== construct this%vel_diff_mat

      !===vel_mass_mat
      CALL MatDuplicate(this%vel_diff_mat, MAT_SHARE_NONZERO_PATTERN, this%vel_mass_mat, ierr)
      !CALL MatSetOption (this%vel_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
      CALL MatSetOption (this%vel_mass_mat, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      CALL qs_mass_block_M (mesh, 1.d0, this%LA_vel, this%vel_mass_mat)

      !===vel_mat
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

END MODULE stokes_parabolic_matrices_module