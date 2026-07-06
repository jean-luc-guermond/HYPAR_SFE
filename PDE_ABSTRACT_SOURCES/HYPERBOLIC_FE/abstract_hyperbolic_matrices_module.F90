MODULE hyperbolic_matrices_module
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE solver_petsc
   USE periodic_data_module
   USE compute_periodic
   USE fem_petsc_matrix_factory_module, &
              ONLY : construct_lumped_mass_vector, construct_cij
   TYPE hyperbolic_matrices_type
      INTEGER                          :: method, which_mass
      !=== START Objects eventually on low order stencil  ===!
      REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: cijL_norm_loc_array
      REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: nijL_loc_array
      Mat                              :: dijL, dijH, stiffL
      Mat, DIMENSION(:,:), ALLOCATABLE :: cij_loc         !=== temporary
      Mat, DIMENSION(:),   POINTER     :: cijL
      Vec,                 POINTER     :: lump_mass_vec_L !=== for dt
      !=== END   Objects eventually on low order stencil  ===!
      !=== START Objects eventually on high order stencil ===!
      Mat, DIMENSION(:), POINTER     :: cij
      Vec,               POINTER     :: lump_mass_vec
      Mat                            :: mass, Al_mass 
      KSP                            :: ksp_consistent_mass
      !=== END   Objects eventually on high order stencil ===!
   CONTAINS
      PROCEDURE, PUBLIC :: construct => construct_hyperbolic_matrices
      PROCEDURE, PRIVATE :: construct_loc_nij, construct_consistent_mass_solver, construct_Al
   END TYPE hyperbolic_matrices_type

   INTEGER, PRIVATE, PARAMETER :: METHOD_VISCOUS=1, METHOD_HIGH=2
   INTEGER, PRIVATE, PARAMETER :: LUMPED_MASS=1, QUASI_CONSISTENT_MASS=2, CONSISTENT_MASS=3

CONTAINS

   SUBROUTINE construct_hyperbolic_matrices(this, communicator, mesh, LA, mesh_L, LA_L)
      USE space_dim
      USE fem_M
      USE st_matrix, ONLY: create_my_ghost
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      TYPE(mesh_type),    INTENT(IN) :: mesh, mesh_L
      type(petsc_csr_LA), INTENT(IN) :: LA,   LA_L
      INTEGER                        :: k, ierr
      INTEGER, DIMENSION(:), POINTER :: ifrom
      MPI_Comm                       :: communicator
      IS,      DIMENSION(1)          :: is !<=== FIXME

      !========================================================!
      !=== START BY CONSTRUCTING LOW ORDER STENCIL MATRICES ===!
      !========================================================!
      
      !=== temporary, to be reassigned to high order stencil later
      CALL create_local_petsc_matrix(mesh%comm, LA_L, this%mass, clean = .FALSE.)
      !=== temporary, to be reassigned to high order stencil later

      ALLOCATE(this%cijL(k_dim))
      ALLOCATE(this%cij_loc(1, k_dim))
      ALLOCATE(this%lump_mass_vec_L)
      DO k = 1, k_dim
         CALL create_local_petsc_matrix(mesh%comm, LA_L, this%cijL(k), clean = .FALSE.)
      END DO
      
      !=== Construct lump_mass for low order stencil (for dt)
      CALL qs_mass_diff_M (mesh_L, 1.d0, 0.d0, LA_L, this%mass)
      CALL periodic_matrix_petsc(mesh_L%per, LA_L, this%mass)
      
      CALL create_my_ghost(mesh_L, LA_L, ifrom)
      CALL VecCreateGhost(mesh%comm, mesh_L%dom_np, &
           PETSC_DETERMINE, SIZE(ifrom), ifrom, this%lump_mass_vec_L, ierr)
      CALL construct_lumped_mass_vector(mesh_L, LA_L, this%mass, this%lump_mass_vec_L)
      !=== END Construct lump_mass for low order stencil (for dt)
      
      !=== Construct cij for low order stencil (for viscous method)
      CALL construct_cij(mesh_L, LA_L, this%cijL)
      CALL ISCreateGeneral(mesh%comm, mesh_L%np, LA_L%loc_to_glob(1, :) - 1, PETSC_COPY_VALUES, is(1), ierr)
      DO k = 1, k_dim
         CALL MatCreateSubMatrices(this%cijL(k), 1, is, is, MAT_INITIAL_MATRIX, this%cij_loc(:, k), ierr)
      END DO
      ALLOCATE(this%cijL_norm_loc_array(mesh_L%gauss%n_w, mesh_L%gauss%n_w, mesh_L%me))
      ALLOCATE(this%nijL_loc_array(mesh_L%gauss%n_w, mesh_L%gauss%n_w, k_dim, mesh_L%me))
      CALL this%construct_loc_nij(mesh_L)
      !=== END Construct cij for low order stencil (for viscous method)

      !=== Cleanup
      DO k=1, k_dim
         CALL MatDestroy(this%cij_loc(1, k), ierr)
      END DO
      DEALLOCATE(this%cij_loc)
      !=== End   Cleanup

      CALL create_local_petsc_matrix(mesh%comm, LA_L, this%dijL, clean = .FALSE.)
      IF (this%method == METHOD_HIGH) THEN
         CALL create_local_petsc_matrix(mesh%comm, LA_L, this%dijH, clean = .FALSE.)
         CALL create_local_petsc_matrix(mesh%comm, LA_L, this%stiffL, clean = .FALSE.)
         CALL qs_mass_diff_M (mesh_L, 0.d0, 1.d0, LA_L, this%stiffL)
      END IF

      !=================================================!
      !=== NOW CONSTRUCT HIGH ORDER STENCIL MATRICES ===!
      !=================================================!

      !=== Only create other matrices in case mesh_L and mesh are actually different ===!
      IF (mesh_L%me /= mesh%me) THEN
         ALLOCATE(this%lump_mass_vec)
         ALLOCATE(this%cij(k_dim))

         !=== Create masses for solver ===!
         CALL MatDestroy(this%mass, ierr)
         CALL create_local_petsc_matrix(mesh%comm, LA, this%mass, clean = .FALSE.)
         CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, this%mass)
         CALL periodic_matrix_petsc(mesh%per, LA, this%mass)
         
         CALL create_my_ghost(mesh, LA, ifrom)
         CALL VecCreateGhost(mesh%comm, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, this%lump_mass_vec, ierr)
         CALL construct_lumped_mass_vector(mesh, LA, this%mass, this%lump_mass_vec)
         !=== Create masses for solver ===!

         !=== Create cij for high order solver ===!
         DO k = 1, k_dim
            CALL create_local_petsc_matrix(mesh%comm, LA, this%cij(k), clean = .FALSE.)
         END DO
         CALL construct_cij(mesh, LA, this%cij)
         !=== Create cij for high order solver ===!

      !=== Build pointers in non-hybrid case ===!
      ELSE
         this%lump_mass_vec => this%lump_mass_vec_L
         this%cij           => this%cijL
      END IF

      !===============================!
      !=== Solvers for mass matrix ===!
      !===============================!

      IF (this%which_mass==CONSISTENT_MASS) THEN
         CALL this%construct_consistent_mass_solver(mesh%comm)
      ELSE IF (this%which_mass==QUASI_CONSISTENT_MASS) THEN
         CALL this%construct_Al
      END IF

   END SUBROUTINE construct_hyperbolic_matrices

   SUBROUTINE construct_loc_nij(this, mesh)
      USE space_dim
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      TYPE(mesh_type), INTENT(IN) :: mesh
      REAL(KIND = 8), DIMENSION(1, k_dim) :: cij_c
      REAL(KIND = 8), DIMENSION(1, 1) :: nij_c, norm
      REAL(KIND =8) :: xx
      INTEGER, DIMENSION(1) :: i, j
      LOGICAL, DIMENSION(mesh%medge) :: virgin_edge
      INTEGER :: k, m, n, ni, nj, ierr, nw, edge

      nw = mesh%gauss%n_w
      this%cijL_norm_loc_array = 0.d0
      this%nijL_loc_array = 0.d0
      DO m = 1, mesh%me
         DO ni = 1, nw
            DO nj=1, nw
               IF (ni==nj) CYCLE
               i = mesh%jj(ni, m)
               j = mesh%jj(nj, m)

               xx = 0.d0
               DO k = 1, k_dim
                  CALL MatGetValues(this%cij_loc(1, k), 1, i - 1, 1, j - 1, cij_c(:, k), ierr)
                  xx = xx + cij_c(1, k)**2
               END DO
               xx = SQRT(xx)
               this%cijL_norm_loc_array(ni, nj, m) = xx

               DO k = 1, k_dim
                  this%nijL_loc_array(ni, nj, k, m) = cij_c(1, k)/xx
               END DO

            END DO
         END DO
      END DO


   END SUBROUTINE construct_loc_nij

   SUBROUTINE construct_consistent_mass_solver(this, communicator)
      USE solver_petsc
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      TYPE(solver_data_type)          :: my_par_solver
      MPI_Comm                        :: communicator

      !===Create ksp solver
      CALL init_solver(communicator, my_par_solver, this%ksp_consistent_mass, this%mass)      

   END SUBROUTINE construct_consistent_mass_solver

   SUBROUTINE construct_Al(this)
      IMPLICIT NONE
      CLASS(hyperbolic_matrices_type) :: this
      INTEGER :: ierr
      Vec     :: xx
      CALL VecDuplicate(this%lump_mass_vec, xx, ierr)
      CALL VecSet(xx, 1.d0, ierr)
      CALL MatDuplicate(this%mass, MAT_COPY_VALUES, this%Al_mass, ierr)
      CALL VecPointWiseDivide(xx, xx, this%lump_mass_vec, ierr)
      CALL MatDiagonalScale(this%Al_mass, PETSC_NULL_VEC, xx, ierr)
      CALL MatScale(this%Al_mass, -1.d0, ierr)
      CALL MatShift(this%Al_mass, 1.d0, ierr)
      CALL VecDestroy(xx, ierr)
   END SUBROUTINE construct_Al

END MODULE hyperbolic_matrices_module
