MODULE euler_matrices_module
#include "petsc/finclude/petsc.h"
   USE space_dim
   USE petsc
   USE def_type_mesh
   USE space_dim
   USE solver_petsc

   TYPE euler_matrices_type
      Mat :: mass, dij
      Mat, DIMENSION(k_dim) :: cij, nij_loc
      Mat, DIMENSION(1, k_dim) :: cij_loc
      Mat :: cij_norm_loc
      REAL(KIND = 8), DIMENSION(:), POINTER :: lumped_mass
   CONTAINS
      PROCEDURE, PUBLIC :: construct => construct_euler_matrices
      PROCEDURE, PRIVATE :: construct_loc_nij
   END TYPE euler_matrices_type

CONTAINS

   SUBROUTINE construct_euler_matrices(this, communicator, mesh, LA)
      USE fem_M
      IMPLICIT NONE
      CLASS(euler_matrices_type) :: this
      TYPE(mesh_type), INTENT(IN) :: mesh
      type(petsc_csr_LA), INTENT(IN) :: LA
      INTEGER :: k, ierr
      MPI_Comm       :: communicator
      IS, DIMENSION(1) :: is

      !===Mat allocations
      CALL create_local_petsc_matrix(communicator, LA, this%mass, clean = .FALSE.)
      CALL create_local_petsc_matrix(communicator, LA, this%dij, clean = .FALSE.)
      DO k = 1, k_dim
         CALL create_local_petsc_matrix(communicator, LA, this%cij(k), clean = .FALSE.)
      END DO
      ALLOCATE(this%lumped_mass(mesh%np))

      !===Mat construction
      CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, this%mass)
      CALL construct_lumped_mass(mesh, LA, this%mass, this%lumped_mass)
      CALL construct_cij(mesh, LA, this%cij)

      CALL ISCreateGeneral(communicator, mesh%np, LA%loc_to_glob(1, :) - 1, PETSC_COPY_VALUES, is(1), ierr)
      DO k = 1, k_dim
         CALL MatCreateSubMatrices(this%cij(k), 1, is, is, MAT_INITIAL_MATRIX, this%cij_loc(:, k), ierr)
         CALL MatDuplicate(this%cij_loc(1, k), MAT_DO_NOT_COPY_VALUES, this%nij_loc(k), ierr)
      END DO
      CALL MatDuplicate(this%cij_loc(1, 1), MAT_DO_NOT_COPY_VALUES, this%cij_norm_loc, ierr)

      CALL this%construct_loc_nij(mesh)

   END SUBROUTINE construct_euler_matrices

   SUBROUTINE construct_lumped_mass(mesh, LA, mass, lumped_mass)
      USE st_matrix
      IMPLICIT NONE
      TYPE(mesh_type), INTENT(IN) :: mesh
      type(petsc_csr_LA), INTENT(IN) :: LA
      Mat, INTENT(IN) :: mass
      REAL(KIND = 8), DIMENSION(:), POINTER :: lumped_mass
      Vec :: vec_one, xx, x_ghost
      INTEGER, POINTER, DIMENSION(:) :: ifrom  ! for ghost structure
      INTEGER :: ierr

      !===Create ghost structure
      CALL create_my_ghost(mesh, LA, ifrom)
      CALL VecCreateGhost(PETSC_COMM_WORLD, mesh%dom_np, &
           PETSC_DETERMINE, SIZE(ifrom), ifrom, xx, ierr)
      CALL VecDuplicate(xx, vec_one, ierr)

      CALL VecSet(vec_one, 1.d0, ierr)
      CALL MatMult(mass, vec_one, xx, ierr)
      CALL VecGhostGetLocalForm(xx, x_ghost, ierr)
      CALL VecGhostUpdateBegin(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)
      CALL VecGhostUpdateEnd(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)
      CALL extract(x_ghost, 1, 1, LA, lumped_mass)

   END SUBROUTINE construct_lumped_mass


   SUBROUTINE construct_cij(mesh, LA, cij)
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type), INTENT(IN) :: mesh
      type(petsc_csr_LA), INTENT(IN) :: LA
      Mat, DIMENSION(:) :: cij
      REAL(KIND = 8), DIMENSION(mesh%gauss%n_w * mesh%gauss%n_w) :: vv_rowise
      INTEGER, DIMENSION(mesh%gauss%n_w) :: idx
      INTEGER :: m, ni, nj, l, k, ierr

      DO k = 1, k_dim
         CALL MatZeroEntries (cij(k), ierr)
         DO m = 1, mesh%dom_me
            idx = LA%loc_to_glob(1, mesh%jj(:, m)) - 1
            l = 0
            DO ni = 1, mesh%gauss%n_w
               DO nj = 1, mesh%gauss%n_w
                  l = l + 1
                  vv_rowise(l) = -SUM(mesh%gauss%dw(k, nj, :, m) * mesh%gauss%ww(ni, :) * mesh%gauss%rj(:, m))
               ENDDO
            ENDDO
            CALL MatSetValues(cij(k), mesh%gauss%n_w, idx, mesh%gauss%n_w, idx, vv_rowise, ADD_VALUES, ierr)
         ENDDO
         CALL MatAssemblyBegin(cij(k), MAT_FINAL_ASSEMBLY, ierr)
         CALL MatAssemblyEnd  (cij(k), MAT_FINAL_ASSEMBLY, ierr)
      END DO

   END SUBROUTINE construct_cij

   SUBROUTINE construct_loc_nij(this, mesh)
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(euler_matrices_type) :: this
      TYPE(mesh_type), INTENT(IN) :: mesh
      REAL(KIND = 8), DIMENSION(1, k_dim) :: cij_c
      REAL(KIND = 8), DIMENSION(1, 1) :: norm, nij_c
      INTEGER, DIMENSION(1) :: i, j
      LOGICAL, DIMENSION(mesh%medge) :: virgin_edge = .true.
      INTEGER :: k, m, n, ni, nj, ierr

      DO m = 1, mesh%me
         DO n = 1, mesh%n_e
            IF (mesh%attrib_e(mesh%jce(n, m))) THEN
               IF (.not. virgin_edge(mesh%jce_loc(n, m))) CYCLE

               ni = MOD(n, nw) + 1
               nj = MOD(n + 1, nw) + 1
               i = mesh%jj(ni, m)
               j = mesh%jj(nj, m)

               norm = 0.d0
               DO k = 1, k_dim
                  CALL MatGetValues(this%cij_loc(1, k), 1, i, 1, j, cij_c(:, k), ierr)
                  norm = norm + cij_c(1, k)**2
               END DO
               norm = SQRT(norm)

               CALL MatSetValues(this%cij_norm_loc, 1, i, 1, j, norm, ADD_VALUES, ierr)

               DO k = 1, k_dim
                  nij_c = cij_c(1, k) / norm
                  CALL MatSetValues(this%nij_loc(k), 1, i, 1, j, nij_c, ADD_VALUES, ierr)
               END DO
            END IF
         END DO
      END DO

   END SUBROUTINE construct_loc_nij


END MODULE euler_matrices_module