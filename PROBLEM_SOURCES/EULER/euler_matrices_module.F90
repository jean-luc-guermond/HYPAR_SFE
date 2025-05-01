MODULE euler_matrices_module
#include "petsc/finclude/petsc.h"
   USE space_dim
   USE petsc
   USE def_type_mesh
   USE space_dim
   USE solver_petsc

   TYPE euler_matrices_type
      Mat :: mass, dij
      Mat, DIMENSION(k_dim) :: cij
      Mat, DIMENSION(k_dim) :: cij_loc
      REAL(KIND = 8), DIMENSION(:), POINTER :: lumped_mass
   CONTAINS
      PROCEDURE, PUBLIC :: construct => construct_euler_matrices
      PROCEDURE, PUBLIC :: compute_dij
   END TYPE euler_matrices_type

CONTAINS

   SUBROUTINE construct_euler_matrices(this, communicator, mesh, LA)
      USE fem_M
      CLASS(euler_matrices_type) :: this
      TYPE(mesh_type), INTENT(IN) :: mesh
      type(petsc_csr_LA), INTENT(IN) :: LA
      INTEGER :: k, ierr
      MPI_Comm       :: communicator

      !===Mat allocations
      CALL create_local_petsc_matrix(communicator, LA, this%mass, clean = .FALSE.)
      CALL create_local_petsc_matrix(communicator, LA, this%dij, clean = .FALSE.)
      DO k = 1, k_dim
         CALL create_local_petsc_matrix(communicator, LA, this%cij(k), clean = .FALSE.)
      END DO
      ALLOCATE(this%lumped_mass(mesh%np))
      write(*, *) 'alloc mat ok'
      !===Mat construction
      CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, this%mass)
      CALL construct_lumped_mass(mesh, LA, this%mass, this%lumped_mass)
      CALL construct_cij(mesh, LA, this%cij)
      write(*, *) 'const base mat ok'
      DO k = 1, k_dim
         CALL MatCreateSubMatrices(this%cij(k), 0, LA%loc_to_glob(1, :) - 1, &
              LA%loc_to_glob(1, :) - 1, MAT_INITIAL_MATRIX, this%cij_loc(k), ierr)
         write(*, *) 'a', k
      END DO

   END SUBROUTINE construct_euler_matrices

   SUBROUTINE construct_lumped_mass(mesh, LA, mass, lumped_mass)
      USE st_matrix
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

   SUBROUTINE compute_dij(this, mesh, LA)
      USE space_dim
      USE petsc
      USE my_util
      USE def_type_mesh
      CLASS(euler_matrices_type) :: this
      TYPE(mesh_type), INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), INTENT(IN) :: LA
      REAL(KIND = 8), DIMENSION(:), POINTER :: lumped_mass
      INTEGER :: m, ni, nj, e, nw, n, i, j, k, ierr

      nw = mesh%gauss%n_w

      DO m = 1, mesh%me
         IF (MINVAL(mesh%jj(:, m))>mesh%dom_np) CALL error_petsc('Cell with no vertices own by processor. Fix mesh distribution.')

         DO n = 1, nw
            IF (mesh%neigh(n, m) < m .AND. mesh%neigh(n, m) > 0) CYCLE !==cycle if neighbour is a cell already done

            ni = MOD(n, nw) + 1
            nj = MOD(n + 1, nw) + 1
            i = mesh%jj(ni, m)
            j = mesh%jj(nj, m)

            k = MIN(i, j)
            j = MAX(i, j)
            i = k

            IF (i > mesh%dom_np) CYCLE !===Verify that at least one point belong to processor

         END DO

      END DO

   END SUBROUTINE compute_dij

END MODULE euler_matrices_module