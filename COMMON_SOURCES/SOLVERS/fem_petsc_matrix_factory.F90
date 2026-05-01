MODULE fem_petsc_matrix_factory_module
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
CONTAINS
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

   SUBROUTINE construct_lumped_mass_vector(mesh, LA, mass, lumped_mass)
      USE st_matrix, ONLY : create_my_ghost
      USE compute_periodic, ONLY : periodic_vector_petsc
      IMPLICIT NONE
      TYPE(mesh_type), INTENT(IN) :: mesh
      type(petsc_csr_LA), INTENT(IN) :: LA
      Mat, INTENT(IN) :: mass
      Vec :: vec_one, lumped_mass
      INTEGER, POINTER, DIMENSION(:) :: ifrom  ! for ghost structure
      INTEGER :: ierr

      !===Create ghost structure
      CALL create_my_ghost(mesh, LA, ifrom)
      CALL VecCreateGhost(PETSC_COMM_WORLD, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, lumped_mass, ierr)
      CALL VecDuplicate(lumped_mass, vec_one, ierr)

      CALL VecSet(vec_one, 1.d0, ierr)
      CALL MatMult(mass, vec_one, lumped_mass, ierr)
      CALL periodic_vector_petsc(mesh%per%nb_bords, mesh%per%list, mesh%per%perlist, lumped_mass, LA)
   END SUBROUTINE construct_lumped_mass_vector


   SUBROUTINE construct_cij(mesh, LA, cij)
      USE space_dim
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
                  vv_rowise(l) = SUM(mesh%gauss%dw(k, nj, :, m) * mesh%gauss%ww(ni, :) * mesh%gauss%rj(:, m))
               ENDDO
            ENDDO
            CALL MatSetValues(cij(k), mesh%gauss%n_w, idx, mesh%gauss%n_w, idx, vv_rowise, ADD_VALUES, ierr)
         ENDDO
         CALL MatAssemblyBegin(cij(k), MAT_FINAL_ASSEMBLY, ierr)
         CALL MatAssemblyEnd  (cij(k), MAT_FINAL_ASSEMBLY, ierr)
      END DO


   END SUBROUTINE construct_cij
END MODULE fem_petsc_matrix_factory_module
