#include "petsc/finclude/petsc.h"
MODULE fem_petsc_matrix_factory_module
   USE petsc
   USE def_type_mesh
   USE petsc_csr_LA_module
   PUBLIC :: construct_lumped_mass, construct_lumped_mass_vector, construct_cij, &
             construct_elasticity_M
   PRIVATE

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

   SUBROUTINE construct_lumped_mass_vector(mesh, LA, mass, lumped_mass, opt_per)
      USE st_matrix, ONLY : create_my_ghost
      USE compute_periodic, ONLY : periodic_vector_petsc
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type), INTENT(IN) :: mesh
      type(petsc_csr_LA), INTENT(IN) :: LA
      Mat, INTENT(IN) :: mass
      Vec :: vec_one, lumped_mass
      INTEGER, POINTER, DIMENSION(:) :: ifrom  ! for ghost structure
      INTEGER :: ierr
      LOGICAL, OPTIONAL :: opt_per
      LOGICAL :: per

      CALL pack_opt(per, .TRUE., opt_per)

      !===Create ghost structure
      CALL create_my_ghost(mesh, LA, ifrom)
      CALL VecCreateGhost(PETSC_COMM_WORLD, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, lumped_mass, ierr)
      CALL VecDuplicate(lumped_mass, vec_one, ierr)

      CALL VecSet(vec_one, 1.d0, ierr)
      CALL MatMult(mass, vec_one, lumped_mass, ierr)
      IF (per) THEN
         CALL periodic_vector_petsc(mesh%per%nb_bords, mesh%per%list, mesh%per%perlist, lumped_mass, LA)
      END IF
      CALL VecDestroy(vec_one, ierr)
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
         DO m = 1, mesh%me
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


   SUBROUTINE construct_elasticity_M (mesh, LA, elasticity_M, lambda_viscosity, mu_viscosity)
      !> subroutine to build elasticity matrix of size k_dim
      !! (IN): lambda_viscosity, mu_viscosity
      !! (IN): mesh, LA
      !! (OUT): elasticity_M
      !! WARNING: possible optimizations:
      !!      - fill with petsc_csr_LA_enhanced
      !!      - hardcode type of finite element for loops
      !=================================================
      USE space_dim
      USE def_type_mesh,       ONLY: mesh_type
      USE petsc_csr_la_module, ONLY: petsc_csr_LA
      IMPLICIT NONE
      TYPE(mesh_type),    INTENT(IN) :: mesh
      type(petsc_csr_LA), INTENT(IN) :: LA
      REAL(KIND=8),       INTENT(IN) :: lambda_viscosity, mu_viscosity
      INTEGER, DIMENSION(mesh%gauss%n_w) :: jj_loc
      INTEGER, DIMENSION(k_dim*mesh%gauss%n_w) :: idxm, idxn
      INTEGER :: m, mi, i, ki, iglob, ix, nj, j, kj, jx, l, n_w, ni, k1, jglob
      REAL(KIND=8) :: x, y, lambda_elast
      REAL(KIND = 8), DIMENSION(k_dim*mesh%gauss%n_w, k_dim*mesh%gauss%n_w) :: mat_loc
      TYPE(tMat) :: elasticity_M
      PetscErrorCode                     :: ierr


      lambda_elast = lambda_viscosity-2.d0/3.d0*mu_viscosity

      CALL MatZeroEntries (elasticity_M, ierr)
      CALL MatSetOption (elasticity_M, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
      CALL MatSetOption (elasticity_M, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)
      n_w = mesh%gauss%n_w
      DO m=1, mesh%me
         jj_loc = mesh%jj(:, m)
         mat_loc = 0.d0
         DO ni = 1, n_w
            i = jj_loc(ni)
            DO ki = 1, k_dim
               iglob = LA%loc_to_glob(ki, i)
               ix = (ki - 1) * n_w + ni
               idxm(ix) = iglob - 1
               DO nj = 1, n_w
                  j = jj_loc(nj)
                  DO kj = 1, k_dim
                     jglob = LA%loc_to_glob(kj, j)
                     jx = (kj - 1) * n_w + nj
                     idxn(jx) = jglob - 1
                     x = 0
                     DO l = 1, mesh%gauss%l_G
                        y =  mu_viscosity*mesh%gauss%dw(kj,ni,l,m)*mesh%gauss%dw(ki,nj,l,m) &
                       + lambda_elast*mesh%gauss%dw(ki,ni,l,m)*mesh%gauss%dw(kj,nj,l,m)
                        IF (kj.EQ.ki) THEN
                           DO k1 = 1, k_dim
                              y = y + mu_viscosity*mesh%gauss%dw(k1,ni,l,m)*mesh%gauss%dw(k1,nj,l,m)
                           END DO
                        END IF
                        x = x + y * mesh%gauss%rj(l,m)
                     END DO
                     mat_loc(ix,jx) = x
                  END DO
               END DO
            END DO
         END DO
         CALL MatSetValues(elasticity_M, k_dim * n_w, idxm, k_dim * n_w, idxn, mat_loc, ADD_VALUES, ierr)
      ENDDO

      CALL MatAssemblyBegin(elasticity_M, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd(elasticity_M, MAT_FINAL_ASSEMBLY, ierr)

   END SUBROUTINE construct_elasticity_M

END MODULE fem_petsc_matrix_factory_module
