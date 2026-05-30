MODULE mesh_refinement_1d
#include "petsc/finclude/petsc.h"
   USE petsc
   USE mesh_tools
   USE mesh_distribution_1d
   PUBLIC :: refinement_P1_mesh_1D, build_jce_1D
   PRIVATE
CONTAINS

   SUBROUTINE refinement_P1_mesh_1D_arbitrary_order(mesh_Pk, mesh_refined)
      USE def_type_mesh
      USE my_util, ONLY: to_str, local_error_petsc
      IMPLICIT NONE
      TYPE(mesh_type)       :: mesh_Pk, mesh_refined
      INTEGER :: n, m, i, p, k, refinement_order, old_me, old_cell_g, old_cell_l, mextra, new_cell_l
      INTEGER, DIMENSION(2) :: num_jj
      LOGICAL               :: test

      refinement_order = SIZE(mesh_Pk%jj, 1) - 1
      CALL copy_mesh(mesh_Pk, mesh_refined)

      DEALLOCATE(mesh_refined%jj, mesh_refined%jj_extra)

      mesh_refined%me     = refinement_order * mesh_Pk%me
      mesh_refined%medge  = refinement_order * mesh_Pk%medge
      mesh_refined%mextra = refinement_order * mesh_Pk%mextra

      ALLOCATE(mesh_refined%jj(2, mesh_refined%me))
      ALLOCATE(mesh_refined%jj_extra(2, mesh_refined%mextra))

      CALL mesh_refined%gather_dom_np
      CALL mesh_refined%gather_me
      CALL mesh_refined%gather_medge

      !=== Bulk
      DO m=1, mesh_Pk%me
         num_jj(1) = 1
         IF (refinement_order > 1) THEN
            num_jj(2) = 3     
         ELSE
            num_jj(2) = 2     
         END IF
         mesh_refined%jj(:, m) = mesh_Pk%jj(num_jj, m)
         DO k=1, refinement_order - 2
            num_jj(1) = 2 + k
            num_jj(2) = num_jj(1) + 1
            mesh_refined%jj(:, m+(k*mesh_Pk%me)) = mesh_Pk%jj(num_jj, m)
         END DO
         IF (refinement_order > 1) THEN
            num_jj(1) = refinement_order + 1
            num_jj(2) = 2     
            mesh_refined%jj(:, m+(refinement_order - 1)*mesh_Pk%me) = mesh_Pk%jj(num_jj, m)
         END IF
      END DO

      ! === Boundary
      DO m=1, mesh_refined%mes
         old_me = mesh_Pk%neighs(m)
         DO n=1, 2
            IF (mesh_Pk%neigh(n, m) <= 0) THEN
               test = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. test) THEN
            CALL local_error_petsc('BUG in refinement 1D => did not find neigh corresponding to neighs')
         END IF
         mesh_refined%neighs(m) = old_me + (refinement_order - 1)*mesh_Pk%me*(n-1)
      END DO

      ! === Extra layer
      mextra = 0
      DO m=1, mesh_Pk%mextra
         mextra = mextra + 1
         num_jj(1) = 1
         IF (refinement_order > 1) THEN
            num_jj(2) = 3
         ELSE
            num_jj(2) = 2     
         END IF
         mesh_refined%jj_extra(:, mextra) = mesh_Pk%jj_extra(num_jj, m)
         DO k=1, refinement_order - 2
            mextra = mextra + 1
            num_jj(1) = 2 + k
            num_jj(2) = num_jj(1) + 1
            mesh_refined%jj_extra(:, mextra) = mesh_Pk%jj_extra(num_jj, m)
         END DO
         IF (refinement_order > 1) THEN
            mextra = mextra + 1
            num_jj(1) = refinement_order + 1
            num_jj(2) = 2     
            mesh_refined%jj_extra(:, mextra) = mesh_Pk%jj_extra(num_jj, m)
         END IF
      END DO

      DEALLOCATE(mesh_refined%jcc_extra)
      ALLOCATE(mesh_refined%jcc_extra(mesh_refined%mextra))

      mextra = 0
      DO m=1, mesh_Pk%mextra
         old_cell_g = mesh_Pk%jcc_extra(m)
         p = mesh_Pk%get_proc(old_cell_g, 'me')
         old_cell_l = old_cell_g - (mesh_Pk%discell(p) - 1)
         DO i=1, refinement_order
            mextra = mextra + 1
            new_cell_l = old_cell_l + (i-1)*mesh_Pk%domcell(p)
            mesh_refined%jcc_extra(mextra) = new_cell_l + (mesh_refined%discell(p) - 1)
         END DO
      END DO


!=== DEBUGGING
      ! write(*,*) "loc_to_glob refined", size(mesh_refined%loc_to_glob), MAXVAL(mesh_refined%jj)
      ! write(*,*) "loc_to_glob Pk", size(mesh_pk%loc_to_glob), MAXVAL(mesh_pk%jj)
      ! ! DO m=1, mesh_Pk%me
         
      ! !    write(*,*) "Pk element ", m, ' on proc ', mesh_Pk%proc, mesh_Pk%rr(1, mesh_Pk%jj(:, m)), mesh_Pk%jj(:, m)
      ! !    write(*,*) mesh_Pk%loc_to_glob(mesh_Pk%jj(:, m))
      ! !    DO k=0, refinement_order-1
      ! !        n = m + k*mesh_Pk%me
      ! !        write(*,*) "refined element ", n, ' on proc ', mesh_refined%proc, mesh_refined%rr(1, mesh_refined%jj(:, n)),&
      ! !        mesh_refined%jj(:, n)
      ! !        write(*,*) mesh_refined%loc_to_glob(mesh_refined%jj(:, n))
      ! ! !    END DO

      ! ! ! END DO
      ! ! write(*,*) "Pk coords dom_np on proc ", mesh_Pk%proc, mesh_Pk%rr(1,:mesh_Pk%dom_np)
      ! ! write(*,*) "Pk coords np on proc ", mesh_Pk%proc, mesh_Pk%rr(1,mesh_Pk%dom_np+1:)
      ! ! write(*,*) "refined coords dom_np on proc ", mesh_refined%proc, mesh_refined%rr(1,:mesh_refined%dom_np)
      ! ! write(*,*) "refined coords np on proc ", mesh_refined%proc, mesh_refined%rr(1,mesh_refined%dom_np+1:)
      ! ! write(*,*) "loc_to_glob: ", mesh_refined%loc_to_glob
      ! DO m=1, SIZE(mesh_Pk%jj_extra, 2)
      !    write(*,*) m, "Pk extra layer ", mesh_Pk%proc, mesh_Pk%jj_extra(:, m)!, mesh_Pk%rr(1,:mesh_Pk%dom_np)
      ! END DO
      ! DO m=1, SIZE(mesh_Pk%jcc_extra)
      !    write(*,*) "Pk jcc_extra", mesh_Pk%proc, mesh_Pk%jcc_extra(m)!, mesh_Pk%rr(1,:mesh_Pk%dom_np)
      ! END DO
      ! DO m=1, SIZE(mesh_refined%jcc_extra)
      !    write(*,*) "refined jcc_extra", mesh_refined%proc, mesh_refined%jcc_extra(m)!, mesh_Pk%rr(1,:mesh_Pk%dom_np)
      ! END DO
      ! ! write(*,*) "Pk coords np on proc ", mesh_Pk%proc, mesh_Pk%rr(1,mesh_Pk%dom_np+1:)
      ! ! write(*,*) "refined coords dom_np on proc ", mesh_refined%proc, mesh_refined%rr(1,:mesh_refined%dom_np)
      ! ! write(*,*) "refined coords np on proc ", mesh_refined%proc, mesh_refined%rr(1,mesh_refined%dom_np+1:)
      ! ! write(*,*) "loc_to_glob: ", mesh_refined%loc_to_glob
!=== DEBUGGING
      CALL free_mesh(mesh_Pk)
      CALL build_jce_1D(mesh_refined)

   END SUBROUTINE refinement_P1_mesh_1D_arbitrary_order

   SUBROUTINE refinement_P1_mesh_1D(mesh_P1, mesh_refined, refinement_order)
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type)     :: mesh_P1, mesh_Pk, mesh_refined
      INTEGER, INTENT(IN) :: refinement_order

      CALL create_Pk_mesh_1D(mesh_P1, mesh_Pk, refinement_order)
      CALL refinement_P1_mesh_1D_arbitrary_order(mesh_Pk, mesh_refined)
   END SUBROUTINE refinement_P1_mesh_1D

   SUBROUTINE build_jce_1D(mesh)
      USE def_type_mesh
      USE my_util, ONLY: local_error_petsc, to_str
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      INTEGER :: m

      IF (ASSOCIATED(mesh%jce)) NULLIFY(mesh%jce)
      ALLOCATE(mesh%jce(1, mesh%medge))
      
      IF (mesh%me /= mesh%medge) CALL local_error_petsc(&
      &'BUG in 1D: me(='//to_str(mesh%me)//') /= &
      &medge(='//to_str(mesh%medge)//')??')
      
      mesh%jce(1, :) = [(m, m=1, mesh%me)] + mesh%disedge(mesh%proc) - 1

   END SUBROUTINE build_jce_1D

END MODULE mesh_refinement_1d
