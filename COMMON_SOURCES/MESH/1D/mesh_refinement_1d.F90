MODULE mesh_refinement_1d
#include "petsc/finclude/petsc.h"
   USE petsc
   USE mesh_tools
   USE mesh_distribution_1d
   PUBLIC :: refinement_P1_mesh_1D, build_jce_1D
   PRIVATE
   LOGICAL, PRIVATE :: per_bool ! <== FIXME
CONTAINS

   SUBROUTINE refinement_P1_mesh_1D(communicator, mesh_P1, mesh_refined, refinement_order)
      USE def_type_mesh
      USE my_util, ONLY: error_petsc, to_str
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_P1, mesh_Pk, mesh_refined
      REAL(KIND=8) :: dx, x0
      INTEGER :: n, m, i, n_shift, i_shift, n_loc, other_proc, other_m_loc, p, k, jj_shift, dom_np_diff, np_diff
      INTEGER, INTENT(IN) :: refinement_order
      INTEGER, DIMENSION(2) :: num_jj
      integer :: ierr
      MPI_Comm :: communicator

      CALL create_Pk_mesh_1D(communicator, mesh_P1, mesh_Pk, refinement_order)
      CALL copy_mesh(mesh_Pk, mesh_refined)

      DEALLOCATE(mesh_refined%jj, mesh_refined%jj_extra)
      mesh_refined%me =  refinement_order * mesh_P1%me
      mesh_refined%medge = mesh_P1%medge * refinement_order
      mesh_refined%dom_me =  refinement_order * mesh_P1%dom_me
      mesh_refined%mextra =  refinement_order * mesh_P1%mextra
      ALLOCATE(mesh_refined%jj(2, mesh_refined%me))
      ALLOCATE(mesh_refined%jj_extra(2, mesh_refined%me))

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
      mesh_refined%neighs(2) = mesh_refined%me


      ! === Extra layer
      DO m=1, mesh_Pk%mextra
         num_jj(1) = 1
         IF (refinement_order > 1) THEN
            num_jj(2) = 3     
         ELSE
            num_jj(2) = 2     
         END IF
         mesh_refined%jj_extra(:, m) = mesh_Pk%jj_extra(num_jj, m)
         DO k=1, refinement_order - 2
            num_jj(1) = 2 + k
            num_jj(2) = num_jj(1) + 1
            mesh_refined%jj_extra(:, m+(k*mesh_Pk%mextra)) = mesh_Pk%jj_extra(num_jj, m)
         END DO
         IF (refinement_order > 1) THEN
            num_jj(1) = refinement_order + 1
            num_jj(2) = 2     
            mesh_refined%jj_extra(:, m+(refinement_order - 1)*mesh_Pk%mextra) = mesh_Pk%jj_extra(num_jj, m)
         END IF
      END DO

!=== DEBUGGING
      ! write(*,*) "loc_to_glob refined", size(mesh_refined%loc_to_glob), MAXVAL(mesh_refined%jj)
      ! write(*,*) "loc_to_glob Pk", size(mesh_pk%loc_to_glob), MAXVAL(mesh_pk%jj)
      ! DO m=1, mesh_Pk%me
         
      !    write(*,*) "Pk element ", m, ' on proc ', mesh_Pk%proc, mesh_Pk%rr(1, mesh_Pk%jj(:, m)), mesh_Pk%jj(:, m)
      !    write(*,*) mesh_Pk%loc_to_glob(mesh_Pk%jj(:, m))
      !    DO k=0, refinement_order-1
      !        n = m + k*mesh_Pk%me
      !        write(*,*) "refined element ", n, ' on proc ', mesh_refined%proc, mesh_refined%rr(1, mesh_refined%jj(:, n)),&
      !        mesh_refined%jj(:, n)
      !        write(*,*) mesh_refined%loc_to_glob(mesh_refined%jj(:, n))
      !    END DO

      ! END DO
      ! write(*,*) "Pk coords dom_np on proc ", mesh_Pk%proc, mesh_Pk%rr(1,:mesh_Pk%dom_np)
      ! write(*,*) "Pk coords np on proc ", mesh_Pk%proc, mesh_Pk%rr(1,mesh_Pk%dom_np+1:)
      ! write(*,*) "refined coords dom_np on proc ", mesh_refined%proc, mesh_refined%rr(1,:mesh_refined%dom_np)
      ! write(*,*) "refined coords np on proc ", mesh_refined%proc, mesh_refined%rr(1,mesh_refined%dom_np+1:)
      ! write(*,*) mesh_refined%loc_to_glob
!=== DEBUGGING
      CALL free_mesh(mesh_Pk)

   END SUBROUTINE refinement_P1_mesh_1D

   SUBROUTINE build_jce_1D(mesh)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      INTEGER :: m, idx
      INTEGER, DIMENSION(:), ALLOCATABLE :: i

      IF (ASSOCIATED(mesh%jce)) NULLIFY(mesh%jce)
      ALLOCATE(mesh%jce(1, mesh%medge))
      
      IF (mesh%me /= mesh%medge) CALL error_petsc('BUG in 1D: me(='//to_str(mesh%me)//') /= &
      medge(='//to_str(mesh%medge)//')??')
      
      mesh%jce(1, :) = [(m, m=1, mesh%me)] + mesh%disedge(mesh%proc) - 1

   END SUBROUTINE build_jce_1D

END MODULE mesh_refinement_1d
