MODULE hyperbolic_bc_tools

CONTAINS

   SUBROUTINE construct_udotn(mesh, LA, udotn_bc, udotn_normal_vtx)
      USE petsc
#include "petsc/finclude/petsc.h"
      !  USE sub_plot
      USE dirichlet_type_module, ONLY : dirichlet_bc
      USE space_dim
      USE def_type_mesh, ONLY: mesh_type, petsc_csr_LA
      USE hyperbolic_matrices_module, ONLY : x1vec, x2vec, x2_ghost
      USE st_matrix,                  ONLY : extract_through_ghost
      IMPLICIT NONE
      TYPE(mesh_type),                 INTENT(IN) :: mesh
      TYPE(petsc_csr_LA),              INTENT(IN) :: LA
      TYPE(dirichlet_bc),                           INTENT(INOUT):: udotn_bc
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :), INTENT(OUT):: udotn_normal_vtx
      LOGICAL,        DIMENSION(mesh%nps)        :: virgin
      REAL(KIND = 8), DIMENSION(mesh%nps, k_dim) :: normal_vtx
      REAL(KIND = 8), DIMENSION(mesh%np)         :: dummy_normal_vtx
      INTEGER,        DIMENSION(SIZE(mesh%jjs,1)):: idxms
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: stuff
      REAL(KIND = 8) :: norm
      INTEGER :: ms, ns, js, n, ierr

      !===Normal at vertices
      SELECT CASE(k_dim)
      CASE(2)

         CALL VecZeroEntries(x1vec, ierr)
         CALL VecZeroEntries(x2vec, ierr)

         virgin = .TRUE.
         n = 0
         DO ms = 1, mesh%mes

            IF (.NOT. ANY(mesh%sides(ms)==udotn_bc%list_sides)) CYCLE
            DO ns = 1, mesh%gauss%n_ws
               js = mesh%iis(ns,ms)
               IF (virgin(js)) THEN
                  virgin(js) =.FALSE.
                  n = n+1
                  udotn_bc%jsd(n) = mesh%jjs(ns,ms)
               END IF
               normal_vtx(js,:) = mesh%gauss%rnorms_v(:,ns,ms)
            END DO
            idxms = LA%loc_to_glob(1, mesh%jjs(:,ms)) -1

            CALL VecSetValues(x1vec, mesh%gauss%n_ws, idxms, normal_vtx(mesh%iis(:,ms),1), ADD_VALUES, ierr)
            CALL VecSetValues(x2vec, mesh%gauss%n_ws, idxms, normal_vtx(mesh%iis(:,ms),2), ADD_VALUES, ierr)
         END DO

         CALL extract_through_ghost(x1vec, x2_ghost, 1, 1, LA, dummy_normal_vtx(:), 'insert', opt_assemble=.TRUE.)
         DO ms = 1, mesh%mes
            IF (.NOT. ANY(mesh%sides(ms) == udotn_bc%list_sides)) CYCLE
            normal_vtx(mesh%iis(:,ms), 1) = dummy_normal_vtx(mesh%jjs(:,ms))
         END DO

         CALL extract_through_ghost(x2vec, x2_ghost, 1, 1, LA, dummy_normal_vtx(:), 'insert', opt_assemble=.TRUE.)
         DO ms = 1, mesh%mes
            IF (.NOT. ANY(mesh%sides(ms) == udotn_bc%list_sides)) CYCLE
            normal_vtx(mesh%iis(:,ms), 2) = dummy_normal_vtx(mesh%jjs(:,ms))
         END DO

         ALLOCATE(udotn_normal_vtx(SIZE(udotn_bc%jsd),k_dim))
         n = 0
         virgin = .TRUE.
         DO ms = 1, mesh%mes
            IF (MINVAL(ABS(mesh%sides(ms) - udotn_bc%list_sides)).NE.0) CYCLE
            DO ns = 1, mesh%gauss%n_ws
               js = mesh%iis(ns,ms)
               IF (virgin(js)) THEN
                  virgin(js) =.FALSE.
                  n = n+1
                  udotn_normal_vtx(n,:) = normal_vtx(js,:)/SQRT(SUM(normal_vtx(js,:)**2))
               END IF
            END DO
         END DO

!!$       !===Check normal vector
!!$       ALLOCATE(stuff(mesh_data_info%k_dim,mesh%np))
!!$       stuff = 0.d0
!!$       stuff(1,this%udotn_bc%jsd) = this%udotn_normal_vtx(:,1)
!!$       stuff(2,this%udotn_bc%jsd) = this%udotn_normal_vtx(:,2)
!!$       WRITE(char, '(I5)') mesh%rank    
!!$       CALL plot_arrow_label(mesh%jj, mesh%rr, stuff, 'normal'//trim(adjustl(char))//'.plt')
!!$       DEALLOCATE(stuff)
      CASE(1)
         ALLOCATE(udotn_bc%jsd(0))
         ALLOCATE(udotn_normal_vtx(0,1))
      END SELECT

   END SUBROUTINE construct_udotn
   
END MODULE hyperbolic_bc_tools