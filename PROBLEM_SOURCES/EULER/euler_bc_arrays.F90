MODULE euler_bc_arrays
   USE dirichlet_type_module, ONLY : dirichlet_bc
   USE def_type_mesh,         ONLY : mesh_type, petsc_csr_LA

   TYPE euler_bc_type
      INTEGER :: syst_dim
      TYPE(dirichlet_bc) :: rho_bc, ux_bc, uy_bc, whole_bdy_bc, udotn_bc
      REAL(KIND = 8), POINTER, DIMENSION(:, :) :: udotn_normal_vtx
   CONTAINS
      PROCEDURE :: construct_euler_bc
   END TYPE euler_bc_type
CONTAINS
   SUBROUTINE construct_euler_bc(this, mesh, LA)
      USE petsc
#include "petsc/finclude/petsc.h"
      !  USE sub_plot
      USE space_dim
      USE euler_matrices_module, ONLY : x1vec, x2vec, x2_ghost
      USE st_matrix,             ONLY : extract_through_ghost
      IMPLICIT NONE
      CLASS(euler_bc_type)                       :: this
      TYPE(mesh_type)                            :: mesh
      TYPE(petsc_csr_LA)                         :: LA
      LOGICAL,        DIMENSION(mesh%nps)        :: virgin
      REAL(KIND = 8), DIMENSION(mesh%nps, k_dim) :: normal_vtx
      REAL(KIND = 8), DIMENSION(mesh%np)         :: dummy_normal_vtx
      INTEGER,        DIMENSION(SIZE(mesh%jjs,1))        :: idxms
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: stuff
      REAL(KIND = 8) :: norm
      INTEGER :: ms, ns, js, n, ierr

      CALL this%rho_bc%set(mesh, "density", "DIRICHLET BC PARAMETERS")

      CALL this%ux_bc%set(mesh, "ux")
      
      IF (k_dim>1) THEN
         CALL this%uy_bc%set(mesh, "uy")
         CALL this%whole_bdy_bc%set(mesh, "whole boundary")
         CALL this%udotn_bc%set(mesh, "u.n=0")
      END IF

      !===Normal at vertices
      SELECT CASE(k_dim)
      CASE(2)

         CALL VecZeroEntries(x1vec, ierr)
         CALL VecZeroEntries(x2vec, ierr)

         virgin = .TRUE.
         n = 0
         DO ms = 1, mesh%mes

            IF (MINVAL(ABS(mesh%sides(ms) - this%udotn_bc%list_sides)).NE.0) CYCLE
            DO ns = 1, mesh%gauss%n_ws
               js = mesh%iis(ns,ms)
               IF (virgin(js)) THEN
                  virgin(js) =.FALSE.
                  n = n+1
                  this%udotn_bc%jsd(n) = mesh%jjs(ns,ms)
               END IF
               normal_vtx(js,:) = mesh%gauss%rnorms_v(:,ns,ms)
            END DO
            idxms = LA%loc_to_glob(1, mesh%jjs(:,ms)) -1

            CALL VecSetValues(x1vec, mesh%gauss%n_ws, idxms, normal_vtx(mesh%iis(:,ms),1), ADD_VALUES, ierr)
            CALL VecSetValues(x2vec, mesh%gauss%n_ws, idxms, normal_vtx(mesh%iis(:,ms),2), ADD_VALUES, ierr)
         END DO

         CALL extract_through_ghost(x1vec, x2_ghost, 1, 1, LA, dummy_normal_vtx(:), 'insert', opt_assemble=.TRUE.)
         DO ms = 1, mesh%mes
            IF (.NOT. ANY(mesh%sides(ms) == this%udotn_bc%list_sides)) CYCLE
            normal_vtx(mesh%iis(:,ms), 1) = dummy_normal_vtx(mesh%jjs(:,ms))
         END DO

         CALL extract_through_ghost(x2vec, x2_ghost, 1, 1, LA, dummy_normal_vtx(:), 'insert', opt_assemble=.TRUE.)
         DO ms = 1, mesh%mes
            IF (.NOT. ANY(mesh%sides(ms) == this%udotn_bc%list_sides)) CYCLE
            normal_vtx(mesh%iis(:,ms), 2) = dummy_normal_vtx(mesh%jjs(:,ms))
         END DO

         ALLOCATE(this%udotn_normal_vtx(SIZE(this%udotn_bc%jsd),k_dim))
         n = 0
         virgin = .TRUE.
         DO ms = 1, mesh%mes
            IF (MINVAL(ABS(mesh%sides(ms) - this%udotn_bc%list_sides)).NE.0) CYCLE
            DO ns = 1, mesh%gauss%n_ws
               js = mesh%iis(ns,ms)
               IF (virgin(js)) THEN
                  virgin(js) =.FALSE.
                  n = n+1
                  this%udotn_normal_vtx(n,:) = normal_vtx(js,:)/SQRT(SUM(normal_vtx(js,:)**2))
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
         ALLOCATE(this%udotn_bc%jsd(0))
         ALLOCATE(this%udotn_normal_vtx(0,1))
      END SELECT

   END SUBROUTINE construct_euler_bc

END MODULE euler_bc_arrays
