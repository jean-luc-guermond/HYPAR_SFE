MODULE euler_bc_arrays
   USE dirichlet_type_module

   TYPE euler_bc_type
      INTEGER :: syst_dim
      TYPE(dirichlet_bc) :: rho_bc, ux_bc, uy_bc, whole_bdy_bc, udotn_bc
      REAL(KIND = 8), POINTER, DIMENSION(:, :) :: udotn_normal_vtx
   CONTAINS
      PROCEDURE, PUBLIC :: construct_euler_bc
   END TYPE euler_bc_type
CONTAINS
  SUBROUTINE construct_euler_bc(this, mesh)
    USE sub_plot
    USE dir_nodes_petsc
    USE def_type_mesh
    USE mesh_parameters
   !  USE space_dim
    IMPLICIT NONE
    CLASS(euler_bc_type) :: this
    TYPE(mesh_type)      :: mesh
    LOGICAL,        DIMENSION(mesh%nps)        :: virgin
    REAL(KIND = 8), DIMENSION(mesh%nps, mesh_data_info%k_dim) :: normal_vtx
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: stuff
    INTEGER :: ms, ns, js, n
!!$    CHARACTER(LEN=5) :: char

    CALL this%rho_bc%set(mesh, "density", "DIRICHLET BC PARAMETERS")

    CALL this%ux_bc%set(mesh, "ux")
 
    IF (mesh_data_info%k_dim>1) THEN
       CALL this%uy_bc%set(mesh, "uy")
       CALL this%whole_bdy_bc%set(mesh, "whole boundary")
       CALL this%udotn_bc%set(mesh, "u.n=0")
    END IF

    !===Normal at vertices
    SELECT CASE(mesh_data_info%k_dim)
    CASE(2)
       normal_vtx = 0.d0
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
             normal_vtx(js,:) = normal_vtx(js,:) + mesh%gauss%rnorms_v(:,ns,ms)
          END DO
       END DO
       ALLOCATE(this%udotn_normal_vtx(SIZE(this%udotn_bc%jsd),mesh_data_info%k_dim))
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
