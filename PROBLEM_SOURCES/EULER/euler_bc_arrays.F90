MODULE euler_bc_arrays
   USE dirichlet_type_module

   TYPE euler_bc_type
      INTEGER :: syst_dim
      TYPE(dirichlet_bc) :: rho_bc, ux_bc, uy_bc, DIR_bc, whole_bdy_bc
      TYPE(dirichlet_bc) :: udotn_bc, surf_udotn_bc
      REAL(KIND = 8), POINTER, DIMENSION(:, :) :: surf_normal_vtx
      REAL(KIND = 8), POINTER, DIMENSION(:, :) :: DIR_normal_vtx
   CONTAINS
      PROCEDURE, PUBLIC :: construct_euler_bc
   END TYPE euler_bc_type
CONTAINS
  SUBROUTINE construct_euler_bc(this, mesh)
    USE sub_plot
    USE dir_nodes_petsc
    USE def_type_mesh
    USE space_dim
    IMPLICIT NONE
    CLASS(euler_bc_type) :: this
    TYPE(mesh_type) :: mesh
    LOGICAL, DIMENSION(mesh%nps) :: virgin
    REAL(KIND = 8), DIMENSION(k_dim, mesh%nps) :: normal_vtx
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: stuff
    INTEGER :: ms, ns, js, n

    !CALL this%euler_bc%whole_bdy_js_D%set(mesh, "")

    CALL this%rho_bc%set(mesh, "density")

    CALL this%ux_bc%set(mesh, "ux")

    IF (k_dim>1) THEN
       CALL this%uy_bc%set(mesh, "uy")
       CALL this%whole_bdy_bc%set(mesh, "whole boundary")
       CALL this%udotn_bc%set(mesh, "u.n=0")
    END IF
    return
    !CALL this%DIR_js_D%set(mesh, "")


    !===Normal at vertices
    SELECT CASE(k_dim)
    CASE(2)
       ALLOCATE(this%surf_normal_vtx(k_dim,mesh%nps))
       ALLOCATE(this%DIR_normal_vtx(k_dim,SIZE(DIR_js_D)))
       surf_normal_vtx = 0.d0
       normal_vtx = 0.d0
       DO ms = 1, mesh%mes
          DO ns = 1, mesh%gauss%n_ws
             js = mesh%iis(ns,ms)
             normal_vtx(:,js) = normal_vtx(:,js) +  mesh%gauss%rnorms_v(:,ns,ms)
          END DO
          !IF (MINVAL(ABS(mesh%sides(ms) - inputs%udotn_zero_list)).NE.0) CYCLE
          IF (MINVAL(ABS(mesh%sides(ms) - this%udotn_bc%list_sides)).NE.0) CYCLE
          DO ns = 1, mesh%gauss%n_ws
             js = mesh%iis(ns,ms)
             this%surf_normal_vtx(:,js) = this%surf_normal_vtx(:,js) + mesh%gauss%rnorms_v(:,ns,ms)
          END DO
       END DO
       DO ns = 1, mesh%nps
          this%surf_normal_vtx(:,ns) = this%surf_normal_vtx(:,ns)/SQRT(SUM(this%surf_normal_vtx(:,ns)**2))
          normal_vtx(:,ns)      = normal_vtx(:,ns)/SQRT(SUM(normal_vtx(:,ns)**2))
       END DO

       dir = .FALSE.
       dir(this%udotn_bc%list_sides) = .TRUE.
       CALL dirichlet_nodes(mesh%iis, mesh%sides, Dir, surf_udotn_js_D)
       ALLOCATE(udotn_js_D(SIZE(surf_udotn_js_D)))
       udotn_js_D = mesh%j_s(surf_udotn_js_D)

       virgin = .TRUE.
       virgin(DIR_js_D) = .FALSE.
       n = 0
       DO ns = 1, mesh%nps
          IF (virgin(ns)) CYCLE
          n = n+1
          DIR_js_D(n) = mesh%j_s(ns)
          DIR_normal_vtx(:,n) = normal_vtx(:,ns)
       END DO
       !===Check normal vector
       ALLOCATE(stuff(k_dim,mesh%np))
       stuff = 0.d0
       stuff(1,DIR_js_D)  = DIR_normal_vtx(1,:)
       stuff(2,DIR_js_D)  = DIR_normal_vtx(2,:)
       CALL plot_arrow_label(mesh%jj, mesh%rr, stuff, 'normal.plt')
       DEALLOCATE(stuff)
    CASE(1)
       ALLOCATE(udotn_js_D(0))
       ALLOCATE(surf_normal_vtx(1,2))
       surf_normal_vtx(1,1) =-1.d0
       surf_normal_vtx(1,2) = 1.d0
       normal_vtx(1,1) =-1.d0
       normal_vtx(1,2) = 1.d0
    END SELECT

  END SUBROUTINE construct_euler_bc

END MODULE euler_bc_arrays
