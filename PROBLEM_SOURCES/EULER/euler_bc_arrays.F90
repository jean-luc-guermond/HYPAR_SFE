MODULE euler_bc_arrays
  TYPE euler_bc_type
  INTEGER, POINTER, DIMENSION(:) :: rho_js_D, ux_js_D, uy_js_D, DIR_js_D, whole_bdy_js_D
  INTEGER, POINTER, DIMENSION(:) :: udotn_js_D, surf_udotn_js_D
  REAL(KIND=8), POINTER, DIMENSION(:,:):: surf_normal_vtx
  REAL(KIND=8), POINTER, DIMENSION(:,:):: DIR_normal_vtx
  CONTAINS
    PROCEDURE, PUBLIC :: construct_euler_bc
 END TYPE euler_bc_type
  PRIVATE
CONTAINS 
  SUBROUTINE construct_euler_bc(this)
    USE sub_plot
    USE dir_nodes_petsc
    IMPLICIT NONE
    TYPE(euler_type) :: this
    LOGICAL, POINTER, DIMENSION(:) :: dir
    LOGICAL, DIMENSION(this%mesh%nps) :: virgin
    REAL(KIND=8), DIMENSION(k_dim,this%mesh%nps) :: normal_vtx
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: stuff
    INTEGER :: ms, ns, js, n
    ALLOCATE (dir(MAXVAL(this%mesh%sides)))
    dir = .TRUE.
    CALL dirichlet_nodes(this%mesh%jjs, this%mesh%sides, Dir, this%euler_bc%whole_bdy_js_D)
    dir = .FALSE.
    dir(inputs%rho_Dir_list) = .TRUE.
    CALL dirichlet_nodes(this%mesh%jjs, this%mesh%sides, Dir,  this%euler_bc%rho_js_D)
    dir = .FALSE.
    dir(inputs%ux_Dir_list) = .TRUE.
    CALL dirichlet_nodes(this%mesh%jjs, this%mesh%sides, Dir,  this%euler_bc%ux_js_D)
    dir = .FALSE.
    dir(inputs%uy_Dir_list) = .TRUE.
    CALL dirichlet_nodes(this%mesh%jjs, this%mesh%sides, Dir,  this%euler_bc%uy_js_D)
    dir = .FALSE.
    dir(inputs%Dir_list) = .TRUE.
    CALL dirichlet_nodes(this%mesh%iis, this%mesh%sides, Dir,  this%euler_bc%DIR_js_D)


!!$    !===Normal at vertices
!!$    SELECT CASE(k_dim)
!!$    CASE(2)
!!$       ALLOCATE(surf_normal_vtx(k_dim,mesh%nps))
!!$       ALLOCATE(DIR_normal_vtx(k_dim,SIZE(DIR_js_D)))
!!$       surf_normal_vtx = 0.d0
!!$       normal_vtx = 0.d0
!!$       DO ms = 1, mesh%mes
!!$          DO ns = 1, mesh%gauss%n_ws
!!$             js = mesh%iis(ns,ms)
!!$             normal_vtx(:,js) = normal_vtx(:,js) +  mesh%gauss%rnorms_v(:,ns,ms)
!!$          END DO
!!$          IF (MINVAL(ABS(mesh%sides(ms) - inputs%udotn_zero_list)).NE.0) CYCLE 
!!$          DO ns = 1, mesh%gauss%n_ws
!!$             js = mesh%iis(ns,ms)
!!$             surf_normal_vtx(:,js) = surf_normal_vtx(:,js) +  mesh%gauss%rnorms_v(:,ns,ms)
!!$          END DO
!!$       END DO
!!$       DO ns = 1, mesh%nps
!!$          surf_normal_vtx(:,ns) = surf_normal_vtx(:,ns)/SQRT(SUM(surf_normal_vtx(:,ns)**2))
!!$          normal_vtx(:,ns)      = normal_vtx(:,ns)/SQRT(SUM(normal_vtx(:,ns)**2))
!!$       END DO
!!$
!!$       dir = .FALSE.
!!$       dir(inputs%udotn_zero_list) = .TRUE.
!!$       CALL dirichlet_nodes(mesh%iis, mesh%sides, Dir, surf_udotn_js_D)
!!$       ALLOCATE(udotn_js_D(SIZE(surf_udotn_js_D)))
!!$       udotn_js_D = mesh%j_s(surf_udotn_js_D)
!!$
!!$       virgin = .TRUE. 
!!$       virgin(DIR_js_D) = .FALSE.
!!$       n = 0
!!$       DO ns = 1, mesh%nps
!!$          IF (virgin(ns)) CYCLE
!!$          n = n+1
!!$          DIR_js_D(n) = mesh%j_s(ns)
!!$          DIR_normal_vtx(:,n) = normal_vtx(:,ns)
!!$       END DO
!!$       !===Check normal vector
!!$       ALLOCATE(stuff(k_dim,mesh%np))
!!$       stuff = 0.d0
!!$       stuff(1,DIR_js_D)  = DIR_normal_vtx(1,:)
!!$       stuff(2,DIR_js_D)  = DIR_normal_vtx(2,:)
!!$       CALL plot_arrow_label(mesh%jj, mesh%rr, stuff, 'normal.plt')
!!$       DEALLOCATE(stuff)
!!$    CASE(1)
!!$       ALLOCATE(udotn_js_D(0))
!!$       ALLOCATE(surf_normal_vtx(1,2))
!!$       surf_normal_vtx(1,1) =-1.d0
!!$       surf_normal_vtx(1,2) = 1.d0
!!$       normal_vtx(1,1) =-1.d0
!!$       normal_vtx(1,2) = 1.d0
!!$    END SELECT
    DEALLOCATE(dir)

  END SUBROUTINE construct_euler_bc

END MODULE euler_bc_arrays
