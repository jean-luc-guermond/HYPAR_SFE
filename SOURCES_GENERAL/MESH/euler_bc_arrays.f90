MODULE euler_bc_arrays
  USE mesh_handling
  USE dir_nodes
  USE space_dim
  USE input_data
  INTEGER, POINTER, DIMENSION(:), PUBLIC :: rho_js_D, ux_js_D, uy_js_D, DIR_js_D, whole_bdy_js_D
  INTEGER, POINTER, DIMENSION(:), PUBLIC :: udotn_js_D, surf_udotn_js_D
  REAL(KIND=8), POINTER, DIMENSION(:,:), PUBLIC:: surf_normal_vtx
  REAL(KIND=8), POINTER, DIMENSION(:,:), PUBLIC:: DIR_normal_vtx
  PUBLIC :: construct_euler_bc
  PRIVATE
CONTAINS 
  SUBROUTINE construct_euler_bc
    USE sub_plot
    IMPLICIT NONE
    LOGICAL, POINTER, DIMENSION(:) :: dir
    LOGICAL, DIMENSION(mesh%nps) :: virgin
    REAL(KIND=8), DIMENSION(k_dim,mesh%nps) :: normal_vtx
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: stuff
    INTEGER :: ms, ns, js, n
    ALLOCATE (dir(MAXVAL(mesh%sides)))
    dir = .TRUE.
    CALL dirichlet_nodes(mesh%jjs, mesh%sides, Dir, whole_bdy_js_D)
    dir = .FALSE.
    dir(inputs%rho_Dir_list) = .TRUE.
    CALL dirichlet_nodes(mesh%jjs, mesh%sides, Dir, rho_js_D)
    dir = .FALSE.
    dir(inputs%ux_Dir_list) = .TRUE.
    CALL dirichlet_nodes(mesh%jjs, mesh%sides, Dir, ux_js_D)
    dir = .FALSE.
    dir(inputs%uy_Dir_list) = .TRUE.
    CALL dirichlet_nodes(mesh%jjs, mesh%sides, Dir, uy_js_D)
    dir = .FALSE.
    dir(inputs%Dir_list) = .TRUE.
    CALL dirichlet_nodes(mesh%iis, mesh%sides, Dir, DIR_js_D)


    !===Normal at vertices
    SELECT CASE(k_dim)
    CASE(2)
       ALLOCATE(surf_normal_vtx(k_dim,mesh%nps))
       ALLOCATE(DIR_normal_vtx(k_dim,SIZE(DIR_js_D)))
       surf_normal_vtx = 0.d0
       normal_vtx = 0.d0
       DO ms = 1, mesh%mes
          DO ns = 1, mesh%gauss%n_ws
             js = mesh%iis(ns,ms)
             normal_vtx(:,js) = normal_vtx(:,js) +  mesh%gauss%rnorms_v(:,ns,ms)
          END DO
          IF (MINVAL(ABS(mesh%sides(ms) - inputs%udotn_zero_list)).NE.0) CYCLE 
          DO ns = 1, mesh%gauss%n_ws
             js = mesh%iis(ns,ms)
             surf_normal_vtx(:,js) = surf_normal_vtx(:,js) +  mesh%gauss%rnorms_v(:,ns,ms)
          END DO
       END DO
       DO ns = 1, mesh%nps
          surf_normal_vtx(:,ns) = surf_normal_vtx(:,ns)/SQRT(SUM(surf_normal_vtx(:,ns)**2))
          normal_vtx(:,ns)      = normal_vtx(:,ns)/SQRT(SUM(normal_vtx(:,ns)**2))
       END DO

       dir = .FALSE.
       dir(inputs%udotn_zero_list) = .TRUE.
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
       ALLOCATE(surf_normal_vtx(1,2))
       surf_normal_vtx(1,1) =-1.d0
       surf_normal_vtx(1,2) = 1.d0
       normal_vtx(1,1) =-1.d0
       normal_vtx(1,2) = 1.d0
    END SELECT
    DEALLOCATE(dir)

  END SUBROUTINE construct_euler_bc

END MODULE euler_bc_arrays
