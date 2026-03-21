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
    USE space_dim
    IMPLICIT NONE
    CLASS(euler_bc_type) :: this
    TYPE(mesh_type)      :: mesh
    LOGICAL,        DIMENSION(mesh%nps)        :: virgin
    REAL(KIND = 8), DIMENSION(k_dim, mesh%nps) :: normal_vtx
    REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: stuff
    INTEGER :: ms, ns, js, n

    !CALL this%euler_bc%whole_bdy_js_D%set(mesh, "")

    CALL this%rho_bc%set(mesh, "density")

    CALL this%ux_bc%set(mesh, "ux")
  write(*,*) ' OOO'
    IF (k_dim>1) THEN
       CALL this%uy_bc%set(mesh, "uy")
       CALL this%whole_bdy_bc%set(mesh, "whole boundary")
       CALL this%udotn_bc%set(mesh, "u.n=0")
    END IF
 
    !CALL this%DIR_js_D%set(mesh, "")

    write(*,*) ' OOO'
    STOP
    !===Normal at vertices
    SELECT CASE(k_dim)
    CASE(2)
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
             normal_vtx(:,js) = normal_vtx(:,js) + mesh%gauss%rnorms_v(:,ns,ms)
          END DO
       END DO
       DO ns = 1, mesh%nps
          normal_vtx(:,ns) = &
               normal_vtx(:,ns)/SQRT(SUM(normal_vtx(:,ns)**2))
       END DO
       ALLOCATE(this%udotn_normal_vtx(k_dim,SIZE(this%udotn_bc%jsd)))
       DO ms = 1, mesh%mes
          IF (MINVAL(ABS(mesh%sides(ms) - this%udotn_bc%list_sides)).NE.0) CYCLE
          DO ns = 1, mesh%gauss%n_ws
             js = mesh%iis(ns,ms)
             IF (virgin(js)) THEN
                virgin(js) =.FALSE.
                n = n+1
                this%udotn_normal_vtx(:,n) = normal_vtx(:,js)
             END IF
          END DO
       END DO
    write(*,*) ' OOOJJJJJJ'   
       !===Check normal vector
       ALLOCATE(stuff(k_dim,mesh%np))
       stuff = 0.d0
       stuff(1,this%udotn_bc%jsd) = this%udotn_normal_vtx(1,:)
       stuff(2,this%udotn_bc%jsd) = this%udotn_normal_vtx(2,:)
       CALL plot_arrow_label(mesh%jj, mesh%rr, stuff, 'normal.plt')
       DEALLOCATE(stuff)
    CASE(1)
       ALLOCATE(this%udotn_bc%jsd(0))
       ALLOCATE(this%udotn_normal_vtx(1,0))
    END SELECT

  END SUBROUTINE construct_euler_bc

END MODULE euler_bc_arrays
