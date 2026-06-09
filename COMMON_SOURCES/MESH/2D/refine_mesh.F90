MODULE refine_mesh
   USE mesh_tools
   PUBLIC :: create_iso_grid_distributed, refinement_iso_grid_distributed, general_refinement_iso_grid_distributed
   PRIVATE
CONTAINS

   SUBROUTINE create_iso_grid_distributed(mesh_p1, mesh, type_fe)
      !===jj(:, :)    nodes of the  volume_elements of the input grid
      !===jjs(:, :)    nodes of the surface_elements of the input grid
      !===rr(:, :)    cartesian coordinates of the nodes of the input grid
      !===m_op(:,:)   volume element opposite to each node
      !===neigh_el(:) volume element ajacent to the surface element
      !===jj_f(:, :)  nodes of the  volume_elements of the output p2 grid
      !===jjs_f(:, :)  nodes of the surface_elements of the output p2 grid
      !===rr_f(:, :)  cartesian coordinates of the nodes of the output p2 grid
      USE def_type_mesh
      USE my_util, ONLY : error_petsc, to_str, local_error_petsc
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_p1, mesh
      INTEGER, INTENT(IN) :: type_fe
      LOGICAL, DIMENSION(:), ALLOCATABLE :: virgin
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: j_edge, jjs_edge
      INTEGER :: np, me, mes, mes_int, nw, nws, kd, n, m, k, l, n_dof, dom_np, ne
      INTEGER :: n1, n2, n3, ms, n_start, n_end, n_loc, n_loc_dof
      INTEGER :: n_k1, n_k2, m_op_k, kk, i, mm, p_e, p_c
      INTEGER, DIMENSION(2) :: n_ks
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: r_mid
      INTEGER :: nb_angle, f_dof_edge, f_dof_cell,edge_g, edge_l, n_new_start, proc, nb_proc, edges, p, cell_g, cell_l
      INTEGER :: interface
      LOGICAL :: iso


      !=== dummy to avoid warnings (TODO => add curved interfaces)
      interface = -1
      !===

      CALL mesh%info%copy(mesh_p1%info)

      nw = SIZE(mesh_p1%jj, 1)   !===nodes in each volume element (3 in 2D)
      ne = 3                     !===Number of edges in a cell (3 for a triangle)
      me = SIZE(mesh_p1%jj, 2)   !===number of cells
      kd = SIZE(mesh_p1%rr, 1)   !===space dimensions
      np = mesh_p1%np            !===number of P1 vertices connected to grid
      dom_np = mesh_p1%dom_np    !===number of P1 vertices attributed to proc
      mes = SIZE(mesh_p1%jjs, 2)
      mes_int = SIZE(mesh_p1%jjs_int, 2)
      nws = SIZE(mesh_p1%jjs, 1)
      f_dof_edge = type_fe - 1
      f_dof_cell = (type_fe - 1) * (type_fe - 2) / 2
      nb_proc = SIZE(mesh_p1%domnp)

      ! Fix VB 27/05
      ! ALLOCATE(mesh%jjs_int(SIZE(mesh_p1%jjs_int, 1), SIZE(mesh_p1%jjs_int, 2)))
      ! mesh%jjs_int = mesh_p1%jjs_int
      ALLOCATE(mesh%sides_int(SIZE(mesh_p1%sides_int)))
      mesh%sides_int = mesh_p1%sides_int
      ALLOCATE(mesh%neighs_int(2, SIZE(mesh_p1%neighs_int, 2)))
      mesh%neighs_int = mesh_p1%neighs_int
      ! Fix VB 27/05

!=== Unchanged number of edges and elements
      mesh%me = mesh_p1%me
      mesh%mes = mesh_p1%mes
      mesh%medge = mesh_p1%medge
      mesh%medges = mesh_p1%medges
      mesh%mextra = mesh_p1%mextra
      mesh%mes_extra = mesh_p1%mes_extra
      mesh%mes_int = mesh_p1%mes_int

!=== New nodes
      mesh%np = mesh_p1%np !=== previous P1 points
      mesh%np = mesh%np + mesh_p1%medge * (type_fe - 1) !=== New points on each edge
      mesh%np = mesh%np + mesh_p1%me * (type_fe - 2)*(type_fe - 1)/2 !=== New points inside each element
      mesh%np = mesh%np + mesh_p1%medges*(type_fe - 1) !=== New ghost points on boundaries

      mesh%dom_np = mesh_p1%dom_np !=== previous P1 points
      mesh%dom_np = mesh%dom_np + mesh_p1%medge * (type_fe - 1) !=== New points on each edge
      mesh%dom_np = mesh%dom_np + mesh_p1%me * (type_fe - 2)*(type_fe - 1)/2 !=== New points inside each element

      mesh%gauss%n_w = (type_fe + 1) * (type_fe + 2) / 2
      mesh%gauss%n_ws = type_fe + 1

!=== Unchanged cell/edges structures
      ALLOCATE(mesh%jce(ne, me))
      mesh%jce = mesh_p1%jce
      ALLOCATE(mesh%neigh(ne, mesh%me))
      mesh%neigh = mesh_p1%neigh
      ALLOCATE(mesh%sides(mesh%mes))
      mesh%sides = mesh_p1%sides
      ALLOCATE(mesh%neighs(mesh%mes))
      mesh%neighs = mesh_p1%neighs
      ALLOCATE(mesh%sides_extra(mesh%mes_extra))
      mesh%sides_extra = mesh_p1%sides_extra
      ALLOCATE(mesh%neighs_extra(mesh%mes_extra))
      mesh%neighs_extra = mesh_p1%neighs_extra
      ALLOCATE(mesh%i_d(mesh%me))
      mesh%i_d = mesh_p1%i_d
      ALLOCATE(mesh%jcc_extra(mesh%mextra))
      mesh%jcc_extra = mesh_p1%jcc_extra


      mesh%nis = mesh_p1%nis
      ALLOCATE(mesh%isolated_interfaces(mesh_p1%nis, 2))
      mesh%isolated_interfaces = mesh_p1%isolated_interfaces

      ALLOCATE(mesh%jce_extra(SIZE(mesh_p1%jce_extra, 1), SIZE(mesh_p1%jce_extra, 2)))
      mesh%jce_extra = mesh_p1%jce_extra
      ALLOCATE(mesh%jees(SIZE(mesh_p1%jees)))
      mesh%jees = mesh_p1%jees
      ALLOCATE(mesh%jecs(SIZE(mesh_p1%jecs)))
      mesh%jecs = mesh_p1%jecs


      CALL mesh%create_comm(mesh_p1%comm)
      CALL mesh%gather_dom_np
      CALL mesh%gather_me
      CALL mesh%gather_medge

!=== New nodes and new global numbering
      ALLOCATE(mesh%isolated_jjs(mesh_p1%nis))
      ALLOCATE(mesh%jjs_int(nws+f_dof_edge, mes_int))

      ALLOCATE(mesh%jj(mesh%gauss%n_w, me), SOURCE=-1)
      ALLOCATE(mesh%jjs(mesh%gauss%n_ws, mes))   

      ALLOCATE(mesh%jj_extra(mesh%gauss%n_w, mesh%mextra)) 
      ALLOCATE(mesh%jjs_extra(mesh%gauss%n_ws, mesh%mes_extra))

      ALLOCATE(mesh%rr(kd, mesh%np))    
      ALLOCATE(mesh%rrs_extra(kd, mesh%gauss%n_w, mesh%mes_extra))

      ALLOCATE(mesh%loc_to_glob(mesh%np)) 
      ALLOCATE(mesh%proc_np_loc(2, mesh%np-mesh%dom_np))


      IF (type_fe==1 .OR. mesh_p1%me == 0) THEN

         mesh%jj = mesh_p1%jj
         mesh%jjs = mesh_p1%jjs
         mesh%jjs_extra = mesh_p1%jjs_extra
         mesh%jj_extra = mesh_p1%jj_extra
         mesh%rr = mesh_p1%rr
         mesh%rrs_extra = mesh_p1%rrs_extra
         mesh%loc_to_glob = mesh_p1%loc_to_glob

         mesh%isolated_jjs = mesh_p1%isolated_jjs
         mesh%proc_np_loc(:,:) = mesh_P1%proc_np_loc(:,:)

         RETURN
      END IF
      
      nb_angle = 0
      ALLOCATE(virgin(mesh_p1%medge), j_edge(nw * f_dof_edge, me), jjs_edge(f_dof_edge, mes), r_mid(kd))

      IF (kd == 3) THEN
         CALL error_petsc('create_iso_grid_distributed: 3D case not programmed yet !')
      END IF

      proc = mesh_p1%proc
      mesh%proc_np_loc = 0
      mesh%proc_np_loc(:,1:SIZE(mesh_P1%proc_np_loc,2)) = mesh_P1%proc_np_loc(:,:)

      !===GENERATION OF THE Pk GRID

!=== Copy the existing P1 grid ===!
      mesh%rr(:, 1:dom_np) = mesh_p1%rr(:, 1:dom_np)
      mesh%rr(:, mesh%dom_np + 1:mesh%dom_np + np - dom_np) = mesh_p1%rr(:, dom_np + 1:)
      mesh%jj(1:nw, :) = mesh_p1%jj
      mesh%jj_extra(1:nw, :) = mesh_p1%jj_extra

      IF (size(mesh_p1%rrs_extra)/=0) THEN
         mesh%rrs_extra(:, 1:nw, 1:mesh_p1%mes_extra) = mesh_p1%rrs_extra !<== JLG April 16, 2026
         !mesh%rrs_extra(:, 1:nw, :) = mesh_p1%rrs_extra
      END IF

      mesh%isolated_jjs = mesh_p1%isolated_jjs &
           + (mesh_p1%disedge(proc) - 1) * (f_dof_edge) + (mesh_p1%discell(proc) - 1) * (f_dof_cell)

      DO m = 1, mesh%me
         DO n = 1, nw
            IF (mesh_p1%jj(n, m) > mesh_p1%dom_np) THEN
               n_loc = mesh_p1%jj(n, m) + mesh_p1%medge * f_dof_edge + mesh_p1%me * f_dof_cell
               mesh%jj(n, m) = n_loc
            END IF
         END DO
      END DO

      mesh%proc_np_loc(:, 1:mesh_p1%np - mesh_p1%dom_np) = mesh_p1%proc_np_loc(:, :)

      DO m = 1, mesh_p1%mextra
         DO n = 1, nw
            p = mesh_p1%get_proc(mesh_p1%jj_extra(n, m), 'np')
            mesh%jj_extra(n, m) = mesh_p1%jj_extra(n, m) &
                 + (mesh_p1%disedge(p) - 1) * f_dof_edge + (mesh_p1%discell(p) - 1) * f_dof_cell
         END DO
      END DO

      DO m = 1, mesh_p1%mes_extra
         DO n = 1, nws
            p = mesh_p1%get_proc(mesh_p1%jjs_extra(n, m), 'np')
            mesh%jjs_extra(n, m) = mesh_p1%jjs_extra(n, m) &
                 + (mesh_p1%disedge(p) - 1) * f_dof_edge + (mesh_p1%discell(p) - 1) * f_dof_cell
         END DO
      END DO
!=== End copy the existing P1 grid ===!


!=== Loop on inner edges, owned by this proc ===!
      virgin = .TRUE.
      n_dof = 0
      DO m = 1, me !===loop on the elements
         DO k = 1, nw !===loop on the nodes (sides) of the element
            edge_g = mesh_p1%jce(k, m)
            edge_l = edge_g - mesh_p1%disedge(proc) + 1

            IF (edge_l <= 0) CYCLE

            n_new_start = (edge_l - 1) * f_dof_edge + mesh_p1%dom_np

            m_op_k = mesh_p1%neigh(k, m)
            n_k1 = MODULO(k, nw) + 1
            n_k2 = MODULO(k + 1, nw) + 1
            n1 = mesh_p1%jj(n_k1, m)
            n2 = mesh_p1%jj(n_k2, m)

            IF (n_k1<n_k2) THEN !===Go from lowest global index to highest global index
               n_start = n1
               n_end = n2
            ELSE
               n_start = n2
               n_end = n1
            END IF

            iso = .FALSE.
            IF (m_op_k == 0) THEN  !===the side is on the boundary
               DO ms = 1, SIZE(mesh_p1%neighs) + 1
                  IF (ms == SIZE(mesh_p1%neighs) + 1) THEN
                     CALL error_petsc("BUG in create_iso_grid_distributed: cell near boundary isnt in neighs "&
                     //to_str(m_op_k)//', '//to_str(m))
                  END IF
                  IF (mesh_p1%neighs(ms) == m) EXIT
               END DO
               !               CALL is_on_curved_interface(mesh_p1%sides(ms), iso, interface)

            END IF

            IF (virgin(edge_l)) THEN !===This side is new
               DO l = 1, f_dof_edge
                  n_dof = n_dof + 1 !===New index created
                  j_edge((k - 1) * f_dof_edge + l, m) = l + n_new_start
                  mesh%rr(:, l + n_new_start) = mesh_p1%rr(:, n_start) &
                       + l * (mesh_p1%rr(:, n_end) - mesh_p1%rr(:, n_start)) / type_fe
                  !                  IF (iso) THEN
                  !                     CALL rescale_to_curved_boundary(mesh%rr(:, l + n_new_start), interface)
                  !                  END IF
               END DO
            ELSE !===the side has been already considered
               mm = m_op_k
               DO i = 1, nw
                  IF (mesh_p1%neigh(i, mm) == m) THEN
                     kk = i
                     EXIT
                  END IF
               ENDDO
               DO l = 1, f_dof_edge
                  j_edge((k - 1) * f_dof_edge + l, m) = j_edge((kk - 1) * f_dof_edge + l, mm) !===New index created
               END DO
            ENDIF
            virgin(edge_l) = .FALSE.
         ENDDO
      ENDDO

      IF (n_dof /= mesh_p1%medge * f_dof_edge) THEN
         CALL error_petsc('create_iso_grid_distributed: n_dof /= mesh_p1%medge * f_dof_edge (I)')
      END IF
      
      
!=== Loop on outer edges, not owned by this proc ===!
      n_dof = 0
      DO edges = 1, mesh_p1%medges
         edge_g = mesh_p1%jees(edges)
         m = mesh_p1%jecs(edges)
         DO k = 1, nw
            IF (mesh_p1%jce(k, m) == edge_g) EXIT
         ENDDO

         p = mesh_p1%get_proc(edge_g, 'medge')

         n_k1 = MODULO(k, nw) + 1
         n_k2 = MODULO(k + 1, nw) + 1
         n1 = mesh_p1%jj(n_k1, m)
         n2 = mesh_p1%jj(n_k2, m)
         IF (n_k1<n_k2) THEN !===Go from lowest global index to highest global index
            n_start = n1
            n_end = n2
         ELSE
            n_start = n2
            n_end = n1
         END IF

         edge_l = edge_g - mesh_p1%disedge(p) + 1
         n_new_start = (edges - 1) * f_dof_edge + mesh_p1%me * f_dof_cell + mesh_p1%medge * f_dof_edge + mesh_p1%np

         DO l = 1, f_dof_edge
            n_dof = n_dof + 1 !===New index created
            j_edge((k - 1) * f_dof_edge + l, m) = l + n_new_start
            mesh%rr(:, n_new_start + l) = mesh_p1%rr(:, n_start) &
                 + l * (mesh_p1%rr(:, n_end) - mesh_p1%rr(:, n_start)) / type_fe
            
            IF ((p==proc) .OR. (n_new_start + l <= mesh%dom_np)) THEN
               CALL error_petsc("BUG in create_iso_grid_distributed: edge on proc boundary should not be owned by this proc: "&
               //to_str(p)//" == "//to_str(proc)//" OR dom_np = "//to_str(mesh%dom_np)//&
               " >= n_new_start + l = "//to_str(n_new_start + l))
            END IF
            mesh%proc_np_loc(1, n_new_start + l - mesh%dom_np) = p
            mesh%proc_np_loc(2, n_new_start + l - mesh%dom_np) = l + (edge_l - 1) * f_dof_edge + mesh_p1%domnp(p)
         END DO

      END DO

      IF (n_dof /= mesh_p1%medges * f_dof_edge) THEN
         CALL error_petsc('create_iso_grid_distributed: n_dof /= mesh_p1%medges * f_dof_edge (II)')
      END IF
      n_dof = 0

      !===connectivity array for iso grid deduced from j_edge, + additional point in P3
      DO m = 1, me
         DO n = 1, f_dof_edge * ne
            mesh%jj(nw + n, m) = j_edge(n, m)
         END DO

         ! IF (type_fe==3) THEN
         !    n_dof = n_dof + 1
         !    n_new_start = m + mesh_p1%dom_np + mesh_p1%medge * f_dof_edge
         !    mesh%jj(10, m) = n_new_start
         !    mesh%rr(:, n_new_start) = &
         !         (mesh_p1%rr(:, mesh_p1%jj(1, m)) + mesh_p1%rr(:, mesh_p1%jj(2, m)) + mesh_p1%rr(:, mesh_p1%jj(3, m))) / type_fe
         ! END IF
         n_loc_dof = 0
         DO n1=1, f_dof_edge
            DO n2=1, f_dof_edge
               DO n3=1, f_dof_edge
                  IF (n1+n2+n3/=type_fe) CYCLE !only retain barycentric coordinates inside the triangle (hence loops starting from 1 and not 0)
                  n_loc_dof = n_loc_dof + 1
                  n_dof = n_dof + 1
                  n_new_start = m + me*(n_loc_dof-1) + mesh_p1%dom_np + mesh_p1%medge * f_dof_edge
            
                  mesh%jj(nw+f_dof_edge*ne + n_loc_dof, m) = n_new_start
                  mesh%rr(:, n_new_start) = &
                     (n1*mesh_p1%rr(:, mesh_p1%jj(1, m)) &
                  +   n2*mesh_p1%rr(:, mesh_p1%jj(2, m)) &
                  +   n3*mesh_p1%rr(:, mesh_p1%jj(3, m))) / type_fe
               END DO
            END DO
         END DO
         IF (n_loc_dof /= f_dof_cell) THEN
            CALL error_petsc('&
            &BUG in create_iso_grid_distributed: type_fe '//to_str(type_fe)//", f_dof_cell="//to_str(f_dof_cell)//&
            &" and yet could find n_loc_dof = "//to_str(n_loc_dof))
         END IF
      END DO

      IF (type_fe >= 3 .AND. n_dof /= mesh_p1%me*f_dof_cell) THEN
         CALL error_petsc('BUG in create_iso_grid_distributed, type_fe >= 3 .and. n_dof /= mesh_p1%me*f_dof_cell')
      END IF
      IF (MINVAL(mesh%jj) < 0) THEN
         CALL local_error_petsc("BUG => apparently all jj were not attributed")
      END IF

      !==connectivity array the surface elements of the iso grid
      DO ms = 1, mes
         m = mesh_p1%neighs(ms)
         DO n = 1, kd + 1 !===Simplices
            IF (.NOT. (ANY(mesh_p1%jj(n, m)==mesh_p1%jjs(:, ms)))) THEN
               kk = n
               EXIT
            END IF
         ENDDO
         DO l = 1, f_dof_edge
            jjs_edge(l, ms) = j_edge((kk - 1) * f_dof_edge + l, m) !===New index created
         END DO
      ENDDO
      mesh%jjs(1:nws, :) = mesh_p1%jjs
      DO i = 1, SIZE(mesh%jjs, 2)
         DO n = 1, nws
            IF (mesh%jjs(n, i) > mesh_p1%dom_np) THEN
               mesh%jjs(n, i) = mesh_p1%jjs(n, i) + mesh_p1%medge * f_dof_edge + mesh_p1%me * f_dof_cell
            END IF
         END DO
      END DO

      mesh%jjs(nws + 1:, :) = jjs_edge

      DEALLOCATE(virgin, j_edge, jjs_edge, r_mid)


      !===Internal surface elements
      DO ms = 1, mes_int
         m = mesh%neighs_int(1, ms)
         !===Finding the corresponding side in the cell
         DO k = 1, nw
            IF (MINVAL(ABS(mesh_p1%jj(k, m) - mesh_p1%jjs_int(:, ms)))/=0) EXIT
         ENDDO

         n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
         IF (n_ks(1)>n_ks(2)) THEN
            n_ks = (/n_ks(2), n_ks(1)/)
         END IF

         mesh%jjs_int(1, ms) = mesh_p1%jj(n_ks(1), m)
         IF (mesh%jjs_int(1, ms) > mesh_p1%dom_np) mesh%jjs_int(1, ms) = mesh%jjs_int(1, ms) + mesh%dom_np - mesh_p1%dom_np
         mesh%jjs_int(3, ms) = mesh_p1%jj(n_ks(2), m)
         IF (mesh%jjs_int(3, ms) > mesh_p1%dom_np) mesh%jjs_int(3, ms) = mesh%jjs_int(3, ms) &
              + mesh%dom_np - mesh_p1%dom_np
         mesh%jjs_int(2, ms) = mesh_p1%jj(k, 4 * (m - 1) + 1)

         !         CALL is_on_curved_interface(mesh_p2%sides_int(ms), iso, interface)
         !         IF (iso) THEN
         !            CALL rescale_to_curved_boundary(mesh%rr(:, mesh%jj(k, 4 * (m - 1) + 1)), interface)
         !         END IF
      ENDDO

      !===new vertices on extra cells
      DO m = 1, mesh%mextra
         DO k = 1, nw !===loop on the nodes (sides) of the element
            edge_g = mesh_p1%jce_extra(k, m)
            p_e = mesh_p1%get_proc(edge_g, 'medge')
            edge_l = edge_g - mesh_p1%disedge(p_e) + 1

            cell_g = mesh_p1%jcc_extra(m)
            p_c = mesh_p1%get_proc(cell_g, 'me')
            cell_l = cell_g - mesh_p1%discell(p_c) + 1

            DO l = 1, f_dof_edge
               mesh%jj_extra(nw + (k - 1) * f_dof_edge + l, m) = l &
                    + (edge_l - 1) * f_dof_edge + mesh_p1%domnp(p_e) + mesh%disp(p_e) - 1
            END DO

            ! IF (type_fe==3) THEN
            !    mesh%jj_extra(10, m) = cell_l + mesh_p1%domedge(p_c) * 2 + mesh_p1%domnp(p_c) + mesh%disp(p_c) - 1
            ! END IF

            n_loc_dof = 0
            DO n1=1, f_dof_edge
               DO n2=1, f_dof_edge
                  DO n3=1, f_dof_edge
                     IF (n1+n2+n3/=type_fe) CYCLE !only retain barycentric coordinates inside the triangle (hence loops starting from 1 and not 0)
                     n_loc_dof = n_loc_dof + 1
                     n_new_start = cell_l + (n_loc_dof-1)*mesh_p1%domcell(p_c) + &
                                   mesh_p1%domedge(p_c) * f_dof_edge + mesh_p1%domnp(p_c) + (mesh%disp(p_c) - 1)
               
                     mesh%jj_extra(nw+f_dof_edge*ne + n_loc_dof, m) = n_new_start
                  END DO
               END DO
            END DO

         END DO
      END DO
      !==connectivity array the surface elements of the iso grid for extras
      DO ms = 1, mesh%mes_extra
         iso = .FALSE.
         !         CALL is_on_curved_interface(mesh%sides_extra(ms), iso, interface)

         cell_g = mesh%neighs_extra(ms)
         DO m = 1, mesh%mextra !find associated extra cell
            IF (mesh_p1%jcc_extra(m) == cell_g) EXIT
         END DO

         kk = 0
         DO n = 1, kd + 1 !===find side in cell
            IF (.NOT. (ANY(mesh_p1%jj_extra(n, m)==mesh_p1%jjs_extra(:, ms)))) THEN
               kk = n
               EXIT
            END IF
         ENDDO
         IF (kk==0) THEN
            CALL local_error_petsc("BUG in create_iso_grid_distributed: did not find extra cell")
         END IF
         DO l = 1, f_dof_edge
            mesh%jjs_extra(nws + l, ms) = mesh%jj_extra(nw + (kk - 1) * f_dof_edge + l, m)
         END DO

         DO k = 1, kd + 1
            n_k1 = MODULO(k, nw) + 1
            n_k2 = MODULO(k + 1, nw) + 1
            IF (n_k1<n_k2) THEN !===Go from lowest index to highest index
               n_start = n_k1
               n_end = n_k2
            ELSE
               n_start = n_k2
               n_end = n_k1
            END IF

            DO l = 1, f_dof_edge
               mesh%rrs_extra(:, nw + (k - 1) * f_dof_edge + l, ms) = mesh_p1%rrs_extra(:, n_start, ms) &
                    + l * (mesh_p1%rrs_extra(:, n_end, ms) - mesh_p1%rrs_extra(:, n_start, ms)) / type_fe
               !               IF (iso) THEN
               !                  CALL rescale_to_curved_boundary(mesh%rrs_extra(:, nw + (k - 1) * f_dof + l, ms), interface)
               !               END IF
            END DO
         END DO

         ! IF (type_fe==3) THEN
         !    mesh%rrs_extra(:, 10, ms) = &
         !         (mesh_p1%rrs_extra(:, 1, ms) + mesh_p1%rrs_extra(:, 2, ms) + mesh_p1%rrs_extra(:, 3, ms)) / 3
         ! END IF

         n_loc_dof = 0
         DO n1=1, f_dof_edge
            DO n2=1, f_dof_edge
               DO n3=1, f_dof_edge
                  IF (n1+n2+n3/=type_fe) CYCLE !only retain barycentric coordinates inside the triangle (hence loops starting from 1 and not 0)
                  n_loc_dof = n_loc_dof + 1
            
                  mesh%rrs_extra(:,nw+f_dof_edge*ne + n_loc_dof, ms) = &
                     (n1*mesh_p1%rrs_extra(:, 1, ms) &
                  +   n2*mesh_p1%rrs_extra(:, 2, ms) &
                  +   n3*mesh_p1%rrs_extra(:, 3, ms)) / type_fe
               END DO
            END DO
         END DO

      ENDDO

      CALL mesh%build_loc_to_glob

   END SUBROUTINE  create_iso_grid_distributed

   SUBROUTINE refinement_iso_grid_distributed(mesh_p1_in, mesh_p1_out, refine_factor)
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type)     :: mesh_p1_in, mesh_pk, mesh_p1_out
      INTEGER, INTENT(IN) :: refine_factor

      IF (refine_factor<2) THEN
         CALL copy_mesh(mesh_p1_in, mesh_p1_out)
         RETURN
      END IF

      CALL create_iso_grid_distributed(mesh_p1_in, mesh_pk, refine_factor)
      CALL general_refinement_iso_grid_distributed(mesh_pk, mesh_p1_out)
      CALL free_mesh(mesh_pk)
   END SUBROUTINE refinement_iso_grid_distributed

   SUBROUTINE general_refinement_iso_grid_distributed(mesh_pk, mesh)
      ! subroutine that takes iso_grid_distributed mesh_pk as argument (meant for Pk)
      ! returns mesh which is the P1 mesh having same nodes as Pk but connected (many more triangles)
      
      !===jj(:, :)    nodes of the  volume_elements of the input grid
      !===jjs(:, :)    nodes of the surface_elements of the input grid
      !===rr(:, :)    cartesian coordinates of the nodes of the input grid
      !===m_op(:,:)   volume element opposite to each node
      !===neigh_el(:) volume element ajacent to the surface element
      !===jj_f(:, :)  nodes of the  volume_elements of the output p2 grid
      !===jjs_f(:, :)  nodes of the surface_elements of the output p2 grid
      !===rr_f(:, :)  cartesian coordinates of the nodes of the output p2 grid
      USE def_type_mesh
      USE my_util, ONLY: error_petsc, to_str, local_error_petsc
      USE tab_subcell_iso_refinement_module
      IMPLICIT NONE
      TYPE(mesh_type)                        :: mesh_pk, mesh
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: raw_jj_HL, raw_jjs_HL
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: jj_HL, jjs_HL
      INTEGER :: ms, i, p_c, m_new, ne, nc, me, mes, nw, nws, kd, n, m, k, mextra
      INTEGER, DIMENSION(2) :: n_ks, pts_jjs
      INTEGER, DIMENSION(3) :: edges_g, p_es, nodes_new
      INTEGER, DIMENSION(4) :: neigh_coarse_triangles
      INTEGER :: proc, nb_proc, p, cell_g, cell_l, nb_subcell, nb_subcells, nc_neigh, cur_local_edge, refine_factor
      INTEGER :: m1, interface, m_center, mes_int, n_loc_1, n_loc_2, offset, idx
      INTEGER :: sub_edge, old_e_g, idx1, idx2, cur_local_edge_offset
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: outer_edges
      INTEGER, DIMENSION(:,:),   ALLOCATABLE :: inner_edges, inner_edges_to_subcell, outer_edges_coarse
      INTEGER :: nc_complementary, m_new_complementary, n_complementary, idx_m_new_complementary
      INTEGER, DIMENSION(:),     ALLOCATABLE :: idx_n_loc_1, idx_n_loc_2, tab_subcell_order, tab_subcells_order
      LOGICAL :: iso, is_inner_edge

      !=== dummy to avoid warning (TODO => add curved boundary)
      iso = .FALSE.; interface = -1
      IF (iso) write(*,*) sgn(1.d0)
      !=== end dummy

      refine_factor = SIZE(mesh_pk%jjs, 1)-1

      nb_subcell = refine_factor**2
      nb_subcells = refine_factor
      ne = 3
      ALLOCATE(jj_HL(nb_subcell, ne))
      ALLOCATE(jjs_HL(nb_subcells, ne-1))
      ALLOCATE(tab_subcell_order(refine_factor**2))
      ALLOCATE(tab_subcells_order(refine_factor))

      tab_subcell_order = [(n, n=1, nb_subcell)]
      tab_subcells_order = [(n, n=1, nb_subcells)]

      SELECT CASE(refine_factor)
      CASE(2)
         raw_jj_HL  = jj_HL_2
         raw_jjs_HL = jjs_HL_2

      CASE(3)
         raw_jj_HL  = jj_HL_3
         raw_jjs_HL = jjs_HL_3

      CASE(4)
         raw_jj_HL  = jj_HL_4
         raw_jjs_HL = jjs_HL_4

      CASE DEFAULT
         CALL build_jj_HL (refine_factor, raw_jj_HL)
         CALL build_jjs_HL(refine_factor, raw_jjs_HL)
         ! CALL error_petsc("BUG in refinement_iso_grid_distributed => order = "//to_str(refine_factor)//&
         ! &" does not have its subcell order defined")
      END SELECT

      IF (SIZE(jj_HL, 1)/=SIZE(raw_jj_HL, 1)) THEN
         CALL error_petsc("BUG in refinement_iso_grid_distributed: references jj_HL seem wrong for "//to_str(refine_factor)&
         &//": should have "//to_str(nb_subcell)//" subcell, not "//to_str(SIZE(raw_jj_HL, 1)))
      ELSE IF (SIZE(jjs_HL, 1)/=SIZE(raw_jjs_HL, 1)) THEN
         CALL error_petsc("BUG in refinement_iso_grid_distributed: references jjs_HL seem wrong for "//to_str(refine_factor)&
         &//": should have "//to_str(nb_subcells)//" subcell, not "//to_str(SIZE(raw_jjs_HL, 1)))
      END IF

      DO n=1, nb_subcell
         jj_HL(tab_subcell_order(n), :) = raw_jj_HL(n, :)
      END DO
      DO n=1, nb_subcells
         jjs_HL(tab_subcells_order(n), :) = raw_jjs_HL(n, :)
      END DO

      CALL build_outer_inner_edges(jj_HL, jjs_HL)

      IF (mesh_pk%me == 0) RETURN

      CALL mesh%info%copy(mesh_pk%info)
      nw  = 3   !===nodes in each volume element (3 in 2D)
      nws = 2  !===nodes in each edge element (2 in 2D)
      me = SIZE(mesh_pk%jj, 2)   !===number of cells
      kd = SIZE(mesh_pk%rr, 1)   !===space dimensions
      mes     = SIZE(mesh_pk%jjs, 2)
      mes_int = SIZE(mesh_pk%jjs_int, 2)
      nb_proc = SIZE(mesh_pk%domnp)

      mesh%me      = refine_factor**2 * mesh_pk%me
      mesh%mes     = refine_factor * mesh_pk%mes !---> something with news to take into account
      mesh%mes_int = refine_factor * mesh_pk%mes_int
      mesh%np      = mesh_pk%np
      mesh%medge   = refine_factor * mesh_pk%medge + 3*(refine_factor-1)*refine_factor/2*mesh_pk%me
      mesh%medges  = refine_factor * mesh_pk%medges
      !naive mextra, can definitely be reduced with some work
      mesh%mextra = refine_factor**2 * mesh_pk%mextra
      !however mes_extra cannot be reduced
      mesh%mes_extra = refine_factor * mesh_pk%mes_extra

      ALLOCATE(mesh%jj(nw, mesh%me)) !--->done
      ALLOCATE(mesh%jjs(nws, mesh%mes))  !--->done
      ALLOCATE(mesh%rr(kd, mesh%np)) !--->done
      ALLOCATE(mesh%loc_to_glob(mesh%np)) !--->done


      ALLOCATE(mesh%jce(ne, mesh%me), SOURCE=-1)  !--->done
      !ALLOCATE(mesh%jev(nw - 1, mesh%medge)) !----->never needed so not constructed
      ALLOCATE(mesh%jees(mesh%medges), SOURCE=0)  !--->done
      ALLOCATE(mesh%jecs(mesh%medges), SOURCE=0)  !--->done

      ALLOCATE(mesh%neigh(ne, mesh%me), SOURCE=-10000) !--->done
      ALLOCATE(mesh%sides(mesh%mes))   !---> still don't know what this is ?? done ?
      ALLOCATE(mesh%neighs(mesh%mes))  !--->done
      ALLOCATE(mesh%i_d(mesh%me)) !--->done

      ALLOCATE(mesh%jjs_int(nws, mesh%mes_int))  !--->done
      ALLOCATE(mesh%neighs_int(2, mesh%mes_int)) !--->done
      ALLOCATE(mesh%sides_int(mesh%mes_int))

      mesh%dom_np = mesh_pk%dom_np
      CALL mesh%create_comm(mesh_pk%comm)
      CALL mesh%gather_dom_np
      CALL mesh%gather_me
      CALL mesh%gather_medge

      mesh%nis = mesh_pk%nis
      ALLOCATE(mesh%isolated_interfaces(mesh_pk%nis, 2))
      mesh%isolated_interfaces = mesh_pk%isolated_interfaces
      ALLOCATE(mesh%isolated_jjs(mesh_pk%nis))

      IF (kd == 3) THEN
         CALL error_petsc('refinement_iso_grid_distributed: 3D case not programmed yet !')
      END IF   

      proc = mesh_pk%get_proc(mesh_pk%loc_to_glob(1), 'np')

      !===GENERATION OF THE Pk GRID
      mesh%rr = mesh_pk%rr

      ALLOCATE(mesh%proc_np_loc(2, mesh%np-mesh%dom_np))
      mesh%proc_np_loc = mesh_pk%proc_np_loc
      mesh%isolated_jjs = mesh_pk%isolated_jjs

      DO m=1, me
         m_center = nb_subcell * (m - 1) + 1
         DO i=1, nb_subcell
            mesh%jj(:, m_center + i - 1) = mesh_pk%jj(jj_HL(i, :), m)
         END DO
      END DO

      ! Loop on coarse triangles
      DO m = 1, me
         
         cur_local_edge_offset = 3*(refine_factor-1)*refine_factor/2*(m-1) + 1
         cur_local_edge = cur_local_edge_offset

         neigh_coarse_triangles(1) = m
         neigh_coarse_triangles(2:4) = mesh_pk%neigh(:, m)

         !define i_d of all subcells of the coarse triangle as the i_d of the coarse triangle
         mesh%i_d(refine_factor**2 * (m - 1) + [(nc, nc=1, nb_subcell)]) = mesh_pk%i_d(m)         

         !Loop on subcells of the coarse triangle   
         DO nc = 1, nb_subcell
            m_new = refine_factor**2 * (m - 1) + nc
            nodes_new = mesh_pk%jj(jj_HL(nc, :), m)
               

            ! Loop on summits/vertices of the subcell
            DO n=1, ne
               ! n_ks(dim=2)      => complementary of n in [1,2,3]
               ! n_loc_1, n_loc_2 => numbering of n_ks(1), n_ks(2) at the coarse cell scale
               ! idx1, idx2       => see subroutine find_edge
               CALL get_edge(n, nc, n_ks, n_loc_1, n_loc_2, idx1, idx2, is_inner_edge)
               ! Find the neighbour of the coarse triangle sharing the edge of the subcell
               IF (is_inner_edge) THEN
                  i = 1
               ELSE
                  DO i = 2, 4
                     IF (neigh_coarse_triangles(i) <= 0) CYCLE
                     IF (ANY(nodes_new(n_ks(1))==mesh_pk%jj(:,neigh_coarse_triangles(i)))) THEN
                        IF (ANY(nodes_new(n_ks(2))==mesh_pk%jj(:,neigh_coarse_triangles(i)))) EXIT
                     END IF
                  END DO
               END IF
               ! Find the subcell in the neighbour coarse triangle sharing the edge of the subcell
               IF (i==5) THEN
                  !Case where neigh is <= 0
                  mesh%neigh(n, m_new) = mesh_pk%neigh(idx2, m)
               ELSE
                  DO nc_neigh=1, nb_subcell
                     IF (ANY(nodes_new(n_ks(1))==mesh_pk%jj(jj_HL(nc_neigh, :), neigh_coarse_triangles(i)))) THEN
                        IF (ANY(nodes_new(n_ks(2))==mesh_pk%jj(jj_HL(nc_neigh, :), neigh_coarse_triangles(i)))) THEN
                           IF (i>1) THEN
                              ! neighboor is outer cell => found our guy
                              EXIT
                           ELSE
                              ! special case where neighboor is the same cell => make sure we don't take the same subcell
                              IF (.NOT. ANY(nodes_new(n)==mesh_pk%jj(jj_HL(nc_neigh, :), neigh_coarse_triangles(i)))) THEN
                                 EXIT
                              END IF
                           END IF
                        END IF
                     END IF
                  END DO
                  !found our neigh
                  mesh%neigh(n, m_new) = refine_factor**2 * (neigh_coarse_triangles(i) - 1) + nc_neigh
               END IF

               !build our edge connectivity
               IF (is_inner_edge) THEN !---> new inner edge inside the coarse cell
                  IF (inner_edges_to_subcell(idx1, 1) == nc) THEN
                     nc_complementary = inner_edges_to_subcell(idx1, 2)
                  ELSE
                     nc_complementary = inner_edges_to_subcell(idx1, 1)
                  END IF
                  m_new_complementary = m_new - nc + nc_complementary
                  idx_n_loc_1 = MINLOC(ABS(jj_HL(nc_complementary,:)-n_loc_1))
                  idx_n_loc_2 = MINLOC(ABS(jj_HL(nc_complementary,:)-n_loc_2))
                  DO n_complementary=1, ne
                     IF (n_complementary /= idx_n_loc_1(1)) THEN
                        IF (n_complementary /= idx_n_loc_2(1)) EXIT                        
                     END IF
                  END DO
!====== DEBUGGING ======!
! write(*,*) 'the two adjacent subcells', nc_complementary, nc
! write(*,*) 'idx of inner triangle:', n_loc_1, n_loc_2, jj_HL(nc, n), idx1!, m_new - nc
! write(*,*) 'n_complementary vs th = ', n_complementary, idx_n_loc_1(1), idx_n_loc_2(1), MOD(m_new, refine_factor**2)
! write(*,*) mesh%jce(n, m_new), mesh%jce(n, m_new), mesh%jce(n_complementary, m_new_complementary)
!====== DEBUGGING ======!
                  DEALLOCATE(idx_n_loc_1, idx_n_loc_2)
                  IF (mesh%jce(n_complementary, m_new_complementary) /= -1) THEN
                     mesh%jce(n, m_new) = mesh%jce(n_complementary, m_new_complementary)
                  ELSE
                     mesh%jce(n, m_new) = cur_local_edge + (mesh%disedge(proc) - 1)     & ! offset for edges on previous procs
                                        + refine_factor * mesh_pk%medge  ! offset for edges resulting from previously existing coarse edges
                     cur_local_edge = cur_local_edge + 1
                  END IF
!====== DEBUGGING ======!
! write(*,*) mesh%jce(n, m_new), mesh%jce(n, m_new)
! IF (mesh%jce(n, m_new) /= mesh%jce(n, m_new)) STOP
! if (mesh%jce(n, m_new) /= mesh%jce(n, m_new)) THEN
!    write(*,*) 'm = ', m, 'n = ', n
!    write(*,*) 'mesh%jce', mesh%jce(:, m_new)
!    write(*,*) ' mesh%jce',  mesh%jce(:, m_new)
!    call error_petsc('wrong in inner edge')
! ! else
!    ! write(*,*) 'inner jce passed'
! end if
!=========== DEBUGGING =========! 
               ELSE                  !---> new outer cell, i.e resulting from partitioning of older triangle edges
                  old_e_g = mesh_pk%jce(idx2, m)
                  IF (old_e_g < mesh_pk%disedge(proc)) THEN
                     p = mesh_pk%get_proc(old_e_g, 'medge')
                  ELSE
                     p = proc
                  END IF
                  IF ((n_loc_1==outer_edges_coarse(idx2, 1))) THEN ! i.e one of the point is one of the old coarse triangle summits
                     sub_edge = 0
                  ELSE IF ((n_loc_1==outer_edges_coarse(idx2, 2))) THEN
                     sub_edge = refine_factor - 1
                  !dummy if refinement of order 2
                  ELSE
                     sub_edge = n_loc_2 - outer_edges_coarse(idx2, 3)
                  END IF
                  mesh%jce(n, m_new) = old_e_g - mesh_pk%disedge(p) + 1 & !old local numbering
                                     + mesh%disedge(p) - 1 & ! offset by all previous edges for new global numbering
                                     + sub_edge * mesh_pk%domedge(p) ! offset related to local numbering at current coarse triangle scale

                  ! IF (old_e_g < mesh_pk%disedge(proc)) THEN 
                  IF (p < proc) THEN
                     DO ms=1, mesh_pk%medges
                        IF (mesh_pk%jees(ms) == old_e_g) EXIT
                     END DO
                     mesh%jees(ms + mesh_pk%medges*sub_edge) = mesh%jce(n, m_new)
                     mesh%jecs(ms + mesh_pk%medges*sub_edge) = m_new
                  END IF

!=========== DEBUGGING =========! 
! if (mesh%jce(n, m_new) /= mesh%jce(n, m_new)) THEN
!    write(*,*) 'm = ', m, 'n = ', n, 'domedge = ', mesh_pk%domedge(p) 
!    write(*,*) 'old values = ', mesh_pk%jce(idx2, m), mesh_pk%jce(:, m)
!    write(*,*) 'mesh%jce', mesh%jce(:, m_new)
!    write(*,*) ' mesh%jce',  mesh%jce(:, m_new)
!    write(*,*) 'sub_edge = ', sub_edge, ', th_sub_edge = ', &
!    (mesh%jce(n, m_new)-(old_e_g - mesh_pk%disedge(p) + 1 + mesh%disedge(p) - 1))/mesh_pk%domedge(p)
!    write(*,*) 'MOD == 0: ', MOD(mesh%jce(n, m_new)-(old_e_g - mesh_pk%disedge(p) + 1 + mesh%disedge(p) - 1),&
!     mesh_pk%domedge(p))
!    write(*,*) 'n_loc = ', n_loc_1, n_loc_2, outer_edges_coarse(idx2, 3)
!    call error_petsc('wrong in outer edge')
! else
!    write(*,*) 'outer jce passed'
! end if
!=========== DEBUGGING =========! 
               END IF
            END DO !end do on vertices of the subcell
         END DO !end do on subcells of the coarse triangle
      END DO !end do on coarse triangles

      DO ms = 1, mes
         DO nc = 1, nb_subcells
            mesh%jjs(:, ms + (nc - 1) * mes) = mesh_pk%jjs(jjs_HL(nc, :), ms)
         ENDDO
      END DO

      DO ms=1, mes
         m = mesh_pk%neighs(ms)
         m_new = refine_factor**2 * (m - 1)
         DO nc = 1, nb_subcells
            mesh%sides(ms + mes*(nc-1)) = mesh_pk%sides(ms)
            pts_jjs = mesh%jjs(:, ms + mes*(nc-1))
            ! Find subcell containing new neighs
            DO idx=1, nb_subcell
               IF (ANY(pts_jjs(1)==mesh%jj(:, m_new+idx))) THEN
                  IF (ANY(pts_jjs(2)==mesh%jj(:, m_new+idx))) EXIT
               END IF
            END DO
            IF (idx == nb_subcell+1) THEN
               CALL local_error_petsc("BUG in surface points => did not find corresponding subcell (I)")
            END IF
            mesh%neighs(ms + mes*(nc-1)) = m_new+idx
         END DO
      END DO

      !===Internal surface elements
      DO ms = 1, mes_int
         m = mesh_pk%neighs_int(1, ms)
         !===Finding the corresponding side in the cell
         DO k = 1, nw
            IF (MINVAL(ABS(mesh_pk%jj(k, m) - mesh_pk%jjs_int(:, ms)))/=0) EXIT
         ENDDO

         n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
         IF (n_ks(1)>n_ks(2)) THEN
            n_ks = (/n_ks(2), n_ks(1)/)
         END IF

         ! VB 27/05/2026 => no test or guarantee that it works, i also don't know where these are used
         DO i=1, refine_factor
            mesh%jjs_int(:, ms + (i - 1) * mes) = mesh_pk%jjs_int(jjs_HL(i, :), ms)
         ENDDO
         ! mesh%jjs_int(1, ms) = mesh_p1%jj(n_ks(1), m)
         ! IF (mesh%jjs_int(1, ms) > mesh_p1%dom_np) mesh%jjs_int(1, ms) = mesh%jjs_int(1, ms) + mesh%dom_np - mesh_p1%dom_np
         ! mesh%jjs_int(1, mes_int + ms) = mesh_p1%jj(n_ks(2), m)
         ! IF (mesh%jjs_int(1, mes_int + ms) > mesh_p1%dom_np) mesh%jjs_int(1, mes_int + ms) = mesh%jjs_int(1, mes_int + ms) &
         !      + mesh%dom_np - mesh_p1%dom_np
         ! mesh%jjs_int(2, ms) = mesh%jj(k, 4 * (m - 1) + 1)
         ! mesh%jjs_int(2, mes_int + ms) = mesh%jj(k, 4 * (m - 1) + 1)
         mesh%neighs_int(1, ms) = mesh%neigh(n_ks(1), 4 * (m - 1) + 1)
         mesh%neighs_int(1, mes_int + ms) = mesh%neigh(n_ks(2), 4 * (m - 1) + 1)
         mesh%sides_int(ms) = mesh_pk%sides_int(ms)
         mesh%sides_int(mes_int + ms) = mesh_pk%sides_int(ms)
         ! VB 27/05/2026 => see also create_iso_grid_distributed, no test or guarantee that it works, i also don't know where these are used

         !         CALL is_on_curved_interface(mesh_p2%sides_int(ms), iso, interface)
         !         IF (iso) THEN
         !            CALL rescale_to_curved_boundary(mesh%rr(:, mesh%jj(k, 4 * (m - 1) + 1)), interface)
         !         END IF
      ENDDO

      !===Counting number of new extra cells
      ! mesh%mextra = 0
      !===Need to rework that to do it the smart way and the update conditions when constructing cells
      !      DO m = 1, mesh_p1%mextra
      !         a = 0
      !         DO k = 1, nw
      !            IF (mesh_p1%disedge(proc) <= mesh_p1%jce_extra(k, m) &
      !                 .and. mesh_p1%jce_extra(k, m) < mesh_p1%disedge(proc + 1) .and. a == 0) THEN
      !               a = 1
      !               mesh%mextra = mesh%mextra + 1
      !            END IF
      !         END DO
      !         DO k = 1, nw
      !            n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
      !            IF (n_ks(1)>n_ks(2)) THEN
      !               n_ks = (/n_ks(2), n_ks(1)/)
      !            END IF
      !
      !            IF (mesh_p1%jj_extra(k, m) < mesh_p1%disp(proc + 1) .and. &
      !                 (mesh_p1%disp(proc) <= mesh_p1%jj_extra(k, m)  .or. &
      !                      (mesh_p1%disedge(proc) <= mesh_p1%jce_extra(n_ks(1), m) &
      !                           .and. mesh_p1%jce_extra(n_ks(1), m) < mesh_p1%disedge(proc + 1)) .or.  &
      !                      (mesh_p1%disedge(proc) <= mesh_p1%jce_extra(n_ks(2), m) &
      !                           .and. mesh_p1%jce_extra(n_ks(2), m) < mesh_p1%disedge(proc + 1)))) THEN
      !               mesh%mextra = mesh%mextra + 1
      !            END IF
      !         END DO
      !      END DO
      !===In the mean time naive mesh%mextra = mesh_pk%mextra * refine_factor**2
      ALLOCATE(mesh%jj_extra(nw, mesh%mextra)) !---->
      ALLOCATE(mesh%jce_extra(nw, mesh%mextra), SOURCE=-1) !---->
      ALLOCATE(mesh%jcc_extra(mesh%mextra)) !---->

      mextra = 0
      DO m=1, mesh_pk%mextra
         DO i=1, refine_factor**2
            mextra = mextra + 1
            mesh%jj_extra(:, mextra) = mesh_pk%jj_extra(jj_HL(i, :), m)
         END DO
      END DO

      mextra = 0
      DO m = 1, mesh_pk%mextra
         !!! getting global and local numbering of old coarse extra cell
         cell_g = mesh_pk%jcc_extra(m)
         p_c = mesh_pk%get_proc(cell_g, 'me')
         cell_l = cell_g - mesh_pk%discell(p_c) + 1
         !!! getting global and local numbering of old coarse extra cell
         DO nc = 1, nb_subcell
            mextra = mextra + 1
            mesh%jcc_extra(mextra) = mesh%discell(p_c) - 1     & !offset due to global numbering
                                   + nb_subcell*(cell_l-1)+1 & !offset due to cells owned by proc p_c
                                   + nc - 1                    !renumbering at the scale of the new fine cell
         END DO
      END DO

      mextra = 0
      DO m = 1, mesh_pk%mextra
         p_c = mesh_pk%get_proc(mesh_pk%jcc_extra(m), 'me')
         !!! WARNING => jcc_extra is the GLOBAL numbering of the elements: need -mesh_pk%discell(p_c) to make it local
         cur_local_edge_offset = 3*(refine_factor-1)*refine_factor/2*(mesh_pk%jcc_extra(m) - mesh_pk%discell(p_c)) + 1
         cur_local_edge = cur_local_edge_offset
         !!! WARNING => jcc_extra is the GLOBAL numbering of the elements: need -mesh_pk%discell(p_c) to make it local

         edges_g = mesh_pk%jce_extra(:, m)
         DO i = 1, ne
            p_es(i) = mesh_pk%get_proc(edges_g(i), 'medge')
         END DO

         DO nc = 1, nb_subcell
            m_new = refine_factor**2 * (mesh_pk%jcc_extra(m) - 1) + nc

            mextra = mextra + 1
            DO n = 1, ne
               ! n_ks(dim=2)      => complementary of n in [1,2,3]
               ! n_loc_1, n_loc_2 => numbering of n_ks(1), n_ks(2) at the coarse cell scale
               ! idx1, idx2       => see subroutine find_edge
               CALL get_edge(n, nc, n_ks, n_loc_1, n_loc_2, idx1, idx2, is_inner_edge)
               IF (is_inner_edge) THEN !---> new inner edge inside the coarse cell
                  IF (inner_edges_to_subcell(idx1, 1) == nc) THEN
                     nc_complementary = inner_edges_to_subcell(idx1, 2)
                  ELSE
                     nc_complementary = inner_edges_to_subcell(idx1, 1)
                  END IF
                  m_new_complementary = m_new - nc + nc_complementary
                  DO idx_m_new_complementary=1, SIZE(mesh%jcc_extra) 
                     IF (mesh%jcc_extra(idx_m_new_complementary) == m_new_complementary) EXIT
                  END DO
                  idx_n_loc_1 = MINLOC(ABS(jj_HL(nc_complementary,:)-n_loc_1))
                  idx_n_loc_2 = MINLOC(ABS(jj_HL(nc_complementary,:)-n_loc_2))
                  DO n_complementary=1, ne
                     IF (n_complementary /= idx_n_loc_1(1)) THEN
                        IF (n_complementary /= idx_n_loc_2(1)) EXIT                        
                     END IF
                  END DO
!====== DEBUGGING ======!
! write(*,*) 'the two adjacent subcells', nc_complementary, nc
! write(*,*) 'idx of inner triangle:', n_loc_1, n_loc_2, jj_HL(nc, n), idx1!, m_new - nc
! write(*,*) 'n_complementary vs th = ', n_complementary, idx_n_loc_1(1), idx_n_loc_2(1), MOD(m_new, refine_factor**2)
! write(*,*) mesh%jce(n, m_new), mesh%jce(n, m_new), mesh%jce(n_complementary, m_new_complementary)
!====== DEBUGGING ======!
                  DEALLOCATE(idx_n_loc_1, idx_n_loc_2)
                  IF (mesh%jce_extra(n_complementary, idx_m_new_complementary) /= -1) THEN
                     mesh%jce_extra(n, mextra) = mesh%jce_extra(n_complementary, idx_m_new_complementary)
                  ELSE
                     mesh%jce_extra(n, mextra) = cur_local_edge + (mesh%disedge(p_c) - 1)     & ! offset for edges on previous procs
                                        + refine_factor * mesh_pk%domedge(p_c)  ! offset for edges resulting from previously existing coarse edges
                     cur_local_edge = cur_local_edge + 1
                  END IF
!====== DEBUGGING ======!
! write(*,*) mesh%jce(n, m_new), mesh%jce(n, m_new)
! IF (mesh%jce(n, m_new) /= mesh%jce(n, m_new)) STOP
! if (mesh%jce(n, m_new) /= mesh%jce(n, m_new)) THEN
!    write(*,*) 'm = ', m, 'n = ', n
!    write(*,*) 'mesh%jce', mesh%jce(:, m_new)
!    write(*,*) ' mesh%jce',  mesh%jce(:, m_new)
!    call error_petsc('wrong in inner edge')
! ! else
!    ! write(*,*) 'inner jce passed'
! end if
! write(*,*) 'bug in inner edges', mesh%jce_extra(n, mextra), mesh%jce_extra(n, mextra)
! write(*,*) 'domedge P2 = ', mesh_pk%domedge
! write(*,*) 'disedge refined = ', mesh%disedge
! write(*,*) 'p_c = ', p_c
! write(*,*) 'cur_local_edge = ', cur_local_edge
! write(*,*) 'already defined before = ', mesh%jce_extra(n_complementary, idx_m_new_complementary) /= -1
! write(*,*) 'cur triangle  = ', n, n_loc_1, n_loc_2, nc
! write(*,*) 'complementary triangle = ', nc_complementary
! write(*,*) 'complementary edge = ', n_complementary, n
! if (mesh%jce_extra(n, mextra) /= mesh%jce_extra(n, mextra)) THEN
!    stop
! end if
!=========== DEBUGGING =========! 
               ELSE                  !---> new outer cell, i.e resulting from partitioning of older triangle edges
                  old_e_g = mesh_pk%jce_extra(idx2, m)
                  IF ((n_loc_1==outer_edges_coarse(idx2, 1))) THEN ! i.e one of the point is one of the old coarse triangle summits
                     sub_edge = 0
                  ELSE IF ((n_loc_1==outer_edges_coarse(idx2, 2))) THEN
                     sub_edge = refine_factor - 1
                  !dummy if refinement of order 2
                  ELSE
                     sub_edge = n_loc_2 - outer_edges_coarse(idx2, 3)
                  END IF
                  mesh%jce_extra(n, mextra) = old_e_g - mesh_pk%disedge(p_es(idx2)) + 1 & !old local numbering
                                     + mesh%disedge(p_es(idx2)) - 1 & ! offset by all previous edges for new global numbering
                                     + sub_edge * mesh_pk%domedge(p_es(idx2)) ! offset related to local numbering at current coarse triangle scale

!=========== DEBUGGING =========! 
! if (mesh%jce_extra(n, mextra) /= mesh%jce_extra(n, mextra)) THEN
!    write(*,*) 'bug in outer edges', mesh%jce_extra(n, mextra), mesh%jce_extra(n, mextra)
! end if
!=========== DEBUGGING =========! 
               END IF
            END DO !end do edge of subcell
         END DO !end do subcell of coarse triangle
      END DO !end do coarse extra triangles

      !===Constructing the extra cells at interfaces
      ALLOCATE(mesh%jjs_extra(nws, mesh%mes_extra))
      ALLOCATE(mesh%rrs_extra(2, nw, mesh%mes_extra))
      ALLOCATE(mesh%sides_extra(mesh%mes_extra), mesh%neighs_extra(mesh%mes_extra))

      mextra = 0
      DO m=1, mesh_pk%mes_extra
         DO i=1, refine_factor
            mextra = mextra + 1
            mesh%jjs_extra(:, mextra) = mesh_pk%jjs_extra(jjs_HL(i, :), m)
         END DO
      END DO

      mextra = 0
      DO m = 1, mesh_pk%mes_extra
         !global and local numbering of coarse extra cell
         cell_g = mesh_pk%neighs_extra(m)
         p_c = mesh_pk%get_proc(cell_g, 'me')
         cell_l = cell_g - (mesh_pk%discell(p_c) - 1)

         !find associated idx of cell_g/cell_l inside jcc_extra
         DO m1 = 1, mesh_pk%mextra 
            IF (mesh_pk%jcc_extra(m1) == cell_g) EXIT
         END DO

         ! Loop on new subcells
         DO nc = 1, nb_subcells
            mextra = mextra + 1
            mesh%sides_extra(mextra) = mesh_pk%sides_extra(m)
            pts_jjs = mesh%jjs_extra(:, mextra)
            ! Find subcell containing new neighs
            DO idx=1, nb_subcell
               IF (ANY(pts_jjs(1)==mesh%jj_extra(:, refine_factor**2*(m1-1)+idx))) THEN
                  IF (ANY(pts_jjs(2)==mesh%jj_extra(:, refine_factor**2*(m1-1)+idx))) EXIT
               END IF
            END DO
            IF (idx == nb_subcell+1) THEN
               CALL local_error_petsc("BUG in surface points => did not find corresponding subcell (II)")
            END IF
            mesh%neighs_extra(mextra) = refine_factor**2*(cell_l-1)+idx & !local numbering on proc containing this point
                                       + (mesh%discell(p_c) - 1)! shift for global numbering

            mesh%rrs_extra(:, :, mextra) = mesh_pk%rrs_extra(:, jj_HL(idx,:), m)
!========= DEBUGGING ========!
! if (mesh%neighs_extra(mextra) /= mesh%neighs_extra(mextra)) THEN
   ! write(*,*) mextra, mesh%neighs_extra(mextra), mesh%neighs_extra(mextra)
   ! call error_petsc('wrong neighs_extra')
! end if
!========= DEBUGGING ========!
         END DO !end do nb_subcells
      END DO !end do coarse mes_extra

      !jjs_extra !(extra layer of cells not own by proc but with dofs own by proc)
      !rrs_extra  ! coordinates for cells at interfaces
      !sides_extra, neighs_extra !interfaces
      !mes_extra
      CALL mesh%build_loc_to_glob

   CONTAINS

      SUBROUTINE build_outer_inner_edges(jj_HL, jjs_HL)
         USE my_util, ONLY: error_petsc, to_str, local_error_petsc
         IMPLICIT NONE
         INTEGER, DIMENSION(:,:), INTENT(IN) :: jj_HL
         INTEGER, DIMENSION(:,:), INTENT(IN) :: jjs_HL
         INTEGER :: out_new_node_min, out_new_node_max, i, num_edge, old_nc
         INTEGER, DIMENSION(2)    :: cur_edge
         INTEGER, DIMENSION(3, 2) :: tab_edges

         ALLOCATE(outer_edges_coarse    (3, refine_factor+1                   ))
         ALLOCATE(outer_edges           (3, refine_factor,                   2))
         ALLOCATE(inner_edges           (3*(refine_factor-1)*refine_factor/2, 2))
         ALLOCATE(inner_edges_to_subcell(3*(refine_factor-1)*refine_factor/2, 2))


!========================= EXAMPLE OF THE ARRAYS THAT ARE CONSTRUCTED CONSIDERING THE FOLLOWING JJ_HL AND JJS_HL
         ! jj_HL(1, :) = [4, 5, 6]
         ! jj_HL(2, :) = [1, 5, 6]
         ! jj_HL(3, :) = [2, 4, 6]
         ! jj_HL(4, :) = [3, 4, 5]

         ! jjs_HL(1, :) = [1, 3]
         ! jjs_HL(2, :) = [2, 3]

!       outer_edges(1, 1, :) = [2, 4]
!       outer_edges(1, 2, :) = [3, 4]

!       outer_edges(2, 1, :) = [1, 5]
!       outer_edges(2, 2, :) = [3, 5]

!       outer_edges(3, 1, :) = [1, 6]
!       outer_edges(3, 2, :) = [2, 6]

!       offset = 0
!       outer_edges_coarse(1, 1:2) = [2,3]
!       outer_edges_coarse(1, 3:) = [(3+offset+n, n=1, refine_factor - 1)]

!       offset = refine_factor - 1
!       outer_edges_coarse(2, 1:2) = [1,3]
!       outer_edges_coarse(2, 3:) = [(3+offset+n, n=1, refine_factor - 1)]

!       offset = 2 * (refine_factor - 1)
!       outer_edges_coarse(3, 1:2) = [1,2]
!       outer_edges_coarse(3, 3:) = [(3+offset+n, n=1, refine_factor - 1)]

!       inner_edges(1,:) = [4, 5]
!       inner_edges(2,:) = [5, 6]
!       inner_edges(3,:) = [4, 6]

!       inner_edges_to_subcell(1,:) = [1, 4]
!       inner_edges_to_subcell(2,:) = [1, 2]
!       inner_edges_to_subcell(3,:) = [1, 3]

!========================= EXAMPLE OF THE ARRAYS THAT ARE CONSTRUCTED CONSIDERING THE FOLLOWING JJ_HL AND JJS_HL

!=== Build arrays containing new fine outer edges
         out_new_node_min = 3 + 1
         out_new_node_max = 3 + 1 + (refine_factor - 1) - 1
         outer_edges(1, 1, :) = [2, out_new_node_min]
         DO i=1, refine_factor-2
            outer_edges(1, i+1, :) = [out_new_node_min+i-1, out_new_node_min+i]
         END DO
         outer_edges(1, refine_factor, :) = [3, out_new_node_max]

         out_new_node_min = 3 + 1 + (refine_factor - 1)
         out_new_node_max = 3 + 1 + 2*(refine_factor - 1) - 1
         outer_edges(2, 1, :) = [1, out_new_node_min]
         DO i=1, refine_factor-2
            outer_edges(2, i+1, :) = [out_new_node_min+i-1, out_new_node_min+i]
         END DO
         outer_edges(2, refine_factor, :) = [3, out_new_node_max]

         out_new_node_min = 3 + 1 + 2*(refine_factor - 1)
         out_new_node_max = 3 + 1 + 3*(refine_factor - 1) - 1
         outer_edges(3, 1, :) = [1, out_new_node_min]
         DO i=1, refine_factor-2
            outer_edges(3, i+1, :) = [out_new_node_min+i-1, out_new_node_min+i]
         END DO
         outer_edges(3, refine_factor, :) = [2, out_new_node_max]


!=== Build arrays containing old coarse (outer) edges 
!=== NAME IS MISLEADING AS IT CONTAINS ALL POINTS ON THAT OUTER EDGE
         offset = 0
         outer_edges_coarse(1, 1:2) = [2,3]
         outer_edges_coarse(1, 3:) = [(3+offset+n, n=1, refine_factor - 1)]

         offset = refine_factor - 1
         outer_edges_coarse(2, 1:2) = [1,3]
         outer_edges_coarse(2, 3:) = [(3+offset+n, n=1, refine_factor - 1)]

         offset = 2 * (refine_factor - 1)
         outer_edges_coarse(3, 1:2) = [1,2]
         outer_edges_coarse(3, 3:) = [(3+offset+n, n=1, refine_factor - 1)]

!=== Build arrays containing new inner edges and the two subcells containing each of them
         num_edge = 0
         inner_edges_to_subcell = -1
         tab_edges(1,:) = [1, 2]
         tab_edges(2,:) = [2, 3]
         tab_edges(3,:) = [1, 3]
         DO nc=1, SIZE(jj_HL, 1) !loop on new subcell
edge_loop:  DO n=1, ne !loop on edges of each subcell
               cur_edge = jj_HL(nc, tab_edges(n, :))
               ! First making sure that the considered edge is indeed an inner edge
               DO i=1, ne
                  IF (ANY(outer_edges_coarse(i, :)==cur_edge(1))) THEN
                     IF (ANY(outer_edges_coarse(i, :)==cur_edge(2))) CYCLE edge_loop
                  END IF
               END DO
               i = 1 !to handle case num_edge = 0
               DO i=1, num_edge
                  IF (MAXVAL(ABS(cur_edge-inner_edges(i, :)))==0) THEN 
                     !cur_edge was already explored, but here we found the second cell containing it, so we update inner_edges_to_subcell
                     IF (MINVAL(inner_edges_to_subcell(i,:)) /= -1) THEN
                        CALL error_petsc("BUG in refine_iso_grid/build_edges => how can an edge belong to three different cells???")
                     END IF
                     old_nc = inner_edges_to_subcell(i,1)
                     inner_edges_to_subcell(i, 1) = MIN(nc, old_nc)
                     inner_edges_to_subcell(i, 2) = MAX(nc, old_nc)
                     EXIT
                  END IF
               END DO
               IF (i==num_edge+1) THEN !means cur_edge was never explored in previous loops
                  num_edge = num_edge + 1
                  inner_edges(num_edge, :) = cur_edge
                  inner_edges_to_subcell(num_edge, 1) = nc
               END IF
            END DO edge_loop
         END DO

         IF (num_edge /= 3*(refine_factor-1)*refine_factor/2) THEN
            CALL error_petsc(&
            &"BUG in refine_iso/build_edges => did not find all inner edges: found a total of "//&
            &to_str(num_edge)//" while there should be "//to_str(3*(refine_factor-1)*refine_factor/2)//". List of edges &
            &that were found pts1 = "//to_str(inner_edges(:,1))//", pts2 = "//to_str(inner_edges(:,2)))
         END IF

         IF (MINVAL(inner_edges_to_subcell) == -1) THEN
            CALL error_petsc(&
            &"BUG in refine_iso/build_edges => all inner edges were not correctly attributed their two subcells.&
            & Subcell 1 = "//to_str(inner_edges_to_subcell(:,1))//", Subcell 2 = "//to_str(inner_edges_to_subcell(:,2)))
         END IF
         
      !=== dummy to avoid warning unused variable
         RETURN
         i = MAXVAL(JJs_HL)
      !=== dummy to avoid warning unused variable
      END SUBROUTINE build_outer_inner_edges

      SUBROUTINE get_edge(n, nc, n_ks, n_loc_1, n_loc_2, idx1, idx2, is_inner_edge)
         USE my_util, ONLY: error_petsc, to_str
         IMPLICIT NONE
         INTEGER,               INTENT(IN)  :: n, nc
         INTEGER, DIMENSION(2), INTENT(OUT) :: n_ks
         INTEGER, INTENT(OUT) :: n_loc_1, n_loc_2
         INTEGER, INTENT(OUT) :: idx1, idx2
         LOGICAL, INTENT(OUT) :: is_inner_edge

         n_ks(1) = MODULO(n, ne) + 1
         n_ks(2) = MODULO(n + 1, ne) + 1
         IF (n_ks(1)>n_ks(2)) THEN
            n_ks = (/n_ks(2), n_ks(1)/)
         END IF
         ! Determine if edge [n_ks(1), nk_s(2)] is the refinement of an old coarse edge or a brand new inner edge
         n_loc_1 = jj_HL(nc, n_ks(1))
         n_loc_2 = jj_HL(nc, n_ks(2))
         IF (n_loc_1 > n_loc_2) THEN
            CALL error_petsc("BUG in refinement_iso_grid => your jj_HL might not have &
            &increasing numbering of local nodes: "//to_str(jj_HL(nc,:)))
         END IF

         CALL find_edge(n_loc_1, n_loc_2, idx1, idx2, is_inner_edge)
      END SUBROUTINE get_edge

      SUBROUTINE find_edge(n_loc_1, n_loc_2, idx1, idx2, is_inner_edge)
         !> subroutine based on correct construction of outer_edges(:,:,:)/inner_edges(:,:)/inner_edges_to_subcell
         !! Two cases:
         !!
         !! is_inner_edge = .TRUE. ==> the edge is inner_edge(idx1, :)
         !!                            idx2 is dummy
         !! is_inner_edge = .FALSE. ==> idx1 is the index of the old outer coarse edge, idx2 the local number of the new edge
         !!                             the edge is outer_edge(idx2, idx1)
         USE my_util, ONLY: error_petsc, to_str
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: n_loc_1, n_loc_2
         INTEGER, INTENT(OUT) :: idx1, idx2
         LOGICAL, INTENT(OUT) :: is_inner_edge
         INTEGER :: n1, n2

         IF (n_loc_1 >= n_loc_2) THEN
            CALL error_petsc(&
            &"BUG in refinement_iso_grid/find_edge => should have n_loc_1 < n_loc_2, values found are: "&
            &//to_str(n_loc_1)//', '//to_str(n_loc_2))
         END IF

         !Check inside inner_edges first
         DO n1=1, SIZE(inner_edges,1)
            IF (n_loc_1==inner_edges(n1,1) .AND. n_loc_2==inner_edges(n1,2)) THEN !this assumes n_loc_1<n_loc_2 and inner_edges(n1,1)<inner_edges(n1,2)
               is_inner_edge = .TRUE.
               idx1 = n1
               idx2 = -1
               RETURN
            END IF
         END DO

         !Finally check outer edges
         is_inner_edge = .FALSE.
         DO n2=1, SIZE(outer_edges,1)
            DO n1=1, SIZE(outer_edges,2)
               IF (n_loc_1==outer_edges(n2,n1,1) .AND. n_loc_2==outer_edges(n2,n1,2)) THEN !this assumes n_loc_1<n_loc_2 and inner_edges(n1,1)<inner_edges(n1,2)
                  idx1 = n1
                  idx2 = n2
                  RETURN
               END IF
            END DO
         END DO

         CALL error_petsc(&
         &"BUG in refinement_iso_grid/find_edge => &
         &did not find an edge for (n_loc_1,n_loc_2)=("&
         &//to_str(n_loc_1)//','//to_str(n_loc_2)//')')

      END SUBROUTINE find_edge

   END SUBROUTINE general_refinement_iso_grid_distributed


   ! ! TODO add back curved boundaries
   ! !   SUBROUTINE is_on_curved_interface(side, iso, interface)
   ! !      USE input_data
   ! !      INTEGER :: side, interface
   ! !      LOGICAL :: iso
   ! !      interface = -1
   ! !      iso = .FALSE.
   ! !
   ! !      IF (inputs%nb_spherical + inputs%nb_curved > 0) THEN
   ! !         IF (MINVAL(ABS(side - inputs%list_spherical)) == 0 .OR. &
   ! !              MINVAL(ABS(side - inputs%list_curved)) == 0) THEN
   ! !            DO interface = 1, inputs%nb_spherical + inputs%nb_curved
   ! !               IF (interface <= inputs%nb_spherical) THEN
   ! !                  IF (side - inputs%list_spherical(interface) == 0) EXIT
   ! !               ELSE
   ! !                  IF (side - inputs%list_curved(interface - inputs%nb_spherical) == 0) EXIT
   ! !               END IF
   ! !            END DO
   ! !            iso = .TRUE.
   ! !         ELSE
   ! !            iso = .FALSE.
   ! !         END IF
   ! !      ELSE
   ! !         iso = .FALSE.
   ! !      END IF
   ! !
   ! !   END SUBROUTINE is_on_curved_interface
   ! !
   ! !   SUBROUTINE rescale_to_curved_boundary(rr, interface)
   ! !      USE input_data
   ! !      USE boundary
   ! !      REAL(KIND = 8), DIMENSION(2) :: rr, rr_ref
   ! !      INTEGER :: interface
   ! !      REAL(KIND = 8) :: rescale, pi = ACOS(-1.d0), theta
   ! !      IF (interface <= inputs%nb_spherical) THEN
   ! !         rr_ref = rr - inputs%origin_spherical(:, interface)
   ! !         rescale = inputs%radius_spherical(interface) / SQRT(SUM(rr_ref * rr_ref))
   ! !         rr = rr_ref * rescale + inputs%origin_spherical(:, interface)
   ! !      ELSE
   ! !         rr_ref = rr - inputs%origin_curved(:, interface - inputs%nb_spherical)
   ! !         theta = pi - pi / 2 * (1 + sgn(rr_ref(1))) * (1 - sgn(rr_ref(2) * rr_ref(2))) &
   ! !              - pi / 4 * (2 + sgn(rr_ref(1))) * sgn(rr_ref(2)) &
   ! !              - sgn(rr_ref(1) * rr_ref(2)) * ATAN((ABS(rr_ref(1)) - ABS(rr_ref(2))) / (ABS(rr_ref(1)) + ABS(rr_ref(2))))
   ! !         rescale = curved_boundary_radius(inputs%list_curved(interface - inputs%nb_spherical), theta) &
   ! !              / SQRT(SUM(rr_ref * rr_ref))
   ! !         rr = rr_ref * rescale + inputs%origin_curved(:, interface - inputs%nb_spherical)
   ! !      END IF
   ! !
   ! !   END SUBROUTINE rescale_to_curved_boundary

   FUNCTION sgn(x) RESULT(out)
      REAL(KIND = 8) :: x, out
      IF (x > 0.d0) THEN
         out = 1.d0
      ELSE IF (x < 0.d0) THEN
         out = -1.d0
      ELSE
         out = 0.d0
      ENDIF
   END FUNCTION sgn

END MODULE refine_mesh
