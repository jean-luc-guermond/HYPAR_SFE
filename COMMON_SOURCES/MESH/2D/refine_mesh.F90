MODULE refine_mesh
   USE mesh_tools
   PUBLIC :: create_iso_grid_distributed, refinement_iso_grid_distributed
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
      USE my_util, ONLY : error_petsc, to_str
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_p1, mesh
      INTEGER, INTENT(IN) :: type_fe
      LOGICAL, DIMENSION(:), ALLOCATABLE :: virgin
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: j_edge, jjs_edge
      INTEGER :: np, me, mes, mes_int, nw, nws, kd, n, m, k, l, n_dof, dom_np, ne
      INTEGER :: n1, n2, ms, n_start, n_end, n_loc
      INTEGER :: n_k1, n_k2, m_op_k, kk, i, mm, p_e, p_c
      INTEGER, DIMENSION(2) :: n_ks
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: r_mid
      INTEGER :: nb_angle, f_dof, f_dof_edge, f_dof_cell,edge_g, edge_l, n_new_start, proc, nb_proc, edges, p, cell_g, cell_l
      INTEGER :: interface
      LOGICAL :: iso

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
      f_dof = type_fe - 1
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

      ALLOCATE(mesh%jj(mesh%gauss%n_w, me))   
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
      ALLOCATE(virgin(mesh_p1%medge), j_edge(nw * f_dof, me), jjs_edge(f_dof, mes), r_mid(kd))

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
         ! n_loc = f_dof_edge * ne + 3 !offset by summits and edge points
         ! DO n = 1, f_dof_cell !TODO AFTER VALIDATING REFINEMENT ORDER 3
         IF (type_fe==3) THEN
            n_dof = n_dof + 1
            n_new_start = m + mesh_p1%dom_np + mesh_p1%medge * f_dof_edge
            mesh%jj(10, m) = n_new_start
            ! mesh%jj(n_loc+n, m) = n_new_start
            mesh%rr(:, n_new_start) = &
                 (mesh_p1%rr(:, mesh_p1%jj(1, m)) + mesh_p1%rr(:, mesh_p1%jj(2, m)) + mesh_p1%rr(:, mesh_p1%jj(3, m))) / type_fe
         END IF
         ! END DO
      END DO

      IF (type_fe == 3 .AND. n_dof /= mesh_p1%me) THEN
         CALL error_petsc('BUG in create_iso_grid_distributed, type_fe == 3 .and. n_dof /= mesh_p1%me')
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
         DO l = 1, f_dof
            jjs_edge(l, ms) = j_edge((kk - 1) * f_dof + l, m) !===New index created
         END DO
      ENDDO
      mesh%jjs(1:nws, :) = mesh_p1%jjs
      DO i = 1, SIZE(mesh%jjs, 2)
         DO n = 1, nws
            IF (mesh%jjs(n, i) > mesh_p1%dom_np) THEN
               mesh%jjs(n, i) = mesh_p1%jjs(n, i) + mesh_p1%medge * f_dof + mesh_p1%me * (f_dof - 1)
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

            DO l = 1, f_dof
               mesh%jj_extra(nw + (k - 1) * f_dof + l, m) = l &
                    + (edge_l - 1) * f_dof + mesh_p1%domnp(p_e) + mesh%disp(p_e) - 1
            END DO

            IF (type_fe==3) THEN
               mesh%jj_extra(10, m) = cell_l + mesh_p1%domedge(p_c) * 2 + mesh_p1%domnp(p_c) + mesh%disp(p_c) - 1
            END IF
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
            CALL error_petsc("BUG in create_iso_grid_distributed: did not find extra cell")
         END IF
         DO l = 1, f_dof
            mesh%jjs_extra(nws + l, ms) = mesh%jj_extra(nw + (kk - 1) * f_dof + l, m)
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

            DO l = 1, f_dof
               mesh%rrs_extra(:, nw + (k - 1) * f_dof + l, ms) = mesh_p1%rrs_extra(:, n_start, ms) &
                    + l * (mesh_p1%rrs_extra(:, n_end, ms) - mesh_p1%rrs_extra(:, n_start, ms)) / type_fe
               !               IF (iso) THEN
               !                  CALL rescale_to_curved_boundary(mesh%rrs_extra(:, nw + (k - 1) * f_dof + l, ms), interface)
               !               END IF
            END DO
         END DO

         IF (type_fe==3) THEN
            mesh%rrs_extra(:, 10, ms) = &
                 (mesh_p1%rrs_extra(:, 1, ms) + mesh_p1%rrs_extra(:, 2, ms) + mesh_p1%rrs_extra(:, 3, ms)) / 3
         END IF
      ENDDO

      CALL mesh%build_loc_to_glob

   END SUBROUTINE  create_iso_grid_distributed

   SUBROUTINE refinement_iso_grid_distributed(mesh_p1)
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_p1, mesh_p2
      INTEGER :: refine_order

      refine_order = 2

      CALL create_iso_grid_distributed(mesh_p1, mesh_p2, refine_order)
      CALL free_mesh(mesh_p1)
      CALL refinement_iso_grid_distributed_order_2(mesh_p2, mesh_p1)
      CALL free_mesh(mesh_p2)
   END SUBROUTINE refinement_iso_grid_distributed

   SUBROUTINE refinement_iso_grid_distributed_order_2(mesh_p2, mesh)
      !===jj(:, :)    nodes of the  volume_elements of the input grid
      !===jjs(:, :)    nodes of the surface_elements of the input grid
      !===rr(:, :)    cartesian coordinates of the nodes of the input grid
      !===m_op(:,:)   volume element opposite to each node
      !===neigh_el(:) volume element ajacent to the surface element
      !===jj_f(:, :)  nodes of the  volume_elements of the output p2 grid
      !===jjs_f(:, :)  nodes of the surface_elements of the output p2 grid
      !===rr_f(:, :)  cartesian coordinates of the nodes of the output p2 grid
      USE def_type_mesh
      USE my_util, ONLY: error_petsc
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_p2, mesh
      INTEGER :: np, me, mes, nw, nws, kd, n, m, k, n_dof, dom_np, e_g, mextra
      INTEGER :: n1, n2, ms, neigh, k_neigh, n_kneigh1, n_kneigh2, swap, n_loc
      INTEGER :: i, p_c, m_new, e_k, p_j, ne, refine_order
      INTEGER, DIMENSION(3) :: edges_g, edges_l, p_es, sub_cell
      INTEGER, DIMENSION(2) :: n_ks
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: jj_HL, jjs_HL
      INTEGER :: proc, nb_proc, p, cell_g, cell_l
      INTEGER :: m1, m2, interface, m_center, tab1, tab2, mes_int
      LOGICAL :: iso


      refine_order = 2
      ne = 3

      ALLOCATE(jj_HL(refine_order**2, ne))
      ALLOCATE(jjs_HL(refine_order, ne-1))

      jj_HL(1, :) = [4, 5, 6]
      jj_HL(2, :) = [1, 5, 6]
      jj_HL(3, :) = [2, 4, 6]
      jj_HL(4, :) = [3, 4, 5]

      jjs_HL(1, :) = [1, 3]
      jjs_HL(2, :) = [2, 3]

      IF (mesh_p2%me == 0) RETURN

      CALL mesh%info%copy(mesh_p2%info)
      nw  = 3   !===nodes in each volume element (3 in 2D)
      nws = 2  !===nodes in each edge element (2 in 2D)
      me = SIZE(mesh_p2%jj, 2)   !===number of cells
      kd = SIZE(mesh_p2%rr, 1)   !===space dimensions
      mes     = SIZE(mesh_p2%jjs, 2)
      mes_int = SIZE(mesh_p2%jjs_int, 2)
      nb_proc = SIZE(mesh_p2%domnp)

      mesh%me      = refine_order**2 * mesh_p2%me
      mesh%mes     = refine_order * mesh_p2%mes !---> something with news to take into account
      mesh%mes_int = refine_order * mesh_p2%mes_int
      mesh%np      = mesh_p2%np
      mesh%medge   = refine_order * mesh_p2%medge + 3*(refine_order-1)*refine_order/2*mesh_p2%me
      mesh%medges  = refine_order * mesh_p2%medges

      ALLOCATE(mesh%jj(nw, mesh%me)) !--->done
      ALLOCATE(mesh%jjs(nws, mesh%mes))  !--->done
      ALLOCATE(mesh%rr(kd, mesh%np)) !--->done
      ALLOCATE(mesh%loc_to_glob(mesh%np)) !--->done

      ALLOCATE(mesh%jce(nw, mesh%me))  !--->done
      !ALLOCATE(mesh%jev(nw - 1, mesh%medge)) !----->never needed so not constructed
      ALLOCATE(mesh%jees(mesh%medges))  !--->done
      ALLOCATE(mesh%jecs(mesh%medges))  !--->done

      ALLOCATE(mesh%neigh(nw, mesh%me)) !--->done
      ALLOCATE(mesh%sides(mesh%mes))   !---> still don't know what this is ?? done ?
      ALLOCATE(mesh%neighs(mesh%mes))  !--->done
      ALLOCATE(mesh%i_d(mesh%me)) !--->done

      ALLOCATE(mesh%jjs_int(nws, mesh%mes_int))  !--->done
      ALLOCATE(mesh%neighs_int(2, mesh%mes_int)) !--->done
      ALLOCATE(mesh%sides_int(mesh%mes_int))

      mesh%dom_np = mesh_p2%dom_np
      CALL mesh%create_comm(mesh_p2%comm)
      CALL mesh%gather_dom_np
      CALL mesh%gather_me
      CALL mesh%gather_medge

      mesh%nis = mesh_p2%nis
      ALLOCATE(mesh%isolated_interfaces(mesh_p2%nis, 2))
      mesh%isolated_interfaces = mesh_p2%isolated_interfaces
      ALLOCATE(mesh%isolated_jjs(mesh_p2%nis))

      IF (kd == 3) THEN
         CALL error_petsc('refinement_iso_grid_distributed: 3D case not programmed yet !')
      END IF

      

      proc = mesh_p2%get_proc(mesh_p2%loc_to_glob(1), 'np')

      !===GENERATION OF THE Pk GRID
      mesh%rr = mesh_p2%rr

      ALLOCATE(mesh%proc_np_loc(2, mesh%np-mesh%dom_np))
      mesh%proc_np_loc = mesh_p2%proc_np_loc
      mesh%isolated_jjs = mesh_p2%isolated_jjs

      DO m=1, me
         m_center = refine_order**2 * (m - 1) + 1
         DO i=1, refine_order**2
            mesh%jj(:, m_center + i - 1) = mesh_p2%jj(jj_HL(i, :), m)
         END DO
      END DO


      n_dof = 0
      DO m = 1, me !===loop on the elements unrefined mesh
         edges_g = mesh_p2%jce(:, m)
         edges_l = edges_g - mesh_p2%disedge(proc) + 1
         m_center = refine_order**2 * (m - 1) + 1
         !===Center cell
         mesh%i_d(m_center) = mesh_p2%i_d(m)
         DO i = 1, nw
            !===Setting up neighbours and edges
            mesh%neigh(i, m_center) = m_center + i
            mesh%i_d(mesh%neigh(i, m_center)) = mesh_p2%i_d(m)
            mesh%jce(i, m_center) = mesh%disedge(proc) - 1 + refine_order * mesh_p2%medge &
                                    + 3*(refine_order-1)*refine_order/2*(m - 1) + i
         END DO


         !===Corner cells
         DO k = 1, nw
            m_new = mesh%neigh(k, m_center)

            n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
            IF (n_ks(1)>n_ks(2)) THEN
               n_ks = (/n_ks(2), n_ks(1)/)
            END IF

            !===Setting the center cell as a neighbour
            mesh%neigh(1, m_new) = m_center
            mesh%jce(1, m_new) = mesh%jce(k, m_center)

            DO e_k = 1, 2

               !===Adding neighbours
               neigh = mesh_p2%neigh(n_ks(MODULO(e_k, 2) + 1), m)

               IF (neigh <= 0) THEN
                  mesh%neigh(e_k + 1, m_new) = neigh
               ELSE
                  DO k_neigh = 1, nw
                     IF (mesh_p2%neigh(k_neigh, neigh) == m) exit
                  END DO
                  n_kneigh1 = MODULO(k_neigh, nw) + 1
                  n_kneigh2 = MODULO(k_neigh + 1, nw) + 1
                  IF (n_kneigh1>n_kneigh2) THEN
                     swap = n_kneigh1
                     n_kneigh1 = n_kneigh2
                     n_kneigh2 = swap
                  END IF
                  IF (k < n_ks(e_k)) THEN
                     mesh%neigh(e_k + 1, m_new) = 4 * (neigh - 1) + 1 + n_kneigh1
                  ELSE
                     mesh%neigh(e_k + 1, m_new) = 4 * (neigh - 1) + 1 + n_kneigh2
                  END IF
               END IF

               !===Adding edge
               e_g = mesh_p2%jce(n_ks(MODULO(e_k, 2) + 1), m)
               IF (e_g < mesh_p2%disedge(proc)) THEN
                  p = mesh_p2%get_proc(e_g, 'medge')
               ELSE
                  p = proc
               END IF

               IF (k < n_ks(e_k)) THEN
                  mesh%jce(e_k + 1, m_new) = mesh%disedge(p) - 1 + e_g - mesh_p2%disedge(p) + 1
               ELSE
                  mesh%jce(e_k + 1, m_new) = mesh%disedge(p) - 1 + mesh_p2%domedge(p) + e_g - mesh_p2%disedge(p) + 1
               END IF

               IF (e_g < mesh_p2%disedge(proc)) THEN
                  DO ms = 1, mesh_p2%medges
                     IF (mesh_p2%jees(ms) == e_g) EXIT
                  END DO
                  IF (k < n_ks(e_k)) THEN
                     mesh%jees(ms) = mesh%jce(e_k + 1, m_new)
                     mesh%jecs(ms) = m_new
                  ELSE
                     mesh%jees(mesh_p2%medges + ms) = mesh%jce(e_k + 1, m_new)
                     mesh%jecs(mesh_p2%medges + ms) = m_new
                  END IF
               END IF
            END DO
         END DO
      END DO

      DO ms = 1, mes
         DO i=1, refine_order
            mesh%jjs(:, ms + (i - 1) * mes) = mesh_p2%jjs(jjs_HL(i, :), ms)
         ENDDO
      END DO

      !===Surface elements
      DO ms = 1, mes
         m = mesh_p2%neighs(ms)
         !===Finding the corresponding side in the cell
         DO k = 1, nw
            IF (MINVAL(ABS(mesh_p2%jj(k, m) - mesh_p2%jjs(:, ms)))/=0) EXIT
         ENDDO

         n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
         IF (n_ks(1)>n_ks(2)) THEN
            n_ks = (/n_ks(2), n_ks(1)/)
         END IF

         mesh%neighs(ms) = mesh%neigh(n_ks(1), 4 * (m - 1) + 1)
         mesh%neighs(mes + ms) = mesh%neigh(n_ks(2), 4 * (m - 1) + 1)
         mesh%sides(ms) = mesh_p2%sides(ms)
         mesh%sides(mes + ms) = mesh_p2%sides(ms)

         !         CALL is_on_curved_interface(mesh_p2%sides(ms), iso, interface)
         !         IF (iso) THEN
         !            CALL rescale_to_curved_boundary(mesh%rr(:, mesh%jj(k, 4 * (m - 1) + 1)), interface)
         !         END IF
      ENDDO

      !===Internal surface elements
      DO ms = 1, mes_int
         m = mesh_p2%neighs_int(1, ms)
         !===Finding the corresponding side in the cell
         DO k = 1, nw
            IF (MINVAL(ABS(mesh_p2%jj(k, m) - mesh_p2%jjs_int(:, ms)))/=0) EXIT
         ENDDO

         n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
         IF (n_ks(1)>n_ks(2)) THEN
            n_ks = (/n_ks(2), n_ks(1)/)
         END IF

         ! VB 27/05/2026 => no test or guarantee that it works, i also don't know where these are used
         DO i=1, refine_order
            mesh%jjs_int(:, ms + (i - 1) * mes) = mesh_p2%jjs_int(jjs_HL(i, :), ms)
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
         mesh%sides_int(ms) = mesh_p2%sides_int(ms)
         mesh%sides_int(mes_int + ms) = mesh_p2%sides_int(ms)
         ! VB 27/05/2026 => see also create_iso_grid_distributed, no test or guarantee that it works, i also don't know where these are used

         !         CALL is_on_curved_interface(mesh_p2%sides_int(ms), iso, interface)
         !         IF (iso) THEN
         !            CALL rescale_to_curved_boundary(mesh%rr(:, mesh%jj(k, 4 * (m - 1) + 1)), interface)
         !         END IF
      ENDDO

      !===Counting number of new extra cells
      mesh%mextra = 0
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
      !===In the mean time
      mesh%mextra = mesh_p2%mextra * refine_order**2
      ALLOCATE(mesh%jj_extra(nw, mesh%mextra)) !---->
      ALLOCATE(mesh%jce_extra(nw, mesh%mextra)) !---->
      ALLOCATE(mesh%jcc_extra(mesh%mextra)) !---->

      mextra = 0
      DO m=1, mesh_p2%mextra
         DO i=1, refine_order**2
            mextra = mextra + 1
            mesh%jj_extra(:, mextra) = mesh_p2%jj_extra(jj_HL(i, :), m)
         END DO
      END DO


      !===Constructing the extra cells
      mextra = 0
      DO m = 1, mesh_p2%mextra

         !center cell
         mextra = mextra + 1

         cell_g = mesh_p2%jcc_extra(m)
         p_c = mesh_p2%get_proc(cell_g, 'me')
         cell_l = cell_g - mesh_p2%discell(p_c) + 1
         mesh%jcc_extra(mextra) = 4 * (cell_l - 1) + 1 + mesh%discell(p_c) - 1

         DO i = 1, 3
            mesh%jce_extra(i, mextra) = mesh%disedge(p_c) - 1 + 2 * mesh_p2%domedge(p_c) + 3 * (cell_l - 1) + i
         END DO

         !corner cells
         DO k = 1, nw
            n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
            IF (n_ks(1)>n_ks(2)) THEN
               n_ks = (/n_ks(2), n_ks(1)/)
            END IF

            mextra = mextra + 1
            cell_g = mesh_p2%jcc_extra(m)
            p_c = mesh_p2%get_proc(cell_g, 'me')

            cell_l = cell_g - mesh_p2%discell(p_c) + 1
            mesh%jcc_extra(mextra) = mesh%discell(p_c) - 1 + 4 * (cell_l - 1) + 1 + k

            edges_g = mesh_p2%jce_extra(:, m)

            DO i = 1, nw
               p_es(i) = mesh_p2%get_proc(edges_g(i), 'medge')
            END DO

            !===Adding the edges
            mesh%jce_extra(1, mextra) = mesh%disedge(p_c) - 1 + 2 * mesh_p2%domedge(p_c) + 3 * (cell_l - 1) + k

            DO e_k = 1, 2
               IF (k < n_ks(e_k)) THEN
                  mesh%jce_extra(e_k + 1, mextra) = mesh%disedge(p_es(n_ks(MODULO(e_k, 2) + 1))) - 1 + &
                       edges_g(n_ks(MODULO(e_k, 2) + 1)) - mesh_p2%disedge(p_es(n_ks(MODULO(e_k, 2) + 1))) + 1
               ELSE
                  mesh%jce_extra(e_k + 1, mextra) = mesh%disedge(p_es(n_ks(MODULO(e_k, 2) + 1))) - 1 + &
                       mesh_p2%domedge(p_es(n_ks(MODULO(e_k, 2) + 1))) + edges_g(n_ks(MODULO(e_k, 2) + 1)) &
                       - mesh_p2%disedge(p_es(n_ks(MODULO(e_k, 2) + 1))) + 1
               END IF
            END DO
         END DO
      END DO

      !===Constructing the extra cells at interfaces
      mesh%mes_extra = 2 * mesh_p2%mes_extra
      ALLOCATE(mesh%jjs_extra(nws, mesh%mes_extra))
      ALLOCATE(mesh%rrs_extra(2, nw, mesh%mes_extra))
      ALLOCATE(mesh%sides_extra(mesh%mes_extra), mesh%neighs_extra(mesh%mes_extra))

      mextra = 0
      DO m = 1, mesh_p2%mes_extra
         !CALL is_on_curved_interface(mesh_p2%sides_extra(m), iso, interface)

         cell_g = mesh_p2%neighs_extra(m)
         DO m1 = 1, mesh_p2%mextra !find associated extra cell
            IF (mesh_p2%jcc_extra(m1) == cell_g) EXIT
         END DO

         p_c = mesh_p2%get_proc(cell_g, 'me')
         cell_l = cell_g - mesh_p2%discell(p_c) + 1

         DO n = 1, 3 !===find side in cell
            IF (.NOT.(ANY(mesh_p2%jj_extra(n, m1) == mesh_p2%jjs_extra(:, m)))) THEN
               EXIT
            END IF
            IF (n == 3) CALL error_petsc('BUG in refinement : didnt find face in extra cell for extra edge')
         ENDDO

         !==cell index of edge
         n_ks = (/MODULO(n, nw) + 1, MODULO(n + 1, nw) + 1/)
         IF (n_ks(1)>n_ks(2)) THEN
            n_ks = (/n_ks(2), n_ks(1)/)
         END IF

         DO k = 1, nws
            mextra = mextra + 1
            mesh%sides_extra(mextra) = mesh_p2%sides_extra(m)
            mesh%neighs_extra(mextra) = mesh%discell(p_c) - 1 + 4 * (cell_l - 1) + 1 + n_ks(k)

            !!$   VB comment
            !!$   !! example on n = 2
            !!$   3
            !!$   5  4
            !!$   1  6  2
            !!$
            !!$   !! n = 2 ==> rrs receives triangles [1, 5, 6] and [3, 4, 5] (opposite to point nw = n = 2)
            !!$   !! k = 1 (resp 2) ==> second point considered is nw = 1 (resp nw = 3) => rrs receives triangle ecompassing
            !!$      edge [1, 5] (resp [3, 5]) => triangle is therefore [1, 5, 6] (resp [3, 4, 5])
            !!$   VB comment

            IF ((n==1 .AND. k==1) .OR. (n==3 .AND. k==2)) THEN
               sub_cell = [2, 4, 6]
            ELSE IF ((n==1 .AND. k==2) .OR. (n==2 .AND. k==2)) THEN
               sub_cell = [3, 4, 5]
            ELSE IF ((n==2 .AND. k==1) .OR. (n==3 .AND. k==1)) THEN
               sub_cell = [1, 5, 6]
            END IF
            mesh%rrs_extra(:, :, mextra) = mesh_p2%rrs_extra(:,sub_cell, (mextra+1)/2)

            ! IF (n == 1) THEN
            !    tab1 = 2
            !    tab2 = 3
            ! ELSE IF (n == 3) THEN
            !    tab1 = 3
            !    tab2 = 2
            ! ELSE
            !    IF (k == 1) THEN
            !       tab1 = 2
            !       tab2 = 3
            !    ELSE
            !       tab1 = 3
            !       tab2 = 2
            !    END IF
            ! END IF
            ! mesh%rrs_extra(:, tab1, mextra) = (mesh_p1%rrs_extra(:, n_ks(1), m) + mesh_p1%rrs_extra(:, n_ks(2), m)) / 2
            ! !            IF (iso) THEN
            ! !               CALL rescale_to_curved_boundary(mesh%rrs_extra(:, tab1, mextra), interface)
            ! !            END IF
            ! mesh%rrs_extra(:, tab2, mextra) = (mesh_p1%rrs_extra(:, n_ks(k), m) + mesh_p1%rrs_extra(:, n, m)) / 2

         END DO
      END DO

      
      mextra = 0
      DO m=1, mesh_p2%mes_extra
         DO i=1, refine_order
            mextra = mextra + 1
            mesh%jjs_extra(:, mextra) = mesh_p2%jjs_extra(jjs_HL(i, :), (mextra+1)/2)
         END DO
      END DO


      !jjs_extra !(extra layer of cells not own by proc but with dofs own by proc)
      !rrs_extra  ! coordinates for cells at interfaces
      !sides_extra, neighs_extra !interfaces
      !mes_extra
      CALL mesh%build_loc_to_glob

   END SUBROUTINE refinement_iso_grid_distributed_order_2


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
