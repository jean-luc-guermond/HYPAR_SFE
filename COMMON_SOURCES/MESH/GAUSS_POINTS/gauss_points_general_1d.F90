!===
!Author: Jean-Luc Guermond, Copyright December 24th 2018
!===
MODULE gauss_points_1d
   USE space_dim, ONLY: k_dim
   PRIVATE
   PUBLIC create_gauss_points_1d
   INTEGER, PARAMETER :: k_d = k_dim
   ABSTRACT INTERFACE
      SUBROUTINE template_element_1d(w, d, p, n_ws, l_Gs)
         IMPLICIT NONE
         INTEGER,                                INTENT(IN)  :: n_ws, l_Gs
         REAL(KIND=8), DIMENSION(   n_ws, l_Gs), INTENT(OUT) :: w
         REAL(KIND=8), DIMENSION(1, n_ws, l_Gs), INTENT(OUT) :: d
         REAL(KIND=8), DIMENSION(l_Gs),          INTENT(OUT) :: p
      END SUBROUTINE template_element_1d

      SUBROUTINE template_element_2d(w, d, p, n_w, l_G)
         IMPLICIT NONE
         INTEGER,                              INTENT(IN)  :: n_w, l_G
         REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
         REAL(KIND=8), DIMENSION(2, n_w, l_G), INTENT(OUT) :: d
         REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p
      END SUBROUTINE template_element_2d

      SUBROUTINE template_element_2d_boundary(face, d, w, n_w, l_Gs)
         IMPLICIT NONE
         INTEGER,                                INTENT(IN)  :: n_w, l_Gs
         INTEGER,                                INTENT(IN)  :: face
         REAL(KIND=8), DIMENSION(2, n_w,  l_Gs), INTENT(OUT) :: d
         REAL(KIND=8), DIMENSION(   n_w, l_Gs),  INTENT(OUT) :: w
      END SUBROUTINE template_element_2d_boundary

      SUBROUTINE template_element_1d_at_nodes(d, n_ws)
         IMPLICIT NONE
         INTEGER,                             INTENT(IN)  :: n_ws
         REAL(KIND=8), DIMENSION(n_ws, n_ws), INTENT(OUT) :: d
      END SUBROUTINE template_element_1d_at_nodes
   END INTERFACE

   PROCEDURE(template_element_1d),          POINTER :: element_1d_Pk
   PROCEDURE(template_element_2d),          POINTER :: element_2d_Pk
   PROCEDURE(template_element_2d_boundary), POINTER :: element_2d_Pk_boundary
   PROCEDURE(template_element_1d_at_nodes), POINTER :: element_1d_Pk_at_nodes

CONTAINS
   SUBROUTINE create_gauss_points_1d(mesh, type_fe)
      USE def_type_mesh
      USE GP_2d_p1
      USE GP_2d_p2
      USE GP_2d_p3
      USE sub_plot
      IMPLICIT NONE
      TYPE(mesh_type), TARGET :: mesh
      INTEGER :: type_fe
      INTEGER, POINTER :: me, mes
      INTEGER, DIMENSION(:, :), POINTER :: jj
      INTEGER, DIMENSION(:, :), POINTER :: js
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: rr
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: ww
      REAL(KIND = 8), DIMENSION(:, :, :, :), POINTER :: dw
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: rj
      REAL(KIND = 8), DIMENSION(:, :, :), ALLOCATABLE :: dd
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: pp
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: r
      REAL(KIND = 8), DIMENSION(k_d, k_d) :: dr
      REAL(KIND = 8) :: rjac, rjacs, x
      INTEGER :: m, l, k, k1, n, n1, n2, ms, ns, ls, face, cote, orient
      INTEGER :: n_w, n_ws, l_G, l_Gs
      REAL(KIND = 8), DIMENSION(k_d) :: rnor, rsd

      SELECT CASE(type_fe)
      CASE(1)
         n_w = 2;  n_ws = 1; l_G = 2; l_Gs = 0
         element_1d_Pk => element_1d_P1
      CASE(2)
         n_w = 3;  n_ws = 1; l_G = 3; l_Gs = 0
         element_1d_Pk => element_1d_P2
      CASE(3)
         n_w = 4; n_ws = 1; l_G = 4; l_Gs = 0
         element_1d_Pk => element_1d_P3
      CASE DEFAULT
         WRITE(*, *) ' FE not programmed yet', type_fe
         STOP
      END SELECT
      ALLOCATE(dd(k_d, n_w, l_G), pp(l_G), r(n_w))

      me => mesh%me
      jj => mesh%jj
      rr => mesh%rr

      ALLOCATE(mesh%gauss%ww(n_w, l_g))
      ALLOCATE(mesh%gauss%dw(k_d, n_w, l_G, me))
      ALLOCATE(mesh%gauss%rj(l_G, me))

      ww => mesh%gauss%ww
      dw => mesh%gauss%dw
      rj => mesh%gauss%rj

      mesh%gauss%k_d = k_d
      mesh%gauss%n_w = n_w
      mesh%gauss%n_e = 1
      mesh%gauss%l_G = l_G
      mesh%gauss%n_ws = n_ws
      mesh%gauss%l_Gs = l_Gs

      
      CALL element_1d_Pk(ww, dd, pp, n_w, l_G)
      !===create jacobian elements
      DO m = 1, me

         DO l = 1, l_G
            DO k = 1, k_d
               r = rr(k, jj(:, m))
               DO k1 = 1, k_d
                  dr(k, k1) = SUM(r * dd(k1, :, l))
               ENDDO
            ENDDO
            rjac = dr(1, 1)
            rj(l, m) = ABS(dr(1, 1)) * pp(l)
            DO n = 1, n_w
               dw(1, n, l, m) = dd(1, n, l) / rjac
            ENDDO
         ENDDO
      ENDDO

      DEALLOCATE(dd, pp, r)

   END SUBROUTINE create_gauss_points_1d
END MODULE gauss_points_1d
