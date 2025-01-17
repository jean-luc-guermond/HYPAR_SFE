MODULE mesh_1d
   USE def_type_mesh
   USE input_data
   USE space_dim
   uSE mesh_tools
   PUBLIC :: load_mesh_1d
   PRIVATE
CONTAINS

   SUBROUTINE load_mesh_1d(mesh)
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      INTEGER, PARAMETER :: in_unit = 30
      INTEGER :: type_fe
      INTEGER :: i, n, m
      REAL(KIND = 8) :: x0, x1, dx
      OPEN(in_unit, FILE = TRIM(ADJUSTL(inputs%directory)) // &
           '/' // TRIM(ADJUSTL(inputs%file_name)), FORM = 'formatted')
      READ(in_unit, *) mesh%me, inputs%type_fe
      READ(in_unit, *) x0, x1
      type_fe = inputs%type_fe
      ALLOCATE(mesh%jj(type_fe + 1, mesh%me))
      ALLOCATE(mesh%neigh(2, mesh%me))
      DO m = 1, mesh%me
         DO n = 1, type_fe + 1
            mesh%jj(n, m) = type_fe * (m - 1) + n
         END DO
         mesh%neigh(1, m) = m + 1
         mesh%neigh(2, m) = m - 1
      END DO
      mesh%np = type_fe * mesh%me + 1

      mesh%neigh(2, 1) = 0
      mesh%neigh(1, mesh%me) = 0
      ALLOCATE(mesh%i_d(mesh%me))
      mesh%i_d = 1
   write(*,*) 'ok111'

      ALLOCATE(mesh%rr(1, mesh%np))
      dx = (x1 - x0) / (mesh%np - 1)
      DO i = 1, mesh%np
         mesh%rr(1, i) = x0 + (i - 1) * dx
      END DO
   write(*,*) 'ok112'

      mesh%nps = 2
      mesh%mes = 2
      ALLOCATE(mesh%jjs(1, mesh%mes))
      mesh%jjs(1, 1) = 1
      mesh%jjs(1, mesh%mes) = mesh%np
      ALLOCATE(mesh%sides(mesh%mes))
      READ(in_unit, *) mesh%sides
   write(*,*) 'ok113'

      IF (inputs%type_fe==1) THEN
         CALL GAUSS_POINT_1d(mesh)
      ELSE
         WRITE(*, *) ' BUG load_mesh_1d: FE not programmed yet'
         STOP
      END IF
   write(*,*) 'ok114'

   END SUBROUTINE load_mesh_1d

   SUBROUTINE GAUSS_POINT_1d(mesh)
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      REAL(KIND = 8) :: one = 1.d0, two = 2.d0, three = 3.d0
      REAL(KIND = 8) :: f1, f2, x, dhatxdx
      REAL(KIND = 8), DIMENSION(2) :: xx
      INTEGER :: l, m
      f1(x) = (one - x) / two
      f2(x) = (x + one) / two

      mesh%gauss%k_d = 1
      mesh%gauss%n_w = 2
      mesh%gauss%l_G = 2
      ALLOCATE(mesh%gauss%ww(mesh%gauss%n_w, mesh%gauss%l_G))
      ALLOCATE(mesh%gauss%dw(k_dim, mesh%gauss%n_w, mesh%gauss%l_G, mesh%me))
      ALLOCATE(mesh%gauss%rj(mesh%gauss%l_G, mesh%me))

      xx(1) = - one / SQRT(three)
      xx(2) = + one / SQRT(three)

      DO l = 1, mesh%gauss%l_G
         mesh%gauss%ww(1, l) = f1(xx(l))
         mesh%gauss%ww(2, l) = f2(xx(l))
      ENDDO

      DO m = 1, mesh%me
         dhatxdx = 2 / ABS(mesh%rr(1, mesh%jj(1, m)) - mesh%rr(1, mesh%jj(2, m)))
         DO l = 1, mesh%gauss%l_G
            mesh%gauss%dw(1, 1, l, m) = - one / two * dhatxdx
            mesh%gauss%dw(1, 2, l, m) = one / two * dhatxdx
            mesh%gauss%rj(l, m) = 1 / dhatxdx
         END DO
      END DO

   END SUBROUTINE GAUSS_POINT_1d

END MODULE mesh_1d
