MODULE mesh_1d
   USE def_type_mesh
   USE input_mesh_data
   USE space_dim
   uSE mesh_tools
   PUBLIC :: load_mesh_1d, GAUSS_POINT_1d
   PRIVATE
CONTAINS

   SUBROUTINE load_mesh_1d(mesh)
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      INTEGER, PARAMETER :: in_unit = 30
      INTEGER :: type_fe
      INTEGER :: i, n, m, nb_procs
      REAL(KIND = 8) :: x0, x1, dx
      nb_procs = 1

      OPEN(in_unit, FILE = TRIM(ADJUSTL(mesh_data%directory)) // &
           '/' // TRIM(ADJUSTL(mesh_data%file_name)), FORM = 'formatted')
      READ(in_unit, *) mesh%me, mesh_data%type_fe
      READ(in_unit, *) x0, x1
      type_fe = mesh_data%type_fe
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

      ALLOCATE(mesh%rr(1, mesh%np))
      dx = (x1 - x0) / (mesh%np - 1)
      DO i = 1, mesh%np
         mesh%rr(1, i) = x0 + (i - 1) * dx
      END DO

      mesh%nps = 2
      mesh%mes = 2
      ALLOCATE(mesh%jjs(1, mesh%mes))
      mesh%jjs(1, 1) = 1
      mesh%jjs(1, mesh%mes) = mesh%np
      ALLOCATE(mesh%neighs(mesh%mes))
      mesh%neighs(1) = 1
      mesh%neighs(2) = mesh%np - 1
      ALLOCATE(mesh%sides(mesh%mes))
      READ(in_unit, *) mesh%sides

      ALLOCATE(mesh%loc_to_glob(mesh%np))
      DO i = 1, mesh%np
         mesh%loc_to_glob(i) = i
      END DO

      ALLOCATE(mesh%i_d(mesh%me))
      mesh%i_d = 1

      mesh%nis = 2
      ALLOCATE(mesh%isolated_jjs(mesh%nis), mesh%isolated_interfaces(mesh%nis, 1))
      mesh%isolated_jjs(1) = 1
      mesh%isolated_jjs(2) = mesh%np
      mesh%isolated_interfaces(1, 1) = 1
      mesh%isolated_interfaces(2, 1) = 2

      mesh%mi = 0
      mesh%medge = mesh%me
      mesh%medges = 0
      mesh%mes_extra = 0
      mesh%mes_int = 0
      mesh%dom_np = mesh%np
      mesh%dom_me = mesh%me
      mesh%dom_mes = mesh%mes
      mesh%mextra = 0

      ALLOCATE(mesh%iis(0, 0))
      ALLOCATE(mesh%jj_extra(2, mesh%mextra), mesh%jce_extra(0, mesh%medge), &
           mesh%jjs_extra(0, mesh%mes_extra))
      ALLOCATE(mesh%jjs_int(0, 0), mesh%jcc_extra(mesh%mextra))
      ALLOCATE(mesh%jees(0), mesh%jecs(0))
      ALLOCATE(mesh%jji(0, 0, 0), mesh%jjsi(0, 0), mesh%j_s(0))
      ALLOCATE(mesh%rrs_extra(1, 2, 0))
      ALLOCATE(mesh%neighi(0, 0))
      ALLOCATE(mesh%sides_extra(mesh%mes_extra), mesh%neighs_extra(mesh%mes_extra))
      ALLOCATE(mesh%sides_int(mesh%mes_int), mesh%neighs_int(2, mesh%mes_int))
      ALLOCATE(mesh%disp(nb_procs + 1), mesh%disedge(nb_procs + 1), mesh%discell(nb_procs + 1))
      ALLOCATE(mesh%domnp(nb_procs), mesh%domedge(nb_procs), mesh%domcell(nb_procs))

      ALLOCATE(mesh%jce(1, mesh%medge))
      DO i = 1, mesh%medge
         mesh%jce(1, i) = i
      END DO
      mesh%disp(1) = 1
      mesh%disp(2) = mesh%dom_np + 1
      mesh%discell(1) = 1
      mesh%discell(2) = mesh%dom_me + 1
      mesh%disedge(1) = 1
      mesh%disedge(2) = mesh%medge + 1
      mesh%domnp(1) = mesh%dom_np
      mesh%domcell(1) = mesh%dom_me
      mesh%domedge(1) = mesh%medge
      IF (mesh_data%type_fe==1) THEN
         CALL GAUSS_POINT_1d(mesh)
      ELSE
         WRITE(*, *) ' BUG load_mesh_1d: FE not programmed yet'
         STOP
      END IF

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
      mesh%gauss%n_ws = 1
      mesh%gauss%l_Gs = 0
      mesh%gauss%n_e = 1
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
