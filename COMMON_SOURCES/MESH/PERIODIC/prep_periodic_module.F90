!
!Authors Jean-Luc Guermond, Copyrights 1996
!
MODULE prep_periodic_module
   USE def_type_mesh
   USE def_type_periodic
   USE dyn_line_type
   USE input_periodic_data
   USE periodic_data_module
   IMPLICIT NONE

   PUBLIC :: prep_periodic
   PRIVATE

CONTAINS
   SUBROUTINE prep_periodic(mesh, periodic, opt_nb_bloc)
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      TYPE(periodic_type) :: periodic
      INTEGER, OPTIONAL :: opt_nb_bloc

      CALL read_periodic_data('data')

      IF (PRESENT(opt_nb_bloc)) THEN
         CALL prep_periodic_bloc(periodic_data, mesh, periodic, opt_nb_bloc)
      ELSE
         CALL prep_periodic_scal(periodic_data, mesh, periodic)
      END IF

   END SUBROUTINE prep_periodic

   SUBROUTINE prep_periodic_scal(my_periodic, mesh, periodic)
      !=========================================
      USE character_strings
      IMPLICIT NONE
      TYPE(periodic_data_type), INTENT(IN) :: my_periodic
      TYPE(mesh_type) :: mesh
      TYPE(periodic_type) :: periodic
      INTEGER, DIMENSION(:), POINTER :: list_loc, perlist_loc, list_dom, perlist_dom
      INTEGER :: n, side1, side2, n_b, nx, i
      REAL(KIND = 8), DIMENSION(:), POINTER :: e

      WRITE (*, *) 'Loading periodic-data file ...'

      IF (mesh%np == 0) THEN
         WRITE(*, *) 'no mesh on this proc'
         RETURN
      END IF

      ALLOCATE(e(SIZE(my_periodic%vect_e, 1)))

      periodic%n_bord = my_periodic%nb_periodic_pairs
      IF (periodic%n_bord .GT. 20) THEN
         WRITE(*, *) 'PREP_MESH_PERIODIC: trop de bords periodiques'
         STOP
      END IF

      DO n = 1, periodic%n_bord

         side1 = my_periodic%list_periodic(1, n)
         side2 = my_periodic%list_periodic(2, n)
         e = my_periodic%vect_e(:, n)

         CALL list_periodic(mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
              list_loc, perlist_loc)

         n_b = SIZE(list_loc)
         ALLOCATE(list_dom(n_b), perlist_dom(n_b))
         nx = 0
         DO i = 1, n_b
            IF (MAX(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
               nx = nx + 1
               list_dom(nx) = list_loc(i)
               perlist_dom(nx) = perlist_loc(i)
            ELSE IF (MIN(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
               WRITE(*, *) 'BUG in prep_periodic_scal'
               STOP
            END IF
         END DO
         IF (n_b /= nx) WRITE(*, *) 'WARNING, I have removed', n_b - nx, ' periodic pairs in prep_periodic_scal'
         n_b = nx

         ALLOCATE (periodic%list(n)%DIL(n_b), periodic%perlist(n)%DIL(n_b))
         periodic%list(n)%DIL = list_dom(1:n_b)
         periodic%perlist(n)%DIL = perlist_dom(1:n_b)

         DEALLOCATE(list_loc, perlist_loc, list_dom, perlist_dom)
      END DO

      DEALLOCATE(e)

      WRITE (*, *) 'Treatment of periodic-data done'

   END SUBROUTINE prep_periodic_scal

   SUBROUTINE prep_periodic_bloc(my_periodic, mesh, periodic, nb_bloc)
      !=========================================
      USE character_strings
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(periodic_data_type), INTENT(IN) :: my_periodic
      TYPE(mesh_type) :: mesh
      TYPE(periodic_type) :: periodic
      INTEGER, INTENT(IN) :: nb_bloc
      INTEGER, DIMENSION(:), POINTER :: list_loc, perlist_loc, list_dom, perlist_dom
      INTEGER :: n, side1, side2, nsize, n_b
      INTEGER :: k, k_deb, k_fin, nx, i
      REAL(KIND = 8), DIMENSION(2) :: e

      WRITE (*, *) 'Loading periodic-data file ...'

      IF (mesh%np == 0) THEN
         WRITE(*, *) 'no mesh on this proc'
         RETURN
      END IF

      periodic%n_bord = my_periodic%nb_periodic_pairs
      IF (periodic%n_bord .GT. 20) THEN
         WRITE(*, *) 'PREP_MESH_PERIODIC: trop de bords periodiques'
         STOP
      END IF

      DO n = 1, periodic%n_bord

         side1 = my_periodic%list_periodic(1, n)
         side2 = my_periodic%list_periodic(2, n)
         e = my_periodic%vect_e(:, n)

         CALL list_periodic(mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
              list_loc, perlist_loc)

         !n_b = SIZE(list_loc)
         n_b = SIZE(perlist_loc)
         ALLOCATE(list_dom(n_b), perlist_dom(n_b))
         nx = 0
         DO i = 1, n_b
            IF (MAX(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
               nx = nx + 1
               list_dom(nx) = list_loc(i)
               perlist_dom(nx) = perlist_loc(i)
            ELSE IF (MIN(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
               WRITE(*, *) 'BUG in prep_periodic_bloc'
               STOP
            END IF
         END DO
         IF (n_b /= nx) WRITE(*, *) 'WARNING, I have removed', n_b - nx, ' periodic pairs in prep_periodic_bloc'
         n_b = nx

         nsize = nb_bloc * n_b !SIZE(list_loc) !n_b

         ALLOCATE(periodic%list(n)%DIL(nsize), periodic%perlist(n)%DIL(nsize))

         DO k = 1, nb_bloc
            k_deb = (k - 1) * n_b + 1
            k_fin = k * n_b
            periodic%list(n)%DIL(k_deb:k_fin) = list_dom(1:n_b) + (k - 1) * mesh%dom_np ! First bloc
            periodic%perlist(n)%DIL(k_deb:k_fin) = perlist_dom(1:n_b) + (k - 1) * mesh%dom_np ! First bloc
         END DO

         DEALLOCATE(list_loc, perlist_loc, list_dom, perlist_dom)

      END DO

      WRITE (*, *) 'Treatment of periodic-data done'

   END SUBROUTINE prep_periodic_bloc
   !jan 29 2007

   SUBROUTINE list_periodic(np, jjs, sides, rr, side1, side2, e, list_out, perlist_out)
      !============================================================================
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: np
      INTEGER, DIMENSION(:, :), INTENT(IN) :: jjs
      INTEGER, DIMENSION(:), INTENT(IN) :: sides
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      INTEGER, INTENT(IN) :: side1, side2
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: e
      INTEGER, DIMENSION(:), POINTER :: list_out, perlist_out
      INTEGER, DIMENSION(:), ALLOCATABLE :: list, perlist
      LOGICAL, DIMENSION(np) :: virgin
      REAL(KIND = 8), DIMENSION(SIZE(rr, 1)) :: ri
      INTEGER :: ms, ns, i, j, long, inter
      REAL(KIND = 8) :: r, epsilon = 1.d-9
      LOGICAL :: verif

      IF (ALLOCATED(list))    DEALLOCATE(list)
      IF (ALLOCATED(perlist)) DEALLOCATE(perlist)

      ALLOCATE (list(np), perlist(np))
      virgin = .TRUE.

      i = 0; j = 0
      DO ms = 1, SIZE(sides)

         IF (sides(ms) .EQ. side1) THEN
            DO ns = 1, SIZE(jjs, 1)
               IF (virgin(jjs(ns, ms))) THEN
                  i = i + 1
                  list(i) = jjs(ns, ms)
                  virgin(jjs(ns, ms)) = .FALSE.
               END IF
            END DO
         ELSE IF (sides(ms) .EQ. side2) THEN
            DO ns = 1, SIZE(jjs, 1)
               IF (virgin(jjs(ns, ms))) THEN
                  j = j + 1
                  perlist(j) = jjs(ns, ms)
                  virgin(jjs(ns, ms)) = .FALSE.
               END IF
            END DO

         END IF

      END DO

      IF (i .NE. j) THEN
         WRITE(*, *) ' FEM_PERIODIC: side1 and side2 have', &
              ' different numbers of points'
         STOP
      END IF
      long = i

      DO i = 1, long
         ri = rr(:, list(i)) + e(:)
         verif = .FALSE.
         !if (i==2) stop
         DO j = i, long
            r = SUM(ABS(ri - rr(:, perlist(j))))
            !if (i==1) write(*,*) ' r',r,'j',  j
            IF (r .LE. epsilon) THEN
               inter = perlist(i)
               perlist(i) = perlist(j)
               perlist(j) = inter
               verif = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.verif) THEN
            WRITE(*, *) ' BUG dans  data_periodic ou le maillage:', &
                 ' side1 + e /= side2'
            WRITE(*, *) ' i = ', i
            !         STOP
         END IF
      END DO

      ALLOCATE (list_out(long))
      list_out(1:long) = list(1:long)
      ALLOCATE (perlist_out(long))
      perlist_out(1:long) = perlist(1:long)

   END SUBROUTINE list_periodic

END MODULE prep_periodic_module
