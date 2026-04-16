MODULE periodic_data_module
   USE dyn_line_type
    
   IMPLICIT NONE
   
   TYPE periodic_type
      CHARACTER(100)                           :: name
      INTEGER                                  :: nb_bords      = 0
      INTEGER, DIMENSION(:, :), POINTER        :: list_periodic
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: vect_e
      TYPE(dyn_int_line), DIMENSION(20)        :: list
      TYPE(dyn_int_line), DIMENSION(20)        :: perlist
      INTEGER, POINTER, DIMENSION(:)           :: pnt
   END TYPE periodic_type

CONTAINS
!!$ 
!!$   SUBROUTINE init_periodic_bc_data(this, name, section_name)
!!$      CLASS(periodic_type), INTENT(INOUT) :: this
!!$      CHARACTER(LEN=*), INTENT(IN) :: name
!!$      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name
!!$      IF (PRESENT(section_name)) THEN
!!$         CALL this%read(name, section_name)
!!$      ELSE
!!$         CALL this%read(name)
!!$      END IF
!!$   END SUBROUTINE init_periodic_bc_data
!!$
!!$   SUBROUTINE read_periodic_bc_data(this, name, section_name)
!!$      USE petsc
!!$      USE character_strings
!!$!VB 2/04/2026
!!$      ! USE space_dim
!!$      USE mesh_parameters
!!$!VB 2/04/2026
!!$      IMPLICIT NONE
!!$      CHARACTER(*), INTENT(IN)            :: name
!!$      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name
!!$
!!$      
!!$      CLASS(periodic_type), INTENT(INOUT) :: this
!!$      TYPE(argument_periodic_type)        :: argument_data
!!$
!!$      this%name = name
!!$      argument_data%nb_bords = '=== How many pieces of periodic boundary on ' // trim(adjustl(this%name)) // '? ==='
!!$      argument_data%list_periodic = '=== Indices of periodic boundaries and corresponding vectors on ' // trim(adjustl(this%name)) // '? ==='
!!$
!!$!================
!!$!=== MANDATORY Reading all data file
!!$!================
!!$      IF (PRESENT(section_name)) THEN
!!$         CALL read_data_init_list(section_name)
!!$      ELSE
!!$         CALL read_data_init_list()
!!$      END IF
!!$
!!$!================
!!$!=== We now find the relevant information for this specific periodic BC
!!$!================
!!$
!!$   !=== Number of periodic boundaries 
!!$      CALL read_data(argument_data%nb_bords, this%nb_bords)
!!$
!!$   !=== List of periodic boundaries (has its specific subroutine, see character_strings.F90 module) 
!!$      ALLOCATE(this%list_periodic(2, this%nb_bords))
!!$      ALLOCATE(this%vect_e(mesh_data_info%k_dim, this%nb_bords))
!!$      CALL read_periodic_data(argument_data%list_periodic, this%nb_bords, this%list_periodic, this%vect_e)
!!$
!!$!================
!!$!=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
!!$!================
!!$     CALL finalize_rewrite_data
!!$
!!$   END SUBROUTINE read_periodic_bc_data
!!$
!!$   SUBROUTINE prep_periodic(this, mesh, opt_nb_bloc)
!!$      IMPLICIT NONE
!!$      TYPE(mesh_type) :: mesh
!!$      CLASS(periodic_type), INTENT(INOUT)  :: this
!!$      INTEGER, OPTIONAL :: opt_nb_bloc
!!$
!!$      IF (PRESENT(opt_nb_bloc)) THEN
!!$         CALL prep_periodic_bloc(this, mesh, opt_nb_bloc)
!!$      ELSE
!!$         CALL prep_periodic_scal(this, mesh)
!!$      END IF
!!$
!!$   END SUBROUTINE prep_periodic
!!$
!!$   SUBROUTINE prep_periodic_scal(periodic, mesh)
!!$      !=========================================
!!$      USE character_strings
!!$      IMPLICIT NONE
!!$      TYPE(periodic_type) :: periodic
!!$      TYPE(mesh_type) :: mesh
!!$      INTEGER, DIMENSION(:), POINTER :: list_loc, perlist_loc, list_dom, perlist_dom
!!$      INTEGER :: n, side1, side2, n_b, nx, i
!!$      REAL(KIND = 8), DIMENSION(:), POINTER :: e
!!$
!!$      IF (mesh%np == 0) THEN
!!$         !WRITE(*, *) 'no mesh on this proc'
!!$         RETURN
!!$      END IF
!!$
!!$      ALLOCATE(e(SIZE(periodic%vect_e, 1)))
!!$
!!$      IF (periodic%nb_bords .GT. 20) THEN
!!$         WRITE(*, *) 'PREP_MESH_PERIODIC: too many periodic pieces'
!!$         STOP
!!$      END IF
!!$
!!$      DO n = 1, periodic%nb_bords
!!$
!!$         side1 = periodic%list_periodic(1, n)
!!$         side2 = periodic%list_periodic(2, n)
!!$         e = periodic%vect_e(:, n)
!!$
!!$         CALL list_periodic(mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
!!$              list_loc, perlist_loc)
!!$
!!$         n_b = SIZE(list_loc)
!!$         ALLOCATE(list_dom(n_b), perlist_dom(n_b))
!!$         nx = 0
!!$         DO i = 1, n_b
!!$            IF (MAX(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
!!$               nx = nx + 1
!!$               list_dom(nx) = list_loc(i)
!!$               perlist_dom(nx) = perlist_loc(i)
!!$            ELSE IF (MIN(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
!!$               WRITE(*, *) 'BUG in prep_periodic_scal, one of the boundary point is not attributed the same processor.'
!!$               STOP
!!$            END IF
!!$         END DO
!!$         IF (n_b /= nx) WRITE(*, *) 'WARNING, I have removed', n_b - nx, ' periodic pairs in prep_periodic_scal'
!!$         n_b = nx
!!$
!!$         ALLOCATE (periodic%list(n)%DIL(n_b), periodic%perlist(n)%DIL(n_b))
!!$         periodic%list(n)%DIL = list_dom(1:n_b)
!!$         periodic%perlist(n)%DIL = perlist_dom(1:n_b)
!!$
!!$         DEALLOCATE(list_loc, perlist_loc, list_dom, perlist_dom)
!!$      END DO
!!$
!!$      DEALLOCATE(e)
!!$
!!$   END SUBROUTINE prep_periodic_scal
!!$
!!$   SUBROUTINE prep_periodic_bloc(periodic, mesh, nb_bloc)
!!$      !=========================================
!!$      USE character_strings
!!$      USE def_type_mesh
!!$      IMPLICIT NONE
!!$      TYPE(mesh_type) :: mesh
!!$      TYPE(periodic_type) :: periodic
!!$      INTEGER, INTENT(IN) :: nb_bloc
!!$      INTEGER, DIMENSION(:), POINTER :: list_loc, perlist_loc, list_dom, perlist_dom
!!$      INTEGER :: n, side1, side2, nsize, n_b
!!$      INTEGER :: k, k_deb, k_fin, nx, i
!!$      REAL(KIND = 8), DIMENSION(2) :: e
!!$
!!$      WRITE (*, *) 'Loading periodic-data file ...'
!!$
!!$      IF (mesh%np == 0) THEN
!!$         WRITE(*, *) 'no mesh on this proc'
!!$         RETURN
!!$      END IF
!!$
!!$      IF (periodic%nb_bords .GT. 20) THEN
!!$         WRITE(*, *) 'PREP_MESH_PERIODIC: trop de bords periodiques'
!!$         STOP
!!$      END IF
!!$
!!$      DO n = 1, periodic%nb_bords
!!$
!!$         side1 = periodic%list_periodic(1, n)
!!$         side2 = periodic%list_periodic(2, n)
!!$         e = periodic%vect_e(:, n)
!!$
!!$         CALL list_periodic(mesh%np, mesh%jjs, mesh%sides, mesh%rr, side1, side2, e, &
!!$              list_loc, perlist_loc)
!!$
!!$         !n_b = SIZE(list_loc)
!!$         n_b = SIZE(perlist_loc)
!!$         ALLOCATE(list_dom(n_b), perlist_dom(n_b))
!!$         nx = 0
!!$         DO i = 1, n_b
!!$            IF (MAX(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
!!$               nx = nx + 1
!!$               list_dom(nx) = list_loc(i)
!!$               perlist_dom(nx) = perlist_loc(i)
!!$            ELSE IF (MIN(list_loc(i), perlist_loc(i)) .LE. mesh%dom_np) THEN
!!$               WRITE(*, *) 'BUG in prep_periodic_bloc'
!!$               STOP
!!$            END IF
!!$         END DO
!!$         IF (n_b /= nx) WRITE(*, *) 'WARNING, I have removed', n_b - nx, ' periodic pairs in prep_periodic_bloc'
!!$         n_b = nx
!!$
!!$         nsize = nb_bloc * n_b !SIZE(list_loc) !n_b
!!$
!!$         ALLOCATE(periodic%list(n)%DIL(nsize), periodic%perlist(n)%DIL(nsize))
!!$
!!$         DO k = 1, nb_bloc
!!$            k_deb = (k - 1) * n_b + 1
!!$            k_fin = k * n_b
!!$            periodic%list(n)%DIL(k_deb:k_fin) = list_dom(1:n_b) + (k - 1) * mesh%dom_np ! First bloc
!!$            periodic%perlist(n)%DIL(k_deb:k_fin) = perlist_dom(1:n_b) + (k - 1) * mesh%dom_np ! First bloc
!!$         END DO
!!$
!!$         DEALLOCATE(list_loc, perlist_loc, list_dom, perlist_dom)
!!$
!!$      END DO
!!$
!!$      WRITE (*, *) 'Treatment of periodic-data done'
!!$
!!$   END SUBROUTINE prep_periodic_bloc
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
            WRITE(*, *) ' BUG in data_periodic or in mesh:', &
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

END MODULE periodic_data_module
