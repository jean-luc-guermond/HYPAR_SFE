MODULE periodic_data_module
   USE dyn_line_type
   USE def_type_mesh
    
   IMPLICIT NONE
   INTEGER, PARAMETER, PRIVATE :: rec_length=200

   TYPE argument_periodic_type
      CHARACTER(rec_length)                    :: nb_bords
      CHARACTER(rec_length)                    :: list_periodic
   END TYPE argument_periodic_type
   
   TYPE periodic_type
      CHARACTER(100)                           :: name
      INTEGER                                  :: nb_bords      = 0
      INTEGER, DIMENSION(:, :), POINTER        :: list_periodic
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: vect_e
      TYPE(dyn_int_line), DIMENSION(20)        :: list
      TYPE(dyn_int_line), DIMENSION(20)        :: perlist
      INTEGER, POINTER, DIMENSION(:)           :: pnt
   CONTAINS
      PROCEDURE, PUBLIC :: read => read_periodic_data
      PROCEDURE, PUBLIC :: set => prep_periodic
   END TYPE periodic_type

CONTAINS
   SUBROUTINE read_periodic_data(this, name)
      USE character_strings
      USE space_dim
      USE petsc
      IMPLICIT NONE
      INTEGER, PARAMETER :: in_unit = 21, list_length=200, length_template_begin=28, length_template_end=26
      CHARACTER(LEN=length_template_begin), PARAMETER :: template_begin_section = '%%% BEGIN SECTION: PERIODIC '  
      CHARACTER(LEN=length_template_end),   PARAMETER :: template_end_section   = '%%% END SECTION: PERIODIC '
      CHARACTER(LEN=4),                     PARAMETER :: template_tip           = ' %%%'
      CHARACTER(LEN=:), ALLOCATABLE                   :: begin_section, end_section
      
      CHARACTER(LEN=rec_length), DIMENSION(list_length):: list, record
      CHARACTER(*) :: name
      CHARACTER(LEN=rec_length)       :: string, string_default
      
      CLASS(periodic_type), INTENT(INOUT) :: this
      
      LOGICAL :: okay
      INTEGER :: k, length_begin, length_end, length_name
      INTEGER :: rank, ierr, record_size, i_list, j, i

      TYPE(argument_periodic_type)        :: argument_data

      !===Initialize data to zero and false by default
      list = ""
      record = ""
      CALL MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

      this%name = name
      argument_data%nb_bords = '=== How many pieces of periodic boundary on ' // trim(adjustl(this%name)) // '? ==='
      argument_data%list_periodic = '=== Indices of periodic boundaries and corresponding vectors on ' // trim(adjustl(this%name)) // '? ==='

      !=== dynamic BEGIN/END in data
      length_name = LEN(this%name)
      length_begin = length_template_begin + length_name + LEN(template_tip) 
      length_end   = length_template_end   + length_name + LEN(template_tip)
      ALLOCATE(CHARACTER(LEN=length_begin) :: begin_section)
      ALLOCATE(CHARACTER(LEN=length_end  ) :: end_section)
      begin_section = template_begin_section // TRIM(ADJUSTL(this%name)) // template_tip
      end_section = template_end_section // TRIM(ADJUSTL(this%name)) // template_tip
      
      !===Initializing record
      CALL read_data_in_record(record_size, record, begin_section, end_section)

      !===Now we reorganize record
  
      i_list = 1
      list(i_list) = begin_section

      !===Initialize data to zero and false by default
      this%nb_bords = 0
      WRITE(string_default,*) this%nb_bords
      CALL compare_string(record, list, argument_data%nb_bords, string_default, okay, i_list, j)
      IF (okay) THEN
         READ(list(i_list),*) this%nb_bords  
      END IF

      ALLOCATE(this%list_periodic(2, this%nb_bords))
      ALLOCATE(this%vect_e(k_dim, this%nb_bords))
      string = argument_data%list_periodic
      string_default = "0 0 0.d0 0.d0"
      IF (this%nb_bords > 0) THEN
         okay = .FALSE.
         i_list = i_list+1
         list(i_list) = string
         DO i = 1, SIZE(record)
            !=== detecting if there is Periodic BC
            IF (TRIM(ADJUSTL(record(i)))==list(i_list)) THEN
               j = i
               record(j) = ''
               okay = .TRUE.
               EXIT
            END IF
         END DO
         IF (okay) THEN
            !=== reading all Periodic BC if detected
            DO k=1, this%nb_bords
               i_list = i_list + 1
               list(i_list) = record(j+k)
               record(j+k) = ''
               READ(list(i_list), *) this%list_periodic(:, k), this%vect_e(:, k)
            END DO
         ELSE
            !=== default value if no Periodic BC detected
            i_list = i_list+1
            list(i_list) = string_default
         END IF
      ELSE
         i_list = i_list+1
         list(i_list) = string
         i_list = i_list+1
         list(i_list) = string_default
      END IF

      i_list = i_list+1
      list(i_list) = end_section
     

      !===Closing unit
      CALL rewrite_data_from_list_record(rank, list, record, i_list, record_size)

      DEALLOCATE(begin_section)
      DEALLOCATE(end_section)

   END SUBROUTINE read_periodic_data

   SUBROUTINE prep_periodic(this, mesh, opt_nb_bloc)
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      CLASS(periodic_type), INTENT(INOUT)  :: this
      INTEGER, OPTIONAL :: opt_nb_bloc

      IF (PRESENT(opt_nb_bloc)) THEN
         CALL prep_periodic_bloc(this, mesh, opt_nb_bloc)
      ELSE
         CALL prep_periodic_scal(this, mesh)
      END IF

   END SUBROUTINE prep_periodic

   SUBROUTINE prep_periodic_scal(periodic, mesh)
      !=========================================
      USE character_strings
      IMPLICIT NONE
      TYPE(periodic_type) :: periodic
      TYPE(mesh_type) :: mesh
      INTEGER, DIMENSION(:), POINTER :: list_loc, perlist_loc, list_dom, perlist_dom
      INTEGER :: n, side1, side2, n_b, nx, i
      REAL(KIND = 8), DIMENSION(:), POINTER :: e

      IF (mesh%np == 0) THEN
         !WRITE(*, *) 'no mesh on this proc'
         RETURN
      END IF

      ALLOCATE(e(SIZE(periodic%vect_e, 1)))

      IF (periodic%nb_bords .GT. 20) THEN
         WRITE(*, *) 'PREP_MESH_PERIODIC: trop de bords periodiques'
         STOP
      END IF

      DO n = 1, periodic%nb_bords

         side1 = periodic%list_periodic(1, n)
         side2 = periodic%list_periodic(2, n)
         e = periodic%vect_e(:, n)

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
               WRITE(*, *) 'BUG in prep_periodic_scal, one of the boundary point is not attributed the same processor.'
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

   END SUBROUTINE prep_periodic_scal

   SUBROUTINE prep_periodic_bloc(periodic, mesh, nb_bloc)
      !=========================================
      USE character_strings
      USE def_type_mesh
      IMPLICIT NONE
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

      IF (periodic%nb_bords .GT. 20) THEN
         WRITE(*, *) 'PREP_MESH_PERIODIC: trop de bords periodiques'
         STOP
      END IF

      DO n = 1, periodic%nb_bords

         side1 = periodic%list_periodic(1, n)
         side2 = periodic%list_periodic(2, n)
         e = periodic%vect_e(:, n)

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

END MODULE periodic_data_module
