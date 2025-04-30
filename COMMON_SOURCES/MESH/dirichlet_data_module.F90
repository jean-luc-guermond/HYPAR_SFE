MODULE dirichlet_type_module
   IMPLICIT NONE
   TYPE dirichlet_bc
      CHARACTER(100) :: name
      INTEGER :: nb_sides
      INTEGER, DIMENSION(:), POINTER :: list_sides
      INTEGER, DIMENSION(:), POINTER :: jsd
   CONTAINS
      PROCEDURE, PUBLIC :: set => dirichlet_nodes_parallel
      PROCEDURE, PRIVATE :: read_dirichlet_data
   END type dirichlet_bc
CONTAINS

 SUBROUTINE read_dirichlet_data(this)
      IMPLICIT NONE
      CLASS(dirichlet_bc) :: this
      INTEGER, PARAMETER :: in_unit = 21
      CHARACTER(LEN = 100) :: argument
      LOGICAL :: test
      !===Initialize data to zero and false by default

      OPEN(UNIT = in_unit, FILE = "data", FORM = 'formatted', STATUS = 'unknown')

      CALL find_string(21, '===How many pieces of boundaries for bcs on ' // trim(adjustl(this%name)) // '?===', test)
      IF (test) THEN
         READ (21, *) this%nb_sides
      ELSE
         this%nb_sides = 0
         WRITE(*,*) 'Boundaries for '//  trim(adjustl(this%name)) // ' not found. Set to none by default.'
      END IF
      ALLOCATE(this%list_sides(this%nb_sides))

      IF (this%nb_sides > 0) THEN
         CALL read_until(21, '===List of boundaries for bcs on ' // trim(adjustl(this%name)) // '===')
         READ(21, *) this%list_sides
      END IF

      CLOSE(in_unit)
   END SUBROUTINE read_dirichlet_data

   SUBROUTINE dirichlet_nodes_parallel(this, mesh, name)
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(dirichlet_bc) :: this
      TYPE(mesh_type) :: mesh
      CHARACTER(*) :: name
      LOGICAL, DIMENSION(:), POINTER :: virgin
      INTEGER :: nn, ms, n, p, n_D, nws, n_D_me, k
      LOGICAL :: test

      this%name = name

      CALL this%read_dirichlet_data

      IF (SIZE(this%list_sides)==0) THEN
         ALLOCATE(this%jsd(0))
         RETURN
      END IF

      nws = SIZE(mesh%jjs, 1)
      nn = 0
      ALLOCATE(virgin(mesh%dom_np))
      virgin = .TRUE.
      DO ms = 1, mesh%dom_mes
         IF (MINVAL(ABS(mesh%sides(ms) - this%list_sides))/=0) CYCLE
         DO n = 1, nws
            p = mesh%jjs(n, ms)
            IF (p>mesh%dom_np) CYCLE
            IF (virgin(p)) THEN
               virgin(p) = .FALSE.
               nn = nn + 1
            END IF
         END DO
      END DO
      n_D_me = nn
      DO ms = 1, mesh%nis
         test = .false.
         DO k = 1, mesh%gauss%n_ws
            test = test .OR. MINVAL(ABS(mesh%isolated_interfaces(ms, k) - this%list_sides))==0
         END DO
         IF (test) THEN
            nn = nn + 1
         END IF
      END DO
      n_D = nn
      ALLOCATE(this%jsd(n_D))
      nn = 0
      virgin = .TRUE.
      DO ms = 1, mesh%dom_mes
         IF (MINVAL(ABS(mesh%sides(ms) - this%list_sides))/=0) CYCLE
         DO n = 1, nws
            p = mesh%jjs(n, ms)
            IF (p>mesh%dom_np) CYCLE
            IF (virgin(p)) THEN
               virgin(p) = .FALSE.
               nn = nn + 1
               this%jsd(nn) = mesh%jjs(n, ms)
            END IF
         END DO
      END DO
      DO ms = 1, mesh%nis
         test = .false.
         DO k = 1, mesh%gauss%n_ws
            test = test .OR. MINVAL(ABS(mesh%isolated_interfaces(ms, k) - this%list_sides))==0
         END DO
         IF (test) THEN
            nn = nn + 1
            this%jsd(nn) = mesh%isolated_jjs(nn - n_D_me) - mesh%loc_to_glob(1) + 1
         END IF
      END DO
      DEALLOCATE(virgin)
   END SUBROUTINE dirichlet_nodes_parallel

END MODULE dirichlet_type_module
