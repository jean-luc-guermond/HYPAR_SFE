MODULE dirichlet_type_module
   IMPLICIT NONE

   INTEGER, PARAMETER, PRIVATE :: rec_length=200

   TYPE argument_dirichlet_bc
      CHARACTER(rec_length) :: nb_sides
      CHARACTER(rec_length) :: list_sides
   END TYPE argument_dirichlet_bc

   TYPE dirichlet_bc
      CHARACTER(rec_length)          :: name
      INTEGER                        :: nb_sides   = 0
      INTEGER, DIMENSION(:), POINTER :: list_sides
      INTEGER, DIMENSION(:), POINTER :: jsd
    CONTAINS
      PROCEDURE, PUBLIC :: set => dirichlet_nodes_local !===With ghosted nodes
      PROCEDURE, PUBLIC :: read => read_dirichlet_data
      PROCEDURE, PUBLIC :: init => init_dirichlet_data
      !!PROCEDURE, PUBLIC :: set => dirichlet_nodes_parallel !===Without ghosted nodes
   END type dirichlet_bc
CONTAINS

   SUBROUTINE init_dirichlet_data(this, section_name)
      CLASS(dirichlet_bc), INTENT(INOUT) :: this
      CHARACTER(LEN=rec_length), OPTIONAL, INTENT(IN) :: section_name
      IF (PRESENT(section_name)) THEN
         CALL this%read(section_name)
      ELSE
         CALL this%read
      END IF
   END SUBROUTINE init_dirichlet_data

   SUBROUTINE read_dirichlet_data(this, section_name)
      USE character_strings
      IMPLICIT NONE

      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name

      CLASS(dirichlet_bc), INTENT(INOUT)              :: this
      TYPE(argument_dirichlet_bc)                     :: argument_data

      !===Initialize data arguments (depends on the name of the BC)
      argument_data%nb_sides = '=== How many pieces of Dirichlet boundaries for bcs on ' // trim(adjustl(this%name)) // '? ==='
      argument_data%list_sides = '=== List of Dirichlet boundaries for bcs on ' // trim(adjustl(this%name)) // ' ==='

!================
!=== MANDATORY Reading all data file
!================
      IF (PRESENT(section_name)) THEN
         CALL read_data_init_list(section_name)
      ELSE
         CALL read_data_init_list()
      END IF

!================
!=== We now find the relevant information for this specific DIRICHLET BC
!================

    !=== number of sides where to impose the DIRICHLET BC 
      CALL read_data(argument_data%nb_sides, this%nb_sides)

    !=== list of sides where to impose the DIRICHLET BC (special treatment if the list is empty)
      ALLOCATE(this%list_sides(this%nb_sides))

      IF (this%nb_sides>0) THEN
         this%list_sides = 0
      END IF
      CALL read_data(argument_data%list_sides, this%list_sides, opt_skip_data=this%nb_sides==0)
   
!================
!=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
!================
      CALL finalize_rewrite_data

   END SUBROUTINE read_dirichlet_data

   SUBROUTINE dirichlet_nodes_local(this, mesh, name, section_name)
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(dirichlet_bc) :: this
      TYPE(mesh_type) :: mesh
      CHARACTER(*) :: name
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name
      LOGICAL, DIMENSION(:), POINTER :: virgin
      INTEGER :: nn, ms, n, p, n_D, nws, n_D_me, k
      LOGICAL :: test

      this%name = name

      IF (this%name.NE.'whole boundary') THEN
            IF (PRESENT(section_name)) THEN
               CALL this%read(section_name)
            ELSE
               CALL this%read
            END IF
         !  CALL this%read_dirichlet_data
      ELSE
         n = MINVAL(mesh%sides)
         p = MAXVAL(mesh%sides)
         this%nb_sides = p-n+1
         ALLOCATE(this%list_sides(this%nb_sides))
         DO k = n, p
            this%list_sides(k-n+1) = k
         END DO
      END IF

      IF (SIZE(this%list_sides)==0) THEN
         ALLOCATE(this%jsd(0))
         RETURN
      END IF

      nws = SIZE(mesh%jjs, 1)
      nn = 0
      ALLOCATE(virgin(mesh%np)) !===mesh%np instead of mesh%dom_np
      virgin = .TRUE.
      DO ms = 1, mesh%dom_mes
         IF (MINVAL(ABS(mesh%sides(ms) - this%list_sides))/=0) CYCLE
         DO n = 1, nws
            p = mesh%jjs(n, ms)
            IF (p>mesh%np) CYCLE! ===mesh%np instead of mesh%dom_np
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
            IF (p>mesh%np) CYCLE !===mesh%np instead of mesh%dom_np
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
   END SUBROUTINE dirichlet_nodes_local

   SUBROUTINE dirichlet_nodes_parallel(this, mesh, name, section_name)
      USE def_type_mesh
      IMPLICIT NONE
      CLASS(dirichlet_bc) :: this
      TYPE(mesh_type) :: mesh
      CHARACTER(*) :: name
         CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name
      LOGICAL, DIMENSION(:), POINTER :: virgin
      INTEGER :: nn, ms, n, p, n_D, nws, n_D_me, k
      LOGICAL :: test

      this%name = name

      IF (this%name.NE.'whole boundary') THEN
         IF (PRESENT(section_name)) THEN
            CALL this%read(section_name)
         ELSE  
            CALL this%read
         END IF
         !  CALL this%read_dirichlet_data
      ELSE
         n = MINVAL(mesh%sides)
         p = MAXVAL(mesh%sides)
         this%nb_sides = p-n+1
         ALLOCATE(this%list_sides(this%nb_sides))
         DO k = n, p
            this%list_sides(k-n+1) = k
         END DO
      END IF

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
