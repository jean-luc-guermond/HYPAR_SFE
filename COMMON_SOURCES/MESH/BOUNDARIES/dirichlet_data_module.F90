MODULE dirichlet_type_module
   IMPLICIT NONE

   INTEGER, PARAMETER, PRIVATE :: rec_length=200

   TYPE argument_dirichlet_bc
      CHARACTER(200) :: nb_sides
      CHARACTER(200) :: list_sides
   END TYPE argument_dirichlet_bc

   TYPE dirichlet_bc
      CHARACTER(100)                 :: name
      INTEGER                        :: nb_sides   = 0
      INTEGER, DIMENSION(:), POINTER :: list_sides
      INTEGER, DIMENSION(:), POINTER :: jsd
   CONTAINS
      PROCEDURE, PUBLIC       :: set => dirichlet_nodes_parallel
      PROCEDURE, PRIVATE      :: read_dirichlet_data
   END type dirichlet_bc

CONTAINS

  SUBROUTINE read_dirichlet_data(this)
    use petsc
    USE character_strings
    IMPLICIT NONE
    CLASS(dirichlet_bc) :: this

    INTEGER, PARAMETER :: in_unit = 21, list_length=200, length_template_begin=32, length_template_end=30
    CHARACTER(LEN=length_template_begin), PARAMETER :: template_begin_section = '%%% BEGIN SECTION: DIRICHLET BC '  
    CHARACTER(LEN=length_template_end),   PARAMETER :: template_end_section   = '%%% END SECTION: DIRICHLET BC '
    CHARACTER(LEN=4),                     PARAMETER :: template_tip             = ' %%%'
    CHARACTER(LEN=:), ALLOCATABLE   :: begin_section, end_section
    INTEGER            :: length_begin, length_end, length_name
    CHARACTER(LEN=rec_length), DIMENSION(list_length) :: list, record
    CHARACTER(LEN=rec_length)       :: string_default
    LOGICAL :: okay
    INTEGER :: rank, ierr, record_size, i_list, j

    TYPE(argument_dirichlet_bc)        :: argument_data

    !===Initialize data to zero and false by default
    list = ""
    record = ""
    CALL MPI_Comm_rank(petsc_comm_world, rank, ierr)
    
    argument_data%nb_sides = '=== How many pieces of boundaries for bcs on ' // trim(adjustl(this%name)) // '? ==='
    argument_data%list_sides = '=== List of boundaries for bcs on ' // trim(adjustl(this%name)) // ' ==='
    
    !=== dynamic BEGIN/END in data
    length_name = LEN(trim(adjustl(this%name)))
    length_begin = length_template_begin + length_name + len(template_tip) 
    length_end   = length_template_end   + length_name + len(template_tip)
    ALLOCATE(CHARACTER(LEN=length_begin) :: begin_section)
    ALLOCATE(CHARACTER(LEN=length_end  ) :: end_section)
    begin_section = template_begin_section // trim(adjustl(this%name)) // template_tip
    end_section = template_end_section // trim(adjustl(this%name)) // template_tip
     
    CALL read_data_in_record(record_size, record, begin_section, end_section)

    !===Now we reorganize record
    i_list = 1
    list(i_list) = begin_section


    WRITE(string_default,*) this%nb_sides
    CALL compare_string(record, list, argument_data%nb_sides, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%nb_sides  
    END IF

    ALLOCATE(this%list_sides(this%nb_sides))
    WRITE(string_default,*) this%list_sides
    CALL compare_string(record, list, argument_data%list_sides, string_default, okay, i_list, j)
    IF (okay) THEN
       READ(list(i_list),*) this%list_sides
    ELSE
       list(i_list) = "0" 
    END IF

    i_list = i_list+1
    list(i_list) = end_section

    !===Closing unit
    CALL rewrite_data_from_list_record(rank, list, record, i_list, record_size)!!, .TRUE.)

    DEALLOCATE(begin_section)
    DEALLOCATE(end_section)
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
