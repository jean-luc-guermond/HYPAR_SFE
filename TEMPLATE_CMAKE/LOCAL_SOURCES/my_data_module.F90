   TYPE argument_limiting_type
      CHARACTER(LEN=rec_length) :: if_limiting       = '=== Apply cell-limiting (T/F) ? ==='
      CHARACTER(LEN=rec_length) :: limit_max         = '=== How many limiting iterations? ==='
      CHARACTER(LEN=rec_length) :: if_relax_bounds   = '=== Apply bound relaxation for limiting (T/F) ? ==='
      CHARACTER(LEN=rec_length) :: relaxation_method = '=== Relaxation method (avg/minmod) ==='
   END TYPE argument_limiting_type

   TYPE limiting_type
      CHARACTER(100) :: name
      LOGICAL                   :: if_limiting       = .False.
      INTEGER                   :: limit_max         = 2
      LOGICAL                   :: if_relax_bounds   = .False.
      CHARACTER(len=rec_length) :: relaxation_method ='minmod'
   CONTAINS
      PROCEDURE, PUBLIC  :: init => init_limiting
      PROCEDURE, PUBLIC  :: read => read_limiting_data
   END TYPE limiting_type

   SUBROUTINE read_limiting_data(this, section_name)
      USE read_inputs_module
      IMPLICIT NONE
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: section_name
      CLASS(limiting_type),    INTENT(INOUT) :: this
      TYPE(argument_limiting_type)           :: argument_data
      !================
      !=== MANDATORY Reading all data file
      !================
      IF (PRESENT(section_name)) THEN
         CALL read_data_init_list(section_name)
      ELSE
         CALL read_data_init_list()
      END IF

      !================
      !=== We now find the relevant information for this specific limiting data
      !================
      !===if_limiting
      CALL read_data(argument_data%if_limiting , this%if_limiting, &
      opt_name=this%name)

      !===Number of limiting iterations
      CALL read_data(argument_data%limit_max , this%limit_max, &
      opt_name=this%name, opt_add=this%if_limiting)

      !===if_relax_bounds
      CALL read_data(argument_data%if_relax_bounds, this%if_relax_bounds, &
      opt_name=this%name, opt_add=this%if_limiting)

      !===relaxation_method
      CALL read_data(argument_data%relaxation_method, this%relaxation_method, &
      opt_name=this%name, opt_add=(this%if_limiting .AND. this%if_relax_bounds))

      !================
      !=== MANDATORY to close data for the current section and
      !=== rewrite it with new information for the next sections
      !================
      CALL finalize_rewrite_data
   END SUBROUTINE read_limiting_data
