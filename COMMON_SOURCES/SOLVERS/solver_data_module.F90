MODULE solver_data_module
   USE read_inputs_module, ONLY: rec_length
   IMPLICIT NONE
  !===chain of characters that should appear in data file

   TYPE argument_solver_data_type
      CHARACTER(len = rec_length) :: it_max          = '=== Maximum number of iterations ==='
      CHARACTER(len = rec_length) :: rel_tol         = '=== Relative tolerance ==='
      CHARACTER(len = rec_length) :: abs_tol         = '=== Absolute tolerance ==='
      CHARACTER(len = rec_length) :: solver          = '=== Which solver? (CG,GMRES,...) ==='
      CHARACTER(len = rec_length) :: precond         = '=== Which precond? (MUMPS,HYPRE,...) ==='
      CHARACTER(len = rec_length) :: boomeramg_strong_threshold = '=== BoomerAMG (HYPRE): Strong threshold ==='
      CHARACTER(len = rec_length) :: boomeramg_coarsen_type     = '=== BoomerAMG : coarsening type (Falgout,...) ==='
      CHARACTER(len = rec_length) :: boomeramg_relax_type_all   = '=== BoomerAMG : relax type (Chebyshev,...) ==='
      CHARACTER(len = rec_length) :: if_fixed_v_cycle           = '=== BoomerAMG : fixed v cycle? ==='
      CHARACTER(len = rec_length) :: number_v_cycle             = '=== BoomerAMG : numer of v cycle (if fixed) ==='
      CHARACTER(len = rec_length) :: if_verbose      = '=== verbose solver? (True/False) ==='
      CHARACTER(len = rec_length) :: if_residual     = '=== monitor residuals? (True/False) ==='
   END TYPE argument_solver_data_type

   !===default value in simulation
   TYPE solver_data_type
   !=== parameters from data file
      INTEGER         :: it_max          = 5000
      REAL(KIND=8)    :: rel_tol         = 1.d-10
      REAL(KIND=8)    :: abs_tol         = 1.d-18
      LOGICAL                        :: if_verbose  = .False.
      LOGICAL                        :: if_residual = .False.
      CHARACTER(len = rec_length)    :: solver = "GMRES"
      CHARACTER(len = rec_length)    :: precond = "MUMPS"
      CHARACTER(LEN = rec_length)    :: boomeramg_strong_threshold = '0.1' 
      CHARACTER(LEN = rec_length)    :: boomeramg_coarsen_type     = 'Falgout'  
      CHARACTER(LEN = rec_length)    :: boomeramg_relax_type_all   = 'Chebyshev'    
      CHARACTER(LEN = rec_length)    :: number_v_cycle             = '2'     
      LOGICAL                        :: if_fixed_v_cycle           = .FALSE.  
      CHARACTER(100)                 :: name = ''
   !=== parameters throughout resolutions
      INTEGER                       :: count
      REAL(KIND = 8)                :: tps_solver_ref
      REAL(KIND = 8)                :: tps_ratio      = 1000.d0

   CONTAINS
      PROCEDURE, PUBLIC              :: read => read_solver_data
      PROCEDURE, PUBLIC              :: init => init_solver_data
   END TYPE solver_data_type
CONTAINS

   SUBROUTINE init_solver_data(this, opt_name)
      IMPLICIT NONE
      CLASS(solver_data_type), INTENT(INOUT) :: this
      CHARACTER(LEN=*), OPTIONAL               :: opt_name
      CALL this%read(opt_name)
   END SUBROUTINE init_solver_data

   SUBROUTINE read_solver_data(this, opt_name)
      USE space_dim
      USE read_inputs_module
      IMPLICIT NONE

      CHARACTER(LEN=rec_length)              :: raw_section_name='SOLVER PARAMETERS '
      CHARACTER(LEN=:), ALLOCATABLE          :: section_name

      CLASS(solver_data_type), INTENT(INOUT) :: this
      TYPE(argument_solver_data_type)        :: argument_data
      CHARACTER(LEN=*), OPTIONAL               :: opt_name
      LOGICAL :: if_hypre

      !=== Reading all data file
      IF (PRESENT(opt_name)) THEN
         this%name = TRIM(ADJUSTL(opt_name)) 
         section_name = TRIM(ADJUSTL(raw_section_name)) // ' ' // this%name
      ELSE
         section_name = TRIM(ADJUSTL(raw_section_name))
      END IF
      CALL read_data_init_list(section_name)

      !================
      !=== We now find the relevant information for the mesh
      !================

      !=== it_max
      CALL read_data(argument_data%it_max, this%it_max, opt_name=this%name)
      !=== rel_tol
      CALL read_data(argument_data%rel_tol, this%rel_tol, opt_name=this%name)
      !=== abs_tol
      CALL read_data(argument_data%abs_tol, this%abs_tol, opt_name=this%name)
      !=== solver
      CALL read_data(argument_data%solver, this%solver, opt_name=this%name)
      !=== precond
      CALL read_data(argument_data%precond, this%precond, opt_name=this%name)

      if_hypre=TRIM(ADJUSTL(this%precond))=='HYPRE'
      !=== strong_thr
      CALL read_data(argument_data%boomeramg_strong_threshold, this%boomeramg_strong_threshold, opt_name=this%name,opt_add=if_hypre)

      !=== boomeramg_coarsen_type
      CALL read_data(argument_data%boomeramg_coarsen_type, this%boomeramg_coarsen_type, opt_name=this%name,opt_add=if_hypre)

      !=== boomeramg_relax_type_all
      CALL read_data(argument_data%boomeramg_relax_type_all, this%boomeramg_relax_type_all, opt_name=this%name,opt_add=if_hypre)

      !=== if_fixed_v_cycle
      CALL read_data(argument_data%if_fixed_v_cycle, this%if_fixed_v_cycle, opt_name=this%name,opt_add=if_hypre)

      !=== number_v_cycle
      CALL read_data(argument_data%number_v_cycle, this%number_v_cycle, opt_name=this%name, opt_add=this%if_fixed_v_cycle.AND.if_hypre)

      !=== if_verbose      
      CALL read_data(argument_data%if_verbose, this%if_verbose, opt_name=this%name)

      !=== if_residual      
      CALL read_data(argument_data%if_residual, this%if_residual, opt_name=this%name)
      
      !================
      !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
      !================
      CALL finalize_rewrite_data

   END SUBROUTINE read_solver_data

END MODULE solver_data_module