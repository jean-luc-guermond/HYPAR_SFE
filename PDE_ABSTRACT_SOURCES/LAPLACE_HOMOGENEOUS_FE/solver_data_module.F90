MODULE solver_data_module
   USE read_inputs_module, ONLY: rec_length
   IMPLICIT NONE
  !===chain of characters that should appear in data file

   TYPE argument_solver_data_type
      CHARACTER(len = rec_length) :: it_max          = '=== Maximum number of iterations ==='
      CHARACTER(len = rec_length) :: rel_tol         = '=== Relative tolerance ==='
      CHARACTER(len = rec_length) :: abs_tol         = '=== Absolute tolerance ==='
      CHARACTER(len = rec_length) :: verbose         = '=== verbose solver? (True/False) ==='
      CHARACTER(len = rec_length) :: solver          = '=== Which solver? (MUMPS) ==='
      CHARACTER(len = rec_length) :: precond         = '=== Which precond? (MUMPS) ==='
      ! CHARACTER(len = rec_length) :: strong_thr         = '=== Strong threshold for HYPRE preconditioner ==='
   END TYPE argument_solver_data_type
   !===default value in simulation
   TYPE solver_data_type
      INTEGER         :: it_max          = 5000
      REAL(KIND=8)    :: rel_tol         = 1.d-10
      REAL(KIND=8)    :: abs_tol         = 1.d-18
      LOGICAL                        :: verbose = .False.
      CHARACTER(len = rec_length)    :: solver = "GMRES"
      CHARACTER(len = rec_length)    :: precond = "MUMPS"
      ! CHARACTER(len = rec_length)    :: strong_thr = '0.7'
      CHARACTER(100)                 :: name = ''
   CONTAINS
      PROCEDURE, PUBLIC              :: read => read_solver_data
      PROCEDURE, PUBLIC              :: init => init_solver_data
      PROCEDURE, PUBLIC              :: set  => set_my_par_solver
   END TYPE solver_data_type
CONTAINS

   SUBROUTINE init_solver_data(this, opt_name)
      IMPLICIT NONE
      CLASS(solver_data_type), INTENT(INOUT) :: this
      CHARACTER(100), OPTIONAL               :: opt_name
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
      CHARACTER(100), OPTIONAL               :: opt_name

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
      !=== verbose      
      CALL read_data(argument_data%verbose, this%verbose, opt_name=this%name)
      !=== solver
      CALL read_data(argument_data%solver, this%solver, opt_name=this%name)
      !=== precond
      CALL read_data(argument_data%precond, this%precond, opt_name=this%name)
      ! !=== strong_thr
      ! CALL read_data(argument_data%strong_thr, this%strong_thr)
      !================
      !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
      !================
      CALL finalize_rewrite_data

   END SUBROUTINE read_solver_data


   SUBROUTINE set_my_par_solver(this, my_par)
      USE solver_petsc
      CLASS(solver_data_type), INTENT(INOUT) :: this
      TYPE(solver_param),      INTENT(INOUT) :: my_par

      my_par%it_max  = this%it_max
      my_par%rel_tol = this%rel_tol
      my_par%abs_tol = this%abs_tol
      my_par%verbose = this%verbose
      my_par%solver  = TRIM(ADJUSTL(this%solver))
      my_par%precond = TRIM(ADJUSTL(this%precond))

   END SUBROUTINE set_my_par_solver

END MODULE solver_data_module