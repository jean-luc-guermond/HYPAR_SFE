MODULE create_laplace_solver_ksp_module
#include "petsc/finclude/petsc.h"
    USE petsc
    USE def_type_mesh,         ONLY: mesh_type
    USE petsc_csr_LA_module,   ONLY: petsc_csr_LA

    USE dirichlet_type_module, ONLY: dirichlet_bc
    USE solver_data_module,    ONLY: solver_data_type
    USE read_inputs_module,    ONLY: rec_length

    TYPE argument_abstract_laplace_solver_type
        CHARACTER(LEN=rec_length) :: viscosity         = '=== Viscosity for parabolic ==='
        CHARACTER(LEN=rec_length) :: mass              = '=== mass for parabolic ==='
    END TYPE argument_abstract_laplace_solver_type

    TYPE, ABSTRACT :: abstract_laplace_solver_type
        REAL(KIND=8)   :: viscosity = 1.d0
        REAL(KIND=8)   :: mass      = 1.d0
        CHARACTER(LEN=:), ALLOCATABLE :: name
        TYPE(mesh_type),    POINTER   :: mesh
        TYPE(petsc_csr_LA), POINTER   :: LA
        TYPE(dirichlet_bc)            :: dir
        TYPE(solver_data_type)        :: solver_data

        Vec :: vec, vec_ghost, rhs
        Mat :: laplace_operator
        KSP :: my_ksp
        MPI_Comm, POINTER :: communicator
    CONTAINS
        PROCEDURE, PUBLIC  :: init     => init_laplace
        PROCEDURE, PRIVATE :: read     => read_laplace
        PROCEDURE, PRIVATE :: init_lhs => init_laplace_lhs
        PROCEDURE, PRIVATE :: init_vectors_rhs => init_vectors_rhs
        PROCEDURE, PUBLIC  :: solve    => solve_laplace
        PROCEDURE(template_dir_bc), DEFERRED :: dir_bc
    END TYPE abstract_laplace_solver_type

   ABSTRACT INTERFACE
        FUNCTION template_dir_bc(this, rr) RESULT(uu)
            USE def_type_mesh
            IMPORT :: abstract_laplace_solver_type
            IMPLICIT NONE
            CLASS(abstract_laplace_solver_type),    INTENT(INOUT) :: this
            REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
            REAL(KIND=8), DIMENSION(SIZE(rr,2))      :: uu
        END FUNCTION template_dir_bc
    END INTERFACE

CONTAINS

!========================================
!======= INIT PARAMETERS ================
!========================================
    SUBROUTINE init_laplace(this, mesh, LA, communicator, opt_name)
        IMPLICIT NONE
        CLASS(abstract_laplace_solver_type), INTENT(INOUT) :: this
        CHARACTER(LEN=100), OPTIONAL :: opt_name
        TYPE(mesh_type),    TARGET :: mesh
        TYPE(petsc_csr_LA), TARGET :: LA
        MPI_Comm,           TARGET :: communicator

        CALL this%read(opt_name)
        CALL this%init_lhs(mesh, LA, communicator, opt_name=opt_name)
        CALL this%init_vectors_rhs

    END SUBROUTINE init_laplace

    SUBROUTINE read_laplace(this, opt_name)
        USE space_dim
        USE read_inputs_module
        IMPLICIT NONE

        CHARACTER(LEN=rec_length)              :: raw_section_name='LAPLACE PARAMETERS '
        CHARACTER(LEN=:), ALLOCATABLE          :: section_name

        CLASS(abstract_laplace_solver_type), INTENT(INOUT) :: this
        TYPE(argument_abstract_laplace_solver_type)        :: argument_data
        CHARACTER(100), OPTIONAL                  :: opt_name

        !=== Reading all data file
        IF (PRESENT(opt_name)) THEN
            this%name = TRIM(ADJUSTL(opt_name)) 
            section_name = TRIM(ADJUSTL(raw_section_name)) // ' ' // this%name
        ELSE
            section_name = TRIM(ADJUSTL(raw_section_name))
        END IF
        CALL read_data_init_list(section_name)

        !================
        !=== We now find the relevant information for laplace
        !================

        !=== viscosity
        CALL read_data(argument_data%viscosity, this%viscosity, opt_name=this%name)
        !=== mass
        CALL read_data(argument_data%mass, this%mass, opt_name=this%name)

        !================
        !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
        !================
        CALL finalize_rewrite_data

    END SUBROUTINE read_laplace

!====================================
!======= INIT SOLVER ================
!====================================

    SUBROUTINE init_laplace_lhs(this, mesh, LA, communicator, opt_name)
        USE fem_M,            ONLY: qs_mass_diff_M
        USE solver_petsc,     ONLY: init_solver, create_local_petsc_matrix
        USE dir_nodes_petsc,  ONLY: Dirichlet_M_parallel
        USE compute_periodic, ONLY: periodic_matrix_petsc
        USE my_util,          ONLY: pack_opt
        IMPLICIT NONE
        CLASS(abstract_laplace_solver_type), INTENT(INOUT) :: this
        CHARACTER(LEN=100), OPTIONAL :: opt_name
        TYPE(mesh_type),    TARGET :: mesh
        TYPE(petsc_csr_LA), TARGET :: LA
        MPI_Comm,           TARGET :: communicator

    
        CALL pack_opt(this%name, "Laplace", opt_name)
        
        this%mesh => mesh
        this%LA   => LA
        CALL this%dir%set(mesh, this%name)
        this%communicator => communicator
        !===create petsc matrix ============================================
        CALL create_local_petsc_matrix(this%communicator, this%LA, this%laplace_operator, clean = .FALSE.)
        CALL qs_mass_diff_M (this%mesh, this%mass, this%viscosity, this%LA, this%laplace_operator)
        CALL periodic_matrix_petsc(this%mesh%per, this%LA, this%laplace_operator)
        CALL Dirichlet_M_parallel(this%laplace_operator, this%LA%loc_to_glob(1,this%dir%jsd))

        !===Create ksp solver
        CALL this%solver_data%init(opt_name)
        CALL init_solver(this%communicator, this%solver_data, this%my_ksp, this%laplace_operator)

        IF (this%mesh%rank == 0) THEN
            WRITE(*,*) "finished initializing solver for Laplace equation " // this%name
        END IF
    END SUBROUTINE init_laplace_lhs

    SUBROUTINE init_vectors_rhs(this)
        USE st_matrix, ONLY: create_my_ghost

        IMPLICIT NONE
        CLASS(abstract_laplace_solver_type), INTENT(INOUT) :: this
        INTEGER, POINTER, DIMENSION(:)            :: ifrom
        INTEGER :: ierr

        CALL create_my_ghost(this%mesh, this%LA, ifrom)
        CALL VecCreateGhost(this%communicator, this%mesh%dom_np, PETSC_DETERMINE, SIZE(ifrom), ifrom, this%vec, ierr)
        CALL VecDuplicate(this%vec, this%rhs, ierr)
        CALL VecGhostGetLocalForm(this%vec, this%vec_ghost, ierr)

    END SUBROUTINE init_vectors_rhs

!====================================
!======= APPLY SOLVER ===============
!====================================

    SUBROUTINE solve_laplace(this, u_in, u_out)
        USE fem_rhs,          ONLY: qs_00
        USE compute_periodic, ONLY: periodic_rhs_petsc
        USE dir_nodes_petsc,  ONLY: dirichlet_rhs
        USE st_matrix,        ONLY: extract_through_ghost
        USE solver_petsc,     ONLY: solver
        IMPLICIT NONE
        CLASS(abstract_laplace_solver_type),            INTENT(INOUT) :: this
        REAL(KIND=8), DIMENSION(:),            INTENT(IN)    :: u_in
        REAL(KIND=8), DIMENSION(SIZE(u_in,1)), INTENT(OUT)   :: u_out

        !=== Building rhs vector
        CALL qs_00 (this%mesh, this%LA, u_in, this%rhs)
        CALL periodic_rhs_petsc(this%mesh%per%nb_bords, this%mesh%per%list, this%mesh%per%perlist, this%rhs, this%LA)
        CALL dirichlet_rhs(this%LA%loc_to_glob(1, this%dir%jsd) - 1, this%dir_bc(this%mesh%rr(:, this%dir%jsd)), this%rhs)

        !=== Solving the linear system
        CALL solver(this%my_ksp, this%rhs, this%vec, reinit = .FALSE., verbose = this%solver_data%if_verbose)

        !=== Extract from Petsc vector to Fortran array
        CALL extract_through_ghost(this%vec, 1, 1, this%LA, u_out, opt_assemble=.FALSE.)

    END SUBROUTINE solve_laplace

END MODULE create_laplace_solver_ksp_module