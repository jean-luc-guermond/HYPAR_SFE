MODULE navier_stokes_module
!>> limited global uses to avoid unexpected behaviors
#include "petsc/finclude/petsc.h"
   USE petsc
   USE petsc_tools,                          ONLY: array_to_petsc_vec
   USE euler_type_module,                    ONLY: euler_type
   USE stokes_parabolic_module,              ONLY: stokes_parabolic_type
   USE cell_limiting_engine_parallel_module, ONLY: limiting_type, limiting_functional_type, limiting_all_functional_type
   USE def_type_mesh,                        ONLY: mesh_type
   USE read_inputs_module,                   ONLY: rec_length
   USE space_dim,                            ONLY: k_dim
   USE Butcher_tableau
   USE Implicit_Butcher_tableau
!>> limited global uses to avoid unexpected behaviors

 IMPLICIT NONE
   TYPE argument_navier_stokes_type
      CHARACTER(LEN=rec_length) :: erk_sv               = '=== ERK ? ==='
      CHARACTER(LEN=rec_length) :: mu_viscosity         = '=== Value of mu viscosity ==='
      CHARACTER(LEN=rec_length) :: lambda_viscosity     = '=== Value of lambda viscosity ==='
      CHARACTER(LEN=rec_length) :: thermal_diffusivity  = '=== Value of thermal diffusivity ==='
      CHARACTER(LEN=rec_length) :: cv                   = '=== Value of thermal capacity at constant volume ==='
   END TYPE argument_navier_stokes_type

   TYPE :: navier_stokes_type
        REAL(KIND = 8)                :: thermal_diffusivity = 0.d0
        REAL(KIND = 8)                :: mu_viscosity = 0.d0
        REAL(KIND = 8)                :: lambda_viscosity =0.d0
        REAL(KIND = 8)                :: cv = 0.d0
        !===Parameters built along way
        MPI_Comm :: communicator
        Vec          :: x1vec, x2vec, x2_ghost, vec_loc, x3vec
        Vec          :: x4vec, x5vec !!!!!Conveniance vectors to be used only inside procedures!!!!
        CHARACTER(LEN=:), ALLOCATABLE :: name
        INTEGER                       :: syst_dim = k_dim + 2
        REAL(KIND = 8) :: dt, time, final_time
        TYPE(mesh_type),     POINTER :: mesh
        INTEGER                      :: imex_sv
        TYPE(BT)                     :: ERK
        TYPE(IBT)                    :: IRK
        TYPE(euler_type)             :: euler
        TYPE(stokes_parabolic_type)  :: stokes
    CONTAINS
        PROCEDURE, PUBLIC   :: init => init_navier_stokes
        PROCEDURE, PUBLIC   :: update
        PROCEDURE, PUBLIC   :: set_times
        PROCEDURE, PRIVATE  :: read => read_navier_stokes_data!, init_vectors
   END TYPE navier_stokes_type

CONTAINS

SUBROUTINE init_navier_stokes(this, communicator, name, mesh, limiting_functionals_euler) 
    USE space_dim
    USE my_util,            ONLY: error_petsc, to_str
    USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
    USE construct_mesh, ONLY: generate_boundary_structure, refine_mesh_1D2D_Pk2P1, create_gauss_points_1D_2D
    IMPLICIT NONE
    CLASS(navier_stokes_type), INTENT(INOUT) :: this
    MPI_Comm,                   INTENT(IN) :: communicator
    CHARACTER(LEN=*),             INTENT(IN) :: name
    TYPE(mesh_type), TARGET,    INTENT(IN) :: mesh
    TYPE(limiting_functional_type), DIMENSION(:), TARGET :: limiting_functionals_euler

    this%name = name
    this%mesh => mesh
    this%communicator = communicator

    CALL this%read("NAVIER STOKES PARAMETERS FOR "//trim(adjustl(this%name)))

    !=== Define Navier-Stokes IMEX tableaux
    CALL this%ERK%init(this%imex_sv)
    CALL this%IRK%init(this%imex_sv)
    !=== Define Navier-Stokes IMEX tableaux
    
    !===Start Euler & ERK
    CALL this%euler%init_euler(TRIM(ADJUSTL(name)) // '/Euler')
    this%euler%erk_sv = this%imex_sv
    CALL this%euler%init_hyperbolic(communicator, TRIM(ADJUSTL(name)) // '/Euler', mesh, limiting_functionals_euler)

    !=== Start Stokes & IRK
    this%stokes%thermal_diffusivity = this%thermal_diffusivity
    this%stokes%mu_viscosity        = this%mu_viscosity
    this%stokes%lambda_viscosity    = this%lambda_viscosity
    this%stokes%cv                  = this%cv
    this%stokes%irk_sv              = this%imex_sv
    CALL this%stokes%init(communicator, TRIM(ADJUSTL(name)) // '/Stokes', mesh)

END SUBROUTINE init_navier_stokes

SUBROUTINE read_navier_stokes_data(this, section_name)
    USE read_inputs_module
    USE my_util, ONLY: get_tab_idx_char
    IMPLICIT NONE
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: section_name
    CLASS(navier_stokes_type), INTENT(INOUT)  :: this
    TYPE(argument_navier_stokes_type)         :: argument_data

    !================
    !=== MANDATORY Reading all data file
    !================
    IF (PRESENT(section_name)) THEN
        CALL read_data_init_list(section_name)
    ELSE
        CALL read_data_init_list()
    END IF

    !================
    !=== We now find the relevant information for this specific Euler data
    !================
    !===mu
    CALL read_data(argument_data%mu_viscosity, this%mu_viscosity, opt_name=this%name)

    !===lambda
    CALL read_data(argument_data%lambda_viscosity, this%lambda_viscosity, opt_name=this%name)

    !===lambda
    CALL read_data(argument_data%thermal_diffusivity, this%thermal_diffusivity, opt_name=this%name)
    
    !===cv
    CALL read_data(argument_data%cv, this%cv, opt_name=this%name)

    !================
    !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
    !================
    CALL finalize_rewrite_data
END SUBROUTINE read_navier_stokes_data

SUBROUTINE set_times(this, times)
    IMPLICIT NONE
    CLASS(navier_stokes_type)              :: this
    REAL(KIND=8), DIMENSION(2), INTENT(IN) :: times
    this%time       = times(1) !<==initial_time
    this%final_time = times(2) !<==final_time
    CALL this%euler%set_times(times)
    CALL this%stokes%set_times(times)
END SUBROUTINE set_times

SUBROUTINE update(this,un_in)
    IMPLICIT NONE
    CLASS(navier_stokes_type)                                          :: this
    REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim)               :: un_in
    REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s+1) :: urk
    REAL(KIND=8), DIMENSION(this%mesh%np, this%syst_dim, this%ERK%s)   :: euler_flux_rk_at_dof

    INTEGER  :: stage
    urk(:,:,1) = un_in
    DO stage = 2, this%ERK%s+1
        CALL this%euler%one_step_ERK(stage,urk,euler_flux_rk_at_dof)
        this%dt = this%euler%dt
        this%stokes%dt = this%euler%dt
        CALL this%stokes%one_step_IRK(stage,urk)
    END DO
    un_in = urk(:,:,this%ERK%s+1)
    this%time = this%time + this%dt
    this%euler%time = this%euler%time + this%dt
    this%stokes%time = this%stokes%time + this%dt
END SUBROUTINE update
END MODULE navier_stokes_module