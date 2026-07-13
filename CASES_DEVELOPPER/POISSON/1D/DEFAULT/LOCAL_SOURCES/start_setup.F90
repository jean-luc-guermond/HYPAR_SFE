MODULE start_setup_MODULE

#include "petsc/finclude/petsc.h"
    USE petsc
    USE def_type_mesh
    USE read_inputs_module
    USE setup
    USE petsc_csr_LA_module, ONLY: petsc_csr_LA

    TYPE argument_setup_data_type
        CHARACTER(LEN=rec_length) :: if_analytical_ref  = '=== Do we compare with analytical reference? (true/false) ==='
    END TYPE argument_setup_data_type

    TYPE setup_data_type
        LOGICAL        :: if_analytical_ref   = .FALSE.
        INTEGER        :: syst_size
    CONTAINS
        PROCEDURE, PUBLIC :: read => read_setup_data
        PROCEDURE, PUBLIC :: init => init_setup_data
    END TYPE setup_data_type

    TYPE(mesh_type),                   PUBLIC :: mesh
    TYPE(petsc_csr_LA),               PRIVATE :: LA
    TYPE(setup_data_type),             PUBLIC :: setup_data
    TYPE(periodic_type), DIMENSION(1), PUBLIC :: per
    TYPE(laplace_solver_type),         PUBLIC :: Laplace

    MPI_Comm :: communicator
    PUBLIC :: start_setup
    PRIVATE

CONTAINS

  SUBROUTINE start_setup
    ! use periodic_data_module
    USE construct_mesh,     ONLY: get_mesh
    USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
    USE setup
    USE options_module
    IMPLICIT NONE
    CHARACTER(100) :: name
    INTEGER        :: rank
    PetscErrorCode :: ierr

    !===Start PETSC and MPI (mandatory)

    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    communicator = PETSC_COMM_WORLD
    CALL MPI_Comm_rank(communicator, rank, ierr)

    !===Read executable arguments
    CALL read_all_arguments

    !===Clean data once
    CALL clean_data_once

    !===Construct mesh
    CALL get_mesh(communicator, mesh)

    !===Construct LA
    CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA)
    
    !===Read
    CALL setup_data%init

    !===Start Laplace
    name = 'Laplace 1'
    CALL Laplace%init(mesh, LA, communicator, opt_name=name)
    
    setup_data%syst_size = 1

  END SUBROUTINE start_setup

  SUBROUTINE init_setup_data(this)
    CLASS(setup_data_type), INTENT(INOUT) :: this
    CALL this%read
  END SUBROUTINE init_setup_data

  SUBROUTINE read_setup_data(this)
    IMPLICIT NONE

    CHARACTER(LEN=rec_length) :: section_name='SETUP PARAMETERS'
    CLASS(setup_data_type)             :: this
    TYPE(argument_setup_data_type)     :: argument_data
    CHARACTER(LEN=rec_length)          :: string

    !================
    !=== MANDATORY Reading all data file
    !================
    CALL read_data_init_list(section_name)

    !================
    !=== We now find the relevant information for this setup
    !================

    !===Analytical reference
    CALL read_data(argument_data%if_analytical_ref, this%if_analytical_ref)

    !================
    !=== MANDATORY to close data for the current section and rewrite it with new information for the next sections
    !================
    CALL finalize_rewrite_data

  END SUBROUTINE read_setup_data


END MODULE start_setup_MODULE
