MODULE start_setup_MODULE
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE eos
   USE euler_type_module
   MPI_Comm        :: communicator
   INTEGER :: rank
   TYPE(mesh_type) :: mesh
   TYPE(petsc_csr_LA) :: LA
   TYPE(euler_type) :: euler
CONTAINS
   SUBROUTINE start_setup
      use def_type_periodic
      USE construct_mesh
      USE st_matrix
      USE prep_periodic_module
      PetscErrorCode :: ierr
      TYPE(periodic_type) :: per
      !===Start PETSC and MPI (mandatory)
      communicator = PETSC_COMM_WORLD
      CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      CALL MPI_Comm_rank(communicator, rank, ierr)

      !===Construct mesh
      CALL get_mesh(communicator, mesh)
      CALL prep_periodic(mesh, per)

      !===Construct LA
      CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA, opt_per = per)

      !===Start Euler
      !FIXE ME ERK -21
      CALL euler%init(communicator, mesh, LA, per, pressure, -21, impose_bc_euler)

      !===Read data setup
   END SUBROUTINE start_setup
END MODULE start_setup_MODULE
