MODULE start_setup_MODULE
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE eos
   USE euler_type_module
   MPI_Comm        :: communicator
   TYPE setup_data_type
      LOGICAL :: if_regression_test
      LOGICAL :: if_restart
      REAL(KIND = 8) :: checkpointing_freq
      REAL(KIND = 8) :: final_time
      INTEGER :: syst_size
   CONTAINS
      PROCEDURE, PUBLIC :: init
   END TYPE setup_data_type

   TYPE(mesh_type), PUBLIC :: mesh
   TYPE(petsc_csr_LA), PRIVATE :: LA
   TYPE(euler_type), PUBLIC :: euler
   TYPE(setup_data_type), PUBLIC :: setup_data
   PUBLIC :: start_setup
   PRIVATE

CONTAINS

   SUBROUTINE start_setup
      use def_type_periodic
      USE construct_mesh
      USE st_matrix
      USE prep_periodic_module
      USE setup
      IMPLICIT NONE
      PetscErrorCode :: ierr
      TYPE(periodic_type) :: per
      REAL(KIND = 8) :: init_time = 0.d0
      CHARACTER(100) :: name = 'Euler 1'
      INTEGER :: ier, rank

      !===Start PETSC and MPI (mandatory)
      CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      communicator = PETSC_COMM_WORLD
      CALL MPI_Comm_rank(communicator, rank, ierr)
      !===Construct mesh
      CALL get_mesh(communicator, mesh)
      CALL prep_periodic(mesh, per)
      !===Construct LA
      CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA, opt_per = per)
      !===Read
      CALL read_setup_data(rank)

      !===Start Euler
      !FIXE ME init_time too
      CALL euler%init(communicator, name, mesh, LA, per, pressure, impose_bc_euler, init_time)

      !===Read data setup
   END SUBROUTINE start_setup

   SUBROUTINE init(this)
      CLASS(setup_data_type), INTENT(INOUT) :: this
      !===Logicals
      this%if_regression_test = .FALSE.
      this%if_restart = .FALSE.
      !===Reals
      this%checkpointing_freq = 1.d20
      this%final_time = 0.d0
      !===Characters
      !===Integers
   END SUBROUTINE init

   SUBROUTINE read_setup_data(rank)
      USE character_strings
      IMPLICIT NONE
      INTEGER, PARAMETER :: in_unit = 21
      CHARACTER(LEN = 100) :: argument
      INTEGER :: rank
      LOGICAL :: okay
      !===Initialize data to zero and false by default
      CALL setup_data%init

      OPEN(UNIT = in_unit, FILE = "data", FORM = 'formatted', STATUS = 'unknown')

      !===Restart
      argument = '===Restart (true/false)==='
      CALL find_string(in_unit, argument, okay)
      IF (okay) THEN
         READ (in_unit, *) setup_data%if_restart
      ELSE
         CALL default_data(rank, in_unit, argument, '.F.')
      END IF

      !===Checkpointing
      argument = '===Checkpointing frequency==='
      CALL find_string(in_unit, argument, okay)
      IF (okay) THEN
         READ (in_unit, *) setup_data%checkpointing_freq
      ELSE
         CALL default_data(rank, in_unit, argument, '1.d20')
      END IF

      !===Checkpointing
      argument = '===Final time==='
      CALL find_string(in_unit, argument, okay)
      IF (okay) THEN
         READ (in_unit, *) setup_data%final_time
      ELSE
         CALL default_data(rank, in_unit, argument, '0.1d0')
      END IF

      !===Regression test
      CALL getarg(1, argument)
      IF (trim(adjustl(argument))=='regression') THEN
         setup_data%if_regression_test = .true.
      ELSE
         setup_data%if_regression_test = .false.
      END IF

      CLOSE(in_unit)
   END SUBROUTINE read_setup_data

   SUBROUTINE default_data(rank, in_unit, argument, char)
      IMPLICIT NONE
      INTEGER :: rank, in_unit
      CHARACTER(*), INTENT(IN) :: argument, char
      IF (rank==0) WRITE(*, *)  TRIM(ADJUSTL(argument)) // ' not found. Set to ' // TRIM(ADJUSTL(char)) // &
           ' by default and added to data file.'
      CLOSE(in_unit)
      OPEN(UNIT = in_unit, FILE = "data", FORM = 'formatted', STATUS = 'unknown', position = 'append')
      IF (rank==0) THEN
         WRITE(in_unit, '(g0)') TRIM(ADJUSTL(argument))
         WRITE(in_unit, '(g0)') TRIM(ADJUSTL(char))
      END IF
   END SUBROUTINE default_data


END MODULE start_setup_MODULE
