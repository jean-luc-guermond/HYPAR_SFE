MODULE euler_type_MODULE
#include "petsc/finclude/petsc.h"
  USE petsc
  USE def_type_mesh
  USE def_type_periodic
  USE euler_bc_arrays
  USE Butcher_tableau
  IMPLICIT NONE
  ABSTRACT INTERFACE
     FUNCTION function_template(un) RESULT(vv)
       REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
       REAL(KIND=8), DIMENSION(SIZE(un,1))            :: vv
     END FUNCTION function_template
   END INTERFACE
  
  TYPE euler_type
     MPI_Comm                    :: communicator
     TYPE(mesh_type), POINTER    :: mesh
     TYPE(petsc_csr_LA), POINTER :: LA
     TYPE(periodic_type), POINTER:: per
     PROCEDURE(function_template), nopass, POINTER :: pressure
     TYPE(BT), PUBLIC :: ERK
     TYPE(euler_bc_type) :: euler_bc
   CONTAINS
     PROCEDURE, PUBLIC :: init
     PROCEDURE, PRIVATE :: construct_euler_bc
  END TYPE euler_type
  
CONTAINS
  SUBROUTINE init(a, communicator, mesh, LA, per, pressure, erk_sv)
    CLASS(euler_type), INTENT(INOUT) :: a
    MPI_Comm, INTENT(IN) ::  communicator
    TYPE(mesh_type), TARGET, INTENT(IN) :: mesh
    TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
    TYPE(periodic_type), TARGET, INTENT(IN):: per
    INTEGER :: erk_sv
    INTERFACE
       FUNCTION pressure(un) RESULT(vv)
         REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
         REAL(KIND=8), DIMENSION(SIZE(un,1))            :: vv
       END FUNCTION pressure
    END INTERFACE
    a%mesh => mesh
    a%communicator = communicator
    a%LA => LA
    a%per => per
    a%pressure => pressure
    a%ERK%init(erk_sv)
    a%construct_euler_bc
    
  END SUBROUTINE init
END MODULE euler_type_MODULE
