MODULE euler_type_MODULE
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE def_type_periodic
   USE euler_bc_arrays
   USE Butcher_tableau
   USE euler_matrices_module
   USE space_dim
   IMPLICIT NONE

   ABSTRACT INTERFACE
      FUNCTION function_template_pressure(un) RESULT(vv)
         REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
         REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: vv
      END FUNCTION function_template_pressure
   END INTERFACE

   ABSTRACT INTERFACE
      SUBROUTINE function_template_impose_bc(un, euler_bc, mesh, time)
         USE petsc
         USE def_type_mesh
         REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
         TYPE(mesh_type) :: mesh
         TYPE(euler_bc_type) :: euler_bc
         REAL(KIND = 8) :: time
      END SUBROUTINE function_template_impose_bc
   END INTERFACE

   TYPE euler_type
      MPI_Comm :: communicator
      TYPE(mesh_type), POINTER :: mesh
      TYPE(petsc_csr_LA), POINTER :: LA
      TYPE(periodic_type), POINTER :: per
      PROCEDURE(function_template_pressure), NOPASS, POINTER :: pressure
      PROCEDURE(function_template_impose_bc), NOPASS, POINTER :: pressure
      TYPE(BT), PUBLIC :: ERK
      TYPE(euler_bc_type) :: euler_bc
      TYPE(euler_matrices_type) :: matrices
      REAL(KIND = 8) :: dt, time
      INTEGER :: syst_dim = k_dim + 2
   CONTAINS
      PROCEDURE, PUBLIC :: init
      PROCEDURE, PUBLIC :: update
   END TYPE euler_type

CONTAINS
   SUBROUTINE init(a, communicator, mesh, LA, per, pressure, erk_s, impose_bc, time_init)
      CLASS(euler_type), INTENT(INOUT) :: a
      MPI_Comm, INTENT(IN) :: communicator
      TYPE(mesh_type), TARGET, INTENT(IN) :: mesh
      TYPE(petsc_csr_LA), TARGET, INTENT(IN) :: LA
      TYPE(periodic_type), TARGET, INTENT(IN) :: per
      INTEGER :: erk_sv
      REAL(KIND = 8) :: time_init

      INTERFACE
         FUNCTION pressure(un) RESULT(vv)
            REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
            REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: vv
         END FUNCTION pressure
      END INTERFACE
      a%mesh => mesh
      a%communicator = communicator
      a%LA => LA
      a%per => per
      a%pressure => pressure
      a%impose_bc => impose_bc
      a%euler_bc%syst_size = a%syst_size
      a%time = time_init
      CALL a%ERK%init(erk_sv)
      CALL a%euler_bc%construct_euler_bc(a%mesh)
      CALL a%matrices%construct(a%communicator, a%mesh, a%LA)
   END SUBROUTINE init

   SUBROUTINE update(this, un)
      USE petsc_tools
      USE euler_flux
      CLASS(euler_type) :: this
      REAL(KIND = 8), DIMENSION(this%mesh%np, this%syst_dim), INTENT(INOUT) :: un
      REAL(KIND = 8), DIMENSION(this%mesh%np, k_dim) :: ff
      REAL(KIND = 8), DIMENSION(this%mesh%np) :: rk
      INTEGER k

      DO comp = 1, thid%syst_size
         ff = flux(comp, un)

         CALL VecSet(x2vec, 0.d0, ierr)
         DO k = 1, k_dim
            !=== set flux_k in x1vec
            CALL array_to_petsc_vec(ff(:, k), x1vec, this%mesh, this%LA, 'insert')
            !=== compute sum_j cij_k * fluxj_k in x3vec
            CALL MatMult(this%matrices%cij(k), x1vec, x3vec, ierr)
            !=== construct sum_k sum_j cij_k flux_k into x2vec
            CALL VecAXPY(x2vec, 1.d0, x3vec, ierr)
         END DO

         CALL VecGhostGetLocalForm(x2vec, x2_ghost, ierr)
         CALL VecGhostUpdateBegin(x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(x2vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL extract(x2_ghost, 1, 1, LA, rk)

         rk = rk * this%dt / this%matrices%lumped_mass

         this%time = this%time + this%dt

         un(:, comp) = un(:, comp) + rk

         CALL this%impose_bc(un, this%euler_bc, this%mesh, this%time)
      END DO

   END SUBROUTINE update


END MODULE euler_type_MODULE
