PROGRAM test_matrix
#include "petsc/finclude/petsc.h"
   USE construct_mesh
   USE input_data
   USE def_type_mesh
   USE petsc
   USE solver_petsc
   USE fem_M
   IMPLICIT NONE
   TYPE(mesh_type) :: mesh
   TYPE(petsc_csr_LA) :: LA
   Mat :: mass
   Vec :: test_vec, test_vec2
   INTEGER, POINTER, DIMENSION(:) :: js_d_loc
   MPI_Comm       :: communicator
   PetscErrorCode :: ierr
   INTEGER :: rank
   !===Start PETSC and MPI (mandatory)=============================================
   CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   CALL MPI_Comm_rank(communicator, rank, ierr)
   !CALL create_cart_comm(k_dim, comm_cart, comm_one_d, coord_cart)

   !===User reads his/her own data=================================================
   CALL read_my_data('data')
   write(*, *) 'ok1'
   CALL get_mesh(PETSC_COMM_WORLD, mesh, LA, js_d_loc, 1)
   write(*, *) 'ok2'
   write(*, *) rank, mesh%disp, mesh%discell
   write(*, *) rank, mesh%jj
   write(*, *) rank, mesh%jj_extra
   write(*, *) rank, mesh%loc_to_glob

   CALL create_local_petsc_matrix(PETSC_COMM_WORLD, LA, mass, clean = .FALSE.)
   write(*, *) 'ok3'

   CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, mass)

   CALL dirichlet_rhs(LA%loc_to_glob(1, :) - 1, source(mesh%rr), test_vec)

   CALL MatMult(mass, test_vec, test_vec2, ierr)

   write(*, *) 'ok4'

CONTAINS

   FUNCTION source(rr) RESULT(uu)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: uu
      REAL(KIND = 8) :: pi
      pi = ACOS(-1.d0)
      uu = COS(16 * rr(1, :)) * COS(16 * rr(2, :))
   END FUNCTION source

END PROGRAM test_matrix