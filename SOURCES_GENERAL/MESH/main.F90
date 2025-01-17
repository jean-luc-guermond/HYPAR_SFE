PROGRAM test_matrix
#include "petsc/finclude/petsc.h"
   USE construct_mesh
   USE input_data
   USE def_type_mesh
   USE petsc
   USE solver_petsc
   USE fem_M
   USE dir_nodes
   USE dir_nodes_petsc
   USE st_matrix
   IMPLICIT NONE
   TYPE(mesh_type) :: mesh
   TYPE(petsc_csr_LA) :: LA
   Mat :: mass
   Vec :: test_vec, test_vec2, test_vec3
   KSP   :: my_ksp
   INTEGER, POINTER, DIMENSION(:) :: js_d_loc
   INTEGER, POINTER, DIMENSION(:) :: ifrom
   REAL(KIND = 8) :: error
   TYPE(solver_param), PRIVATE :: my_par
   MPI_Comm       :: communicator
   PetscErrorCode :: ierr
   INTEGER :: rank
   !===Start PETSC and MPI (mandatory)=============================================
   CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   CALL MPI_Comm_rank(communicator, rank, ierr)
   !CALL create_cart_comm(k_dim, comm_cart, comm_one_d, coord_cart)

   my_par%it_max = 5000
   my_par%rel_tol = 1.d-10
   my_par%abs_tol = 1.d-18
   my_par%verbose = .FALSE.
   my_par%solver = 'MUMPS'
   my_par%precond = 'MUMPS'

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

   CALL create_my_ghost(mesh, LA, ifrom)
   CALL VecCreateGhost(PETSC_COMM_WORLD, mesh%dom_np, PETSC_DETERMINE, SIZE(ifrom), ifrom, test_vec, ierr)
   CALL VecDuplicate(test_vec, test_vec2, ierr)
   CALL VecDuplicate(test_vec, test_vec3, ierr)

   CALL dirichlet_rhs(LA%loc_to_glob(1, :) - 1, source(mesh%rr), test_vec)

   CALL MatMult(mass, test_vec, test_vec2, ierr)

   write(*, *) 'ok4'

   CALL init_solver(my_par, my_ksp, mass, PETSC_COMM_WORLD, olver = 'MUMPS', precond = 'MUMPS')

   CALL solver(my_ksp, test_vec2, test_vec3, reinit = .FALSE., verbose = .FALSE., type = 'AAAA')

   write(*, *) 'ok5'

   CALL VecAXPY(test_vec, -1.d0, test_vec3, ierr)
   CALL VecNorm(test_vec, NORM_1, error, ierr)

   write(*,*) 'error', error

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