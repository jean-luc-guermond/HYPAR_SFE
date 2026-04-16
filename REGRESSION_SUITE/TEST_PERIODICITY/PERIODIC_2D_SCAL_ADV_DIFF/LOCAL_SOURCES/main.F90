PROGRAM test_matrix
#include "petsc/finclude/petsc.h"
  USE petsc
  USE petsc_tools
  USE construct_mesh
  USE def_type_mesh
  USE periodic_data_module
  USE dirichlet_type_module
  USE compute_periodic
  USE solver_petsc
  USE fem_M
  USE fem_rhs
  USE dir_nodes
  USE dir_nodes_petsc
  USE st_matrix
  USE sub_plot
  USE character_strings
  USE post_processing_debug_MODULE
  IMPLICIT NONE
  TYPE(mesh_type)     :: mesh
  TYPE(petsc_csr_LA)  :: LA
  TYPE(dirichlet_bc)  :: dir
  TYPE(solver_param)  :: my_par
  !TYPE(periodic_type), DIMENSION(1) :: per
  !INTEGER, POINTER, DIMENSION(:) :: js_d_loc
  INTEGER, POINTER, DIMENSION(:) :: ifrom
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: un
  REAL(KIND=8) :: error, norm
  INTEGER :: rank
  CHARACTER(5) :: char
! for regression test
  CHARACTER(100) :: string
  REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: tab_norm
! for regression test
  Mat :: mass
  Vec :: rhs, sol, ghost_sol
  KSP :: my_ksp
  MPI_Comm       :: communicator
  PetscErrorCode :: ierr

  !===Start PETSC and MPI (mandatory)=============================================

  CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  communicator = PETSC_COMM_WORLD
  CALL MPI_Comm_rank(communicator, rank, ierr)

  my_par%it_max = 5000
  my_par%rel_tol = 1.d-10
  my_par%abs_tol = 1.d-18
  my_par%verbose = .FALSE.
  my_par%solver = 'MUMPS'
  my_par%precond = 'MUMPS'

  !===read mesh data =================================================
  CALL clean_data_once
  CALL get_mesh(communicator, mesh)

  CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA)
  CALL dir%set(mesh, "a")

  !===create petsc matrix ============================================
  CALL create_local_petsc_matrix(PETSC_COMM_WORLD, LA, mass, clean = .FALSE.)
  CALL qs_mass_diff_M (mesh, 1.d0, 1.d0, LA, mass)
  CALL periodic_matrix_petsc(mesh%per, LA, mass)
  CALL Dirichlet_M_parallel(mass, LA%loc_to_glob(1,dir%jsd))


  CALL create_my_ghost(mesh, LA, ifrom)
  CALL VecCreateGhost(PETSC_COMM_WORLD, mesh%dom_np, PETSC_DETERMINE, SIZE(ifrom), ifrom, sol, ierr)
  CALL VecDuplicate(sol, rhs, ierr)
  CALL VecGhostGetLocalForm(sol, ghost_sol, ierr)

  !===set PBC + rhs ==================================================
  CALL qs_00 (mesh, LA, source(mesh%rr), rhs)
  CALL periodic_rhs_petsc(mesh%per%nb_bords, mesh%per%list, mesh%per%perlist, rhs, LA)
  write(*,*) mesh%per%nb_bords, mesh%per%list(1)%DIL, mesh%per%perlist(1)%DIL
  CALL dirichlet_rhs(LA%loc_to_glob(1, dir%jsd) - 1, ex_sol(mesh%rr(:, dir%jsd)), rhs)

  !===solving linear system ==========================================
  CALL init_solver(my_par, my_ksp, mass, PETSC_COMM_WORLD, solver = 'MUMPS', precond = 'MUMPS')

  CALL solver(my_ksp, rhs, sol, reinit = .FALSE., verbose = .FALSE.)
  CALL VecGhostUpdateBegin(sol, INSERT_VALUES, SCATTER_FORWARD, ierr)
  CALL VecGhostUpdateEnd(sol, INSERT_VALUES, SCATTER_FORWARD, ierr)
  ALLOCATE(un(mesh%np))
  CALL extract(ghost_sol, 1, 1, LA, un)

  !===post processing =================================================
  WRITE(char, '(I5)') mesh%rank
  CALL plot_scalar_field(mesh%jj, mesh%rr, un, 'SOL' // trim(adjustl(char)) // '.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un-ex_sol(mesh%rr), 'error' // trim(adjustl(char)) // '.plt')

  CALL array_to_petsc_vec(ex_sol(mesh%rr), rhs, mesh, LA, 'insert')
  CALL VecAssemblyBegin(rhs, ierr)
  CALL VecAssemblyEnd(rhs, ierr)
  CALL VecNorm(rhs, NORM_1, norm, ierr)
  CALL VecAXPY(rhs, -1.d0, sol, ierr)
  CALL VecNorm(rhs, NORM_1, error, ierr)
  if(rank==0) WRITE(*, '(A,g12.3)') 'L1 NORM error ', error / norm

  !===regression test =================================================
  CALL getarg(1, string)
  IF (trim(adjustl(string))=='regression') THEN
       ALLOCATE(tab_norm(1))
       tab_norm(1) = error / norm
       CALL regression(tab_norm)
  END IF

CONTAINS

  FUNCTION source(rr) RESULT(uu)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:, :) :: rr
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: uu
    REAL(KIND = 8) :: kxmax=ACOS(-1.d0), kymax=3*ACOS(-1.d0)
    uu = (1+kxmax**2+kymax**2)*SIN(kxmax * rr(1, :))*SIN(kymax * rr(2, :) + .1d0)
  END FUNCTION source
  FUNCTION ex_sol(rr) RESULT(uu)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:, :) :: rr
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: uu
    REAL(KIND = 8) :: kxmax=ACOS(-1.d0), kymax=3*ACOS(-1.d0)
    uu = SIN(kxmax * rr(1, :))*SIN(kymax * rr(2, :) + .1d0)
  END FUNCTION ex_sol
END PROGRAM test_matrix
