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
  USE read_inputs_module
  USE post_processing_debug_MODULE
  IMPLICIT NONE
  TYPE(mesh_type)     :: mesh
  TYPE(petsc_csr_LA)  :: LA
  TYPE(dirichlet_bc)  :: dir
  TYPE(solver_param)  :: my_par
  !TYPE(periodic_type), DIMENSION(1) :: per
  !INTEGER, POINTER, DIMENSION(:) :: js_d_loc
  INTEGER, POINTER, DIMENSION(:) :: ifrom
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: un, tab_norm
  REAL(KIND=8) :: error, norm
  INTEGER :: rank
  CHARACTER(5) :: char
! for regression test
  CHARACTER(100) :: string
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
call write_l1_petsc_mat(mass, mesh, LA, 'qs_mass_diff_M')
  CALL periodic_matrix_petsc(mesh%per, LA, mass)
call write_l1_petsc_mat(mass, mesh, LA, 'periodic mat')
  CALL Dirichlet_M_parallel(mass, LA%loc_to_glob(1,dir%jsd))
call write_l1_petsc_mat(mass, mesh, LA, 'dirichlet mat')


  CALL create_my_ghost(mesh, LA, ifrom) !=== creating ghost structures
  CALL VecCreateGhost(PETSC_COMM_WORLD, mesh%dom_np, PETSC_DETERMINE, SIZE(ifrom), ifrom, sol, ierr)
  CALL VecDuplicate(sol, rhs, ierr)
  CALL VecGhostGetLocalForm(sol, ghost_sol, ierr) !=== creating pointer ghost_sol => sol used after solving lin syst

  !===set PBC + rhs ==================================================
  CALL qs_00 (mesh, LA, source(mesh%rr), rhs) !=== create rhs with a scalar source term
call write_l1_petsc_vec(rhs, ghost_sol, mesh, LA, 'qs_00')
  CALL periodic_rhs_petsc(mesh%per%nb_bords, mesh%per%list, mesh%per%perlist, rhs, LA) !=== give periodic structure to rhs
call write_l1_petsc_vec(rhs, ghost_sol, mesh, LA, 'periodic')
  ! write(*,*) mesh%per%nb_bords, mesh%per%list(1)%DIL, mesh%per%perlist(1)%DIL
  !=== setting rhs: LA%... - 1 ==> global indexing starts at 0); ex_sol... => associated values
  CALL dirichlet_rhs(LA%loc_to_glob(1, dir%jsd) - 1, ex_sol(mesh%rr(:, dir%jsd)), rhs)
call write_l1_petsc_vec(rhs, ghost_sol, mesh, LA, 'dirichlet')

  !===solving linear system ==========================================
  CALL init_solver(my_par, my_ksp, mass, PETSC_COMM_WORLD, solver = 'MUMPS', precond = 'MUMPS')

  CALL solver(my_ksp, rhs, sol, reinit = .FALSE., verbose = .FALSE.)
  CALL VecGhostUpdateBegin(sol, INSERT_VALUES, SCATTER_FORWARD, ierr)
  CALL VecGhostUpdateEnd(sol, INSERT_VALUES, SCATTER_FORWARD, ierr)
  ALLOCATE(un(mesh%np))
  CALL extract(ghost_sol, 1, 1, LA, un) !=== extracting solution from pointer ghost_sol

  !===post processing =================================================
  WRITE(char, '(I5)') mesh%rank
  CALL plot_scalar_field(mesh%jj, mesh%rr, un, 'SOL' // trim(adjustl(char)) // '.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un-ex_sol(mesh%rr), 'error' // trim(adjustl(char)) // '.plt')

  CALL array_to_petsc_vec(ex_sol(mesh%rr), rhs, mesh, LA, 'insert') !=== HYPAR subroutine
  CALL VecAssemblyBegin(rhs, ierr)
  CALL VecAssemblyEnd(rhs, ierr)
  CALL VecNorm(rhs, NORM_1, norm, ierr)
  CALL VecAXPY(rhs, -1.d0, sol, ierr)
  CALL VecNorm(rhs, NORM_1, error, ierr)
  ! IF (rank==0) WRITE(*, '(A,g12.9)') 'L1 NORM error ', error / norm
  IF (rank==0) WRITE(*, *) 'L1 NORM error ', error / norm

  ! !===regression test =================================================
  ! CALL getarg(1, string)
  ! IF (trim(adjustl(string))=='regression') THEN
  !   IF (error / norm < 5.d-4) THEN
  !     IF (rank == 0) WRITE(*, '(A,A)') 'Regression test passed', '1234567891'
  !   ELSE
  !     IF (rank == 0) WRITE(*, '(A)') 'Regression test failed'
  !   END IF
  ! END IF

  !===regression test =================================================
  CALL getarg(1, string)
  IF (trim(adjustl(string))=='regression') THEN
       ALLOCATE(tab_norm(1))
       tab_norm(1) = error / norm
      !  CALL get_num_test(num_test)
       CALL regression(tab_norm)!, opt_num_test=num_test)
  END IF


  CALL error_petsc('End of test')
  CALL PetscFinalize(ierr)

CONTAINS

  FUNCTION source(rr) RESULT(uu)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:, :) :: rr
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: uu
    REAL(KIND = 8) :: kmax=4*ACOS(-1.d0)
    uu = (1+kmax**2)*COS(kmax * rr(1, :) +  .7d0)
  END FUNCTION source
  FUNCTION ex_sol(rr) RESULT(uu)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:, :) :: rr
    REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: uu
    REAL(KIND = 8) :: kmax=4*ACOS(-1.d0)
    uu = COS(kmax * rr(1, :) +  .7d0 )
  END FUNCTION ex_sol


    SUBROUTINE write_l1(field_in, mesh, in_char)
      USE fem_tn, ONLY : ns_l1

      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: field_in
      REAL(KIND=8) :: norm_loc, norm
      CHARACTER(LEN=*), INTENT(IN) :: in_char
      INTEGER :: ierr

      CALL ns_l1(mesh, field_in(:), norm_loc)
      CALL MPI_ALLREDUCE(norm_loc,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,mesh%comm,ierr)
      IF(mesh%rank==0) WRITE(*, *) in_char, '=>  L1 Norm = ', norm

   END SUBROUTINE write_l1

   SUBROUTINE write_l1_petsc_vec(vec_in, vec_ghost, mesh, LA, in_char)

      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      TYPE(petsc_csr_la) :: la
      REAL(KIND=8), DIMENSION(mesh%np) :: field_in
      REAL(KIND=8) :: norm_loc, norm
      CHARACTER(LEN=*), INTENT(IN) :: in_char
      INTEGER :: ierr
      Vec :: vec_in, vec_ghost

      CALL extract_through_ghost(vec_in, vec_ghost, 1, 1, LA, field_in, 'insert',&
      opt_assemble=.FALSE.)

      CALL write_l1(field_in, mesh, in_char)

   END SUBROUTINE write_l1_petsc_vec

   SUBROUTINE write_l1_petsc_mat(mat_in, mesh, LA, in_char)

      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      TYPE(petsc_csr_la) :: la
      REAL(KIND=8), DIMENSION(mesh%np) :: field_in
      REAL(KIND=8) :: norm_loc, norm
      CHARACTER(LEN=*), INTENT(IN) :: in_char
      INTEGER :: ierr
      Mat :: Mat_in

      CALL MatNorm(Mat_in, NORM_FROBENIUS, norm, ierr)

      IF (mesh%rank==0) write(*,*) in_char, norm

   END SUBROUTINE write_l1_petsc_mat

END PROGRAM test_matrix
