#=========================================#
#==== This is layer III of CTEST tree ====#
#=========================================#

# Project name and version
get_filename_component(
    CURRENT_CASE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    NAME
)
project(CTEST_${CURRENT_CASE} VERSION 2.0)
MESSAGE(STATUS "======> ${CURRENT_CASE}")

# List of tests
#=== 2D_P1_REFINEMENT_FACTOR_2 ===#
add_regression_test(
  TEST_DIR "P1_REFINEMENT_FACTOR_2"
  TEST_NAME "2D_P1_REFINEMENT_FACTOR_2"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;refinement;poisson;dirichlet;GMRES;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== 2D_P1_REFINEMENT_FACTOR_4 ===#
add_regression_test(
  TEST_DIR "P1_REFINEMENT_FACTOR_4"
  TEST_NAME "2D_P1_REFINEMENT_FACTOR_4"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;refinement;poisson;dirichlet;GMRES;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)


#=== 2D_P2_REFINEMENT_FACTOR_2 ===#
add_regression_test(
  TEST_DIR "P2_REFINEMENT_FACTOR_2"
  TEST_NAME "2D_P2_REFINEMENT_FACTOR_2"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;P2;refinement;poisson;dirichlet;GMRES;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== 2D_SCAL_ADV_PERIODIC ===#
add_regression_test(
  TEST_DIR "PERIODIC"
  TEST_NAME "2D_SCAL_ADV_PERIODIC"
  FE_DIM 2D
  NB_REGEX 3
  TEST_LABELS "quick"
  CASE_LABELS "2D;P2;refinement;poisson;dirichlet;periodic;GMRES;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)