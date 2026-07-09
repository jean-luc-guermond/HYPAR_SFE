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
#=== 1D_P2_DIRICHLET ===#
add_regression_test(
  TEST_DIR "P2_DIRICHLET"
  TEST_NAME "1D_P2_DIRICHLET"
  FE_DIM 1D
  NB_REGEX 2
  TEST_LABELS "quick"
  CASE_LABELS "1D;P2;poisson;dirichlet;GMRES;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== 1D_P3_DIRICHLET ===#
add_regression_test(
  TEST_DIR "P3_DIRICHLET"
  TEST_NAME "1D_P3_DIRICHLET"
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "1D;P3;poisson;dirichlet;GMRES;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== 1D_P3_REFINEMENT_FACTOR_4 ===#
add_regression_test(
  TEST_DIR "P3_DIRICHLET_REFINEMENT_FACTOR_4"
  TEST_NAME "1D_P3_DIRICHLET_REFINEMENT_FACTOR_4"
  FE_DIM 1D
  NB_REGEX 2
  TEST_LABELS "quick"
  CASE_LABELS "1D;P3;refinement;poisson;periodic;dirichlet;GMRES;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== 1D_P8_DIRICHLET ===#
add_regression_test(
  TEST_DIR "P8_DIRICHLET"
  TEST_NAME "1D_P8_DIRICHLET"
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "1D;P8;refinement;poisson;periodic;dirichlet;GMRES;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== PERIODIC ===#
add_regression_test(
  TEST_DIR "PERIODIC"
  TEST_NAME "1D_SCAL_ADV_PERIODIC"
  FE_DIM 1D
  NB_REGEX 2
  TEST_LABELS "quick"
  CASE_LABELS "1D;P1;poisson;periodic;GMRES;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)