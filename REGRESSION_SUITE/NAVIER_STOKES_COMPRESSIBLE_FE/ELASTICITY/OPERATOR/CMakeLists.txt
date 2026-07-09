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
#=== 2D_ELASTICITY_OPERATOR ===#
add_regression_test(
  TEST_DIR "DEFAULT"
  TEST_NAME "ELASTICITY_OPERATOR"
  FE_DIM 2D
  NB_REGEX 2
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;stokes;dirichlet;GMRES;CG;MUMPS;HYPRE"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)