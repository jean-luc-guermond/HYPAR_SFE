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
#=== FE_2D_LINEAR_TRANSPORT_TRANSLATION_GALERKIN ===#
add_regression_test(
  TEST_DIR "DEFAULT"
  TEST_NAME "FE_2D_LINEAR_TRANSPORT_GALERKIN_QUASI_CONSISTENT_ROTATION"
  FE_DIM 2D
  NB_REGEX 2
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;ERK-31;hyperbolic;linear_transport;galerkin;lumped;quasi_consistent"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)
