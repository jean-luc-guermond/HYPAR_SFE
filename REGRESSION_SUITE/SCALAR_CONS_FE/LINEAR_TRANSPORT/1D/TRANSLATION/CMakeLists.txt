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
#=== FE_1D_LINEAR_TRANSPORT_TRANSLATION_GALERKIN ===#
add_regression_test(
  TEST_DIR "PERIODIC_GALERKIN_QUASI_CONSISTENT"
  TEST_NAME "FE_1D_LINEAR_TRANSPORT_GALERKIN_QUASI_CONSISTENT"
  FE_DIM 1D
  NB_REGEX 2
  TEST_LABELS "quick"
  CASE_LABELS "1D;P1;ERK-31;hyperbolic;linear_transport;galerkin;quasi_consistent;consistent;periodic"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== FE_1D_LINEAR_TRANSPORT_TRANSLATION_HIGH_ORDER_MINMOD ===#
add_regression_test(
  TEST_DIR "HIGH_ORDER_LIMITING_MINMOD"
  TEST_NAME "FE_1D_LINEAR_TRANSPORT_HIGH_ORDER_LIMITING_MINMOD"
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "1D;P1;ERK-31;hyperbolic;linear_transport;high;quasi_consistent;dirichlet;limiting;relax_minmod"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)
