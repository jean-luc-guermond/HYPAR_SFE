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
#=== LINEAR TRANSPORT Default ===#
add_regression_test(
  TEST_DIR "DEFAULT"
  TEST_NAME "FFT_LINEAR_TRANSPORT"
  FE_DIM 1D
  NB_REGEX 2
  TEST_LABELS "quick"
  CASE_LABELS "Fourier;ERK-41;hyperbolic;linear_transport;high;limiting;relax_avg"
  DEFAULT_PROC 1
)
