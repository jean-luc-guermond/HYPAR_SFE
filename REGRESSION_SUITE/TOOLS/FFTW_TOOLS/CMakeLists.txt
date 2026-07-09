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

#List of tests
#=== FFT_1D ===# test can only run on 1 proc
add_regression_test(
  TEST_DIR "FFT_1D"
  TEST_NAME "FFT_1D"
  INCLUDE_TEST_SOURCES ON
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "tools"
  CASE_LABELS "fftw"
  DEFAULT_PROC 1
)

#=== FFT_1D_PLAN ===# test can only run on 1 proc
add_regression_test(
  TEST_DIR "FFT_1D_PLAN"
  TEST_NAME "FFT_1D_PLAN"
  INCLUDE_TEST_SOURCES ON 
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "tools"
  CASE_LABELS "fftw"
  DEFAULT_PROC 1
)
