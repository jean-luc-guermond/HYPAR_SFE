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
#=== 1D_LAGRANGE ===# test can only run on 1 proc
add_regression_test(
  TEST_DIR "1D_LAGRANGE"
  TEST_NAME "1D_LAGRANGE" 
  INCLUDE_TEST_SOURCES ON
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "tools"
  CASE_LABELS "1D;Pk"
  DEFAULT_PROC 1
)

#=== LEGENDRE ===# test can only run on 1 proc
add_regression_test(
  TEST_DIR "LEGENDRE"
  TEST_NAME "LEGENDRE"
  INCLUDE_TEST_SOURCES ON 
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "tools"
  CASE_LABELS "1D;Pk"
  DEFAULT_PROC 1
)
