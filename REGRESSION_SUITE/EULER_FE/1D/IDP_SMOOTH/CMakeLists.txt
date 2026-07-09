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
#=== FE_1D_SMOOTH ===#
add_regression_test(
  TEST_DIR "VISCOUS_DEFAULT"
  TEST_NAME "FE_1D_SMOOTH_VISCOUS"
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "1D;ERK11;P1;hyperbolic;euler;viscous;lumped"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== PERIODIC_FE_1D_EULER_SMOOTH ===#
add_regression_test(
  TEST_DIR "VISCOUS_PERIODIC"
  TEST_NAME "PERIODIC_FE_1D_EULER_SMOOTH" 
  FE_DIM 1D 
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "1D;ERK11;P1;hyperbolic;euler;viscous;lumped;periodic"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)