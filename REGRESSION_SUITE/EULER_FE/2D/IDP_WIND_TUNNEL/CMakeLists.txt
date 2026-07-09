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
#=== Wind Tunnel Viscous ===#
add_regression_test(
  TEST_DIR "VISCOUS_DEFAULT"
  TEST_NAME "FE_2D_WIND_TUNNEL_VISCOUS"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;ERK-11;P1;hyperbolic;euler;viscous;lumped;dirichlet;udotn"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)
