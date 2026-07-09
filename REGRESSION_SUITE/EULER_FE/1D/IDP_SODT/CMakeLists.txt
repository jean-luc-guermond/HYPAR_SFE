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
#=== SODT Viscous ===#
add_regression_test(
  TEST_DIR "VISCOUS_DEFAULT"
  TEST_NAME "FE_1D_SODT_VISCOUS" 
  NB_REGEX 1
  FE_DIM 1D
  TEST_LABELS "quick"
  CASE_LABELS "1D;ERK11;P1;hyperbolic;euler;viscous;lumped;dirichlet"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== SODT High order ===#
add_regression_test(
  TEST_DIR "HIGH_ORDER_DEFAULT"
  TEST_NAME "FE_1D_SODT_HIGH_ORDER"
  NB_REGEX 1
  FE_DIM 1D
  TEST_LABELS "quick"
  CASE_LABELS "1D;ERK-31;P1;hyperbolic;euler;high;quasi_consistent;dirichlet"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== SODT High order with limiting ===#
add_regression_test(
  TEST_DIR "HIGH_ORDER_LIMITING"
  TEST_NAME "FE_1D_SODT_HIGH_ORDER_LIMITING"
  NB_REGEX 1
  FE_DIM 1D
  TEST_LABELS "quick"
  CASE_LABELS "1D;ERK31;P1;hyperbolic;euler;high;quasi_consistent;dirichlet;limiting"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== SODT High order with limiting and minmod ===#
add_regression_test(
  TEST_DIR "HIGH_ORDER_LIMITING_MINMOD"
  TEST_NAME "FE_1D_SODT_HIGH_ORDER_LIMITING_MINMOD"
  NB_REGEX 1
  FE_DIM 1D
  TEST_LABELS "quick"
  CASE_LABELS "1D;ERK31;P1;hyperbolic;euler;high;quasi_consistent;dirichlet;limiting;relax_minmod"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)
