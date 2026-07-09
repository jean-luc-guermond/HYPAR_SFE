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
#=== Isentropic vortex Viscous ===#
add_regression_test(
  TEST_DIR "VISCOUS_DEFAULT"
  TEST_NAME "FE_2D_ISENTROPIC_VORTEX_VISCOUS"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;ERK11;P1;hyperbolic;euler;viscous;lumped"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== Isentropic vortex high order ===#
add_regression_test(
  TEST_DIR "HIGH_ORDER_DEFAULT"
  TEST_NAME "FE_2D_ISENTROPIC_VORTEX_HIGH_ORDER"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;ERK31;P1;hyperbolic;euler;high;quasi_consistent"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== Isentropic vortex high order with limiting ===#
add_regression_test(
  TEST_DIR "HIGH_ORDER_LIMITING"
  TEST_NAME "FE_2D_ISENTROPIC_VORTEX_HIGH_ORDER_LIMITING"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;ERK31;P1;hyperbolic;euler;high;quasi_consistent;limiting"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== Isentropic vortex high order with limiting and relax AVG ===#
add_regression_test(
  TEST_DIR "HIGH_ORDER_LIMITING_AVG"
  TEST_NAME "FE_2D_ISENTROPIC_VORTEX_HIGH_ORDER_LIMITING_AVG"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;ERK31;P1;hyperbolic;euler;high;quasi_consistent;limiting;relax_avg"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== Isentropic vortex high order P2 ===#
add_regression_test(
  TEST_DIR "HIGH_ORDER_P2"
  TEST_NAME "FE_2D_ISENTROPIC_VORTEX_HIGH_ORDER_P2"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;ERK31;P2;hyperbolic;euler;high;consistent"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== Isentropic vortex high order P3 ===#
add_regression_test(
  TEST_DIR "HIGH_ORDER_P3"
  TEST_NAME "FE_2D_ISENTROPIC_VORTEX_HIGH_ORDER_P3"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;ERK31;P3;hyperbolic;euler;high;quasi_consistent"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)