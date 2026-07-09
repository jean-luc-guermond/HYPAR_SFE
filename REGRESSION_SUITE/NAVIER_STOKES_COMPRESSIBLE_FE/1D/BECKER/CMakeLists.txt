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
#=== FE_1D_BECKER ===#
add_regression_test(
  TEST_DIR "DEFAULT"
  TEST_NAME "FE_1D_BECKER_HIGH_ORDER"
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "1D;P1;IMEX;ERK-31;IRK-31;navier-stokes;stokes;hyperbolic;euler;dirichlet;high;quasi_consistent;limiting;lumped_stokes;GMRES;CG;MUMPS;HYPRE"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)