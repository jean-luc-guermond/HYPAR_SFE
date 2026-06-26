# Project name and version
project(CTEST_NAVIER_STOKES_COMPRESSIBLE VERSION 1.0)
message(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

# List of tests
#=== ELASTICITY ===#
generate_random_procs_list(LIST_PROC 
MIN 4 MAX ${MAX_PROC_PARALLEL_TEST} DEFAULT 1
LEVEL ${PARALLEL_TEST_LEVEL} LIST_LEVEL (1 4 7 14))
add_regression_test(
  TEST_NAME "ELASTICITY"  NB_REGEX 2
  INCLUDE_TEMPLATE OFF FE_DIM 2D LIST_N_PROCS ${LIST_PROC}
)
