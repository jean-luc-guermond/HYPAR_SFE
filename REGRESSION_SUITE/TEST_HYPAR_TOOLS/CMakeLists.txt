# Project name and version
project(CTEST_PERIODIC VERSION 1.0)
message(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

# List of tests
#=== REFINEMENT_P1 ===#
generate_random_procs_list(LIST_PROC 
MIN 1 MAX 16 DEFAULT 4
LEVEL ${PARALLEL_TEST_LEVEL} LIST_LEVEL (1 3 7 13))
add_regression_test(
  TEST_NAME "REFINEMENT_P1"  NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 2D LIST_N_PROCS ${LIST_PROC}
)
