# Project name and version
project(CTEST_FFT_1D VERSION 1.0)
message(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

#List of tests
#=== FFT_1D ===# test can only run on 1 proc
set(LIST_PROC 1)
add_regression_test(
  TEST_NAME "FFT_1D" NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 1D LIST_N_PROCS ${LIST_PROC}
)

#=== FFT_1D_PLAN ===# test can only run on 1 proc
set(LIST_PROC 1)
add_regression_test(
  TEST_NAME "FFT_1D_PLAN" NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 1D LIST_N_PROCS ${LIST_PROC}
)
