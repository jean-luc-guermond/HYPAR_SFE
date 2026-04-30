# Project name and version
project(CTEST_FFT_1D VERSION 1.0)
message(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

#List of tests
add_regression_test(
  TEST_NAME "FFT_1D" NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 1D
)

add_regression_test(
  TEST_NAME "FFT_1D_PLAN" NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 1D
)
