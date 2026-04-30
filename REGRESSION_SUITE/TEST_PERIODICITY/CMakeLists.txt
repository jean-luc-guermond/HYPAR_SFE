# Project name and version
project(CTEST_PERIODIC VERSION 1.0)
message(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

# List of tests
add_regression_test(
  TEST_NAME "PERIODIC_FE_1D_EULER_SMOOTH"  NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 1D
)

add_regression_test(
  TEST_NAME "PERIODIC_1D_SCAL_ADV_DIFF"  NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 1D
)

add_regression_test(
  TEST_NAME "PERIODIC_2D_SCAL_ADV_DIFF"  NB_REGEX 3
  INCLUDE_TEMPLATE OFF FE_DIM 2D
)
