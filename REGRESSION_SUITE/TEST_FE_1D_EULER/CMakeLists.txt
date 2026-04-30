# Project name and version
project(CTEST_EULER_FE_1D VERSION 1.0)
message(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

# List of tests
add_regression_test(
  TEST_NAME "FE_1D_SODT" NB_REGEX 1
  INCLUDE_TEMPLATE ON FE_DIM 1D
)

add_regression_test(
  TEST_NAME "FE_1D_SMOOTH" NB_REGEX 1
  INCLUDE_TEMPLATE ON FE_DIM 1D
)

add_regression_test(
  TEST_NAME "FE_1D_VANDERWALLS" NB_REGEX 1
  INCLUDE_TEMPLATE ON FE_DIM 1D
)
