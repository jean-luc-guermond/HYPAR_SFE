# Project name and version
project(CTEST_EULER_FE_2D VERSION 1.0)
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

# List of tests
add_regression_test(
  TEST_NAME "FE_2D_ISENTROPIC_VORTEX"  NB_REGEX 1
  INCLUDE_TEMPLATE ON FE_DIM 2D
)

add_regression_test(
  TEST_NAME "FE_2D_WIND_TUNNEL"  NB_REGEX 1
  INCLUDE_TEMPLATE ON FE_DIM 2D
)
