# Project name and version
project(CTEST_NL_scalar_cons_fourier VERSION 1.0)
message(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

# List of tests
add_regression_test(
  TEST_NAME "FOURIER_1D_BURGERS"  NB_REGEX 2
  INCLUDE_TEMPLATE ON FE_DIM 1D
)

add_regression_test(
  TEST_NAME "FOURIER_1D_LINEAR_TRANSPORT"  NB_REGEX 2
  INCLUDE_TEMPLATE ON FE_DIM 1D
)

add_regression_test(
  TEST_NAME "FOURIER_1D_NON_CONVEX"  NB_REGEX 1
  INCLUDE_TEMPLATE ON FE_DIM 1D
)
