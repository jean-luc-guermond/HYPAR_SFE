# Project name and version
project(CTEST_PETSC VERSION 1.0)
message(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

# List of tests
#=== PETSC_MatGet_MatSet ===#
set(LIST_PROC 1)
add_regression_test(
  TEST_NAME "PETSC_MatGet_MatSet"  NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 1D LIST_N_PROCS ${LIST_PROC}
)

#=== PETSC_petscSF ===#
set(LIST_PROC 4)
add_regression_test(
  TEST_NAME "PETSC_petscSF"  NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 2D LIST_N_PROCS ${LIST_PROC}
)

#=== PETSC_MatUpdateMPIAIJWitjArray ===#
set(LIST_PROC 4)
add_regression_test(
  TEST_NAME "PETSC_MatUpdateMPIAIJWitjArray"  NB_REGEX 1
  INCLUDE_TEMPLATE OFF FE_DIM 2D LIST_N_PROCS ${LIST_PROC}
)