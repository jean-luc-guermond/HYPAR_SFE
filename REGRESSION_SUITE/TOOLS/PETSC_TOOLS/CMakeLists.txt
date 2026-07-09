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
#=== PETSC_MatGet_MatSet ===#
add_regression_test(
  TEST_DIR "PETSC_MatGet_MatSet"
  TEST_NAME "PETSC_MatGet_MatSet"
  INCLUDE_TEST_SOURCES ON
  FE_DIM 1D
  NB_REGEX 1
  TEST_LABELS "tools"
  CASE_LABELS "petsc"
  PROC_MIN 1 PROC_MAX 1 DEFAULT_PROC 1
)

#=== PETSC_petscSF ===#
add_regression_test(
  TEST_DIR "PETSC_petscSF"
  TEST_NAME "PETSC_petscSF"
  INCLUDE_TEST_SOURCES ON
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "tools"
  CASE_LABELS "petsc"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)

#=== PETSC_MatUpdateMPIAIJWithArray ===#
add_regression_test(
  TEST_DIR "PETSC_MatUpdateMPIAIJWithArray"
  TEST_NAME "PETSC_MatUpdateMPIAIJWithArray"
  INCLUDE_TEST_SOURCES ON 
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "tools"
  CASE_LABELS "petsc"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
)