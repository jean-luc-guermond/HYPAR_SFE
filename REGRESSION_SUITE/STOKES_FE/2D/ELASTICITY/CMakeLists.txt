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
#=== TEST_1 ===#
add_regression_test(
  TEST_DIR "TEST_1"
  TEST_NAME "FE_2D_STOKES_FORCING_1"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;IRK-21;stokes;forcing;elasticity;dirichlet;lumped_stokes;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 1
)

#=== TEST_2 ===#
add_regression_test(
  TEST_DIR "TEST_2"
  TEST_NAME "FE_2D_STOKES_FORCING_2"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;IRK-21;stokes;forcing;elasticity;dirichlet;lumped_stokes;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 1
)

#=== TEST_3 ===#
add_regression_test(
  TEST_DIR "TEST_3"
  TEST_NAME "FE_2D_STOKES_FORCING_3"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;IRK21;stokes;rho_variable;dirichlet;lumped_stokes;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 1
)

#=== TEST_4 ===#
add_regression_test(
  TEST_DIR "TEST_4"
  TEST_NAME "FE_2D_STOKES_FORCING_4"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;IRK31;stokes;forcing;elasticity;rho_variable;dirichlet;lumped_stokes;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 1
)

#=== TEST_5 ===#
add_regression_test(
  TEST_DIR "TEST_5"
  TEST_NAME "FE_2D_STOKES_FORCING_5"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;IRK-31;stokes;forcing_variable;dirichlet;lumped_stokes;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 1
)

#=== TEST_6 ===#
add_regression_test(
  TEST_DIR "TEST_6"
  TEST_NAME "FE_2D_STOKES_FORCING_6"
  FE_DIM 2D
  NB_REGEX 1
  TEST_LABELS "quick"
  CASE_LABELS "2D;P1;IRK-31;stokes;forcing_variable;elasticity;rho_variable;dirichlet;lumped_stokes;MUMPS"
  PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 1
)