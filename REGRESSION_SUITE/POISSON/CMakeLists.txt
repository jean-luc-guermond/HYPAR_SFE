#========================================#
#==== This is layer II of CTEST tree ====#
#========================================#

# Project name and version
get_filename_component(
    CURRENT_CASE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    NAME
)
project(CTEST_${CURRENT_CASE} VERSION 2.0)
print_case_name(TEST_NAME CTEST_${CURRENT_CASE})


# Generate CASES_USER from CASES_DEVELOPPER
create_cases_user(
  SOURCES_DIR POISSON/1D
  FE_DIM 1
)
add_subdirectory(1D/DEFAULT)

# Generate CASES_USER from CASES_DEVELOPPER
create_cases_user(
  SOURCES_DIR POISSON/2D
  FE_DIM 2
)
add_subdirectory(2D/DEFAULT)