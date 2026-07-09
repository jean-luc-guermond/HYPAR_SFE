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
  SOURCES_DIR SCALAR_CONS_FE/LINEAR_TRANSPORT/1D
  FE_DIM 1
)

add_subdirectory(LINEAR_TRANSPORT/1D/TRANSLATION)

# Generate CASES_USER from CASES_DEVELOPPER
create_cases_user(
  SOURCES_DIR SCALAR_CONS_FE/LINEAR_TRANSPORT/2D
  FE_DIM 2
)

add_subdirectory(LINEAR_TRANSPORT/2D/PERIODIC_TRANSLATION)

add_subdirectory(LINEAR_TRANSPORT/2D/THREE_BODY_ROTATION)

# Generate CASES_USER from CASES_DEVELOPPER
create_cases_user(
  SOURCES_DIR SCALAR_CONS_FE/BURGERS/1D
  FE_DIM 1
)