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

# PETSc tools
add_subdirectory(PETSC_TOOLS)

# fftw3 tools
add_subdirectory(FFTW_TOOLS)

# hypar_sfe FE tools
add_subdirectory(FE_TOOLS)