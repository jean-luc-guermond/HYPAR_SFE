#=======================================#
#==== This is layer I of CTEST tree ====#
#=======================================#

# Project name and version
get_filename_component(
    CURRENT_CASE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    NAME
)
project(CTEST_${CURRENT_CASE} VERSION 2.0)

# POISSON
add_subdirectory(POISSON)

# EULER_FE
add_subdirectory(EULER_FE)

# SCALAR_CONS_FE
add_subdirectory(SCALAR_CONS_FE)

# SCALAR_CONS_FFT
add_subdirectory(SCALAR_CONS_FFT)

# NAVIER_STOKES_COMPRESSIBLE_FE
add_subdirectory(NAVIER_STOKES_COMPRESSIBLE_FE)

# TOOLS
add_subdirectory(TOOLS)