# Project name and version
project(hypar_lib)
MESSAGE(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")
MESSAGE(STATUS "")

include(variables_lib.cmake)

#==========================================================================
#=============== compilation of HYPAR_SFE
#==========================================================================

# Getting mandatory libraries
find_package(MPI REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)

pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)
# Getting mandatory libraries

# Debugging information for libraries
message(STATUS "Working in: ${HYPAR_SFE_DIR}")
message(STATUS "ZLIB_LIBRARIES=${ZLIB_LIBRARIES}")
message(STATUS "MPI_Fortran_LIBRARIES=${MPI_Fortran_LIBRARIES}")
get_target_property(PETSC_INC PkgConfig::PETSC INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "PETSC includes: ${PETSC_INC}")
get_target_property(FFTW_INC PkgConfig::FFTW INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "FFTW includes: ${FFTW_INC}")
# Debugging information for libraries

# Building HYPAR_SFE library
file(GLOB_RECURSE hypar_SFE_sources CONFIGURE_DEPENDS ${HYPAR_SFE_DIR}/COMMON_SOURCES/*.F90
${HYPAR_SFE_DIR}/COMMON_SOURCES/POST_PROCESSING/*.c
${HYPAR_SFE_DIR}/PROBLEM_SOURCES/*.F90)

if(BUILD_ALL_DIMS)
    # Compiling HYPAR_SFE_FE_1D library
   set(DIMS 1 2)
else()
    # Compiling HYPAR_SFE_FE_2D library
   set(DIMS ${FE_dim})
endif()

foreach(dim IN LISTS DIMS)
    set(MY_HYPAR_SFE_LIB hypar_lib_F_${dim}D)
    set(MY_HYPAR_SFE_MODULE_SOURCES ${HYPAR_SFE_DIR}/LIBS/BUILD/${dim}D/mod)

    add_library(${MY_HYPAR_SFE_LIB} SHARED ${hypar_SFE_sources})
    target_link_libraries(${MY_HYPAR_SFE_LIB} PUBLIC
        PkgConfig::PETSC PkgConfig::FFTW
        ${MPI_Fortran_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${ADDITIONAL_LINKS})
    target_compile_definitions(${MY_HYPAR_SFE_LIB} PUBLIC K_DIM_FE=${dim})
    message(STATUS "HYPAR_SFE_FE_${dim} sources: ${MY_HYPAR_SFE_MODULE_SOURCES}")
    set_target_properties(${MY_HYPAR_SFE_LIB} PROPERTIES
        Fortran_MODULE_DIRECTORY ${MY_HYPAR_SFE_MODULE_SOURCES})
endforeach()
