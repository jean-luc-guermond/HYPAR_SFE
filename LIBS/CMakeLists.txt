MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Project name and version
project(hypar_lib)

# Compile C routines and create library
enable_language(C)

# Enable Fortran
enable_language(Fortran)

# Getting useful libraries
find_package(PkgConfig)

pkg_check_modules(PETSC PETSc)
pkg_check_modules(FFTW fftw3)

#=============================
# # Build library hypar_lib_F_1D
# #=============================

# file(GLOB_RECURSE hyper_SFE_sources_1D ${HYPAR_SFE_DIR}/COMMON_SOURCES/*.F90
# ${HYPAR_SFE_DIR}/PROBLEM_SOURCES/EULER/*.F90)
# list(REMOVE_ITEM hyper_SFE_sources_1D ${HYPAR_SFE_DIR}/COMMON_SOURCES/MESH/space_dim_*.F90)
# list(APPEND hyper_SFE_sources_1D ${HYPAR_SFE_DIR}/COMMON_SOURCES/MESH/space_dim_1D.F90)


# add_library(hypar_lib_F_1D SHARED ${hyper_SFE_sources_1D})
# target_include_directories(hypar_lib_F_1D PUBLIC
#     ${FFTW_INCLUDEDIR} ${FFTW_INCLUDE_DIRS}
#     ${PETSC_INCLUDEDIR} ${PETSC_INCLUDE_DIRS}
#   )
# target_link_libraries(hypar_lib_F_1D) 

# #=============================
# # Build library hypar_lib_F_2D
# #=============================

# file(GLOB_RECURSE hyper_SFE_sources_2D ${HYPAR_SFE_DIR}/COMMON_SOURCES/*.F90
# ${HYPAR_SFE_DIR}/PROBLEM_SOURCES/EULER/*.F90)
# list(REMOVE_ITEM hyper_SFE_sources_2D ${HYPAR_SFE_DIR}/COMMON_SOURCES/MESH/space_dim_*.F90)
# list(APPEND hyper_SFE_sources_2D ${HYPAR_SFE_DIR}/COMMON_SOURCES/MESH/space_dim_2D.F90)


# add_library(hypar_lib_F_2D SHARED ${hyper_SFE_sources_2D})
# target_include_directories(hypar_lib_F_2D PUBLIC
#     ${FFTW_INCLUDEDIR} ${FFTW_INCLUDE_DIRS}
#     ${PETSC_INCLUDEDIR} ${PETSC_INCLUDE_DIRS}
#   )
# target_link_libraries(hypar_lib_F_2D) 



file(GLOB_RECURSE hypar_SFE_sources ${HYPAR_SFE_DIR}/COMMON_SOURCES/*.F90
${HYPAR_SFE_DIR}/PROBLEM_SOURCES/EULER/*.F90)


add_library(hypar_lib_F SHARED ${hypar_SFE_sources})
target_include_directories(hypar_lib_F PUBLIC
    ${FFTW_INCLUDEDIR} ${FFTW_INCLUDE_DIRS}
    ${PETSC_INCLUDEDIR} ${PETSC_INCLUDE_DIRS}
  )
target_link_libraries(hypar_lib_F) 