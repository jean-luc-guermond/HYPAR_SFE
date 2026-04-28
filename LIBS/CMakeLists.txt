MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Cmake minimum version
cmake_minimum_required(VERSION 3.29)

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


file(GLOB_RECURSE hypar_SFE_sources ${HYPAR_SFE_DIR}/COMMON_SOURCES/*.F90
${HYPAR_SFE_DIR}/COMMON_SOURCES/POST_PROCESSING/*.c
${HYPAR_SFE_DIR}/PROBLEM_SOURCES/EULER/*.F90
${HYPAR_SFE_DIR}/PROBLEM_SOURCES/IDP_FOURIER_1D/*.F90)


add_library(hypar_lib_F SHARED ${hypar_SFE_sources})
target_include_directories(hypar_lib_F PUBLIC
    ${FFTW_INCLUDEDIR} ${FFTW_INCLUDE_DIRS}
    ${PETSC_INCLUDEDIR} ${PETSC_INCLUDE_DIRS}
  )

target_link_libraries(hypar_lib_F)

message(${ZLIB_LIBRARIES})
