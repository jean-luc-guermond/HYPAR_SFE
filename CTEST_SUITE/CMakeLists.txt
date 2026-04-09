# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

include(variables.cmake)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(regression_test)

include(CTest) # We include CTest which is part of CMake

#==========================================================================
#=============== compilation of HYPAR_SFE 
#==========================================================================

message (STATUS "${HYPAR_SFE_DIR}")


set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug_bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE ${release_bounds} CACHE STRING
	  " Flags used by the compiler during release builds. ")
set(CMAKE_Fortran_FLAGS_NATIVE ${native_bounds} CACHE STRING
	    " Flags used by the compiler during native builds. ")

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug")

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release")

ADD_CUSTOM_TARGET(native
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Native")

# Getting useful libraries
find_package(MPI)
find_package(ZLIB)
find_package(PkgConfig)

pkg_check_modules(PETSC PETSc)
pkg_check_modules(FFTW fftw3)


# Lib hypar_sfe_F
add_subdirectory(${HYPAR_SFE_DIR}/LIBS ${HYPAR_SFE_DIR}/LIBS/BUILD)
# Tests 1--5: Euler FE 1D/2D
add_subdirectory(${HYPAR_SFE_DIR}/REGRESSION_SUITE ${HYPAR_SFE_DIR}/REGRESSION_SUITE/BUILD)
# Test 6: FFT 1D
add_subdirectory(${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST6 ${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST6/BUILD)

enable_testing()
