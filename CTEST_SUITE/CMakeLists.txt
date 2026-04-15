# Cmake minimum version
cmake_minimum_required(VERSION 3.29)

include(variables.cmake)

# Project name and version
project(regression_test)

# Enable Fortran
enable_language(Fortran)

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


#==========================================================================
#=============== compilation of CTESTS
#==========================================================================

set(LIST_REGEX
"1234567891"
"1234567891.*2345678912"
"1234567891.*2345678912.*3456789123"
)

# Lib hypar_sfe_F
add_subdirectory(${HYPAR_SFE_DIR}/LIBS ${HYPAR_SFE_DIR}/LIBS/BUILD)


# TEST_FE_1D_EULER
add_subdirectory(${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_FE_1D_EULER
  ${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_FE_1D_EULER/BUILD)

# TEST_FE_2D_EULER
add_subdirectory(${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_FE_2D_EULER
  ${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_FE_2D_EULER/BUILD)

# TEST_FFT_1D
add_subdirectory(${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_FFT
  ${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_FFT/BUILD)

# TEST_FOURIER_1D_SCALAR
add_subdirectory(${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_FOURIER_1D_SCALAR
  ${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_FOURIER_1D_SCALAR/BUILD)

# TEST_PERIODICITY
add_subdirectory(${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_PERIODICITY/PERIODIC_1D_SCAL_ADV_DIFF
  ${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_PERIODICITY/PERIODIC_1D_SCAL_ADV_DIFF/BUILD)

enable_testing()
