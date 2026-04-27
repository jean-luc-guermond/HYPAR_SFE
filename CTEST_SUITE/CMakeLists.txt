# Cmake minimum version
cmake_minimum_required(VERSION 3.29)
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Project name and version
project(REGRESSION_TEST LANGUAGES C Fortran)
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

include(variables.cmake)
define_compile_mode()
include(function_regression_test.cmake)
include(CTest)

#==========================================================================
#=============== compilation of HYPAR_SFE libs
#==========================================================================

# Lib hypar_sfe_F
set(BUILD_ALL_DIMS ON) 
add_subdirectory(${HYPAR_SFE_DIR}/LIBS ${HYPAR_SFE_DIR}/LIBS/BUILD)

#==========================================================================
#=============== Listing and compilation of CTESTS
#==========================================================================

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
add_subdirectory(${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_PERIODICITY/
  ${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST_PERIODICITY/BUILD)