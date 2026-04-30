# Cmake minimum version
cmake_minimum_required(VERSION 3.29)
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Project name and version
project(REGRESSION_TEST LANGUAGES C Fortran)
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

include(variables.cmake)
# define_compile_mode()

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