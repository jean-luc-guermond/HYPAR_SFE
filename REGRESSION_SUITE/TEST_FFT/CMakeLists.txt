# Cmake minimum version
cmake_minimum_required(VERSION 3.29)

# include(variables.cmake)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(ctest_FFT_1D)


#==========================================================================
#=============== Compilation of test
#==========================================================================
set(TEST_NAME "FFT_1D")
message(${TEST_NAME})
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Create local sources from local dir
file(GLOB_RECURSE local_sources
    ${TEST_DIR}/LOCAL_SOURCES/*.F90)

# Define runtime output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR}/EXECUTABLE)

# Add executable
set(exe "${TEST_NAME}.exe")
add_executable(${exe} ${local_sources})

set(HYPAR_LIB hypar_lib_F)
#===================================================
# Link with the rest
#===================================================
target_link_libraries(${exe} PUBLIC ${HYPAR_LIB}
        ${PETSC_LINK_LIBRARIES}
        ${PETSC_PRIVATE_LIBRARIES}
        ${FFTW_LINK_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${MPI_Fortran_LIBRARIES}
        ${ADDITIONAL_LINKS}
    )
target_include_directories(${exe} PUBLIC ${HYPAR_SFE_DIR}/LIBS/BUILD
    ${FFTW_INCLUDEDIR} ${FFTW_INCLUDE_DIRS}
    ${PETSC_INCLUDEDIR} ${PETSC_INCLUDE_DIRS}
  )
#===================================================
# Define test parameters
#===================================================

#=== chain of characters to define test validation
set(passRegex "1234567891")

#=== test name definition, used for test output
set (nm "${TEST_NAME}")

#=== effective test definition
add_test(NAME ${nm}
    WORKING_DIRECTORY ${TEST_DIR}/REGRESSION_TESTS
    COMMAND /bin/bash/ job.sh ${RUN_PRE_PROC} ${PROC_CALL} ${exe} ${RUN_POST_PROC})
set_property (TEST ${nm} PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

enable_testing()
