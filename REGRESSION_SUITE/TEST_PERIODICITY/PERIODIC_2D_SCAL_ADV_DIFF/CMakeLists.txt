# Cmake minimum version
cmake_minimum_required(VERSION 3.29)

# include(variables.cmake)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(ctest_PERIODIC_2D_SCAL_ADV_DIFF)


#==========================================================================
#=============== Compilation of test
#==========================================================================
set(TEST_NAME "PERIODIC_2D_SCAL_ADV_DIFF")
set(NB_REGEX 3)

message(${TEST_NAME})
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Create local sources from local dir
file(GLOB_RECURSE local_sources
    ${TEST_DIR}/LOCAL_SOURCES/*.F90)

# Define runtime output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR}/EXECUTABLE)

set(HYPAR_LIB hypar_lib_F)
#===================================================
# Link with the rest
#===================================================


set(MY_CMD
    /bin/bash
    job.sh
    ${RUN_PRE_PROC}
    ${PROC_CALL}
    ${RUN_POST_PROC}
    ${NB_REGEX}
)
foreach(i RANGE 1 ${NB_REGEX})

  # Add executable
  set(exe "${TEST_NAME}_${i}.exe")
  add_executable(${exe} ${local_sources})
  set_target_properties(${exe} PROPERTIES
    Fortran_MODULE_DIRECTORY ${TEST_DIR}/BUILD_${i}/mod)
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

  target_include_directories(${exe} PUBLIC ${HYPAR_SFE_DIR}/LIBS/BUILD)
  target_compile_definitions(${exe} PRIVATE regex_num=${i})
  LIST(APPEND MY_CMD ${exe})
endforeach()

#=== chain of characters to define test validation
math(EXPR index_regex "${NB_REGEX} - 1")
list(GET LIST_REGEX ${index_regex} passRegex)

#=== test name definition, used for test output
set (nm "${TEST_NAME}")

#=== effective test definition
add_test(NAME ${nm}
  WORKING_DIRECTORY ${TEST_DIR}/REGRESSION_TESTS
  COMMAND ${MY_CMD})
set_property (TEST ${nm} PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")

enable_testing()
