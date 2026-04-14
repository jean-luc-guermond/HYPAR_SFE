# Cmake minimum version
cmake_minimum_required(VERSION 3.29)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(ctest_EULER_FE VERSION 1.0)


#==========================================================================
#=============== TARGET debug
#==========================================================================

  IF (TARGET debug)
  else()
  ADD_CUSTOM_TARGET(debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all -j20
    COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
    )

  ADD_CUSTOM_TARGET(release
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all -j20
    COMMENT "Switch CMAKE_BUILD_TYPE to Release"
    )

  ADD_CUSTOM_TARGET(native
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all -j20
    COMMENT "Switch CMAKE_BUILD_TYPE to Native"
    )
  endif()


#==========================================================================
#=============== list of tests, comment out whichever you don't want to run
#==========================================================================

set(LIST_TESTS
  1 FE_2D_ISENTROPIC_VORTEX
  1 FE_2D_WIND_TUNNEL
)


#==========================================================================
#=============== function compiling all tests
#==========================================================================

function(add_regression_test TEST_NAME NB_REGEX)

    message(${TEST_NAME})
    set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME})

    # Create local sources from local dir
    file(GLOB_RECURSE local_sources
        ${TEST_DIR}/LOCAL_SOURCES/*.F90
        ${TEST_DIR}/../TEMPLATE/*.F90)

    # Define runtime output directory
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR}/EXECUTABLE)

    #===================================================
    # Define hypar library
    #===================================================

    set(HYPAR_LIB hypar_lib_F)

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

    #===================================================
    # Define test parameters
    #===================================================

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
endfunction()


#==========================================================================
#=============== Compilation of all tests
#==========================================================================

list(LENGTH LIST_TESTS N)
math(EXPR N "${N} - 1")

foreach(i RANGE 0 ${N} 2)
    math(EXPR j "${i} ")
    list(GET LIST_TESTS ${j} NB_REGEX)
    math(EXPR j "${i} + 1")
    list(GET LIST_TESTS ${j} TEST_NAME)

    add_regression_test(${TEST_NAME} ${NB_REGEX})
endforeach()

enable_testing()
