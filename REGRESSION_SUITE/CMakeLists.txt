# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# include(variables.cmake)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(ctest_HYPAR_SFE)


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
#=============== function compiling all tests
#==========================================================================

function(add_regression_test TEST_NAME TEST_NB TEST_DIM)
    
    message(${TEST_NAME})
    set(TEST_DIR ${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST${TEST_NB})

    # Create local sources from local dir
    file(GLOB_RECURSE local_sources 
        ${TEST_DIR}/LOCAL_SOURCES/*.F90
        ${TEST_DIR}/../TEMPLATE_EULER_FE/*.F90)

    # Define runtime output directory
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR}/EXECUTABLE)

    # Add executable
    set(exe "test${TEST_NB}.exe")
    add_executable(${exe} ${local_sources})

    #===================================================
    # Decide which hypar_lib to use based on dimension
    #===================================================
    # if(TEST_DIM STREQUAL "1D")
    #     set(HYPAR_LIB hypar_lib_F_1D)
    # elseif(TEST_DIM STREQUAL "2D")
    #     set(HYPAR_LIB hypar_lib_F_2D)
    # else()
    #     message(FATAL_ERROR "Invalid dimension '${TEST_DIM}'. Allowed: 1D, 2D")
    # endif()
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

    message(${ZLIB_LIBRARIES})
    target_include_directories(${exe} PUBLIC ${HYPAR_SFE_DIR}/LIBS/BUILD)
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
endfunction()

#==========================================================================
#=============== list of tests, comment out whichever you don't want to run 
#==========================================================================

set(LIST_TESTS 
1 1D IDP_SODT_1D
2 1D IDP_SMOOTH_1D
3 1D IDP_VANDERWALLS_1D
4 2D IDP_ISENTROPIC_VORTEX
5 2D IDP_WIND_TUNNEL
)

#==========================================================================
#=============== Compilation of all tests
#==========================================================================

list(LENGTH LIST_TESTS N)
math(EXPR N "${N} - 1")

foreach(i RANGE 0 ${N} 3)
    list(GET LIST_TESTS ${i} TEST_NB)
    math(EXPR j "${i} + 1")
    list(GET LIST_TESTS ${j} TEST_DIM)
    math(EXPR k "${i} + 2")
    list(GET LIST_TESTS ${k} TEST_NAME)

    add_regression_test(${TEST_NAME} ${TEST_NB} ${TEST_DIM})
endforeach()

enable_testing()
