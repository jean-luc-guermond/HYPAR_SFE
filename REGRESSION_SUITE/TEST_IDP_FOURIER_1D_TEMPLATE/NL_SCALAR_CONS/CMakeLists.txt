# Cmake minimum version
cmake_minimum_required(VERSION 3.16)

# include(variables.cmake)

# Enable Fortran
enable_language(Fortran)

# Project name and version
project(ctest_NL_SCALAR_CONS VERSION 1.0)


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

function(add_regression_test TEST_NAME TEST_NB NB_REGEX)
    
    message(${TEST_NAME})
    set(TEST_DIR ${HYPAR_SFE_DIR}/REGRESSION_SUITE/TEST${TEST_NB})

    # Create local sources from local dir
    file(GLOB_RECURSE local_sources 
        ${TEST_DIR}/LOCAL_SOURCES/*.F90
        ${TEST_DIR}/../TEST_IDP_FOURIER_1D_TEMPLATE/NL_SCALAR_CONS/*.F90)

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
    #===================================================
    foreach(i RANGE 1 ${NB_REGEX})
        # Add executable
        set(exe "${TEST_NAME}_${i}.exe")
        # set(exe "test${TEST_NB}_${i}.exe")
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
                ${ADDITIONAL_LINKS})
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
    # set(passRegex "1234567891")
   
    #=== test name definition, used for test output
    set (nm "${TEST_NAME}")

    #=== effective test definition
    add_test(NAME ${nm}
      WORKING_DIRECTORY ${TEST_DIR}/REGRESSION_TESTS
      COMMAND ${MY_CMD})
    set_property (TEST ${nm} PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
endfunction()

#==========================================================================
#=============== list of tests, comment out whichever you don't want to run 
#==========================================================================

set(LIST_TESTS 
7 2 burgers_fft
8 2 linear_transport_fft
9 1 non_convex_fft
)

#==========================================================================
#=============== Compilation of all tests
#==========================================================================

list(LENGTH LIST_TESTS N)
math(EXPR N "${N} - 1")

foreach(i RANGE 0 ${N} 3)
    list(GET LIST_TESTS ${i} TEST_NB)
    math(EXPR j "${i} + 1")
    list(GET LIST_TESTS ${j} NB_REGEX)
    math(EXPR j "${i} + 2")
    list(GET LIST_TESTS ${j} TEST_NAME)

    add_regression_test(${TEST_NAME} ${TEST_NB} ${NB_REGEX})
endforeach()

enable_testing()
