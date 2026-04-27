set(LIST_REGEX
"1234567891"
"1234567891.*2345678912"
"1234567891.*2345678912.*3456789123"
)

function(add_regression_test)
# IN PARAMETERS:
    set(options)
    set(oneValueArgs TEST_NAME NB_REGEX INCLUDE_TEMPLATE FE_DIM)
    set(multiValueArgs)
    cmake_parse_arguments(LOCAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
# IN PARAMETERS:

    message("==> ${LOCAL_TEST_NAME}")
    set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_TEST_NAME})
    
    #===================================================
    # Compiling local sources
    #===================================================
    file(GLOB_RECURSE local_sources CONFIGURE_DEPENDS ${TEST_DIR}/LOCAL_SOURCES/*.F90)   

    IF (LOCAL_INCLUDE_TEMPLATE)
        file(GLOB_RECURSE template_sources CONFIGURE_DEPENDS ${TEST_DIR}/../TEMPLATE/*.F90)
        list(APPEND local_sources ${template_sources})
    ENDIF()

    #===================================================
    # Loop over number of regex to generate executables
    #===================================================
    set(MY_CMD /bin/bash job.sh)
    list(APPEND MY_CMD
        ${RUN_PRE_PROC}
        ${PROC_CALL}
        ${RUN_POST_PROC}
        ${LOCAL_NB_REGEX}
    )

    foreach(i RANGE 1 ${LOCAL_NB_REGEX})
      # Create executable
      set(exe "${LOCAL_TEST_NAME}_${i}.exe")
      add_executable(${exe} ${local_sources})
      
      set_target_properties(${exe} PROPERTIES
      Fortran_MODULE_DIRECTORY ${TEST_DIR}/BUILD/BUILD_${i}/mod RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR}/EXECUTABLE)
      target_compile_definitions(${exe} PRIVATE regex_num=${i})

      # Link with hypar, petsc, fftw, mpi, zlib, and eventually additional libraries
      target_link_libraries(${exe} PRIVATE hypar_lib_F_${LOCAL_FE_DIM})
      target_include_directories(${exe} PRIVATE ${HYPAR_SFE_DIR}/LIBS/BUILD/${LOCAL_FE_DIM}/mod)
      LIST(APPEND MY_CMD ${exe})
    endforeach()

    #===================================================
    # Define test parameters
    #===================================================

    #=== chain of characters to define test validation
    math(EXPR index_regex "${LOCAL_NB_REGEX} - 1")
    list(GET LIST_REGEX ${index_regex} passRegex)

    #=== effective test definition
    add_test(NAME ${LOCAL_TEST_NAME}
      WORKING_DIRECTORY ${TEST_DIR}/REGRESSION_TESTS
      COMMAND ${MY_CMD})
    set_property (TEST ${LOCAL_TEST_NAME} PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
endfunction()