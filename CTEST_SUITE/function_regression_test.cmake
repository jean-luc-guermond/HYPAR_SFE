set(LIST_REGEX
"1234567891"
"1234567891.*2345678912"
"1234567891.*2345678912.*3456789123"
)

function(add_regression_test)
# IN PARAMETERS:
    set(options)
    set(oneValueArgs TEST_NAME NB_REGEX INCLUDE_TEMPLATE FE_DIM)
    set(multiValueArgs LIST_N_PROCS)
    cmake_parse_arguments(LOCAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
# IN PARAMETERS:

    message(STATUS "==> ${LOCAL_TEST_NAME}")
    message(STATUS "np_procs testing = ${LOCAL_LIST_N_PROCS}")
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
    set(MY_CMD /bin/bash "../../../job.sh")
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
    foreach(n_proc IN LISTS LOCAL_LIST_N_PROCS)
        add_test(NAME ${LOCAL_TEST_NAME}_PROC_${n_proc}
        WORKING_DIRECTORY ${TEST_DIR}/REGRESSION_TESTS
        COMMAND ${MY_CMD} ${n_proc})
        set_property (TEST "${LOCAL_TEST_NAME}_PROC_${n_proc}" PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
    endforeach()
endfunction()


#===========================================================================================
#======== FUNCTIONS GENERATING A RANDOM LIST TO TEST NUMBER OF PROCESSES RANDOMLY ==========
#===========================================================================================
function(generate_random_number OUT_INT MIN MAX)
    math(EXPR range_size "${MAX} - ${MIN} + 1")    
    string(RANDOM LENGTH 5 ALPHABET "0123456789" raw_val)
    math(EXPR final_val "(${raw_val} % ${range_size} + ${MIN})")
    set(${OUT_INT} "${final_val}" PARENT_SCOPE)

endfunction()

function(generate_random_procs_list OUT_VAR)

# IN PARAMETERS:
    set(options)
    set(oneValueArgs MIN MAX LEVEL DEFAULT)
    set(multiValueArgs LIST_LEVEL)
    cmake_parse_arguments(LOCAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
# IN PARAMETERS:

    if(LOCAL_LEVEL EQUAL 0)
        # Return result
        set(result ${LOCAL_DEFAULT})
        set(${OUT_VAR} "${result}" PARENT_SCOPE)
    else()
        # Get number of different tests to perform
        list(GET LOCAL_LIST_LEVEL ${LOCAL_LEVEL} N)

        # Safety check
        math(EXPR range_size "${LOCAL_MAX} - ${LOCAL_MIN} + 1")
        if(N GREATER range_size)
        message(FATAL_ERROR
            "Cannot generate ${N} unique numbers in range [${LOCAL_MIN}, ${LOCAL_MAX}]")
        endif()

        # Decide strategy
        math(EXPR range_size "${LOCAL_MAX} - ${LOCAL_MIN} + 1")    
        math(EXPR half_range "${range_size} / 2")
        if(N GREATER half_range)
            math(EXPR small_N "${range_size} - ${N}")
        else()
            math(EXPR small_N "${N}")
        endif()

        # List hosting random numbers
        set(tmp_list)

        while(TRUE)
            list(LENGTH tmp_list len)
            # Exit condition if found all random integers
            if(len EQUAL small_N)
                break()
            endif()

            # Generate random integer, check it was already generated before
            generate_random_number(RAND_INT ${LOCAL_MIN} ${LOCAL_MAX})
            # math(RANDOM OUTPUT_VARIABLE r RANGE ${MIN} ${MAX})
            list(FIND tmp_list ${RAND_INT} idx)
            if(idx EQUAL -1)
                list(APPEND tmp_list ${RAND_INT})
            endif()
        endwhile()

        # Assemble final list
        if(N GREATER half_range)
            # Take the complementary in this case
            set(result)
            foreach(v RANGE ${LOCAL_MIN} ${LOCAL_MAX})
                list(FIND tmp_list ${v} idx)
                if(idx EQUAL -1)
                    list(APPEND result ${v})
                endif()
            endforeach()
        else()
            set(result ${tmp_list})
            list(SORT result COMPARE NATURAL)
        endif()

        # Return result
        set(${OUT_VAR} "${result}" PARENT_SCOPE)
    endif()

endfunction()