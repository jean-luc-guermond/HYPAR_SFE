set(LIST_REGEX
"1234567891"
"1234567891.*2345678912"
"1234567891.*2345678912.*3456789123"
"1234567891.*2345678912.*3456789123*45678901234"
"1234567891.*2345678912.*3456789123*4567891234*5678912345"
"1234567891.*2345678912.*3456789123*4567891234*5678912345*6789123456"
"1234567891.*2345678912.*3456789123*4567891234*5678912345*6789123456*78901234567"
)

function(add_regression_test)
# IN PARAMETERS:
    set(options)
    set(oneValueArgs TEST_DIR 
                    TEST_NAME
                    INCLUDE_TEST_SOURCES
                    NB_REGEX
                    FE_DIM
                    PROC_MIN PROC_MAX DEFAULT_PROC)
    set(multiValueArgs TEST_LABELS CASE_LABELS)
    cmake_parse_arguments(LOCAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
# IN PARAMETERS:

    #===================================================
    # Generate CPU configurations for MPI tests
    #===================================================

    if((NOT DEFINED LOCAL_PROC_MIN) OR (NOT DEFINED LOCAL_PROC_MAX))
        set(LOCAL_PROC_MIN 1)
        set(LOCAL_PROC_MAX 0)
        set(NB_PARALLEL_TEST 0)
    endif()
    generate_random_procs_list(
        MPI_LIST_N_PROCS 
        MIN ${LOCAL_PROC_MIN}
        MAX ${LOCAL_PROC_MAX}
        NB_TESTS ${NB_PARALLEL_TEST})
    set(LOCAL_LIST_N_PROCS ${LOCAL_DEFAULT_PROC} ${MPI_LIST_N_PROCS})
    set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_TEST_DIR})
    
    #===================================================
    # Compiling local sources
    #===================================================

    set(local_sources_dir ${CMAKE_CURRENT_SOURCE_DIR})
    string(REPLACE "REGRESSION_SUITE" "CASES_USER" local_sources_dir "${local_sources_dir}")
    file(GLOB_RECURSE local_sources CONFIGURE_DEPENDS ${local_sources_dir}/LOCAL_SOURCES/*.F90)   
    if(LOCAL_INCLUDE_TEST_SOURCES)
        file(GLOB_RECURSE test_sources CONFIGURE_DEPENDS ${TEST_DIR}/LOCAL_SOURCES/*.F90)   
        list(APPEND local_sources ${test_sources})
    endif()

    #===================================================
    # Selecting Hypar library configuration
    #===================================================

    if(NOT DEFINED static_library)
        set(static_library OFF)
    endif()
    if(static_library)
        set(MY_HYPAR_SFE_LIB hypar_lib_F_${LOCAL_FE_DIM}_static)
    else()
        set(MY_HYPAR_SFE_LIB hypar_lib_F_${LOCAL_FE_DIM}_dynamic)
    endif()

    #===================================================
    # Loop over number of regex to generate executables
    #===================================================
    set(EXEC_CMD /bin/bash "${HYPAR_SFE_DIR}/REGRESSION_SUITE/job.sh")
    set(MY_CMD 
        "${RUN_PRE_PROC} ${RUN_POST_PROC}"
        "${PROC_CALL}"
        "${LOCAL_NB_REGEX}"
        )

    foreach(i RANGE 1 ${LOCAL_NB_REGEX})
      # Create executable
      set(exe "${LOCAL_TEST_NAME}_${i}.exe")
      add_executable(${exe} ${local_sources})
      
      set_target_properties(${exe} PROPERTIES
      Fortran_MODULE_DIRECTORY ${TEST_DIR}/BUILD/BUILD_${i}/mod RUNTIME_OUTPUT_DIRECTORY ${TEST_DIR}/EXECUTABLE)
      target_compile_definitions(${exe} PRIVATE regex_num=${i})

      # Link with hypar, petsc, fftw, mpi, zlib, and eventually additional libraries
      target_link_libraries(${exe} PRIVATE ${MY_HYPAR_SFE_LIB})
      target_include_directories(${exe} PRIVATE ${HYPAR_SFE_DIR}/LIBS/BUILD/${LOCAL_FE_DIM}/mod)
      LIST(APPEND MY_CMD ${exe})
    endforeach()

    #===================================================
    # Define test parameters
    #===================================================

    #=== chain of characters to define test validation
    math(EXPR index_regex "${LOCAL_NB_REGEX} - 1")
    list(GET LIST_REGEX ${index_regex} passRegex)

    #=== Finally the test definition (default + potential mpi configurations)
    set(n_test 0)
    foreach(n_proc IN LISTS LOCAL_LIST_N_PROCS)
        math(EXPR n_test "${n_test} + 1")
        if(n_test EQUAL 1)
            set(TEST_LABEL "${LOCAL_TEST_LABELS};${LOCAL_CASE_LABELS}")
        elseif(NOT n_proc EQUAL LOCAL_DEFAULT_PROC)
            set(TEST_LABEL "${LOCAL_TEST_LABELS};mpi;${LOCAL_CASE_LABELS}")
        else()
            continue()
        endif()
        add_test(
            NAME ${LOCAL_TEST_NAME}_PROC_${n_proc}
            WORKING_DIRECTORY ${TEST_DIR}/REGRESSION_TESTS
            COMMAND ${EXEC_CMD} ${MY_CMD} ${n_proc}
            )
        set_tests_properties(
                ${LOCAL_TEST_NAME}_PROC_${n_proc}
                PROPERTIES
                    LABELS "${TEST_LABEL}"
                    PASS_REGULAR_EXPRESSION "${passRegex}"
            )       
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
    set(oneValueArgs MIN MAX NB_TESTS)
    set(multiValueArgs)
    cmake_parse_arguments(LOCAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
# IN PARAMETERS:

    if(LOCAL_NB_TESTS LESS_EQUAL 0)
        # Return result
        set(result)
        set(${OUT_VAR} "${result}" PARENT_SCOPE)
    else()
        # Get number of different tests to perform
        set(N ${LOCAL_NB_TESTS})
        # Safety check
        math(EXPR range_size "${LOCAL_MAX} - ${LOCAL_MIN} + 1")
        if(N GREATER range_size)
            set(result "")
            foreach(i RANGE ${LOCAL_MIN} ${LOCAL_MAX})
                list(APPEND result ${i})
            endforeach()

            set(${OUT_VAR} "${result}" PARENT_SCOPE)
            return()
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