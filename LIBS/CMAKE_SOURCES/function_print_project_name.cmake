function(print_case_name)
# IN PARAMETERS:
    set(options)
    set(oneValueArgs 
                    TEST_NAME
                    )
    set(multiValueArgs)
    cmake_parse_arguments(LOCAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
# IN PARAMETERS:

    string(LENGTH "${LOCAL_TEST_NAME}" len)
    string(REPEAT "=" 5 partial_eq_char)
    math(EXPR n_equal "5 + 1 + ${len} + 1 + 5")
    string(REPEAT "=" ${n_equal} full_eq_char)
    MESSAGE(STATUS ${full_eq_char})
    MESSAGE(STATUS "${partial_eq_char} ${LOCAL_TEST_NAME} ${partial_eq_char}")
    MESSAGE(STATUS ${full_eq_char})
endfunction()
