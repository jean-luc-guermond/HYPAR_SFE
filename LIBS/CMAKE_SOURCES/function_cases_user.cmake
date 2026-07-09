function (create_cases_user)

# IN PARAMETERS:
    set(options)
    set(oneValueArgs SOURCES_DIR FE_DIM)
    set(multiValueArgs)
    cmake_parse_arguments(LOCAL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
# IN PARAMETERS:

#====================================================================================================
#=============== generating python environment using pip HYPAR_SFE for sources pre-generation
#====================================================================================================

    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    execute_process(
        COMMAND ${Python3_EXECUTABLE} ${HYPAR_SFE_DIR}/LIBS/PYTHON_SOURCES/generate_user_cases.py 
        ${HYPAR_SFE_DIR} ${LOCAL_FE_DIM} ${LOCAL_SOURCES_DIR}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        RESULT_VARIABLE result
    )
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "generate_user_cases.py failed")
    endif()

endfunction()