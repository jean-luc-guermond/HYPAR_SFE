#====================================================================================================
#=============== generating python environment using pip HYPAR_SFE for sources pre-generation
#====================================================================================================

find_package(Python3 REQUIRED COMPONENTS Interpreter)

execute_process(
    COMMAND
        ${Python3_EXECUTABLE} -c
        "import fypp"
    RESULT_VARIABLE HAVE_FYPP
    OUTPUT_QUIET
    ERROR_QUIET
)

if(HAVE_FYPP EQUAL 0)
    message(STATUS "Using system Python with FYPP.")
else()
    set(VENV_DIR ${HYPAR_SFE_DIR}/LIBS/DEPENDENCIES/.venv)
    set(VENV_PYTHON ${VENV_DIR}/bin/python)
    if(EXISTS "${VENV_DIR}")
        message(STATUS "Using existing virtual environment: ${VENV_DIR}")
    else()
        message(STATUS "FYPP not found. Creating virtual environment...")
        # Create the virtual environment
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -m venv ${VENV_DIR}
            RESULT_VARIABLE VENV_RESULT
        )
        if(NOT VENV_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to create Python virtual environment")
        endif()
        # Install dependencies
        execute_process(
            COMMAND ${VENV_PYTHON} -m pip install fypp
            RESULT_VARIABLE PIP_RESULT
        )
        if(NOT PIP_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to install fypp")
        endif()
    endif()
    set(Python3_EXECUTABLE ${VENV_PYTHON})
endif()
