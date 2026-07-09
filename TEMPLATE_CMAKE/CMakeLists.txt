#=================================================
#============ INIT CMAKE ENVIRONMENT =============
#=================================================

# Cmake minimum version
cmake_minimum_required(VERSION 3.29)
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Project name and version
# Define variables for compilation
set(FE_dim "0")
include(variables.cmake)
if(FE_dim EQUAL 0)
    MESSAGE(FATAL_ERROR "BUG IN MY_APP/CMakeLists.txt: FE_dim = 0 ==> very likely means you forgot to define it in your variable.cmake")
elseif(FE_dim GREATER 2)
    MESSAGE(FATAL_ERROR "BUG IN MY_APP/CMakeLists.txt: FE_dim = ${FE_dim} unavailable, only 1 or 2")
endif()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}" "${HYPAR_SFE_DIR}/LIBS/CMAKE_SOURCES")

#=============================================
#============ PROJECT DEFINITION =============
#=============================================

project(BURGERS_1D VERSION 1.0 LANGUAGES Fortran C)
MESSAGE(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")
MESSAGE(STATUS "")

#=============================================================
#============ SETTING VARIABLES FOR COMPILATION  =============
#=============================================================
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  add_compile_options($<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug_bounds}" CACHE STRING
        " Flags used by the compiler during debug builds. " FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${release_bounds}" CACHE STRING
        " Flags used by the compiler during release builds. " FORCE)
set(CMAKE_Fortran_FLAGS_NATIVE "${native_bounds}" CACHE STRING
        " Flags used by the compiler during native builds. " FORCE)

ADD_CUSTOM_TARGET(debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to Debug")

ADD_CUSTOM_TARGET(release
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to Release")

ADD_CUSTOM_TARGET(native
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
    COMMENT "Switch CMAKE_BUILD_TYPE to Native")


#==========================================================================
#=============== HYPAR_SFE sources pre-generation (using python and fypp)
#==========================================================================

if(python_precompile)
    include(CMake_get_python)

    execute_process(
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/python_precompile.py ${HYPAR_SFE_DIR}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        RESULT_VARIABLE result
    )

    if(NOT result EQUAL 0)
        message(FATAL_ERROR "BUG in ${CMAKE_SOURCE_DIR}/CMakeLists.txt: python_precompile.py failed")
    endif()

    file(GLOB_RECURSE fppFiles CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/LOCAL_SOURCES/*.fpp)

    # Pre-process with fpp files
    foreach(infile IN LISTS fppFiles)

        # Generate output file name
        string(REGEX REPLACE ".fpp\$" ".F90" outfile "${infile}")

        # Custom command to do the processing
        add_custom_command(
            OUTPUT "${outfile}"
            COMMAND ${Python3_EXECUTABLE} -m fypp "${infile}" "${outfile}"
            MAIN_DEPENDENCY "${infile}"
            VERBATIM)
        # Finally add output file to a list
        set(outFiles ${outFiles} "${outfile}")
    endforeach(infile)
endif()

#=============================================================
#============ COMPILING AND LINKING WITH HYPAR   =============
#=============================================================

set(exe "a.exe")


# Building HYPAR_SFE library
add_subdirectory(${HYPAR_SFE_DIR}/LIBS ${HYPAR_SFE_DIR}/LIBS/BUILD)

if(NOT DEFINED static_library)
    set(static_library OFF)
endif()
if(static_library)
    set(MY_HYPAR_SFE_LIB hypar_lib_F_${FE_dim}D_static)
    message(STATUS "Linking with STATIC library")
else()
    set(MY_HYPAR_SFE_LIB hypar_lib_F_${FE_dim}D_dynamic)
    message(STATUS "Linking with DYNAMIC library")
endif()

# set(MY_HYPAR_SFE_LIB hypar_lib_F_${FE_dim}D)
set(MY_HYPAR_SFE_MODULE_SOURCES ${HYPAR_SFE_DIR}/LIBS/BUILD/${FE_dim}D/mod)

# Building local sources
file(GLOB_RECURSE local_sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/LOCAL_SOURCES/*.F90)   
add_executable(${exe} ${local_sources})

if(python_precompile)
    add_custom_target(fypp_generated DEPENDS ${outFiles})
    add_dependencies(${exe} fypp_generated)
endif()

set_target_properties(${exe} PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/EXECUTABLE)
# target_compile_definitions(${exe} PRIVATE k_dim=${i})
# target_compile_definitions(${exe} PUBLIC USE_LIGHTKRYLOV)

# Link local_sources with correct hypar lib
target_link_libraries(${exe} PRIVATE ${MY_HYPAR_SFE_LIB})
target_include_directories(${exe} PRIVATE ${MY_HYPAR_SFE_MODULE_SOURCES})
