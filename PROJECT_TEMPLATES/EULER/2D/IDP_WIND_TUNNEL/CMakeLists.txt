# Cmake minimum version
cmake_minimum_required(VERSION 3.29)
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Project name and version
project(MY_PROJECT VERSION 1.0 LANGUAGES Fortran C)
MESSAGE(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")
MESSAGE(STATUS "")

# Define variables for compilation
set(FE_dim "0")
include(variables.cmake)
if(FE_dim EQUAL 0)
    MESSAGE(FATAL_ERROR "BUG IN MY_APP/CMakeLists.txt: FE_dim = 0 ==> very likely means you forgot to define it in your variable.cmake")
elseif(FE_dim GREATER 2)
    MESSAGE(FATAL_ERROR "BUG IN MY_APP/CMakeLists.txt: FE_dim = ${FE_dim} unavailable, only 1 or 2")
endif()

#=============================================================
#============ SETTING VARIABLES FOR COMPILATION  =============
#=============================================================

set(CMAKE_Fortran_FLAGS_DEBUG "${debug_bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE "${release_bounds}")
set(CMAKE_Fortran_FLAGS_NATIVE "${native_bounds}")# CACHE STRING
# set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug_bounds}")
# set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${release_bounds}")
# set(CMAKE_Fortran_FLAGS_NATIVE "${CMAKE_Fortran_FLAGS_NATIVE} ${native_bounds}")# CACHE STRING
        # " Flags used by the compiler during native builds. ")

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


# Create a list of the files to be preprocessed
# set(fppFiles limiter.fpp)
file(GLOB_RECURSE fppFiles CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/LOCAL_SOURCES/TEMPLATED/*.fpp)

# Pre-process
foreach(infile IN LISTS fppFiles)

    # Generate output file name
    string(REGEX REPLACE ".fpp\$" ".F90" outfile "${infile}")

    # Custom command to do the processing
    add_custom_command(
        OUTPUT "${outfile}"
        COMMAND fypp "${infile}" "${outfile}"
        MAIN_DEPENDENCY "${infile}"
        VERBATIM)

    # Finally add output file to a list
    set(outFiles ${outFiles} "${outfile}")
    # message(STATUS "${outfile}")
    # message(STATUS "${outFiles}")
endforeach(infile)


# Building local sources
file(GLOB_RECURSE local_sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/LOCAL_SOURCES/*.F90)
list(APPEND local_sources ${outFiles})
add_executable(${exe} ${local_sources})

set_target_properties(${exe} PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/EXECUTABLE)
# target_compile_definitions(${exe} PRIVATE k_dim=${i})
# target_compile_definitions(${exe} PUBLIC USE_LIGHTKRYLOV)

# Link local_sources with correct hypar lib
target_link_libraries(${exe} PRIVATE ${MY_HYPAR_SFE_LIB})
target_include_directories(${exe} PRIVATE ${MY_HYPAR_SFE_MODULE_SOURCES})
