# Cmake minimum version
cmake_minimum_required(VERSION 3.29)
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Project name and version
list(APPEND CMAKE_MODULE_PATH
     "${CMAKE_SOURCE_DIR}"
)
set(VAR_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/variables_template.cmake")
set(VAR "${CMAKE_CURRENT_SOURCE_DIR}/variables.cmake")
if(NOT EXISTS "${VAR}")
    file(COPY_FILE "${VAR_TEMPLATE}" "${VAR}")
endif()

include(variables.cmake)
get_filename_component(
    HYPAR_SFE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}"
    DIRECTORY
)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}" "${HYPAR_SFE_DIR}/LIBS/CMAKE_SOURCES")
include(function_regression_test)
include(function_cases_user)
include(function_print_project_name)

project(REGRESSION_TEST LANGUAGES C Fortran)
MESSAGE(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")
MESSAGE(STATUS "")

include(CTest)

#==========================================================================
#=============== messaging concerning dynamic/static
#==========================================================================

if(${static_library})
    MESSAGE(WARNING "You selected the static library => be aware that all executables might be very heavy")
endif()
MESSAGE(STATUS "")

#==========================================================================
#=============== setting compile flags
#==========================================================================
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
#=============== compilation of HYPAR_SFE libs
#==========================================================================

# Lib hypar_sfe_F
set(is_ctest True)
set(BUILD_ALL_DIMS ON) 
add_subdirectory(${HYPAR_SFE_DIR}/LIBS ${HYPAR_SFE_DIR}/LIBS/BUILD)

#==========================================================================
#=============== Point towards actual CTESTS
#==========================================================================

set(DIR_TESTS ${HYPAR_SFE_DIR}/REGRESSION_SUITE)
add_subdirectory(${DIR_TESTS}
                 ${DIR_TESTS}/BUILD)