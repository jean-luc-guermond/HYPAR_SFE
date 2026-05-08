# Cmake minimum version
cmake_minimum_required(VERSION 3.29)
MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Project name and version
project(IDP_SODT_LIMITING VERSION 1.0 LANGUAGES Fortran C)
message(STATUS "===== compiling project ${PROJECT_NAME} =====")

# Define variables for compilation
include(variables.cmake)
define_compile_mode()
set(exe "a.exe")

# Building HYPAR_SFE library
add_subdirectory(${HYPAR_SFE_DIR}/LIBS ${HYPAR_SFE_DIR}/LIBS/BUILD)
set(MY_HYPAR_SFE_LIB hypar_lib_F_${FE_dim}D)
set(MY_HYPAR_SFE_MODULE_SOURCES ${HYPAR_SFE_DIR}/LIBS/BUILD/${FE_dim}D/mod)

# Building local sources
file(GLOB_RECURSE local_sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/LOCAL_SOURCES/*.F90)   
add_executable(${exe} ${local_sources})

set_target_properties(${exe} PROPERTIES
RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/EXECUTABLE)
# target_compile_definitions(${exe} PRIVATE k_dim=${i})
# target_compile_definitions(${exe} PUBLIC USE_LIGHTKRYLOV)

# Link local_sources with correct hypar lib
target_link_libraries(${exe} PRIVATE ${MY_HYPAR_SFE_LIB})
target_include_directories(${exe} PRIVATE ${MY_HYPAR_SFE_MODULE_SOURCES})
