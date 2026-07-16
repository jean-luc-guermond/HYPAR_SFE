# Project name and version
project(HYPAR_SFE_LIBS)
MESSAGE(STATUS "")
message(STATUS "===== compiling project ${PROJECT_NAME} =====")
MESSAGE(STATUS "")

include(variables_lib.cmake)

if(FE_dim EQUAL 0)
    MESSAGE(FATAL_ERROR "BUG IN HYPAR_SFE/LIBS/CMakeLists.txt: FE_dim = 0 ==> very likely means you forgot to define it in your variable.cmake")
elseif(FE_dim GREATER 2)
    MESSAGE(FATAL_ERROR "BUG IN HYPAR_SFE/LIBS/CMakeLists.txt: FE_dim = ${FE_dim} unavailable, only 1 or 2")
endif()
#==========================================================================
#=============== Loading of packages and checking versions
#==========================================================================

# Getting mandatory libraries
find_package(MPI REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)
# Getting mandatory libraries

# Debugging information for libraries
message(STATUS "Working in: ${HYPAR_SFE_DIR}")
message(STATUS "ZLIB_LIBRARIES=${ZLIB_LIBRARIES}")
message(STATUS "MPI_Fortran_LIBRARIES=${MPI_Fortran_LIBRARIES}")

#======================================#
#==== load and check PETSc version ====#
#======================================#

#=== Get Petsc
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

#=== Verbose version and paths
get_target_property(PETSC_INC PkgConfig::PETSC INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "PETSC includes: ${PETSC_INC}")
message(STATUS "PETSC version: ${PETSC_VERSION}")

#=== Version checks
if(${PETSC_VERSION} VERSION_LESS ${PETSC_MIN_VERSION})
    MESSAGE(FATAL_ERROR "PETSC version: must be at least ${PETSC_MIN_VERSION}")
elseif(${PETSC_VERSION} VERSION_LESS ${PETSC_MIN_VERSION_WARNING})
    MESSAGE(WARNING "It is strongly encouraged to used at least PETSC version ${PETSC_MIN_VERSION_WARNING}")
elseif((${PETSC_VERSION} VERSION_GREATER ${PETSC_MAX_VERSION}) OR (${PETSC_VERSION} VERSION_EQUAL ${PETSC_MAX_VERSION}))
    MESSAGE(FATAL_ERROR "PETSC version: must be at most ${PETSC_MAX_VERSION}")
endif()

#=== Prepare macros for compilation
string(REPLACE "." ";" PETSC_VERSION_LIST "${PETSC_VERSION}")
# list(GET PETSC_VERSION_LIST 0 PETSC_VERSION_MAJOR)
list(GET PETSC_VERSION_LIST 1 PETSC_VERSION_MINOR)
# list(GET PETSC_VERSION_LIST 2 PETSC_VERSION_SUBMINOR)
if(${PETSC_VERSION} VERSION_LESS "3.23")
    set(PETSC_INCLUDE_MACRO "#include <petsc/finclude/petsc.h>")
else()
    set(PETSC_INCLUDE_MACRO "")
endif()
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/CMAKE_SOURCES/petsc_include.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/petsc_include.h
    @ONLY
)
#======================================#
#==== load and check PETSc version ====#
#======================================#

#===================#
#==== load FFTw ====#
#===================#
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)

get_target_property(FFTW_INC PkgConfig::FFTW INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "FFTW includes: ${FFTW_INC}")
message(STATUS "FFTW version: ${FFTW_VERSION}")
#===================#
#==== load FFTw ====#
#===================#

#==========================================================================
#=============== HYPAR_SFE sources pre-generation (using python and fypp)
#==========================================================================

#====================================
#======= python generation of .fpp
#====================================

include(CMake_get_python)

execute_process(
    COMMAND ${Python3_EXECUTABLE} ${HYPAR_SFE_DIR}/LIBS/PYTHON_SOURCES/python_precompile.py ${HYPAR_SFE_DIR}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "BUG in HYPAR_SFE/LIBS/CMakeLists.txt: python_precompile.py failed")
endif()

#========================================================================================
#======= Manual handling of dependencies on files generated through "execute_process"
#========================================================================================

if(is_ctest)
    file(GLOB_RECURSE commonFppFiles 
    CONFIGURE_DEPENDS "${HYPAR_SFE_DIR}/COMMON_SOURCES/*fpp" "${HYPAR_SFE_DIR}/PDE_ABSTRACT_SOURCES/*fpp" "${HYPAR_SFE_DIR}/PROBLEM_SOURCES/*inc"
    "${HYPAR_SFE_DIR}/CASES_DEVELOPPER/*F90" "${HYPAR_SFE_DIR}/CASES_DEVELOPPER/*fpp")
else()
    file(GLOB_RECURSE commonFppFiles 
    CONFIGURE_DEPENDS "${HYPAR_SFE_DIR}/COMMON_SOURCES/*fpp" "${HYPAR_SFE_DIR}/PDE_ABSTRACT_SOURCES/*fpp" "${HYPAR_SFE_DIR}/PROBLEM_SOURCES/*inc")
endif()
set_property(DIRECTORY APPEND PROPERTY
    CMAKE_CONFIGURE_DEPENDS ${commonFppFiles})
     

#===============================
#======= .fpp => .F90
#===============================

file(GLOB_RECURSE fppFiles CONFIGURE_DEPENDS ${HYPAR_SFE_DIR}/PROBLEM_SOURCES/*.fpp)

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

#==========================================================================
#=============== Building fortran HYPAR_SFE library
#==========================================================================

file(GLOB_RECURSE hypar_SFE_sources CONFIGURE_DEPENDS
${HYPAR_SFE_DIR}/COMMON_SOURCES/*.F90
${HYPAR_SFE_DIR}/PDE_ABSTRACT_SOURCES/*.F90
${HYPAR_SFE_DIR}/COMMON_SOURCES/POST_PROCESSING/*.c
${HYPAR_SFE_DIR}/PROBLEM_SOURCES/*.F90)
list(APPEND hypar_SFE_sources ${outFiles})

if(BUILD_ALL_DIMS)
    # Compiling HYPAR_SFE_FE_1D library
   set(DIMS 1 2)
else()
    # Compiling HYPAR_SFE_FE_2D library
   set(DIMS ${FE_dim})
endif()

foreach(dim IN LISTS DIMS)

    set(MY_HYPAR_SFE_LIB_OBJECT hypar_lib_F_${dim}D_obj)
    set(MY_HYPAR_SFE_LIB_STATIC hypar_lib_F_${dim}D_static)
    set(MY_HYPAR_SFE_LIB_DYNAMIC hypar_lib_F_${dim}D_dynamic)

    set(MY_HYPAR_SFE_MODULE_SOURCES ${HYPAR_SFE_DIR}/LIBS/BUILD/${dim}D/mod)

    message(STATUS "HYPAR_SFE_FE_${dim} sources: ${MY_HYPAR_SFE_MODULE_SOURCES}")

    # Define the three libraries
    set(LIBRARIES "obj" "static" "dynamic")
    add_library(${MY_HYPAR_SFE_LIB_OBJECT} OBJECT ${hypar_SFE_sources})
    add_library(${MY_HYPAR_SFE_LIB_STATIC} STATIC $<TARGET_OBJECTS:${MY_HYPAR_SFE_LIB_OBJECT}>)
    add_library(${MY_HYPAR_SFE_LIB_DYNAMIC} SHARED $<TARGET_OBJECTS:${MY_HYPAR_SFE_LIB_OBJECT}>)

    # Configuring libraries
    foreach(type_lib IN LISTS LIBRARIES)
        set(CUR_HYPAR_SFE_LIB hypar_lib_F_${dim}D_${type_lib})

        set_target_properties(${CUR_HYPAR_SFE_LIB} PROPERTIES
            Fortran_MODULE_DIRECTORY ${MY_HYPAR_SFE_MODULE_SOURCES}
            POSITION_INDEPENDENT_CODE ON)
        
        target_link_libraries(${CUR_HYPAR_SFE_LIB} PUBLIC
            PkgConfig::PETSC PkgConfig::FFTW
            ${MPI_Fortran_LIBRARIES}
            ${ZLIB_LIBRARIES}
            ${ADDITIONAL_LINKS})

        target_include_directories(${CUR_HYPAR_SFE_LIB} PUBLIC ${PETSC_INC} ${CMAKE_CURRENT_BINARY_DIR})

        target_compile_definitions(${CUR_HYPAR_SFE_LIB} PUBLIC 
        K_DIM_FE=${dim}
        PETSC_VERSION_MINOR=${PETSC_VERSION_MINOR}
        )

    endforeach()

endforeach()
