# Project name and version
project(hypar_lib)
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
#=============== compilation of HYPAR_SFE
#==========================================================================

# Getting mandatory libraries
find_package(MPI REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(Python3 REQUIRED)

pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)
pkg_check_modules(FFTW REQUIRED IMPORTED_TARGET fftw3)
# Getting mandatory libraries

# Debugging information for libraries
message(STATUS "Working in: ${HYPAR_SFE_DIR}")
message(STATUS "ZLIB_LIBRARIES=${ZLIB_LIBRARIES}")
message(STATUS "MPI_Fortran_LIBRARIES=${MPI_Fortran_LIBRARIES}")
get_target_property(PETSC_INC PkgConfig::PETSC INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "PETSC includes: ${PETSC_INC}")
get_target_property(FFTW_INC PkgConfig::FFTW INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "FFTW includes: ${FFTW_INC}")
# Debugging information for libraries


execute_process(
    COMMAND ${Python3_EXECUTABLE} ${HYPAR_SFE_DIR}/LIBS/python_precompile.py ${HYPAR_SFE_DIR}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE result
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "python_precompile.py failed")
endif()

file(GLOB_RECURSE fppFiles CONFIGURE_DEPENDS ${HYPAR_SFE_DIR}/PROBLEM_SOURCES/*.fpp)

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

# Building HYPAR_SFE library
# file(GLOB_RECURSE hypar_SFE_sources CONFIGURE_DEPENDS
file(GLOB_RECURSE hypar_SFE_sources
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

        target_compile_definitions(${CUR_HYPAR_SFE_LIB} PUBLIC K_DIM_FE=${dim})
        
        set_target_properties(${CUR_HYPAR_SFE_LIB} PROPERTIES
            Fortran_MODULE_DIRECTORY ${MY_HYPAR_SFE_MODULE_SOURCES}
            POSITION_INDEPENDENT_CODE ON)
        
        target_link_libraries(${CUR_HYPAR_SFE_LIB} PUBLIC
            PkgConfig::PETSC PkgConfig::FFTW
            ${MPI_Fortran_LIBRARIES}
            ${ZLIB_LIBRARIES}
            ${ADDITIONAL_LINKS})

        target_include_directories(${CUR_HYPAR_SFE_LIB} PUBLIC ${PETSC_INC})

    endforeach()

endforeach()
