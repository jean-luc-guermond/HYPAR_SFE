MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Includes various variables depending on native
include(variables.cmake)

# Cmake minimum version
cmake_minimum_required(VERSION 3.29)

# Project name and version
project(sfemans)

# Compile C routines and create library
enable_language(C)

#file(GLOB_RECURSE sfemans_sources_C
#  ${SFEMaNS_DIR}/SOURCES/WRITE_BLOCK_MM/write_block_c.c)

#add_library(sfemans_lib_C SHARED ${sfemans_sources_C})


# Enable Fortran
enable_language(Fortran)

# Getting useful libraries
find_package(MPI)
find_package(ZLIB)
find_package(PkgConfig)

pkg_check_modules(PETSC PETSc)
pkg_check_modules(FFTW fftw3)

# Searching all petsc linked libraries
execute_process (COMMAND pkg-config PETSc --libs --static OUTPUT_VARIABLE PETSC_PRIVATE_LIBRARIES)
string(STRIP "${PETSC_PRIVATE_LIBRARIES}" PETSC_PRIVATE_LIBRARIES)

# Setting for different builds
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug_bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE ${release_bounds} CACHE STRING
  " Flags used by the compiler during release builds. ")
set(CMAKE_Fortran_FLAGS_NATIVE ${native_bounds} CACHE STRING
  " Flags used by the compiler during native builds. ")

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(native
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Native ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Native"
  )

 # Create local sources from local dir
file(GLOB_RECURSE sfemans_sources
	${SFEMaNS_DIR}/COMMON_SOURCES/*.F90
	${SFEMaNS_DIR}/PROJECT_TEMPLATES/EXAMPLE_MESH_CREATION/*.F90
)
# Define runtime output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../EXECUTABLE)

message(${sfemans_sources})

# Add executable
add_executable(a.exe ${sfemans_sources})

# Include and link library
target_link_libraries(a.exe PUBLIC
  ${PETSC_LINK_LIBRARIES}
  ${PETSC_PRIVATE_LIBRARIES}
  ${FFTW_LINK_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${MPI_Fortran_LIBRARIES}
  ${ADDITIONAL_LINKS}
  )
target_include_directories(a.exe PUBLIC
    ${FFTW_INCLUDEDIR} ${FFTW_INCLUDE_DIRS}
    ${PETSC_INCLUDEDIR} ${PETSC_INCLUDE_DIRS}
  )
#target_include_directories(a.exe )



# Adding directory to build lib_f

