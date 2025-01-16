MESSAGE(STATUS "This is CMake ${CMAKE_VERSION}")
MESSAGE(STATUS "")

# Cmake minimum version
cmake_minimum_required(VERSION 3.5)

# Project name and version
project(test_mat)

# Enable Fortran
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 ")
    set(bounds "-Wall -fbounds-check -fimplicit-none -fcheck=all")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS_RELEASE " -O3 " CACHE STRING
  " Flags used by the compiler during release builds. ")
set(CMAKE_Fortran_FLAGS_NATIVE " -O4 -march=native -mtune=native -Ofast " CACHE STRING
  " Flags used by the compiler during release builds. ")

IF (TARGET debug)
else()
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
endif()
  
# Define runtime output directory
file(GLOB_RECURSE global_sources ${EULER_DIR}/SOURCES/COMMON_SOURCES/*.f90 ${EULER_DIR}/SOURCES/IDP_EULER/*.f90
  ${EULER_DIR}/SOURCES/SHOCKTUBE/*.f90 )

# add executable
set(exe "test.exe")
add_executable(${exe} ${CMAKE_CURRENT_SOURCE_DIR}/LOCAL_SOURCES/main.f90 ${global_sources} )

# find Intel MKL
set(BLA_VENDOR Intel10_64lp_seq)
find_package(LAPACK REQUIRED)

# include and link library
target_include_directories(${exe} PUBLIC $ENV{MKLROOT}/include)
target_link_libraries(${exe} ${SMALL_SFEMaNS_LIB} LAPACK::LAPACK)

