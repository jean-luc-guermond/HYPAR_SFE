# Variables for personal laptop
#==== Mandatory for any application ====
set(HYPAR_SFE_DIR "/home/victor.botez/HYPAR_SFE_VIRGIN/HYPAR_SFE")
set(debug_bounds "$ENV(FFLAGS) -O2 -g -traceback -heap-arrays")
set(release_bounds "-O3")
set(FE_dim "1") #1 or 2

#==== Specific to CTEST ====
set(RUN_PRE_PROC "mpirun")
set(PROC_CALL "-np ")
set(RUN_POST_PROC "")



#============ SETTING VARIABLES FOR COMPILATION  =============
function(define_compile_mode)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug_bounds}")
    set(CMAKE_Fortran_FLAGS_RELEASE ${release_bounds} CACHE STRING
        " Flags used by the compiler during release builds. ")
    set(CMAKE_Fortran_FLAGS_NATIVE ${native_bounds} CACHE STRING
            " Flags used by the compiler during native builds. ")

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
endfunction()
