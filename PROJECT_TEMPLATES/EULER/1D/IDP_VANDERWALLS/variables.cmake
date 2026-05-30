# Variables for Whistler
#==== Mandatory for any application ====
set(HYPAR_SFE_DIR "/home/botez/HYPAR_SFE")
set(debug_bounds "-Wall -Wextra -fimplicit-none -finit-real=snan -finit-integer=-999999 -fbounds-check")
#set(debug_bounds "-Wall -fimplicit-none -fbounds-check")
set(release_bounds "-O3")
set(FE_dim "1") #1 or 2
set(static_library ON)

#==== Specific to CTEST ====
set(RUN_PRE_PROC "mpirun")
set(PROC_CALL "-np ")
set(RUN_POST_PROC "")
set(PARALLEL_TEST_LEVEL 0)
