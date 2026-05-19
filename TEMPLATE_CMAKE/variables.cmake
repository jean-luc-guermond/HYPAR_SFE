# Variables for Whistler
#==== Mandatory for any application ====
set(HYPAR_SFE_DIR "/home/guermond/HYPAR_SFE")
set(debug_bounds "-Wall -fimplicit-none -fbounds-check")
set(release_bounds "-O3")
set(FE_dim "1") #1 or 2

#==== Specific to CTEST ====
set(RUN_PRE_PROC "mpirun")
set(PROC_CALL "-np ")
set(RUN_POST_PROC "")
set(PARALLEL_TEST_LEVEL 0)
