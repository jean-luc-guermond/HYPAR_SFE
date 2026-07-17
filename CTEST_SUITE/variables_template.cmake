#==== Mandatory for any application ====
# set(HYPAR_SFE_DIR "my_path_to_hypar/HYPAR_SFE")

#===compile bounds for GNU compiler
set(debug_bounds "-Wall -fbounds-check -fimplicit-none -O0 -ffree-line-length-none")
set(release_bounds "-O3 -ffree-line-length-none")

#===compile bounds for Intel compiler (incoming)
# set(debug_bounds "-O0 -g -traceback -heap-arrays -check bounds -warn all")
# set(release_bounds "-O2") <== FIXME (-O3 fails on Stokes)

#==== Specific to CTEST ====
#set(RUN_PRE_PROC "mpirun")
#set(PROC_CALL "-np ")
#set(RUN_POST_PROC "")
set(RUN_PRE_PROC "srun")
set(PROC_CALL "-n")
set(RUN_POST_PROC "")

#==== Developper settings ====
# How many MPI configurations, and what upper bound for number of CPUS?
set(NB_PARALLEL_TEST 5)
set(MAX_PROC_MPI_TEST 16)
# Do we use static library? (STRONGLY UNRECOMMENDED FOR CTEST AS THERE CAN BE HUNDREDS OF HEAVY EXECUTABLES)
set(static_library OFF)

# GNU compiler debug bounds:
# set(debug_bounds "-Wall \
# -fbounds-check \
# -fimplicit-none \
# -finit-real=snan \
# -finit-integer=-999999 \
# -Wextra \
# -Wconversion \
# -fcheck=all \
# -ffpe-trap=invalid,zero,overflow \
# -fbacktrace \
# -fno-omit-frame-pointer \
# -g3 \
# -O0")