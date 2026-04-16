set(HYPAR_SFE_DIR "/home/victor.botez/HYPAR_SFE_VIRGIN/HYPAR_SFE")
set(ADDITIONAL_LINKS "-lmetis -lz -L /usr/lib/x86_64-linux-gnu/hdf5/serial")
set(debug_bounds "-Wall -fimplicit-none -fbounds-check")
set(release_bounds "-O4")
set(native_bounds "-march=native -mtune=native -Ofast")

set(RUN_PRE_PROC "srun")
set(PROC_CALL "-n")
set(RUN_POST_PROC "")
