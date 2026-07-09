# Variables for Whistler
#==== Mandatory for any application ====
set(HYPAR_SFE_DIR "my_path_to_hypar_sfe/HYPAR_SFE")
set(debug_bounds "-Wall -Wextra -fimplicit-none -finit-real=snan -finit-integer=-999999 -fbounds-check -ffree-line-length-none")
set(release_bounds "-O3 -ffree-line-length-none -flto")
set(FE_dim "1") #1 or 2
set(python_precompile OFF)
set(static_library OFF)
