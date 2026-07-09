import os, sys
from pathlib import Path
dir_top = "my_path_to_hypar_sfe/HYPAR_SFE"
sys.path.append(str(Path(dir_top) / "LIBS/PYTHON_SOURCES"))
from precompile_sources import define_paths, precompile_hyperbolic

define_paths(dir_top)
path_pb = os.getcwd()

#=== Some hyperbolic problem (uncomment and fill as many blocks as needed)
# syst_dim = 1
# n_lim = 2
# list_lim = ['scalar_rho_min', 'scalar_rho_max']
# path_lim = os.path.join(dir_top, "relative_path_to_limiting_functionals") #this path must contain .inc files
# path_pb = os.path.join(dir_top, "relative_path_to_pb") #a TEMPLATED folder containing the generated sources will be created here
# name = 'name_of_object'

# precompile_hyperbolic(syst_dim, n_lim, list_lim, path_lim, path_pb, name)
