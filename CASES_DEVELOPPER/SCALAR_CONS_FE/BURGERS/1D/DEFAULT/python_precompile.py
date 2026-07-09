import os, sys
from pathlib import Path
dir_top = "my_path_to_hypar_sfe/HYPAR_SFE"
sys.path.append(str(Path(dir_top) / "LIBS/PYTHON_SOURCES"))
from precompile_sources import define_paths, precompile_hyperbolic

define_paths(dir_top)
path_pb = os.getcwd()

#=== Definitions for BURGERS problem in 1D
syst_dim = 1
n_lim = 2
list_lim = ['scalar_rho_min', 'scalar_rho_max']
path_lim = str(Path("path_to_case") / Path("LOCAL_SOURCES"))
path_pb = str(Path("path_to_case") / Path("LOCAL_SOURCES"))
name = 'burgers'

precompile_hyperbolic(syst_dim, n_lim, list_lim, path_lim, path_pb, name)
