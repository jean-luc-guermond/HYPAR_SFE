import os, sys
from precompile_sources import *

dir_top = sys.argv[1]
define_paths(dir_top)
path_pb = os.getcwd()

#=== Linear transport
syst_dim = 1
n_lim = 2
list_lim = ['scalar_rho_min', 'scalar_rho_max']
path_lim = os.path.join(dir_top, "PROBLEM_SOURCES/SCALAR_CONSERVATION_FE/LINEAR_TRANSPORT")
path_pb = os.path.join(dir_top, "PROBLEM_SOURCES/SCALAR_CONSERVATION_FE/LINEAR_TRANSPORT")
name = 'linear_transport'

precompile_hyperbolic(syst_dim, n_lim, list_lim, path_lim, path_pb, name)

#=== Euler
syst_dim = "'k_dim+2'"
n_lim = 2
list_lim = ['euler_rho_min', 'euler_rho_max']
path_lim = os.path.join(dir_top, "PROBLEM_SOURCES/EULER")
path_pb = os.path.join(dir_top, "PROBLEM_SOURCES/EULER")
name = 'euler'

precompile_hyperbolic(syst_dim, n_lim, list_lim, path_lim, path_pb, name)

#=== Rk: precompilation of Navier-Stokes relies on precompilation of Euler only, for now Stokes isn't precompiled