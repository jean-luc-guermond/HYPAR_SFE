import os, sys
import shutil

#================================#
#======= PATH DEFINITIONS =======#
#================================#


def define_paths(dir_top):
    global dir_HY, dir_limiter
    dir_top = sys.argv[1]
    dir_HY = os.path.join(dir_top, "PDE_ABSTRACT_SOURCES/HYPERBOLIC_FE/")
    dir_limiter = os.path.join(dir_top, "COMMON_SOURCES/POST_PROCESSING")

#=========================================================#
#======= precompiling list of limiting functionals =======#
#=========================================================#

def precompile_limiting_func(n_lim, list_lim, path_lim, path_pb, name):
    #=== precompile limiter with cell = elements
    os.makedirs(os.path.join(path_pb, "TEMPLATED"), exist_ok=True)
    path_out = os.path.join(path_pb, f"TEMPLATED/limiting_{name}.fpp")

    lines_include = ""
    for i in range(n_lim):
        path_new_lim = os.path.join(path_lim, f'limiting_{list_lim[i]}.inc')
        lines_include = "\n".join([f"#:include '{path_new_lim}'", 
                                   "",
                                lines_include])
    lines = "\n".join([f"MODULE limiting_functionals_{name}_module",
                       "   CONTAINS",
                       lines_include,
                       f"END MODULE limiting_functionals_{name}_module"])

    with open(path_out,'w') as f:
        f.write(lines)

#=============================================================================================================================#
#======= precompiling limiting with cell averaging over elements, using uk_plus & uk_minus (for syst_dim and inlining) =======#
#=============================================================================================================================#

def precompile_limiter_cell_elt(syst_dim, n_lim, list_lim, path_lim, path_pb, name):
    #=== precompile limiter with cell = elements
    path_cell_elts = os.path.join(dir_limiter, "limiter_cell_elt.fpp")
    os.makedirs(os.path.join(path_pb, "TEMPLATED"), exist_ok=True)
    path_out = os.path.join(path_pb, f"TEMPLATED/limiter_cell_elt_{name}.fpp")
    shutil.copy(path_cell_elts, path_out)
    
    with open(path_out,'r') as f:
        lines = f.read()
    lines = lines.replace("list_limiters_input", str(list_lim))
    lines = lines.replace("syst_dim_input", str(syst_dim))
    lines = lines.replace("n_lim_input", str(n_lim))
    lines = lines.replace("template", f"{name}")

    lines = "\n".join(["#:enddef", lines])
    for i in range(n_lim):
        path_new_lim = os.path.join(path_lim, f'limiting_{list_lim[i]}.inc')
        lines = "\n".join([f"#:if name=='{list_lim[i]}'", 
                           f"#:include '{path_new_lim}'", 
                           "#:endif", 
                           lines])
    lines = "\n".join(["#:def limiter_func(name)", lines])

    with open(path_out,'w') as f:
        f.write(lines)

    #=== precompile limiting class (FIXME?)
    path_limiting_class = os.path.join(dir_limiter, "cell_limiting_engine_parallel.fpp")
    os.makedirs(os.path.join(path_pb, "TEMPLATED"), exist_ok=True)
    path_out = os.path.join(path_pb, f"TEMPLATED/cell_limiting_engine_parallel_{name}.fpp")
    shutil.copy(path_limiting_class, path_out)
    
    with open(path_out,'r') as f:
        lines = f.read()
    lines = lines.replace("syst_dim_input", str(syst_dim))
    lines = lines.replace("template", f"{name}")

    with open(path_out,'w') as f:
        f.write(lines)

#==========================================================================#
#======= precompiling bounds computation with uijbar (for inlining) =======#
#==========================================================================#

def precompile_uijbar_bounds(syst_dim, n_lim, list_lim, path_lim, path_pb, name):
    #=== precompile uij_bar bounds computation
    path_uijbar_bounds = os.path.join(dir_HY, "uij_bar_bounds.fpp")
    os.makedirs(os.path.join(path_pb, "TEMPLATED"), exist_ok=True)
    path_out = os.path.join(path_pb, f"TEMPLATED/uij_bar_bounds_{name}.fpp")
    shutil.copy(path_uijbar_bounds, path_out)
    
    with open(path_out,'r') as f:
        lines = f.read()
    lines = lines.replace("list_limiters_input", str(list_lim))
    lines = lines.replace("syst_dim_input", str(syst_dim))
    lines = lines.replace("n_lim_input", str(n_lim))
    lines = lines.replace("template", f"{name}")

    lines = "\n".join(["#:enddef", lines])
    for i in range(n_lim):
        path_new_lim = os.path.join(path_lim, f'limiting_{list_lim[i]}.inc')
        lines = "\n".join([f"#:if name=='{list_lim[i]}'", 
                           f"#:include '{path_new_lim}'", 
                           "#:endif", 
                           lines])
    lines = "\n".join(["#:def limiter_func(name)", lines])

    with open(path_out,'w') as f:
        f.write(lines)

#=============================================================================#
#======= precompiling abstract_hyperbolic.F90 (for syst_dim and n_lim) =======#
#=============================================================================#

def precompile_abstract_hyperbolic_module(syst_dim, n_lim, list_lim, path_lim, path_pb, name):
    #=== precompile limiter with cell = elements
    path_cell_elts = os.path.join(dir_HY, "abstract_hyperbolic.fpp")
    os.makedirs(os.path.join(path_pb, "TEMPLATED"), exist_ok=True)
    path_out = os.path.join(path_pb, f"TEMPLATED/abstract_hyperbolic_{name}.fpp")
    shutil.copy(path_cell_elts, path_out)
    
    with open(path_out,'r') as f:
        lines = f.read()
    lines = lines.replace("list_limiters_input", str(list_lim))
    lines = lines.replace("syst_dim_input", str(syst_dim))
    lines = lines.replace("n_lim_input", str(n_lim))
    lines = lines.replace("template", f"{name}")

    lines = "\n".join(["#:enddef", lines])
    for i in range(n_lim):
        path_new_lim = os.path.join(path_lim, f'limiting_{list_lim[i]}.inc')
        lines = "\n".join([f"#:if name=='{list_lim[i]}'", 
                           f"#:include '{path_new_lim}'", 
                           "#:endif", 
                           lines])
    lines = "\n".join(["#:def limiter_func(name)", lines])

    with open(path_out,'w') as f:
        f.write(lines)

#===============================================#
#======= precompiling hyperbolic problem =======#
#===============================================#

def precompile_hyperbolic(syst_dim, n_lim, list_lim, path_lim, path_pb, name):

    precompile_limiting_func(n_lim, list_lim, path_lim, path_pb, name)
    precompile_limiter_cell_elt(syst_dim, n_lim, list_lim, path_lim, path_pb, name)
    precompile_uijbar_bounds(syst_dim, n_lim, list_lim, path_lim, path_pb, name)
    precompile_abstract_hyperbolic_module(syst_dim, n_lim, list_lim, path_lim, path_pb, name)
