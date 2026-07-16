import sys
from pathlib import Path

#================================#
#======= PATH DEFINITIONS =======#
#================================#

def overwrite_file(lines, file_out, opt_message="file"):
    if (Path(file_out).is_file()):
        with open(file_out,'r') as f:
            old_lines = f.read()

        if (old_lines != lines):
            print("file "+opt_message+" was modified, rebuilding...")
            files_differ = True
        else:
            files_differ = False
    else:
        print("Building "+opt_message+" for first time...")
        files_differ = True

    if (files_differ):
        with open(file_out,'w') as f:
            f.write(lines)

def define_paths(dir_top):
    global dir_HY, dir_limiter
    dir_top = sys.argv[1]
    dir_HY = str(Path(dir_top) / "PDE_ABSTRACT_SOURCES" / "HYPERBOLIC_FE")
    dir_limiter = str(Path(dir_top) / 'COMMON_SOURCES' / 'POST_PROCESSING')

#=========================================================#
#======= precompiling list of limiting functionals =======#
#=========================================================#

def precompile_limiting_func(n_lim, list_lim, path_lim, path_pb, name):
    #=== precompile limiter with cell = elements
    (Path(path_pb) / 'TEMPLATED').mkdir(parents=True, exist_ok=True)
    path_out = Path(path_pb) / "TEMPLATED" / f"limiting_{name}.fpp"

    lines_include = ""
    for i in range(n_lim):
        path_new_lim = str(Path(path_lim) / f"limiting_{list_lim[i]}.inc")
        lines_include = "\n".join([f"#:include '{path_new_lim}'", 
                                   "",
                                lines_include])
    lines = "\n".join([f"MODULE limiting_functionals_{name}_module",
                       "   CONTAINS",
                       lines_include,
                       f"END MODULE limiting_functionals_{name}_module"])

    overwrite_file(lines, path_out, opt_message=f"limiting_{name}.fpp")

#=============================================================================================================================#
#======= precompiling limiting with cell averaging over elements, using uk_plus & uk_minus (for syst_dim and inlining) =======#
#=============================================================================================================================#

def precompile_limiter_cell_elt(syst_dim, n_lim, list_lim, path_lim, path_pb, name):
    #=== precompile limiter with cell = elements
    path_cell_elts = Path(dir_limiter) / "limiter_cell_elt.fpp"
    (Path(path_pb) / "TEMPLATED").mkdir(parents=True, exist_ok=True)
    path_out = Path(path_pb) / "TEMPLATED" / f"limiter_cell_elt_{name}.fpp"
    
    with open(path_cell_elts,'r') as f:
        lines = f.read()
    lines = lines.replace("list_limiters_input", str(list_lim))
    lines = lines.replace("syst_dim_input", str(syst_dim))
    lines = lines.replace("n_lim_input", str(n_lim))
    lines = lines.replace("template", f"{name}")

    lines = "\n".join(["#:enddef", lines])
    for i in range(n_lim):
        path_new_lim = Path(path_lim) / f'limiting_{list_lim[i]}.inc'
        lines = "\n".join([f"#:if name=='{list_lim[i]}'", 
                           f"#:include '{path_new_lim}'", 
                           "#:endif", 
                           lines])
    lines = "\n".join(["#:def limiter_func(name)", lines])

    overwrite_file(lines, path_out, opt_message="limiter_cell_elt")

    #=== precompile limiting class (FIXME?)
    path_limiting_class = Path(dir_limiter) / "cell_limiting_engine_parallel.fpp"
    (Path(path_pb) / "TEMPLATED").mkdir(parents=True, exist_ok=True)
    path_out = Path(path_pb) / "TEMPLATED" / f"cell_limiting_engine_parallel_{name}.fpp"
    
    with open(path_limiting_class,'r') as f:
        lines = f.read()
    lines = lines.replace("syst_dim_input", str(syst_dim))
    lines = lines.replace("template", f"{name}")

    overwrite_file(lines, path_out, opt_message="cell_limiting_engine_parallel")

#==========================================================================#
#======= precompiling bounds computation with uijbar (for inlining) =======#
#==========================================================================#

def precompile_uijbar_bounds(syst_dim, n_lim, list_lim, path_lim, path_pb, name):
    #=== precompile uij_bar bounds computation
    path_uijbar_bounds = Path(dir_HY) / "uij_bar_bounds.fpp"
    (Path(path_pb) / "TEMPLATED").mkdir(parents=True, exist_ok=True)
    path_out = Path(path_pb) / "TEMPLATED" / f"uij_bar_bounds_{name}.fpp"
    
    with open(path_uijbar_bounds,'r') as f:
        lines = f.read()
    lines = lines.replace("list_limiters_input", str(list_lim))
    lines = lines.replace("syst_dim_input", str(syst_dim))
    lines = lines.replace("n_lim_input", str(n_lim))
    lines = lines.replace("template", f"{name}")

    lines = "\n".join(["#:enddef", lines])
    for i in range(n_lim):
        path_new_lim = Path(path_lim) / f'limiting_{list_lim[i]}.inc'
        lines = "\n".join([f"#:if name=='{list_lim[i]}'", 
                           f"#:include '{path_new_lim}'", 
                           "#:endif", 
                           lines])
    lines = "\n".join(["#:def limiter_func(name)", lines])

    overwrite_file(lines, path_out, opt_message="inside uij_bar_bounds")

#=============================================================================#
#======= precompiling abstract_hyperbolic.F90 (for syst_dim and n_lim) =======#
#=============================================================================#

def precompile_abstract_hyperbolic_module(syst_dim, n_lim, list_lim, path_lim, path_pb, name):
    #=== precompile limiter with cell = elements
    path_abstract_hyperbolic = Path(dir_HY) / "abstract_hyperbolic.fpp"
    (Path(path_pb) / "TEMPLATED").mkdir(parents=True, exist_ok=True)
    path_out = Path(path_pb) / "TEMPLATED" / f"abstract_hyperbolic_{name}.fpp"
    
    with open(path_abstract_hyperbolic,'r') as f:
        lines = f.read()
    lines = lines.replace("list_limiters_input", str(list_lim))
    lines = lines.replace("syst_dim_input", str(syst_dim))
    lines = lines.replace("n_lim_input", str(n_lim))
    lines = lines.replace("template", f"{name}")

    lines = "\n".join(["#:enddef", lines])
    for i in range(n_lim):
        path_new_lim = Path(path_lim) / f"limiting_{list_lim[i]}.inc"
        lines = "\n".join([f"#:if name=='{list_lim[i]}'", 
                           f"#:include '{path_new_lim}'", 
                           "#:endif", 
                           lines])
    lines = "\n".join(["#:def limiter_func(name)", lines])

    overwrite_file(lines, path_out, opt_message="abstract_hyperbolic")

#===============================================#
#======= precompiling hyperbolic problem =======#
#===============================================#

def precompile_hyperbolic(syst_dim, n_lim, list_lim, path_lim, path_pb, name):

    precompile_limiting_func(n_lim, list_lim, path_lim, path_pb, name)
    precompile_limiter_cell_elt(syst_dim, n_lim, list_lim, path_lim, path_pb, name)
    precompile_uijbar_bounds(syst_dim, n_lim, list_lim, path_lim, path_pb, name)
    precompile_abstract_hyperbolic_module(syst_dim, n_lim, list_lim, path_lim, path_pb, name)
