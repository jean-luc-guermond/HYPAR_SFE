import shutil
from pathlib import Path
import sys, os

path_hypar = sys.argv[1]
fe_dim = sys.argv[2]
name_developper_case = sys.argv[3]

def copy_if_files_differ(file_1, file_2):
    with open(file_1, 'r') as f:
        lines_1 = f.read()
    if Path(file_2).is_file():
        with open(file_2, 'r') as f:
            lines_2 = f.read()
        files_differ = lines_1 != lines_2
    else:
        files_differ = True
    if files_differ:
        shutil.copy(file_1, file_2)
    return files_differ


def generate_user_case(path_hypar, name_developper_case, fe_dim, exist_ok=True):
    
    path_user_case = Path(path_hypar) / "CASES_USER" / name_developper_case
    path_user_case.mkdir(parents=True, exist_ok=True)

    path_dvp_case = Path(path_hypar) / "CASES_DEVELOPPER" / name_developper_case
    list_spe_case = [file.name for file in Path(path_dvp_case).iterdir() if file.name != 'TEMPLATE']

    list_already_exists = []
    for spe_case in list_spe_case:
        path_spe_case = path_user_case / spe_case
        try:
            (path_spe_case / "LOCAL_SOURCES").mkdir(parents=True, exist_ok=exist_ok)
            (path_spe_case / "BUILD").mkdir(parents=True, exist_ok=exist_ok)

            #=== Copy local_sources (*.F90, *.inc)===#
            for file_F90 in (path_dvp_case / spe_case / "LOCAL_SOURCES").glob('*.F90'):
                _ = copy_if_files_differ(file_F90, path_user_case / spe_case / "LOCAL_SOURCES" / file_F90.name)
            for file_inc in (path_dvp_case / spe_case / "LOCAL_SOURCES").glob('*.inc'):
                _ = copy_if_files_differ(file_inc, path_user_case / spe_case / "LOCAL_SOURCES" / file_inc.name)
            #=== Copy template_sources ===#
            for file_F90 in (path_dvp_case / "TEMPLATE").glob("*.F90"):
                _ = copy_if_files_differ(file_F90, path_user_case / spe_case / "LOCAL_SOURCES" / file_F90.name)

            #=== copy cmake files from template_cmake ===#
            path_in = Path(path_hypar) / "TEMPLATE_CMAKE" / "CMakeLists.txt"
            path_out = path_user_case / spe_case / "CMakeLists.txt"
            _ = copy_if_files_differ(path_in, path_out)

            path_in = Path(path_hypar) / "TEMPLATE_CMAKE" / "variables.cmake"
            path_out_var_cmake = path_user_case / spe_case / "variables.cmake"
            var_cmake_differ = copy_if_files_differ(path_in, path_out_var_cmake)

            # path_in = Path(path_hypar) / "TEMPLATE_CMAKE" / "python_precompile.py"
            path_py_precompile_in = path_dvp_case / spe_case / "python_precompile.py"
            path_py_precompile_out = path_user_case / spe_case / "python_precompile.py"
            if (os.path.exists(path_py_precompile_in)):
                py_precompile_differ = copy_if_files_differ(path_py_precompile_in, path_py_precompile_out)

            #=== Configure variables.cmake ===#
            if (var_cmake_differ):
                with open(path_out_var_cmake, 'r') as f:
                    lines = f.read()
                lines = lines.replace("my_path_to_hypar_sfe/HYPAR_SFE", path_hypar)
                lines = lines.replace("my_FE_dim", str(fe_dim))

                with open(path_out_var_cmake, 'w') as f:
                    f.write(lines)
            
            #=== Configure python_precompile.py ===#
            if (os.path.exists(path_py_precompile_in)):
                if (py_precompile_differ):
                    with open(path_py_precompile_out, 'r') as f:
                        lines = f.read()
                    lines = lines.replace("my_path_to_hypar_sfe/HYPAR_SFE", path_hypar)
                    lines = lines.replace("path_to_case", str(path_spe_case))

                    with open(path_py_precompile_out, 'w') as f:
                        f.write(lines)        
        
        except FileExistsError:
            list_already_exists.append(spe_case)
    # print("Already existing, therefore not overwritten: "+ ", ".join(list_already_exists))
generate_user_case(path_hypar, name_developper_case, fe_dim)