Getting started: installation instructions and CTEST
===================================

List of mandatory dependencies
--------------------

0. **cmake**: must be at least 3.29

1. **PETSc**: version must be at least 3.21 (3.25 strongly encouraged), and installation made by the user

> wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.25.3.tar.gz & tar -xvf petsc-3.25.3.tar.gz & cd petsc-3.25.3

> export PETSC_DIR=$PWD

> export PETSC_ARCH=linux-gnu-c (or whatever other name)

> ./configure --configModules=PETSc.Configure --optionsModule=config.compilerOptions --download-fblaslapack=1 --with-shared-libraries=1 --download-hypre=1 --download-mumps=1 --download-scalapack=1 --download-metis=1 --download-parmetis=1 --download-blas=1 --with-debugging=0 --with-x=0

> make all

> make check

2. **fftw3**: see https://www.fftw.org/download.html, installation made by the user

3. **Python3 and fypp**: fypp is taken care of manually by the CTEST compilation, so the user should only care about the installation of python3

compiling and running the ctests:
---------------------------------

This step serves two purposes:

### Compilation

> cd HYPAR_SFE/CTEST

> mkdir BUILD

> cd BUILD

> cmake ..

> make

> ctest -LE mpi -L quick

More information on ctest options is given below. Also you might need to modify HYPAR_SFE/CTEST/variables.cmake especially if you do not use srun or need to specify additional options. 

The *cmake ..* step also makes sure all dependencies are installed with compatible versions. 

### Generation of folder **CASES_USER**

Compiling CTEST will also generate a bunch of applications ready to copy/paste within HYPAR_SFE/CASES_USER. The generation of those ready to compile applications is based on:

1. HYPAR_SFE/TEMPLATE_CMAKE (for cmake sources)

2. HYPAR_SFE/CASES_DEVELOPPER (for Fortran sources)

**IMPORTANT**: do not modify directly the sources from CASES_USER. They are solely intended for copy and pasting and are overwritten when building again CTEST. 

        
Optional dependencies:
----------------------------

### Post-processing 
=> HYPAR_SFE provides subroutines which makes visualization handy with plotmtv and paraview
### Mesh Generator  
=> download through "wget https://people.tamu.edu/~guermond/DOWNLOADS/MESH_GENERATOR_SFEMaNS.tar.gz". 
                      documentation on: https://people.tamu.edu/~guermond//SFEMaNS/html/doc_mesh_generator.html

Building one's own application
==============================

Compilation
--------------

After having ran the CTEST, the user can basically just copy/paste whatever application generated within the **CASES_USER** folder. For example

        cp -r HYPAR_SFE/CASES_USER/EULER_FE/2D/IDP_WIND_TUNNEL MY_APPLICATIONS_HYPAR_SFE/IDP_WIND_TUNNEL
        cd MY_APPLICATIONS_HYPAR_SFE/IDP_WIND_TUNNEL

The application will already contain the cmake requirements (variables.cmake and CMakeLists.txt) which are already correctly filled if the ctests were passed successfully, as well as the LOCAL_SOURCES/*F90 which are required for compiling the specific application. Exception made for the applications which contain a "python_precompile.py" file, the only things left to do is:

        mkdir BUILD
        cd BUILD
        cmake ..
        make 

Now the executable can be found in MY_APPLICATIONS_HYPAR_SFE/IDP_WIND_TUNNEL/EXECUTABLE and be run directly. The data file to be filled by the user will be generated automatically everytime the executable is run, until the data file contains all information.

For more information regarding "python_precompile.py" see [this link.](./DOCS/my_hyperbolic_class.md)

Adding one's own arguments to the data file
-------------------------------------------

A built-in user data framework is provided with HYPAR_SFE. Further documentation can be found on [this link.](./DOCS/data_reading.md)


Defining one's own hyperbolic conservation equation
---------------------------------------------------

Currently hyperbolic conservation equations are codded through Fortran abstract classes and are compiled thanks to python/fypp code generation. This combination allows the code to be easily maintainable as well as optimized for every given application through hard-coding dimensions and inlining. Documentation on how to program one's own hyperbolic class is provided on [this link.](./DOCS/my_hyperbolic_class.md)


# Developping a new CTEST in HYPAR_SFE

A whole cmake environment is provided for coding ctests. More information can be found on [this link](./DOCS/ctest.md)



TODO LIST & General info
========================

Clusters & compilation information
----------------------------------

1. JeanZay & TGCC => GCC: internal compiler error, submit a full report

2. JeanZay => intel environment available with petsc/3.24

3. TGCC => intel environment possible with petsc/3.15 or petsc/3.17, ask to install petsc/3.25

CTEST currently failing
-----------------------

1. Periodic mesh partitionning after some point (PERIODIC_2D_SCAL_ADV after ~ 14 procs)

2. Wind tunnel with many processes (bug on udotn)

3. Vanderwalls test cases > 1

4. GNU Ctest with MINMOD with debug bound -finit-real=snan -finit-integer=-999999 -ffpe-trap=invalid,zero,overflow

5. Intel On JeanZay: compilation -O3, ctest -L stokes

CTEST to add
------------

1. Euler

- Test limiting with Periodicity

2. Navier-Stokes

- Periodic conditions

3. Convergence tests on Poisson:

- 1D mesh refinements

- 1D data refinements

- 2D with refinement factor

4. Convergence tests on Euler

- Convergence of Sodt & Smooth, P1 & Pk?

- Convergence of Isentropic vortex, P2, P3

5. Convergence tests on Linear transport

- Check already existing ones

6. Convergence tests on Navier-Stokes

- Convergence of Becker test

7. Throughput of Hypre & MUMPS (see already existing Elasticity)

- Against refinements

- Against MPI procs

8. Benchmarks

- Wind Tunnel

- CTVD

- Noh

- Kelvin-Helmholtz


Dependencies
-------------

- modify MESH_GENERATOR user API

- complete doc for fftw

HYPAR_SFE functionalities
-------------------------

1. Entropy/Energy? limiting

2. Limiting method using u_L

3. Bernstein polynomials

4. Euler Fourier with padding

5. Abstract class for parabolic problems?

6. Factorize data in ctest

7. change commutator so that the code runs even though eta_commute might be zero at some points

8. By default, associate the uijbar and limiting calculations outside of the start_setup.F90?

9. Subroutines to Factorize? iterative_LA, st_aij_csr_glob_block_with_extra_layer, reorganize SOLVERS, MESH ?

10. Delete prep_periodic_bloc?