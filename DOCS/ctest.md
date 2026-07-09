The CTEST tree in HYPAR_SFE is based on several layers of CMakeLists.txt which will successively call each other. It starts with

CTEST/CMakeLists.txt
--------------------

Initializes cmake and calls successively Libs and Layer I


LIBS/CMakeLists.txt
-------------------

Builds the HYPAR_SFE libraries, 1D and 2D, and takes care of the automatic generation of .fpp files for all problems implemented inside PROBLEM_SOURCES. This file is called for any application, except a given application specifies the problem dimension in order to compile just the one of interest.


Then there are three cmake layers within REGRESSION_SUITE:

Layer I: REGRESSION_SUITE/CMakeLists.txt
----------------------------------------

This layer is just a succession of add_subdirectory(...) meant to call Layer II.

Layer II: e.g REGRESSION_SUITE/EULER_FE/CMakeLists.txt
-----------------------------------------------------

This layer specializes on a given case problem, e.g EULER_FE, NAVIER_STOKES_COMPRESSIBLE_FE, etc... CASES_USER are generated at that point through the command:

    create_cases_user(
        SOURCES_DIR EULER_FE/2D
        FE_DIM      2
    )

This call generates all user cases it can based on what it finds within CASES_DEVELOPPER/EULER_FE/2D. In this case, it will detect (July 9th 2026):

> CASES_DEVELOPPER/EULER_FE/2D/IDP_WIND_TUNNEL

> CASES_DEVELOPPER/EULER_FE/2D/IDP_ISENTROPIC_VORTEX

> CASES_DEVELOPPER/EULER_FE/2D/TEMPLATE

Therefore, based on the common sources present within TEMPLATE, the code will create the corresponding folders:

> CASES_USER/EULER_FE/2D/IDP_WIND_TUNNEL/LOCAL_SOURCES

> CASES_USER/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/LOCAL_SOURCES

Layer II then calls Layer III through the new add_subdirectory.

Layer III: e.g REGRESSION_SUITE/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/CMakeLists.txt
-----------------------------------------------------------------

This layer specializes in a given physical problem to a given case/boundary conditions. The whole point of layer III is to be able to compile the sources located at the corresponding CASES_USER location (e.g all tests within REGRESSION_SUITE/EULER_FE/2D/IDP_ISENTROPIC_VORTEX will compile the sources located at CASES_USER/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/LOCAL_SOURCES). To add new test:

    #=== Isentropic vortex high order P2 ===#
    add_regression_test(
        TEST_DIR "HIGH_ORDER_P2"
        TEST_NAME "FE_2D_ISENTROPIC_VORTEX_HIGH_ORDER_P2"
        INCLUDE_TEST_SOURCES OFF
        FE_DIM 2D
        NB_REGEX 1
        TEST_LABELS "quick"
        CASE_LABELS "2D;ERK31;P2;hyperbolic;euler;high;consistent"
        PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4
    )

This command will add a new test located at REGRESSION_SUITE/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/HIGH_ORDER_P2:

- TEST_DIR: location at "REGRESSION_SUITE/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/HIGH_ORDER_P2"

- with name "FE_2D_ISENTROPIC_VORTEX_HIGH_ORDER_P2" (usefull for ctest -N some_string)

- INCLUDE_TEST_SOURCES

- FE_DIM: Compiled with HYPAR_SFE library of dimension 2

- NB_REGEX: Only 1 test is generated

- TEST_LABELS: "quick" (usefull for ctest -L quick to run only tests labeled with quick; ctest -LE quick to run tests not labeled with quick)

- CASE_LABLES: "2D;ERK31;P2;hyperbolic;euler;high;consistent" (ctest -L "2D|periodic" e.g to run tests which labels contain at least either 2D or periodic)

- PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4


### Source compilation:

The ctest will compile all sources found within CASES_DEVELOPPER/loc_layer_III/LOCAL_SOURCES. In this example, it means compiling CASES_USER/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/LOCAL_SOURCES. The ctests layer are therefore entirely based on the same file structure as the one that can be found within CASES_DEVELOPPER/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/LOCAL_SOURCES:

- it is required to stick to the same structure when testing on a physical problem coded within CASES_DEVELOPPER

- the corresponding folders must be found inside CASES_USER, i.e created within layer II

**INCLUDE_TEST_SOURCES**: It is possible also to compile the sources within given local_sources through this option. If set to ON, the test will compile both whatever .F90 sources found within 

- CASES_USER/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/LOCAL_SOURCES

- REGRESSION_SUITE/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/HIGH_ORDER_P2/LOCAL_SOURCES

**NB_REGEX**: useful to test convergence, different boundary conditions, different solvers, etc... see below on how to use it

**PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4**: this line designs tests with label 'mpi' and can be configured within CTEST/variables.cmake:

    set(NB_PARALLEL_TEST 5)
    set(MAX_PROC_MPI_TEST 16)

- The label "mpi" is added in cmake, no need to add it at the "TEST_LABELS" level

- **WARNING**: the tests "mpi" are designed to test parallelization and can therefore be particularly long. Therefore most tests should be run with ctest -LE mpi to exclude these tests

**NB_REGEX**: the number of tests to run simultaneously, see below on how to actually implement a hypar_sfe test

### Implement a new test

A new test is explicitly added at Layer III. The creation of layers I and/or II might be required if the test is based on a new physical problem. The newly created directory containing the new test, say REGRESSION_SUITE/EULER_FE/2D/IDP_ISENTROPIC_VORTEX/HIGH_ORDER_P2, then contains:

- Optionally LOCAL_SOURCES: if INCLUDE_TEST_SOURCES is set to ON

- REGRESSION_TESTS/data_n: as many files data_1, data_2, data_3, etc... as the value of NB_REGEX. The test will run the application n times, and check every single of these data files match the corresponding ctest. Success requires all tests to pass

For test number k in [1, n], ctest will look for the corresponding chain of characters within the input [see here](../LIBS/CMAKE_SOURCES/function_regression_test.cmake). This is automatically handled by HYPAR with:

    INTEGER :: num_test
    REAL(KIND=8), DIMENSION(syst_dim) :: tab_norm


    ...
    ...

    CALL get_num_test(num_test)
    CALL regression(tab_norm, opt_num_test=num_test, opt_tol=1d-5)

Where tab_norm contains the values to be tested and can be of whatever dimension.

When running the test for the first time, several files will be generated: 
current_regression_reference_1, regression_reference_1, current_regression_reference_2, regression_reference_2, ..., current_regression_reference_n, regression_reference_n. CTEST tests the value within current_regression_reference_k against regression_reference_k within some relative tolerance (1.d-7 by default if opt_tol not specified). Therefore it is only left to cp current_regression_reference_k regression_reference_k, and the ctest should be successful!



### Remark on tools

The folder REGRESSION_SUITE/TOOLS is particular in that it is not based on a similar folder within CASES_DEVELOPPER. It is meant to test basic features of HYPAR_SFE and its dependencies, such as PETSc and FFTW. *INCLUDE_TEST_SOURCES* is always set to ON.

### CTEST features

All tests can serve a whole a lot of different purposes which can be split into two categories: test/physical test. The user can decide which tests he wants to run according to the label he gives in. 
ctest -L "type_1|type_2":  will only run tests which include labels either type_1 or type_2
ctest -LE "type_3":        will run all runs not labeled type_3
ctest -N:                  gives the list of all available tests (one can for instance for check ctest -N -L "type_1|quick" -LE "type_3" to know how many tests will be run before effectively starting)
ctest -N -VV:              gives more detailed information about all tests
ctest --print-labels:      to obtain the full list of labels


1) Speed test: quick, tool, mpi, convergence, benchmark, performance

The label mpi should not be added manually but instead is taken care of by cmake.

2) Physical test: here is a non-exhaustive list of labels that are implemented:

- Type of test:             quick, mpi, tools

- Problem dimension:        1D, 2D, Fourier

- Matrix solver & precond:  CG, MUMPS, HYPRE, GMRES

- Temporal scheme:          ERK-11, ERK-21, ERK-31, ERK-41, ERK11, ERK31, IMEX, IRK-31

- Spatial scheme:           P1, P2, P3, P8, Pk

- Boundary conditions:      dirichlet, periodic, udotn

- Matrix inversion methods:  lumped, consistent, quasi_consistent, lumped_stokes
  
- Method for conservation eq: high, viscous, galerkin

- Problem solved:  hyperbolic, euler, stokes, navier-stokes, poisson, burgers, non_convex, linear_transport
  
- Solver post-processing: limiting, relax_minmod, relax_avg
  
- HYPAR_SFE tools:       petsc, fftw, refinement  
  
  
A given test often has several labels for its physics, and one or two for its speed. 

The mpi label is only given to tests with parameters:

- PROC_MIN 1 PROC_MAX ${MAX_PROC_MPI_TEST} DEFAULT_PROC 4

set to allow an mpi test. To disable mpi test, it is possible just to set:

- PROC_MIN 1 PROC_MAX 0 DEFAULT_PROC 4
