Hyperbolic class
================

There are currently two classes provided by HYPAR_SFE, namely:

- Linear transport

- Euler

Additional examples not present in PROBLEM_SOURCES are present within CASES_USER:

- Burgers

The abstract class provides the general numerical scheme for solving conservation equations, both spatial (FEM with dij method) and temporal (RK methods). The user needs to program the followings for the application to be complete:

1) flux
2) compute_lambda
3) construct_bc
4) impose_bc

The state vector must be un(mesh%np, syst_dim) where syst_dim is the problem dimension. For example for Euler:

- Linear transport: scalar conservation equation
- 1D Euler: U = (density, ux, E_tot)
- 2D Euler: U = (density, ux, uy, E_tot)

### Flux:
Simply the flux function:

    FUNCTION interface_flux(this, comp, un) RESULT(vv)
        USE space_dim, ONLY: k_dim
        IMPLICIT NONE
        CLASS(my_hyperbolic_type),          INTENT(INOUT) :: this
        REAL(KIND = 8), DIMENSION(:, :),       INTENT(IN) :: un
        INTEGER,                               INTENT(IN) :: comp
        REAL(KIND = 8), DIMENSION(SIZE(un, 1), k_dim)     :: vv
    END FUNCTION interface_flux

which defines the flux for the conservation equation which is considered. comp ranges from 1 to syst_dim.

### compute_lambda:
An upper bound of the speed of propagation of the Riemann fan of the Riemann problem associated to the conservation equation considered:

    SUBROUTINE interface_lambda(this, un, i, j, nij, lambda_max)
        USE space_dim
        IMPLICIT NONE
        CLASS(my_hyperbolic_type),                        INTENT(INOUT) :: this
        REAL(KIND=8), DIMENSION(this%mesh%np, ${syst_dim}$), INTENT(IN) :: un
        INTEGER,                                             INTENT(IN) :: i, j
        REAL(KIND=8), DIMENSION(k_dim),                      INTENT(IN) :: nij
        REAL(KIND=8), DIMENSION(2),                   INTENT(OUT) :: lambda_max
    END SUBROUTINE interface_lambda

The Riemann problem is defined here with U_L = un(i, :) and U_R = un(j, :)

### construct_bc:
A subroutine which will construct the objects meant to impose boundary conditions:

    SUBROUTINE interface_construct_bc(this, mesh, LA)
        USE def_type_mesh
        USE petsc_csr_LA_module
        IMPLICIT NONE
        CLASS(my_hyperbolic_type), INTENT(INOUT) :: this
        TYPE(mesh_type)           :: mesh
        TYPE(petsc_csr_LA)        :: LA
    END SUBROUTINE interface_construct_bc

Usually this subroutine is just a succession of *CALL this%bc%var%set(mesh, "var_name")*. This line will add a line to the data file asking the user for the boundaries where Dirichlet boundary conditions should be imposed to the quantity *var*. So the user should just copy/paste this line for every variable on which he wishes to impose Dirichlet boundary conditions. 
**Remark**: the u.n = f boundary condition can also be enforced. For this the user writes *CALL this%bc%udotn_bc%set(mesh, "var_name")* as usual, and needs to add the following call provided by the hyperbolic framework: *CALL this%construct_udotn(mesh, LA, this%bc%udotn_bc, this%bc%udotn_normal_vtx)*

### impose_bc:
A subroutine meant to impose the boundary conditions explicitly:

    SUBROUTINE interface_impose_bc(this, un, mesh, time)
        USE def_type_mesh
        IMPLICIT NONE
        CLASS(my_hyperbolic_type),       INTENT(INOUT) :: this
        TYPE(mesh_type)                                :: mesh
        REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
        REAL(KIND = 8), INTENT(IN)                     :: time
    END SUBROUTINE interface_impose_bc

This subroutine is a succession of calls *un(this%bc%var%jsd, comp_sys) = some_function(time, mesh%rr(:, this%bc%val%jsd))* where un is the state vector and the user defines some_function which contains the formula to be applied to un(:,comp_sys) on the boundaries specified within the data file.

We refer to *HYPAR_SFE/PROBLEM_SOURCES/EULER/euler_type_module.F90* for an example on how to impose u.n boundary conditions.

These four subroutines are mandatory for the user to define in his own type extending the abstract one *hyperbolic_type*. Optionally he can define some other functions that will allow him to fully use the high order methods capabilities.

### eta_commute:
      FUNCTION interface_eta_commute(un) RESULT(eta)
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(:, :), INTENT(IN)    :: un
         REAL(KIND=8), DIMENSION(SIZE(un,1))         :: eta
      END FUNCTION interface_eta_commute

Takes the state vector as argument, returns a scalar which represents the smoothness indicator for computing the commutator (typically the density or the pressure in Euler)

### Limiting:
The user can define as many limiting functionals as wanted. Here is the way to do it for a given functional
1) Write a file named *limiting_my_hyperbolic_type_limiter_name.inc* (e.g *limiting_euler_rho_min.inc*) which contains three subroutines: 

- **psi_my_hyperbolic_type_limiter_name**: computes psi(U) - psi_min. e.g its definition is straightforward for imposing a minimum density in Euler or any scalar conservation equation where U(:,1) represents the density within the state vector:

        PURE FUNCTION psi_euler_rho_min(x,psi_m) RESULT(v)
            IMPLICIT NONE
            REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: x
            REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: psi_m
            REAL(KIND=8), DIMENSION(SIZE(psi_m))     :: v
            INTEGER :: n
            DO n=1, SIZE(x, 1)
                v(n) = x(n,1)-psi_m(n)
            END DO
        END FUNCTION psi_euler_rho_min

- **zero_of_psi_my_hyperbolic_type_limiter_name**: given a direction vector P, the way of computing lambda such that psi(U + lambda*P) = psi_min:

        PURE FUNCTION zero_of_psi_euler_rho_min(psi_m,u0,P) RESULT(v)
            IMPLICIT NONE
            REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: u0, P
            REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: psi_m
            REAL(KIND=8), DIMENSION(SIZE(psi_m)) :: v
            INTEGER :: n
            DO n=1, SIZE(u0, 1)
                v(n) = (psi_m(n)-u0(n,1))/P(n,1)
            END DO
        END FUNCTION zero_of_psi_euler_rho_min

Notice that the two above limiters take x(:,:) and u0(:,:) as arguments. These are state vectors where the first dimension is the number of state vectors and the second one is the dimension of a state vector. 

- **scal_zero_of_psi_my_hyperbolic_type_limiter_name**: the exact same as  *zero_of_psi_my_hyperbolic_type_limiter_name* but taking a single state vector: 

        PURE FUNCTION scal_zero_of_psi_euler_rho_min(psi_m,u0,P) RESULT(v)
            IMPLICIT NONE
            REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u0, P
            REAL(KIND=8), INTENT(IN)  :: psi_m
            REAL(KIND=8) :: v
            v = (psi_m-u0(1))/P(1)
        END FUNCTION scal_zero_of_psi_euler_rho_min


**Remark 1**: It is very important that those files only contain these subroutines, and nothing else. These files are never compiled as such but rather copied/pasted into modules, or into CONTAINS blocks of some subroutines which are templated within .fpp files in order to force inlining.

**Remark 2**: these subroutines must absolutely be called **psi_my_hyperbolic_type_limiter_name**, **zero_of_psi_my_hyperbolic_type_limiter_name**, **scal_zero_of_psi_my_hyperbolic_type_limiter_name**

2) The user must modify the python file *python_precompile.py* accordingly. An example on how to fill it can be found inside CASES_USER/SCALAR_CONS_FE/BURGERS/1D/DEFAULT. The user needs to fill:

- path_to_hypar

- path_to_application

- syst_dim

- n_lim (number of limiting functionals to compile)

- list_lim (list of limiting functionals' names)

With the example from below, the precompilation will look for two limiting functionals, inside the path "path_to_my_application/LIMITING_FUNCTIONALS" under the form:

1) limiting_scalar_rho_min.inc

2) limiting_scalar_rho_max.inc

They will be compiled based on "path_pb", i.e within "path_to_my_application/LOCAL_SOURCES/TEMPLATED" and with the name burgers

        import os, sys
        from pathlib import Path
        dir_top = "path_to_hypar/HYPAR_SFE"
        path_to_application = "path_to_my_application"
        sys.path.append(str(Path(dir_top) / "LIBS/PYTHON_SOURCES"))
        from precompile_sources import define_paths, precompile_hyperbolic

        define_paths(dir_top)
        path_pb = os.getcwd()

        #=== Definitions for BURGERS problem in 1D
        syst_dim = 1
        n_lim = 2
        list_lim = ['scalar_rho_min', 'scalar_rho_max']
        path_lim = str(Path(path_to_application) / Path("LIMITING_FUNCTIONALS"))
        path_pb = str(Path(path_to_application) / Path("LOCAL_SOURCES"))
        name = 'burgers'

        #=== The line performing precompilation of all source files discussed above
        precompile_hyperbolic(syst_dim, n_lim, list_lim, path_lim, path_pb, name)


3) The execution of the python source generation files will generate specialized modules for every given limiter which allows inlining and therefore a significant performance optimization. These modules can be visualized at the *cmake ..* level within a newly created folder called *TEMPLATED* (located inside *LOCAL_SOURCES*). The list of modules is:
- **abstract_hyperbolic_my_hyperbolic.F90 (MODULE my_hyperbolic_abstract_hyperbolic_module)** (e.g *euler_abstract_hyperbolic_module*) from which the user gets *hyperbolic_type*.
- **uij_bar_bounds_my_hyperbolic.F90 (MODULE my_hyperbolic_uij_bar_bounds)**: specialized subroutines for computing limiting bounds based on $\bar{u}_{ij}$ from which to get *my_hyperbolic_compute_bounds_uijbar*
- **limiting_my_hyperbolic.F90 (MODULE limiting_functionals_my_hyperbolic_module)**: copy and paste of limiting functionals defined by the user within his .inc files
- **limiter_cell_elt_my_hyperbolic.F90 (MODULE my_hyperbolic_limiter_cell_elt)**: contains the specialized subroutines performing convex limiting based on the elements.
- **cell_limiting_engine_parallel_my_hyperbolic.F90**: defines the limiting type specialized for the given problem. This module is called only used by my_hyperbolic_abstract_hyperbolic_module, the shouldn't need to use it.

With these newly created modules, the *start_setup.F90* should look like:


    MODULE start_setup_MODULE
        USE my_hyperbolic_cell_limiting_engine_parallel_module
        USE limiting_functionals_my_hyperbolic_module, ONLY: psi_my_hyperbolic_limiter_name_1, zero_of_psi_my_hyperbolic_limiter_name_1, &
                                                             psi_my_hyperbolic_limiter_name_2, zero_of_psi_my_hyperbolic_limiter_name_2
        USE my_hyperbolic_limiter_cell_elt
        USE my_hyperbolic_uij_bar_bounds

    ...
    ...
    ...

        SUBROUTINE start_setup
            ...
            ...
            ...
            !=== Define my_hyperbolic limiting bounds (should we put this in PROBLEM_SOURCES instead?)
            my_hyperbolic%compute_bounds_uijbar => my_hyperbolic_compute_bounds_uijbar

            ALLOCATE(limiting_functionals_my_hyperbolic(2))
            limiting_functionals_my_hyperbolic(1)%psi => psi_my_hyperbolic_limiter_name_1
            limiting_functionals_my_hyperbolic(1)%zero_of_psi => zero_of_psi_my_hyperbolic_limiter_name_1
            limiting_functionals_my_hyperbolic(1)%name = limiter_name_1
            limiting_functionals_my_hyperbolic(2)%psi => psi_my_hyperbolic_limiter_name_2
            limiting_functionals_my_hyperbolic(2)%zero_of_psi => zero_of_psi_my_hyperbolic_limiter_name_2
            limiting_functionals_my_hyperbolic(2)%name = limiter_name_2

            limiting_functionals_my_hyperbolic(1)%spe_iterative_cell_limiting_procedure => iterative_cell_limiting_procedure_my_hyperbolic_limiter_name
            limiting_functionals_my_hyperbolic(2)%spe_iterative_cell_limiting_procedure => iterative_cell_limiting_procedure_my_hyperbolic_limiter_name
            !=== Define my_hyperbolic limiting bounds
            
            my_hyperbolic%erk_sv = setup_data%erk_sv
            CALL my_hyperbolic%init_hyperbolic(communicator, name, mesh, limiting_functionals_my_hyperbolic)
            CALL init_state_functions(my_hyperbolic, setup_data%which_init)

            CALL my_hyperbolic%set_times(times)

        END SUBROUTINE start_setup
    END MODULE start_setup_MODULE


The simplest example is the application of this logic on Burgers' equation in 1D.