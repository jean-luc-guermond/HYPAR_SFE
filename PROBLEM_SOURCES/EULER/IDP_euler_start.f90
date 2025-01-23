MODULE IDP_euler_start
  USE euler_bc_arrays
  USE boundary_conditions
  USE pardiso_solve
  USE Butcher_tableau
  USE space_dim
  USE mesh_handling
  USE input_data
  USE IDP_update_euler
CONTAINS
  SUBROUTINE IDP_start_euler
    CALL read_my_data('data')
    inputs%syst_size = k_dim+2
    ERK%sv = inputs%RKsv
    CALL ERK%init
    CALL construct_mesh
    CALL construct_euler_bc
    CALL IDP_construct_euler_matrices
    !===Pardiso parameters
    isolve_euler_pardiso = -1 !===
    CALL allocate_pardiso_parameters(1)
    pardiso_param(1)%mtype = 1    !===real and structurally symmetric matrix
    pardiso_param(1)%phase = 33   !===Direct solve
    pardiso_param(1)%parm(4) = 0  !===Direct solve
  END SUBROUTINE IDP_start_euler
END MODULE IDP_euler_start
