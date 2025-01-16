PROGRAM test_matrix
   USE construct_mesh
   USE input_data
   USE timing_tools
   USE def_type_mesh

   TYPE(mesh_type) :: mesh
   TYPE(petsc_csr_LA) :: LA
   Mat :: mass
   INTEGER, POINTER, DIMENSION(:) :: js_d_loc

   INTEGER :: rank
   !===Start PETSC and MPI (mandatory)=============================================
   CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
   CALL create_cart_comm(inputs%ndim, comm_cart, comm_one_d, coord_cart)

   !===User reads his/her own data=================================================
   CALL read_user_data('data')

   CALL construct_mesh(communicator, mesh, LA, js_d_loc, 1)
   CALL create_local_petsc_matrix(PETSC_COMM_WORLD, LA, mass, clean = .FALSE.)
   CALL qs_mass_diff_M (mesh, 1.d0, 0.d0, LA, mass)

END PROGRAM test_matrix