MODULE construct_mesh
   USE def_type_mesh
   USE space_dim
   USE input_data
   TYPE(mesh_type), PUBLIC :: mesh
   PUBLIC :: construct_mesh
   PRIVATE
CONTAINS
   SUBROUTINE construct_mesh(opt_edge_stab)
      USE input_data
      USE prep_maill
      USE HCT_mesh
      USE Powell_Sabin_mesh
      IMPLICIT NONE

      LOGICAL, OPTIONAL :: opt_edge_stab
      LOGICAL :: edge_stab
      TYPE(mesh_type) :: p1_mesh
      IF (.NOT.PRESENT(opt_edge_stab)) THEN
         edge_stab = .FALSE.
      ELSE
         edge_stab = opt_edge_stab
      END IF
      SELECT CASE(k_dim)
      CASE(2)
         CALL mesh_2d(mesh)
      CASE(1)
         CALL load_mesh_1d(mesh)
      CASE DEFAULT
         write(*, *) ' BUG in construct_mesh, k_dim not correct'
         STOP
      END SELECT
   END SUBROUTINE construct_mesh
END MODULE  construct_mesh