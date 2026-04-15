MODULE fem_tn
  PUBLIC :: ns_l1, ns_l2
  PRIVATE
CONTAINS
  SUBROUTINE ns_l1 (mesh, ff, t)
   !===============================================
   !  < |f| >)   ===>   t
   USE def_type_mesh
   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t
   INTEGER ::  m, l, n
   REAL(KIND=8) :: fl
   TYPE(mesh_type), TARGET                     :: mesh
   t = 0
   DO m = 1, mesh%me
      DO l = 1, mesh%gauss%l_G
         fl = 0
         DO n = 1,  mesh%gauss%n_w
            fl = fl + ff(mesh%jj(n,m)) *  mesh%gauss%ww(n,l)
         END DO
         t = t + ABS(fl) *  mesh%gauss%rj(l,m)
      ENDDO
   ENDDO
 END SUBROUTINE ns_l1

 SUBROUTINE ns_l2 (mesh, ff, t)
   !===============================================
   !  < |f| >)   ===>   t
   USE def_type_mesh
   IMPLICIT NONE
   REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
   REAL(KIND=8),                 INTENT(OUT) :: t
   INTEGER ::  m, l, n
   REAL(KIND=8) :: fl
   TYPE(mesh_type), TARGET                     :: mesh
   t = 0
   DO m = 1, mesh%me
      DO l = 1, mesh%gauss%l_G
         fl = 0
         DO n = 1,  mesh%gauss%n_w
            fl = fl + ff(mesh%jj(n,m)) *  mesh%gauss%ww(n,l)
         END DO
         t = t + fl**2 * mesh%gauss%rj(l,m)
      ENDDO
   ENDDO
   t = sqrt(t)
 END SUBROUTINE ns_l2
END MODULE fem_tn
