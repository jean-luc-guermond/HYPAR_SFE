MODULE fem_tn
   PUBLIC :: ns_l1, ns_l1_PAR, ns_l2
   PRIVATE
CONTAINS
   SUBROUTINE ns_l1 (mesh, ff, t)
      !===============================================
      !  < |f| >)   ===>   t
      USE def_type_mesh
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
      REAL(KIND=8),                 INTENT(OUT) :: t
      TYPE(mesh_type)                           :: mesh
      INTEGER      ::  m, l, n
      REAL(KIND=8) :: fl
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

   SUBROUTINE ns_l1_PAR(mesh, f_in, norm_glob, communicator)
      USE def_type_mesh
      USE petsc
#include "petsc/finclude/petsc.h"
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f_in
      REAL(KIND=8),                 INTENT(OUT) :: norm_glob
      TYPE(mesh_type)                           :: mesh
      INTEGER      :: ierr
      REAL(KIND=8) :: norm_loc
      MPI_Comm :: communicator

      CALL ns_l1(mesh, f_in, norm_loc)
      CALL MPI_ALLREDUCE(norm_loc,norm_glob,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

   END SUBROUTINE ns_l1_PAR

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
