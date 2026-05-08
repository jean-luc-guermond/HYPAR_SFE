MODULE euler_post_proc_module
#include "petsc/finclude/petsc.h"
  USE petsc
  USE euler_type_MODULE, ONLY : euler_type
  USE petsc_tools,       ONLY : array_to_petsc_vec
  USE st_matrix,         ONLY : extract_through_ghost
  CONTAINS
    SUBROUTINE schlieren(euler,un,grad)
      IMPLICIT NONE
      CLASS(euler_type) :: euler
      REAL(KIND=8), DIMENSION(euler%mesh%np) :: un
      REAL(KIND=8), DIMENSION(euler%mesh%np) :: grad
      REAL(KIND=8) :: mxmn(2), mxmn_glob(2), mx, mn
      integer ::  k, k_dim, ierr

      !===Gradient
      k_dim = euler%mesh%gauss%k_d
      
      CALL VecSet(euler%x2vec, 0.d0, ierr)
      CALL array_to_petsc_vec(un, euler%x1vec, euler%mesh, euler%LA, 'insert')
      DO k = 1, k_dim
         CALL MatMult(euler%matrices%cij(k), euler%x1vec, euler%x3vec, ierr)
         CALL VecPointwiseMult(euler%x4vec, euler%x3vec, euler%x3vec, ierr)
         CALL VecAXPY(euler%x2vec, 1.0d0, euler%x4vec, ierr)
      END DO

      CALL VecSqrtAbs(euler%x2vec, ierr)
      CALL VecPointWiseDivide(euler%x3vec, euler%x2vec, euler%matrices%lump_mass_vec, ierr)
      CALL extract_through_ghost(euler%x3vec, euler%x2_ghost, 1, 1, euler%LA, grad, &
                                'insert', opt_assemble=.FALSE.)
      ! grad = sqrt(grad)/this%matrices%lumped_mass


      mxmn(1)=MAXVAL(grad)
      mxmn(2)=MAXVAL(-grad)
      CALL MPI_ALLREDUCE(mxmn, mxmn_glob, 2, MPI_DOUBLE, MPI_MAX, euler%communicator, ierr)
      mx = mxmn_glob(1)
      mn = -mxmn_glob(2)
      grad = EXP(-10*(grad-mn)/(mx-mn))
    END SUBROUTINE schlieren
  
END MODULE euler_post_proc_module
