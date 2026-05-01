MODULE euler_post_proc_module
#include "petsc/finclude/petsc.h"
  USE petsc
  USE euler_type_MODULE, ONLY : euler_type
  USE petsc_tools,       ONLY : array_to_petsc_vec
  USE st_matrix,         ONLY : extract_through_ghost
  CONTAINS
    SUBROUTINE schlieren(this,un,grad)
      implicit none
      TYPE(euler_type) :: this
      REAL(KIND=8), DIMENSION(this%mesh%np) :: un
      REAL(KIND=8), DIMENSION(this%mesh%np) :: grad
      REAL(KIND=8) :: mxmn(2), mxmn_glob(2), mx, mn
      integer ::  k, k_dim, ierr

      !===Gradient
      k_dim = this%mesh%gauss%k_d
      
      CALL VecSet(this%x2vec, 0.d0, ierr)
      CALL array_to_petsc_vec(un, this%x1vec, this%mesh, this%LA, 'insert')
      DO k = 1, k_dim
         CALL MatMult(this%matrices%cij(k), this%x1vec, this%x3vec, ierr)
         CALL VecPointwiseMult(this%x4vec, this%x3vec, this%x3vec, ierr)
         CALL VecAXPY(this%x2vec, 1.0d0, this%x4vec, ierr)
      END DO

      CALL VecSqrtAbs(this%x2vec, ierr)
      CALL VecPointWiseDivide(this%x3vec, this%x2vec, this%matrices%lump_mass_vec, ierr)
      CALL extract_through_ghost(this%x3vec, this%x2_ghost, 1, 1, this%LA, grad, &
                                'insert', opt_assemble=.FALSE.)
      ! grad = sqrt(grad)/this%matrices%lumped_mass


      mxmn(1)=MAXVAL(grad)
      mxmn(2)=MAXVAL(-grad)
      CALL MPI_ALLREDUCE(mxmn, mxmn_glob, 2, MPI_DOUBLE, MPI_MAX, this%communicator, ierr)
      mx = mxmn_glob(1)
      mn = -mxmn_glob(2)
      grad = EXP(-10*(grad-mn)/(mx-mn))
    END SUBROUTINE schlieren
  
END MODULE euler_post_proc_module
