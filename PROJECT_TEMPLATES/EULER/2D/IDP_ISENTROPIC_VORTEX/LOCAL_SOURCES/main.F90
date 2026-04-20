PROGRAM prog
#include "petsc/finclude/petsc.h"
  USE petsc
  USE start_setup_MODULE
  USE setup
  USE sub_plot
  USE my_util
  IMPLICIT NONE
  REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: un
  REAL(KIND = 8) :: tps
  CHARACTER(5) :: char
  INTEGER :: n, tot_np, code

  CALL start_setup

  ALLOCATE(un(mesh%np, euler%syst_dim ))
  CALL init(un, 0.d0, euler%mesh%rr)

  WRITE(char, '(I5)') euler%mesh%rank
  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:,1), 'initrho'//trim(adjustl(char))//'.plt')

  tps = user_time()
  n = 0
  DO WHILE(euler%time < setup_data%final_time)
     CALL euler%update(un)
     n = n + 1
     IF (euler%mesh%rank==0) write(*,*) euler%time, euler%dt
  END DO
  tps = user_time() - tps

  CALL MPI_ALLREDUCE(euler%mesh%dom_np,tot_np,1,MPI_INTEGER,MPI_SUM,euler%communicator,code)
  IF(euler%mesh%rank==0) THEN
     WRITE(*,*) ' tot_np', tot_np
     WRITE(*,*) ' Time per time step per dof times proc', tps/(tot_np*n), tps
  END IF

  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:,1), 'rho' // trim(adjustl(char)) // '.plt')
  CALL errors
CONTAINS
  SUBROUTINE errors
    USE fem_tn
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    REAL(KIND=8) :: error_loc, norm_loc, error, norm
    INTEGER :: code
    CALL ns_l1(mesh, un(:,1)-rho_anal(euler%time,mesh%rr), error_loc)
    CALL MPI_ALLREDUCE(error_loc,error,1,MPI_DOUBLE_PRECISION,MPI_SUM,euler%communicator,code)
    CALL ns_l1(mesh, rho_anal(euler%time,mesh%rr), norm_loc)
    CALL MPI_ALLREDUCE(norm_loc,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,euler%communicator,code)
    if(euler%mesh%rank==0) WRITE(*, '(A,g12.3)') 'L1 NORM error ', error/norm
  END SUBROUTINE errors
END PROGRAM prog
