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
  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:,1), 'initrho'//TRIM(ADJUSTL(char))//'.plt')

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
     WRITE(*,*) ' Time per time step per dof times proc', tps/(tot_np*n), tps, n
  END IF

  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:, 1), 'rho' // TRIM(ADJUSTL(char)) // '.plt')
  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:, 2)/un(:, 1), 'ux' // TRIM(ADJUSTL(char)) // '.plt')
  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:, 3)/un(:, 1), 'uy' // TRIM(ADJUSTL(char)) // '.plt')
END PROGRAM prog
