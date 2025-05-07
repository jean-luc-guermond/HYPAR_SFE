PROGRAM prog
#include "petsc/finclude/petsc.h"
  USE start_setup_MODULE
  USE setup
  USE sub_plot
  IMPLICIT NONE
  REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: un
  CHARACTER(5) :: char

  CALL start_setup

  ALLOCATE(un(mesh%np, euler%syst_dim ))
  CALL init(un, 0.d0, euler%mesh%rr)

  WRITE(char, '(I5)') euler%mesh%rank
  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:,1), 'initrho' // trim(adjustl(char)) // '.plt')

  DO WHILE(euler%time < setup_data%final_time)
     CALL euler%update(un)
     IF (euler%mesh%rank==0) write(*,*) euler%time, euler%dt
  END DO

  CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:,1), 'rho' // trim(adjustl(char)) // '.plt')
END PROGRAM prog
