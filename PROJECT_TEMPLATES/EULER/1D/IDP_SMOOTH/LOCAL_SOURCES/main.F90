PROGRAM prog
#include "petsc/finclude/petsc.h"
  USE start_setup_MODULE
  USE setup
  USE sub_plot
  IMPLICIT NONE
  REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: un
  CHARACTER(5) :: car
  INTEGER :: n
  CALL start_setup

  WRITE(car, '(I5)') euler%mesh%rank
  ALLOCATE(un(mesh%np, euler%syst_dim ))
  CALL init(un, 0.d0, euler%mesh%rr)

  CALL plot_1d(euler%mesh%rr(1, :), un(:,1), 'initrho' // trim(adjustl(car)) // '.plt')

  DO WHILE(euler%time < setup_data%final_time)
     IF (euler%mesh%rank==0) write(*,*) euler%time, euler%dt
     CALL euler%update(un)
  END DO

  CALL plot_1d(euler%mesh%rr(1, :), un(:,1), 'rho' // trim(adjustl(car)) // '.plt')
END PROGRAM prog
