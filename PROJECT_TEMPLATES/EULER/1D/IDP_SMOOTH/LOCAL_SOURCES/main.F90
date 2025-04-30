PROGRAM prog
#include "petsc/finclude/petsc.h"
  USE start_setup_MODULE
  USE sub_plot
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un
  CHARACTER(5) :: car
  
  CALL start_setup
 
  WRITE(car,'(I5)') rank
  ALLOCATE(un(mesh%np, 3))
  un(:,1) = euler%mesh%rr(1,:)**2 + 1
  un(:,2) = SIN(euler%mesh%rr(1,:))
  un(:,3) = euler%mesh%rr(1,:)**2 + 2
  !inputs%time =0.d0
  !CALL init(un,mesh%rr)
  CALL plot_1d(euler%mesh%rr(1,:),euler%pressure(un),'p' // trim(adjustl(car)) // '.plt')
END PROGRAM prog
