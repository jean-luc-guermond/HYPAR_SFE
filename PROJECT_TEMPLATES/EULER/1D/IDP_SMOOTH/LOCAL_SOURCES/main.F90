PROGRAM prog
#include "petsc/finclude/petsc.h"
   USE start_setup_MODULE
   USE sub_plot
   IMPLICIT NONE
   REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: un
   CHARACTER(5) :: car
   INTEGER :: n

   CALL start_setup

   WRITE(car, '(I5)') rank
   ALLOCATE(un(mesh%np, 3))
   CALL init(un, 0.d0, euler%mesh%rr)

   CALL plot_1d(euler%mesh%rr(1, :), un(:,1), 'initrho' // trim(adjustl(car)) // '.plt')

   euler%dt = 1.d-2 / real(SUM(euler%mesh%domnp))

   DO n = 1, 1000
      CALL euler%update(un)
   END DO

   CALL plot_1d(euler%mesh%rr(1, :), un(:,1), 'rho' // trim(adjustl(car)) // '.plt')
END PROGRAM prog
