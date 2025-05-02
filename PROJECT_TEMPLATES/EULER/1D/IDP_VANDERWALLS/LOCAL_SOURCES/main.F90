PROGRAM prog
#include "petsc/finclude/petsc.h"
   USE start_setup_MODULE
   USE sub_plot
   USE my_util
   IMPLICIT NONE
   REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: un
   CHARACTER(5) :: car
   REAL(KIND = 8) :: t0
   INTEGER :: n, i

   CALL start_setup

   WRITE(car, '(I5)') rank
   ALLOCATE(un(mesh%np, 3))
   CALL init(un, 0.d0, euler%mesh%rr)

   CALL plot_1d(euler%mesh%rr(1, :), un(:,1), 'initrho' // trim(adjustl(car)) // '.plt')
   CALL plot_1d(euler%mesh%rr(1, :), un(:,2)/un(:,1), 'initvt' // trim(adjustl(car)) // '.plt')
   CALL plot_1d(euler%mesh%rr(1, :), un(:,3), 'initE' // trim(adjustl(car)) // '.plt')

   euler%cfl = 0.5

   t0 = user_time()
   i = 0
   DO WHILE (euler%time < .1 )
      CALL euler%update(un)
      i = i + 1
   END DO
   write(*,*)  'th', i*SUM(euler%mesh%domnp)/(user_time() - t0)
   CALL plot_1d(euler%mesh%rr(1, :), un(:,1), 'rho' // trim(adjustl(car)) // '.plt')
   CALL plot_1d(euler%mesh%rr(1, :), un(:,2), 'mt' // trim(adjustl(car)) // '.plt')
END PROGRAM prog
