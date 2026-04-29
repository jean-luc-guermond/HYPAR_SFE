PROGRAM prog
#include "petsc/finclude/petsc.h"
  USE petsc
  USE start_setup_MODULE
  USE setup
  USE sub_plot
  USE my_util
  USE plot_vtu_module
  USE euler_post_proc_module
  IMPLICIT NONE
  REAL(KIND = 8), DIMENSION(:, :), ALLOCATABLE :: un
  REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: grad
  REAL(KIND = 8) :: tps, t_plot
  CHARACTER(5) :: char
  INTEGER :: n, tot_np, code, it, it_max, it_plot

  CALL start_setup

  !===Allocate
  ALLOCATE(un(mesh%np, euler%syst_dim ))
  ALLOCATE(grad(mesh%np))

  !===Restart
  IF (setup_data%if_restart) THEN
     euler%time = -1.d0
     CALL init(un, euler%time, euler%mesh%rr)
     it_plot = euler%time/setup_data%checkpointing_freq                                        
     t_plot = (it_plot+1)*setup_data%checkpointing_freq
  ELSE
     CALL init(un, 0.d0, euler%mesh%rr)
     it_plot = 0                                                                 
     t_plot = setup_data%checkpointing_freq  
  END IF

  grad =0.d0
  CALL make_vtu_file_2D(euler%communicator, euler%mesh, 'test', un(:,1), 'Density', 'new', opt_it=0)
  CALL make_vtu_file_2D(euler%communicator, euler%mesh, 'test', grad, 'schlieren', 'old', opt_it=0)
  tps = user_time()
  n = 0
  DO WHILE(euler%time < setup_data%final_time)
     CALL euler%update(un)
     n = n + 1
     IF (euler%mesh%rank==0) write(*,*) euler%time, euler%dt

     IF ((euler%time<t_plot .AND. euler%time+euler%dt/2>t_plot) .OR. &        
          (euler%time.GE.t_plot .AND. euler%time-euler%dt/2<t_plot)) THEN        
        it_plot = it_plot + 1                                                                                            
        t_plot = t_plot+setup_data%checkpointing_freq                                          
        CALL schlieren(euler,un(:,1),grad)
        !WRITE(char,'(I5)') euler%mesh%rank
        !CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, grad, &
        !     'grad' // TRIM(ADJUSTL(char)) // '.plt')
        CALL make_vtu_file_2D(euler%communicator, euler%mesh, 'test', un(:,1), &
             'Density', 'old', opt_it=it_plot)
        CALL make_vtu_file_2D(euler%communicator, euler%mesh, 'test', grad, &
             'schlieren', 'old', opt_it=it_plot)
     END IF
  END DO
  tps = user_time() - tps

  CALL MPI_ALLREDUCE(euler%mesh%dom_np,tot_np,1,MPI_INTEGER,MPI_SUM,euler%communicator,code)
  IF(euler%mesh%rank==0) THEN
     WRITE(*,*) ' tot_np', tot_np
     WRITE(*,*) ' Time per time step per dof times proc', tps/(tot_np*n), tps, n
  END IF
  !WRITE(char,'(I5)') euler%mesh%rank
  !CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:,1), &
  !           'density' // TRIM(ADJUSTL(char)) // '.plt')
  !CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:, 2)/un(:, 1), 'ux' // TRIM(ADJUSTL(char)) // '.plt')
  !CALL plot_scalar_field(euler%mesh%jj, euler%mesh%rr, un(:, 3)/un(:, 1), 'uy' // TRIM(ADJUSTL(char)) // '.plt')

  CALL write_restart(un)

CONTAINS
  SUBROUTINE write_restart(un)
    USE mesh_parameters
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: un
    CHARACTER(len=5) :: char
    WRITE(char, '(I5)') euler%mesh%rank
    OPEN(unit = 10, &
         file = 'restart_'//trim(adjustl(char))//'_'//trim(adjustl(mesh_data_info%file_name)),&
         form = 'unformatted', status = 'unknown')
    WRITE(10) euler%mesh%rank, euler%time, un
    WRITE(*,*) ' time at checkpoint', euler%time
    CLOSE(10)
  END SUBROUTINE write_restart
END PROGRAM prog
