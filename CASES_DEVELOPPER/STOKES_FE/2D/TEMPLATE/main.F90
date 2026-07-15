PROGRAM test_matrix
  USE fem_tn, ONLY: ns_l1_par
  USE fem_rhs
  USE start_setup_MODULE

  USE sub_plot
  USE post_processing_debug_MODULE
  USE my_util, ONLY: user_time
  USE space_dim

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: u_out
  REAL(KIND=8) :: tps
  INTEGER      :: k, n, tot_np
  CHARACTER(5) :: char
  CHARACTER(1), DIMENSION(2), PARAMETER :: char_comp = ["x", "y"]

  INTEGER :: ierr

!========================!
!==== INITIALIZATION ====!
!========================!

  CALL start_setup

!=====================!
!==== UPDATE STOKES ==!
!=====================!

  tps = user_time() 

  DO n=1, setup_data%max_it
    !=== Update ===!
    CALL stokes%update(un)

    !=== Verbose ===!
    IF (mod(n, setup_data%verbose_freq)==0) THEN
      IF (stokes%mesh%rank==0) write(*,*) "stokes ", n
    END IF
  END DO

  tps = user_time() - tps

!=========================!
!==== POST-PROCESSING ====!
!=========================!

  ALLOCATE(u_out(mesh%np, k_dim))
  DO k=1, k_dim
    u_out(:, k) = un(:, k+1) / un(:, 1)
  END DO

  WRITE(char, '(I5)') mesh%rank
  IF (mesh%rank==0) THEN
    tot_np = mesh%disp(mesh%nb_proc+1)-1
    WRITE(*,*) 'tot_np = ', tot_np
    WRITE(*,*) ' Time per time step per dof times proc', stokes%mesh%nb_proc*tps/(tot_np*n*(stokes%IRK%s)), tps, n*(stokes%IRK%s)
  END IF

  DO k=1, k_dim
    CALL plot_scalar_field(mesh%jj, mesh%rr, u_out(:,k), 'u' // trim(adjustl(char_comp(k)))// "_" // TRIM(ADJUSTL(char)) // '.plt')
    
    IF (setup_data%if_analytical_ref) THEN
      CALL plot_scalar_field(mesh%jj, mesh%rr, u_out(:,k)-stokes%bc%vit_anal(k, stokes%time,mesh%rr), 'error' // trim(adjustl(char_comp(k)))// "_" // trim(adjustl(char)) // '.plt')
      CALL plot_scalar_field(mesh%jj, mesh%rr, stokes%bc%vit_anal(k, stokes%time,mesh%rr), 'ref' // trim(adjustl(char_comp(k)))// "_" // trim(adjustl(char)) // '.plt')
    END IF
  END DO

!=========================!
!==== REGRESSION TEST ====!
!=========================!

  CALL errors

!=====================!
!==== END PROGRAM ====!
!=====================!
  CALL PetscFinalize(ierr)

CONTAINS

  SUBROUTINE errors
    USE fem_tn
    USE post_processing_debug_MODULE
    USE options_module
    IMPLICIT NONE

    REAL(KIND=8) :: error, norm, norm_anal
    REAL(KIND=8), DIMENSION(k_dim) :: tab_norm

!==== Put final processing stuff here ====!
    DO k=1, k_dim
      IF (setup_data%if_analytical_ref) THEN
        CALL ns_l1_PAR(mesh, u_out(:,k)-stokes%bc%vit_anal(k,stokes%time,mesh%rr), error, communicator)
        CALL ns_l1_PAR(mesh, stokes%bc%vit_anal(k,stokes%time,mesh%rr), norm_anal, communicator)
        norm = error/norm_anal
        IF(mesh%rank==0) WRITE(*, *) 'Comp = ',k,'; Relative error, L1-norm = ', error/norm_anal
      ELSE
        CALL ns_l1_PAR(mesh, u_out(:,k), norm, communicator)
        IF(mesh%rank==0) WRITE(*, *) 'Comp = ',k,'; no analytical ref, L1-norm = ', norm
      END IF
    END DO

!==== For regression tests ====!
    IF (options%if_regression) THEN
      DO k=1, k_dim
        IF (setup_data%if_analytical_ref) THEN
          CALL ns_l1_PAR(mesh, u_out(:,k)-stokes%bc%vit_anal(k,stokes%time,mesh%rr), error, communicator)
          CALL ns_l1_PAR(mesh, stokes%bc%vit_anal(k,stokes%time,mesh%rr), norm_anal, communicator)
          norm = error/norm_anal
        ELSE
          CALL ns_l1_PAR(mesh, u_out(:,k), norm, communicator)
        END IF
        tab_norm(k) = norm
      END DO
      CALL regression(tab_norm, opt_num_test=options%num_regex)
    END IF

  END SUBROUTINE errors

END PROGRAM test_matrix
