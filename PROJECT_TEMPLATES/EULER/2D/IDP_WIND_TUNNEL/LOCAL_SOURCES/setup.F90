MODULE setup
   USE mesh_parameters
   PUBLIC :: sol_anal, init, rho_anal, press_anal, mt_anal, E_anal, impose_bc_euler, pressure
   PRIVATE
   REAL(KIND=8), PARAMETER :: r0=0.15d0, x0=0d0, y0=0.0d0
   REAL(KIND=8), PARAMETER :: u_infty=0.d0, rho_infty=1.d0, p_infty=1.d0, beta0=5.d0, gamma = 1.4d0
   CONTAINS
   
!==========================================================================
!================= DEF PRESSURE FOR SETUP =================================
!==========================================================================

   FUNCTION pressure(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      REAL(KIND = 8) :: gamma
      gamma = 7.0 / 5.0
      vv = rho * e * (gamma - 1)
   END FUNCTION pressure

!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================
   
   SUBROUTINE impose_bc_euler(un, euler_bc, mesh, time)
      USE euler_bc_arrays
      USE def_type_mesh
      TYPE(mesh_type) :: mesh
      TYPE(euler_bc_type) :: euler_bc
      REAL(KIND = 8) :: time
      REAL(KIND = 8), DIMENSION(:, :), INTENT(INOUT) :: un
      REAL(KIND=8), DIMENSION(SIZE(euler_bc%udotn_bc%jsd)) :: mdotn
      INTEGER :: comp
      
      DO comp = 1, euler_bc%syst_dim
         SELECT CASE(comp)
         CASE(1)
            un(euler_bc%rho_bc%jsd, comp) = rho_anal(time, mesh%rr(:, euler_bc%rho_bc%jsd))
         CASE(2)
            un(euler_bc%ux_bc%jsd, comp) = mt_anal(comp - 1, time, mesh%rr(:, euler_bc%ux_bc%jsd))
         CASE(3)
            un(euler_bc%uy_bc%jsd, comp) = mt_anal(comp - 1, time, mesh%rr(:, euler_bc%uy_bc%jsd))
         END SELECT
      END DO
      IF (SIZE(euler_bc%udotn_bc%jsd).NE.0) THEN
         mdotn = euler_bc%udotn_normal_vtx(:,1)*un(euler_bc%udotn_bc%jsd,2) &
         +  euler_bc%udotn_normal_vtx(:,2)*un(euler_bc%udotn_bc%jsd,3)
         un(euler_bc%udotn_bc%jsd,2) = un(euler_bc%udotn_bc%jsd,2) - mdotn*euler_bc%udotn_normal_vtx(:,1)
         un(euler_bc%udotn_bc%jsd,3) = un(euler_bc%udotn_bc%jsd,3) - mdotn*euler_bc%udotn_normal_vtx(:,2)
         
         mdotn = euler_bc%udotn_normal_vtx(:,1)*un(euler_bc%udotn_bc%jsd,2) &
         +  euler_bc%udotn_normal_vtx(:,2)*un(euler_bc%udotn_bc%jsd,3)
      END IF
      
   END SUBROUTINE impose_bc_euler
   
   SUBROUTINE init(un, time, rr)
     USE def_of_gamma
     USE lambda_module
     USE mesh_parameters
     USE petsc
     USE my_util
     IMPLICIT NONE
     REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
     REAL(KIND = 8), DIMENSION(SIZE(rr, 2), mesh_data_info%k_dim + 2), INTENT(OUT) :: un
     REAL(KIND = 8)   :: time, my_time
     CHARACTER(len=5) :: char
     INTEGER :: ierr, rank, my_rank
     CALL MPI_Comm_rank(PETSC_COMM_WORLD, my_rank, ierr)
     IF (time<0.d0) THEN !<==Restart
        WRITE(char, '(I5)') my_rank
        OPEN(unit = 10, &
             file = 'restart_'//TRIM(ADJUSTL(char))//'_'//TRIM(ADJUSTL(mesh_data_info%file_name)),&
             form = 'unformatted', status = 'unknown', err=100)
        READ(10,err=101,END=101) rank, my_time, un
        time = my_time
        IF (rank/=my_rank) THEN
           CALL error_petsc('Error in setup: Wrong processor mapping')
        END IF
        IF (my_rank==0) WRITE(*,*) ' time at checkpoint restart: ', time
        CLOSE(10)
        RETURN
100     CONTINUE
        CALL error_petsc('Error in setup: error opening restart files. Wrong number of procs?')
101     CONTINUE
        CALL error_petsc('Error in setup: error reading restart files.')      
     ELSE
        un(:, 1) = rho_anal(time, rr)
        un(:, 2) = mt_anal(1, time, rr)
        un(:, 3) = mt_anal(2, time, rr)
        un(:, 4) = E_anal(time, rr)
     END IF
   END SUBROUTINE init
   
   FUNCTION rho_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      IF (SIZE(vv)==0) RETURN
      vv = gamma
   END FUNCTION rho_anal
   
   FUNCTION press_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      IF (SIZE(vv)==0) RETURN
      vv = 1.d0
   END FUNCTION press_anal
   
   FUNCTION vit_anal(comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      IF (SIZE(vv)==0) RETURN
      IF (comp==1) THEN
         IF (time<1.d-8) THEN
            vv = 3.d0
            RETURN
         END IF
         DO n = 1, SIZE(vv)
            IF (rr(1, n)<1.d-8) THEN
               vv(n) = 3.0
            ELSE
               vv(n) = 0.d0
            END IF
         END DO
      ELSE IF (comp==2) THEN
         vv = 0.d0
      ELSE
         WRITE(*, *) ' BUG '
         STOP
      END IF
   END FUNCTION vit_anal
   
   FUNCTION E_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = press_anal(time, rr) / (gamma - 1.d0) &
      + rho_anal(time, rr) * (vit_anal(1, time, rr)**2 + vit_anal(2, time, rr)**2) / 2
   END FUNCTION E_anal
   
   FUNCTION mt_anal(comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = rho_anal(time, rr) * vit_anal(comp, time, rr)
   END FUNCTION mt_anal
   
   
   FUNCTION sol_anal(comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      SELECT CASE(comp)
      CASE(1)
         vv = rho_anal(time, rr)
      CASE(2)
         vv = mt_anal(1, time, rr)
      CASE(3)
         vv = mt_anal(2, time, rr)
      CASE(4)
         vv = E_anal(time, rr)
      CASE DEFAULT
         WRITE(*, *) ' BUG in sol_anal'
         STOP
      END SELECT
   END FUNCTION sol_anal
END MODULE setup
