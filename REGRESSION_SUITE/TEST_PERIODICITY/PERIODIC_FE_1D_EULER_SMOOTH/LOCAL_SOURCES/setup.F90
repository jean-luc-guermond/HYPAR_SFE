MODULE setup
   USE space_dim
   PUBLIC :: sol_anal, init, rho_anal, press_anal, mt_anal, E_anal, impose_bc_euler, pressure
   PRIVATE
   REAL(KIND = 8), PARAMETER, PRIVATE :: x0=0.1d0, x1=0.3d0, gamma=1.4d0

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
      INTEGER :: comp

      DO comp = 1, euler_bc%syst_dim
         SELECT CASE(comp)
         CASE(1)
            un(euler_bc%rho_bc%jsd, comp) = rho_anal(time, mesh%rr(:, euler_bc%rho_bc%jsd))
         CASE(2:k_dim + 1)
            un(euler_bc%rho_bc%jsd, comp) = mt_anal(comp - 1, time, mesh%rr(:, euler_bc%rho_bc%jsd))
         CASE(k_dim + 2)
            un(euler_bc%rho_bc%jsd, comp) = E_anal(time, mesh%rr(:, euler_bc%rho_bc%jsd))
         END SELECT
      END DO

   END SUBROUTINE impose_bc_euler

   SUBROUTINE init(un, time, rr)
      USE def_of_gamma
      USE lambda_module
      USE my_util, ONLY : error_petsc, to_str
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), k_dim + 2), INTENT(OUT) :: un
      REAL(KIND = 8), INTENT(IN) :: time
      INTEGER :: comp
      ! gamma = 1.4d0
      ! x0 = 0.1
      ! x1 = 0.3
      DO comp=1, SIZE(un, 2)
         SELECT CASE(comp)
         CASE(1)
            un(:, comp) = rho_anal(time, rr)
         CASE(2:k_dim+1)
            un(:, comp) = mt_anal(comp-1, time, rr)
         CASE(k_dim+2)
            un(:, comp) = E_anal(time, rr)
         CASE DEFAULT
            CALL error_petsc("BUG in init setup, wrong component "//to_str(comp)//&
                             " Max authorized is "//to_str(k_dim+2))
         END SELECT
      END DO
      CALL set_gamma_for_riemann_solver(gamma)
   END SUBROUTINE init

   FUNCTION rho_anal(time, rr) RESULT(vv)
     USE mesh_parameters
     IMPLICIT NONE
     REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
     REAL(KIND = 8), INTENT(IN) :: time
     REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
     INTEGER :: n, k
     REAL(KIND = 8) :: length, x
     IF (mesh_data_info%nb_bords==0) THEN
        length=1.d30
     ELSE
        length = abs(mesh_data_info%vect_e(1,1))
     END IF
     IF (SIZE(vv)==0) RETURN
     DO n = 1, SIZE(vv)
        k = floor((rr(1, n) - time)/length)
        x = rr(1, n) - time -k*length  
        IF (x<x0 .OR. x>x1) THEN
           vv(n) = 1.d0
        ELSE
           vv(n) = 1 + (2 / (x1 - x0))**6 * (x - x0)**3 * (x1 - x)**3
        END IF
     END DO
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
      IF (SIZE(vv)==0) RETURN
      vv = 1.d0
   END FUNCTION vit_anal

   FUNCTION E_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = press_anal(time, rr) / (gamma - 1.d0) &
           + rho_anal(time, rr) * (vit_anal(1, time,  rr)**2) / 2
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
         vv = E_anal(time, rr)
      CASE DEFAULT
         WRITE(*, *) ' BUG in sol_anal'
         STOP
      END SELECT
   END FUNCTION sol_anal
END MODULE setup
