MODULE setup
   USE mesh_parameters
   PUBLIC :: sol_anal, init, rho_anal, press_anal, mt_anal, E_anal, impose_bc_euler, init_eos_for_setup
   PRIVATE
   REAL(KIND = 8), PARAMETER :: x0 = 0.1d0, x1 = 0.3d0, gamma = 1.4d0
CONTAINS

!==========================================================================
!================= INIT EOS FOR SETUP  ====================================
!==========================================================================

   SUBROUTINE init_eos_for_setup
      USE eos_examples
      USE eos
      IMPLICIT NONE
      TYPE(eos_pointer_type) :: eos_type

      eos_type%pressure => pressure_ideal_diatomic_gas
      
      CALL assign_eos(eos_type)
   END SUBROUTINE init_eos_for_setup

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
         ! SELECT CASE(comp)
         ! CASE(1)
         IF (comp == 1) THEN
            un(euler_bc%rho_bc%jsd, comp) = rho_anal(time, mesh%rr(:, euler_bc%rho_bc%jsd))
         ! CASE(2:k_dim + 1)
         ELSE IF ((2<=comp) .AND. (comp<=mesh_data_info%k_dim + 1)) THEN
            un(euler_bc%rho_bc%jsd, comp) = mt_anal(comp - 1, time, mesh%rr(:, euler_bc%rho_bc%jsd))
         ! CASE(k_dim + 2)
         ELSE IF (comp == mesh_data_info%k_dim + 2) THEN
            un(euler_bc%rho_bc%jsd, comp) = E_anal(time, mesh%rr(:, euler_bc%rho_bc%jsd))
         ! END SELECT
         END IF
      END DO

   END SUBROUTINE impose_bc_euler

   SUBROUTINE init(un, time, rr)
      USE def_of_gamma
      USE lambda_module
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), mesh_data_info%k_dim + 2), INTENT(OUT) :: un
      ! REAL(KIND = 8), DIMENSION(SIZE(rr, 2), k_dim + 2), INTENT(OUT) :: un
      REAL(KIND = 8), INTENT(IN) :: time
      un(:, 1) = rho_anal(time, rr)
      un(:, 2) = mt_anal(1, time, rr)
      un(:, 3) = E_anal(time, rr)
      CALL set_gamma_for_riemann_solver(gamma)
   END SUBROUTINE init

   FUNCTION rho_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF ((rr(1, n) - time)<x0 .OR. (rr(1, n) - time)>x1) THEN
            vv(n) = 1.d0
         ELSE
            vv(n) = 1 + (2 / (x1 - x0))**6 * (rr(1, n) - time - x0)**3 * (x1 - rr(1, n) + time)**3
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
           + rho_anal(time, rr) * (vit_anal(1, time, rr)**2) / 2
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
