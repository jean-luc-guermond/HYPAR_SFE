MODULE setup
   USE space_dim
   USE eos
   PUBLIC :: sol_anal, init, rho_anal, press_anal, mt_anal, E_anal, impose_bc_euler
   PRIVATE
   REAL(KIND = 8) :: x0, x1
CONTAINS

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

   SUBROUTINE init(un, rr)
      USE def_of_gamma
      USE lambda_module
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), k_dim + 2), INTENT(OUT) :: un
      REAL(KIND = 8) :: cl, cr, rho_plus
      REAL(KIND = 8) :: in_state(2), in_data(3), out_state(4)

      gamma = 1.4d0
      long = 1.d0
      x0 = long * 0.5d0
      rhol = 1.d0
      rhor = 0.125d0
      pl = 1.d0
      pr = 0.1d0
      ul = 0.d0
      ur = 0.d0

      un(:, 1) = rho_anal(rr)
      un(:, 2) = mt_anal(1, rr)
      un(:, 3) = E_anal(rr)
      CALL set_gamma_for_riemann_solver(gamma)
   END SUBROUTINE init

   FUNCTION rho_anal(rr) RESULT(vv)
      USE input_data
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      REAL(KIND = 8) :: xi, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
      INTEGER :: n
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (rr(1, n)<x0) THEN
            vv(n) = rhol
         ELSE
            vv(n) = rhor
         END IF
      END DO
   END FUNCTION rho_anal

   FUNCTION press_anal(rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      REAL(KIND = 8) :: xi, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
      INTEGER :: n
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (rr(1, n)<x0) THEN
            vv(n) = pl
         ELSE
            vv(n) = pr
         END IF
      END DO
   END FUNCTION press_anal

   FUNCTION vit_anal(comp, rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      REAL(KIND = 8) :: xi, rhostarL, rhostarR, vstar, pstar, lambda1, lambda3
      INTEGER :: n
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (rr(1, n)<x0) THEN
            vv(n) = ul
         ELSE
            vv(n) = ur
         END IF
      END DO
   END FUNCTION vit_anal

   FUNCTION E_anal(rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = press_anal(rr) / (gamma - 1.d0) &
           + rho_anal(rr) * (vit_anal(1, rr)**2) / 2
   END FUNCTION E_anal

   FUNCTION mt_anal(comp, rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = rho_anal(rr) * vit_anal(comp, rr)
   END FUNCTION mt_anal

   FUNCTION pressure(un) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
      REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: vv
      vv = (un(:, 3) - 0.5d0 * (un(:, 2)**2) / un(:, 1)) * (gamma - 1)
   END FUNCTION pressure

   FUNCTION sol_anal(comp, rr, time) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      REAL(KIND = 8) :: time
      SELECT CASE(comp)
      CASE(1)
         vv = rho_anal(rr)
      CASE(2)
         vv = mt_anal(1, rr)
      CASE(3)
         vv = E_anal(rr)
      CASE DEFAULT
         WRITE(*, *) ' BUG in sol_anal'
         STOP
      END SELECT
   END FUNCTION sol_anal
END MODULE setup
