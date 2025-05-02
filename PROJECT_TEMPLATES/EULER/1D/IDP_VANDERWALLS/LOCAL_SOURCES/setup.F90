MODULE setup
   USE space_dim
   USE eos
   USE vdw
   PUBLIC :: sol_anal, init, rho_anal, press_anal, mt_anal, E_anal, impose_bc_euler
   PRIVATE
   REAL(KIND = 8) :: x0, x1
   INTEGER :: VdW_test_case = 0
   REAL(KIND = 8) :: rhol, pl, ul
   REAL(KIND = 8) :: rhor, pr, ur
   REAL(KIND = 8) :: long
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

   SUBROUTINE init(un, time, rr)
      USE def_of_gamma
      USE lambda_module
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), k_dim + 2), INTENT(OUT) :: un
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8) :: cl, cr, rho_plus
      REAL(KIND = 8) :: in_state(2), in_data(3), out_state(4)
      SELECT CASE(VdW_test_case)
      CASE(0) !===Convergence test Section 6.1, SISC Vol. 44, No. 1, pp. A444-A470
         rho_plus = 0.35d0
         rhol = 0.1d0
         rhor = 0.39d0
         in_state(1) = rhol
         in_state(2) = rhor
         in_data(1) = 1.d0
         in_data(2) = 1.d0
         in_data(3) = 1.02d0
         CALL initialize_vdw(rho_plus, in_state, in_data, out_state)
         uL = out_state(1)
         uR = out_state(2)
         pL = out_state(3)
         pR = out_state(4)
      CASE(1) !===Stability test 1, Eq. (6.4) section 6.2, SISC Vol. 44, No. 1, pp. A444-A470
         avdw = 1.d0
         bvdw = 1.d0
         gamma_vdw = 1.02d0
         rhol = (0.5d0) * (2 - gamma_vdw) / (2 * bvdw)  !===rho must be smaller than (2-gamma_vdw)/(2*bvdw)
         rhor = (0.25d0) * (2 - gamma_vdw) / (2 * bvdw)
         pl = 1.01 * avdw * rhol**2 * (2 - gamma_vdw - 2 * bvdw * rhol) / (gamma_vdw)
         pr = 1.913 * avdw * rhor**2 * (2 - gamma_vdw - 2 * bvdw * rhor) / (gamma_vdw)
         ul = 0.
         ur = 0.
      CASE(2) !===Stability test 1, Eq. (6.5) section 6.2, SISC Vol. 44, No. 1, pp. A444-A470
         avdw = 1.d0
         bvdw = 1.d0
         gamma_vdw = 1.02d0
         rhol = 2.5d-1
         rhor = 4.9d-5
         ul = 0
         ur = 0
         pl = 3.d-2
         pr = 5.d-8
      CASE(3) !===Stability test 1, Eq. (6.6) section 6.2, SISC Vol. 44, No. 1, pp. A444-A470
         avdw = 1.d0
         bvdw = 1.d0
         gamma_vdw = 1.02d0
         rhol = 0.9932
         rhor = 0.95
         ul = 3
         ur = -3
         pl = 2
         pr = 2
      END SELECT
      WRITE(*, *) avdw, bvdw, gamma_vdw
      WRITE(*, *) 'rhol', rhol, 'rhor', rhor
      write(*, *) 'vl', ul, 'vr', ur
      write(*, *) 'pl', pl, 'pr', pr
      WRITE (*, *) 'cl', SQRT(gamma_vdw * (pl + avdw * rhol**2) / (rhol * (1 - bvdw * rhol)) - 2 * avdw * rhol)
      WRITE (*, *) 'cR', SQRT(gamma_vdw * (pr + avdw * rhor**2) / (rhor * (1 - bvdw * rhor)) - 2 * avdw * rhor)

      un(:, 1) = rho_anal(time, rr)
      un(:, 2) = mt_anal(1, time, rr)
      un(:, 3) = E_anal(time, rr)
   END SUBROUTINE init

   FUNCTION rho_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
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

   FUNCTION press_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
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

   FUNCTION vit_anal(comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
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

   FUNCTION E_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = (press_anal(time, rr) + avdw * rho_anal(time, rr)**2) * (1 - bvdw * rho_anal(time, rr)) / (gamma_vdw - 1.d0) &
           - avdw * rho_anal(time, rr)**2 + rho_anal(time, rr) * (vit_anal(1, time, rr)**2) / 2
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
