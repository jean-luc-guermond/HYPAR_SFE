MODULE setup
   USE space_dim, ONLY : k_dim
   USE euler_bc_arrays, ONLY: euler_bc_type, mt_anal_rho_times_vit

   PUBLIC :: pressure, init_state_functions

   PRIVATE
   REAL(KIND = 8) :: x0, x1
   INTEGER :: VdW_test_case = 0
   REAL(KIND = 8) :: rhol, pl, ul
   REAL(KIND = 8) :: rhor, pr, ur
   REAL(KIND = 8) :: long

CONTAINS

!==========================================================================
!================= DEF PRESSURE FOR SETUP =================================
!==========================================================================

   FUNCTION pressure(rho, e) RESULT(vv)
      USE vdw
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      vv = pressure_vdw(rho, e)
   END FUNCTION pressure

!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================

   SUBROUTINE init_state_functions(bc)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: bc
      CALL init_vdw

      bc%mt_anal    => mt_anal_rho_times_vit

      bc%E_anal     => E_anal_vdw
      bc%rho_anal   => rho_anal_vdw
      bc%press_anal => press_anal_vdw
      bc%vit_anal   => vit_anal_vdw
   END SUBROUTINE init_state_functions

   SUBROUTINE init_vdw
      USE vdw
      IMPLICIT NONE
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
      ! WRITE(*, *) avdw, bvdw, gamma_vdw
      ! WRITE(*, *) 'rhol', rhol, 'rhor', rhor
      ! write(*, *) 'vl', ul, 'vr', ur
      ! write(*, *) 'pl', pl, 'pr', pr
      ! WRITE (*, *) 'cl', SQRT(gamma_vdw * (pl + avdw * rhol**2) / (rhol * (1 - bvdw * rhol)) - 2 * avdw * rhol)
      ! WRITE (*, *) 'cR', SQRT(gamma_vdw * (pr + avdw * rhor**2) / (rhor * (1 - bvdw * rhor)) - 2 * avdw * rhor)
   END SUBROUTINE init_vdw

   FUNCTION rho_anal_vdw(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
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
   END FUNCTION rho_anal_vdw

   FUNCTION press_anal_vdw(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
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
   END FUNCTION press_anal_vdw

   FUNCTION vit_anal_vdw(this, comp, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
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
   END FUNCTION vit_anal_vdw

   FUNCTION E_anal_vdw(this, time, rr) RESULT(vv)
      USE vdw
      IMPLICIT NONE
      CLASS(euler_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      vv = (this%press_anal(time, rr) + avdw * this%rho_anal(time, rr)**2) * &
                (1 - bvdw * this%rho_anal(time, rr)) / (gamma_vdw - 1.d0) &
           - avdw * this%rho_anal(time, rr)**2 + this%rho_anal(time, rr) * (this%vit_anal(1, time, rr)**2) / 2
   END FUNCTION E_anal_vdw


END MODULE setup
