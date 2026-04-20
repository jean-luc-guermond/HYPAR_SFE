MODULE setup
   USE mesh_parameters
   PUBLIC :: sol_anal, init, rho_anal, press_anal, mt_anal, E_anal, impose_bc_euler, pressure
   PRIVATE
   REAL(KIND = 8) :: x0, x1
   REAL(KIND = 8) :: long
   REAL(KIND = 8), PARAMETER :: gamma = 1.4d0
   REAL(KIND = 8), PARAMETER :: rhoL=1.d0, rhor=0.125d0, pl=1.d0, pr=0.1d0, ul=0.d0, ur=0.d0,&
        l1m=-1.183215956619923d0, l1p=-0.07027281256118334d0, &
        l3=1.7521557320301779, ustar=0.92745262004894991d0, rhoLstar=0.4263194281784952d0, &
        rhoRstar=0.26557371170530708d0, pstar=0.3031301780506468, cL=SQRT(gamma*pL/rhoL) 
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
         ! SELECT CASE(comp)
         ! CASE(1)
         IF (comp == 1) THEN
            un(euler_bc%rho_bc%jsd, comp) = rho_anal(time, mesh%rr(:, euler_bc%rho_bc%jsd))
         ELSE IF ((2<=comp) .AND. (comp<=mesh_data_info%k_dim + 1)) THEN
         ! CASE(2:mesh_data%k_dim + 1)
            un(euler_bc%rho_bc%jsd, comp) = mt_anal(comp - 1, time, mesh%rr(:, euler_bc%rho_bc%jsd))
         ELSE IF (comp == mesh_data_info%k_dim + 2) THEN
         ! CASE(mesh_data%k_dim + 2)
            un(euler_bc%rho_bc%jsd, comp) = E_anal(time, mesh%rr(:, euler_bc%rho_bc%jsd))
         END IF
         ! END SELECT
      END DO

   END SUBROUTINE impose_bc_euler

   SUBROUTINE init(un, time, rr)
      USE def_of_gamma
      USE lambda_module
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), 3), INTENT(OUT) :: un
      ! REAL(KIND = 8), DIMENSION(SIZE(rr, 2), mesh_data%k_dim + 2), INTENT(OUT) :: un
      REAL(KIND = 8), INTENT(IN) :: time
      long = 1.d0
      x0 = long * 0.5d0
!!$      rhol = 1.d0
!!$      rhor = 0.125d0
!!$      pl = 1.d0
!!$      pr = 0.1d0
!!$      ul = 0.d0
!!$      ur = 0.d0
!!$
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
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (time.LE.1.d-12) THEN
            IF (rr(1, n)<x0) THEN
               vv(n) = rhol
            ELSE
               vv(n) = rhor
            END IF
         ELSE
            xi = (rr(1,n)-x0)/time
            IF (xi.LE.l1m) THEN
               vv(n) = rhoL
            ELSE IF (xi.LE.l1p) THEN
               vv(n) = rhoL*(2/(gamma+1) + (uL-xi)*(gamma-1)/((gamma+1)*cL))**(2/(gamma-1))
            ELSE IF (xi.LE.ustar) THEN
               vv(n) = rhoLstar
            ELSE IF (xi.LE.l3) THEN
               vv(n) = rhoRstar
            ELSE
               vv(n) = rhoR
            END IF
         END IF
      END DO
   END FUNCTION rho_anal

   FUNCTION press_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (time.LE.1.d-12) THEN
            IF (rr(1, n)<x0) THEN
               vv(n) = pl
            ELSE
               vv(n) = pr
            END IF
         ELSE
            xi = (rr(1,n)-x0)/time
            IF (xi.LE.l1m) THEN
               vv(n) = pL
            ELSE IF (xi.LE.l1p) THEN
               vv(n) = pL*(2/(gamma+1) + (uL-xi)*(gamma-1)/((gamma+1)*cL))**(2*gamma/(gamma-1))
            ELSE IF (xi.LE.l3) THEN
               vv(n) = pstar
            ELSE
               vv(n) = pR
            END IF
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
      REAL(KIND = 8) :: xi
      IF (SIZE(vv)==0) RETURN
      DO n = 1, SIZE(vv)
         IF (time.LE.1.d-12) THEN
            IF (rr(1, n)<x0) THEN
               vv(n) = ul
            ELSE
               vv(n) = ur
            END IF
         ELSE
            xi = (rr(1,n)-x0)/time
            IF (xi.LE.l1m) THEN
               vv(n) = uL
            ELSE IF (xi.LE.l1p) THEN
               vv(n) = (2/(gamma+1))*(cL + uL*(gamma-1)/2+xi)
            ELSE IF (xi.LE.l3) THEN
               vv(n) = ustar
            ELSE
               vv(n) = uR
            END IF
         END IF
      END DO
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
