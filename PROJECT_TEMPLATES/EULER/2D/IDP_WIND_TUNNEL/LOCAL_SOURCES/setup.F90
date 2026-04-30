MODULE setup
   USE space_dim, ONLY : k_dim
!VB TEST
   USE my_util, ONLY: WRITE_rank_0
!VB TEST
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
      IMPLICIT NONE
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
      IF (size(euler_bc%udotn_bc%jsd).NE.0) THEN
         mdotn = euler_bc%udotn_normal_vtx(:,1)*un(euler_bc%udotn_bc%jsd,2) &
         +  euler_bc%udotn_normal_vtx(:,2)*un(euler_bc%udotn_bc%jsd,3)
         un(euler_bc%udotn_bc%jsd,2) = un(euler_bc%udotn_bc%jsd,2) - mdotn*euler_bc%udotn_normal_vtx(:,1)
         un(euler_bc%udotn_bc%jsd,3) = un(euler_bc%udotn_bc%jsd,3) - mdotn*euler_bc%udotn_normal_vtx(:,2)
      END IF
      
   END SUBROUTINE impose_bc_euler
   
   SUBROUTINE init(un, time, rr)
      USE def_type_mesh
      USE def_of_gamma
      USE lambda_module
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2), k_dim + 2), INTENT(OUT) :: un
      REAL(KIND = 8), INTENT(IN) :: time
      un(:, 1) = rho_anal(time, rr)
      un(:, 2) = mt_anal(1, time, rr)
      un(:, 3) = mt_anal(2, time, rr)
      un(:, 4) = E_anal(time, rr)
      CALL set_gamma_for_riemann_solver(gamma)
   END SUBROUTINE init
   
   FUNCTION rho_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      IF (SIZE(vv)==0) RETURN
      vv = gamma
      ! CALL WRITE_rank_0("(IV) FLAG EULER_BC")

   END FUNCTION rho_anal
   
   FUNCTION press_anal(time, rr) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      IF (SIZE(vv)==0) RETURN
      vv = 1.d0
   END FUNCTION press_anal
   
   FUNCTION vit_anal(comp, time, rr, mesh) RESULT(vv)
      USE def_type_mesh
      USE petsc
#include "petsc/finclude/petsc.h"
      IMPLICIT NONE
      TYPE(mesh_type), optional :: mesh
      INTEGER, INTENT(IN) :: comp
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8), INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2)) :: vv
      INTEGER :: n, rank, ierr

      CALL MPI_Comm_Rank(PETSC_COMM_WORLD, rank, ierr)
      IF (SIZE(vv)==0) RETURN
      ! CALL WRITE_rank_0("(V) FLAG EULER_BC")
      vv = 0.d0
      IF (comp==1) THEN
         IF (time<1.d-8) THEN
            vv = 3.d0
            RETURN
         END IF

         WHERE(rr(1, :)<1.d-8)
            vv(:) = 3.d0
         ELSEWHERE
            vv(:) = 0.d0
         END WHERE
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
      USE def_type_mesh
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
      CASE(2:k_dim+1)
         vv = mt_anal(comp-1, time, rr)
      CASE(k_dim+2)
         vv = E_anal(time, rr)
      CASE DEFAULT
         WRITE(*, *) ' BUG in sol_anal, comp=', comp, 'should be <=', k_dim+2
         STOP
      END SELECT
   END FUNCTION sol_anal






!VB TEST

   SUBROUTINE write_l1_mesh(field_in, mesh, in_char, comp)
      USE fem_tn, ONLY : ns_l1
      USE def_type_mesh
      USE petsc
#include "petsc/finclude/petsc.h"

      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: field_in
      REAL(KIND=8) :: norm_loc, norm
      CHARACTER(LEN=*), INTENT(IN) :: in_char
      INTEGER, INTENT(IN) :: comp
      INTEGER :: ierr

      CALL ns_l1(mesh, field_in(:), norm_loc)
      CALL MPI_ALLREDUCE(norm_loc,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,PETSC_COMM_WORLD,ierr)
      IF(mesh%rank==0) WRITE(*, *) in_char, ' comp= ', comp, '=>  L1 Norm = ', norm

   END SUBROUTINE write_l1_mesh

!VB TEST


END MODULE setup
