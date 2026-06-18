MODULE setup
   USE space_dim, ONLY : k_dim
   USE linear_transport_bc_arrays, ONLY: linear_transport_bc_type
   USE linear_transport_type_module, ONLY: linear_transport_type
   USE linear_transport_eta_commute, ONLY: default_eta_commute

   PUBLIC :: init_state_functions, my_linear_transport

   TYPE, EXTENDS(linear_transport_type) :: my_linear_transport
   CONTAINS
      PROCEDURE :: transport => sbr_transport 
   END TYPE my_linear_transport

   PRIVATE
   REAL(KIND=8), PARAMETER :: x0=0.4d0, a = 0.3d0
CONTAINS

!==========================================================================
!================= ANALYTICAL SOLUTIONS ===================================
!==========================================================================
   
   SUBROUTINE init_state_functions(linear_transport)
      IMPLICIT NONE
      CLASS(linear_transport_type), INTENT(INOUT) :: linear_transport
      linear_transport%bc%rho_anal => three_body
      ! linear_transport%bc%rho_anal => rho_circle
      ! linear_transport%bc%rho_anal   => sine_rho_anal
      ! linear_transport%bc%rho_anal   => step_rho_anal
      linear_transport%eta_commute   => eta_commute
      
   END SUBROUTINE init_state_functions

    FUNCTION eta_commute(un) RESULT(eta)
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(:, :), INTENT(IN)    :: un
        REAL(KIND=8), DIMENSION(SIZE(un,1))         :: eta
        eta = ABS(un(:,1)) + 1.d0
    END FUNCTION eta_commute

   FUNCTION rho_circle(this, time, rr) RESULT(vv)
      IMPLICIT NONE
      CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
      REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
      REAL(KIND = 8),                  INTENT(IN) :: time
      REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
      INTEGER :: n
      REAL(KIND=8) :: xt, yt, pi=ACOS(-1.d0)

      xt = COS(2*pi*time)*x0
      yt = SIN(2*pi*time)*x0

      DO n=1, SIZE(rr,2)
         vv(n) = 0.5d0*(1.d0-tanh(((rr(1,n)-xt)**2+(rr(2,n)-yt)**2)/a**2 - 1.d0))
      END DO
      !=== Dummy to avoid warning
      RETURN
      vv = SUM(this%rho_bc%jsd)*1.d0
      !=== Dummy to avoid warning
   END FUNCTION rho_circle

  FUNCTION three_body(this, time, rr) RESULT(vv)                                             
      IMPLICIT NONE                                                                
      CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
      REAL(KIND=8), DIMENSION(:,:),        INTENT(IN) :: rr                        
      REAL(KIND=8),                        INTENT(IN) :: time                   
      REAL(KIND=8), DIMENSION(SIZE(rr,2))             :: vv                        
      REAL(KIND=8), DIMENSION(2) :: r0 = (/0.5d0, 0.d0/), rc, rh, rd               
      REAL(KINd=8) ::  theta, radius, a, dc, dh, dd                                
      REAL(KIND=8) :: pi = ACOS(-1.d0)                                             
      INTEGER :: i                                                                 
                                                                                 
      !===Three body                                                               
      theta = 2*pi*time                                                            
      rh=(/-0.5d0,0.d0/)                                                           
      rc=(/0.d0,-0.5d0/)                                                           
      rd=(/0.d0,0.5d0/)                                                            
      a =.3d0                                                                      
      radius =0.5d0                                                                
      rh(1) = radius*COS(theta-pi/2)                                               
      rh(2) = radius*SIN(theta-pi/2)                                               
      rc(1) = radius*COS(theta+pi)                                                 
      rc(2) = radius*SIN(theta+pi)                                                 
      rd(1) = radius*COS(theta+pi/2)                                               
      rd(2) = radius*SIN(theta+pi/2)                                               
      DO i = 1, SIZE(rr,2)                                                         
         dh = SQRT((rr(1,i)-rh(1))**2+(rr(2,i)-rh(2))**2)                          
         dc = SQRT((rr(1,i)-rc(1))**2+(rr(2,i)-rc(2))**2)                          
         dd = SQRT((rr(1,i)-rd(1))**2+(rr(2,i)-rd(2))**2)                          
         IF (dh<a) THEN                                                            
            vv(i) = (1.d0+COS(pi*dh/a))/4                                          
         ELSE IF(dc<a) THEN                                                        
            vv(i) = 1.d0 - dc/a                                                    
         ELSE IF(dd<a .AND. (ABS(rr(1,i)).GE.0.05d0 .OR. rr(2,i).GE.0.7d0)) THEN   
            vv(i) = 1.d0                                                           
         ELSE                                                                      
            vv(i) = 0.d0                                                           
         END IF                                                                    
      END DO                                                                       
   END FUNCTION three_body

   ! FUNCTION step_rho_anal(this, time, rr) RESULT(vv)
   !    IMPLICIT NONE
   !    CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
   !    REAL(KIND = 8),                  INTENT(IN) :: time
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
   !    INTEGER :: n

   !    DO n=1, SIZE(rr,2)
   !       IF (rr(1, n)-cx*time < x0) THEN
   !          vv(n) = rhol
   !       ELSE
   !          vv(n) = rhor
   !       END IF
   !    END DO
      ! !=== Dummy to avoid warning
      ! RETURN
      ! vv = SUM(this%rho_bc%jsd)*1.d0
      ! !=== Dummy to avoid warning
   ! END FUNCTION step_rho_anal

   ! FUNCTION sine_rho_anal(this, time, rr) RESULT(vv)
   !    IMPLICIT NONE
   !    CLASS(linear_transport_bc_type), INTENT(INOUT) :: this
   !    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: rr
   !    REAL(KIND = 8),                  INTENT(IN) :: time
   !    REAL(KIND = 8), DIMENSION(SIZE(rr, 2))      :: vv
   !    REAL(KIND = 8),                   PARAMETER :: pi=ACOS(-1.d0)
   !    INTEGER                                     :: n

   !    DO n=1, SIZE(rr,2)
   !       vv(n) = SIN((rr(1,n)-cx*time)*2*pi) + 2.d0
   !    END DO
      ! !=== Dummy to avoid warning
      ! RETURN
      ! vv = SUM(this%rho_bc%jsd)*1.d0
      ! !=== Dummy to avoid warning
   ! END FUNCTION sine_rho_anal


   ! FUNCTION cst_transport(this, rr, comp) RESULT(vv)
   !    USE my_util, ONLY: error_petsc, to_str
   !    IMPLICIT NONE
   !    CLASS(my_linear_transport)                  :: this
   !    REAL(KIND = 8), DIMENSION(k_dim), INTENT(IN) :: rr
   !    INTEGER,                          INTENT(IN) :: comp
   !    REAL(KIND = 8)                               :: vv

   !    SELECT CASE(comp)
   !    CASE(1)
   !       vv = cx
   !    CASE(2)
   !       vv = cy
   !    CASE DEFAULT
   !       CALL error_petsc("BUG in cst_transport => wrong comp = "//to_str(comp))
   !    END SELECT
      ! !=== Dummy to avoid warning
      ! RETURN
      ! vv = SUM(this%rho_bc%jsd)*1.d0
      ! !=== Dummy to avoid warning
   ! END FUNCTION cst_transport

   FUNCTION sbr_transport(this, rr, comp) RESULT(vv)
      USE my_util, ONLY: error_petsc, to_str
      USE space_dim
      IMPLICIT NONE
      CLASS(my_linear_transport)                  :: this
      REAL(KIND = 8), DIMENSION(k_dim), INTENT(IN) :: rr
      INTEGER,                          INTENT(IN) :: comp
      REAL(KIND = 8)                               :: vv
      REAL(KIND = 8), PARAMETER :: pi=ACOS(-1.d0)

      SELECT CASE(comp)
      CASE(1)
         vv = -2*pi*rr(2)
      CASE(2)
         vv = +2*pi*rr(1)
      CASE DEFAULT
         CALL error_petsc("BUG in cst_transport => wrong comp = "//to_str(comp))
      END SELECT
      !=== Dummy to avoid warning
      RETURN
      vv = this%syst_dim*1.d0
      !=== Dummy to avoid warning
   END FUNCTION sbr_transport

END MODULE setup
