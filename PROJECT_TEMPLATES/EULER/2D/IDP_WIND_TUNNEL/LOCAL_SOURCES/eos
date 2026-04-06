MODULE eos
   IMPLICIT NONE
   REAL(KIND=8), PARAMETER :: gamma=1.4d0
CONTAINS
   FUNCTION pressure(rho, e) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
      REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      vv = rho * e * (gamma - 1)
    END FUNCTION pressure

    FUNCTION sound_speed(un) RESULT(vv)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:),       INTENT(IN) :: un
      REAL(KIND=8), DIMENSION(SIZE(un,2))            :: rho, e, vv
      rho = un(1,:)
      e=(un(4,:)-0.5d0*(un(2,:)**2+un(3,:)**2)/rho)*(gamma-1) !===Dim = 2
      vv= SQRT(gamma*pressure(rho,e)/rho)
    END FUNCTION sound_speed
END MODULE eos

