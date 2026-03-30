MODULE euler_flux
   USE eos !<===FIX ME
   PUBLIC :: flux
CONTAINS
  FUNCTION flux(comp, un) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND = 8), DIMENSION(:, :), INTENT(IN) :: un
    INTEGER, INTENT(IN) :: comp
    REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: vv
    REAL(KIND = 8), DIMENSION(SIZE(un, 1)) :: u, ie
    SELECT CASE(comp)
    CASE(1)
       vv(:) = un(:,2)
    CASE(2)
       u = un(:,2)/un(:,1)
       vv(:) = un(:,2)*u
       ie = un(:,3)/un(:,1)
       ie = ie - 0.5d0*u**2
       vv(:) = vv(:) + pressure(un(:,1), ie)
    CASE(3)
       ie = un(:,3)/un(:,1)
       ie = ie - 0.5d0 *(un(:,2)/un(:,1))**2
       vv(:) = (un(:,2)/un(:,1)) * (un(:,3) + pressure(un(:,1), ie))
    CASE DEFAULT
       WRITE(*, *) ' BUG in flux'
       STOP
    END SELECT
  END FUNCTION flux
END MODULE euler_flux
