MODULE eos_interface
   ABSTRACT INTERFACE
      FUNCTION pressure_template(rho, e) RESULT(vv)
        IMPLICIT NONE
        REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
        REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
      END FUNCTION pressure_template
   END INTERFACE
END MODULE eos_interface
   
MODULE eos
    USE eos_interface
    IMPLICIT NONE

    TYPE eos_pointer_type
       PROCEDURE(pressure_template), NOPASS, POINTER :: pressure => null()
    END TYPE eos_pointer_type

    PROCEDURE(pressure_template), POINTER :: pressure => null()

  CONTAINS 
    
    SUBROUTINE assign_eos(eos_type)
        IMPLICIT NONE
        TYPE(eos_pointer_type), INTENT(IN) :: eos_type

        !=== Associate pressure pointer
        IF (ASSOCIATED(eos_type%pressure)) THEN
           pressure => eos_type%pressure
        ELSE
            WRITE(*,*) "ERROR in assign_eos: pressure pointer not associated in eos_type"
            STOP
        END IF

    END SUBROUTINE assign_eos

END MODULE eos

MODULE eos_examples
    USE eos_interface
    IMPLICIT NONE

  CONTAINS

    FUNCTION pressure_ideal_diatomic_gas(rho, e) RESULT(vv)
       IMPLICIT NONE
       REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
       REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
       REAL(KIND = 8) :: gamma
       gamma = 7.0 / 5.0
       vv = pressure_ideal_gas(rho, e, gamma)
    END FUNCTION pressure_ideal_diatomic_gas

    FUNCTION pressure_ideal_monoatomic_gas(rho, e) RESULT(vv)
       IMPLICIT NONE
       REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
       REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
       REAL(KIND = 8) :: gamma
       gamma = 5.0 / 3.0
       vv = pressure_ideal_gas(rho, e, gamma)
    END FUNCTION pressure_ideal_monoatomic_gas

    FUNCTION pressure_ideal_gas(rho, e, gamma) RESULT(vv)
       IMPLICIT NONE
       REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: rho, e
       REAL(KIND = 8), INTENT(IN) :: gamma
       REAL(KIND = 8), DIMENSION(SIZE(rho)) :: vv
       vv = rho * e * (gamma - 1)
    END FUNCTION pressure_ideal_gas

END MODULE eos_examples