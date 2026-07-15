MODULE petsc_tools
#include "petsc/finclude/petsc.h"
   USE petsc
   USE def_type_mesh
   USE petsc_csr_LA_module

CONTAINS

   ! SUBROUTINE array_to_petsc_vec(uu, xx, LA, operation)!, opt_include_ghost)
   !    !> convert uu(HYPAR vector) to xx(petsc vec)
   !    !! LA
   !    !! operation = 'insert'/'add'/'min'/'max' (for ghost values)
   !    USE my_util
   !    IMPLICIT NONE
   !    TYPE(petsc_csr_LA),           INTENT(IN) :: LA
   !    REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: uu
   !    REAL(KIND = 8), DIMENSION(:), POINTER    :: x_loc
   !    CHARACTER(LEN = *),           INTENT(IN) :: operation
   !    INTEGER, DIMENSION(SIZE(uu)) :: idxm
   !    INTEGER :: i, np, k, kmax, ierr
   !    ! LOGICAL, OPTIONAL :: opt_include_ghost
   !    Vec     :: xx, xx_ghost

   !    kmax = SIZE(LA%loc_to_glob,1)
   !    np = SIZE(uu)/kmax

   !    IF (np /= SIZE(LA%loc_to_glob, 2)) THEN
   !       CALL error_petsc("BUG in array to petsc, wrong np")
   !    END IF

   !    DO k=1, kmax
   !       DO i = 1, np
   !          idxm(i+(k-1)*np) = LA%loc_to_glob(k, i) - 1
   !       END DO
   !    END DO
   !    SELECT CASE (operation)
   !    CASE('insert')
   !       CALL VecSetValues(xx, np*kmax, idxm, uu, INSERT_VALUES, ierr)
   !    CASE('add')
   !       CALL VecSetValues(xx, np*kmax, idxm, uu, ADD_VALUES, ierr)
   !    CASE('min')
   !       !=== Define ghost vectors LOCALLY
   !       CALL VecGhostGetLocalForm(xx, xx_ghost, ierr)
   !       CALL VecGetArrayF90(xx_ghost, x_loc, ierr)
   !       x_loc(:) = uu(:)
   !       CALL VecRestoreArrayF90(xx_ghost, x_loc, ierr)
   !       CALL VecGhostRestoreLocalForm(xx, xx_ghost, ierr)

   !       !=== At global level, computing min values on Ghost points across processes and:
   !          !=== 1) Update dom_np+1:np (compute min between overlapping values, update ghost values)
   !       CALL VecGhostUpdateBegin(xx, MIN_VALUES, SCATTER_REVERSE, ierr)
   !       CALL VecGhostUpdateEnd(xx, MIN_VALUES, SCATTER_REVERSE, ierr)
   !          !=== 2) Update 1:dom_np (update owned dofs with ghost values computed with min)
   !       CALL VecGhostUpdateBegin(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)
   !       CALL VecGhostUpdateEnd(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)

   !    CASE('max')
   !       !=== Define ghost vectors LOCALLY
   !       CALL VecGhostGetLocalForm(xx, xx_ghost, ierr)
   !       CALL VecGetArrayF90(xx_ghost, x_loc, ierr)
   !       x_loc(:) = uu(:)
   !       CALL VecRestoreArrayF90(xx_ghost, x_loc, ierr)
   !       CALL VecGhostRestoreLocalForm(xx, xx_ghost, ierr)
   !       !=== At global level, computing max values on Ghost points across processes and:
   !          !=== 1) Update dom_np+1:np (compute max between overlapping values, update ghost values)
   !       CALL VecGhostUpdateBegin(xx, MAX_VALUES, SCATTER_REVERSE, ierr)
   !       CALL VecGhostUpdateEnd(xx, MAX_VALUES, SCATTER_REVERSE, ierr)
   !          !=== 2) Update 1:dom_np (update owned dofs with ghost values computed with max)
   !       CALL VecGhostUpdateBegin(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)
   !       CALL VecGhostUpdateEnd(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)


   !    CASE DEFAULT
   !       CALL error_petsc(&
   !       &'Wrong option operation '//operation//' in array_to_petsc_vec ==> should be: add/insert/min/max')
   !    END SELECT

   !    CALL VecAssemblyBegin(xx, ierr)
   !    CALL VecAssemblyEnd(xx, ierr)
   ! END SUBROUTINE array_to_petsc_vec

   SUBROUTINE array_to_petsc_vec(uu, xx, LA, operation, opt_include_ghost)
      !> convert uu(HYPAR vector) to xx(petsc vec)
      !! LA
      !! operation = 'insert'/'add'/'min'/'max' (for ghost values)
      USE my_util, ONLY: error_petsc, pack_opt, to_str
      IMPLICIT NONE
      TYPE(petsc_csr_LA),           INTENT(IN) :: LA
      REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: uu
      REAL(KIND = 8), DIMENSION(:), POINTER    :: x_loc
      CHARACTER(LEN = *),           INTENT(IN) :: operation
      ! INTEGER,      DIMENSION(SIZE(uu)) :: idxm
      ! REAL(KIND=8), DIMENSION(SIZE(uu)) :: values
      INTEGER,      DIMENSION(:), ALLOCATABLE :: idxm
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: values
      INTEGER :: np, dom_np, k, kmax, ierr
      INTEGER :: np_eff
      LOGICAL, OPTIONAL :: opt_include_ghost
      LOGICAL :: include_ghost
      Vec     :: xx, xx_ghost

      kmax = SIZE(LA%loc_to_glob,1)
      np = SIZE(uu)/kmax
      dom_np = LA%dom_np(1)
      IF (MOD(SIZE(uu), kmax) /= 0) THEN
         CALL error_petsc("BUG in array to petsc: SIZE(uu)="//to_str(SIZE(uu))//" incompatible with kmax="//to_str(kmax))
      ELSE IF (np /= SIZE(LA%loc_to_glob, 2)) THEN
         CALL error_petsc("BUG in array to petsc: np="//to_str(np)//" incompatible with LA%loc_to_glob="//to_str(SIZE(LA%loc_to_glob, 2)))
      ELSE IF (np < dom_np) THEN
         CALL error_petsc("BUG in array to petsc: np<dom_np ("//to_str(np)//"<"//to_str(dom_np)//"); how is that possible??")
      END IF

      ! !=== By default, fortran => petsc operation operates on ghosts as well
      CALL pack_opt(include_ghost, .TRUE., opt_include_ghost)
      IF (include_ghost) THEN
         np_eff = np
      ELSE
         np_eff = dom_np
      END IF
      ! np_eff = np
      ALLOCATE(idxm(kmax*np_eff))
      ALLOCATE(values(kmax*np_eff))

      DO k=1, kmax
         values((k-1)*np_eff+1 : k*np_eff) = uu((k-1)*np + 1 : (k-1)*np + np_eff)
         idxm((k-1)*np_eff+1 : k*np_eff) = LA%loc_to_glob(k, 1:np_eff) - 1
      END DO
      ! write(*,*) sum(values), sum(idxm)
      ! values = uu
      ! DO k=1, kmax
      !    idxm((k-1)*np_eff+1 : k*np_eff) = LA%loc_to_glob(k, 1:np_eff) - 1
      ! END DO

      SELECT CASE (operation)
      CASE('insert')
         CALL VecSetValues(xx, np_eff*kmax, idxm, values, INSERT_VALUES, ierr)
      CASE('add')
         CALL VecSetValues(xx, np_eff*kmax, idxm, values, ADD_VALUES, ierr)
      CASE('min')
         IF (.NOT. include_ghost) THEN
            CALL error_petsc("BUG in array_to_petsc_vec: min incompatible with .NOT. include_ghost")
         END IF
         !=== Define ghost vectors LOCALLY
         CALL VecGhostGetLocalForm(xx, xx_ghost, ierr)
         CALL VecGetArrayF90(xx_ghost, x_loc, ierr)
         x_loc(:) = uu(:)
         CALL VecRestoreArrayF90(xx_ghost, x_loc, ierr)
         CALL VecGhostRestoreLocalForm(xx, xx_ghost, ierr)

         !=== At global level, computing min values on Ghost points across processes and:
            !=== 1) Update dom_np+1:np (compute min between overlapping values, update ghost values)
         CALL VecGhostUpdateBegin(xx, MIN_VALUES, SCATTER_REVERSE, ierr)
         CALL VecGhostUpdateEnd(xx, MIN_VALUES, SCATTER_REVERSE, ierr)
            !=== 2) Update 1:dom_np (update owned dofs with ghost values computed with min)
         CALL VecGhostUpdateBegin(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)

      CASE('max')
         IF (.NOT. include_ghost) THEN
            CALL error_petsc("BUG in array_to_petsc_vec: max incompatible with .NOT. include_ghost")
         END IF
         !=== Define ghost vectors LOCALLY
         CALL VecGhostGetLocalForm(xx, xx_ghost, ierr)
         CALL VecGetArrayF90(xx_ghost, x_loc, ierr)
         x_loc(:) = uu(:)
         CALL VecRestoreArrayF90(xx_ghost, x_loc, ierr)
         CALL VecGhostRestoreLocalForm(xx, xx_ghost, ierr)
         !=== At global level, computing max values on Ghost points across processes and:
            !=== 1) Update dom_np+1:np (compute max between overlapping values, update ghost values)
         CALL VecGhostUpdateBegin(xx, MAX_VALUES, SCATTER_REVERSE, ierr)
         CALL VecGhostUpdateEnd(xx, MAX_VALUES, SCATTER_REVERSE, ierr)
            !=== 2) Update 1:dom_np (update owned dofs with ghost values computed with max)
         CALL VecGhostUpdateBegin(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)
         CALL VecGhostUpdateEnd(xx, INSERT_VALUES, SCATTER_FORWARD, ierr)


      CASE DEFAULT
         CALL error_petsc(&
         &'Wrong option operation '//operation//' in array_to_petsc_vec ==> should be: add/insert/min/max')
      END SELECT
      CALL VecAssemblyBegin(xx, ierr)
      CALL VecAssemblyEnd(xx, ierr)
      DEALLOCATE(idxm, values)

   END SUBROUTINE array_to_petsc_vec

END MODULE petsc_tools