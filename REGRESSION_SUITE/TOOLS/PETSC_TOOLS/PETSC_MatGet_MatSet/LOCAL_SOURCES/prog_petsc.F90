PROGRAM test_matrix
#include "petsc/finclude/petsc.h"
    USE petsc

    INTEGER :: rank, ierr, ni, m, nj
    INTEGER, DIMENSION(2) :: idx, jdx
    INTEGER, DIMENSION(0:2) :: ia
    INTEGER, DIMENSION(0:3) :: ja
    REAL(KIND=8), DIMENSION(0:3) :: aa
    REAL(KIND=8), DIMENSION(2,2) :: block_1, block_2, mat_loc
    REAL(KIND=8), DIMENSION(2,2) :: mat_ref, mat_ref_transpose 
    LOGICAL :: test

    Mat :: mass_block, mass_seq

    !===== CTEST
    mat_ref(1,1) = 0.d0
    mat_ref(1,2) = 1.d0
    mat_ref(2,1) = 0.d0
    mat_ref(2,2) = 1.d0

    mat_ref_transpose = TRANSPOSE(mat_ref)

    !===== INIT
    ia(0) = 0
    ia(1) = 2
    ia(2) = 4
    ja(0) = 0
    ja(1) = 1
    ja(2) = 0
    ja(3) = 1
    aa = 0.d0

    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

    CALL MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD,2,2,PETSC_DECIDE, &
         PETSC_DECIDE, ia, ja, aa, mass_block, ierr)
    CALL MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD,2,2,PETSC_DECIDE, &
         PETSC_DECIDE, ia, ja, aa, mass_seq, ierr)
    

    !===== Fill either in block or with loops
    DO ni=1, 2
        idx(ni) = ni - 1
    END DO   
    !===== Fill with loops
    DO ni=1, 2
        DO nj=1, 2
           CALL MatSetValues(mass_seq, 1, idx(ni:ni), 1, idx(nj:nj), idx(nj:nj)*1.d0, INSERT_VALUES, ierr)
        END DO
    END DO 
    CALL MatAssemblyBegin(mass_seq, MAT_FINAL_ASSEMBLY, ierr)
    CALL MatAssemblyEnd  (mass_seq, MAT_FINAL_ASSEMBLY, ierr)   

    !===== Fill in block
    DO ni=1, 2
        mat_loc(ni, :) = idx(:)*1.d0
    END DO
    CALL MatSetValues(mass_block, 2, idx(:), 2, idx(:), mat_loc, INSERT_VALUES, ierr)
    CALL MatAssemblyBegin(mass_block, MAT_FINAL_ASSEMBLY, ierr)
    CALL MatAssemblyEnd  (mass_block, MAT_FINAL_ASSEMBLY, ierr)


    test = .TRUE.

   !=== Extract results from mass_seq
    block_1 = -1.d0
    block_2 = -1.d0
    !=== Extract in block
    CALL MatGetValues(mass_seq, 2, idx(:), 2, idx(:), block_1, ierr)
    !=== Extract by loops
    DO ni=1, 2
        DO nj=1, 2
            CALL MatGetValues(mass_seq, 1, idx(ni:ni), 1, idx(nj:nj), block_2(ni:ni, nj:nj), ierr)
        END DO
    END DO  
    CALL MatAssemblyBegin(mass_seq, MAT_FINAL_ASSEMBLY, ierr)
    CALL MatAssemblyEnd  (mass_seq, MAT_FINAL_ASSEMBLY, ierr)
    !=== Write results from mass_seq
    IF (MAXVAL(ABS(mat_ref_transpose - block_1)) > 1.d-10) THEN
        test = .FALSE.
        WRITE(*,*) "test 1 wrong"
    END IF

    IF (MAXVAL(ABS(mat_ref - block_2)) > 1.d-10) THEN
        test = .FALSE.
        WRITE(*,*) "test 2 wrong"
    END IF
    ! WRITE(*,*) "======== Extraction from MASS_SEQ ========"
    ! WRITE(*,*) "==== Extraction in block (is wrongly transposed) ==="
    ! DO ni=1, 2
    !     WRITE(*,*) 'row ', ni, block_1(ni, :)
    ! END DO
    ! WRITE(*,*) "==== Extraction by loops (is correct) ==="
    ! DO ni=1, 2
    !     WRITE(*,*) 'row ', ni, block_2(ni, :)
    ! END DO

   !=== Extract results from mass_block
    block_1 = -1.d0
    block_2 = -1.d0  
    !=== Extract in block
    CALL MatGetValues(mass_block, 2, idx(:), 2, idx(:), block_1, ierr)
    !=== Extract by loops
    DO ni=1, 2
        DO nj=1, 2
            CALL MatGetValues(mass_block, 1, idx(ni:ni), 1, idx(nj:nj), block_2(ni:ni, nj:nj), ierr)
        END DO
    END DO  
    CALL MatAssemblyBegin(mass_block, MAT_FINAL_ASSEMBLY, ierr)
    CALL MatAssemblyEnd  (mass_block, MAT_FINAL_ASSEMBLY, ierr)
    !=== Write results from mass_block
    IF (MAXVAL(ABS(mat_ref - block_1)) > 1.d-10) THEN
        test = .FALSE.
        WRITE(*,*) "test 3 wrong"
    END IF

    IF (MAXVAL(ABS(mat_ref_transpose - block_2)) > 1.d-10) THEN
        test = .FALSE.
        WRITE(*,*) "test 4 wrong"
    END IF
    ! WRITE(*,*)
    ! WRITE(*,*)
    ! WRITE(*,*) "======== Extraction from MASS_BLOCK ========"
    ! WRITE(*,*) "==== Extraction in block (is correct) ==="    
    ! DO ni=1, 2 
    !     WRITE(*,*) 'row ', ni, block_1(ni, :)
    ! END DO
    ! WRITE(*,*) "==== Extraction by loops (is wrongly transposed) ==="
    ! DO ni=1, 2
    !     WRITE(*,*) 'row ', ni, block_2(ni, :)
    ! END DO

    IF (test) THEN
        WRITE(*,*) 'test passed ', '1234567891'
    ELSE
        WRITE(*,*) 'test failed'
    END IF

END PROGRAM test_matrix