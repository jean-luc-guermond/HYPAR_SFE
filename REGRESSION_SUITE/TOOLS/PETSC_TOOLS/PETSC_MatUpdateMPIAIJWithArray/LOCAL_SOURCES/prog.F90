PROGRAM prog_petscSF

#include "petsc/finclude/petsc.h"
    USE petsc
    USE def_type_mesh
    USE read_inputs_module
    USE my_util
    USE st_matrix
    USE solver_petsc
    USE petsc_csr_LA_module

    IMPLICIT NONE
    TYPE(mesh_type)                    :: mesh
    TYPE(petsc_csr_LA)        :: LA
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: un_block
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: un_small_block
    REAL(KIND=8), DIMENSION(3,3) :: given_un_small_block
    REAL(KIND=8) :: tps_update, tps_set, tps_set_given, norm_mat1, norm_mat2, norm_diff
    INTEGER :: nw, ni, nj, m, i, j, n, rank, ierr
    INTEGER, DIMENSION(:), ALLOCATABLE :: idxm
    
    !=== TEST CASES
    !test=1 (CTEST) => test that results are the same with MatSetValues and MatUpdateMPIAIJWithArray
    !test=2         => test efficiency of MatUpdateMPIAIJWithArray VS old MatSetValues method 
    INTEGER, PARAMETER :: test=1
    !=== TEST CASES
    INTEGER :: n_max
    Mat :: matrix_1, matrix_2
    MPI_Comm :: communicator
    
    SELECT CASE(test)
    CASE(1)
        n_max = 1
    CASE(2)
        n_max = 1000
    CASE DEFAULT
        CALL error_petsc("BUG in select test => wrong value "//to_str(test))
    END SELECT

    !====== INITIALIZE PETSC AND MPI
    CALL start_setup
    rank = mesh%rank
    WRITE(*,'(A,I0,A,I0)') 'rank: ', rank, '; np: ', mesh%np
    nw = mesh%gauss%n_w

    !=== init fields
    ALLOCATE(un_block(mesh%me, nw, nw))
    ALLOCATE(un_small_block(nw, nw))
    ALLOCATE(idxm(nw))

    !=== preliminary First loop defining values inside un_block
        DO m=1, mesh%me
            DO ni=1, nw
                i = mesh%jj(ni, m)
                DO nj=1, nw
                    j = mesh%jj(nj, m)
                    un_block(m,nj,ni) = 10.*mesh%loc_to_glob(i) + mesh%loc_to_glob(j) !=== WARNING ROW ORIENTATION
                END DO
            END DO
        END DO
    !=== preliminary First loop defining values inside un_block

    tps_set = user_time()
    DO n=1, n_max    
    !=== Filling with MatSetValues approach
        CALL MatZeroEntries(matrix_2, ierr)
        DO m=1, mesh%me
            DO ni=1, nw
                i = mesh%jj(ni, m)
                DO nj=1, nw
                    j = mesh%jj(nj, m)
                    un_small_block(nj,ni) = 10.*mesh%loc_to_glob(i) + mesh%loc_to_glob(j) !=== WARNING ROW ORIENTATION
                END DO
            END DO
            idxm = LA%loc_to_glob(1, mesh%jj(:, m)) - 1
            CALL MatSetValues(matrix_2, nw, idxm, nw, idxm, un_small_block(:,:), ADD_VALUES, ierr)
        END DO

        CALL MatAssemblyBegin(matrix_2, MAT_FINAL_ASSEMBLY, ierr)
        CALL MatAssemblyEnd(matrix_2, MAT_FINAL_ASSEMBLY, ierr)
    END DO
    !=== MatSetValues approach
    tps_set = user_time() - tps_set


    tps_set_given = user_time()
    DO n=1, n_max    
    !=== Filling with MatSetValues approach
        CALL MatZeroEntries(matrix_2, ierr)
        DO m=1, mesh%me
            idxm = LA%loc_to_glob(1, mesh%jj(:, m)) - 1
            CALL MatSetValues(matrix_2, nw, idxm, nw, idxm, un_block(m,:,:), ADD_VALUES, ierr)
        END DO

        CALL MatAssemblyBegin(matrix_2, MAT_FINAL_ASSEMBLY, ierr)
        CALL MatAssemblyEnd(matrix_2, MAT_FINAL_ASSEMBLY, ierr)
    END DO
    !=== MatSetValues approach
    tps_set_given = user_time() - tps_set_given

    !======= CREATE petscSF STRUCTURE REQUIRED FOR COMMUNICATION
    CALL LA%init_mat_loc_to_glob(mesh)
    tps_update = user_time()
    DO n=1, n_max
        !=== Second loop defining un_contiguous
        !======= GATHER OPERATION ==========!
        ASSOCIATE(arr => LA%zz_contig_1, mat_loc_to_glob => LA%mat_loc_to_glob)
        arr(:) = 0.d0
        DO m=1, mesh%me
            DO ni=1, nw
                i = mesh%jj(ni, m)
                DO nj=1, nw
                    j = mesh%jj(nj, m)
                    arr(mat_loc_to_glob(ni, nj, m)) = arr(mat_loc_to_glob(ni, nj, m)) + 10.*mesh%loc_to_glob(i) + mesh%loc_to_glob(j) !=== WARNING ROW ORIENTATION
                END DO
            END DO
        END DO
        !======= COMMUNICATE & (instantanous) MATRIX FILL ==========!
        CALL LA%fill_mat(matrix_1, arr)
        END ASSOCIATE
        !=== Second loop defining un_contiguous

! !=================================!
! !=========== DEBUGGING ===========!
! !=================================!
!     ! DO m=1, mesh%me
!     !     DO ni=1, nw
!     !         WRITE(*,*) '(0) rank = ', mesh%rank, '; m, ni = ', m, ni, mesh%jj(ni, m), '; un = ', un_block(m,:,ni)
!     !     END DO
!     ! END DO

!     ! DO m=1, mesh%me
!     !     DO ni=1, nw
!     !         WRITE(*,*) '(I) rank = ', mesh%rank, '; m, ni = ', m, ni, mesh%jj(ni, m), '; un = ', un_contiguous(mat_loc_to_glob(m, ni, :))
!     !     END DO
!     ! END DO
!     ! WRITE(*,*) '(I) ghosts ', mesh%rank, un_contiguous(LA%ia(mesh%dom_np):)
!     ! WRITE(*,*) '(I) shift_ini ', mesh%rank, un_contiguous(shift_ini:)
!     ! WRITE(*,*) '(I) owned ', mesh%rank, un_contiguous(:nroots-1)
!     ! WRITE(*,*) '(I) ghosts ', mesh%rank, un_contiguous(nroots:)
! !=================================!
! !=========== DEBUGGING ===========!
! !=================================!

!         !=== Communicate ghosts
!         CALL PetscSFReduceBegin(sf_node, MPI_DOUBLE_PRECISION, un_contiguous(nroots: ), un_contiguous(0 : nroots - 1), MPI_SUM, ierr)
!         CALL PetscSFReduceEnd(sf_node, MPI_DOUBLE_PRECISION, un_contiguous(nroots: ), un_contiguous(0 : nroots - 1), MPI_SUM, ierr)
!         !=== Communicate ghosts

! !=================================!
! !=========== DEBUGGING ===========!
! !=================================!
!     ! DO m=1, mesh%me
!     !     DO ni=1, nw
!     !         WRITE(*,*) '(II) rank = ', mesh%rank, '; m, ni = ', m, ni, mesh%jj(ni, m), '; un = ', un_contiguous(mat_loc_to_glob(m, ni, :))
!     !     END DO
!     ! END DO
!     ! WRITE(*,*) '(II) ghosts ', mesh%rank, un_contiguous(LA%ia(mesh%dom_np):)
!     ! WRITE(*,*) '(II) shift_ini ', mesh%rank, un_contiguous(shift_ini:)
!     ! WRITE(*,*) '(II) owned ', mesh%rank, un_contiguous(:nroots-1)
!     ! WRITE(*,*) '(II) ghosts ', mesh%rank, un_contiguous(nroots:)
! !=================================!
! !=========== DEBUGGING ===========!
! !=================================!

    END DO
    tps_update = user_time() - tps_update


    SELECT CASE(test)
    CASE(1)
        CALL MatAXPY(matrix_1, -1.d0, matrix_2, SAME_NONZERO_PATTERN, ierr)
        CALL MatNorm(matrix_1, NORM_1, norm_diff, ierr)
        IF (norm_diff > 1.d-10) THEN
            IF (mesh%rank==0) WRITE(*,*) "test MatSetValues vs MatUpdateMPIAIJWithArray failed", norm_diff
        ELSE
            IF (mesh%rank==0) WRITE(*,*) "test MatSetValues vs MatUpdateMPIAIJWithArray successfull    ", '1234567891'
        END IF
    CASE(2)
        CALL MatNorm(matrix_1, NORM_1, norm_mat1, ierr)
        CALL MatNorm(matrix_2, NORM_1, norm_mat2, ierr)
        CALL MatAXPY(matrix_1, -1.d0, matrix_2, SAME_NONZERO_PATTERN, ierr)
        CALL MatNorm(matrix_1, NORM_1, norm_diff, ierr)

        IF (mesh%rank==0) THEN
            WRITE(*,*) 'tps_update = ', tps_update/n_max
            WRITE(*,*) 'tps_set = ', tps_set/n_max
            WRITE(*,*) 'tps_set_given = ', tps_set_given/n_max

            write(*,*) 'norm matrix_1 = ', norm_mat1
            write(*,*) 'norm matrix_2 = ', norm_mat2
            write(*,*) 'norm diff = ', norm_diff
        END IF
    END SELECT

CONTAINS

    SUBROUTINE start_setup
        USE construct_mesh,     ONLY: get_mesh
        USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
        USE options_module
        IMPLICIT NONE
        INTEGER :: rank
        PetscErrorCode :: ierr
        
        !===Start PETSC and MPI (mandatory)
        CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        communicator = PETSC_COMM_WORLD
        CALL MPI_Comm_rank(communicator, rank, ierr)

        !===Read executable arguments
        CALL read_all_arguments

        !===Clean data once
        CALL clean_data_once

        !===Construct mesh
        CALL get_mesh(communicator, mesh)

        !===Construct LA
        CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA)
        !===Construct matrix
        CALL create_local_petsc_matrix(communicator, LA, matrix_1, clean = .FALSE.)
        CALL create_local_petsc_matrix(communicator, LA, matrix_2, clean = .FALSE.)
    END SUBROUTINE start_setup

END PROGRAM prog_petscSF