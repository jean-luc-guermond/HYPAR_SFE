PROGRAM prog_petscSF

#include "petsc/finclude/petsc.h"
    USE petsc
    USE def_type_mesh
    USE read_inputs_module
    USE my_util
    USE st_matrix

    IMPLICIT NONE
    TYPE(mesh_type)                    :: mesh
    TYPE(petsc_csr_LA)                 :: LA
    TYPE(periodic_type), DIMENSION(1)  :: per
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: un, un_extr, un_sf
    INTEGER :: n, rank, ierr, k
    INTEGER, DIMENSION(:), POINTER :: ifrom

    !=== TEST CASES
    !test=1         => 1D small scale test to understand the communication of ghost values
    !test=2 (CTEST) => test the min/max communication of ghost values with petscSF VS old extraction method
    !test=3         => test efficiency of the min/max communication of ghost values with petscSF VS old extraction method 
    INTEGER, PARAMETER :: test=2
    INTEGER :: n_test_2 = 20
    LOGICAL :: test_passed=.TRUE.
    !=== TEST CASES
    INTEGER, PARAMETER :: n_max=10000
    REAL(KIND=8) :: tps
    MPI_Comm :: communicator
    Vec :: xx
    !====== INITIALIZE PETSC AND MPI
    CALL start_setup
    rank = mesh%rank
    ALLOCATE(un(mesh%np))

    WRITE(*,'(A,I0,A,I0)') 'rank: ', rank, '; np: ', mesh%np

    SELECT CASE (test)
    CASE(1)
        IF (mesh%disp(mesh%nb_proc+1) /= 22) THEN
            CALL error_petsc("BUG in prog_petscSF, please use np = 10 with 1 refinement level and PBC")
        END IF
        !=== INITIALIZE un AND WRITE OUTPUT
        CALL init_un(un, mesh)
        CALL write_output(un, mesh, "INITIAL")

        !=== TEST COMMUNICATIONS
        CALL init_un(un, mesh)
        CALL mesh%comm_ghost(un, MPI_SUM)
        CALL write_output(un, mesh, "SUM")

        CALL init_un(un, mesh)
        CALL mesh%comm_ghost(un, MPI_MAX)
        CALL write_output(un, mesh, "MAX")

        CALL init_un(un, mesh)
        CALL mesh%comm_ghost(un, MPI_MIN)
        CALL write_output(un, mesh, "MIN")

        CALL init_un(un, mesh)
        CALL mesh%comm_ghost(un, MPI_PROD)
        CALL write_output(un, mesh, "PROD")
    CASE(2)
        !=== initialize vectors
        ALLOCATE(un_extr(mesh%np), un_sf(mesh%np))
        CALL create_my_ghost(mesh, LA, ifrom)
        CALL VecCreateGhost(mesh%comm, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, xx, ierr)
        !=== initialize vectors

        !=== TEST COMMUNICATIONS several times
        DO k=1, n_test_2
            CALL RANDOM_NUMBER(un)
            CALL array_to_petsc_vec_min_max(un, xx, LA, 'min')
            CALL extract_through_ghost(xx, 1, 1, LA, un_extr, opt_assemble=.FALSE.)
            un_sf = un
            CALL mesh%comm_ghost(un_sf, MPI_MIN)
            IF (SUM(ABS(un_extr-un_sf)) > 1.d-10) THEN
                test_passed = .FALSE.
                ! CALL error_petsc("BUG in prog_petscSF, un_extr and un_sf not consistent")
            END IF

            CALL RANDOM_NUMBER(un)
            CALL array_to_petsc_vec_min_max(un, xx, LA, 'max')
            CALL extract_through_ghost(xx, 1, 1, LA, un_extr, opt_assemble=.FALSE.)
            un_sf = un
            CALL mesh%comm_ghost(un_sf, MPI_MAX)
            IF (SUM(ABS(un_extr-un_sf)) > 1.d-10) THEN
                test_passed = .FALSE.
                ! CALL error_petsc("BUG in prog_petscSF, un_extr and un_sf not consistent")
            END IF
        END DO
        !=== TEST COMMUNICATIONS several times

        IF (test_passed) THEN
            IF (mesh%rank==0) WRITE(*,*) "TEST PETSC_SF PASSED", 1234567891
        ELSE
            IF (mesh%rank==0) WRITE(*,*) "TEST PETSC_SF FAILED"
        END IF

    CASE(3)
        ALLOCATE(un_extr(mesh%np), un_sf(mesh%np))
        CALL create_my_ghost(mesh, LA, ifrom)
        CALL VecCreateGhost(mesh%comm, mesh%dom_np, &
            PETSC_DETERMINE, SIZE(ifrom), ifrom, xx, ierr)
        CALL RANDOM_NUMBER(un)
        tps = user_time()
        DO n=1, n_max
            CALL array_to_petsc_vec_min_max(un, xx, LA, 'min')
            CALL extract_through_ghost(xx, 1, 1, LA, un_extr, opt_assemble=.FALSE.)
        END DO
        tps = user_time() - tps
        IF (mesh%rank==0) WRITE(*,*) 'Time for array_to_petsc_vec + extract_through_ghost: ', tps/n_max, ' s', ' n_max: ', n_max, " np: ", mesh%np, " dom_np: ", mesh%dom_np, "ghost_np: ", mesh%np-mesh%dom_np
        tps = user_time()
        DO n=1, n_max
            un_sf = un
            CALL mesh%comm_ghost(un_sf, MPI_MIN)
        END DO
        tps = user_time() - tps
        IF (mesh%rank==0) WRITE(*,*) 'Time for mesh%comm_ghost: ', tps/n_max, ' s', 'n_max: ', n_max
    CASE DEFAULT
        WRITE(*,*) ' BUG in prog_petscSF, test not implemented: ', test
        STOP
    END SELECT


CONTAINS

    SUBROUTINE start_setup
        USE construct_mesh,     ONLY: get_mesh
        USE st_matrix,          ONLY: st_aij_csr_glob_block_with_extra_layer
        IMPLICIT NONE
        PetscErrorCode :: ierr
        INTEGER :: rank

        !===Start PETSC and MPI (mandatory)

        CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        communicator = PETSC_COMM_WORLD
        CALL MPI_Comm_rank(communicator, rank, ierr)

        !===Clean data once
        CALL clean_data_once

        !===Construct mesh
        CALL get_mesh(communicator, mesh)

        !===Construct LA
        CALL st_aij_csr_glob_block_with_extra_layer(communicator, 1, mesh, LA)
    END SUBROUTINE start_setup

    SUBROUTINE init_un(un, mesh)
        USE def_type_mesh, ONLY: mesh_type
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN) :: mesh
        REAL(KIND=8), DIMENSION(mesh%np), INTENT(OUT) :: un
        INTEGER :: n

        un = 1.d0*mesh%proc
        DO n=mesh%dom_np, mesh%np
            un(n) = n
        END DO
    END SUBROUTINE init_un

    SUBROUTINE write_output(un, mesh, char)
        USE def_type_mesh, ONLY: mesh_type
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN) :: mesh
        REAL(KIND=8), DIMENSION(:), INTENT(IN) :: un
        INTEGER :: n, rank
        CHARACTER(LEN=*), INTENT(IN) :: char

        rank = mesh%rank
        IF (mesh%rank==0) THEN
            WRITE(*,*) "================= un "//char//" ================="
            WRITE(*,*) ""
        END IF
        DO n=1, mesh%nb_proc
            IF (n==mesh%rank+1) THEN
                DO k=1, mesh%np
                    WRITE(*,*) 'rank: ', rank, mesh%loc_to_glob(k), mesh%rr(1,k), '; CAST un=', un(k)
                END DO
            END IF
            CALL MPI_BARRIER(communicator, ierr)
        END DO
    END SUBROUTINE write_output

    SUBROUTINE array_to_petsc_vec_min_max(uu, xx, LA, operation)
        !> convert uu(HYPAR vector) to xx(petsc vec)
        !! LA
        !! operation = 'insert'/'add'/'min'/'max' (for ghost values)
        USE my_util
        IMPLICIT NONE
        TYPE(petsc_csr_LA),           INTENT(IN) :: LA
        REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: uu
        REAL(KIND = 8), DIMENSION(:), POINTER    :: x_loc
        CHARACTER(LEN = *),           INTENT(IN) :: operation
        INTEGER, DIMENSION(SIZE(uu)) :: idxm
        INTEGER :: i, np, k, kmax, ierr
        Vec     :: xx, xx_ghost
        SELECT CASE (operation)
        CASE('min')
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
            CALL error_petsc("BUG in array_to_petsc_vec_min, operation not implemented: "//operation)
        END SELECT

    END SUBROUTINE array_to_petsc_vec_min_max


END PROGRAM prog_petscSF