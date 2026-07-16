MODULE petsc_csr_LA_module
    !> Update VB 01/07/2026
    !! module containing object petsc_csr_LA
    !! fill_mat: offers a faster way to fill petsc matrices than with MatSetValues
    !! WARNING: a matrix provided by fill_mat must have been generated with petsc_create_local_matrix, not MatDuplicate!!
    !! TODO => generalize to LA(dim>1)
    !!
    !!    MatSetValues:
    !! 
    !!    CALL MatZeroEntries(matrix_2, ierr)
    !!    DO m=1, mesh%me
    !!        idxm = LA%loc_to_glob(1, mesh%jj(:, m)) - 1
    !!        DO ni=1, nw
    !!            DO nj=1, nw
    !!                CALL my_compute(un_block(nj, ni), m, ni, nj)  !=== WARNING ROW ORIENTATION
    !!            END DO
    !!        END DO
    !!        CALL MatSetValues(matrix_2, nw, idxm, nw, idxm, un_block(:,:), ADD_VALUES, ierr)
    !!    END DO
    !!    CALL MatAssemblyBegin(matrix_2, MAT_FINAL_ASSEMBLY, ierr)
    !!    CALL MatAssemblyEnd(matrix_2, MAT_FINAL_ASSEMBLY, ierr)
    !!
    !!  
    !!    LA%fill_mat:
    !!
    !!    !======= GATHER OPERATION ==========!
    !!    ASSOCIATE(arr => LA%zz_contig_1, mat_loc_to_glob => LA%mat_loc_to_glob)
    !!    arr(:) = 0.d0
    !!    DO m=1, mesh%me
    !!        DO ni=1, nw
    !!            DO nj=1, nw
    !!                CALL my_compute(un_block(nj, ni), m, ni, nj)
    !!                arr(mat_loc_to_glob(ni, nj, m)) = arr(mat_loc_to_glob(ni, nj, m)) + un_block(nj, ni) !=== WARNING ROW ORIENTATION
    !!            END DO
    !!        END DO
    !!    END DO
    !!    !======= COMMUNICATE & (instantanous) MATRIX FILL ==========!
    !!    CALL LA%fill_mat(matrix_1, arr)
    !!    END ASSOCIATE

    USE def_type_mesh
    USE petscmat
    USE petscvec

    TYPE :: petsc_csr_LA
    !=== Elementary objects
        INTEGER, DIMENSION(:), POINTER :: ia, ja
        INTEGER, DIMENSION(:, :), POINTER :: loc_to_glob
        INTEGER :: kmax
        INTEGER, DIMENSION(:), POINTER :: np
        INTEGER, DIMENSION(:), POINTER :: dom_np
    !=== Added objects & procedures for performance
        TYPE(tpetscSF)                         :: sf_node
        INTEGER                                :: nroots, n_ghosts
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mat_loc_to_glob
        INTEGER, DIMENSION(:,:),   ALLOCATABLE :: mat_proc_np_loc
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: zz_contig_1
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: zz_contig_2
    CONTAINS
        PROCEDURE :: init_mat_loc_to_glob
        PROCEDURE :: fill_mat
    END TYPE petsc_csr_LA


CONTAINS

    SUBROUTINE init_mat_loc_to_glob(this, mesh)
        USE my_util, ONLY: error_petsc, local_error_petsc, to_str
        IMPLICIT NONE
        CLASS(petsc_csr_LA), INTENT(INOUT) :: this
        TYPE(mesh_type),              INTENT(IN)    :: mesh
        INTEGER :: m, nw, ni, i, im1, nj, j, jm1, ja_p_loc, p, proc_ja_p, n
        INTEGER :: shift_ini, idx, nleaves, ierr
        LOGICAL :: if_test, local_if_found_all_pairs, local_saw_all_stencils, saw_all_stencils, if_found_all_pairs
        
        INTEGER, DIMENSION(:),           ALLOCATABLE :: other_j_loc, idx_loc_stencil
#if (PETSC_VERSION_MINOR < 23)
        TYPE(PetscSFNode), DIMENSION(:), ALLOCATABLE :: roots
#else
        TYPE(sPetscSFNode), DIMENSION(:), ALLOCATABLE :: roots
#endif

        !=== Consistency checks
        IF (.NOT. ASSOCIATED(this%loc_to_glob)) THEN
            CALL error_petsc('BUG in init_mat_loc_to_glob: start by calling st_csr_aij_glob_block_with_extra_layer')
        ELSE IF (SIZE(this%loc_to_glob, 1) > 1) THEN
            CALL error_petsc('BUG in init_mat_loc_to_glob: SIZE(LA%loc_to_glob, 1) > 1 not programmed yet - you tried with size = '&
            &//to_str(SIZE(this%loc_to_glob, 1)))
        END IF
        !=== Consistency checks

        !=== Necessary and useful allocations
        nw = mesh%gauss%n_w
        this%n_ghosts = 0
        ALLOCATE(this%mat_loc_to_glob(nw, nw, mesh%me), SOURCE=-1)
        !=== Necessary and useful allocations

        !=========================!
        !===== STEP 1 ============!  fill every entry of mat_loc_to_glob that 
        !=========================!  can be filled easily with only the information from that proc
        DO m=1, mesh%me
            DO ni=1, nw
                i = mesh%jj(ni, m)
                im1 = i - 1
                DO nj=1, nw
                    j = mesh%jj(nj, m)
                    jm1 = j - 1
                    if_test = .FALSE.
                    !=== stencil owned locally ===!
                    IF (i <= mesh%dom_np) THEN
                        DO p = this%ia(im1), this%ia(im1+1)-1
                            proc_ja_p = mesh%get_proc(this%ja(p)+1, 'np') !=== WARNING: LA%ja is global and p is local
                            IF (proc_ja_p /= mesh%proc) THEN
                                DO n=mesh%dom_np+1, mesh%np
                                    IF (mesh%loc_to_glob(n) == this%ja(p)+1) EXIT
                                END DO
                                ja_p_loc = n - 1
                            ELSE
                                ja_p_loc = this%ja(p) - (mesh%disp(mesh%proc)-1)
                            END IF
                            IF (ja_p_loc==jm1) THEN
                                this%mat_loc_to_glob(ni,nj,m) = p !=== WARNING ROW ORIENTATION
                                if_test = .TRUE.
                                EXIT
                            END IF
                        END DO
                        IF (.NOT.if_test) THEN
                            CALL local_error_petsc("BUG in Init mat_loc_to_glob, cannot find the value for i,j: "//to_str(i)//','//to_str(j)//" in the CSR structure")
                        END IF
                    !=== stencil owned by another proc ===!
                    ELSE
                        this%n_ghosts = this%n_ghosts + 1
                    END IF
                END DO
            END DO
        END DO

        !=========================!
        !===== END STEP 1 ========!
        !=========================!

        !=========================!
        !===== STEP 2 ============!  Now fill the parts of mat_loc_to_glob that require communicating
        !=========================!

        ALLOCATE(this%mat_proc_np_loc(2, this%n_ghosts))
        ALLOCATE(other_j_loc(mesh%np), SOURCE=0)

        if_found_all_pairs = .FALSE.
        saw_all_stencils   = .FALSE.
        
        ALLOCATE(idx_loc_stencil(mesh%np))
        idx_loc_stencil(1:mesh%dom_np) = this%ia(0:mesh%dom_np-1)

        !=== Go through all stencils ===!
        !=== Fails because mat_loc_to_glob only contains (i,j) \in (dom_np, np) and not \in (dom_np, n_extra)
        ! shift_ini = MAXVAL(mat_loc_to_glob)
        !=== Fails because mat_loc_to_glob only contains (i,j) \in (dom_np, np) and not \in (dom_np, n_extra)

        shift_ini = this%ia(mesh%dom_np) - 1
        idx = 0

        !=== Loop that will check all stencils of all nodes, and communicate the relevant ones ===!
        DO WHILE (.NOT. (if_found_all_pairs .OR. saw_all_stencils)) !=== Loop that will check all stencils of all nodes, and communicate the relevant ones
            local_saw_all_stencils = .TRUE.
            !===============================================================!
            !=== SUBSTEP 2.1: move to next node in line for all stencils ===!
            !===============================================================!
            DO i=1, mesh%dom_np
                IF (idx_loc_stencil(i) < this%ia(i)) THEN
                    p = idx_loc_stencil(i)
                    other_j_loc(i) = this%ja(p) !=== WARNING: THIS IS TRICKY BECAUSE LA%JA IS ACTUALLY A GLOBAL NUMBERING
                    local_saw_all_stencils = .FALSE.
                    idx_loc_stencil(i) = idx_loc_stencil(i) + 1 !=== SEE WARNING BELOW
                END IF
            END DO

            !====================================================================================================!
            !=== SUBSTEP 2.2: communicate node i owned by proc n to proc np if i \in [dom_np+1,np] of proc np ===!
            !====================================================================================================!
            CALL mesh%bulk_to_ghost_int(other_j_loc, MPI_REPLACE)
            CALL mesh%bulk_to_ghost_int(idx_loc_stencil, MPI_REPLACE) !=== BEWARE FOR idx_loc_stencil IS INCREMENTED ONE TIME TOO MUCH AT THAT POINT

            !====================================================================================================!
            !=== SUBSTEP 2.3: checking if node i owned by proc n is in a stencil owned by proc np ===============!
            !====================================================================================================!
            DO m=1, mesh%me
                DO ni=1, nw
                    i = mesh%jj(ni, m)
                    IF (i <= mesh%dom_np) CYCLE
                    im1 = i - 1
                    DO nj=1, nw
                        j = mesh%jj(nj, m)
                        !=== stencil owned locally ===!
                        IF (other_j_loc(i) == mesh%loc_to_glob(j)-1 .AND. this%mat_loc_to_glob(ni,nj,m) == -1) THEN
                            idx = idx + 1
                            this%mat_loc_to_glob(ni,nj,m) = shift_ini + idx
                            this%mat_proc_np_loc(1, idx)  = mesh%proc_np_loc(1, i-mesh%dom_np)
                            this%mat_proc_np_loc(2, idx)  = idx_loc_stencil(i) !=== ends up in Fortran convention (starting from 1, not 0)
                        END IF
                        !=== stencil owned by another proc ===!
                    END DO
                END DO
            END DO

            !===============================================================================================================!
            !=== SUBSTEP 2.4: checking if there are still stencils owned by np but computed by n which are still unknown ===!
            !===============================================================================================================!
            local_if_found_all_pairs = (.NOT. ANY(this%mat_loc_to_glob == -1))

            CALL MPI_ALLREDUCE(local_saw_all_stencils,   saw_all_stencils,   1, MPI_LOGICAL, MPI_LAND, mesh%comm, ierr)
            CALL MPI_ALLREDUCE(local_if_found_all_pairs, if_found_all_pairs, 1, MPI_LOGICAL, MPI_LAND, mesh%comm, ierr)

        END DO

        !=== one last consistency check
        IF (saw_all_stencils .AND. ANY(this%mat_loc_to_glob == -1)) THEN
            CALL local_error_petsc('BUG when looking through stencils => saw all stencils, and yet could not find all matrices pairs.')
        END IF
        !=== one last consistency check
        !=========================!
        !===== END STEP 2 ========!  Now mat_loc_to_glob is complete, ghosts included
        !=========================!  mat_proc_np_loc contains respectively the rank and local index of all ghosted elements of mat_loc_to_glob

        
        !=========================!
        !===== STEP 3 ============!  Create petscSF for communicating ghosts
        !=========================!

        CALL PetscSFCreate(mesh%comm, this%sf_node, ierr)
        CALL PetscSFSetFromOptions(this%sf_node, ierr)

        nleaves = this%n_ghosts
        this%nroots  = this%ia(mesh%dom_np)
        ALLOCATE(roots(nleaves))
        DO n=1, nleaves
            !=== WARNING: petsc conventions indexing starts from 0, not 1
            roots(n)%rank  = this%mat_proc_np_loc(1, n) - 1
            roots(n)%index = this%mat_proc_np_loc(2, n) - 1
            !=== WARNING: petsc conventions indexing starts from 0, not 1
        END DO
#if (PETSC_VERSION_MINOR <= 21)
    CALL PetscSFSetGraph(this%sf_node, this%nroots, nleaves, PETSC_NULL_INTEGER, PETSC_COPY_VALUES, roots, PETSC_COPY_VALUES, ierr)
#else
    CALL PetscSFSetGraph(this%sf_node, this%nroots, nleaves, PETSC_NULL_INTEGER_ARRAY, PETSC_COPY_VALUES, roots, PETSC_COPY_VALUES, ierr)
#endif

        !=========================!
        !===== END STEP 3 ========!  Now everything is known to use MatUpdateMPIAIJ
        !=========================!

        ALLOCATE(this%zz_contig_1(0 : (this%ia(mesh%dom_np)-1) + this%n_ghosts))
        ALLOCATE(this%zz_contig_2(0 : (this%ia(mesh%dom_np)-1) + this%n_ghosts))
        DEALLOCATE(roots, other_j_loc, idx_loc_stencil)
    END SUBROUTINE

    SUBROUTINE fill_mat(this, matrix, zz_contig)
        !=== WARNING 02/07/2026 MATRIX MUST NOT HAVE BEEN GENERATED USING MATDUPLICATE
        USE petscmat
        USE petscvec
        USE my_util
        IMPLICIT NONE
        CLASS(petsc_csr_LA) :: this
        TYPE(tMat)                   :: matrix
        INTEGER                      :: ierr
        REAL(KIND=8), DIMENSION(0:SIZE(this%zz_contig_1)-1), INTENT(INOUT) :: zz_contig

        CALL PetscSFReduceBegin(this%sf_node, MPI_DOUBLE_PRECISION, zz_contig(this%nroots: ), zz_contig(0 : this%nroots - 1), MPI_SUM, ierr)
        CALL PetscSFReduceEnd  (this%sf_node, MPI_DOUBLE_PRECISION, zz_contig(this%nroots: ), zz_contig(0 : this%nroots - 1), MPI_SUM, ierr)
        CALL MatUpdateMPIAIJWithArray(matrix, zz_contig(:this%nroots - 1), ierr)

    END SUBROUTINE fill_mat

END MODULE petsc_csr_LA_module