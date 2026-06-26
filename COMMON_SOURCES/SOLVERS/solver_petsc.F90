MODULE solver_petsc
    USE my_util
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    USE solver_data_module
   
CONTAINS

    SUBROUTINE init_solver(communicator, solver_paramater, my_ksp, matrix, opt_mat_pre, opt_re_init)          
        USE character_strings
        USE petsc
        USE petscmat
        USE petscsys
        USE my_util
        IMPLICIT NONE
        LOGICAL, INTENT(IN), OPTIONAL :: opt_re_init
        TYPE(solver_data_type) :: solver_paramater
        CHARACTER(LEN=:), ALLOCATABLE :: solver, precond
        LOGICAL :: re_init
        INTEGER :: deb, fin

        Mat            :: matrix, matrix_shell
        Mat, OPTIONAL :: opt_mat_pre
        KSP            :: my_ksp
        PC             :: prec
        PetscErrorCode :: ierr
        MPI_Comm       :: communicator
        PetscViewerAndFormat :: vf

        solver = TRIM(ADJUSTL(solver_paramater%solver))
        precond = TRIM(ADJUSTL(solver_paramater%precond))

        IF (.NOT.PRESENT(opt_re_init)) THEN
            re_init = .FALSE.
        ELSE
            re_init = opt_re_init
        END IF

        IF (solver_paramater%it_max.LE.0) THEN
            !   solver_paramater%it_max = 100
        END IF
        IF (solver_paramater%rel_tol.LE.0.d0) THEN
            solver_paramater%rel_tol = 1.d-8
        END IF
        IF (solver_paramater%abs_tol.LE.0.d0) THEN
            solver_paramater%abs_tol = 1.d-14
        END IF
        IF (.NOT.re_init) CALL KSPCreate(communicator, my_ksp, ierr)

        IF (PRESENT(opt_mat_pre)) THEN
            CALL KSPSetOperators(my_ksp, matrix, opt_mat_pre, ierr) !Petsc 3.7.2
        ELSE
            CALL KSPSetOperators(my_ksp, matrix, matrix, ierr) !Petsc 3.7.2
        END IF

        IF (solver_paramater%if_residual) CALL KSPMonitorSet(my_ksp, MyKSPMonitor, vf, PetscViewerAndFormatDestroy, ierr)

        SELECT CASE(solver)
        CASE('BCGS')
            CALL KSPSetType(my_ksp, KSPBCGS, ierr)
        CASE('GMRES')
            CALL KSPSetType(my_ksp, KSPGMRES, ierr)
        CASE('FGMRES')
            CALL KSPSetType(my_ksp, KSPFGMRES, ierr)
        CASE('PCR')
            CALL KSPSetType(my_ksp, KSPCR, ierr)
        CASE('CHEBYCHEV')
            CALL KSPSetType(my_ksp, KSPCHEBYSHEV, ierr)
        CASE('RICHARDSON')
            CALL KSPSetType(my_ksp, KSPRICHARDSON, ierr)
            CALL KSPSetNormType(my_ksp, KSP_NORM_NONE, ierr)
        CASE('CG')
            CALL KSPSetType(my_ksp, KSPCG, ierr)
        CASE DEFAULT
            CALL KSPSetType(my_ksp, KSPFGMRES, ierr)
        END SELECT

        CALL KSPGetPC(my_ksp, prec, ierr)
        CALL KSPSetTolerances(my_ksp, solver_paramater%rel_tol, solver_paramater%abs_tol, &
            PETSC_DEFAULT_REAL, solver_paramater%it_max, ierr)

        !  CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-ksp_initial_guess_nonzero', 'true', ierr)
        SELECT CASE(precond)
        CASE('JACOBI')
            !CALL PCSetType(prec, PCJACOBI, ierr)
            CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_type', 'jacobi', ierr)
        CASE('HYPRE')
            IF (solver_paramater%if_fixed_v_cycle) THEN
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_type', 'hypre', ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_type', 'boomeramg', ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_strong_threshold', &
                    solver_paramater%boomeramg_strong_threshold, ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_coarsen_type', &
                    solver_paramater%boomeramg_coarsen_type, ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_relax_type_all', &
                    solver_paramater%boomeramg_relax_type_all, ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_max_iter', &
                    solver_paramater%number_v_cycle, ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_tol', &
                    '0.0', ierr)
            ELSE
               !my_par%rel_tol = 1.d-4
               !CALL PCSetType(prec, PCHYPRE, ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_type', 'hypre', ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_type', 'boomeramg', ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_nodal_coarsen', '1', ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_coarsen_type', '0', ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_strong_threshold', &
                    solver_paramater%boomeramg_strong_threshold, ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_coarsen_type', &
                    solver_paramater%boomeramg_coarsen_type, ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_agg_nl', '2', ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_agg_num_paths', '4', ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_truncfactor','.05', ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_interp_type', 'multipass', ierr)
               CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_relax_type_all', &
                    solver_paramater%boomeramg_relax_type_all, ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_relax_type_down', 'Chebyshev', ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_relax_type_up', 'Chebyshev', ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_relax_type_coarse', 'Gaussian-elimination', ierr)
               !!$HYPRE preconditioner options
               !!$  -pc_hypre_type <boomeramg> (choose one of) pilut parasails boomeramg ams (PCHYPRESetType)
               !!$HYPRE BoomerAMG Options
               !!$  -pc_hypre_boomeramg_cycle_type <V> (choose one of) V W (None)
               !!$  -pc_hypre_boomeramg_max_levels <25>: Number of levels (of grids) allowed (None)
               !!$  -pc_hypre_boomeramg_max_iter <1>: Maximum iterations used PER hypre call (None)
               !!$  -pc_hypre_boomeramg_tol <0.>: Convergence tolerance PER hypre call (0.0 = use a fixed number of iterations) (None)
               !!$  -pc_hypre_boomeramg_truncfactor <0.>: Truncation factor for interpolation (0=no truncation) (None)
               !!$  -pc_hypre_boomeramg_P_max <0>: Max elements per row for interpolation operator (0=unlimited) (None)
               !!$  -pc_hypre_boomeramg_agg_nl <0>: Number of levels of aggressive coarsening (None)
               !!$  -pc_hypre_boomeramg_agg_num_paths <1>: Number of paths for aggressive coarsening (None)
               !!$  -pc_hypre_boomeramg_strong_threshold <0.25>: Threshold for being strongly connected (None)
               !!$  -pc_hypre_boomeramg_max_row_sum <0.9>: Maximum row sum (None)
               !!$  -pc_hypre_boomeramg_grid_sweeps_all <1>: Number of sweeps for the up and down grid levels (None)
               !!$  -pc_hypre_boomeramg_nodal_coarsen <0>: Use a nodal based coarsening 1-6 (HYPRE_BoomerAMGSetNodal)
               !!$  -pc_hypre_boomeramg_vec_interp_variant <0>: Variant of algorithm 1-3 (HYPRE_BoomerAMGSetInterpVecVariant)
               !!$  -pc_hypre_boomeramg_grid_sweeps_down <1>: Number of sweeps for the down cycles (None)
               !!$  -pc_hypre_boomeramg_grid_sweeps_up <1>: Number of sweeps for the up cycles (None)
               !!$  -pc_hypre_boomeramg_grid_sweeps_coarse <1>: Number of sweeps for the coarse level (None)
               !!$  -pc_hypre_boomeramg_smooth_type <Schwarz-smoothers> (choose one of) Schwarz-smoothers Pilut ParaSails Euclid (None)
               !!$  -pc_hypre_boomeramg_smooth_num_levels <25>: Number of levels on which more complex smoothers are used (None)
               !!$  -pc_hypre_boomeramg_eu_level <0>: Number of levels for ILU(k) in Euclid smoother (None)
               !!$  -pc_hypre_boomeramg_eu_droptolerance <0.>: Drop tolerance for ILU(k) in Euclid smoother (None)
               !!$  -pc_hypre_boomeramg_eu_bj: <FALSE> Use Block Jacobi for ILU in Euclid smoother? (None)
               !!$  -pc_hypre_boomeramg_relax_type_all <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
               !!$  -pc_hypre_boomeramg_relax_type_down <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
               !!$  -pc_hypre_boomeramg_relax_type_up <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
               !!$  -pc_hypre_boomeramg_relax_type_coarse <Gaussian-elimination> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination      CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
               !!$  -pc_hypre_boomeramg_relax_weight_all <1.>: Relaxation weight for all levels (0 = hypre estimates, -k = determined with k CG steps) (None)
               !!$  -pc_hypre_boomeramg_relax_weight_level <1.>: Set the relaxation weight for a particular level (weight,level) (None)
               !!$  -pc_hypre_boomeramg_outer_relax_weight_all <1.>: Outer relaxation weight for all levels (-k = determined with k CG steps) (None)
               !!$  -pc_hypre_boomeramg_outer_relax_weight_level <1.>: Set the outer relaxation weight for a particular level (weight,level) (None)
               !!$  -pc_hypre_boomeramg_no_CF: <FALSE> Do not use CF-relaxation (None)
               !!$  -pc_hypre_boomeramg_measure_type <local> (choose one of) local global (None)
               !!$  -pc_hypre_boomeramg_coarsen_type <Falgout> (choose one of) CLJP Ruge-Stueben  modifiedRuge-Stueben   Falgout  PMIS  HMIS (None)
               !!$  -pc_hypre_boomeramg_interp_type <classical> (choose one of) classical   direct multipass multipass-wts ext+i ext+i-cc standard standard-wts   FF FF1 (None)
               !!$  -pc_hypre_boomeramg_print_statistics: Print statistics (None)
               !!$  -pc_hypre_boomeramg_print_statistics <3>: Print statistics (None)
               !!$  -pc_hypre_boomeramg_print_debug: Print debug information (None)
               !!$  -pc_hypre_boomeramg_nodal_relaxation: <FALSE> Nodal relaxation via Schwarz (None)

               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_boomeramg_vec_interp_variant', '1', ierr)
               !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_factor_levels','100',ierr)
               !CALL PCSetType(prec, PCGAMG, ierr)
            END IF
        CASE('SSOR')
            !CALL PCSetType(prec, PCSOR, ierr)
            CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_type', 'none', ierr)
        CASE('MUMPS')
            !CALL PCSetType(prec, PCLU, ierr)
            CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_type', 'lu', ierr)
            !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mat_mumps_icntl_35', '1', ierr)
            !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mat_mumps_icntl_36', '1', ierr)
            !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mat_mumps_cntl_7', '0.000001', ierr)
            !CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-mat_mumps_cntl_1', '0.0', ierr)
            IF (.NOT. PRESENT(opt_mat_pre)) THEN
               CALL KSPSetType(my_ksp, KSPPREONLY, ierr)
            END IF
            !CALL PCFactorSetMatSolverType(prec, MATSOLVERMUMPS, ierr) !
        CASE DEFAULT
            !CALL PCSetType(prec, PCHYPRE, ierr)
            CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_type', 'hypre', ierr)
            CALL PetscOptionsSetValue(PETSC_NULL_OPTIONS, '-pc_hypre_type', 'boomeramg', ierr)
        END SELECT

        CALL KSPSetFromOptions(my_ksp, ierr)
    END SUBROUTINE init_solver

    SUBROUTINE solver(my_ksp, b, x, reinit, verbose, type)
        use petsc
        IMPLICIT NONE
        LOGICAL, OPTIONAL :: reinit, verbose
        CHARACTER(*), OPTIONAL, INTENT(IN) :: type
        INTEGER :: its, rank
        KSP            :: my_ksp
        PetscErrorCode :: ierr
        Vec            :: x, b
        KSPConvergedReason :: reason
        IF (.NOT.PRESENT(reinit)) reinit = .TRUE.
        IF (.NOT.PRESENT(verbose)) verbose = .FALSE.

        IF (reinit) CALL VecZeroEntries (x, ierr)
        CALL KSPSolve(my_ksp, b, x, ierr)
        IF(PRESENT(type)) THEN
            CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
            CALL KSPGetIterationNumber(my_ksp, its, ierr)
            IF (rank == 0)  WRITE(*, *) "Nb of GMRES iterations for ", type, its
        END IF
        IF (verbose) THEN
            CALL KSPGetIterationNumber(my_ksp, its, ierr)
            CALL KSPGetConvergedReason(my_ksp, reason, ierr)
            SELECT CASE(reason)
            CASE(2)
                WRITE(*, *) "KSP_CONVERGED_RTOL, Nb of iterations", its
            CASE(3)
                WRITE(*, *) "KSP_CONVERGED_ATOL, Nb of iterations", its
            CASE(4)
                WRITE(*, *) "Converged after one single iteration of the preconditioner is applied"
            CASE(5, 6, 7, 8)
                WRITE(*, *) "Converge for strange reason:", reason
            CASE(-2)
                WRITE(*, *) "KSP_DIVERGED_NULL"
            CASE(-3)
                WRITE(*, *) "Not converged after it_max", its
            CASE(-4)
                WRITE(*, *) "Not converged: explosion"
            CASE(-5, -6, -7)
                WRITE(*, *) "Not converged for strange reasons", reason
            CASE(-8)
                WRITE(*, *) "Not converged: Indefinite preconditioner"
            CASE(-9)
                WRITE(*, *) "Not converged: NAN"
            CASE(-10)
                WRITE(*, *) "Not converged: Indefinite matrix"
            CASE DEFAULT
                WRITE(*, *) "Something strange happened", reason
            END SELECT
        END IF

    END SUBROUTINE solver

   SUBROUTINE create_local_petsc_matrix(communicator, LA, matrix, clean)
      USE def_type_mesh
      use petsc
      IMPLICIT NONE
      TYPE(petsc_csr_LA) :: LA
      LOGICAL, OPTIONAL :: clean
      REAL(KIND = 8), DIMENSION(:), POINTER :: aa
      INTEGER :: nnzm1, dom_np
      LOGICAL :: test_clean
      !!$  INTEGER, DIMENSION(:), POINTER :: ia, ja
      MPI_Comm       :: communicator
      Mat            :: matrix
      PetscErrorCode :: ierr
      !------------------------------------------------------------------------------
      dom_np = SIZE(LA%ia) - 1
      nnzm1 = LA%ia(dom_np) - LA%ia(0) - 1
      ALLOCATE(aa(0:nnzm1))
      aa = 0.d0


      !!$  ALLOCATE(ia(0:dom_np),ja(0:nnzm1))
      !!$  ia = LA%ia
      !!$  ja = LA%ja
      !!$  CALL MatCreateMPIAIJWithArrays(communicator,dom_np,dom_np,PETSC_DECIDE, &
      !!$       PETSC_DECIDE, ia, ja, aa, matrix, ierr)
      !!$DEALLOCATE(ia,ja)

      CALL MatCreateMPIAIJWithArrays(communicator, dom_np, dom_np, PETSC_DECIDE, &
           PETSC_DECIDE, LA%ia, LA%ja, aa, matrix, ierr)

      DEALLOCATE(aa)
      IF (PRESENT(clean)) THEN
         test_clean = clean
      ELSE
         test_clean = .TRUE.
      END IF
      IF (test_clean) THEN
         IF (ASSOCIATED(LA%ia)) DEALLOCATE(LA%ia)
         IF (ASSOCIATED(LA%ja)) DEALLOCATE(LA%ja)
      END IF
   END SUBROUTINE create_local_petsc_matrix

   SUBROUTINE create_local_petsc_matrix_a_detruire(communicator, aij, i_loc, matrix)
      USE def_type_mesh
      use petsc
      IMPLICIT NONE
      TYPE(aij_type), INTENT(IN) :: aij
      INTEGER, DIMENSION(2) :: i_loc
      INTEGER, DIMENSION(:), POINTER :: ia, ja
      REAL(KIND = 8), DIMENSION(:), POINTER :: aa
      INTEGER :: nnzm1, dom_np, p, i, n

      MPI_Comm       :: communicator
      Mat            :: matrix
      PetscErrorCode :: ierr
      !------------------------------------------------------------------------------
      dom_np = i_loc(2) - i_loc(1) + 1
      nnzm1 = aij%ia(i_loc(2) + 1) - aij%ia(i_loc(1)) - 1
      ALLOCATE(ia(0:dom_np), ja(0:nnzm1))
      ia(0) = 0
      DO i = 1, dom_np
         n = i_loc(1) + i - 1
         ia(i) = aij%ia(n + 1) - aij%ia(i_loc(1))
         DO p = aij%ia(n), aij%ia(n + 1) - 1
            ja(p - aij%ia(i_loc(1))) = aij%ja(p) - 1
         END DO
      END DO
      !------------------------------------------------------------------------------
      ALLOCATE(aa(0:nnzm1))
      aa = 0
      CALL MatCreateMPIAIJWithArrays(communicator, dom_np, dom_np, PETSC_DECIDE, &
           PETSC_DECIDE, ia, ja, aa, matrix, ierr)

      DEALLOCATE(ia, ja, aa)
   END SUBROUTINE create_local_petsc_matrix_a_detruire

   SUBROUTINE create_local_petsc_block_matrix(communicator, n_b, aij, i_loc, matrix)
      USE def_type_mesh
      use petsc
      IMPLICIT NONE
      TYPE(aij_type), INTENT(IN) :: aij
      INTEGER, DIMENSION(2) :: i_loc
      INTEGER :: n_b
      INTEGER, DIMENSION(:), POINTER :: ia, ja
      REAL(KIND = 8), DIMENSION(:), POINTER :: aa
      INTEGER :: nnzm1, dom_np, p, i, n, ib, k

      MPI_Comm       :: communicator
      Mat            :: matrix
      PetscErrorCode :: ierr

      dom_np = i_loc(2) - i_loc(1) + 1
      nnzm1 = n_b * (aij%ia(i_loc(2) + 1) - aij%ia(i_loc(1)) - 1)
      ALLOCATE(ia(0:n_b * dom_np), ja(0:nnzm1))
      ia(0) = 0
      DO k = 1, n_b
         DO i = 1, dom_np
            ib = i + (k - 1) * dom_np
            n = i_loc(1) + i - 1
            ia(i) = n_b * (aij%ia(n + 1) - aij%ia(i_loc(1)))
            DO p = aij%ia(n), aij%ia(n + 1) - 1
               ja(p - aij%ia(i_loc(1))) = aij%ja(p) - 1
            END DO
         END DO
      END DO
      !------------------------------------------------------------------------------
      ALLOCATE(aa(0:nnzm1))
      aa = 0
      CALL MatCreateMPIAIJWithArrays(communicator, dom_np, dom_np, PETSC_DECIDE, &
           PETSC_DECIDE, ia, ja, aa, matrix, ierr)

      DEALLOCATE(ia, ja, aa)
   END SUBROUTINE create_local_petsc_block_matrix

   subroutine MyKSPMonitor(ksp, n, rnorm, dummy, ierr)
      use petscksp

      KSP         ::     ksp
      Vec        ::      x
      PetscErrorCode :: ierr
      PetscInt :: n, dummy
      PetscMPIInt :: rank
      PetscReal :: rnorm
      CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
      IF (rank == 0)  WRITE(*, *) 'Residual', rnorm

   end subroutine MyKSPMonitor


END MODULE solver_petsc
