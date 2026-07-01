MODULE profiler_module

    TYPE :: profiler_type
        REAL(KIND=8),        DIMENSION(:), ALLOCATABLE :: tot_time
        REAL(KIND=8),        DIMENSION(:), ALLOCATABLE :: cur_time
        INTEGER,             DIMENSION(:), ALLOCATABLE :: counter
        CHARACTER(LEN=16),   DIMENSION(:), ALLOCATABLE :: name
        INTEGER :: nb_profilings
    CONTAINS
        PROCEDURE :: init   => init_profiling
        PROCEDURE :: start  => start_profiling
        PROCEDURE :: end    => end_profiling
        PROCEDURE :: output => output_profiling
    END TYPE profiler_type

CONTAINS

    SUBROUTINE init_profiling(this, name_in)
        IMPLICIT NONE
        CLASS(profiler_type), INTENT(INOUT) :: this
        CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: name_in

        this%name = name_in
        this%nb_profilings = SIZE(name_in)
        ALLOCATE(this%tot_time(this%nb_profilings), SOURCE=0.d0)
        ALLOCATE(this%cur_time(this%nb_profilings), SOURCE=0.d0)
        ALLOCATE(this%counter(this%nb_profilings), SOURCE=0)
    END SUBROUTINE init_profiling

    SUBROUTINE start_profiling(this, idx)
        USE my_util, ONLY: user_time
        IMPLICIT NONE
        CLASS(profiler_type), INTENT(INOUT) :: this
        INTEGER,         INTENT(IN)    :: idx

        this%cur_time(idx) = user_time()
    END SUBROUTINE start_profiling

    SUBROUTINE end_profiling(this, idx)
        USE my_util, ONLY: user_time
        IMPLICIT NONE
        CLASS(profiler_type), INTENT(INOUT) :: this
        INTEGER,         INTENT(IN)    :: idx

        this%tot_time(idx) = this%tot_time(idx) + (user_time()-this%cur_time(idx))
        this%cur_time(idx) = 0.d0
        this%counter(idx) = this%counter(idx) + 1
    END SUBROUTINE end_profiling

    SUBROUTINE output_profiling(this)
        USE my_util, ONLY: user_time
        USE petscmpi
        IMPLICIT NONE
        CLASS(profiler_type), INTENT(INOUT) :: this
        INTEGER                        :: idx, rank, ierr

        CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
        DO idx=1, this%nb_profilings
            WRITE(*,*) idx, 'RANK ', rank, this%name(idx), ", avg time = ", this%tot_time(idx)/this%counter(idx), &
            'percentage = ', (this%tot_time(idx)) / (this%tot_time(1)) * 100
            ! 'percentage = ', (this%tot_time(idx)/this%counter(idx)) / (this%tot_time(1)/this%counter(1)) * 100
        END DO


    END SUBROUTINE output_profiling

END MODULE profiler_module