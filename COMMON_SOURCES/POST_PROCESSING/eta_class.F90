MODULE ETA_module
    USE my_util

    TYPE ETA_type
        REAL(KIND=8) :: dt_avg, tps_eta, ETA, init_time
        INTEGER :: counter, rank
    CONTAINS
        PROCEDURE :: init
        PROCEDURE :: print
        PROCEDURE :: update
    END TYPE ETA_type

CONTAINS

    SUBROUTINE init(this, rank, init_time)
        IMPLICIT NONE
        CLASS(ETA_type) :: this
        INTEGER      :: rank
        REAL(KIND=8) :: init_time
        this%counter = 0
        this%dt_avg  = 0.d0
        this%tps_eta = USER_time()
        this%rank = rank
        this%init_time = init_time
    END SUBROUTINE init

    SUBROUTINE update(this, dt)
        IMPLICIT NONE
        CLASS(ETA_type) :: this
        REAL(KIND=8)    :: dt
        this%dt_avg = this%dt_avg + dt
        this%counter = this%counter + 1
    END SUBROUTINE update

    SUBROUTINE print(this, cur_time, final_time)
        IMPLICIT NONE
        CLASS(ETA_type) :: this
        REAL(KIND=8), INTENT(IN) :: cur_time, final_time
        !== compute averages
        this%tps_eta = USER_time() - this%tps_eta
        this%tps_eta = this%tps_eta/this%counter
        this%dt_avg = this%dt_avg/this%counter
        this%ETA = (final_time-cur_time)*this%tps_eta/this%dt_avg
        !=== Print averages
        IF (this%rank==0) THEN
            IF (this%ETA<60.d0) THEN
                write(*, '(A,f4.1,A,f5.1,A,e9.3,A,e9.3,A,e9.3)') 'ETA: ', this%ETA, ' secs; execution ratio: ', &
                (cur_time-this%init_time)/(final_time-this%init_time)*100, '%, timestep: ', this%dt_avg, &
                ', t_sim: ', cur_time, ', t_fin: ', final_time
            ELSE IF (this%ETA<3600.d0) THEN
                write(*, '(A,f4.1,A,f5.1,A,e9.3,A,e9.3,A,e9.3)') 'ETA: ', this%ETA/60, ' mins; execution ratio: ', &
                (cur_time-this%init_time)/(final_time-this%init_time)*100, '%, timestep: ', this%dt_avg, &
                ', t_sim: ', cur_time, ', t_fin: ', final_time
            ELSE IF (this%ETA<3600.d0*24.d0) THEN
                write(*, '(A,f4.1,A,f5.1,A,e9.3,A,e9.3,A,e9.3)') 'ETA: ', this%ETA/3600, ' hours; execution ratio: ', &
                (cur_time-this%init_time)/(final_time-this%init_time)*100, '%, timestep: ', this%dt_avg, &
                ', t_sim: ', cur_time, ', t_fin: ', final_time
            ELSE
                write(*, '(A,f4.1,A,f5.1,A,e9.3,A,e9.3,A,e9.3)') 'ETA: ', this%ETA/(3600*24), ' days; execution ratio: ', &
                (cur_time-this%init_time)/(final_time-this%init_time)*100, '%, timestep: ', this%dt_avg, &
                ', t_sim: ', cur_time, ', t_fin: ', final_time
            END IF
        END IF
        !=== Reinit for next it
        CALL this%init(this%rank, this%init_time)
    END SUBROUTINE print

END MODULE ETA_module