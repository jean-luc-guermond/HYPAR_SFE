MODULE restart_module

    TYPE read_write_type
        INTEGER :: counter
    CONTAINS
        PROCEDURE :: write_restart
        PROCEDURE :: read_restart    
    END TYPE read_write_type

CONTAINS

    SUBROUTINE write_restart(this, mesh, time, un, filename_in, opt_if_series)
        USE def_type_mesh
        USE mesh_parameters
        IMPLICIT NONE
        CLASS(read_write_type)                   :: this
        TYPE(mesh_type)                          :: mesh
        REAL(KIND=8),                 INTENT(IN) :: time
        REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: un
        CHARACTER(len=*),             INTENT(IN) :: filename_in !'EULER/NS'
        LOGICAL,  OPTIONAL,           INTENT(IN) :: opt_if_series
        CHARACTER(LEN=200)                       :: filename, char_S
        CHARACTER(LEN=4)                         :: char_I
        INTEGER :: n_blank, n
        
        WRITE(char_S, '(I0)') mesh%rank

        IF (PRESENT(opt_if_series)) THEN
            IF (opt_if_series) THEN
                this%counter = this%counter + 1
                WRITE(char_I, '(I4)') this%counter
                n_blank = len(char_I) - len(TRIM(ADJUSTL(char_I)))
                char_I(1:n_blank) = REPEAT('0',n_blank)
                filename = 'series_'//TRIM(ADJUSTL(filename_in))//'_I'//TRIM(ADJUSTL(char_I))//'_S'//TRIM(ADJUSTL(char_S))//'.'//mesh_data_info%file_name
            ELSE
                filename = 'backup_'//TRIM(ADJUSTL(filename_in))//'_S'//TRIM(ADJUSTL(char_S))//'.'//mesh_data_info%file_name
            END IF
        ELSE
            filename = 'backup_'//TRIM(ADJUSTL(filename_in))//'_S'//TRIM(ADJUSTL(char_S))//'.'//mesh_data_info%file_name
        END IF

        OPEN(UNIT = 10, FILE = filename, FORM = 'unformatted', STATUS = 'unknown')
        WRITE(10) time, mesh%np, SIZE(un,2)
        WRITE(10) un
        CLOSE(10)
    END SUBROUTINE write_restart

    SUBROUTINE read_restart(this, mesh, time, un, filename_in, opt_it)
        USE def_type_mesh
        USE my_util
        USE mesh_parameters
        IMPLICIT NONE
        CLASS(read_write_type)                   :: this
        TYPE(mesh_type)                          :: mesh
        REAL(KIND=8),                 INTENT(OUT):: time
        REAL(KIND=8), DIMENSION(:,:), INTENT(OUT):: un
        CHARACTER(len=*),             INTENT(IN) :: filename_in !'EULER/NS'
        INTEGER,  OPTIONAL,           INTENT(IN) :: opt_it 
        INTEGER :: io_error, np, syst_dim
        CHARACTER(LEN=200)                       :: filename, char_S
        CHARACTER(LEN=4)                         :: char_I
        INTEGER :: n_blank

        WRITE(char_S, '(I0)') mesh%rank

        IF (PRESENT(opt_it)) THEN
            IF (opt_it >= 0) THEN
                this%counter = opt_it
                WRITE(char_I, '(I4)') this%counter
                n_blank = len(char_I) - len(TRIM(ADJUSTL(char_I)))
                char_I(1:n_blank) = REPEAT('0',n_blank)

                filename = 'series_'//TRIM(ADJUSTL(filename_in))//'_I'//TRIM(ADJUSTL(char_I))//'_S'//TRIM(ADJUSTL(char_S))//'.'//mesh_data_info%file_name
            ELSE
                filename = 'backup_'//TRIM(ADJUSTL(filename_in))//'_S'//TRIM(ADJUSTL(char_S))//'.'//mesh_data_info%file_name
                this%counter = 0
            END IF
        ELSE
            filename = 'backup_'//TRIM(ADJUSTL(filename_in))//'_S'//TRIM(ADJUSTL(char_S))//'.'//mesh_data_info%file_name
            this%counter = 0
        END IF

        OPEN(UNIT = 10, FILE = filename, FORM = 'unformatted', STATUS = 'old', iostat=io_error)
        IF (io_error > 0) THEN
            CALL error_petsc('BUG in read_restart: could not find file '//filename)
        END IF
        READ(10) time, np, syst_dim
        IF (np /= mesh%np) THEN
            CALL error_petsc('BUG in read_restart: wrong mesh%np '// to_str(np)//'/='//to_str(mesh%np))
        ELSE IF (syst_dim /= SIZE(un, 2)) THEN
            CALL error_petsc('BUG in read_restart: wrong syst_dim '// to_str(syst_dim)//'/='//to_str(SIZE(un, 2)))
        END IF
        IF (mesh%rank==0) WRITE(*,*) 'Reading series '//TRIM(ADJUSTL(filename_in))//' it= '//to_str(this%counter)
        READ(10) un
        CLOSE(10)
    END SUBROUTINE read_restart

END MODULE restart_module