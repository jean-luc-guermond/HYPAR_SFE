MODULE plot_vtu_module
#include "petsc/finclude/petsc.h"
  USE petsc
  LOGICAL, PARAMETER, PRIVATE :: if_xml=.true.
  CONTAINS
SUBROUTINE make_vtu_file_2D(communicator, mesh_in, header, field_in, field_name, what, opt_it)
    USE def_type_mesh
    USE my_util
    !USE plot_vtk
    IMPLICIT NONE
    TYPE(mesh_type),       INTENT(IN), TARGET :: mesh_in
    TYPE(mesh_type),               POINTER    :: mesh
    CHARACTER(*),                  INTENT(IN) :: header
    CHARACTER(*),                  INTENT(IN) :: field_name, what
    INTEGER, OPTIONAL,             INTENT(IN) :: opt_it
    REAL(KIND=8), DIMENSION(:),  INTENT(IN), TARGET :: field_in !<==Scalar fields
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE, TARGET :: field_tmp
    REAL(KIND=8), DIMENSION(:),  POINTER    :: field
    INTEGER                                   :: j, it
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: file_list
    CHARACTER(LEN=3)                          :: st_rank, st_it
    PetscErrorCode                            :: ierr
    PetscMPIInt                               :: rank, nb_procs
    MPI_Comm                                  :: communicator
    CALL MPI_Comm_rank(communicator, rank, ierr)
    CALL MPI_Comm_size(communicator, nb_procs, ierr)
    ALLOCATE(file_list(nb_procs))

    IF (mesh_in%gauss%n_w==10 .AND. mesh_in%gauss%k_d==2) THEN
       CALL error_petsc('Error in make_vtu_file_2D: P3 elements not programmed')
    ELSE
       field => field_in
       mesh => mesh_in
    END IF

    IF (PRESENT(opt_it)) THEN
       it = opt_it
       WRITE(st_it,'(I3)') it
       DO j = 1, nb_procs
          WRITE(st_rank,'(I3)') j
          file_list(j) = TRIM(header)//'_proc_'//TRIM(ADJUSTL(st_rank))//&
               '_it_'//TRIM(ADJUSTL(st_it))
       END DO
    ELSE
       DO j = 1, nb_procs
          WRITE(st_rank,'(I3)') j
          file_list(j) = TRIM(header)//'_proc_'//TRIM(ADJUSTL(st_rank))
       END DO
    END IF

    CALL check_list(communicator, file_list, mesh%np)
    IF (rank==0) THEN
       IF (PRESENT(opt_it)) THEN
          it = opt_it
       ELSE
          it = 1
       END IF
       CALL create_pvd_file(file_list, TRIM(header), it, TRIM(what))
    END IF

    CALL create_xml_vtu_scal_file(field, mesh, TRIM(ADJUSTL(file_list(rank+1))), field_name)

    IF (ALLOCATED(field_tmp)) DEALLOCATE(field_tmp)
  END SUBROUTINE make_vtu_file_2D

  SUBROUTINE check_list(communicator, file_list, check)
    IMPLICIT NONE
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: file_list
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: dummy_list
    INTEGER, DIMENSION(SIZE(file_list)) :: check_mylist
    INTEGER                             :: check, n, count
    !#include "petsc/finclude/petsc.h"
    MPI_Comm                            :: communicator
    PetscMPIInt                         :: rank, nb_procs
    PetscErrorCode                      :: ierr
    CALL MPI_Comm_rank(communicator, rank, ierr)
    CALL MPI_Comm_size(communicator, nb_procs, ierr)

    CALL MPI_ALLGATHER(check, 1, MPI_INTEGER, check_mylist, 1, &
         MPI_INTEGER, communicator, ierr)

    count = 0
    DO n = 1, SIZE(file_list)
       IF (check_mylist(n)==0) CYCLE
       count = count + 1
    END DO
    ALLOCATE(dummy_list(count))
    count = 0
    DO n = 1, SIZE(file_list)
       IF (check_mylist(n)==0) CYCLE
       count = count + 1
       dummy_list(count) = file_list(n)
    END DO
    DEALLOCATE(file_list)
    ALLOCATE(file_list(count))
    file_list = dummy_list
  END SUBROUTINE check_list

  SUBROUTINE create_pvd_file(file_list, file_header, time_step, what)
    IMPLICIT NONE
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: file_list
    CHARACTER(*),               INTENT(IN) :: file_header, what
    INTEGER,                    INTENT(IN) :: time_step
    INTEGER                                :: unit_file=789, j
    CHARACTER(len=5)                       :: tit, tit_part

    IF (what=='new') THEN
       OPEN (UNIT=unit_file, FILE=file_header//'.pvd', FORM = 'formatted', &
            ACCESS = 'append', STATUS = 'replace')
       WRITE(unit_file, '(A)') '<?xml version="1.0"?>'
       WRITE(unit_file, '(A)') '<VTKFile type="Collection" version="0.1"'// &
            ' byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
       WRITE(unit_file, '(A)') '<Collection>'
    ELSE
       OPEN (UNIT=unit_file, FILE=file_header//'.pvd', FORM = 'formatted', &
            ACCESS = 'append', STATUS = 'old')
       BACKSPACE(unit_file)
       BACKSPACE(unit_file)
    END IF
    WRITE(tit,'(I5)') time_step
    DO j = 1, SIZE(file_list)
       WRITE(tit_part,'(I5)') j
       WRITE(unit_file,'(A)') '<DataSet timestep="'//TRIM(ADJUSTL(tit))//&
            '" group="" part="'// &
            TRIM(ADJUSTL(tit_part))//'" file="./'//TRIM(ADJUSTL(file_list(j)))//&
            '.vtu'//'"/>'
    END DO
    WRITE(unit_file, '(A)') '</Collection>'
    WRITE(unit_file, '(A)') '</VTKFile>'

    CLOSE(unit_file) !Ecriture pour paraview
  END SUBROUTINE create_pvd_file

    SUBROUTINE create_xml_vtu_scal_file(field, mesh, file_name, field_name)
    USE def_type_mesh
    !USE input_data
    USE zlib_base64
    IMPLICIT NONE
    TYPE(mesh_type)                        :: mesh
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: field
    CHARACTER(*),               INTENT(IN) :: file_name, field_name
    INTEGER                                :: unit_file=789, m, n, type_cell
    REAL(KIND=4),    DIMENSION(3*mesh%np)              :: r4_threed_xml_field
    INTEGER(KIND=4), DIMENSION(mesh%gauss%n_w*mesh%me) :: i4_threed_xml_field
    INTEGER(KIND=4), DIMENSION(mesh%me)                :: i4_xml_field
    INTEGER(KIND=1), DIMENSION(mesh%me)                :: i1_xml_field
    CHARACTER(LEN=200)                         :: ascii_or_binary

    IF (SIZE(field)==0) RETURN

    IF (if_xml) THEN
       ascii_or_binary = 'binary'
       OPEN (UNIT=unit_file, FILE=file_name//'.vtu', STATUS = 'unknown')
       WRITE(unit_file,'(A)') '<?xml version="1.0" ?>'
       WRITE(unit_file,'(A)', advance="no") '<VTKFile type="UnstructuredGrid" version="0.1" '
       WRITE(unit_file,'(A)') 'compressor="vtkZLibDataCompressor" byte_order="LittleEndian">'
    ELSE
       ascii_or_binary = 'ascii'
       OPEN (UNIT=unit_file, FILE=file_name//'.vtu',&
            FORM = 'formatted', STATUS = 'unknown')
       WRITE(unit_file,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1"'// &
            ' byte_order="LittleEndian">'
    END IF

    WRITE(unit_file,'(A)') '<UnstructuredGrid>'
    WRITE(unit_file,'(A,I9,A,I9,A)') '<Piece NumberOfPoints="', mesh%np, &
         '" NumberOfCells="', mesh%me, '">'

    !===PointData Block
    WRITE(unit_file,'(A)') '<PointData Scalars="truc">'
    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'&
         //TRIM(ADJUSTL(field_name))//'" format="'//TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (if_xml) THEN
       CALL write_compressed_block(unit_file,  REAL(field(1:mesh%np),4))
    ELSE
       DO n = 1, mesh%np
          WRITE(unit_file,'(e14.7)') field(n)
       ENDDO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</PointData>'
    !===End of PointData Block

    !===CellData Block
    WRITE(unit_file,'(A)') '<CellData>'
    WRITE(unit_file,'(A)') '</CellData>'
    !===End of CellData Block

    !===Points Block
    WRITE(unit_file,'(A)') '<Points>'
    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="Points" '//&
         'NumberOfComponents="3" format="'//TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (if_xml) THEN
       DO n = 1, mesh%np
          r4_threed_xml_field(3*(n-1)+1) = REAL(mesh%rr(1,n),4)
          r4_threed_xml_field(3*(n-1)+2) = REAL(0.,4)
          r4_threed_xml_field(3*(n-1)+3) = REAL(mesh%rr(2,n),4)
       END DO
       CALL write_compressed_block(unit_file, r4_threed_xml_field(1:3*mesh%np))
    ELSE
       DO n = 1, mesh%np
          WRITE(unit_file,'(3(e14.7,x))') mesh%rr(1,n), 0.d0 , mesh%rr(3,n)
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Points>'
    !===End of Points Block

    !===Cells Block
    IF (mesh%gauss%n_w==3) THEN
       type_cell = 5
    ELSE IF (mesh%gauss%n_w==6) THEN
       type_cell = 22
    END IF
    WRITE(unit_file,'(A)') '<Cells>'
    WRITE(unit_file,'(A)') '<DataArray type="Int32" Name="connectivity" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (if_xml) THEN
       DO m = 1, mesh%me
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+1) = INT(mesh%jj(1,m)-1,4)
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+2) = INT(mesh%jj(2,m)-1,4)
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+3) = INT(mesh%jj(3,m)-1,4)
          IF (type_cell==22) THEN
             i4_threed_xml_field(mesh%gauss%n_w*(m-1)+4) = INT(mesh%jj(6,m)-1,4)
             i4_threed_xml_field(mesh%gauss%n_w*(m-1)+5) = INT(mesh%jj(4,m)-1,4)
             i4_threed_xml_field(mesh%gauss%n_w*(m-1)+6) = INT(mesh%jj(5,m)-1,4)
          END IF
       END DO
       CALL write_compressed_block(unit_file, i4_threed_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(3(I8,1x))') mesh%jj(1:3,m)-1
          IF (type_cell==22) THEN
             WRITE(unit_file,'(3(I8,1x))') mesh%jj(6,m)-1 , mesh%jj(4,m)-1 , mesh%jj(5,m)-1
          END IF
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'

    WRITE(unit_file,'(A)') '<DataArray type="Int32" Name="offsets" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (if_xml) THEN
       DO m = 1, mesh%me
          i4_xml_field(m) = INT(mesh%gauss%n_w*m,4)
       END DO
       CALL write_compressed_block(unit_file, i4_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(I8)') m*mesh%gauss%n_w
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'

    WRITE(unit_file,'(A)') '<DataArray type="UInt8" Name="types" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (if_xml) THEN
       DO m = 1, mesh%me
          i1_xml_field(m) = INT(type_cell,1)
       END DO
       CALL write_compressed_block(unit_file, i1_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(I8)') type_cell
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Cells>'
    !===End of Cells Block

    WRITE(unit_file,'(A)') '</Piece>'
    WRITE(unit_file,'(A)') '</UnstructuredGrid>'
    WRITE(unit_file,'(A)') '</VTKFile>'

    CLOSE(unit_file)
  END SUBROUTINE create_xml_vtu_scal_file
  
END MODULE plot_vtu_module
