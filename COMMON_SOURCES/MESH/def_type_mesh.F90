!
!Authors: Jean-Luc Guermond, Copyright 2000
!
MODULE def_type_mesh
   USE dyn_line_type
   USE space_dim
   USE periodic_data_module
   USE mesh_data_module, ONLY: mesh_info_type
#include "petsc/finclude/petsc.h"
   USE petsc
   IMPLICIT NONE

   TYPE aij_type
      INTEGER, POINTER, DIMENSION(:) :: ia, ja
   END TYPE aij_type

   TYPE petsc_csr_LA
      INTEGER, DIMENSION(:), POINTER :: ia, ja
      INTEGER, DIMENSION(:, :), POINTER :: loc_to_glob
      INTEGER :: kmax
      INTEGER, DIMENSION(:), POINTER :: np
      INTEGER, DIMENSION(:), POINTER :: dom_np
   END TYPE petsc_csr_LA

   !------------------------------------------------------------------------------
   !  REAL(KIND=8), DIMENSION(n_w,  l_G),  PUBLIC :: ww
   !  REAL(KIND=8), DIMENSION(n_ws, l_Gs), PUBLIC :: wws
   !  REAL(KIND=8), DIMENSION(k_d,  n_w,  l_G,   me),   PUBLIC :: dw
   !  REAL(KIND=8), DIMENSION(n_w,l_G,1:2,me),          PUBLIC :: dwni !d/dn, interface (JLG, April 2009)
   !  REAL(KIND=8), DIMENSION(k_d,        l_Gs,  mes),  PUBLIC :: rnorms
   !  REAL(KIND=8), DIMENSION(l_G,   me),               PUBLIC :: rj
   !  REAL(KIND=8), DIMENSION(l_Gs,  mes),              PUBLIC :: rjs
   !  REAL(KIND=8), DIMENSION(k_d, n_w, l_Gs, mes)             :: dw_s !gradient at the boundary
   !------------------------------------------------------------------------------

   TYPE gauss_type
      INTEGER :: k_d, n_w, l_G, n_ws, l_Gs, n_e
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: ww   !Shape functions on cells
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: wws  !Shape function at boundary or DG interface
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: wwsi !Interface shape function (JLG, June 4 2012)
      REAL(KIND = 8), DIMENSION(:, :, :, :), POINTER :: dw !Gradient on cells
      REAL(KIND = 8), DIMENSION(:, :, :, :), POINTER :: dw_s !Gradient at the boundary
      REAL(KIND = 8), DIMENSION(:, :, :, :), POINTER :: dwni !Interface gradient (JLG, April 2009)
      REAL(KIND = 8), DIMENSION(:, :, :), POINTER :: rnorms
      REAL(KIND = 8), DIMENSION(:, :, :), POINTER :: rnorms_v !(JLG Aug 31, 2017)
      REAL(KIND = 8), DIMENSION(:, :, :, :), POINTER :: rnormsi !Interface normal (JLG, June 4 2012)
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: rj   !Interface weight (JLG, April 2009)
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: rji
      REAL(KIND = 8), DIMENSION(:, :), POINTER :: rjs
      REAL(KIND = 8), DIMENSION(:, :, :, :), POINTER :: dwps !special!
      REAL(KIND = 8), DIMENSION(:, :, :, :), POINTER :: dws  !SPECIAL!
   END TYPE gauss_type

   !------------------------------------------------------------------------------
   !  loc_to_glob(np)   gives global numbering from local numbering on current processor
   !  jj(n_w,   me)     nodes of the  volume_elements
   !  jji(n_w, 1:2, mi) edge to node conectivity array --> volume numbering (JLG April 2009)
   !  neighi(1:2, mi)   interfaces to volume elements --> cell 1 has lowest cell number
   !  jjsi(n_ws, mi)    nodes of the interface elements --> volume numbering (JLG April 2009)
   !  jjs(n_ws, mes)    nodes of the surface_elements --> volume numbering
   !  iis(n_ws, mes)    nodes of the surface_elements --> surface numbering
   !  mm(me)           (possibly sub) set of elements for quadrature
   ! mms(mes)          (possibly sub) set of surface_elements for surf_quadrature
   !------------------------------------------------------------------------------

   TYPE mesh_type
      INTEGER, POINTER, DIMENSION(:, :) :: jj, jjs, iis
      INTEGER, POINTER, DIMENSION(:, :) :: jj_extra, jce_extra, jjs_extra !(extra layer of cells not own by proc but with dofs own by proc)
      INTEGER, POINTER, DIMENSION(:, :) :: jjs_int
      INTEGER, POINTER, DIMENSION(:) :: jcc_extra
      !=== SIZE(nt, mesh%me) => !!!GLOBAL!!! numbering of edges (use jce_loc to get local numbering)
      INTEGER, POINTER, DIMENSION(:, :) :: jce! cell-> edge (JLG+MC Sept 2022)
      INTEGER, POINTER, DIMENSION(:) :: jees, jecs !edges belonging to another proc (MC Sept 2022)
      INTEGER, POINTER, DIMENSION(:, :, :) :: jji  ! (JLG April 2009)
      INTEGER, POINTER, DIMENSION(:, :) :: jjsi ! (JLG April 2009)
      INTEGER, POINTER, DIMENSION(:) :: j_s  ! boundary nodes --> volume numbering
      REAL(KIND = 8), POINTER, DIMENSION(:, :) :: rr
      REAL(KIND = 8), POINTER, DIMENSION(:, :, :) :: rrs_extra  ! coordinates for cells at interfaces
      INTEGER, POINTER, DIMENSION(:, :) :: neigh
      INTEGER, POINTER, DIMENSION(:, :) :: neighi ! (JLG April 2009)
      INTEGER, POINTER, DIMENSION(:) :: sides, neighs, sides_extra, neighs_extra !interfaces
      INTEGER, POINTER, DIMENSION(:) :: sides_int
      INTEGER, POINTER, DIMENSION(:, :) :: neighs_int
      INTEGER, POINTER, DIMENSION(:) :: i_d
      !==Parallel structure
      !=== SIZE(mesh%np)
      INTEGER, POINTER, DIMENSION(:) :: loc_to_glob ! (JLG+FL, January 2011)
      !=== SIZE(nb_proc), SIZE(nb_proc + 1)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: domnp, disp ! (JLG+FL, January 2011) resp. ALLGATHER(mesh%dom_np) and cumsum starting from 1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: domcell, discell ! (MC, Sept 2022) resp. ALLGATHER(mesh%me) and cumsum starting from 1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: disedge, domedge ! (MC, Sept 2022) resp. ALLGATHER(mesh%medge) and cumsum starting from 1
      INTEGER :: dom_me, dom_np, dom_mes ! (JLG+FL, January 2011)
      !==Isolated nodes at interfaces
      INTEGER, POINTER, DIMENSION(:) :: isolated_jjs !give glob index of isolated point
      INTEGER, POINTER, DIMENSION(:, :) :: isolated_interfaces !give the number of the interfaces
      INTEGER :: nis !number of isolated points
      ! dom_me and dom_mes are obsolete structures.
      ! dom_np is the number of nodes owned by the processor: dom_np .LE. mesh%np
      !==End parallel structure
      INTEGER :: me, mes, np, nps, mi, medge, medges, mextra, mes_extra, mes_int
      LOGICAL :: edge_stab ! edge stab, yes/no, (JLG April 2009)
      TYPE(gauss_type) :: gauss
      REAL(KIND = 8), POINTER, DIMENSION(:) :: hloc ! local mesh size (JLG+LC January, 21, 2015)
      REAL(KIND = 8), POINTER, DIMENSION(:) :: hloc_gauss ! local mesh size (JLG+LC January, 21, 2015)
      REAL(KIND = 8) :: global_diameter !diameter of domain (LC 2017/01/27)
      REAL(KIND = 8), POINTER, DIMENSION(:) :: hm !local meshsize in azimuth (JLG April 7, 2017)
      TYPE(periodic_type) :: per !<==Periodic structure is attached to the mesh

      ! !VB 14/05/2026 => array o size np - dom_np which value is the proc to which the node belongs
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: proc_np_loc 
      INTEGER                       :: rank, proc, nb_proc !VB 11/05/2026
      CHARACTER(LEN=:), ALLOCATABLE :: name
      TYPE(mesh_info_type)          :: info
      MPI_Comm, POINTER             :: comm !VB 11/05/2026
   CONTAINS
      PROCEDURE :: jj_glob
      PROCEDURE :: jce_loc
      PROCEDURE :: attr_e
      PROCEDURE :: side_edge
      PROCEDURE :: create_comm, gather_dom_np, gather_me, gather_medge
      PROCEDURE :: get_proc, global_numbering, build_loc_to_glob
   END TYPE mesh_type

   TYPE mesh_type_interface
      INTEGER, POINTER, DIMENSION(:) :: slave_elem ! list slave elemt in interface
      INTEGER, POINTER, DIMENSION(:) :: list_slave_node ! list of slave nodes on interface
      INTEGER, POINTER, DIMENSION(:, :) :: master_node ! local --> global numbering; master nodes
      INTEGER, POINTER, DIMENSION(:, :) :: slave_node  ! local --> global numbering; slave nodes
      INTEGER :: me ! nb of slave elemt in interface
   END TYPE mesh_type_interface

   TYPE mesh_type_boundary
      INTEGER, POINTER, DIMENSION(:) :: master ! list master boundary elemts
      INTEGER, POINTER, DIMENSION(:) :: slave ! list slave boundary elemts not in interface
      INTEGER, POINTER, DIMENSION(:) :: INTERFACE ! list slave boundary elemts in the interface
      INTEGER, POINTER, DIMENSION(:, :) :: master_node ! local --> global numbering; master nodes
   END TYPE mesh_type_boundary

   TYPE interface_type
      INTEGER :: mes ! number of interface elements
      INTEGER :: mes_extra ! number of interface elements on extra cells
      INTEGER, POINTER, DIMENSION(:) :: mesh1 ! list slave interface elements
      INTEGER, POINTER, DIMENSION(:) :: mesh2 ! list master interface elements
      INTEGER, POINTER, DIMENSION(:, :) :: jjs1 ! list of slave node on interface elements
      INTEGER, POINTER, DIMENSION(:, :) :: jjs2 ! list of master nodes on interface elements
      INTEGER, POINTER, DIMENSION(:) :: mesh1_extra ! list slave interface elements on extra cells
      INTEGER, POINTER, DIMENSION(:) :: mesh2_extra! list master interface elements on extra cells
      INTEGER, POINTER, DIMENSION(:, :) :: jjs1_extra ! list of slave node on interface elements on extra cells
      INTEGER, POINTER, DIMENSION(:, :) :: jjs2_extra ! list of master nodes on interface elements on extra cells
   END TYPE interface_type

CONTAINS

   FUNCTION jce_loc(this, n, m) RESULT(out)
      CLASS(mesh_type) :: this
      INTEGER :: n, m, out
      out = this%jce(n, m) - this%disedge(this%rank + 1) + 1
   END FUNCTION jce_loc

   FUNCTION side_edge(this, n, m) RESULT(out)
      CLASS(mesh_type) :: this
      INTEGER :: n, m
      LOGICAL :: out
      SELECT CASE(k_dim)
      CASE(1)
         out = this%neigh(1, m) == 0 .or. this%neigh(2, m) == 0
      CASE(2)
         out = this%neigh(n, m) == 0
      CASE DEFAULT
         WRITE(*, *) 'space dimension is not implemented for this k_dim in side_edge in mesh_type'
      END SELECT

   END FUNCTION side_edge


   FUNCTION jj_glob(this, n, m) RESULT(out)
      CLASS(mesh_type) :: this
      INTEGER :: n, m, out
      out = this%loc_to_glob(this%jj(n, m))
   END FUNCTION jj_glob


   FUNCTION attr_e(this, e) RESULT(out)
      CLASS(mesh_type) :: this
      INTEGER :: e
      LOGICAL :: out
      out = this%disedge(this%rank + 1) <= e .AND. e < this%disedge(this%rank + 2)
   END FUNCTION attr_e

!==========================================
!=== communication subroutines VB 11/05/2026
!==========================================

   SUBROUTINE create_comm(this, communicator)
      IMPLICIT NONE
      CLASS(mesh_type) :: this
      INTEGER          :: ierr
      MPI_Comm, TARGET :: communicator

      this%comm => communicator
      CALL MPI_Comm_rank(this%comm, this%rank, ierr)
      this%proc = this%rank + 1
      CALL MPI_Comm_Size(this%comm, this%nb_proc, ierr)
   END SUBROUTINE create_comm

   SUBROUTINE gather_dom_np(this)
      IMPLICIT NONE
      CLASS(mesh_type) :: this
      INTEGER          :: n, ierr

      IF (.NOT. ALLOCATED(this%disp)) THEN
         ALLOCATE(this%disp(this%nb_proc + 1), this%domnp(this%nb_proc))
      END IF
      CALL MPI_ALLGATHER(this%dom_np, 1, MPI_INTEGER, this%domnp, 1, &
           MPI_INTEGER, this%comm, ierr)
      this%disp(1) = 1
      DO n = 1, this%nb_proc
         this%disp(n + 1) = this%disp(n) + this%domnp(n)
      END DO
   END SUBROUTINE gather_dom_np

   SUBROUTINE gather_me(this)
      IMPLICIT NONE
      CLASS(mesh_type) :: this
      INTEGER          :: n, ierr

      IF (.NOT. ALLOCATED(this%discell)) THEN
         ALLOCATE(this%discell(this%nb_proc + 1), this%domcell(this%nb_proc))
      END IF
      CALL MPI_ALLGATHER(this%me, 1, MPI_INTEGER, this%domcell, 1, &
           MPI_INTEGER, this%comm, ierr)
      this%discell(1) = 1
      DO n = 1, this%nb_proc
         this%discell(n + 1) = this%discell(n) + this%domcell(n)
      END DO
   END SUBROUTINE gather_me

   SUBROUTINE gather_medge(this)
      IMPLICIT NONE
      CLASS(mesh_type) :: this
      INTEGER          :: n, ierr

      IF (.NOT. ALLOCATED(this%disedge)) THEN
         ALLOCATE(this%disedge(this%nb_proc + 1), this%domedge(this%nb_proc))
      END IF
      CALL MPI_ALLGATHER(this%medge, 1, MPI_INTEGER, this%domedge, 1, &
           MPI_INTEGER, this%comm, ierr)
      this%disedge(1) = 1
      DO n = 1, this%nb_proc
         this%disedge(n + 1) = this%disedge(n) + this%domedge(n)
      END DO

   END SUBROUTINE gather_medge

!==========================================
!=== functions to get proc number VB 14/05/2026
!==========================================

   FUNCTION get_proc(this, val_glob, char_in) RESULT(p)
      !> function to get the proc owning node/element/edge val_glob
      !! char_in: np, me, or medge 

      !! np necessarily comes from this%loc_to_glob, this%jj_extra, this%jjs_extra
      !! me necessarily comes from this%jcc_extra
      !! medge necessarily comes from this%jees, this%jce_extra, this%jce 
      USE my_util, ONLY: error_petsc, to_str
      IMPLICIT NONE
      CLASS(mesh_type) :: this
      INTEGER,          INTENT(IN) :: val_glob
      CHARACTER(LEN=*), INTENT(IN) :: char_in
      INTEGER, DIMENSION(:), ALLOCATABLE :: cumul_over_procs
      INTEGER :: p

      SELECT CASE(char_in)
      CASE('np')
         cumul_over_procs = this%disp
      CASE('me')
         cumul_over_procs = this%discell
      CASE('medge')
         cumul_over_procs = this%disedge
      CASE DEFAULT
         CALL error_petsc("BUG in get_proc => wrong char_in "//char_in//".&
         Should be in 'np; me; medge'")
      END SELECT

      DO p = 1, this%nb_proc
         IF (val_glob < cumul_over_procs(p + 1)) RETURN
      END DO

   END FUNCTION get_proc

!=================================================================================
!   Subroutine to create loc_to_glob, KNOWING jj and the cumulative quantities
!=================================================================================

   FUNCTION global_numbering(this, p, n_loc) RESULT(n_glob)
      IMPLICIT NONE
      CLASS(mesh_type)    :: this
      INTEGER, INTENT(IN) :: p, n_loc
      INTEGER             :: n_glob

      n_glob = n_loc + this%disp(p) - 1
   END FUNCTION global_numbering

   SUBROUTINE build_loc_to_glob(this)
      !> subroutine building the loc_to_glob array, provided the construction of:
      !! this%jj, this%proc_np_loc, this%proc, this%np, this%dom_np, cumulative quantities
      IMPLICIT NONE
      CLASS(mesh_type)    :: this
      LOGICAL, DIMENSION(:), ALLOCATABLE :: virgin
      INTEGER :: m, n, n_loc, p

      !=== Create loc_to_glob
      IF (.NOT. ASSOCIATED(this%loc_to_glob)) ALLOCATE(this%loc_to_glob(this%np), source=-1)

      != nodes owned by proc
      ALLOCATE(virgin(this%dom_np), source=.TRUE.)
      DO m=1, this%me
         DO n=1, SIZE(this%jj,1)
            n_loc = this%jj(n, m)
            IF (n_loc > this%dom_np) CYCLE
            this%loc_to_glob(n_loc) = this%global_numbering(this%proc, n_loc)
            virgin(n_loc) = .FALSE.
         END DO
      END DO
      IF (ANY(virgin)) THEN
         DO n=1, this%dom_np
            IF (virgin(n)) WRITE(*,*) n, 'is virgin (out of) ', this%dom_np, this%np
         END DO
         WRITE(*,*) 'BUG in def loc_to_glob: how can jj not have all values between 1 and dom_np??'
         STOP
      END IF

      != nodes owned by other proc
      DEALLOCATE(virgin)
      ALLOCATE(virgin(this%np-this%dom_np), source=.TRUE.)
      IF (SIZE(this%proc_np_loc,2)+this%dom_np /= this%np) THEN
         WRITE(*,*) 'sizes mismatch in proc_np_loc ',SIZE(this%proc_np_loc,2)+this%dom_np , this%np
         STOP
      END IF
      DO n=1, SIZE(this%proc_np_loc,2)
         p = this%proc_np_loc(1, n)
         n_loc = this%proc_np_loc(2, n)
         this%loc_to_glob(n+this%dom_np) = this%global_numbering(p, n_loc)
         virgin(n) = .FALSE.
      END DO

      IF (ANY(virgin)) THEN
         DO n=1, this%np-this%dom_np
            IF (virgin(n)) WRITE(*,*) n, ' is virgin (out of) ', this%dom_np, this%np
         END DO
         WRITE(*,*) 'BUG in def loc_to_glob: proc_np_loc does not seem to contain all ghost points for proc ', this%proc
         STOP
      END IF
   END SUBROUTINE build_loc_to_glob
END MODULE def_type_mesh
