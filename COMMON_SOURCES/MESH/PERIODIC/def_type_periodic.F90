MODULE def_type_periodic
   USE dyn_line_type
   TYPE periodic_type
      TYPE(dyn_int_line), DIMENSION(20) :: list
      TYPE(dyn_int_line), DIMENSION(20) :: perlist
      INTEGER, POINTER, DIMENSION(:) :: pnt
      INTEGER :: n_bord
   END TYPE periodic_type
END MODULE def_type_periodic