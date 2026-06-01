MODULE tab_subcell_iso_refinement_module

    IMPLICIT NONE
!========================!
!=== ORDER 2 SUBCELLS ===!
!========================!
    INTEGER, TARGET :: jj_HL_2(4, 3) = transpose(reshape([ &
        4, 5, 6, &
        1, 5, 6, &
        2, 4, 6, &
        3, 4, 5], [3, 4]))
    INTEGER, TARGET :: jjs_HL_2(2, 2) = transpose(reshape([ &
        1, 3, &
        2, 3], [2, 2]))

!========================!
!=== ORDER 3 SUBCELLS ===!
!========================!
    INTEGER, TARGET :: jj_HL_3(9,3) = transpose(reshape([ &
        1, 6, 8, &
        6, 8,10, &
        8, 9,10, &
        4, 9,10, &
        2, 4, 9, &
        6, 7,10, &
        5, 7,10, &
        4, 5,10, &
        3, 5, 7 ], [3, 9]))
    INTEGER, TARGET :: jjs_HL_3(3,2) = transpose(reshape([ &
        1, 3, &
        3, 4, &
        2, 4], [2, 3]))

!========================!
!=== ORDER 4 SUBCELLS ===!
!========================!
    INTEGER, TARGET :: jj_HL_4(16,3) = transpose(reshape([ &
        1,  7,  10, &
        7,  10, 15, &
        10, 11, 15, &
        11, 14, 15, &
        11, 12, 14, &
        4,  12, 14, &
        2,  4,  12, &
        7,  8,  15, &
        8,  13, 15, &
        13, 14, 15, &
        5,  13, 14, &
        4,  5,  14, &
        8,  9,  13, &
        6,  9,  13, &
        5,  6,  13, &
        3,  6,  9], [3, 16]))
    INTEGER, TARGET :: jjs_HL_4(4,2) = transpose(reshape([ &
        1, 3, &
        3, 4, &
        4, 5, &
        2, 5], [2, 4]))

    PUBLIC :: build_jj_HL, build_jjs_HL, jj_HL_2, jjs_HL_2, jj_HL_3, jjs_HL_3, jj_HL_4, jjs_HL_4

CONTAINS

    SUBROUTINE build_jj_HL(refine_factor, jj_HL)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: refine_factor
        INTEGER :: n3, i1, i2, i3, a, num_cell
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: jj_HL
        INTEGER, DIMENSION(3,3) :: n
        INTEGER, DIMENSION(3)   :: i_tuple

        num_cell = 0
        ALLOCATE(jj_HL(refine_factor**2, 3), SOURCE=-1)
        DO n3=0, refine_factor-1
            IF (n3==0) THEN
                i1 = 1
            ELSE
                i1 = (3 + 1*(refine_factor-1)) + n3
            END IF
            IF (n3 == refine_factor-1) THEN
                i2 = 3
            ELSE
                i2 = (3 + 1*(refine_factor-1)) + n3+1
            END IF
            n(1, :) = get_n(i1)
            n(2, :) = get_n(i2)
            DO a=1, 1+2*(refine_factor-n3-1)
                num_cell = num_cell + 1
                IF (n(2, 3) < n(1, 3)) THEN
                    n(3, 1) = n(2, 1)-1
                    n(3, 2) = n(2, 2)
                    n(3, 3) = n(2, 3)+1
                ELSE IF (n(2, 3) > n(1, 3)) THEN
                    n(3, 1) = n(1, 1)-1
                    n(3, 2) = n(1, 2)+1
                    n(3, 3) = n(1, 3)
                END IF
                i3 = get_i(n(3, :))
                i_tuple = sort_i(i1, i2, i3)
                jj_HL(num_cell, :) = i_tuple
                n(1, :) = n(2, :)
                i1 = i2
                n(2, :) = n(3, :)
                i2 = i3
            END DO
        END DO

    CONTAINS

        !=== get barycentric coordinates knowing local index of node
        FUNCTION get_n(i_in) RESULT(n_array)
            USE my_util, ONLY: error_petsc, to_str
            IMPLICIT NONE
            INTEGER, INTENT(IN)   :: i_in
            INTEGER :: edge, i_loc, n1, n2, n3
            INTEGER, DIMENSION(3) :: n_array, tab_123

            tab_123 = [1, 2, 3]
            n1 = -1
            n2 = -1
            n3 = -1
            ! Case where the point is on summits
            IF (i_in <= 3) THEN
                n_array = 0
                n_array(i_in) = refine_factor
            ! Case where the point is on an external edge
            ELSE IF (i_in <= 3+3*(refine_factor-1)) THEN
                n_array = 0
                edge = (i_in-3 - 1) / (refine_factor-1) + 1
                i_loc = (i_in - 3) - (edge-1)*(refine_factor-1)
                n1 = MINVAL(PACK(tab_123, tab_123/=edge))
                n2 = MAXVAL(PACK(tab_123, tab_123/=edge))
                n_array(n1) = refine_factor - i_loc
                n_array(n2) = i_loc
            ! Case where the point is internal (order is lexicographic)
            ELSE IF (i_in <= (refine_factor + 1) * (refine_factor + 2) / 2) THEN
                i_loc = 3+3*(refine_factor-1)
                DO n1=1, refine_factor - 1
                    DO n2=1, refine_factor - 1
                        DO n3=1, refine_factor - 1
                            IF (n1+n2+n3/=refine_factor) CYCLE !only retain barycentric coordinates inside the triangle (hence loops starting from 1 and not 0)
                            i_loc = i_loc + 1
                            IF (i_loc == i_in) THEN
                                n_array = [n1, n2, n3]
                                RETURN
                            END IF
                        END DO
                    END DO
                END DO
            ELSE
                CALL error_petsc("BUG in get_n, i is too high"//to_str(i_in)//"; should be at most"//&
                &to_str((refine_factor + 1) * (refine_factor + 2) / 2))
            END IF
        END FUNCTION get_n

        !=== get local index of node knowing barycentric coordinates
        FUNCTION get_i(n_in) RESULT(i_out)
            USE my_util, ONLY: error_petsc, to_str
            IMPLICIT NONE
            INTEGER, DIMENSION(3), INTENT(IN) :: n_in
            INTEGER                           :: i, i_out

            DO i_out=1, (refine_factor + 1) * (refine_factor + 2) / 2
                IF (MAXVAL(ABS(get_n(i_out)-n_in))==0) RETURN
            END DO
            CALL error_petsc("BUG in get_i => could not find i for n_in = "//to_str(n_in))
        END FUNCTION get_i

        !=== naive sorting of 3 indices array
        FUNCTION sort_i(i1, i2, i3) RESULT(tab_i)
            USE my_util, ONLY: error_petsc, to_str
            IMPLICIT NONE
            INTEGER, INTENT(IN)   :: i1, i2, i3
            INTEGER, DIMENSION(3) :: tab_i
            IF (i1 > i2) THEN
                IF (i2 > i3) THEN
                    tab_i = [i3, i2, i1]
                ELSE
                    IF (i3 > i1) THEN
                        tab_i = [i2, i1, i3]
                    ELSE
                        tab_i = [i2, i3, i1]
                    END IF
                END IF
            ELSE
                IF (i1 > i3) THEN
                    tab_i = [i3, i1, i2]
                ELSE
                    IF (i3 > i2) THEN
                        tab_i = [i1, i2, i3]
                    ELSE
                        tab_i = [i1, i3, i2]
                    END IF
                END IF
            END IF

            IF (MINVAL(ABS(tab_i(1:2)-tab_i(2:3)))==0) THEN
                CALL error_petsc("BUG in sort_i: some are duplicates "//to_str(tab_i))
            END IF

        END FUNCTION sort_i

    END SUBROUTINE build_jj_HL

    SUBROUTINE build_jjs_HL(refine_factor, jjs_HL)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: refine_factor
        INTEGER :: n
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: jjs_HL

        ALLOCATE(jjs_HL(refine_factor, 2))
        jjs_HL(1, :)            = [1, 3]
        jjs_HL(refine_factor, :) = [2, refine_factor+1]
        DO n=1, refine_factor-2
            jjs_HL(n+1, :) = [n+2, n+3]
        END DO

    END SUBROUTINE build_jjs_HL

END MODULE tab_subcell_iso_refinement_module
