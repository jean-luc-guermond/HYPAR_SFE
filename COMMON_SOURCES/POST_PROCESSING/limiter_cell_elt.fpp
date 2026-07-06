#:set list_name = list_limiters_input
#:set n_lim = n_lim_input
#:set syst_dim = syst_dim_input
#include "petsc/finclude/petsc.h"
MODULE  template_limiter_cell_elt
    USE template_cell_limiting_engine_parallel_module

CONTAINS

#:def limiter_spe(name, syst_dim)
    SUBROUTINE iterative_cell_limiting_procedure_${name}$(this, xx_in, loc_min)  
        USE petsc_tools
        USE petsc 
        USE compute_periodic
        USE my_util, ONLY: error_petsc
        IMPLICIT NONE
        CLASS(limiting_type),             INTENT(IN) :: this
        REAL(KIND=8), DIMENSION(:),                           INTENT(IN) :: loc_min
        REAL(KIND=8), DIMENSION(SIZE(loc_min),${syst_dim}$),  INTENT(INOUT) :: xx_in
        ! REAL(KIND=8), DIMENSION(k_dim+1, ${syst_dim}$, size(this%jj,2))     :: xx
        REAL(KIND=8), DIMENSION(SIZE(loc_min),${syst_dim}$) :: xx_out
        REAL(KIND=8), DIMENSION(${syst_dim}$)               :: uk_minus, uk_plus, uu, pp
        LOGICAL,      DIMENSION(SIZE(loc_min))              :: mask_up_vec, mask_down_vec
        REAL(KIND=8), DIMENSION(k_dim+1,${syst_dim}$)       :: xx_loc
        REAL(KIND=8), DIMENSION(k_dim+1) :: dummy, weight_mass_minus, weight_mass_plus, mass_loc
        INTEGER,      DIMENSION(k_dim+1) :: jloc
        INTEGER :: m, n, me, comp, ierr
        INTEGER, PARAMETER :: syst_size = ${syst_dim}$
        LOGICAL, DIMENSION(k_dim+1) :: mask_up, mask_down
        REAL(KIND=8) :: mass_plus, mass_minus, inv_minus, inv_plus, wm, wp, &
                lambda_K_minus, lambda_K_plus, &
                lambda_star_minus, lambda_star_plus, &
                lambda_minus, lambda_plus

        me = SIZE(this%jj,2)
        xx_out = 0.d0

        mask_up_vec   = psi_${name}$(xx_in,loc_min + this%epsilon*ABS(loc_min(:)))>0
        mask_down_vec = psi_${name}$(xx_in,loc_min - this%epsilon*ABS(loc_min))<0

        DO m = 1, me
            weight_mass_minus = 0.d0
            weight_mass_plus  = 0.d0
            lambda_minus = 1.d0
            lambda_plus = 1.d0
            uk_minus = 0.d0
            uk_plus = 0.d0
            mass_minus = 0.d0
            mass_plus = 0.d0
            DO n=1, k_dim+1
                jloc(n) = this%jj(n,m)
                mass_loc(n) = this%localized_mass(n,m)
                mask_up(n)   = mask_up_vec(jloc(n))
                mask_down(n) = mask_down_vec(jloc(n))
                DO comp=1, syst_size
                    xx_loc(n,comp) = xx_in(jloc(n),comp)
                END DO
                IF ( mask_down(n) ) THEN
                    wp = 0.d0
                    wm = mass_loc(n)
                    weight_mass_minus(n) = wm
                    mass_minus = mass_minus + wm
                    DO comp=1, syst_size
                        uk_minus(comp)=uk_minus(comp) + xx_loc(n,comp)*wm
                    END DO
                ELSE IF (mask_up(n)) THEN
                    wm = 0.d0
                    wp = mass_loc(n)
                    weight_mass_plus(n) = wp
                    mass_plus = mass_plus + wp
                    DO comp=1, syst_size
                        uk_plus(comp)=uk_plus(comp) + xx_loc(n,comp)*wp
                    END DO
                END IF
            END DO
            inv_minus = 1.d0/max(mass_minus,this%mass_eps)
            inv_plus  = 1.d0/max(mass_plus ,this%mass_eps)
            DO comp=1, syst_size
                uk_minus(comp) = uk_minus(comp)/max(mass_minus,this%mass_eps)
                uk_plus(comp) = uk_plus(comp)/max(mass_plus,this%mass_eps)
            END DO

            DO n=1, k_dim+1
                IF ( mask_down(n) ) THEN
                    DO comp=1, syst_size
                        uu(comp) = uk_plus(comp)
                        pp(comp) = xx_loc(n,comp)-uk_plus(comp)
                    END DO
                    lambda_minus = MIN(lambda_minus, scal_zero_of_psi_${name}$(loc_min(jloc(n)),uu,pp))
                ELSE IF (mask_up(n)) THEN
                    DO comp=1, syst_size
                        uu(comp) = xx_loc(n,comp)
                        pp(comp) = uk_minus(comp)-xx_loc(n,comp)
                    END DO
                    lambda_plus = MIN(lambda_plus, scal_zero_of_psi_${name}$(loc_min(jloc(n)),uu,pp))
                END IF
            END DO
            Lambda_star_minus = MAX(lambda_minus, 0.d0)
            Lambda_star_plus = MAX(lambda_plus, 0.d0)
            Lambda_K_minus = MAX(Lambda_star_minus, 1.d0-Lambda_star_plus*mass_plus/mass_minus)
            Lambda_K_plus  = MIN(Lambda_star_plus, (1.d0-Lambda_star_minus)*mass_minus/mass_plus)

            !=== DEBUGGING ===!
            !write(*,*)  'm possible limiting', m
            !write(*,*) lambda_minus,  lambda_plus
            !write(*,*) lambda_star_minus,  lambda_star_plus
            !write(*,*)  Lambda_K_minus, Lambda_K_plus
            !=== DEBUGGING ===!
            !    ! !!$ ===P2 fix
            !    ! IF (ABS(this%lumped_mass(jloc(n))).LE.this%mass_eps) THEN
            !    !    xx(n,m,:) = uk_plus(:)
            !    ! ELSE
            !    ! !!$ ===END fix
            DO comp=1, syst_size
                dummy = xx_loc(:,comp)*mass_loc  &
                        +weight_mass_minus(:)*(1-Lambda_K_minus)*(uk_plus(comp)-xx_loc(:,comp))&
                        +weight_mass_plus(:) *     Lambda_K_plus*(uK_minus(comp)-xx_loc(:,comp)) 
                xx_out(jloc,comp) = xx_out(jloc,comp) + dummy
                ! dummy  = limit_zero(:)*xx_loc(:,k) &
                ! +limit_minus(:)*(xx_loc(:,k)+(1-Lambda_K_minus)*(uk_plus(k)-xx_loc(:,k)))&
                ! +limit_plus(:) *(xx_loc(:,k)+     Lambda_K_plus*(uK_minus(k)-xx_loc(:,k)))
                ! xx_out(jloc,k) = xx_out(jloc,k) + dummy*this%localized_mass(:,m)
            END DO
        END DO

    !===Now we average over the nodes=========

        xx_in = xx_out
        DO comp = 1, syst_size
            CALL this%mesh%reduce_through_ghost(xx_in(:,comp), MPI_SUM)
            xx_in(:,comp) = xx_in(:,comp)/this%lumped_mass

            ! CALL VecSetValues(this%xvect1, np, this%LA%loc_to_glob(1, :) - 1, xx_out(:,k), ADD_VALUES, ierr)
            ! CALL VecAssemblyBegin(this%xvect1, ierr)
            ! CALL VecAssemblyEnd(this%xvect1, ierr)
            ! CALL VecZeroEntries(this%xvect1, ierr)
            ! xx_out(:,k) = xx_out(:,k)/this%lumped_mass
            CALL array_to_petsc_vec(xx_in(:,comp), this%xvect1, this%LA, 'insert')
            CALL periodic_add_vector_petsc(this%per%nb_bords, this%per%list, this%per%perlist, this%xvect1, this%LA)
            CALL extract_through_ghost(this%xvect1, 1, 1, this%LA, xx_in(:,comp), opt_assemble=.TRUE.)
        END DO    
    CONTAINS
$:limiter_func(name)
    END SUBROUTINE iterative_cell_limiting_procedure_${name}$

#:enddef

#:for name in list_name
$:limiter_spe(name, syst_dim)
#:endfor

END MODULE template_limiter_cell_elt
