#include "petsc/finclude/petsc.h"
MODULE  template_uij_bar_bounds
#:set list_name = list_limiters_input
#:set n_lim = n_lim_input
#:set syst_dim = syst_dim_input
    USE template_abstract_hyperbolic_module
    USE petsc
CONTAINS

#:def compute_bounds_uijbar(syst_dim, n_lim)
    SUBROUTINE template_compute_bounds_uijbar(this, flux_array, un, bounds)
        USE space_dim
        IMPLICIT NONE
        CLASS(hyperbolic_type) :: this
        REAL(KIND = 8), DIMENSION(this%mesh_L%np, k_dim, ${syst_dim}$), INTENT(IN) :: flux_array
        REAL(KIND = 8), DIMENSION(this%mesh_L%np, ${syst_dim}$),        INTENT(IN)  :: un
        REAL(KIND = 8), DIMENSION(this%mesh_L%np, ${n_lim}$),        INTENT(OUT) :: bounds
        REAL(KIND = 8), DIMENSION(1, ${syst_dim}$) :: uijbar
        REAL(KIND = 8), DIMENSION(1), PARAMETER   :: zero=0.d0
        REAL(KIND = 8), DIMENSION(k_dim)    :: nij_c_scal
        REAL(KIND = 8)                      :: norm_c_scal, max_lambda
        INTEGER :: nl, m, n, ni, nj, i, j, edge_idx, nw, comp

        ASSOCIATE(mesh => this%mesh_L,&
            lim_func => this%limiting_all_functionals%limiting_functionals, &
            max_lambda_array => this%xxx2, zero_vec => this%x1)

        !===================================================!
        !=== Computation of bounds for limiting method 1 ===!
        !===================================================!
            
        nw = mesh%gauss%n_w
        zero_vec = 0.d0
#:for nl, name in enumerate(list_name)
        bounds(:, ${nl+1}$) = psi_${name}$(un(:, :), zero_vec)
#:endfor

        DO edge_idx=1, this%n_edge_in
            m = this%edge_in(1, edge_idx)
            n = this%edge_in(2, edge_idx)

            ni = MOD(n, nw) + 1
            nj = MOD(n + 1, nw) + 1
            i = mesh%jj(ni, m)
            j = mesh%jj(nj, m)
            nij_c_scal(:) = this%matrices%nijL_loc_array(ni,nj,:,m)
            norm_c_scal = this%matrices%cijL_norm_loc_array(ni,nj,m)

            max_lambda = max_lambda_array(ni, nj, m)
            DO comp=1, ${syst_dim}$
                uijbar(1,comp) =  (un(i, comp)+un(j, comp))/2.d0 - &
                SUM((flux_array(j, :, comp)-flux_array(i, :, comp))*nij_c_scal(:))/(2.d0*max_lambda)
            END DO

#:for nl, name in enumerate(list_name)
            bounds(i:i, ${nl+1}$) = MIN(bounds(i, ${nl+1}$), psi_${name}$(uijbar, zero))
            bounds(j:j, ${nl+1}$) = MIN(bounds(j, ${nl+1}$), psi_${name}$(uijbar, zero))
#:endfor
        END DO

        DO edge_idx=1, this%n_edge_out
            m = this%edge_out(1, edge_idx)
            n = this%edge_out(2, edge_idx)

            ni = MOD(n, nw) + 1
            nj = MOD(n + 1, nw) + 1
            i = mesh%jj(ni, m)
            j = mesh%jj(nj, m)
            nij_c_scal(:) = this%matrices%nijL_loc_array(ni,nj,:,m)
            norm_c_scal = this%matrices%cijL_norm_loc_array(ni,nj,m)

            max_lambda = max_lambda_array(ni, nj, m)
            DO comp=1, ${syst_dim}$
                uijbar(1,comp) =  (un(i, comp)+un(j, comp))/2.d0 - &
                SUM((flux_array(j, :, comp)-flux_array(i, :, comp))*nij_c_scal(:))/(2.d0*max_lambda)
            END DO

#:for nl, name in enumerate(list_name)
            bounds(i:i, ${nl+1}$) = MIN(bounds(i, ${nl+1}$), psi_${name}$(uijbar, zero))
            bounds(j:j, ${nl+1}$) = MIN(bounds(j, ${nl+1}$), psi_${name}$(uijbar, zero))
#:endfor
        END DO
        DO nl=1, this%limiting_all_functionals%nl
            CALL mesh%reduce_through_ghost(bounds(:,nl), MPI_MIN)
        END DO
        END ASSOCIATE

    CONTAINS
!=== Force inlining
#:for name in list_name
$:limiter_func(name)
#:endfor
!=== Force inlining
    END SUBROUTINE template_compute_bounds_uijbar
#:enddef

$:compute_bounds_uijbar(syst_dim, n_lim)

END MODULE template_uij_bar_bounds
