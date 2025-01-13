MODULE IDP_euler_boundary_conditions
  USE input_data
  PUBLIC :: IDP_euler_bc
  PRIVATE
CONTAINS
  SUBROUTINE IDP_euler_bc(unext,time)
    USE mesh_handling
    USE euler_bc_arrays
    USE boundary_conditions
    USE space_dim
    USE riemann_solution_module
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np,inputs%syst_size) :: unext
    REAL(KIND=8), DIMENSION(SIZE(udotn_js_D)) :: mdotn
    REAL(KIND=8), DIMENSION(SIZE(DIR_js_D)) :: rho, p, a, s, udotn, rhoD, pD, aD, sD, udotnD, &
         R1, R1D, R3, R3D, l1, l2, l3
    REAL(KIND=8), DIMENSION(k_dim,SIZE(DIR_js_D)) :: vel, vel_perp, velD, vel_perpD, rrD
    REAL(KIND=8), DIMENSION(k_dim+2,SIZE(DIR_js_D)) :: uD, uo
    REAL(KIND=8):: time, x, vn, tol ,rho0, u0, p0,lambda_max,pstar,v11,v32, ad_small
    INTEGER :: i, n, k, it, it_max=1
    unext(rho_js_D,1) = sol_anal(1,mesh%rr(:,rho_js_D),time)
    unext(rho_js_D,inputs%syst_size) = sol_anal(inputs%syst_size,mesh%rr(:,rho_js_D),time)
    SELECT CASE(k_dim)
    CASE(1)
       unext(ux_js_D,2)  = sol_anal(2,mesh%rr(:,ux_js_D),time)
    CASE(2)
       IF (size(udotn_js_D).NE.0) THEN
          mdotn = surf_normal_vtx(1,surf_udotn_js_D)*unext(udotn_js_D,2) &
               +  surf_normal_vtx(2,surf_udotn_js_D)*unext(udotn_js_D,3)
          unext(udotn_js_D,2)  = unext(udotn_js_D,2) - mdotn*surf_normal_vtx(1,surf_udotn_js_D)
          unext(udotn_js_D,3)  = unext(udotn_js_D,3) - mdotn*surf_normal_vtx(2,surf_udotn_js_D)
       END IF
       IF (SIZE(ux_js_D).NE.0) THEN
          unext(ux_js_D,2)  = sol_anal(2,mesh%rr(:,ux_js_D),time)
       END IF
       IF (SIZE(uy_js_D).NE.0) THEN
          unext(uy_js_D,3)  = sol_anal(3,mesh%rr(:,uy_js_D),time)
       END IF
       IF (SIZE(DIR_js_D).NE.0) THEN
          !===Dirichlet data
          rrD = mesh%rr(:,DIR_js_D)
          DO k = 1, k_dim+2
             !uD(k,:) = sol_anal(k,mesh%rr(:,DIR_js_D),time,opt_out=.true.)
             uD(k,:) = sol_anal(k,mesh%rr(:,DIR_js_D),time)
          END DO
          rhoD = uD(1,:)
          DO k = 1, k_dim
             velD(k,:) = uD(k+1,:)/rhoD
          END DO
          pD = pressure(uD)
          aD = SQRT(gamma*pD/uD(1,:))
          sD = pD/rhoD**gamma
          DO n = 1, SIZE(DIR_js_D)
             udotnD(n) = SUM(velD(1:k_dim,n)*DIR_normal_vtx(:,n))
             vel_perpD(1:k_dim,n) = velD(1:k_dim,n) - udotnD(n)*DIR_normal_vtx(1:k_dim,n)
          END DO
          R1D = udotnD+2*aD/(gamma-1)
          R3D = udotnD-2*aD/(gamma-1)
          ad_small = maxVAL(ad)*1.d-9

          DO it = 1, it_max
             !DO k = 1, k_dim+2
             !   uo(k,:) = unext(DIR_js_D,k)
             !END DO
             !===Subsonic/supersonic BC
             rho = unext(DIR_js_D,1)
             p = pressure(unext(DIR_js_D,:))
             a = SQRT(gamma*p/rho)
             s = p/rho**gamma
             DO k = 1, k_dim
                vel(k,:) = unext(DIR_js_D,k+1)/rho
             END DO
             DO n = 1, SIZE(DIR_js_D)
                udotn(n) = SUM(vel(:,n)*DIR_normal_vtx(:,n))
                vel_perp(1:k_dim,n) = vel(1:k_dim,n) - udotn(n)*DIR_normal_vtx(1:k_dim,n)
             END DO
             R1 = udotn+2*a/(gamma-1)
             R3 = udotn-2*a/(gamma-1)

             l1 = udotn - a
             l2 = udotn
             l3 = udotn + a
             !l1 = udotnD - aD
             !l2 = udotnD
             !l3 = udotnD + aD

             IF (.true.) THEN
                !===Riemann pb
                tol = 1.d-10
                DO n = 1, SIZE(DIR_js_D)
                   i = DIR_js_D(n)
                   IF (l3(n) <  0.d0) THEN !===Supersonic inflow (should be Dirichlet)
                      unext(i,:) = uD(:,i)
                      !write(*,*) 'sup in', mesh%rr(:,i)
                   ELSE IF (l2(n) < ad_small) THEN !===Subsonic inflow
                      CALL lambda(gamma,tol,rho(n),udotn(n),p(n),rhoD(n),udotnD(n),pD(n),lambda_max,pstar,k,v11,v32)
                      CALL riemann_solution_at_zero(pstar,rho(n),udotn(n),p(n),rhoD(n),udotnD(n),pD(n),rho0,u0,p0)
                      unext(i,1) = rho0
                      unext(i,2:k_dim+1) = rho0*(vel_perpD(1:k_dim,n)+u0*DIR_normal_vtx(1:k_dim,n))
                      unext(i,k_dim+2) = p0/(gamma-1) + SUM(unext(i,2:k_dim+1)**2)/(2*rho0)
                      !write(*,*) 'sub in', mesh%rr(:,i)
                   ELSE IF (l1(n) < 0.d0) THEN !===Subsonic outflow
                      !HACK to be removed
                      !rhoD(n) = 0.2*rhoD(n)
                      CALL lambda(gamma,tol,rho(n),udotn(n),p(n),rhoD(n),udotnD(n),pD(n),lambda_max,pstar,k,v11,v32)
                      CALL riemann_solution_at_zero(pstar,rho(n),udotn(n),p(n),rhoD(n),udotnD(n),pD(n),rho0,u0,p0)
                      unext(i,1) = rho0
                      unext(i,2:k_dim+1) = rho0*(vel_perp(1:k_dim,n)+u0*DIR_normal_vtx(1:k_dim,n))
                      unext(i,k_dim+2) = p0/(gamma-1) + SUM(unext(i,2:k_dim+1)**2)/(2*rho0)
                      !write(*,*) 'sub out', mesh%rr(:,i)
                   END IF
                END DO
             ELSE
                !===Postprocessing
                DO n = 1, SIZE(DIR_js_D)
                   i = DIR_js_D(n)
                   IF (l3(n) < 0.d0) THEN !===Supersonic inflow (should be Dirichlet)
                      unext(i,:) = uD(:,i)
                   ELSE IF (l2(n) < ad_small) THEN !===Subsonic inflow
                      x = ((((gamma-1)/4)*(R3D(n)-R1(n)))**2/(gamma*sD(n)))**(1/(gamma-1))
                      unext(i,1) = x
                      p(n) = sD(n)*x**gamma
                      vn = (R3D(n)+R1(n))/2
                      unext(i,2:k_dim+1) = x*(vel_perpD(1:k_dim,n)+vn*DIR_normal_vtx(1:k_dim,n))
                      unext(i,k_dim+2) = p(n)/(gamma-1) + SUM(unext(i,2:k_dim+1)**2)/(2*x)
                   ELSE IF (l1(n) < 0.d0) THEN !===Subsonic outflow
                      x = ((((gamma-1)/4)*(R3D(n)-R1(n)))**2/(gamma*s(n)))**(1/(gamma-1))
                      unext(i,1) = x
                      p(n) = s(n)*x**gamma
                      vn = (R3D(n)+R1(n))/2
                      unext(i,2:k_dim+1) = x*(vel_perp(1:k_dim,n)+vn*DIR_normal_vtx(1:k_dim,n))
                      unext(i,k_dim+2) = p(n)/(gamma-1) + SUM(unext(i,2:k_dim+1)**2)/(2*x)
                      !WRITE(*,*) 'Outlow subsonic', mesh%rr(1,i), mesh%rr(2,i), l1(n), l2(n)
                    ELSE
                      !WRITE(*,*) 'Outlow supersonic', mesh%rr(1,i), mesh%rr(2,i)
                   END IF
                END DO
             END IF
          END DO
       END IF
    END SELECT
  END SUBROUTINE IDP_euler_bc
END MODULE IDP_euler_boundary_conditions
