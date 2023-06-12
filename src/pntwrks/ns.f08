module ns
implicit none
contains
subroutine slv_ns(vel,V,P,ro, pnts, totpnts, krnls, vel1, P0, &
    F, phi, viscosity, density, totphases, segma, bcs, &
     dt, av_prmtr, P_prmtr, lagrangian)

    use krnl_struct
    use bc_struct
    use bndry
    use intrf
    use omp_lib

    real(8), dimension(:,:), intent(IN)                ::  pnts
    integer, intent(IN)                                ::  totpnts
    integer, intent(IN)                                ::  totphases
    logical, intent(IN)                                ::  lagrangian
    real(8), intent(IN)                                ::  dt, P_prmtr, av_prmtr, segma
    real(8), dimension(:), intent(IN)                  ::  P0, phi, viscosity, density
    real(8), dimension(:,:), intent(IN)                ::  vel1, F
    type(bc), dimension(4), intent(IN)                 ::  bcs
    type(kernel), dimension(:), intent(IN)             ::  krnls
    integer                                            ::  iter
    integer                                            ::  vt
    real(8)                                            ::  eps, nu
    real(8), dimension(:), allocatable, intent(OUT)    ::  P, V, ro
    real(8), dimension(:,:), allocatable, intent(OUT)  ::  vel
    real(8), dimension(:,:), allocatable               ::  vel_inter
    real(8), dimension(:), allocatable                 ::  P_fut
    real(8), dimension(:), allocatable                 ::  curv
    real(8)                                            ::  adv_x, adv_y, adv_z, abs_dphi, alpha, c, d0, h0, pii0
    real(8)                                            ::  A_pressure_x, A_pressure_y, A_pressure_z
    real(8)                                            ::  art_visc_x, art_visc_y, art_visc_z, dif_x, dif_y, dif_z
    real(8)                                            ::  dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz
    real(8)                                            ::  nabla2_P, nabla2_vx, nabla2_vy, nabla2_vz
    real(8)                                            ::  tension_x, tension_y, tension_z, trm, trm0, trm_f
    real(8)                                            ::  uij_rij, xij, vxij, yij, vyij, zij, vzij, v0, Vab_x, Vab_y, Vab_z, vmax
    real(8)                                            ::  phi_diff, dphidx, dphidy, dphidz
    real(8)                                            ::  vx_diff, vy_diff, vz_diff
    integer                                            ::  i, j, s, nbr, pntnum

    !SOLVER BASED ON THE ARTIFICIAL COMPRESSIBILITY AND PROJECTION METHOD
    iter = 1;
    eps = 0.00001;

    !do not touch
    allocate(vel(totpnts, 3));
    allocate(P(totpnts));
    vel = vel1
    P = P0
    !do not touch

    allocate(V(totpnts))
    allocate(P_fut(totpnts))
    allocate(ro(totpnts))
    allocate(vel_inter(totpnts,3))
    adv_x = 0
    adv_y = 0
    adv_z = 0

    !============================================================================
    ! Computing intrmediate velocities
    !============================================================================
    !$omp parallel do private(i)
    do i = 1,totpnts
        dvxdx = 0.0
        dvxdy = 0.0
        dvxdz = 0.0
        dvydx = 0.0
        dvydy = 0.0
        dvydz = 0.0
        dvzdx = 0.0
        dvzdy = 0.0
        dvzdz = 0.0
        do j = 1, krnls(i)%totnbrs
            pntnum = krnls(i)%nbrs(j)
            vx_diff = vel(pntnum,1)-vel(i,1)
            dvxdx = dvxdx + krnls(i)%nx(j)*vx_diff
            dvxdy = dvxdy + krnls(i)%ny(j)*vx_diff
            dvxdz = dvxdz + krnls(i)%nz(j)*vx_diff

            vy_diff = vel(pntnum,2)-vel(i,2)
            dvydx = dvydx + krnls(i)%nx(j)*vy_diff
            dvydy = dvydy + krnls(i)%ny(j)*vy_diff
            dvydz = dvydz + krnls(i)%nz(j)*vy_diff

            vz_diff = vel(pntnum,3)-vel(i,3)
            dvzdx = dvzdx + krnls(i)%nx(j)*vz_diff
            dvzdy = dvzdy + krnls(i)%ny(j)*vz_diff
            dvzdz = dvzdz + krnls(i)%nz(j)*vz_diff
        end do
        nabla2_vx = 0.0
        nabla2_vy = 0.0
        nabla2_vz = 0.0
        do j = 1,krnls(i)%totnbrs
            pntnum = krnls(i)%nbrs(j)
            nabla2_vx = nabla2_vx + (vel(pntnum,1)-vel(i,1))*krnls(i)%nabla2(j)
            nabla2_vy = nabla2_vy + (vel(pntnum,2)-vel(i,2))*krnls(i)%nabla2(j)
            nabla2_vz = nabla2_vz + (vel(pntnum,3)-vel(i,3))*krnls(i)%nabla2(j)
        end do
        if (lagrangian .eqv. .false.) then
            !===========================================================
            ! computing acceleration due to advection
            !===========================================================
            adv_x = vel(i,1)*dvxdx + vel(i,2)*dvxdy + vel(i,3)*dvxdz
            adv_y = vel(i,1)*dvydx + vel(i,2)*dvydy + vel(i,3)*dvydz
            adv_z = vel(i,1)*dvzdx + vel(i,2)*dvzdy + vel(i,3)*dvzdz
        end if

        !===========================================================
        ! computing acceleration due to viscosity
        !===========================================================
        if (totphases == 1) then
            nu = viscosity(1)
        else
            nu = phi(i)*viscosity(2) + (1.0-phi(i))*viscosity(1)
        end if

        dif_x = nu*nabla2_vx
        dif_y = nu*nabla2_vy
        dif_z = nu*nabla2_vz
        !===========================================================
        ! computing intrmediate velocities
        !===========================================================
        vel_inter(i,1) = vel(i,1) + dt*(-adv_x + dif_x + F(i,1))
        vel_inter(i,2) = vel(i,2) + dt*(-adv_y + dif_y + F(i,2))
        vel_inter(i,3) = vel(i,3) + dt*(-adv_z + dif_z + F(i,3))
    end do
    !$omp end parallel do

!-------------------------------------------------------------------------
    ! Computing curvature for two-phase flows with surface tension
    if (segma > 0.0) then
        call get_curv(phi, pnts, krnls, totpnts, curv)
        !make sure intrpolants have xy, yx, xz
        
        !$omp parallel do private(i)
        do i = 1,totpnts
            dphidx = 0.0
            dphidy = 0.0
            dphidz = 0.0
            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                phi_diff = phi(nbr) - phi(i)
                dphidx = dphidx + krnls(i)%nx(j) * phi_diff
                dphidy = dphidy + krnls(i)%ny(j) * phi_diff
                dphidz = dphidz + krnls(i)%nz(j) * phi_diff
            end do
            ! computing acceleration due to surface tension
            abs_dphi = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz)
            tension_x = segma * curv(i) * dphidx / (abs_dphi + eps)
            tension_y = segma * curv(i) * dphidy / (abs_dphi + eps)
            tension_z = segma * curv(i) * dphidz / (abs_dphi + eps)
            vel_inter(i, 1) = vel_inter(i, 1) + dt * tension_x
            vel_inter(i, 2) = vel_inter(i, 2) + dt * tension_y
            vel_inter(i, 3) = vel_inter(i, 3) + dt * tension_z
        end do
        !$omp end parallel do
    end if
!-------------------------------------------------------------------------

    !============================================================================
    ! Computing the pressure iteratively
    !============================================================================
    vmax = 0;
    !$omp parallel do private(i)
    do i = 1, totpnts
        v0 = sqrt(vel_inter(i,1) * vel_inter(i,1)  + vel_inter(i,2) * vel_inter(i,2)  + vel_inter(i,3) * vel_inter(i,3))
        if (v0 >= vmax) then
            vmax = v0
        end if
    end do
    !$omp end parallel do
    c = 10*vmax
    
    if (totphases == 1) then
        ro(:) = density(1)
    else
        ro = phi*density(2) + (1.0-1.0*phi)*density(1)
    end if
    
    do vt = 1, iter
    ! computing pressure using the artificial compressibility method
        !$omp parallel do private(i)
        do i = 1, totpnts
            trm = 0
            do s = 1,krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(s)
                Vab_x = vel_inter(i,1)-vel_inter(nbr,1)
                Vab_y = vel_inter(i,2)-vel_inter(nbr,2)
                Vab_z = vel_inter(i,3)-vel_inter(nbr,3)
                trm = trm + (Vab_x*krnls(i)%nx(s)+Vab_y*krnls(i)%ny(s)+Vab_z*krnls(i)%nz(s))
            end do
            P_fut(i)  = P(i) + P_prmtr*dt*ro(i)*trm
        end do
        !$omp end parallel do
        P = P_fut

        ! diffusive pressure trm acts as stabilizer
        !$omp parallel do private(i)
        do i = 1,totpnts
            nabla2_P = 0
            do j = 1,krnls(i)%totnbrs
                pntnum = krnls(i)%nbrs(j)
                nabla2_P = nabla2_P + (P(pntnum)-P(i))*krnls(i)%nabla2(j)
            end do
            P_fut(i) = P(i) + dt*(nabla2_P) ! + 2*(dvydx*dvxdy+dvydz*dvzdy+dvxdz*dvzdx-dvxdx*dvydy*dvzdz));
        end do
        !$omp end parallel do
        P=P_fut
    end do
    call set_var(P, bcs(4))

    !$omp parallel do private(i) !be careful...for some reason it screwed up the artificial viscocity trm
    do i = 1,totpnts
        A_pressure_x = 0
        A_pressure_y = 0
        A_pressure_z = 0
        art_visc_x = 0
        art_visc_y = 0
        art_visc_z = 0

        do s = 1, krnls(i)%totnbrs
            j = krnls(i)%nbrs(s)

            !===========================================================
            ! computing acceleration due to pressure
            !===========================================================
            A_pressure_x = A_pressure_x + (P(i)+P(j))*krnls(i)%nx(s)    !do not touch you idiot
            A_pressure_y = A_pressure_y + (P(i)+P(j))*krnls(i)%ny(s)    !do not touch you idiot
            A_pressure_z = A_pressure_z + (P(i)+P(j))*krnls(i)%nz(s)    !do not touch you idiot

            !===========================================================
            ! computing acceleration due to artificial viscosity
            !===========================================================
            xij = pnts(i,1)-pnts(j,1)
            yij = pnts(i,2)-pnts(j,2)
            zij = pnts(i,3)-pnts(j,3)
            vxij = vel_inter(i,1)-vel_inter(j,1)
            vyij = vel_inter(i,2)-vel_inter(j,2)
            vzij = vel_inter(i,3)-vel_inter(j,3)
            uij_rij = xij*vxij + yij*vyij + zij*vzij

            alpha = av_prmtr
            d0 = sqrt(xij * xij + yij * yij + zij * zij);
            h0 = 2 * d0;
            trm0 = ((alpha * c * h0) / ro(i)) / (d0 * d0 + h0 * 1e-5)
            if (lagrangian .eqv. .true.) then
                trm0 = ((alpha * c * h0))/(d0 * d0 + h0 * 1e-5)
            end if
            trm_f = -uij_rij*trm0
            pii0 = trm_f
            if (trm_f < 0) then
                pii0 = 0
            end if
            art_visc_x = art_visc_x + pii0 * krnls(i)%nx(s)
            art_visc_y = art_visc_y + pii0 * krnls(i)%ny(s)
            art_visc_z = art_visc_z + pii0 * krnls(i)%nz(s)
        end do

        !===========================================================
        ! computing final velocity
        !===========================================================
        vel(i,1) = vel_inter(i,1) + dt * (-(1 / (ro(i) + eps)) * A_pressure_x - art_visc_x)
        vel(i,2) = vel_inter(i,2) + dt * (-(1 / (ro(i) + eps)) * A_pressure_y - art_visc_y)
        vel(i,3) = vel_inter(i,3) + dt * (-(1 / (ro(i) + eps)) * A_pressure_z - art_visc_z)
    end do
    !$omp end parallel do

    !============================================================================
    ! Applying boundary conditions
    !============================================================================
    call set_var(vel(:,1), bcs(1))
    call set_var(vel(:,2), bcs(2))
    call set_var(vel(:,3), bcs(3))

    v = sqrt(vel(:,1)*vel(:,1) + vel(:,2)*vel(:,2) + vel(:,3)*vel(:,3))

end subroutine slv_ns

end module ns
