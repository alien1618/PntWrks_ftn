module intrf
implicit none
contains

pure subroutine add_sphr(c, a, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit sphere to phi
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c
    real(8), intent(in)                    ::  a
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i
!-----------------------------------------------------------
    do i = 1, totpnts
        phi(i) = min(phi(i), ((pnts(i, 1) - c(1))**2) + ((pnts(i, 2) - c(2))**2) &
                     + ((pnts(i, 3) - c(3))**2) - (a**2));
    end do
end subroutine add_sphr

pure subroutine sub_sphr(c, a, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit sphere to phi
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c
    real(8), intent(in)                    ::  a
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i
!-----------------------------------------------------------
    do i = 1, totpnts
        phi(i) = max(phi(i), ((pnts(i, 1) - c(1))**2) + ((pnts(i, 2) - c(2))**2) &
                     + ((pnts(i, 3) - c(3))**2) - (a**2));
    end do
end subroutine sub_sphr

pure subroutine add_cyl_z(c, a, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit cylinder to phi along the z axis
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c
    real(8), intent(in)                    ::  a
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i   
!-----------------------------------------------------------
    do i = 1, totpnts
        phi(i) = min(phi(i), ((pnts(i, 1) - c(1))**2) + ((pnts(i, 2) - c(2))**2) - (a**2))
    end do
end subroutine add_cyl_z

pure subroutine sub_cyl_z(c, a, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit cylinder to phi along the z axis
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c
    real(8), intent(in)                    ::  a
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i   
!-----------------------------------------------------------
    do i = 1, totpnts
        phi(i) = max(phi(i), ((pnts(i, 1) - c(1))**2) + ((pnts(i, 2) - c(2))**2) - (a**2))
    end do
end subroutine sub_cyl_z

pure subroutine add_cube(c, l, pnts, totpnts, phi)
!-----------------------------------------------------------
!  add an implicit cube to phi
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c, l
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i  
    real(8)                                ::  x, y, z, width, height, depth 
!-----------------------------------------------------------

    x = c(1) + l(1) / 2.0
    y = c(2) + l(2) / 2.0
    z = c(3) + l(3) / 2.0
    width = l(1) / 2.0
    height = l(2) / 2.0
    depth = 1.0
    if (l(3) > 0.0) depth = l(3) / 2.0

    do i = 1, totpnts
        phi(i) = min(phi(i), (((pnts(i, 1) - x) / width)**100) + (((pnts(i, 2) - y) / height)**100) &
                    + (((pnts(i,3) - z) / depth)**6) - 1.0)
    end do
end subroutine add_cube

pure subroutine sub_cube(c, l, pnts, totpnts, phi)
!-----------------------------------------------------------
!  subtract an implicit cube from phi
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  pnts
    real(8), dimension(3), intent(in)      ::  c, l
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i  
    real(8)                                ::  x, y, z, width, height, depth
!-----------------------------------------------------------

    x = c(1) + l(1) / 2.0
    y = c(2) + l(2) / 2.0
    z = c(3) + l(3) / 2.0

    width = l(1) / 2.0
    height = l(2) / 2.0
    depth = 1.0
    if (l(3) > 0.0) depth = l(3) / 2.0

    do i = 1, totpnts
        phi(i) = max(phi(i), -((((pnts(i, 1) - x) / width)**100) + (((pnts(i, 2)-y) / height)**100) & 
                    + (((pnts(i, 3) - z) / depth)**6) - 1.0))
    end do
end subroutine sub_cube

pure subroutine get_vof(phi, totpnts, w)
!-----------------------------------------------------------
! convert a signed distance function to a heaveside function
! from 0 to 1
!-----------------------------------------------------------
    real(8), intent(in)                    ::  w
    integer, intent(in)                    ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)   ::  phi
!-----------------------------------------------------------
    integer                                ::  i
!-----------------------------------------------------------

    do i = 1, totpnts
        phi(i) = 0.5 * (1 - tanh(phi(i) / w))
        if (phi(i) > 1) phi(i) = 1
        if (phi(i) < 0) phi(i) = 0
    end do
end subroutine get_vof

pure subroutine set_vrtx_vel(vel, pnts, totpnts)
!--------------------------------------------------------------------
! assign vrtx flow field distribution
!--------------------------------------------------------------------
    real(8), dimension(:,:), intent(in)         ::  pnts
    integer, intent(in)                         ::  totpnts
!--------------------------------------------------------------------
    real(8), dimension(totpnts,3), intent(out)  ::  vel
!--------------------------------------------------------------------
    integer                                     ::  i
!--------------------------------------------------------------------

    do i = 1, totpnts
        vel(i, 1) = -(sin(3.14 * pnts(i, 1))**2) * sin(2 * 3.14 * pnts(i, 2))
        vel(i, 2) =  (sin(3.14 * pnts(i, 2))**2) * sin(2 * 3.14 * pnts(i, 1))
    end do
end subroutine set_vrtx_vel

pure subroutine set_shearing_vel(vel, pnts, totpnts)
!-----------------------------------------------------------
! assign shearing flow field distribution
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)                 ::  pnts
    integer, intent(in)                                 ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(totpnts,3), intent(out)   ::  vel
!-----------------------------------------------------------
    integer                                             ::  i
!-----------------------------------------------------------

    do i = 1, totpnts
        vel(i, 1) = sin(4 * 3.14 * (pnts(i, 1) + 0.5)) * sin(4 * 3.14 * (pnts(i, 2) + 0.5))
        vel(i, 2) = cos(4 * 3.14 * (pnts(i, 1) + 0.5)) * cos(4 * 3.14 * (pnts(i, 2) + 0.5))
    end do
end subroutine set_shearing_vel

pure subroutine set_rot_vel(vel, pnts, totpnts, u, v)
!-----------------------------------------------------------
! assign rotating flow field distribution
!-----------------------------------------------------------
    real(8), intent(in)                         :: u, v
    real(8), dimension(:,:), intent(in)         ::  pnts
    integer, intent(in)                         ::  totpnts
!-----------------------------------------------------------
    real(8), dimension(totpnts,3), intent(out)  ::  vel
!-----------------------------------------------------------
    integer                                     ::  i
!-----------------------------------------------------------
    do i = 1, totpnts
        vel(i, 1) = u * pnts(i, 2)
        vel(i, 2) = -v * pnts(i, 1)
    end do
end subroutine set_rot_vel

pure subroutine set_branch_vn(v_ext, pnts, totpnts, n, b)
!-----------------------------------------------------------
! assign 4 or 6 branch normal vector flow field
!-----------------------------------------------------------
    real(8), intent(in)                      :: b, n
    real(8), dimension(:,:), intent(in)      :: pnts
    integer, intent(in)                      :: totpnts
!-----------------------------------------------------------
    real(8), dimension(totpnts), intent(out) :: v_ext
!-----------------------------------------------------------
    integer                                  :: i
    real(8)                                  :: theta
!-----------------------------------------------------------
    do i = 1, totpnts
        theta = atan(pnts(i, 2) / pnts(i, 1))
        V_ext(i) = b * cos(n * theta)
    end do
end subroutine set_branch_vn

subroutine vof(phi_now, vel, vn, krnls, totpnts, dt, itrs)
!-----------------------------------------------------------
! solve volume of fluid (or level set) equation
!-----------------------------------------------------------
    use krnl_struct
    use omp_lib
!-----------------------------------------------------------    
    integer, intent(in)                     :: totpnts, itrs
    real(8), dimension(:,:), intent(in)     :: vel
    real(8), dimension(:),intent(in)        :: vn
    type (kernel),dimension(:), intent(in)  :: krnls
    real(8), intent(in)                     :: dt
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)    :: phi_now
!-----------------------------------------------------------
    integer                                 :: it, cntr = 0, i,  s
    real(8), dimension(:), allocatable      :: phi_itr 
    real(8)                                 :: dphi, dphidx, dphidy, dphidz, absdphi, adv
    real(8)                                 :: res, sum_old
!-----------------------------------------------------------
    allocate(phi_itr(totpnts))
    phi_itr = phi_now
    do it = 1, itrs
        cntr = cntr + 1
        sum_old = sum(phi_itr)
        !$omp parallel do private(i)
        do i=1,totpnts
            dphidx = 0
            dphidy = 0
            dphidz = 0
            do s = 1, krnls(i)%totnbrs
                dphi = phi_itr(krnls(i)%nbrs(s)) - phi_itr(i)
                dphidx = dphidx + krnls(i)%nx(s) * dphi
                dphidy = dphidy + krnls(i)%ny(s) * dphi
                dphidz = dphidz + krnls(i)%nz(s) * dphi
            end do
            absdphi = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz)
            adv = vel(i, 1) * dphidx + vel(i, 2) * dphidy + vel(i, 3) * dphidz
            phi_itr(i) = phi_now(i) - dt * (adv + vn(i) * absdphi)
        end do
        !$omp end parallel do
        res = abs(sum(phi_itr) - sum_old)
        !write(*,*) "iter = ", it, " , res = ", res
    end do
    phi_now = phi_itr
end subroutine vof

subroutine vof_gfd(phi_now, vel, vn, krnls, totpnts, dt, itrs)
!-----------------------------------------------------------
! solve volume of fluid (or level set) equation
!-----------------------------------------------------------
    use krnl_struct
    use omp_lib
!-----------------------------------------------------------
    integer, intent(in)                     ::  totpnts, itrs
    real(8), dimension(:,:), intent(in)     ::  vel
    real(8), dimension(:),intent(in)        ::  vn
    type (kernel),dimension(:), intent(in)  ::  krnls
    real(8), intent(in)                     ::  dt
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)    ::  phi_now
!-----------------------------------------------------------
    integer                                 ::  it, cntr = 0, i
    real(8), dimension(totpnts)             ::  phi_itr 
    real(8)                                 ::  dphidx, dphidy, dphidz, absdphi, adv
    real(8)                                 ::  res, sum_old
    real(8), dimension(9)                   ::  b, grad
!-----------------------------------------------------------

    phi_itr = phi_now
    do it = 1, itrs
        cntr = cntr + 1
        sum_old = sum(phi_itr)
        ! !$omp parallel do private(i)
        do i=1, totpnts
            b = matmul(krnls(i)%gfd_trm, (phi_itr(krnls(i)%nbrs) - phi_itr(i)))
            grad = matmul(krnls(i)%gfd_inv_a, b)
            dphidx = grad(1)
            dphidy = grad(2)
            dphidz = grad(3)
            absdphi = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz)
            adv = vel(i, 1) * dphidx + vel(i, 2) * dphidy + vel(i, 3) * dphidz
            phi_itr(i) = phi_now(i) - dt * (adv + vn(i) * absdphi)
        end do
        ! !$omp end parallel do
        res = abs(sum(phi_itr) - sum_old)
        !write(*,*) "iter = ", it, " , res = ", res
    end do
    phi_now = phi_itr
end subroutine vof_gfd

subroutine vof2(phi_now, vel, vn, krnls, totpnts, dt, itrs)
!-----------------------------------------------------------
! solves the conservative form of the volume of fluid equation
!-----------------------------------------------------------
    use krnl_struct
    use omp_lib
!-----------------------------------------------------------
    integer, intent(in)                     ::  totpnts, itrs
    real(8), dimension(:,:), intent(in)     ::  vel
    real(8), dimension(:), intent(in)       ::  vn
    type (kernel),dimension(:), intent(in)  ::  krnls
    real(8), intent(in)                     ::  dt
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)    ::  phi_now
!-----------------------------------------------------------
    integer                                 ::  it, nbr, i,  s
    real(8), dimension(totpnts)             ::  phi_itr
    real(8)                                 ::  dphi, dphidx, dphidy, dphidz, absdphi, adv
    real(8)                                 ::  dphidx1, dphidy1, dphidz1
!-----------------------------------------------------------

    phi_itr = phi_now
    do it = 1,itrs
        !$omp parallel do private(i)
        do i=1,totpnts
            dphidx = 0
            dphidy = 0
            dphidz = 0
            dphidx1 = 0
            dphidy1 = 0
            dphidz1 = 0
            do s = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(s)
                dphidx = dphidx + krnls(i)%nx(s) * (phi_itr(nbr) * vel(nbr, 1) - phi_itr(i) * vel(i, 1))
                dphidy = dphidy + krnls(i)%ny(s) * (phi_itr(nbr) * vel(nbr, 2) - phi_itr(i) * vel(i, 2))
                dphidz = dphidz + krnls(i)%nz(s) * (phi_itr(nbr) * vel(nbr, 3) - phi_itr(i) * vel(i, 3))
                dphi = phi_itr(nbr) - phi_itr(i)
                dphidx1 = dphidx1 + krnls(i)%nx(s) * dphi
                dphidy1 = dphidy1 + krnls(i)%ny(s) * dphi
                dphidz1 = dphidz1 + krnls(i)%nz(s) * dphi
            end do
            absdphi = sqrt(dphidx1 * dphidx1 + dphidy1 * dphidy1 + dphidz1 * dphidz1)
            adv = dphidx + dphidy + dphidz
            phi_itr(i) = phi_now(i) - dt * (adv + vn(i) * absdphi)
        end do
        !$omp end parallel do
    end do
    phi_now = phi_itr
end subroutine vof2

subroutine shrpn_intrf(phi, krnls, totpnts, vof, beta, eps, dx, dt, nt)
    !-----------------------------------------------------------
    ! iteratively sharpen a diffuse interface
    !-----------------------------------------------------------
    use krnl_struct
    real(8), intent(IN)                        ::  beta, eps, dx, dt
    integer, intent(IN)                     ::  totpnts, nt
    logical, intent(IN)                     ::  vof
    type(kernel), dimension(:), intent(IN)  ::  krnls
    real(8), dimension(:), intent(INOUT)       ::  phi
    integer                                 ::  nbr, k, i, j
    real(8)                                    ::  dxdx_inv, abs_phi, dphidx, dphidy, dphidz
    real(8)                                    ::  nabla_phi, dtrmxdx, dtrmydy, dtrmzdz
    real(8)                                    ::  nx, ny, nz
    real(8), dimension(:), allocatable         ::  phi_fut, trm_x, trm_y, trm_z

    allocate(phi_fut(totpnts))
    allocate(trm_x(totpnts))
    allocate(trm_y(totpnts))
    allocate(trm_z(totpnts))
    dxdx_inv = 1.0/(dx*dx)
    phi_fut = phi
    do k = 1,nt
        !$omp parallel do private(i)
        do i = 1,totpnts
            dphidx = 0
            dphidy = 0
            dphidz = 0
            do j = 1,krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                dphidx = dphidx + krnls(i)%nx(j)*(phi(nbr)-phi(i))
                dphidy = dphidy + krnls(i)%ny(j)*(phi(nbr)-phi(i))
                dphidz = dphidz + krnls(i)%nz(j)*(phi(nbr)-phi(i))
            end do
            abs_phi = (sqrt(dphidx*dphidx+dphidy*dphidy+dphidz*dphidz)+1e-6)
            nx= dphidx/abs_phi
            ny= dphidy/abs_phi
            nz= dphidz/abs_phi
            if (vof .eqv. .true.) then
                !for phi varying between 0 to 1
                trm_x(i) = phi(i)*(1-phi(i))*nx
                trm_y(i) = phi(i)*(1-phi(i))*ny
                trm_z(i) = phi(i)*(1-phi(i))*nz
            else
                !for phi varying between -1 to 1
                trm_x(i) = 0.25*(phi(i)+1)*(1-phi(i))*nx
                trm_y(i) = 0.25*(phi(i)+1)*(1-phi(i))*ny
                trm_z(i) = 0.25*(phi(i)+1)*(1-phi(i))*nz
            end if
        end do
        !$omp end parallel do

        !$omp parallel do private(i)
        do i = 1, totpnts
            dtrmxdx = 0
            dtrmydy = 0
            dtrmzdz = 0
            nabla_phi = 0
            do j = 1,krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                dtrmxdx = dtrmxdx + krnls(i)%nx(j)*(trm_x(nbr)-trm_x(i))
                dtrmydy = dtrmydy + krnls(i)%ny(j)*(trm_y(nbr)-trm_y(i))
                dtrmzdz = dtrmzdz + krnls(i)%nz(j)*(trm_z(nbr)-trm_z(i))
                nabla_phi = nabla_phi + dxdx_inv*krnls(i)%N(j)*(phi(nbr)-phi(i))
            end do
            phi_fut(i) = phi(i) + dt*beta*(eps*nabla_phi - (dtrmxdx+dtrmydy+dtrmzdz))
        end do
        !$omp end parallel do
        phi = phi_fut;
    end do

end subroutine shrpn_intrf

subroutine get_curv(phi, pnts, krnls, totpnts, curvature)
!-----------------------------------------------------------
! compute the curvature of a diffuse intrf
!-----------------------------------------------------------
    use krnl_struct
    use omp_lib
!-----------------------------------------------------------
    integer, intent(in)                      :: totpnts
    type(kernel), dimension(:), intent(in)   :: krnls
    real(8), dimension(:,:), intent(in)      :: pnts
    real(8), dimension(:), intent(in)        :: phi
!-----------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: curvature
!-----------------------------------------------------------
    real(8) :: dphidx, dphidy, dphidz
    real(8) :: dphidxx, dphidyy, dphidzz
    real(8) :: dphidxy, dphidxz, dphidyz
    real(8) :: absgradphi, dphi
    real(8) :: eps = 1e-6, t1, t2, t3, t4, t5, t6, t7, t8, t9
    real(8) :: d, dx, dy, dz, prmtr = 1.0
    real(8) :: nijx, nijy, nijz
    real(8) :: nxx, nyy, nzz, nxy, nxz, nyz
    integer :: i, j, nbr
!-----------------------------------------------------------
    allocate(curvature(totpnts))
    
    !$omp parallel do private(i)
    do i = 1, totpnts
        dphidx = 0
        dphidy = 0
        dphidz = 0
        dphidxx = 0
        dphidyy = 0
        dphidzz = 0
        dphidxy = 0
        dphidxz = 0
        dphidyz = 0
        do j = 1, krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            dphi = phi(nbr) - phi(i)
            dphidx = dphidx + dphi * krnls(i)%nx(j)
            dphidy = dphidy + dphi * krnls(i)%ny(j)
            dphidz = dphidz + dphi * krnls(i)%nz(j)

            dx = pnts(i, 1) - pnts(nbr, 1)
            dy = pnts(i, 2) - pnts(nbr, 2)
            dz = pnts(i, 3) - pnts(nbr, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz) + eps
            nijx = -dx / d
            nijy = -dy / d
            nijz = -dz / d
            nxx = prmtr * (nijx / d) * krnls(i)%nx(j)
            nyy = prmtr * (nijy / d) * krnls(i)%ny(j)
            nzz = prmtr * (nijz / d) * krnls(i)%nz(j)
            nxy = prmtr * (nijy / d) * krnls(i)%nx(j)
            nxz = prmtr * (nijz / d) * krnls(i)%nx(j)
            nyz = prmtr * (nijz / d) * krnls(i)%ny(j)

            dphidxx = dphidxx + dphi * nxx
            dphidyy = dphidyy + dphi * nyy
            dphidzz = dphidzz + dphi * nzz
            dphidxy = dphidxy + dphi * nxy
            dphidxz = dphidxz + dphi * nxz
            dphidyz = dphidyz + dphi * nyz

            !dphidxx = dphidxx + dphi * krnls(i)%nxx(j)
            !dphidyy = dphidyy + dphi * krnls(i)%nyy(j)
            !dphidzz = dphidzz + dphi * krnls(i)%nzz(j)
            !dphidxy = dphidxy + dphi * krnls(i)%nxy(j)
            !dphidxz = dphidxz + dphi * krnls(i)%nxz(j)
            !dphidyz = dphidyz + dphi * krnls(i)%nyz(j)

        end do
        absgradphi = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz) + eps
        t1 = dphidx * dphidx * dphidyy
        t2 = dphidx * dphidy * dphidxy
        t3 = dphidy * dphidy * dphidxx
        t4 = dphidx * dphidx * dphidzz
        t5 = dphidx * dphidz * dphidxz
        t6 = dphidz * dphidz * dphidxx
        t7 = dphidy * dphidy * dphidzz
        t8 = dphidy * dphidz * dphidyz
        t9 = dphidz * dphidz * dphidyy
        curvature(i) = (t1 - 2 * t2 + t3 + t4 - 2 * t5 + t6 + t7 - 2 * t8 + t9) / (absgradphi**3)
    end do
    !$omp end parallel do
end subroutine get_curv

subroutine get_curv_gfd(phi, krnls, totpnts, curvature)
!-----------------------------------------------------------
! compute the curvature of a diffuse intrf
!-----------------------------------------------------------
    use krnl_struct
    use omp_lib
!-----------------------------------------------------------
    integer, intent(in)                      :: totpnts
    type(kernel), dimension(:), intent(in)   :: krnls
    real(8), dimension(:), intent(in)        :: phi
!-----------------------------------------------------------
    real(8), dimension(totpnts), intent(out) :: curvature
!-----------------------------------------------------------
    real(8)                  :: dphidx, dphidy, dphidz
    real(8)                  :: dphidxx, dphidyy, dphidzz
    real(8)                  :: dphidxy, dphidxz, dphidyz
    real(8)                  :: absgradphi
    real(8)                  :: eps = 1e-6, t1, t2, t3, t4, t5, t6, t7, t8, t9
    integer                  :: i
    real(8), dimension(9)    :: b, grad
!-----------------------------------------------------------
    do i = 1, totpnts
        b = matmul(krnls(i)%gfd_trm,(phi(krnls(i)%nbrs)-phi(i)))
        grad = matmul(krnls(i)%gfd_inv_a,b)
        dphidx = grad(1)
        dphidy = grad(2)
        dphidz = grad(3)
        dphidxx = grad(4)
        dphidyy = grad(5)
        dphidzz = grad(6)
        dphidxy = grad(7)
        dphidxz = grad(8)
        dphidyz = grad(9)

        absgradphi = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz) + eps
        t1 = dphidx * dphidx * dphidyy
        t2 = dphidx * dphidy * dphidxy
        t3 = dphidy * dphidy * dphidxx
        t4 = dphidx * dphidx * dphidzz
        t5 = dphidx * dphidz * dphidxz
        t6 = dphidz * dphidz * dphidxx
        t7 = dphidy * dphidy * dphidzz
        t8 = dphidy * dphidz * dphidyz
        t9 = dphidz * dphidz * dphidyy
        curvature(i) = (t1 - 2 * t2 + t3 + t4 - 2 * t5 + t6 + t7 - 2 * t8 + t9) / (absgradphi**3)
    end do

end subroutine get_curv_gfd

subroutine ch(phi_now, vel, vn, krnls, totpnts, w, m0, dt, itrs)
!-----------------------------------------------------------
! solve cahn-hilliard equation
! phi needs to be 0 to 1
!-----------------------------------------------------------
    use krnl_struct
    use omp_lib
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)    ::  vel
    real(8), dimension(:), intent(in)      ::  vn
    type(kernel), dimension(:), intent(in) ::  krnls
    integer, intent(in)                    ::  totpnts, itrs
    real(8), intent(in)                    ::  w, m0, dt
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout) ::  phi_now
!-----------------------------------------------------------
    integer                              ::  i, j, it, nbr
    real(8), dimension(totpnts)          ::  mu, phi_itr
    real(8)                              ::  dphi, absdphi, adv, dfdcon
    real(8)                              ::  dphidx, dphidy, dphidz, nabla2_phi, nabla2_mu
!-----------------------------------------------------------
    phi_itr = phi_now
    do it = 1, itrs
        !$omp parallel do private(i)
        do i = 1, totpnts
            nabla2_phi = 0
            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                nabla2_phi = nabla2_phi + (phi_itr(nbr) - phi_itr(i)) * krnls(i)%nabla2(j)
            end do
            dfdcon = 2.0 * phi_itr(i) * ((1 - phi_itr(i))**2) - 2.0 * (phi_itr(i)**2) * (1.0 - phi_itr(i))
            mu(i) = dfdcon - w * (nabla2_phi)
        end do
        !$omp end parallel do

        !$omp parallel do private(i)
        do i = 1, totpnts
            dphidx = 0
            dphidy = 0
            dphidz = 0
            nabla2_mu = 0
            do j = 1, krnls(i)%totnbrs   
                nbr = krnls(i)%nbrs(j)
                dphi = phi_itr(nbr) - phi_itr(i)
                dphidx = dphidx + dphi * krnls(i)%nx(j)
                dphidy = dphidy + dphi * krnls(i)%ny(j)
                dphidz = dphidz + dphi * krnls(i)%nz(j)
                nabla2_mu = nabla2_mu + (mu(nbr) - mu(i)) * krnls(i)%nabla2(j)
            end do
            absdphi = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz)
            adv = dphidx * vel(i,1) + dphidy * vel(i,2) + dphidz * vel(i,3)
            phi_itr(i) = phi_now(i) + dt * (m0 * nabla2_mu - adv + absdphi * vn(i))
            if (phi_itr(i) >= 0.999) phi_itr(i)= 1
            if (phi_itr(i) <= 0.0001) phi_itr(i) = 0
        end do
        !$omp end parallel do
    end do
    phi_now = phi_itr
end subroutine ch

subroutine ch_gfd(phi_now, vel, vn, krnls, totpnts, w, m0, dt, itrs)
!-----------------------------------------------------------
! solve cahn-hilliard equation
! phi needs to be 0 to 1
!-----------------------------------------------------------
    use krnl_struct
    use omp_lib
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)        :: vel
    real(8), dimension(:), intent(in)          :: vn
    type(kernel), dimension(:), intent(in)     :: krnls
    integer, intent(in)                        :: totpnts, itrs
    real(8), intent(in)                        :: w, m0, dt
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)       :: phi_now
!-----------------------------------------------------------
    integer                                    :: i, it
    real(8), dimension(totpnts)                :: mu, phi_itr
    real(8)                                    :: absdphi, adv, dfdcon
    real(8)                                    :: dphidx, dphidy, dphidz, nabla2_phi, nabla2_mu
    real(8), dimension(9)                      :: b, grad
!-----------------------------------------------------------

    phi_itr = phi_now
    do it = 1, itrs
        !$omp parallel do private(i)
        do i = 1, totpnts
            b = matmul(krnls(i)%gfd_trm, (phi_itr(krnls(i)%nbrs) - phi_itr(i)))
            grad = matmul(krnls(i)%gfd_inv_a, b)
            dphidx = grad(1)
            dphidy = grad(2)
            dphidz = grad(3)
            nabla2_phi = sum(grad(4 : 6))
            dfdcon =2.0 * phi_itr(i) * ((1 - phi_itr(i))**2) - 2.0 * (phi_itr(i)**2) * (1.0 - phi_itr(i))
            mu(i) = dfdcon - w * (nabla2_phi)
        end do
        !$omp end parallel do

        !$omp parallel do private(i)
        do i = 1,totpnts
            b = matmul(krnls(i)%gfd_trm, (mu(krnls(i)%nbrs) - mu(i)))
            grad = matmul(krnls(i)%gfd_inv_a, b)
            nabla2_mu = grad(4) + grad(5) + grad(6)
            absdphi = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz)
            adv = dphidx * vel(i,1) + dphidy * vel(i,2) + dphidz*vel(i,3)
            phi_itr(i)  = phi_now(i) + dt * (m0 * nabla2_mu - adv + absdphi * vn(i))
            if (phi_itr(i) >= 0.999) phi_itr(i)= 1.0
            if (phi_itr(i) <= 0.0001) phi_itr(i) = 0.0
        end do
        !$omp end parallel do
    end do
    phi_now = phi_itr
end subroutine ch_gfd

subroutine ch2(phi_now, vel, vn, krnls, totpnts, w, m0, dt, itrs)
!-----------------------------------------------------------
! another version to solve cahn-hilliard equation
! phi needs to be -1 to 1
!-----------------------------------------------------------
    use krnl_struct
    use omp_lib
!-----------------------------------------------------------
    real(8), dimension(:,:), intent(in)        ::  vel
    real(8), dimension(:), intent(in)          ::  vn
    integer, intent(in)                        ::  totpnts, itrs
    real(8), intent(in)                        ::  w, m0, dt
    type(kernel), dimension(:), intent(in)     ::  krnls
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)       ::  phi_now
!-----------------------------------------------------------
    integer                                    ::  i, j, it, nbr
    real(8), dimension(totpnts)                ::  mu, phi_itr
    real(8)                                    ::  dphi, absdphi, adv
    real(8)                                    ::  dphidx, dphidy, dphidz, nabla2_phi, nabla2_mu
!-----------------------------------------------------------

    phi_itr = phi_now  
    do it = 1, itrs
        !$omp parallel do private(i)
        do i = 1, totpnts
            nabla2_phi = 0
            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                nabla2_phi = nabla2_phi + (phi_itr(nbr) - phi_itr(i)) * krnls(i)%nabla2(j)
            end do
            mu(i) =  (phi_itr(i) * ((phi_itr(i) * phi_itr(i)) - 1) - w * w * nabla2_phi) !be careful the sign of diff
        end do
        !$omp end parallel do

        !$omp parallel do private(i)
        do i = 1, totpnts
            dphidx = 0
            dphidy = 0
            dphidz = 0
            nabla2_mu = 0
            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                dphi = phi_itr(nbr) - phi_itr(i)
                dphidx = dphidx + dphi * krnls(i)%nx(j)
                dphidy = dphidy + dphi * krnls(i)%ny(j)
                dphidz = dphidz + dphi * krnls(i)%nz(j)
                nabla2_mu = nabla2_mu + (mu(nbr) - mu(i)) * krnls(i)%nabla2(j)
            end do
            absdphi = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz)
            adv = dphidx * vel(i,1) + dphidy * vel(i,2) + dphidz * vel(i,3)
            phi_itr(i) = phi_now(i) + dt * (m0 * nabla2_mu - adv + absdphi * vn(i))
            if(phi_itr(i) > 1) phi_itr(i)= 1
            if(phi_itr(i) < -1) phi_itr(i) = -1
        end do
        !$omp end parallel do
    end do
    phi_now = phi_itr
end subroutine ch2

subroutine ac(phi_now, vel, vn, krnls, totpnts, w, m0, dt, itrs)
!-----------------------------------------------------------
! solve allen-cah equation using the conservative form
! of the advection term
!-----------------------------------------------------------
    use omp_lib
    use krnl_struct
!-----------------------------------------------------------
    real(8), intent(in)                        :: w, m0, dt
    integer, intent(in)                        :: itrs, totpnts
    real(8), dimension(:), intent(in)          :: vn
    real(8), dimension(:,:), intent(in)        :: vel
    type(kernel), dimension(:), intent(in)     :: krnls
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)       :: phi_now
!-----------------------------------------------------------
    real(8), dimension(totpnts)                :: phi_itr, trmx, trmy, trmz, absdphi
    integer                                    :: i, s, nbr, it
    real(8)                                    :: dphi, dphidx, dphidy, dphidz
    real(8)                                    :: dtrmdx, dtrmdy, dtrmdz, rhs, adv
    real(8)                                    :: dtrm2dx, dtrm2dy, dtrm2dz
!-----------------------------------------------------------

    phi_itr = phi_now
    do it = 1, itrs
        !$omp parallel do private(i)
        do i=1, totpnts
            dphidx = 0
            dphidy = 0
            dphidz = 0
            do s = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(s)
                dphi = phi_itr(nbr) - phi_itr(i)
                dphidx = dphidx + krnls(i)%nx(s) * dphi
                dphidy = dphidy + krnls(i)%ny(s) * dphi
                dphidz = dphidz + krnls(i)%nz(s) * dphi
            end do
            absdphi(i) = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz)
            trmx(i) = w * dphidx - phi_itr(i) * (1 - phi_itr(i)) * (dphidx / absdphi(i))
            trmy(i) = w * dphidy - phi_itr(i) * (1 - phi_itr(i)) * (dphidy / absdphi(i))
            trmz(i) = w * dphidz - phi_itr(i) * (1 - phi_itr(i)) * (dphidz / absdphi(i))
        end do
        !$omp end parallel do

        !$omp parallel do private(i)
        do i=1,totpnts
            dtrmdx = 0;
            dtrmdy = 0;
            dtrmdz = 0;

            dtrm2dx = 0;
            dtrm2dy = 0;
            dtrm2dz = 0;
            do s = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(s);
                dtrmdx = dtrmdx + (trmx(nbr) - trmx(i)) * krnls(i)%nx(s)
                dtrmdy = dtrmdy + (trmy(nbr) - trmy(i)) * krnls(i)%ny(s)
                dtrmdz = dtrmdz + (trmz(nbr) - trmz(i)) * krnls(i)%nz(s)

                dtrm2dx = dtrm2dx + ((phi_itr(nbr) * vel(nbr,1)) - (phi_itr(i) * vel(i,1))) * krnls(i)%nx(s)
                dtrm2dy = dtrm2dy + ((phi_itr(nbr) * vel(nbr,2)) - (phi_itr(i) * vel(i,2))) * krnls(i)%ny(s)
                dtrm2dz = dtrm2dz + ((phi_itr(nbr) * vel(nbr,3)) - (phi_itr(i) * vel(i,3))) * krnls(i)%nz(s)
            end do
            rhs = dtrmdx + dtrmdy + dtrmdz
            adv = dtrm2dx + dtrm2dy + dtrm2dz
            phi_itr(i) = phi_now(i) + dt * (m0 * rhs - adv + absdphi(i) * vn(i))
            if (phi_itr(i) < 0) phi_itr(i) = 0
            if (phi_itr(i) > 1) phi_itr(i) = 1
        end do
        !$omp end parallel do
    end do
    phi_now = phi_itr
end subroutine ac

subroutine ac2(phi_now, vel, vn, krnls, totpnts, m0, dt, itrs)
!-----------------------------------------------------------
! another version to solve allen-cahn equation
!-----------------------------------------------------------
    use krnl_struct
    use omp_lib
!-----------------------------------------------------------
    real(8), intent(in)                        :: m0, dt
    integer, intent(in)                        :: itrs, totpnts
    real(8), dimension(:), intent(in)          :: vn
    real(8), dimension(:,:), intent(in)        :: vel
    type(kernel), dimension(:), intent(in)     :: krnls
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)       :: phi_now
!-----------------------------------------------------------
    real(8), dimension(totpnts)                :: phi_itr
    integer                                    :: i, s, nbr, it
    real(8)                                    :: eps = 1e-6, phi2, dphi, dphidx, dphidy, dphidz, diff, adv
    real(8)                                    :: absdphi, sum_F_prime, sum_F, beta, k, nabla2_phi
!-----------------------------------------------------------
    phi_itr = phi_now
    do it = 1, itrs
        sum_F_prime = 0;
        sum_F = 0;
        do i = 1, totpnts
            phi2 = phi_now(i) * phi_now(i)
            sum_F_prime = sum_F_prime + phi_now(i) * (phi2 - 1)
            sum_F = sum_F + 0.25 * ((phi2 - 1)**2)
        end do
        beta = sum_F_prime / (sum_F + eps)

        !$omp parallel do private(i)
        do i = 1, totpnts
            dphidx = 0
            dphidy = 0
            dphidz = 0
            nabla2_phi = 0
            do s = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(s)
                dphi =  phi_itr(nbr) - phi_itr(i)
                dphidx = dphidx + dphi * krnls(i)%nx(s)
                dphidy = dphidy + dphi * krnls(i)%ny(s)
                dphidz = dphidz + dphi * krnls(i)%nz(s)
                nabla2_phi = nabla2_phi + dphi * krnls(i)%nabla2(s)
            end do
            absdphi = sqrt(dphidx * dphidx + dphidy * dphidy + dphidz * dphidz)
            adv = dphidx * vel(i,1) + dphidy * vel(i,2) + dphidz * vel(i,3)
            phi2 = phi_now(i) * phi_now(i)
            k = 2 * sqrt(0.25 * ((phi2 - 1)**2))
            diff =  m0 * nabla2_phi - phi_now(i) * (phi2 - 1)
            phi_itr(i) = phi_now(i) + dt*(diff - adv + absdphi * vn(i) + beta * k)

            if(phi_itr(i) > 1) phi_itr(i) = 1
            if(phi_itr(i) < -1) phi_itr(i) = -1
        end do
        !$omp end parallel do
    end do
    phi_now = phi_itr
end subroutine ac2


subroutine get_intrf_pnts(pnts, phi, krnls, totpnts, dxx, ip, totip)
    use krnl_struct
!----------------------------------------------------------------------    
    real(8), dimension(:,:), intent(in) :: pnts
    real(8), dimension(:), intent(in) :: phi
    real(8), intent(in) :: dxx
    integer, intent(in) :: totpnts
    type (kernel), dimension(:), intent(in)    :: krnls
!----------------------------------------------------------------------
    integer, intent(out) :: totip
    real(8), dimension(:,:), allocatable, intent(out) :: ip
!----------------------------------------------------------------------
    real(8), dimension(:), allocatable :: ipx, ipy, ipz, unq_x, unq_y, unq_z
    integer :: i, j, m, nbr, totunqpnts
    real(8) :: phi1, phi2, x1, x2, y1, y2, z1, z2, px, py, pz, product
    real(8) :: dx, dy, dz, d
    integer, dimension(:), allocatable :: indx
!----------------------------------------------------------------------
    !Function works for a VOF function (0 to 1)
    totip = 0
    allocate(ipx(0))
    allocate(ipy(0))
    allocate(ipz(0))
    do i = 1, totpnts
        phi1 = phi(i)
        do j = 1, krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            phi2 = phi(nbr)
            product = phi1 * phi2
            if (product <= 1e-6) then
                x1 = pnts(i, 1)
                y1 = pnts(i, 2)
                z1 = pnts(i, 3)
                x2 = pnts(nbr, 1)
                y2 = pnts(nbr, 2)
                z2 = pnts(nbr, 3)
                totip = totip + 1
                px = (phi1*x2-phi2*x1)/(phi1-phi2)
                py = (phi1*y2-phi2*y1)/(phi1-phi2)
                pz = (phi1*z2-phi2*z1)/(phi1-phi2)
                ipx = [ipx, px]
                ipy = [ipy, py]
                ipz = [ipz, pz]
            end if
        end do
    end do

    totunqpnts = 0
    allocate(unq_x(0))
    allocate(unq_y(0))
    allocate(unq_z(0))
    
    allocate(indx(totip))
    indx(:) = 0

    do i = 1, totip
        m = 0
        do j = 1, totip
            dx = ipx(i) - ipx(j)
            dy = ipy(i) - ipy(j)
            dz = ipz(i) - ipz(j)
            d = (dx*dx + dy*dy + dz*dz)**0.5
            if (indx(j) == 0 .and. d <= 0.5*dxx) then
                m = m + 1
                indx(j) = m
            end if
        end do
    end do
    do i = 1, totip
        if (indx(i) == 1) then
            totunqpnts = totunqpnts + 1
            unq_x = [unq_x, ipx(i)]
            unq_y = [unq_y, ipy(i)]
            unq_z = [unq_z, ipz(i)]
        end if
    end do
    
    totip = totunqpnts
    allocate(ip(totip,3))
    ip(:,1) = unq_x
    ip(:,2) = unq_y
    ip(:,3) = unq_z

end subroutine get_intrf_pnts

end module intrf