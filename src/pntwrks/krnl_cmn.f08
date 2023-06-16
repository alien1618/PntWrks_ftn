module krnl_cmn
implicit none

contains

pure subroutine get_distance(pnt1, pnt2, d)
    real(8), dimension(3), intent(in) :: pnt1, pnt2
    real(8), intent(out) :: d
    real(8) :: dx, dy, dz
    
    dx = pnt1(1) - pnt2(2)
    dy = pnt1(2) - pnt2(2)
    dz = pnt1(3) - pnt2(2)
    d = sqrt(dx * dx + dy * dy + dz * dz)
end subroutine get_distance 

subroutine get_nbrs_bf(h, pnts, totpnts, sctr, totsctr, krnls)
!----------------------------------------------------------------------
! searches for neightbours using brute force
!----------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!----------------------------------------------------------------------
    real(8), dimension(:,:),intent(in)               :: pnts, sctr
    integer, intent(in)                              :: totpnts, totsctr
    real(8), intent(in)                              :: h
!----------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout)        :: krnls
!----------------------------------------------------------------------
    real(8)                                          :: dx, dy, dz, d
    real(8)                                          :: tol = 1e-6
    integer                                          :: i, j
!----------------------------------------------------------------------
    !DONOT PARALLELIZE. RESULTED IN SLOW DOWNS BY A FACTOR OF 2 !!!
    do i = 1,totpnts
        krnls(i)%totnbrs = 0
        allocate(krnls(i)%nbrs(0))
        do j = 1, totsctr
            dx = pnts(i, 1) - sctr(j, 1)
            dy = pnts(i, 2) - sctr(j, 2)
            dz = pnts(i, 3) - sctr(j, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= h .and. d >= tol) then
                krnls(i)%nbrs = [krnls(i)%nbrs, j]
                krnls(i)%totnbrs = krnls(i)%totnbrs + 1
            end if
        end do
    end do
end subroutine get_nbrs_bf

subroutine get_nbrs_bg(pnts, pntbins, totpnts, bins, nx, nxny, totbins, dim, h, krnls)
!----------------------------------------------------------------------
! searches for neighbours using background grid
!----------------------------------------------------------------------
    use omp_lib
    use krnl_struct
!----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in)         ::  pnts
    integer, intent(in)                        ::  dim, totpnts, nx, nxny, totbins
    real(8), intent(in)                        ::  h
    type(kernel), dimension(:), intent(in)     ::  bins
    integer, dimension(:), intent(in)          ::  pntbins
!----------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout)  ::  krnls
!----------------------------------------------------------------------
    real(8)                                    ::  dx, dy, dz, d, tol = 1.0e-6
    integer                                    ::  i, j, c, nbr, tot
    integer, dimension(:), allocatable         :: cells
!----------------------------------------------------------------------

    tot = 9
    if (dim == 3) tot = 27
    allocate(cells(tot))

    !DONOT PARALLELIZE. RESULTED IN SLOW DOWNS BY A FACTOR OF 2 !!!
    do i = 1,totpnts
        krnls(i)%totnbrs = 0
        allocate(krnls(i)%nbrs(0))

        cells(1) = pntbins(i) !cntr
        cells(2) = cells(1) + 1 !right
        cells(3) = cells(1) - 1 !left
        cells(4) = cells(1) + nx !up
        cells(5) = cells(1) - nx !low
        cells(6) = cells(4) + 1  !up right
        cells(7) = cells(4) - 1  !up left
        cells(8) = cells(5) + 1  !low right
        cells(9) = cells(5) - 1  !low left
        if (dim == 3) then 
            cells(10) = pntbins(i) + nxny
            cells(11) = cells(1) + 1 + nxny
            cells(12) = cells(1) - 1 + nxny
            cells(13) = cells(1) + nx + nxny
            cells(14) = cells(1) - nx + nxny
            cells(15) = cells(4) + 1  + nxny
            cells(16) = cells(4) - 1  + nxny
            cells(17) = cells(5) + 1  + nxny
            cells(18) = cells(5) - 1  + nxny

            cells(19) = pntbins(i) - nxny
            cells(20) = cells(1) + 1 - nxny
            cells(21) = cells(1) - 1 - nxny
            cells(22) = cells(1) + nx - nxny
            cells(23) = cells(1) - nx - nxny
            cells(24) = cells(4) + 1  - nxny
            cells(25) = cells(4) - 1  - nxny
            cells(26) = cells(5) + 1  - nxny
            cells(27) = cells(5) - 1  - nxny
        end if
        do c = 1, tot
            if (cells(c) > 0 .and. cells(c) < totbins) then
                do j = 1, bins(cells(c))%totnbrs
                    nbr = bins(cells(c))%nbrs(j)
                    dx = pnts(i, 1) - pnts(nbr, 1)
                    dy = pnts(i, 2) - pnts(nbr, 2)
                    dz = pnts(i, 3) - pnts(nbr, 3)
                    d = sqrt(dx * dx + dy * dy + dz * dz)
                    if (d <= h .and. d >= tol) then
                        krnls(i)%nbrs = [krnls(i)%nbrs, nbr]
                        krnls(i)%totnbrs = krnls(i)%totnbrs + 1
                    end if
                end do
            end if
        end do
    end do
end subroutine get_nbrs_bg

subroutine get_nbrs_sp(h, phi, pnts, totpnts, krnls)
!----------------------------------------------------------------------
! searches for same phase neighbours
!----------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!----------------------------------------------------------------------
    real(8), dimension(:,:),intent(in)                      ::  pnts
    real(8), dimension(:), intent(in)                       ::  phi
    integer, intent(in)                                     ::  totpnts
    real(8), intent(in)                                     ::  h
!----------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout)               ::  krnls
!----------------------------------------------------------------------
    real(8)                                                 ::  dx, dy, dz, d, tol = 1e-6
    integer                                                 ::  i, j
!----------------------------------------------------------------------

    !DONOT PARALLELIZE. RESULTED IN SLOW DOWNS BY A FACTOR OF 2 !!!
    do i = 1,totpnts
        krnls(i)%totnbrs = 0
        allocate(krnls(i)%nbrs(0))
        do j = 1, totpnts
            dx = pnts(i, 1) - pnts(j, 1)
            dy = pnts(i, 2) - pnts(j, 2)
            dz = pnts(i, 3) - pnts(j, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= h .and. d >= tol) then
                if (phi(i) <= 0.5 .and. phi(j) <= 0.5) then
                    krnls(i)%nbrs = [krnls(i)%nbrs, j]
                    krnls(i)%totnbrs = krnls(i)%totnbrs + 1
                else if (phi(i) >= 0.5 .and. phi(j) >= 0.5) then
                    krnls(i)%nbrs = [krnls(i)%nbrs, j]
                    krnls(i)%totnbrs = krnls(i)%totnbrs + 1
                end if
            end if
        end do
    end do
end subroutine get_nbrs_sp

subroutine get_intrps_o2(pnts, totpnts, krnls)
!----------------------------------------------------------------------
! Compute laplacian using first order interpolants
!----------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in)                 :: totpnts
!----------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) :: krnls
!----------------------------------------------------------------------
    integer :: i, j, nbr, totnbrs
    real(8) :: eps = 1.0e-6, d, dx, dy, dz
    real(8) :: nijx, nijy, nijz, prmtr = 1.0, nxx, nyy, nzz
!----------------------------------------------------------------------

    !$omp parallel do private(i)
    do i = 1,totpnts
        totnbrs = krnls(i)%totnbrs
        allocate(krnls(i)%nabla2(totnbrs))
        do j = 1,totnbrs
            nbr = krnls(i)%nbrs(j)
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
            krnls(i)%nabla2(j) = nxx + nyy + nzz
        end do
    end do
    !$omp end parallel do
end subroutine get_intrps_o2

subroutine get_intrps_o22(ip, totip, pnts, krnls)
!----------------------------------------------------------------------
! Compute laplacian using first order interpolants
!----------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: ip, pnts
    integer, intent(in)                 :: totip
!----------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) :: krnls
!----------------------------------------------------------------------
    integer :: i, j, nbr, totnbrs
    real(8) :: eps = 1.0e-6, d, dx, dy, dz
    real(8) :: nijx, nijy, nijz, prmtr = 1.0, nxx, nyy, nzz
!----------------------------------------------------------------------

!    !$omp parallel do private(i)
    do i = 1,totip
        totnbrs = krnls(i)%totnbrs
       ! allocate(krnls(i)%nabla2(totnbrs))
        do j = 1,totnbrs
            nbr = krnls(i)%nbrs(j)
            dx = ip(i, 1) - pnts(nbr, 1)
            dy = ip(i, 2) - pnts(nbr, 2)
            dz = ip(i, 3) - pnts(nbr, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz) + eps
            nijx = -dx / d
            nijy = -dy / d
            nijz = -dz / d
            nxx = prmtr * (nijx / d) * krnls(i)%nx(j)
            nyy = prmtr * (nijy / d) * krnls(i)%ny(j)
            nzz = prmtr * (nijz / d) * krnls(i)%nz(j)
            krnls(i)%nabla2(j) = nxx + nyy + nzz
        end do
    end do
 !   !$omp end parallel do
end subroutine get_intrps_o22

subroutine get_intrps_o2_v2(pnts, totpnts, krnls)
!----------------------------------------------------------------------
! Compute first order gradients and laplacian using
! only zero order interpolants
!----------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in)                 :: totpnts
!----------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) :: krnls
!----------------------------------------------------------------------
    integer :: i, j, nbr, totnbrs
    real(8) :: d, xj_xi, yj_yi, zj_zi, prmtr = 1.0
!----------------------------------------------------------------------

    !$omp parallel do private(i)
    do i = 1, totpnts
        totnbrs = krnls(i)%totnbrs
        allocate(krnls(i)%nabla2(totnbrs))
        do j = 1, totnbrs
            nbr = krnls(i)%nbrs(j)
            xj_xi = pnts(nbr, 1) - pnts(i, 1)
            yj_yi = pnts(nbr, 2) - pnts(i, 2)
            zj_zi = pnts(nbr, 3) - pnts(i, 3)
            d = sqrt((xj_xi**2) + (yj_yi**2) + (zj_zi**2))

            krnls(i)%nx(j) = prmtr * xj_xi * krnls(i)%n(j)
            krnls(i)%ny(j) = prmtr * yj_yi * krnls(i)%n(j)
            krnls(i)%nz(j) = prmtr * zj_zi * krnls(i)%n(j)
            krnls(i)%nabla2(j) = prmtr * krnls(i)%n(j) / (d * d)
        end do
    end do
    !$omp end parallel do
end subroutine get_intrps_o2_v2

pure subroutine get_grad(ui, u_nbrs, krnl, ux, uy, uz, nabla2_u)          
!----------------------------------------------------------------------
! Compute the gradients of field variable u at point i using
! u at neighbouring points and kernel gradient interpolants
!----------------------------------------------------------------------
    use krnl_struct
    implicit none
!----------------------------------------------------------------------
    real(8), intent(in) :: ui
    type(kernel), intent(in) :: krnl
    real(8), dimension(:), intent(in) :: u_nbrs
!----------------------------------------------------------------------
    real(8), intent(out) :: ux, uy, uz, nabla2_u
!----------------------------------------------------------------------
    integer :: j
    real(8) :: du
!----------------------------------------------------------------------

    ux = 0
    uy = 0
    uz = 0
    nabla2_u = 0
    do j = 1, krnl%totnbrs
        du = u_nbrs(j) - ui
        ux = ux + du * krnl%nx(j)
        uy = uy + du * krnl%ny(j)
        uz = uz + du * krnl%nz(j)
        nabla2_u = nabla2_u + du * krnl%nabla2(j)
    end do
end subroutine get_grad

subroutine get_bin_bnds(pnts, dxx, dim, nx, nxy, totbins, bg_b, bg_l, bg_d)
!-------------------------------------------------------------------------
! Calculate bounds of neighbour search bins
!-------------------------------------------------------------------------
    use krnl_struct
    implicit none
!-------------------------------------------------------------------------   
    integer, intent(in) :: nx
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: dim
    real(8), intent(in) :: dxx
!-------------------------------------------------------------------------
    integer, intent(out) :: totbins, nxy
    real(8), dimension(3), intent(out) :: bg_b, bg_l, bg_d
!-------------------------------------------------------------------------
    real(8) :: dx, dy, dz, lx, ly, lz, minx, miny, minz, maxx, maxy, maxz
    integer :: ny, nz
!-------------------------------------------------------------------------
    minx = minval(pnts(:,1)) - 2 * dxx
    maxx = maxval(pnts(:,1)) + 2 * dxx
    miny = minval(pnts(:,2)) - 2 * dxx
    maxy = maxval(pnts(:,2)) + 2 * dxx
    minz = 0.0
    maxz = 0.0
    if (dim == 3) then
        minz = minval(pnts(:,3)) - 2 * dxx
        maxz = maxval(pnts(:,3)) + 2 * dxx
    end if
    lx = maxx - minx
    ly = maxy - miny
    lz = maxz - minz
    ny = int(ly * nx / lx)
    nz = int(lz * nx / lx)
    if (dim == 2) nz = 1
    dx = (lx-minx) / nx
    dy = (ly-miny) / ny
    dz = (lz-minz) / nz
    nxy = nx * ny
    totbins = nx * ny * nz
    bg_b = [minx, miny, minz]
    bg_l = [lx, ly, lz]
    bg_d = [dx, dy, dz]
end subroutine get_bin_bnds

subroutine get_bins(nx, nxy, totbins, bg_b, bg_d, pnts, totpnts, dim, bins, pntbins)
!------------------------------------------------------------------------------------------
! Calculate bounds of neighbour search bins
!------------------------------------------------------------------------------------------
    use krnl_struct
    implicit none
!------------------------------------------------------------------------------------------
    integer, intent(in) :: totbins, totpnts, dim, nx, nxy
    real(8), dimension(3), intent(in) :: bg_b, bg_d
    real(8), dimension(:,:), intent(in) :: pnts
!------------------------------------------------------------------------------------------
    integer, dimension(:), allocatable, intent(inout)      :: pntbins
    type(kernel), dimension(:), allocatable, intent(inout) :: bins
!------------------------------------------------------------------------------------------
    integer :: i, ix, iy, iz
!------------------------------------------------------------------------------------------
        do i = 1, totbins
            bins(i)%totnbrs = 0
            allocate(bins(i)%nbrs(0))
        end do
        do i = 1, totpnts
            ix = int((pnts(i,1) - bg_b(1)) / bg_d(1))
            iy = int((pnts(i,2) - bg_b(2)) / bg_d(2))
            if (dim == 3) then
                iz = int((pnts(i,3) - bg_b(3)) / bg_d(3))
                pntbins(i) = 1 + ix + nx * iy + nxy * iz
            else
                pntbins(i) = 1 + ix + nx * iy
            end if
            if (pntbins(i) > totbins) then 
                write(*,*) "ERROR: point is out of bounds"
                write(*,*) pnts(i,1), " ", pnts(i,2), " ", pnts(i,3)
                call exit(0)
            end if
            bins(pntbins(i))%nbrs = [bins(pntbins(i))%nbrs, i]
            bins(pntbins(i))%totnbrs = bins(pntbins(i))%totnbrs + 1
        end do
end subroutine get_bins

subroutine gen_bins(pnts, totpnts, dim, dx, nx, nxy, totbins, bins, pntbins)
!-------------------------------------------------------------------------
! Calculate bounds of neighbour search bins
!-------------------------------------------------------------------------
    use krnl_struct
    implicit none
!-------------------------------------------------------------------------   
    integer, intent(in) :: nx
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: totpnts, dim
    real(8), intent(in) :: dx
!-------------------------------------------------------------------------
    integer, intent(out) :: totbins, nxy
    integer, dimension(:), allocatable, intent(out)      :: pntbins
    type(kernel), dimension(:), allocatable, intent(out) :: bins
!-------------------------------------------------------------------------
    real(8), dimension(3) :: bg_b, bg_l, bg_d
!-------------------------------------------------------------------------

    call get_bin_bnds(pnts, dx, dim, nx, nxy, totbins, bg_b, bg_l, bg_d)
    allocate(bins(totbins))
    allocate(pntbins(totpnts)) 
    call get_bins(nx, nxy, totbins, bg_b, bg_d, pnts, totpnts, dim, bins, pntbins)
end subroutine gen_bins

subroutine get_bases(order, pnt, monomials, bases)
!-----------------------------------------------------------
! Compute bases polynomials and its derivatives to 
! second order
!-----------------------------------------------------------
    use krnl_struct
!-----------------------------------------------------------
    integer, intent(in)                :: order
    real(8), dimension(:), intent(in)  :: pnt
!-----------------------------------------------------------
    type(kernel), intent(out) :: bases
    integer, intent(out)      :: monomials
!-----------------------------------------------------------
    real(8)                :: x, y, z
    real(8), dimension(10) :: b, bx, by, bz, bxx, byy, bzz, bxy, bxz, byz
!-----------------------------------------------------------

    x = pnt(1)
    y = pnt(2)
    z = pnt(3)

    if (order == 1) then
        monomials = 4
    else if (order == 2) then
        monomials = 10
    else
        write(*,*) "ERROR: order for basis construction is invalid"
        call exit(0)
    end if 

    b =  [1.0d0, x, y, z, x * x, y * y, z * z, x * y, x * z, y * z]
    bx = [0.0d0, 1.0d0, 0.0d0, 0.0d0, 2.0d0*x, 0.0d0, 0.0d0, y, z, 0.0d0]
    by = [0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 2.0d0*y, 0.0d0, x, 0.0d0, z]
    bz = [0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 2.0d0*z, 0.0d0, x, y]
    bxx = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    byy = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    bzz = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0]
    bxy = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0]
    bxz = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0]
    byz = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0]

    bases%n = b(1:monomials)
    bases%nx = bx(1:monomials)
    bases%ny = by(1:monomials)
    bases%nz = bz(1:monomials)
    bases%nxx = bxx(1:monomials)
    bases%nyy = byy(1:monomials)
    bases%nzz = bzz(1:monomials)
    bases%nxy = bxy(1:monomials)
    bases%nxz = bxz(1:monomials)
    bases%nyz = byz(1:monomials)

end subroutine get_bases

subroutine get_bases_n(order, pnt, n)
!-----------------------------------------------------------
! Compute second order basis polynomials only without
! computing the derivatives
!-----------------------------------------------------------
    use krnl_struct
!-----------------------------------------------------------
    integer, intent(in) :: order
    real(8), dimension(:), intent(in) :: pnt
!-----------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: n
!-----------------------------------------------------------
    real(8) :: x, y, z
    integer :: monomials = 0
    real(8), dimension(10) :: b
!-----------------------------------------------------------
    
    x = pnt(1)
    y = pnt(2)
    z = pnt(3)

    if (order == 1) then
        monomials = 4
    else if (order == 2) then
        monomials = 10
    else
        write(*,*) "ERROR: order for basis construction is invalid"
        call exit (0)
    end if 

    b =  [1.0d0, x, y, z, x * x, y * y, z * z, x * y, x * z, y * z]
    n = b(1:monomials)

end subroutine get_bases_n

end module krnl_cmn
