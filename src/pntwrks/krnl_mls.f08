module krnl_mls
implicit none
contains

subroutine get_mls_krnls(pnts, totpnts, sctr, wt_type, order, h, krnls)
!-------------------------------------------------------------------------
! Compute the moving least squares interpolats at input points 'pnts'
! using neighbouring points from sctr points
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use omp_lib
!-------------------------------------------------------------------------
    integer, intent(in) :: wt_type, order, totpnts
    real(8), intent(in) :: h
    real(8), dimension(:,:), intent(in) :: pnts, sctr
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) :: krnls
!-------------------------------------------------------------------------
    integer :: i
!-------------------------------------------------------------------------
    write(*,'(a,i0)') "MLS weight = " , wt_type
    write(*,'(a,f10.4)') "Kernel radius = " , h
    write(*,'(a,i0)') "Interpolation order = " , order

    !$omp parallel do private(i)
    do i = 1,totpnts
        call get_mls_krnl(pnts(i, :), sctr, wt_type, order, h, krnls(i))
    end do
    !$omp end parallel do

end subroutine get_mls_krnls

subroutine get_mls_krnl(pnt, pnts, wt_type, order, h, krnl)
!-------------------------------------------------------------------------
! Subroutine constructs the zero, first and second order intrpolants
! using moving least squares at an input point 'pnt' and stores 
! the computed data inside the krnl data structure
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use eqslvrs
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: wt_type, order
    real(8), intent(in) :: h
!-------------------------------------------------------------------------
    type(kernel), intent(inout) :: krnl
!-------------------------------------------------------------------------
    type(kernel) :: p, bases
    integer :: monomials
    real(8), dimension(:,:), allocatable :: pp, a, inv_a, dadx, dady, dadz
    real(8), dimension(:,:), allocatable :: dadxx, dadyy, dadzz, dadxy, dadxz, dadyz
    real(8), dimension(:,:), allocatable :: b, dbdx, dbdy, dbdz
    real(8), dimension(:,:), allocatable :: dbdxx, dbdyy, dbdzz, dbdxy, dbdxz, dbdyz
    real(8), dimension(:,:), allocatable :: weights
    real(8), dimension(:), allocatable :: p_dadx_gamma, p_dady_gamma, p_dadz_gamma
    real(8), dimension(:), allocatable :: trm_xx, trm_yy, trm_zz, trm_xy, trm_xz, trm_yz
    real(8), dimension(:), allocatable :: gamma, dgammadx, dgammady, dgammadz
    real(8), dimension(:), allocatable :: dgammadxx, dgammadyy, dgammadzz, dgammadxy, dgammadxz, dgammadyz
    integer :: i, m, n
    real(8) :: sum_n
!-------------------------------------------------------------------------
    ! Calculate bases at point
    call get_bases(order, pnt, monomials, bases)
!-------------------------------------------------------------------------
    ! memory allocation
    allocate(pp(monomials, monomials))
    allocate(A(monomials, monomials))
    allocate(inv_a(monomials, monomials))
    allocate(dAdx(monomials, monomials))
    allocate(dAdy(monomials, monomials))
    allocate(dAdz(monomials, monomials))
    allocate(B(krnl%totnbrs, monomials))
    allocate(dBdx(krnl%totnbrs, monomials))
    allocate(dBdy(krnl%totnbrs, monomials))
    allocate(dBdz(krnl%totnbrs, monomials))   
    allocate(weights(krnl%totnbrs, monomials))
!-------------------------------------------------------------------------
    ! generate weight functions
    select case (wt_type)
        case (1) 
            call mls_spln3(pnt, pnts, krnl%nbrs, krnl%totnbrs, h, weights)
        case (2)
            call mls_spln4(pnt, pnts, krnl%nbrs, krnl%totnbrs, h, weights)
        case (3)
            call mls_spln5(pnt, pnts, krnl%nbrs, krnl%totnbrs, h, weights)
        case (4)
            call mls_reg_spln4(pnt, pnts, krnl%nbrs, krnl%totnbrs, h, weights)
    end select
    
!-------------------------------------------------------------------------
    ! calculate the A and B matrices and its derivatives
    A(:, :) = 0.0
    dAdx(:, :) = 0.0
    dAdy(:, :) = 0.0
    dAdz(:, :) = 0.0
    do i = 1, krnl%totnbrs
        call get_bases(order, pnts(krnl%nbrs(i), :), monomials, p)

        !calculate p*p_transpose
        do n = 1, monomials
            do m = 1, monomials
                pp(m,n) = p%n(m) * p%n(n)
            end do
        end do

        A = A + pp * weights(i, 1)
        dAdx = dAdx + pp * weights(i, 2)
        dAdy = dAdy + pp * weights(i, 3)
        dAdz = dAdz + pp * weights(i, 4)

        B(i,:) = p%n * weights(i, 1)
        dBdx(i,:) = p%n * weights(i, 2)
        dBdy(i,:) = p%n * weights(i, 3)
        dBdz(i,:) = p%n * weights(i, 4)
    end do

!-------------------------------------------------------------------------
    ! calculate the interpolation functions and its derivatives
    call inv(A, monomials, inv_a)
    gamma = matmul(inv_a, bases%n)
    !call slv(A, bases%N, monomials, gamma)
    
    p_dAdx_gamma = bases%nx - matmul(dAdx, gamma)
    p_dAdy_gamma = bases%ny - matmul(dAdy, gamma)
    p_dAdz_gamma = bases%nz - matmul(dAdz, gamma)
       
    dgammadx = matmul(inv_a, p_dAdx_gamma)
    dgammady = matmul(inv_a, p_dAdy_gamma)
    dgammadz = matmul(inv_a, p_dAdz_gamma)

    krnl%n = matmul(B, gamma)
    krnl%nx = matmul(B, dgammadx) + matmul(dBdx, gamma)
    krnl%ny = matmul(B, dgammady) + matmul(dBdy, gamma)
    krnl%nz = matmul(B, dgammadz) + matmul(dBdz, gamma)
!-------------------------------------------------------------------------
    ! calculate the second order intrpolants if needed
    if (order == 2) then
        allocate(dAdxx(monomials, monomials))
        allocate(dAdyy(monomials, monomials))
        allocate(dAdzz(monomials, monomials))
        allocate(dAdxy(monomials, monomials))
        allocate(dAdxz(monomials, monomials))
        allocate(dAdyz(monomials, monomials))

        allocate(dBdxx(krnl%totnbrs, monomials))
        allocate(dBdyy(krnl%totnbrs, monomials))
        allocate(dBdzz(krnl%totnbrs, monomials))
        allocate(dBdxy(krnl%totnbrs, monomials))
        allocate(dBdxz(krnl%totnbrs, monomials))
        allocate(dBdyz(krnl%totnbrs, monomials))

        dAdxx(:, :) = 0.0
        dAdyy(:, :) = 0.0
        dAdzz(:, :) = 0.0
        dAdxy(:, :) = 0.0
        dAdxz(:, :) = 0.0
        dAdyz(:, :) = 0.0
        do i = 1, krnl%totnbrs
            call get_bases(order, pnts(krnl%nbrs(i), :), monomials, p)

            !calculate p*p_transpose
            do n = 1, monomials
                do m = 1, monomials
                    pp(m,n) = p%n(m) * p%n(n)
                end do
            end do
            dAdxx = dAdxx + pp * weights(i, 5)
            dAdyy = dAdyy + pp * weights(i, 6)
            dAdzz = dAdzz + pp * weights(i, 7)
            dAdxy = dAdxy + pp * weights(i, 8)
            dAdxz = dAdxz + pp * weights(i, 9)
            dAdyz = dAdyz + pp * weights(i, 10)

            dBdxx(i,:) = p%n * weights(i, 5)
            dBdyy(i,:) = p%n * weights(i, 6)
            dBdzz(i,:) = p%n * weights(i, 7)
            dBdxy(i,:) = p%n * weights(i, 8)
            dBdxz(i,:) = p%n * weights(i, 9)
            dBdyz(i,:) = p%n * weights(i, 10)
        end do
        trm_xx = bases%nxx - (2 * matmul(dAdx, dgammadx) + matmul(dAdxx, gamma))
        trm_yy = bases%nyy - (2 * matmul(dAdy, dgammady) + matmul(dAdyy, gamma))
        trm_zz = bases%nzz - (2 * matmul(dAdz, dgammadz) + matmul(dAdzz, gamma))
        trm_xy = bases%nxy - (matmul(dAdx, dgammady) + matmul(dAdy, dgammadx) + matmul(dAdxy, gamma))
        trm_xz = bases%nxz - (matmul(dAdx, dgammadz) + matmul(dAdz, dgammadx) + matmul(dAdxz, gamma))
        trm_yz = bases%nyz - (matmul(dAdy, dgammadz) + matmul(dAdz, dgammady) + matmul(dAdyz, gamma))

        dgammadxx = matmul(inv_a, trm_xx)
        dgammadyy = matmul(inv_a, trm_yy)
        dgammadzz = matmul(inv_a, trm_zz)
        dgammadxy = matmul(inv_a, trm_xy)
        dgammadxz = matmul(inv_a, trm_xz)
        dgammadyz = matmul(inv_a, trm_yz)

        krnl%nxx = matmul(B, dgammadxx) + 2 * matmul(dBdx, dgammadx) + matmul(dBdxx, gamma)
        krnl%nyy = matmul(B, dgammadyy) + 2 * matmul(dBdy, dgammady) + matmul(dBdyy, gamma)
        krnl%nzz = matmul(B, dgammadzz) + 2 * matmul(dBdz, dgammadz) + matmul(dBdzz, gamma)
        krnl%nxy = matmul(B, dgammadxy) + matmul(dBdy, dgammadx) + matmul(dBdx, dgammady) + matmul(dBdxy, gamma)
        krnl%nxz = matmul(B, dgammadxz) + matmul(dBdz, dgammadx) + matmul(dBdx, dgammadz) + matmul(dBdxz, gamma)
        krnl%nyz = matmul(B, dgammadyz) + matmul(dBdz, dgammady) + matmul(dBdy, dgammadz) + matmul(dBdyz, gamma)
    end if
!-------------------------------------------------------------------------
    ! check the partition of unity is satisfied
    sum_n = sum(krnl%n)
    if (sum_n <= 0.9 .or. sum_n >= 1.1) then
        write(*,'(a)') "ERROR in constructed MLS intrpolants. Partition of Unity is Not satisfied"
        write(*,'(a,f10.3)') "Sum = ", sum_n 
        !call exit(0)
    end if

end subroutine get_mls_krnl

subroutine mls_spln3(pnt, pnts, nbrs, totnbrs, h, w)
!-------------------------------------------------------------------------
! cbc spln WEIGHT FUNCTION
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in)   :: pnt
    real(8), intent(in)                 :: h
    real(8), dimension(:,:), intent(in) :: pnts
    integer, dimension(:), intent(in)   :: nbrs
    integer, intent(in)                 :: totnbrs
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: w
!-------------------------------------------------------------------------
    real(8) :: wx = 0.0, dwx = 0.0, dwxx = 0.0
    real(8) :: wy = 0.0, dwy = 0.0, dwyy = 0.0
    real(8) :: wz = 0.0, dwz = 0.0, dwzz = 0.0
    real(8) :: dx, dy, dz, drdx, drdy, drdz, rx, ry, rz
    integer :: i
!-------------------------------------------------------------------------

    allocate(w(totnbrs, 10))
    do i = 1, totnbrs
        dx = pnt(1) - pnts(nbrs(i), 1)
        dy = pnt(2) - pnts(nbrs(i), 2)
        dz = pnt(3) - pnts(nbrs(i), 3)

        drdx = sign(1.0d0, dx) / h
        drdy = sign(1.0d0, dy) / h
        drdz = sign(1.0d0, dz) / h
        
        rx = abs(dx) / h
        ry = abs(dy) / h
        rz = abs(dz) / h
        
        if (rx > 0.5) then
            wx = 1.33333333 - 4.0 * rx + 4.0 * rx * rx - 1.333333 * (rx**3)
            dwx = (-4.0 + 8.0 * rx - 4.0 * (rx**2)) * drdx
            dwxx = (8.0 - 8.0 * rx) * drdx * drdx
        else if (rx <= 0.5) then
            wx = 0.6666667 - 4.0 * rx * rx + 4 * (rx**3)
            dwx = (-8.0 * rx + 12.0 * (rx**2)) * drdx
            dwxx = (-8.0 + 24.0 * rx) * drdx * drdx
        end if

        if (ry > 0.5) then
            wy = 1.3333333 - 4.0 * ry + 4.0 * ry * ry -1.3333333 * (ry**3)
            dwy = (-4.0 + 8.0 * ry - 4.0 * (ry**2)) * drdy
            dwyy = (8.0 - 8.0 * ry) * drdy * drdy
        else if (ry <= 0.5) then
            wy = 0.6666667 - 4.0 * ry * ry + 4.0 * (ry**3)
            dwy = (-8.0 * ry + 12.0 * (ry**2)) * drdy
            dwyy = (-8.0 + 24.0 * ry) * drdy * drdy
        end if

        if (rz > 0.5) then
            wz = 1.3333333 - 4.0 * rz + 4.0 * rz * rz - 1.3333333 * (rz**3)
            dwz = (-4.0 + 8.0 * rz - 4 * (rz**2)) * drdz
            dwzz = (8.0 - 8.0 * rz) * drdz * drdz
        else if (rz <= 0.5) then
            wz = 0.6666667 - 4.0 * rz * rz + 4.0 * (rz**3)
            dwz = (-8.0 * rz + 12.0 * (rz**2)) * drdz
            dwzz = (-8.0 + 24.0 * rz) * drdz * drdz
        end if

        w(i, 1) = wx * wy * wz   !w
        w(i, 2) = wy * wz * dwx  !dwdx
        w(i, 3) = wx * wz * dwy  !dwdy
        w(i, 4) = wx * wy * dwz  !dwdz
        w(i, 5) = wy * wz * dwxx   !dwdxx
        w(i, 6) = wx * wz * dwyy   !dwdyy
        w(i, 7) = wx * wy * dwzz   !dwdzz
        w(i, 8) = dwx * dwy * wz   !dwdxdwdy
        w(i, 9) = dwx * wy * dwz   !dwdxdwdz
        w(i, 10) = wx * dwy * dwz  !dwdydwdz
    end do
end subroutine mls_spln3

subroutine mls_spln4(pnt, pnts, nbrs, totnbrs, h, w)
!-------------------------------------------------------------------------
! QUARTIC SPLING WEIGHT FUNCTIONS
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in)   :: pnt
    real(8), intent(in)                 :: h
    real(8), dimension(:,:), intent(in) :: pnts
    integer, dimension(:), intent(in)   :: nbrs
    integer, intent(in)                 :: totnbrs
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: w
!-------------------------------------------------------------------------
    real(8) :: wx, dwx, dwxx, wy, dwy, dwyy, wz, dwz, dwzz
    real(8) :: dx, dy, dz, drdx, drdy, drdz, rx, ry, rz
    integer :: i
!-------------------------------------------------------------------------

    allocate(w(totnbrs, 10))
    do i = 1, totnbrs
        dx = pnt(1) - pnts(nbrs(i), 1)
        dy = pnt(2) - pnts(nbrs(i), 2)
        dz = pnt(3) - pnts(nbrs(i), 3)

        drdx = sign(1.0d0, dx) / h
        drdy = sign(1.0d0, dy) / h
        drdz = sign(1.0d0, dz) / h

        rx = abs(dx) / h
        ry = abs(dy) / h
        rz = abs(dz) / h

        wx = 1.0 - 6.0 * rx * rx + 8.0 * rx * rx * rx - 3.0 * rx * rx * rx * rx
        dwx = (-12.0 * rx + 24.0 * rx * rx - 12.0 * rx * rx * rx) * drdx
        dwxx = (-12.0 + 48.0 * rx - 36.0 * rx * rx) * drdx * drdx

        wy = 1.0 - 6.0 * ry * ry + 8.0 * ry * ry * ry - 3.0 * ry * ry * ry * ry
        dwy = (-12.0 * ry + 24.0 * ry * ry - 12.0 * ry * ry * ry) * drdy
        dwyy = (-12.0 + 48.0 * ry - 36.0 * ry * ry) * drdy * drdy

        wz = 1.0 - 6.0 * rz * rz + 8.0 * rz * rz * rz - 3.0 * rz * rz * rz * rz
        dwz = (-12.0 * rz + 24.0 * rz * rz - 12.0 * rz * rz * rz) * drdz
        dwzz = (-12.0 + 48.0 * rz - 36.0 * rz * rz) * drdz * drdz

        w(i, 1) = wx * wy * wz   !w
        w(i, 2) = wy * wz * dwx  !dwdx
        w(i, 3) = wx * wz * dwy  !dwdy
        w(i, 4) = wx * wy * dwz  !dwdz
        w(i, 5) = wy * wz * dwxx   !dwdxx
        w(i, 6) = wx * wz * dwyy   !dwdyy
        w(i, 7) = wx * wy * dwzz   !dwdzz
        w(i, 8) = dwx * dwy * wz   !dwdxdwdy
        w(i, 9) = dwx * wy * dwz   !dwdxdwdz
        w(i, 10) = wx * dwy * dwz  !dwdydwdz
    end do
end subroutine mls_spln4

subroutine mls_spln5(pnt, pnts, nbrs, totnbrs, h, w)
!-------------------------------------------------------------------------
! FIFTH ORDER splnS WEIGHT FUNCTIONS
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), intent(in)               :: h
    real(8), dimension(:,:), intent(in) :: pnts
    integer, dimension(:), intent(in)   :: nbrs
    integer, intent(in) :: totnbrs
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: w
!-------------------------------------------------------------------------
    real(8) :: wx, dwx, dwxx, wy, dwy, dwyy, wz, dwz, dwzz
    real(8) :: dx, dy, dz, drdx, drdy, drdz, rx, ry, rz
    integer :: i
!-------------------------------------------------------------------------

    allocate(w(totnbrs, 10))
    do i = 1, totnbrs
        dx = pnt(1) - pnts(nbrs(i), 1)
        dy = pnt(2) - pnts(nbrs(i), 2)
        dz = pnt(3) - pnts(nbrs(i), 3)

        drdx = sign(1.0d0, dx) / h
        drdy = sign(1.0d0, dy) / h
        drdz = sign(1.0d0, dz) / h

        rx = abs(dx) / h
        ry = abs(dy) / h
        rz = abs(dz) / h

        wx= 1.0 - 10.0 * rx * rx + 20.0 * rx * rx * rx - 15.0 * rx * rx * rx * rx + 4.0 * rx * rx * rx * rx * rx
        dwx = (-20.0 * rx + 60.0 * rx * rx - 60.0 * rx * rx * rx + 20.0 * rx * rx * rx * rx) * drdx
        dwxx =(-20.0 + 120.0 * rx - 180.0 * rx * rx + 80.0 * rx * rx * rx) * drdx * drdx

        wy = 1.0 - 10.0 * ry * ry + 20.0 * ry * ry * ry - 15.0 * ry * ry * ry * ry + 4.0 * ry * ry * ry * ry * ry
        dwy = (-20.0 * ry + 60.0 * ry * ry - 60.0 * ry * ry * ry + 20.0 * ry * ry * ry * ry) * drdy
        dwyy =(-20.0 + 120.0 * ry - 180.0 * ry * ry + 80.0 * ry * ry * ry) * drdy * drdy

        wz = 1.0 - 10.0 * rz * rz + 20.0 * rz * rz * rz - 15.0 * rz * rz * rz * rz + 4.0 * rz * rz * rz * rz * rz
        dwz = (-20.0 * rz + 60.0 * rz * rz - 60.0 * rz * rz * rz + 20.0 * rz * rz * rz * rz) * drdz
        dwzz =(-20.0 + 120.0 * rz - 180.0 * rz * rz + 80.0 * rz * rz * rz) * drdz * drdz

        w(i, 1) = wx * wy * wz   !w
        w(i, 2) = wy * wz * dwx  !dwdx
        w(i, 3) = wx * wz * dwy  !dwdy
        w(i, 4) = wx * wy * dwz  !dwdz
        w(i, 5) = wy * wz * dwxx   !dwdxx
        w(i, 6) = wx * wz * dwyy   !dwdyy
        w(i, 7) = wx * wy * dwzz   !dwdzz
        w(i, 8) = dwx * dwy * wz   !dwdxdwdy
        w(i, 9) = dwx * wy * dwz   !dwdxdwdz
        w(i, 10) = wx * dwy * dwz  !dwdydwdz
    end do
end subroutine mls_spln5

subroutine mls_reg_spln4(pnt, pnts, nbrs, totnbrs, h, w)
!-------------------------------------------------------------------------
! REGULARIZED QUARTIC splnS WEIGHT FUNCTIONS
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), intent(in) :: h
    real(8), dimension(:,:), intent(in) :: pnts
    integer, dimension(:), intent(in) :: nbrs
    integer, intent(in) :: totnbrs
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: w
!-------------------------------------------------------------------------
    real(8) :: wx, wy, wz
    real(8) :: dwx, dwy, dwz, dwxx, dwyy, dwzz, eps = 1e-3, eps2 = 1e-8
    real(8) :: dx, dy, dz, rx, ry, rz
    real(8) :: trm, trm1_x, trm1_y, trm1_z, trm2_x, trm2_y, trm2_z, trm3
    integer :: i
!-------------------------------------------------------------------------

    allocate(w(totnbrs, 10))
    do i=1, totnbrs
        dx = pnt(1) - pnts(nbrs(i), 1)
        dy = pnt(2) - pnts(nbrs(i), 2)
        dz = pnt(3) - pnts(nbrs(i), 3)

        rx = abs(dx) / h
        ry = abs(dy) / h
        rz = abs(dz) / h

        trm = ((eps**-2) - ((1.0 + eps)**(-2.)))
        
        wx = (((rx**2) + eps)**-2) - ((1.0 + eps)**(-2.))
        wx = wx / trm

        wy = (((ry**2) + eps)**-2) - ((1.0 + eps)**(-2.))
        wy = wy / trm

        wz = (((rz**2) + eps)**-2) - ((1.0 + eps)**(-2.))
        wz = wz / trm

        trm1_x = (((rx**2.) + eps)**(-3.))
        trm1_y = (((ry**2.) + eps)**(-3.))
        trm1_z = (((rz**2.) + eps)**(-3.))

        trm2_x = (((rx**2.) + eps)**(-4.))
        trm2_y = (((ry**2.) + eps)**(-4.))
        trm2_z = (((rz**2.) + eps)**(-4.))

        trm3 = (eps**-2.) - ((1 + eps)**(-2.))

        if (dx >= eps2) then
            dwx = 4.0 * trm1_x * rx / (h * trm3)
            dwxx = -24.0 * ((rx**2.0) / (h**2.0)) * trm2_x - (4.0 / (h**2.0)) * trm1_x
            dwxx = dwxx / trm3
        else
            dwx = 4.0 * (eps**2.0) * (1.0 / h)
            dwxx = 20.0 * (eps**2.0) * (1.0 / (h**2.0))
        end if

        if (dy >= eps2) then
            dwy = 4.0 * trm1_y * ry / (h * trm3)
            dwyy = -24.0 * ((ry**2.0) / (h**2.0)) * trm2_y + (4.0 / (h**2.0)) * trm1_y
            dwyy = dwyy / trm3
        else
            dwy = 4.0 * (eps**2.0) * (1.0 / h)
            dwyy = 20.0 * (eps**2.0) * (1.0 / (h**2.0))
        end if

        if (dz >= eps2) then
            dwz = 4.0 * trm1_z * rz / (h * trm3)
            dwzz = -24.0 * ((rz**2.0) / (h**2.0)) * trm2_z + (4.0 / (h**2.0)) * trm1_z
            dwzz = dwzz / trm3
        else
            dwz = 4.0 * (eps**2.0) * (1.0 / h)
            dwzz = 20.0 * (eps**2.0) * (1.0 / (h**2.0))
        end if

        w(i, 1) = wx * wy * wz   !w
        w(i, 2) = wy * wz * dwx  !dwdx
        w(i, 3) = wx * wz * dwy  !dwdy
        w(i, 4) = wx * wy * dwz  !dwdz
        w(i, 5) = wy * wz * dwxx   !dwdxx
        w(i, 6) = wx * wz * dwyy   !dwdyy
        w(i, 7) = wx * wy * dwzz   !dwdzz
        w(i, 8) = dwx * dwy * wz   !dwdxdwdy
        w(i, 9) = dwx * wy * dwz   !dwdxdwdz
        w(i, 10) = wx * dwy * dwz  !dwdydwdz
    end do
end subroutine mls_reg_spln4

end module krnl_mls
