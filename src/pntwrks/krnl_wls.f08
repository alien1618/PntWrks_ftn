module krnl_wls
implicit none
contains

subroutine get_wls_krnls(pnts, totpnts, sctr, wt_type, order, h, krnls)
!-------------------------------------------------------------------------
! Compute the weighted least squares interpolats at input points 'pnts'
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
    write(*,'(a,i0)') "WLS weight = " , wt_type
    write(*,'(a,f10.4)') "Kernel radius = " , h
    write(*,'(a,i0)') "Interpolation order = " , order

    !$omp parallel do private(i)
    do i = 1,totpnts
        call get_wls_krnl(pnts(i, :), sctr, wt_type, order, h, krnls(i))
    end do
    !$omp end parallel do
end subroutine get_wls_krnls

subroutine get_wls_krnl(pnt, pnts, wt_type, order, h, krnl)
!-------------------------------------------------------------------------
! Subroutine constructs the zero, first and second order interpolants
! using weighted least squares at an input point 'pnt' and stores 
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
    real(8), dimension(:,:), allocatable :: a, inv_a, b, p
    real(8), dimension(:), allocatable :: gamma, dgammadx, dgammady, dgammadz
    real(8), dimension(:), allocatable :: dgammadxx, dgammadyy, dgammadzz
    real(8), dimension(:), allocatable :: dgammadxy, dgammadxz, dgammadyz
    real(8), dimension(:), allocatable :: bn
    real(8), dimension(krnl%totnbrs, krnl%totnbrs) :: w
    real(8) :: dx, dy, dz, r, rx, ry, rz
    real(8) :: wx, wy, wz
    real(8) :: trm, c, w0, sum_n
    integer :: i, j, monomials
    type(kernel):: bases
!-------------------------------------------------------------------------

    call get_bases(order, pnt, monomials, bases)

    allocate(a(monomials, monomials))
    allocate(inv_a(monomials, monomials))
    allocate(b(krnl%totnbrs, monomials))
    allocate(p(monomials,krnl%totnbrs))

    do j = 1, krnl%totnbrs
        call get_bases_n(order, pnts(krnl%nbrs(j), :), bn)
        P(:,j) = bn
    end do
    c = 0.1
    trm = exp(-(1.1 * h / c) * (1.1 * h / c))
    w(:, :) = 0.0

    do i = 1, krnl%totnbrs
        dx = pnt(1) - pnts(krnl%nbrs(i), 1)
        dy = pnt(2) - pnts(krnl%nbrs(i), 2)
        dz = pnt(3) - pnts(krnl%nbrs(i), 3)
        select case (wt_type)
            case (0)
                !standard least squares method
                w(i, i) = 1
            case (1)
                !weighted least squares with gaussian exponential weights 		
                r = sqrt(dx * dx + dy * dy + dz * dz)
                w0 = exp(-(r / c) * (r / c)) - trm
                w0 = w0 / (1 - trm)
                w(i, i) = w0
            case (2)
                !weighted least squares with cbc spln weights
                rx = abs(dx) / h
                ry = abs(dy) / h
                rz = abs(dz) / h
                wx = 1.0
                wy = 1.0
                wz = 1.0
                if (rx > 0.5) then
                    wx = 1.33333333 - 4.0 * rx + 4.0 * rx * rx - 1.333333 * (rx**3.0)
                else 
                    wx = 0.6666667 - 4.0 * rx * rx + 4.0 * (rx**3.0)
                end if

                if (ry > 0.5) then
                    wy = 1.3333333 - 4.0 * ry + 4.0 * ry * ry - 1.3333333 * (ry**3.0)
                else 
                    wy = 0.6666667 - 4.0 * ry * ry + 4.0 * (ry**3.0)
                end if

                if (rz > 0.5) then
                    wz = 1.3333333 - 4.0 * rz + 4.0 * rz * rz - 1.3333333 * (rz**3.0)
                else 
                    wz = 0.6666667 - 4.0 * rz * rz + 4.0 * (rz**3.0)
                end if

                w(i, i) = wx * wy * wz   
            case (3)
                !weighted least squares with quartic spln weights
                rx = abs(dx) / h
                ry = abs(dy) / h
                rz = abs(dz) / h

                wx = 1.0 - 6.0 * rx * rx + 8.0 * rx * rx * rx - 3.0 * rx * rx * rx * rx
                wy = 1.0 - 6.0 * ry * ry + 8.0 * ry * ry * ry - 3.0 * ry * ry * ry * ry
                wz = 1.0 - 6.0 * rz * rz + 8.0 * rz * rz * rz - 3.0 * rz * rz * rz * rz

                w(i, i) = wx * wy * wz
            case (4)
                !weighted least squares with quintic spln weights
                rx = abs(dx) / h
                ry = abs(dy) / h
                rz = abs(dz) / h

                wx = 1.0 - 10.0 * rx * rx + 20.0 * rx * rx * rx - 15.0 * rx * rx * rx * rx + 4.0 * rx * rx * rx * rx * rx
                wy = 1.0 - 10.0 * ry * ry + 20.0 * ry * ry * ry - 15.0 * ry * ry * ry * ry + 4.0 * ry * ry * ry * ry * ry
                wz = 1.0 - 10.0 * rz * rz + 20.0 * rz * rz * rz - 15.0 * rz * rz * rz * rz + 4.0 * rz * rz * rz * rz * rz
                w(i, i) = wx * wy * wz
        end select
    end do

    b = matmul(w, transpose(p))
    a = matmul(p, b)
    call inv(a, monomials, inv_a)
    
    !--------------------------------------------------------------------------
    ! CALCULATE THE INTERPOLATION FUNCTIONS AND THEIR DERIVATIVES
    !--------------------------------------------------------------------------
    if (order >= 1) then
        gamma = matmul(inv_a, bases%n)
        dgammadx = matmul(inv_a, bases%nx)
        dgammady = matmul(inv_a, bases%ny)
        dgammadz = matmul(inv_a, bases%nz)

        !call slv(a, bases%n, monomials, gamma)
        !call slv(a, bases%nx, monomials, dgammadx)
        !call slv(a, bases%ny, monomials, dgammady)
        !call slv(a, bases%nz, monomials, dgammadz)
        
        krnl%n = matmul(b, gamma)
        krnl%nx = matmul(b, dgammadx)
        krnl%ny = matmul(b, dgammady)
        krnl%nz = matmul(b, dgammadz)
    end if

    if (order == 2) then
        dgammadxx = matmul(inv_a, bases%nxx)
        dgammadyy = matmul(inv_a, bases%nyy)
        dgammadzz = matmul(inv_a, bases%nzz)
        dgammadxy = matmul(inv_a, bases%nxy)
        dgammadxz = matmul(inv_a, bases%nxz)
        dgammadyz = matmul(inv_a, bases%nyz)

        !call slv(a, bases%nxx, monomials, dgammadxx)
        !call slv(a, bases%nyy, monomials, dgammadyy)
        !call slv(a, bases%nzz, monomials, dgammadzz)
        !call slv(a, bases%nxy, monomials, dgammadxy)
        !call slv(a, bases%nxz, monomials, dgammadxz)
        !call slv(a, bases%nyz, monomials, dgammadyz)

        krnl%nxx = matmul(b, dgammadxx)
        krnl%nyy = matmul(b, dgammadyy)
        krnl%nzz = matmul(b, dgammadzz)
        krnl%nxy = matmul(b, dgammadxy)
        krnl%nxz = matmul(b, dgammadxz)
        krnl%nyz = matmul(b, dgammadyz)
    end if

    sum_n = sum(krnl%n)
    if (sum_n <= 0.9 .or. sum_n >= 1.1) then
        write(*,'(a)') "ERROR in constructed WLS interpolants. Partition of Unity is Not satisfied"
        write(*,*) "Sum = ", sum_n 
        !call exit(0)
    end if
end subroutine get_wls_krnl

end module krnl_wls
