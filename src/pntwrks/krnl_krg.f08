module krnl_krg
implicit none
contains

subroutine get_krg_krnls(pnts, totpnts, sctr, wt_type, order, h, alfc, krnls)
!-------------------------------------------------------------------------
! Compute the radial basis function interpolats at input points 'pnts'
! using neighbouring points from sctr points
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use omp_lib
!-------------------------------------------------------------------------
    integer, intent(in) :: wt_type, order, totpnts
    real(8), intent(in) :: h, alfc
    real(8), dimension(:,:), intent(in) :: pnts, sctr
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) :: krnls
!-------------------------------------------------------------------------
    integer :: i
!-------------------------------------------------------------------------

    write(*,'(a,i0)') "KRG weight = " , wt_type
    write(*,'(a,f10.4)') "Kernel radius = " , h
    write(*,'(a,i0)') "Interpolation Order = " , order
    
   ! !$omp parallel do private(i)
    do i = 1,totpnts
        call get_krg_krnl(pnts(i, :), sctr, wt_type, order, h, alfc, krnls(i))
    end do
   ! !$omp end parallel do

end subroutine get_krg_krnls

subroutine get_krg_krnl(pnt, pnts, krg_type, krg_order, h, alfc, krnl)
!-------------------------------------------------------------------------
! Subroutine constructs the zero, first and second order intrpolants
! using radial basis functions at an input point 'pnt' and stores 
! the computed data inside the intrp data structure
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use eqslvrs
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: krg_type, krg_order
    real(8), intent(in) :: h, alfc
!-------------------------------------------------------------------------
    type(kernel), intent(inout) ::  krnl
!-------------------------------------------------------------------------
    real(8) :: dx, dy, dz, d, q, rr, rc
    real(8) :: sum_n, rx, ry, rz, wx, wy, wz
    real(8), dimension(:,:), allocatable :: a, b, w, r_mat, inv_r
    real(8), dimension(:,:), allocatable :: trm, inv_trm, p_mat, unt_mat
    real(8), dimension(:), allocatable :: bn
    integer :: i, j, monomials
    type(kernel) :: bases
!-------------------------------------------------------------------------  
    monomials = 4 !for first order interpolants
    if (krg_order == 2) monomials = 10 !for second order interpolants

    allocate(w(krnl%totnbrs, 10))
    allocate(bn(monomials))
    allocate(r_mat(krnl%totnbrs, krnl%totnbrs))
    allocate(inv_r(krnl%totnbrs, krnl%totnbrs))
    allocate(p_mat(krnl%totnbrs, monomials))
    allocate(unt_mat(krnl%totnbrs, krnl%totnbrs))

    unt_mat(:,:) = 0.0
    do i = 1, krnl%totnbrs
        unt_mat(i,i) = 1.0
    end do

    do i = 1, krnl%totnbrs
        do j = 1, krnl%totnbrs
            dx = pnts(krnl%nbrs(i), 1) - pnts(krnl%nbrs(j), 1)
            dy = pnts(krnl%nbrs(i), 2) - pnts(krnl%nbrs(j), 2)
            dz = pnts(krnl%nbrs(i), 3) - pnts(krnl%nbrs(j), 3)
            d = sqrt(dx * dx + dy * dy + dz * dz) 
            rr = d*d
            select case (krg_type)
                case (1)
                    !multi-quartics	
                    q = 0.5
                    rc = alfc * h
                    r_mat(i, j) = ((rr + rc * rc)**q)
                case (2)
                    !inverse multi-quartics	
                    q = -0.5
                    rc = alfc * h
                    r_mat(i, j) = ((rr + rc * rc)**q)
                case (3)
                    !gaussian exponential
                    q = alfc / h / h
                    r_mat(i, j) = exp(-q * rr)
                case (4)
                    !quartic spline
                    !r = d / h
                    !r_mat(i, k) = 1 + r * r * (-6.0 + 8.0 * r - 3 * r * r)
                    rx = abs(dx) / h
                    ry = abs(dy) / h
                    rz = abs(dz) / h

                    wx = 1.0 - 6.0 * rx * rx + 8.0 * rx * rx * rx - 3.0 * rx * rx * rx * rx
                    wy = 1.0 - 6.0 * ry * ry + 8.0 * ry * ry * ry - 3.0 * ry * ry * ry * ry
                    wz = 1.0 - 6.0 * rz * rz + 8.0 * rz * rz * rz - 3.0 * rz * rz * rz * rz
                    r_mat(i, j) = wx * wy * wz 
            end select
        end do
        call get_bases_n(krg_order, pnts(krnl%nbrs(i),:), bn)
        p_mat(i, 1:monomials) = bn
    end do

    select case(krg_type)
        case (1) !multi-quartics
            call krg_mq(pnt, pnts, krnl%nbrs, krnl%totnbrs, rc, q, w)
        case(2) !inverse multi-quartics
            call krg_mq(pnt, pnts, krnl%nbrs, krnl%totnbrs, rc, q, w)
        case (3) !Gaussian Exponential
            !call krg_ge(pnt, pnts, krnl%nbrs, krnl%totnbrs, q, alfc, w)
            call krg_ge(pnt, pnts, krnl%nbrs, krnl%totnbrs, q, w)
        case (4) !Quartic Spline
            call krg_s4(pnt, pnts, krnl%nbrs, krnl%totnbrs, h, w)
    end select

    call inv(r_mat, krnl%totnbrs, inv_r)
    !call gauss_jordan_inv(r_mat, krnl%totnbrs, inv_r)

    trm = matmul(transpose(p_mat), matmul(inv_r, p_mat)) !trm = monomials X monomials
    call inv(trm, monomials, inv_trm) !inv_trm = monomials X monomaials
    !call gauss_jordan_inv(trm, monomials, inv_trm)

    a = matmul(inv_trm, matmul(transpose(p_mat), inv_r)) !a = monomials rows X nbrs clmns
    b = matmul(inv_r, unt_mat - matmul(p_mat, a)) !b = nbrs X nbrs

    call get_bases(krg_order, pnt, monomials, bases)
    if (krg_order >= 1) then
        krnl%n =  matmul(bases%n, a) + matmul(b, w(:, 1))
        krnl%nx = matmul(bases%nx, a) + matmul(b, w(:, 2))
        krnl%ny = matmul(bases%ny, a) + matmul(b, w(:, 3))
        krnl%nz = matmul(bases%nz, a) + matmul(b, w(:, 4))
    end if
    if (krg_order == 2) then
        krnl%nxx =  matmul(bases%nxx, a) + matmul(b, w(:, 5))
        krnl%nyy =  matmul(bases%nyy, a) + matmul(b, w(:, 6))
        krnl%nzz =  matmul(bases%nzz, a) + matmul(b, w(:, 7))
        krnl%nxy =  matmul(bases%nxy, a) + matmul(b, w(:, 8))
        krnl%nxz =  matmul(bases%nxz, a) + matmul(b, w(:, 9))
        krnl%nyz =  matmul(bases%nyz, a) + matmul(b, w(:, 10))
    end if

    sum_n = sum(krnl%n)
    if (sum_n <= 0.9 .or. sum_n >= 1.1) then
        write(*,'(a)') "ERROR in constructed krg intrpolants. Partition of Unity is Not satisfied"
        write(*,'(a,f10.3)') "Sum = ", sum_n 
        !call exit(0)
    end if

end subroutine get_krg_krnl

subroutine krg_s4(pnt, pnts, nbrs, totnbrs, h, w)
!-------------------------------------------------------------------------
! Compute quartic spline weight function at 'pnt'
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), dimension(:,:), intent(in) :: pnts
    integer, dimension(:), intent(in) :: nbrs
    integer, intent(in) :: totnbrs
    real(8), intent(in) :: h
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: w
!-------------------------------------------------------------------------
    integer :: i
    real(8) :: dx, dy, dz, drdx, drdy, drdz, rx, ry, rz
    real(8) :: wx, dwx, dwxx, wy, dwy, dwyy, wz, dwz, dwzz
!-------------------------------------------------------------------------
    allocate(w(totnbrs, 10))
    w(:,:) = 0.0
    do i=1, totnbrs
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
end subroutine krg_s4

subroutine krg_mq(pnt, pnts, nbrs, totnbrs, rc, q, w)
!-------------------------------------------------------------------------
! Compute multi-quartics radial basis functions at 'pnt'
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), dimension(:,:), intent(in) :: pnts
    integer, dimension(:), intent(in) :: nbrs
    integer, intent(in) :: totnbrs
    real(8), intent(in) :: rc, q
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: w
!-------------------------------------------------------------------------
    integer :: i
    real(8) :: dx, dy, dz, rr, rcrc, t
!-------------------------------------------------------------------------
    
    allocate(w(totnbrs, 10))

    w(:,:) = 0.0
    do i=1, totnbrs
        dx = pnt(1) - pnts(nbrs(i), 1)
        dy = pnt(2) - pnts(nbrs(i), 2)
        dz = pnt(3) - pnts(nbrs(i), 3)
        rr = dx * dx + dy * dy + dz * dz

        !multi-quartics
        rcrc = rc * rc
        t = rr + rcrc
        w(i, 1) = (t**q)                                                       !w
        w(i, 2) = 2.0 * q * (t**(q - 1.0)) * dx                                        !dwdx
        w(i, 3) = 2.0 * q * (t**(q - 1.0)) * dy                                        !dwdy
        w(i, 4) = 2.0 * q * (t**(q - 1.0)) * dz                                        !dwdz
        w(i, 5) = 2.0 * q * (t**(q - 1.0)) + 4.0 * (q - 1.0) * q * dx * dx * (t**(q - 2.0))     !dwdxx
        w(i, 6) = 2.0 * q * (t**(q - 1.0)) + 4.0 * (q - 1.0) * q * dy * dy * (t**(q - 2.0))     !dwdyy
        w(i, 7) = 2.0 * q * (t**(q - 1.0)) + 4.0 * (q - 1.0) * q * dz * dz * (t**(q - 2.0))     !dwdzz
        w(i, 8) = 4.0 * q * (q - 1.0) * (t**(q - 2.0)) * dx * dy     !dwdxy
        w(i, 9) = 4.0 * q * (q - 1.0) * (t**(q - 2.0)) * dx * dz     !dwdxz
        w(i, 10) = 4.0 * q * (q - 1.0) * (t**(q - 2.0)) * dy * dz    !dwdyz
    end do
end subroutine krg_mq

subroutine krg_ge(pnt, pnts, nbrs, totnbrs, q, w)
!-------------------------------------------------------------------------
! Compute gaussian exponential radial basis functions at 'pnt'
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), dimension(:,:), intent(in) :: pnts
    integer, dimension(:), intent(in) :: nbrs
    integer, intent(in) :: totnbrs
    real(8), intent(in) :: q
!-------------------------------------------------------------------------
    integer :: i
    real(8) :: dx, dy, dz, rr
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: w
!-------------------------------------------------------------------------

    allocate(w(totnbrs, 10))
    w(:, :) = 0.0
 
    do i=1, totnbrs
        dx = pnt(1) - pnts(nbrs(i), 1)
        dy = pnt(2) - pnts(nbrs(i), 2)
        dz = pnt(3) - pnts(nbrs(i), 3)
        rr = (dx * dx + dy * dy + dz * dz)

        w(i, 1) = exp(-q * rr)                                              !w
        w(i, 2) = -2. * q * exp(-q * rr) * dx                               !dwdx
        w(i, 3) = -2. * q * exp(-q * rr) * dy                               !dwdy
        w(i, 4) = -2. * q * exp(-q * rr) * dz                               !dwdz
        w(i, 5) = -2. * q * exp(-q * rr) + 4 * q * q * dx**2 * exp(-q*rr)   !dwdxx
        w(i, 6) = -2. * q * exp(-q * rr) + 4 * q * q * dy**2 * exp(-q*rr)   !dwdyy
        w(i, 7) = -2. * q * exp(-q * rr) + 4 * q * q * dz**2 * exp(-q*rr)   !dwdzz
        w(i, 8) = 4. * q * q * exp(-q * rr) * dx * dy                       !dwdxy
        w(i, 9) = 4. * q * q * exp(-q * rr) * dx * dz                       !dwdxz
        w(i, 10) = 4. * q * q * exp(-q * rr) * dy * dz                      !dwdyz
    end do
end subroutine krg_ge

end module krnl_krg
