module krnl_rbf
implicit none
contains

subroutine get_rbf_krnls(pnts, totpnts, sctr, wt_type, order, h, alfc, polyex, krnls)
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
    logical, intent(in) :: polyex
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) :: krnls
!-------------------------------------------------------------------------
    integer :: i
!-------------------------------------------------------------------------

    write(*,'(a,i0)') "RBF weight = " , wt_type
    write(*,'(a,f10.4)') "Kernel radius = " , h
    write(*,'(a,i0)') "Interpolation order = " , order
    
    !$omp parallel do private(i)
    do i = 1,totpnts
        call get_rbf_krnl(pnts(i, :), sctr, wt_type, order, h, alfc, polyex, krnls(i))
    end do
    !$omp end parallel do

end subroutine get_rbf_krnls

subroutine get_rbf_krnl(pnt, pnts, rbf_type, rbf_order, h, alfc, polyex, krnl)
!-------------------------------------------------------------------------
! Subroutine constructs the zero, first and second order interpolants
! using radial basis functions at an input point 'pnt' and stores 
! the computed data inside the krnl data structure
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use eqslvrs
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: rbf_type, rbf_order
    real(8), intent(in) :: h, alfc
    logical, intent(in) :: polyex
!-------------------------------------------------------------------------
    type(kernel), intent(inout) ::  krnl
!-------------------------------------------------------------------------
    real(8) :: dx, dy, dz, rr, q, rc
    real(8) :: sum_n
    real(8), dimension(:,:), allocatable :: a, inv_a, w
    real(8), dimension(:), allocatable :: bn, trm
    integer :: i, k, mg, monomials
!-------------------------------------------------------------------------
    
    if (polyex .eqv. .true.) then
        monomials = 4 !for first order interpolants
        if (rbf_order == 2) monomials = 10
    else
        monomials = 0
    end if
    mg = krnl%totnbrs + monomials

    allocate(trm(mg))
    allocate(a(mg, mg))
    allocate(inv_a(mg, mg))
    allocate(w(mg, 10))
    allocate(bn(monomials))

    a(:, :) = 0
    do i = 1, krnl%totnbrs
        do k=1, krnl%totnbrs
            dx = pnts(krnl%nbrs(i), 1) - pnts(krnl%nbrs(k), 1)
            dy = pnts(krnl%nbrs(i), 2) - pnts(krnl%nbrs(k), 2)
            dz = pnts(krnl%nbrs(i), 3) - pnts(krnl%nbrs(k), 3)
            rr = dx*dx+dy*dy+dz*dz
            select case (rbf_type)
                case (1)
                    !multi-quartics
                    q = 0.5
                    rc = alfc * h
                    a(i, k) = ((rr + rc * rc)**q)
                case (2)
                    !inverse multi-quartics	
                    q = -0.5
                    rc = alfc * h
                    a(i, k) = ((rr + rc * rc)**q)
                case (3)
                    !gaussian exponential
                    q = alfc / h / h
                    a(i, k) = exp(-q * rr)
            end select
        end do
        if (monomials > 0) then
            call get_bases_n(rbf_order, pnts(krnl%nbrs(i),:), bn)
            a(i, krnl%totnbrs + 1 : krnl%totnbrs + monomials) = bn
            a(krnl%totnbrs + 1 : krnl%totnbrs + monomials, i) = bn
        end if
    end do

    select case(rbf_type)
        case (1) !MQ
            call rbf_mq(pnt, pnts, krnl%nbrs, krnl%totnbrs, rc, q, polyex, rbf_order, w)
        case (2) !IMQ
            call rbf_mq(pnt, pnts, krnl%nbrs, krnl%totnbrs, rc, q, polyex, rbf_order, w)
        case (3) !GE
            call rbf_ge(pnt, pnts, krnl%nbrs, krnl%totnbrs, q, polyex, rbf_order, w)
    end select

    call inv(a, mg, inv_a)
    !call gauss_jordan_inv(a, mg, inv_a)

    if (rbf_order >= 1) then
        trm = matmul(inv_a, w(:, 1))
        krnl%n = trm(1:krnl%totnbrs)
        
        trm = matmul(inv_a, w(:, 2))
        krnl%nx = trm(1:krnl%totnbrs)

        trm = matmul(inv_a, w(:, 3))
        krnl%ny = trm(1:krnl%totnbrs)

        trm = matmul(inv_a, w(:, 4))
        krnl%nz = trm(1:krnl%totnbrs)
    end if
    if (rbf_order == 2) then
        trm = matmul(inv_a, w(:, 5))
        krnl%nxx = trm(1:krnl%totnbrs)

        trm = matmul(inv_a, w(:, 6))
        krnl%nyy = trm(1:krnl%totnbrs)

        trm = matmul(inv_a, w(:, 7))
        krnl%nzz = trm(1:krnl%totnbrs)

        trm = matmul(inv_a, w(:, 8))
        krnl%nxy = trm(1:krnl%totnbrs)

        trm = matmul(inv_a, w(:, 9))
        krnl%nxz = trm(1:krnl%totnbrs)

        trm = matmul(inv_a, w(:, 10))
        krnl%nyz = trm(1:krnl%totnbrs)
    end if

    sum_n = sum(krnl%n)
    if (sum_n <= 0.9 .or. sum_n >= 1.1) then
        write(*,'(a)') "ERROR in constructed RBF interpolants. Partition of Unity is Not satisfied"
        write(*,'(a,f10.3)') "Sum = ", sum_n 
        !call exit(0)
    end if

end subroutine get_rbf_krnl

subroutine rbf_mq(pnt, pnts, nbrs, totnbrs, rc, q, extended, rbf_order, w)
!-------------------------------------------------------------------------
! Compute multi-quartics radial basis functions at 'pnt'
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), dimension(:,:), intent(in) :: pnts
    integer, dimension(:), intent(in) :: nbrs
    integer, intent(in) :: totnbrs, rbf_order
    real(8), intent(in) :: rc, q
    logical, intent(in) :: extended
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: w
!-------------------------------------------------------------------------
    integer :: monomials, i
    real(8) :: dx, dy, dz, rr, rcrc, t
    type(kernel) :: bases
!-------------------------------------------------------------------------
    
    if (extended .eqv. .true.) then
        call get_bases(rbf_order, pnt, monomials, bases)
    else 
        monomials = 0
    end if
    allocate(w(totnbrs + monomials, 10))

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

    if (monomials > 0) then
        w(totnbrs + 1 : totnbrs + monomials, 1) = bases%n
        w(totnbrs + 1 : totnbrs + monomials, 2) = bases%nx
        w(totnbrs + 1 : totnbrs + monomials, 3) = bases%ny
        w(totnbrs + 1 : totnbrs + monomials, 4) = bases%nz
        w(totnbrs + 1 : totnbrs + monomials, 5) = bases%nxx
        w(totnbrs + 1 : totnbrs + monomials, 6) = bases%nyy
        w(totnbrs + 1 : totnbrs + monomials, 7) = bases%nzz
        w(totnbrs + 1 : totnbrs + monomials, 8) = bases%nxy
        w(totnbrs + 1 : totnbrs + monomials, 9) = bases%nxz
        w(totnbrs + 1 : totnbrs + monomials, 10) = bases%nyz
    end if
end subroutine rbf_mq

subroutine rbf_ge(pnt, pnts, nbrs, totnbrs, q, extended, rbf_order, w)
!-------------------------------------------------------------------------
! Compute gaussian exponential radial basis functions at 'pnt'
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(in) :: pnt
    real(8), dimension(:,:), intent(in) :: pnts
    integer, dimension(:), intent(in) :: nbrs
    integer, intent(in) :: totnbrs, rbf_order
    real(8), intent(in) :: q
    logical, intent(in) :: extended
!-------------------------------------------------------------------------
    integer :: monomials, i
    real(8) :: dx, dy, dz, rr
    type(kernel) :: bases
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: w
!-------------------------------------------------------------------------

    if (extended .eqv. .true.) then
        call get_bases(rbf_order, pnt, monomials, bases)
    else 
        monomials = 0
    end if
    allocate(w(totnbrs + monomials, 10))
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

    if (monomials > 0) then
        w(totnbrs + 1 : totnbrs + monomials, 1) = bases%n
        w(totnbrs + 1 : totnbrs + monomials, 2) = bases%nx
        w(totnbrs + 1 : totnbrs + monomials, 3) = bases%ny
        w(totnbrs + 1 : totnbrs + monomials, 4) = bases%nz
        w(totnbrs + 1 : totnbrs + monomials, 5) = bases%nxx
        w(totnbrs + 1 : totnbrs + monomials, 6) = bases%nyy
        w(totnbrs + 1 : totnbrs + monomials, 7) = bases%nzz
        w(totnbrs + 1 : totnbrs + monomials, 8) = bases%nxy
        w(totnbrs + 1 : totnbrs + monomials, 9) = bases%nxz
        w(totnbrs + 1 : totnbrs + monomials, 10) = bases%nyz
    end if
end subroutine rbf_ge

end module krnl_rbf
