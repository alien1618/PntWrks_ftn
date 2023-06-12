module krnl_gfd
implicit none
contains

subroutine get_gfd_krnls(pnts, totpnts, sctr, wt_type, dm, krnls)
!-------------------------------------------------------------------------
! Subroutine calculates the zero, first, and second order interpolants
! using the generalized finite difference method
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use eqslvrs
    use omp_lib
!-------------------------------------------------------------------------
    integer, intent(in)                     :: wt_type
    real(8), dimension(:,:), intent(in)     :: pnts, sctr
    real(8), intent(in)                     :: dm
    integer, intent(in)                     :: totpnts
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) :: krnls
!-------------------------------------------------------------------------
    integer                                 :: i, j, nbr
    real(8)                                 :: w, dx, dy, dz, d, gamma
    real(8), dimension(:,:), allocatable                :: a, inv_a
    real(8), dimension(:,:), allocatable    :: wmat, h
!-------------------------------------------------------------------------
    allocate(a(9,9))
    allocate(inv_a(9,9))
    do i = 1, totpnts
        allocate(wmat(krnls(i)%totnbrs, krnls(i)%totnbrs))
        allocate(h(krnls(i)%totnbrs, 9))
        wmat(:, :) = 0.0
        do j = 1, krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            dx = sctr(nbr, 1) - pnts(i, 1)
            dy = sctr(nbr, 2) - pnts(i, 2)
            dz = sctr(nbr, 3) - pnts(i, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            w = 0

            select case (wt_type)
            case (1) !cubic spline
                if (d <= 0.5 * dm) then
                    w = 0.66667d0 - 4.0 * ((d / dm)**2) - 4.0 * ((d / dm)**3)
                else if (d > 0.5 * dm .and. d <= dm) then
                    w = 1.333333d0 - 4.0 * (d / dm) + 4.0 * ((d / dm)**2) - 1.33333 * ((d / dm)**3)
                end if
            case (2) !quartic spline
                if (d <= dm) then
                    w = 1.0 - 6 * ((d / dm)**2) + 8 * ((d / dm)**3) - 3 * ((d / dm)**4)
                else
                    w = 0.0
                end if
            case (3) !gaussian exponential
                gamma = 1.0
                if (d / dm <= 1) then
                    w = exp(- gamma * d * d / (dm * dm))
                else
                    w = 0.0
                end if
            end select

            wmat(j, j) = w**2

            h(j, 1) = dx
            h(j, 2) = dy
            h(j, 3) = dz
            h(j, 4) = 0.5d0 * dx * dx
            h(j, 5) = 0.5d0 * dy * dy
            h(j, 6) = 0.5d0 * dz * dz
            h(j, 7) = dx * dy
            h(j, 8) = dx * dz
            h(j, 9) = dy * dz

        end do
        a = matmul(transpose(h), matmul(wmat,h))
        call inv(a, 9, inv_a)
        krnls(i)%gfd_inv_a = inv_a
        krnls(i)%gfd_trm = matmul(transpose(h), wmat)
        deallocate(wmat)
        deallocate(h)
    end do
end subroutine get_gfd_krnls

end module krnl_gfd
