module krnl_sph
implicit none
contains

subroutine get_sph_krnls(pnts, totpnts, sctr, wt_type, dim0, h, krnls)
!-------------------------------------------------------------------------
! Compute smoothed kernel approximants at input points 'pnts'
! using neighbouring points from sctr points
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
!-------------------------------------------------------------------------
    real(8), dimension(:,:),intent(in)                      ::  pnts, sctr
    integer, intent(in)                                     ::  wt_type, dim0, totpnts
    real(8), intent(in)                                     ::  h
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout)  ::  krnls
!-------------------------------------------------------------------------

    !write(*,'(a,i0)') "SPH weight = " , wt_type
    !write(*,'(a,f10.4)') "Kernel radius = " , h

    select case (wt_type)
        case (1)
            call get_sph_exp(dim0, h, pnts, totpnts, sctr, krnls)
        case (2)
            call get_sph_cbc_spln(dim0, h, pnts, totpnts, sctr, krnls)
        case (3)
            call get_sph_wndlnd4(dim0, h, pnts, totpnts, sctr, krnls)
        case (4)
            call get_sph_wndlnd5(dim0, h, pnts, totpnts, sctr, krnls)
        case (5)
            call get_sph_inv_dist(h, pnts, totpnts, sctr, krnls)
        case default
            write(*,'(a)') " Error in get_sph_krnls. Chosen wt_type is invalid"
    end select
end subroutine

subroutine get_sph_exp(dim0, h, pnts, totpnts, sctr, krnls)
!-------------------------------------------------------------------------
! gaussian exponential weight functions
!-------------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!-------------------------------------------------------------------------
    real(8), intent(in)                 :: h
    integer, intent(in)                 :: dim0, totpnts
    real(8), dimension(:,:),intent(in)  :: pnts, sctr
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) ::  krnls
!-------------------------------------------------------------------------
    integer :: i, j, nbr, totnbrs
    real(8) :: pi, eps, alpha, prmtr, dx, dy, dz, d
    real(8) :: q, qx, qy, qz, val, sumn
!-------------------------------------------------------------------------

    pi = 3.14
    eps = 1.0e-6
    prmtr = 1
    alpha = 1.0 / ((pi**(0.5 * dim0)) * (h**dim0))

    !computing intrpolants
    !$omp parallel do private(i)
    do i = 1, totpnts
        totnbrs = krnls(i)%totnbrs
        allocate(krnls(i)%n(totnbrs))
        allocate(krnls(i)%nx(totnbrs))
        allocate(krnls(i)%ny(totnbrs))
        allocate(krnls(i)%nz(totnbrs))
        do j = 1, totnbrs
            nbr = krnls(i)%nbrs(j)
            dx = pnts(i, 1) - sctr(nbr, 1)
            dy = pnts(i, 2) - sctr(nbr, 2)
            dz = pnts(i, 3) - sctr(nbr, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz) + eps
            q = d / h
            qx = dx / (d * h)
            qy = dy / (d * h)
            qz = dz / (d * h)
            if (q .ge. 0 .and. q .le. 3) then
                val = alpha * exp(-q*q)
                krnls(i)%n(j) = val
                krnls(i)%nx(j) = -2.0 * qx * val
                krnls(i)%ny(j) = -2.0 * qy * val
                krnls(i)%nz(j) = -2.0 * qz * val
            else
                krnls(i)%n(j) = 0.0
                krnls(i)%nx(j) = 0.0
                krnls(i)%ny(j) = 0.0
                krnls(i)%nz(j) = 0.0
            end if
        end do
        !normalizing intrpolants
        sumn = sum(krnls(i)%n)
        krnls(i)%n = krnls(i)%n / sumn
        krnls(i)%nx = krnls(i)%nx / sumn
        krnls(i)%ny = krnls(i)%ny / sumn
        krnls(i)%nz = krnls(i)%nz / sumn
    end do
    !$omp end parallel do
end subroutine get_sph_exp

subroutine get_sph_cbc_spln(dim0, h, pnts, totpnts, sctr, krnls)
!-------------------------------------------------------------------------
! cubic spline weight functions
!-------------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!-------------------------------------------------------------------------
    real(8), intent(in)                       :: h
    integer, intent(in)                       :: dim0, totpnts
    real(8), dimension(:,:), intent(in)       :: pnts, sctr
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) :: krnls
!-------------------------------------------------------------------------
    integer :: i, j, nbr, totnbrs
    real(8) :: pi, eps, alpha, dx, dy, dz, d
    real(8) :: q, qx, qy, qz, sumn
!-------------------------------------------------------------------------
    
    pi = 3.14;
    alpha = 1;
    eps = 1.0e-6
    select case (dim0)
        case (1)
            alpha = 1 / (6.0 * h)
        case (2)
            alpha = 15.0 / (14 * pi * h * h)
        case (3)
            alpha = 1.0 / (4 * pi * h * h * h)
        case default
            write(*,'(a)') " Error in get_sph_cbc_spln. Dim is invalid"
    end select
   
    !computing intrpolants
    !$omp parallel do private(i)
    do i = 1, totpnts
        totnbrs = krnls(i)%totnbrs
        allocate(krnls(i)%n(totnbrs))
        allocate(krnls(i)%nx(totnbrs))
        allocate(krnls(i)%ny(totnbrs))
        allocate(krnls(i)%nz(totnbrs))
        do j = 1, totnbrs
            nbr = krnls(i)%nbrs(j)
            dx = pnts(i, 1) - sctr(nbr, 1)
            dy = pnts(i, 2) - sctr(nbr, 2)
            dz = pnts(i, 3) - sctr(nbr, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz) + eps
            q = d / h
            qx = dx / (d * h)
            qy = dy / (d * h)
            qz = dz / (d * h)
            if (q >= 0 .and. q < 1) then
                krnls(i)%n(j) = alpha * (((2 - q)**3) - 4 * ((1.0 - q)**3))
                krnls(i)%nx(j) = alpha * (-3 * qx * ((2 - q)**2) + 12 * qx * ((1 - q)**2))
                krnls(i)%ny(j) = alpha * (-3 * qy * ((2 - q)**2) + 12 * qy * ((1 - q)**2))
                krnls(i)%nz(j) = alpha * (-3 * qz * ((2 - q)**2) + 12 * qz * ((1 - q)**2))
            else if(q >= 1 .and. q < 2) then
                krnls(i)%n(j) = alpha * (((2 - q)**3))
                krnls(i)%nx(j) = alpha * (-3 * qx * ((2 - q)**2))
                krnls(i)%ny(j) = alpha * (-3 * qy * ((2 - q)**2))
                krnls(i)%nz(j) = alpha * (-3 * qz * ((2 - q)**2))
            else if(q >= 2) then
                krnls(i)%n(j) = 0
                krnls(i)%nx(j) = 0
                krnls(i)%ny(j) = 0
                krnls(i)%nz(j) = 0
            end if
        end do
        !normalizing intrpolants
        sumn = sum(krnls(i)%n)
        krnls(i)%n = krnls(i)%n / sumn
        krnls(i)%nx = krnls(i)%nx / sumn
        krnls(i)%ny = krnls(i)%ny / sumn
        krnls(i)%nz = krnls(i)%nz / sumn
    end do
    !$omp end parallel do
end subroutine get_sph_cbc_spln

subroutine get_sph_wndlnd4(dim0, h, pnts, totpnts, sctr, krnls)
!-------------------------------------------------------------------------
! 4th order wendland weight function
!-------------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!-------------------------------------------------------------------------
    real(8), intent(in)                  :: h
    integer, intent(in)                  ::  dim0, totpnts
    real(8), dimension(:,:), intent(in)  ::  pnts, sctr
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) ::  krnls
!-------------------------------------------------------------------------
    integer :: i, j, nbr, totnbrs
    real(8) :: pi, eps, alpha, dx, dy, dz, r
    real(8) :: q, qx, qy, qz, sumn
!-------------------------------------------------------------------------
    
    pi = 3.14;
    alpha = 1;
    eps = 1.0e-6
    select case (dim0)
        case (1)
            alpha = 3.0 / (64.0 * h)
        case (2)
            alpha = 7.0 / (64.0 * pi * h * h)
        case (3)
            alpha = 21.0 / (16.0 * pi * h * h * h)
        case default
            write(*,'(a)') " Error in get_sph_wndlnd4. Dim is invalid"
    end select

   !computing intrpolants
   !$omp parallel do private(i)
   do i = 1, totpnts
        totnbrs = krnls(i)%totnbrs
        allocate(krnls(i)%n(totnbrs))
        allocate(krnls(i)%nx(totnbrs))
        allocate(krnls(i)%ny(totnbrs))
        allocate(krnls(i)%nz(totnbrs))
        do j = 1, totnbrs
            nbr = krnls(i)%nbrs(j)
            dx = pnts(i, 1) - sctr(nbr, 1)
            dy = pnts(i, 2) - sctr(nbr, 2)
            dz = pnts(i, 3) - sctr(nbr, 3)
            r = sqrt(dx * dx + dy * dy + dz * dz) + eps
            q = r / h
            qx = dx / (r * h)
            qy = dy / (r * h)
            qz = dz / (r * h)

            if (q >= 0 .and. q <= 2) then
                krnls(i)%n(j) = alpha * ((2 - q)**4) * (2 * q + 1)
                krnls(i)%nx(j) = -alpha * 4 * ((2 - q)**3) * qx * (2 * q + 1) + alpha * ((2 - q)**4) * (2 * qx)
                krnls(i)%ny(j) = -alpha * 4 * ((2 - q)**3) * qy * (2 * q + 1) + alpha * ((2 - q)**4) * (2 * qy)
                krnls(i)%nz(j) = -alpha * 4 * ((2 - q)**3) * qz * (2 * q + 1) + alpha * ((2 - q)**4) * (2 * qz)
            else
                krnls(i)%n(j) = 0
                krnls(i)%nx(j) = 0
                krnls(i)%ny(j) = 0
                krnls(i)%nz(j) = 0
            end if
        end do
        !normalizing intrpolants
        sumn = sum(krnls(i)%n) + eps
        krnls(i)%n = krnls(i)%n / sumn
        krnls(i)%nx = krnls(i)%nx / sumn
        krnls(i)%ny = krnls(i)%ny / sumn
        krnls(i)%nz = krnls(i)%nz / sumn   
    end do
    !$omp end parallel do
end subroutine get_sph_wndlnd4

subroutine get_sph_wndlnd5(dim0, h, pnts, totpnts, sctr, krnls)
!-------------------------------------------------------------------------
! 5th order wendland weight function
!-------------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!-------------------------------------------------------------------------
    real(8), intent(in)                       :: h
    integer, intent(in)                       :: dim0, totpnts
    real(8), dimension(:,:), intent(in)       :: pnts, sctr
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) ::  krnls
!-------------------------------------------------------------------------
    integer :: i, j, nbr, totnbrs
    real(8) :: pi, eps, alpha, dx, dy, dz, r
    real(8) :: q, qx, qy, qz, sumn
!-------------------------------------------------------------------------
    
    pi = 3.14;
    alpha = 1;
    eps = 1.0e-6
    select case (dim0)
        case (2)
            alpha = 9.0 / (4.0 * pi * h * h)
        case (3)
            alpha = 495.0 / (256.0 * pi * h * h * h)
        case default
            write(*,'(a)') " Error in get_sph_wndlnd5. Dim is invalid"
    end select

   !computing intrpolants
   !$omp parallel do private(i)
    do i = 1, totpnts
        totnbrs = krnls(i)%totnbrs
        allocate(krnls(i)%n(totnbrs))
        allocate(krnls(i)%nx(totnbrs))
        allocate(krnls(i)%ny(totnbrs))
        allocate(krnls(i)%nz(totnbrs))
        do j = 1,totnbrs
            nbr = krnls(i)%nbrs(j)
            dx = pnts(i, 1) - sctr(nbr, 1)
            dy = pnts(i, 2) - sctr(nbr, 2)
            dz = pnts(i, 3) - sctr(nbr, 3)
            r = sqrt(dx * dx + dy * dy + dz * dz) + eps
            q = r / h
            qx = dx / (r * h)
            qy = dy / (r * h)
            qz = dz / (r * h)
            if (q >= 0 .and. q <= 2) then
                krnls(i)%n(j) = alpha * ((2 - q)**6) * ((35.0 / 12.0) * q * q + 3 * q + 1)
                krnls(i)%nx(j) = alpha * 6 * ((2 - q)**5) * (-qx) * ((35.0 / 12.0) * q * q + 3 * q + 1) & 
                           + alpha * ((2 - q)**6) * ((35.0 / 12.0) * 2 * qx + 3 * qx)
                krnls(i)%ny(j) = alpha * 6 * ((2 - q)**5) * (-qy) * ((35.0 / 12.0) * q * q + 3 * q + 1) & 
                           + alpha * ((2 - q)**6) * ((35.0 / 12.0) * 2 * qy + 3 * qy)
                krnls(i)%nz(j) = alpha * 6 * ((2 - q)**5) * (-qz) * ((35.0 / 12.0) * q * q + 3 * q + 1) &
                           + alpha * ((2 - q)**6) * ((35.0 / 12.0) * 2 * qz + 3 * qz)
            else
                krnls(i)%n(j) = 0
                krnls(i)%nx(j) = 0
                krnls(i)%ny(j) = 0
                krnls(i)%nz(j) = 0
            end if
        end do
        !normalizing intrpolants
        sumn = sum(krnls(i)%n)
        krnls(i)%n = krnls(i)%n / sumn
        krnls(i)%nx = krnls(i)%nx / sumn
        krnls(i)%ny = krnls(i)%ny / sumn
        krnls(i)%nz = krnls(i)%nz / sumn
    end do
    !$omp end parallel do
end subroutine get_sph_wndlnd5

subroutine get_sph_inv_dist(h, pnts, totpnts, sctr, krnls)
!-------------------------------------------------------------------------
! inverse distance weight function
! function does not work very well. Maybe should delete...
!-------------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!-------------------------------------------------------------------------
    real(8), intent(in)                       :: h
    integer, intent(in)                       :: totpnts
    real(8), dimension(:,:), intent(in)       :: pnts, sctr
!-------------------------------------------------------------------------
    type(kernel), dimension(:), intent(inout) :: krnls
!-------------------------------------------------------------------------
    integer :: i, j, nbr, totnbrs
    real(8) :: eps, dx, dy, dz, d, sumn
!-------------------------------------------------------------------------

    eps = 1e-6

   !computing intrpolants
   !$omp parallel do private(i) 
   do i = 1, totpnts
        totnbrs = krnls(i)%totnbrs
        allocate(krnls(i)%n(totnbrs))
        allocate(krnls(i)%nx(totnbrs))
        allocate(krnls(i)%ny(totnbrs))
        allocate(krnls(i)%nz(totnbrs))
        do j = 1, totnbrs
            nbr = krnls(i)%nbrs(j)
            dx = pnts(i, 1) - sctr(nbr, 1)
            dy = pnts(i, 2) - sctr(nbr, 2)
            dz = pnts(i, 3) - sctr(nbr, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz) + eps
            krnls(i)%n(j) = ((h / d) - 1)
            krnls(i)%nx(j) = (1 / (d * d)) * krnls(i)%n(j) * dx
            krnls(i)%ny(j) = (1 / (d * d)) * krnls(i)%n(j) * dy
            krnls(i)%nz(j) = (1 / (d * d)) * krnls(i)%n(j) * dz
        end do
        !normalizing intrpolants
        sumn = sum(krnls(i)%n)
        krnls(i)%n = krnls(i)%n / sumn
        krnls(i)%nx = krnls(i)%nx / sumn
        krnls(i)%ny = krnls(i)%ny / sumn
        krnls(i)%nz = krnls(i)%nz / sumn
    end do
    !$omp end parallel do
end subroutine get_sph_inv_dist

end module krnl_sph
