module sdf
implicit none
contains

pure subroutine get_df(sctr, totsctr, pnts, totpnts, df)
!-----------------------------------------------------------------------
! Compute a distance function at 'pnts' relative to a set of 
! scattered points 'sctr'
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: sctr, pnts
    integer, intent(in) :: totsctr, totpnts
!-------------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: df
!-------------------------------------------------------------------------
    integer :: i, j
    real(8) :: dx, dy, dz, d, min_d
!-------------------------------------------------------------------------
    allocate(df(totpnts))
    do i = 1, totpnts
        min_d = 1e6
        do j = 1, totsctr
            dx = pnts(i, 1) - sctr(j, 1)
            dy = pnts(i, 2) - sctr(j, 2)
            dz = pnts(i, 3) - sctr(j, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= min_d) min_d = d
        end do
        df(i) = min_d
    end do
end subroutine get_df

subroutine get_sdf_vec(ip, totip, pnts, totpnts, dim, dxx, sp, phi)
!-------------------------------------------------------------------------
! Constructs an implicit signed distance function by checking the normal
! vector at an interface pointset (ip) compared with the vector from the 
! interface pointset (ip) to the base points set (pnts)
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_wls
    use krnl_mls
    use slvr_prmtrs_struct
    implicit none
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: ip, pnts
    integer, intent(in) :: dim, totip, totpnts
    real(8), intent(in) :: dxx
    type(slvr_prmtrs), intent(in) :: sp       !solver parameters data structure
!-------------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: phi
!-------------------------------------------------------------------------
    type (kernel), dimension(:), allocatable     :: krnls
    !real(8) :: radius, delta
    real(8), dimension(:), allocatable :: df, dfx, dfy, dfz, absdf
    real(8) :: delta
    real(8) :: dx, dy, dz, d
    real(8) :: min_d, min_d_x, min_d_y, min_d_z
    real(8) :: nx, ny, nz
    integer :: i, j, indx, nbr
!-------------------------------------------------------------------------
    
    delta = 10 * dxx !must be a large number

    allocate(dfx(totip))
    allocate(dfy(totip))
    allocate(dfz(totip))
    allocate(absdf(totip))
    allocate(krnls(totpnts))
    
    !compute kernel interpolation functions at interface points
    call get_nbrs_bf(sp%h, ip, totip, pnts, totpnts, krnls)
    select case(sp%krnl)
        case (1)
            call get_rbf_krnls(ip, totip, pnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, sp%rbf_polyex, krnls)
        case (2)
            call get_mls_krnls(ip, totip, pnts, sp%mls, sp%order, sp%h, krnls)
        case (3)
            call get_wls_krnls(ip, totip, pnts, sp%wls, sp%order, sp%h, krnls)
        case default
            call get_sph_krnls(ip, totip, pnts, sp%sph, dim, sp%h,  krnls)
    end select

    !calculate the distance function
    call get_df(ip, totip, pnts, totpnts, df)
    
    !compute the gradients of the distance function at interface points
    do i = 1, totip
        dfx(i) = 0;
        dfy(i) = 0;
        dfz(i) = 0;
        do j = 1, krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            dfx(i) = dfx(i) + krnls(i)%nx(j) * df(nbr)
            dfy(i) = dfy(i) + krnls(i)%ny(j) * df(nbr)
            dfz(i) = dfz(i) + krnls(i)%nz(j) * df(nbr)
        end do
        absdf(i) = sqrt(dfx(i) * dfx(i) + dfy(i) * dfy(i) + dfz(i) * dfz(i)) + 1e-6
    end do

    !computing signed distance function at pointset
    allocate(phi(totpnts))
    do j = 1, totpnts
        !find the closest interface point to point j
        min_d = 1.0e6
        indx = 1;
        do i = 1, totip
            dx = pnts(j,1) - ip(i,1)
            dy = pnts(j,2) - ip(i,2)
            dz = pnts(j,3) - ip(i,3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= min_d) then
                min_d = d
                indx = i
            end if
        end do

        !construct vector between closest interface point and point j
        min_d_x = pnts(j,1) - ip(indx,1)
        min_d_y = pnts(j,2) - ip(indx,2)
        min_d_z = pnts(j,3) - ip(indx,3)
        
        !construct an outer normal vector at closest interface pnt to point j
        nx = ip(indx,1) + delta * (dfx(indx) / absdf(indx))
        ny = ip(indx,2) + delta * (dfy(indx) / absdf(indx))
        nz = ip(indx,3) + delta * (dfz(indx) / absdf(indx))
        
        !the sign of the distance function is the cross product between both vectors
        phi(j) = sign(min_d, min_d_x * nx + min_d_y * ny + min_d_z * nz);
    end do
end subroutine get_sdf_vec

subroutine get_sdf_norm(ip, totip, pnts, totpnts, dim, dxx, sp, phi)
!-------------------------------------------------------------------------
! Constructs an implicit signed distance function by checking the normal
! vector at an interface pointset (ip) compared with the vector from the
! interface pointset (ip) to the base points set (pnts).
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_wls
    use krnl_mls
    use slvr_prmtrs_struct
    implicit none
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: ip, pnts
    integer, intent(in) :: dim, totip, totpnts
    real(8), intent(in) :: dxx
    type(slvr_prmtrs), intent(in) :: sp       !solver parameters data structure
!-------------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: phi
!-------------------------------------------------------------------------
    type (kernel), dimension(:), allocatable     :: krnls
    !real(8) :: radius, delta
    real(8), dimension(:), allocatable :: df, dfx, dfy, dfz, absdf
    real(8) :: delta, s, min_d
    real(8) :: dx, dy, dz, d, d1, d2
    real(8) :: nx_m, ny_m, nz_m, nx_p, ny_p, nz_p
    integer :: i, j, indx, nbr
!-------------------------------------------------------------------------
    !radius = 10*dx !this variable was found critical for successfully constructing the implicit function
    delta = 0.005 * dxx !must be a large number
 
    allocate(dfx(totip))
    allocate(dfy(totip))
    allocate(dfz(totip))
    allocate(absdf(totip))
    allocate(krnls(totpnts))
    
    !compute kernel interpolation functions at interface points
    call get_nbrs_bf(sp%h, ip, totip, pnts, totpnts, krnls)
    select case(sp%krnl)
        case (1)
            call get_rbf_krnls(ip, totip, pnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, sp%rbf_polyex, krnls)
        case (2)
            call get_mls_krnls(ip, totip, pnts, sp%mls, sp%order, sp%h, krnls)
        case (3)
            call get_wls_krnls(ip, totip, pnts, sp%wls, sp%order, sp%h, krnls)
        case default
            call get_sph_krnls(ip, totip, pnts, sp%sph, dim, sp%h,  krnls)
    end select

    !calculate the distance function
    call get_df(ip, totip, pnts, totpnts, df)
    
    !compute the gradients of the distance function at interface points
    do i = 1, totip
        dfx(i) = 0;
        dfy(i) = 0;
        dfz(i) = 0;
        do j = 1, krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            dfx(i) = dfx(i) + krnls(i)%nx(j) * df(nbr)
            dfy(i) = dfy(i) + krnls(i)%ny(j) * df(nbr)
            dfz(i) = dfz(i) + krnls(i)%nz(j) * df(nbr)
        end do
        absdf(i) = sqrt(dfx(i) * dfx(i) + dfy(i) * dfy(i) + dfz(i) * dfz(i)) + 1e-6
    end do

    ! computing signed distance function
    allocate(phi(totpnts))
    do j = 1, totpnts
        !find the closest interface point to point j
        min_d = 1.0e6
        indx = 1;
        do i = 1, totip
            dx = pnts(j,1) - ip(i,1)
            dy = pnts(j,2) - ip(i,2)
            dz = pnts(j,3) - ip(i,3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= min_d) then
                min_d = d
                indx = i
            end if
        end do

        !construct a normal vector on each side of the closest interface point
        nx_p = ip(indx,1) + delta * (dfx(indx) / absdf(indx))
        ny_p = ip(indx,2) + delta * (dfy(indx) / absdf(indx))
        nz_p = ip(indx,3) + delta * (dfz(indx) / absdf(indx))

        nx_m = ip(indx,1) - delta * (dfx(indx) / absdf(indx))
        ny_m = ip(indx,2) - delta * (dfy(indx) / absdf(indx))
        nz_m = ip(indx,3) - delta * (dfz(indx) / absdf(indx))

        !calculate the distance between point j to the both normal vectors
        dx = pnts(j,1) - nx_m
        dy = pnts(j,2) - ny_m
        dz = pnts(j,3) - nz_m
        d1 = sqrt(dx * dx + dy * dy + dz * dz)

        dx = pnts(j,1) - nx_p
        dy = pnts(j,2) - ny_p
        dz = pnts(j,3) - nz_p
        d2 = sqrt(dx * dx + dy * dy + dz * dz)
        
        !the signed distance function at j is that of the closest vector
        s = 0.0
        !if (d1 < d2) then
        !    s = -sign(1.0d, d1 - delta)
        !else
        !    s = -sign(1.0d0, delta - d2)
        !end if

        if (d1 <= d2) then
            s = -1.0d0
        else
            s = 1.0d0
        end if

        phi(j) = s * min_d
    end do
end subroutine get_sdf_norm

subroutine get_sdf_corrected(phi, pnts, totpnts, dxx, ip, totip, phi_corrected)
!-------------------------------------------------------------------------
! Constructs an implicit function by correcting the distance function as
! the interface pointset are advected
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_wls
    use krnl_mls
    use slvr_prmtrs_struct
    implicit none
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: ip, pnts
    real(8), dimension(:), intent(in) :: phi
    real(8), intent(in) :: dxx
    integer, intent(in) :: totip, totpnts
!-------------------------------------------------------------------------
    real(8), dimension(:), allocatable, intent(out) :: phi_corrected
!-------------------------------------------------------------------------
    real(8), dimension(:), allocatable :: phi_ip
    real(8) :: phi_correction, dx, dy, dz, d
    real(8) :: w, sum_w, sum_phi, radius
    type (kernel), dimension(:), allocatable     :: krnls
    integer :: i, j, nbr
!-------------------------------------------------------------------------
    radius = 4 * dxx !best results obtained with 2 to 3 times min_length
    allocate(phi_corrected(totpnts))
    allocate(phi_ip(totip))
    
    !interpolate the level set values from fixed points to interface points
    do i = 1, totip
        sum_w = 0.0
        sum_phi = 0.0
        do j = 1, totpnts
            dx = ip(i,1)-pnts(j,1)
            dy = ip(i,2)-pnts(j,2)
            dz = ip(i,3)-pnts(j,3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= radius) then
                w = 1/d
                !r = 1.1*min_length;
                !q = d/r;
                !w = exp(-pow((q/c),2))-exp(-pow((1/c),2))
                !w= w / (1-exp(-pow((1/c),2)));
                sum_w = sum_w + w
                sum_phi = sum_phi + w * phi(j)
            end if
        end do
        phi_ip(i) = sum_phi / sum_w
    end do

    !but interface points are supposed to have phi = 0. So now we apply the correction
    allocate(krnls(totpnts))
    call get_nbrs_bf(radius, pnts, totpnts, ip, totip, krnls)
   
    do j = 1, totpnts
        phi_correction = 0.0
        sum_w = 0.0
        sum_phi = 0.0
        if (krnls(i)%totnbrs > 0) then
            do i = 1, krnls(i)%totnbrs
                nbr = krnls(j)%nbrs(i)
                dx = pnts(j,1) - ip(nbr,1)
                dy = pnts(j,2) - ip(nbr,2)
                dz = pnts(j,3) - ip(nbr,3)
                d = sqrt(dx * dx + dy * dy + dz * dz) + 1e-6
                w = 1 / d;
                sum_w = sum_w + w
                sum_phi = sum_phi + w * phi_ip(nbr)
            end do
            phi_correction = sum_phi / sum_w
        end if
        phi_corrected(j) = phi(j) - phi_correction
    end do
end subroutine get_sdf_corrected

subroutine get_sdf_grad(ip, totip, phi, pnts, totpnts, krnls, vdt, vnt)
!-------------------------------------------------------------------------
! Constructs an implicit function by using the gradient
! of the distance function from the interface pointset as a driving force
!-------------------------------------------------------------------------
    use intrf
    use krnl_struct
    implicit none
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: ip, pnts
    real(8), intent(in) :: vdt
    integer, intent(in) :: totip, totpnts, vnt
    type(kernel), dimension(:), intent(in) :: krnls
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(inout) :: phi
!-------------------------------------------------------------------------
    integer :: i, j, t, nbr
    real(8) :: dfx, dfy, dfz, phi_x, phi_y, phi_z, absdphi
    real(8), dimension(:), allocatable :: vn, df
    real(8), dimension(:,:), allocatable :: vel
!-------------------------------------------------------------------------

    allocate(vn(totpnts))
    allocate(vel(totpnts,3))
    vel(:,:) = 0.0

    call get_df(ip, totip, pnts, totpnts, df)

    do i = 1, totpnts
        dfx = 0.0
        dfy = 0.0
        dfz = 0.0
        phi_x = 0.0
        phi_y = 0.0
        phi_z = 0.0
        do j = 1, krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            dfx = dfx + krnls(i)%nx(j) * df(nbr)
            dfy = dfy + krnls(i)%ny(j) * df(nbr)
            dfz = dfz + krnls(i)%nz(j) * df(nbr)

            phi_x = phi_x + krnls(i)%nx(j) * phi(nbr)
            phi_y = phi_y + krnls(i)%ny(j) * phi(nbr)
            phi_z = phi_z + krnls(i)%nz(j) * phi(nbr)
        end do
        absdphi = sqrt(phi_x * phi_x + phi_y * phi_y + phi_z * phi_z) + 1e-6
        vn(i) = -(dfx * phi_x + dfy * phi_y + dfz * phi_z) / absdphi
        !vel(i,1) = dfx
        !vel(i,2) = dfy
        !vel(i,3) = dfz
    end do

    do t = 1, vnt
        call vof(phi, vel, vn, krnls, totpnts, vdt, 1)
    end do

end subroutine get_sdf_grad

end module sdf
