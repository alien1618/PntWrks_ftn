module pntst
implicit none
contains

subroutine get_avg_dx(pnts, avg_dx)
!-------------------------------------------------------------------------
    use omp_lib
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
!-------------------------------------------------------------------------
    real(8), intent(out) :: avg_dx
!-------------------------------------------------------------------------
    real(8) :: d_min, dminf = 1e6, dx, dy, dz, d
    integer :: i, j
    real(8) :: sum_d = 0
    integer :: totsamples
!-------------------------------------------------------------------------
    totsamples = int(0.2*size(pnts))
    do i=1, totsamples
        d_min = 1e6
        do j = 1, TotSamples
            dx = pnts(i, 1) - pnts(j, 1)
            dy = pnts(i, 2) - pnts(j, 2)
            dz = pnts(i, 3) - pnts(j, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= d_min .and. d > 1.0e-6) d_min = d
        end do
        if (d_min <= dminf) dminf = d_min
        sum_d = sum_d + d_min
    end do
    !avg_dx = sum_d / TotSamples;
    avg_dx = dminf;
    write(*,'(a,f10.3)') "Average dx = ", avg_dx
end subroutine get_avg_dx

subroutine gen_sctr_pnts(p, l, totsctr, pnts)
!-------------------------------------------------------------------------
! Generate scattered points in a specified region
!-------------------------------------------------------------------------
    real(8), dimension(3), intent(in) :: p, l
    integer, intent(in) :: totsctr
!-------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: pnts
!-------------------------------------------------------------------------
    real(8) :: sctrx, sctry, dx, dy, d, num
    real(8), dimension(:), allocatable :: sx, sy
    integer :: tot, i, j, cntr
!-------------------------------------------------------------------------

    allocate(pnts(totsctr,3))
    pnts(:,:) = 0.0

   allocate(sx(0))
   allocate(sy(0))

   call random_number(num)
   sctrx = p(1) + num * (l(1))
   call random_number(num)
   sctry = p(2) + num * (l(2))
   sx = [sx, sctrx]
   sy = [sy, sctry]
   tot = 1
   do i = 1, 100*totsctr
        call random_number(num)
        sctrx = p(1) + num * (l(1))
        call random_number(num)
        sctry = p(2) + num * (l(2))
        cntr = 0
        do j = 1,tot
            dx = sx(j) - sctrx
            dy = sy(j) - sctry
            d = dx * dx + dy * dy
            if (d <= 0.0001) cntr = cntr + 1
        end do
        if (cntr == 0) then
            sx = [sx, sctrx]
            sy = [sy, sctry]
            tot = tot + 1
        end if
        if (tot >= totsctr) exit
    end do 
    write (*,*) "totsctr = ", tot
    pnts(:, 1) = sx
    pnts(:, 2) = sy

end subroutine gen_sctr_pnts

pure subroutine gen_grdpnts1d(p0, l0, n, ps)
!-------------------------------------------------------------------------
! Generate 1D grid of points
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    integer, dimension(3), intent(in)  ::  n
    real(8), dimension(3), intent(in)  ::  p0, l0
!-------------------------------------------------------------------------
    type(pointset), intent(out)        :: ps
!-------------------------------------------------------------------------
    integer                            ::  i, l
    real(8)                            ::  dx
!-------------------------------------------------------------------------
    dx = l0(1) / (n(1) - 1)
    ps%totpnts = n(1)
    ps%dim = 1
    ps%order = 1
    ps%surfshape = 0
    ps%surfpnts= 0
    ps%totsurfs = 0
    ps%totvols = 0

    allocate(ps%pnts(ps%totpnts + 1, 3))   
    ps%pnts(1, 1) = p0(1)
    ps%pnts(1, 2) = 0.0
    ps%pnts(1, 3) = 0.0
    l = 1;
    do i = 1, n(1)
        ps%pnts(l + 1, 1) = ps%pnts(l, 1) + dx
        ps%pnts(l + 1, 2) = 0.0
        ps%pnts(l + 1, 3) = 0
        l = l + 1
    end do
end subroutine gen_grdpnts1d

pure subroutine gen_quadgrd(p0, l0, n, ps)
!-------------------------------------------------------------------------
! Generate 2D quadrilateral grid
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    integer, dimension(3), intent(in)  ::  n
    real(8), dimension(3), intent(in)  ::  p0, l0
!-------------------------------------------------------------------------
    type(pointset), intent(out)        :: ps
!-------------------------------------------------------------------------
    integer                            ::  elems_x, elems_y, i, j, k, l
    integer, dimension(:), allocatable ::  nd
    real(8)                            ::  dx, dy
!-------------------------------------------------------------------------
    dx = l0(1) / (n(1) - 1)
    dy = l0(2) / (n(2) - 1)
    ps%totpnts = n(1) * n(2);
    elems_x = n(1) - 1;
    elems_y = n(2) - 1;
    ps%dim = 2;
    ps%order = 1;
    ps%surfshape = 2;
    ps%surfpnts=4;
    ps%totsurfs = elems_x * elems_y

    allocate(ps%pnts(ps%totpnts + 1, 3))
    allocate(nd(ps%totpnts))
    allocate(ps%surfs(ps%totsurfs + 1, ps%surfpnts))
    
    ps%pnts(1, 1) = p0(1)
    ps%pnts(1, 2) = p0(2)
    ps%pnts(1, 3) = 0
    l = 1;
    do j = 1, n(2)
        do i = 1, n(1)
            ps%pnts(l + 1, 1) = ps%pnts(l, 1) + dx
            ps%pnts(l + 1, 2) = ps%pnts(l, 2)
            ps%pnts(l + 1, 3) = 0
            l = l + 1
        end do
        ps%pnts(l, 1) = ps%pnts(1, 1)
        ps%pnts(l, 2) = ps%pnts(l, 2) + dy
        ps%pnts(l, 3) = 0
    end do

    do i=1,ps%totpnts
       nd(i) = i
    end do
    k = 1;
    l = 1;
    do j = 1, elems_y
        do i = 1, elems_x
            ps%surfs(l, 1) = nd(k)
            ps%surfs(l, 2) = nd(k + 1)
            ps%surfs(l, 3) = nd(k + elems_x + 2)
            ps%surfs(l, 4) = nd(k + elems_x + 1)
            k = k + 1
            l = l + 1
        end do
        k = k+1
    end do
    !write(*,1) "totpnts = ", ps%totpnts
    !write(*,1) "totelems = ", ps%totsurfs
    !write(*,1) "elempnts = ", ps%surfpnts
    !write(*,1) "order = ", ps%order
    !write(*,1) "elemshape = ", ps%surfshape
    !1 format(a,i0)
end subroutine gen_quadgrd

pure subroutine gen_trigrd(p0,l0,n,msh)
!-------------------------------------------------------------------------
! Generate 2D triangular grid
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    integer, dimension(3), intent(in)  :: n
    real(8),dimension(3), intent(in)   :: p0, l0
!-------------------------------------------------------------------------
    type(pointset), intent(out)        :: msh
!-------------------------------------------------------------------------
    integer                            :: elems_x, elems_y, i, j, k, l
    integer,dimension(:), allocatable  :: nd
    real(8)                            :: dx, dy
!-------------------------------------------------------------------------

    dx = l0(1) / (n(1) - 1)
    dy = l0(2) / (n(2) - 1)
    elems_x = n(1) - 1
    elems_y = n(2) - 1
    msh%order = 1
    msh%surfshape = 1
    msh%surfpnts = 3
    msh%totsurfs = 2 * elems_x * elems_y
    msh%totvols = 0
    allocate(msh%surfs(msh%totsurfs, msh%surfpnts))

    msh%totpnts = (elems_x + 1) * (elems_y + 1)
    allocate(msh%pnts(msh%totpnts+1, 3))
    msh%pnts(1,1) = p0(1)
    msh%pnts(1,2) = p0(2)
    msh%pnts(1,3) = 0
    l = 1;
    do j=1,n(2)
        do i=1,n(1)
            msh%pnts(l + 1, 1) = msh%pnts(l, 1) + dx
            msh%pnts(l + 1, 2) = msh%pnts(l, 2)
            msh%pnts(l + 1, 3) = 0
            l = l+1;
        end do
        msh%pnts(l, 1) = msh%pnts(1, 1)
        msh%pnts(l, 2) = msh%pnts(l, 2) + dy
        msh%pnts(l, 3) = 0
    end do

    allocate(nd(msh%totpnts))
    do i = 1, msh%totpnts
        nd(i) = i
    end do

    k = 1
    l = 1
    do j = 1, elems_y
        do i = 1, elems_x
            if ((mod(i,2) == 0 .and. mod(j,2) /= 0) .or. (mod(i,2) /= 0 .and. mod(j,2) == 0)) then
                msh%surfs(l,1) = nd(k)
                msh%surfs(l,2) = nd(k+1)
                msh%surfs(l,3) = nd(k+elems_x+1)

                msh%surfs(l+1,1) = nd(k+1)
                msh%surfs(l+1,2) = nd(k+elems_x+2)
                msh%surfs(l+1,3) = nd(k+elems_x+1)
            else  if ((mod(i,2) /= 0 .and. mod(j,2) /= 0) .or. (mod(i,2) == 0 .and. mod(j,2)== 0)) then
                msh%surfs(l,1) = nd(k)
                msh%surfs(l,2) = nd(k+1)
                msh%surfs(l,3) = nd(k+elems_x+2)

                msh%surfs(l+1,1) = nd(k)
                msh%surfs(l+1,2) = nd(k+elems_x+2)
                msh%surfs(l+1,3) = nd(k+elems_x+1)
            end if
            k = k+1
            l = l+2
        end do
        k = k+1
    end do
end subroutine

pure subroutine gen_trigrd2(p0,l0,n,msh)
!-------------------------------------------------------------------------
! Generate 2D triangular grid version 2
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    integer, dimension(3), intent(in)   :: n
    real(8), dimension(3), intent(in)   :: p0, l0
!-------------------------------------------------------------------------
    type(pointset), intent(out)         :: msh
!-------------------------------------------------------------------------
    integer                             :: elems_x, elems_y, i, j, k, l
    integer,dimension(:), allocatable    :: nd
    real(8)                             :: dx, dy
!-------------------------------------------------------------------------

    dx = l0(1) / (n(1) - 1)
    dy = l0(2) / (n(2) - 1)
    elems_x = n(1) - 1
    elems_y = n(2) - 1
    msh%order = 1
    msh%surfshape = 1
    msh%surfpnts = 3
    msh%totsurfs = 2 * elems_x * elems_y
    msh%totvols = 0
    allocate(msh%surfs(msh%totsurfs, msh%surfpnts))

    msh%totpnts = (elems_x + 1) * (elems_y + 1)
    allocate(msh%pnts(msh%totpnts + 1, 3))
    msh%pnts(1,1) = p0(1)
    msh%pnts(1,2) = p0(2)
    msh%pnts(1,3) = 0
    l = 1;
    do j=1,n(2)
        do i=1,n(1)
            msh%pnts(l+1,1) = msh%pnts(l,1) + dx;
            msh%pnts(l+1,2) = msh%pnts(l,2);
            msh%pnts(l+1,3) = 0
            l = l+1;
        end do
        msh%pnts(l,1) = msh%pnts(1,1)
        msh%pnts(l,2) = msh%pnts(l,2) + dy
        msh%pnts(l,3) = 0
    end do

    allocate(nd(msh%totpnts))
    do i = 1, msh%totpnts
        nd(i) = i
    end do

    k = 1
    l = 1
    do j = 1, elems_y
        do i = 1, elems_x
            msh%surfs(l,1) = nd(k)
            msh%surfs(l,2) = nd(k+1)
            msh%surfs(l,3) = nd(k+elems_x+1)

            msh%surfs(l+1,1) = nd(k+1)
            msh%surfs(l+1,2) = nd(k+elems_x+2)
            msh%surfs(l+1,3) = nd(k+elems_x+1)
            l = l+2;
            k = k+1;
        end do
        k = k+1;
    end do
end subroutine gen_trigrd2

pure subroutine gen_grdpnts3d(p0, l0, n, ps)
!-------------------------------------------------------------------------
! Generate 3D grid points
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    integer, dimension(3), intent(in) ::  n
    real(8), dimension(3), intent(in) ::  p0, l0
!-------------------------------------------------------------------------
    type(pointset), intent(out)       ::  ps
!-------------------------------------------------------------------------
    integer                           ::  i, j, k, l
    real(8)                           ::  dx, dy, dz
!-------------------------------------------------------------------------
    dx = l0(1) / (n(1) - 1)
    dy = l0(2) / (n(2) - 1)
    dz = l0(3) / (n(3) - 1)
    ps%totpnts = n(1) * n(2) * n(3);
    allocate(ps%pnts(ps%totpnts + 1, 3))
    ps%dim =  3
    ps%pnts(1,1) = p0(1)
    ps%pnts(1,2) = p0(2)
    ps%pnts(1,3) = p0(3)
    l = 1;
    do k = 1, n(3)
        do j = 1, n(2)
            do i = 1, n(1)
                ps%pnts(l + 1, 1) = ps%pnts(l, 1) + dx;
                ps%pnts(l + 1, 2) = ps%pnts(l, 2);
                ps%pnts(l + 1, 3) = ps%pnts(l, 3)
                l = l+1;
            end do
            ps%pnts(l, 1) = ps%pnts(1, 1)
            ps%pnts(l, 2) = ps%pnts(l, 2) + dy
            ps%pnts(l, 3) = ps%pnts(l, 3)
        end do
        ps%pnts(l, 1) = p0(1)
        ps%pnts(l, 2) = p0(2)
        ps%pnts(l, 3) = ps%pnts(l, 3) + dz
    end do
    ps%totsurfs = 0
    ps%totvols = 0
    !write(*,'(a,i0)') "totpnts = ", ps%totpnts
end subroutine gen_grdpnts3d

pure subroutine gen_hexgrd(p0,l0,n,ps)
!-------------------------------------------------------------------------
! Generate 3D hexahedral grid
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    real(8), dimension(3), intent(in)      :: p0, l0
    integer, dimension(3), intent(in)      ::  n
!-------------------------------------------------------------------------
    type(pointset), intent(out)            ::  ps
!-------------------------------------------------------------------------
    integer, dimension(:,:,:), allocatable ::  ndnum
    integer                                ::  i, j, k, l, e, s
!-------------------------------------------------------------------------

    call gen_grdpnts3d(p0, l0, n, ps)
    ps%dim = 3
    ps%order = 1
    ps%volshape = 2
    ps%surfshape = 2
    ps%volpnts = 8
    ps%surfpnts = 4

    ps%totvols = (n(1) - 1) * (n(2) - 1) * (n(3) - 1)
    ps%totsurfs = ps%totvols * 6

    allocate(ndnum(n(1), n(2), n(3)))
    allocate(ps%surfs(ps%totsurfs, ps%surfpnts))
    allocate(ps%vols(ps%totvols, ps%volpnts))

    l = 1
    do k = 1, n(3)
        do j = 1, n(2)
            do i = 1, n(1)
                ndnum(i, j, k) = l
                l = l + 1
            end do
        end do
    end do

    e = 1
    do k = 1, n(3) - 1
        do j = 1, n(2) - 1
            do i = 1, n(1) - 1
                ps%vols(e, 1) = ndnum(i, j, k)
                ps%vols(e, 2) = ndnum(i + 1, j, k)
                ps%vols(e, 3) = ndnum(i + 1, j + 1, k)
                ps%vols(e, 4) = ndnum(i, j + 1, k)
                ps%vols(e, 5) = ndnum(i, j, k + 1)
                ps%vols(e, 6) = ndnum(i + 1, j, k + 1)
                ps%vols(e, 7) = ndnum(i + 1, j + 1, k+1)
                ps%vols(e, 8) = ndnum(i, j + 1, k + 1)
                e = e + 1
            end do
        end do
    end do

    s = 1
    do i = 1,ps%totvols
        ps%surfs(s,1) = ps%vols(i,1)
        ps%surfs(s,2) = ps%vols(i,2)
        ps%surfs(s,3) = ps%vols(i,3)
        ps%surfs(s,4) = ps%vols(i,4)

        ps%surfs(s+1,1) = ps%vols(i,5)
        ps%surfs(s+1,2) = ps%vols(i,6)
        ps%surfs(s+1,3) = ps%vols(i,7)
        ps%surfs(s+1,4) = ps%vols(i,8)

        ps%surfs(s+2,1) = ps%vols(i,2)
        ps%surfs(s+2,2) = ps%vols(i,6)
        ps%surfs(s+2,3) = ps%vols(i,7)
        ps%surfs(s+2,4) = ps%vols(i,3)

        ps%surfs(s+3,1) = ps%vols(i,1)
        ps%surfs(s+3,2) = ps%vols(i,5)
        ps%surfs(s+3,3) = ps%vols(i,8)
        ps%surfs(s+3,4) = ps%vols(i,4)

        ps%surfs(s+4,1) = ps%vols(i,1)
        ps%surfs(s+4,2) = ps%vols(i,2)
        ps%surfs(s+4,3) = ps%vols(i,6)
        ps%surfs(s+4,4) = ps%vols(i,5)

        ps%surfs(s+5,1) = ps%vols(i,4)
        ps%surfs(s+5,2) = ps%vols(i,3)
        ps%surfs(s+5,3) = ps%vols(i,7)
        ps%surfs(s+5,4) = ps%vols(i,8)

        s = s+6
    end do
end subroutine gen_hexgrd

pure subroutine gen_tetgrd(p0,l0,n,msh)
!-------------------------------------------------------------------------
! Generate 3D tetrahedral grid
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    integer, dimension(3), intent(in)      :: n
    real(8), dimension(3), intent(in)      :: p0, l0
!-------------------------------------------------------------------------
    type(pointset), intent(out)            :: msh
!-------------------------------------------------------------------------
    integer,dimension(:,:,:), allocatable  :: ndnum
    integer                                :: i, j, k, l, e, s
    integer                                :: tets_per_cube
!-------------------------------------------------------------------------

    call gen_grdpnts3d(p0,l0,n,msh)
    
    msh%order = 1
    msh%volshape = 1
    msh%volpnts = 4
    msh%surfshape = 1
    msh%surfpnts = 3
    tets_per_cube = 5
    select case (tets_per_cube)
        case (5)
            msh%totvols =  5 * (n(1) - 1) * (n(2) - 1) * (n(3) - 1)
        case (6)
            msh%totvols =  6 * (n(1) - 1) * (n(2) - 1) * (n(3) - 1)
    end select
 
    msh%totsurfs = msh%totvols * 4;

    allocate(ndnum(n(1), n(2), n(3)))
    allocate(msh%surfs(msh%totsurfs, msh%surfpnts))
    allocate(msh%vols(msh%totvols, msh%volpnts))

    l = 1
    do k = 1, n(3)
        do j = 1, n(2)
            do i = 1, n(1)
                ndnum(i, j, k) = l
                l = l + 1
            end do
        end do
    end do
    
    e = 1
    select case (tets_per_cube)
        case (5)
            do k = 1, n(3) - 1
                do j = 1, n(2) - 1
                    do i = 1, n(1) - 1
                        msh%vols(e,1) = ndnum(i,j,k)
                        msh%vols(e,2) = ndnum(i+1,j,k)
                        msh%vols(e,3) = ndnum(i,j+1,k)
                        msh%vols(e,4) = ndnum(i,j,k+1)

                        msh%vols(e+1,1) = ndnum(i+1,j,k)
                        msh%vols(e+1,2) = ndnum(i+1,j+1,k)
                        msh%vols(e+1,3) = ndnum(i,j+1,k)
                        msh%vols(e+1,4) = ndnum(i+1,j+1,k+1)

                        msh%vols(e+2,1) = ndnum(i,j,k+1)
                        msh%vols(e+2,2) = ndnum(i,j+1,k)
                        msh%vols(e+2,3) = ndnum(i,j+1,k+1)
                        msh%vols(e+2,4) = ndnum(i+1,j+1,k+1)

                        msh%vols(e+3,1) = ndnum(i,j,k+1)
                        msh%vols(e+3,2) = ndnum(i+1,j,k)
                        msh%vols(e+3,3) = ndnum(i+1,j+1,k+1)
                        msh%vols(e+3,4) = ndnum(i+1,j,k+1)

                        msh%vols(e+4,1) = ndnum(i,j+1,k)
                        msh%vols(e+4,2) = ndnum(i+1,j,k)
                        msh%vols(e+4,3) = ndnum(i,j,k+1)
                        msh%vols(e+4,4) = ndnum(i+1,j+1,k+1)

                        e = e + 5
                    end do
                end do
            end do
        case (6)
            do k = 1, n(3) - 1
                do j = 1, n(2) - 1
                    do i = 1, n(1) - 1 
                        msh%vols(e,1) = ndnum(i,j,k);
                        msh%vols(e,2) = ndnum(i+1,j,k);
                        msh%vols(e,3) = ndnum(i+1,j+1,k);
                        msh%vols(e,4) = ndnum(i+1,j,k+1);

                        msh%vols(e+1,1) = ndnum(i,j,k);
                        msh%vols(e+1,2) = ndnum(i+1,j,k+1);
                        msh%vols(e+1,3) = ndnum(i,j,k+1);
                        msh%vols(e+1,4) = ndnum(i,j+1,k+1);

                        msh%vols(e+2,1) = ndnum(i,j+1,k);
                        msh%vols(e+2,2) = ndnum(i+1,j+1,k);
                        msh%vols(e+2,3) = ndnum(i+1,j+1,k+1);
                        msh%vols(e+2,4) = ndnum(i+1,j,k+1);

                        msh%vols(e+3,1) = ndnum(i,j,k);
                        msh%vols(e+3,2) = ndnum(i+1,j+1,k);
                        msh%vols(e+3,3) = ndnum(i,j+1,k);
                        msh%vols(e+3,4) = ndnum(i+1,j,k+1);

                        msh%vols(e+4,1) = ndnum(i,j,k);
                        msh%vols(e+4,2) = ndnum(i,j+1,k);
                        msh%vols(e+4,3) = ndnum(i,j+1,k+1);
                        msh%vols(e+4,4) = ndnum(i+1,j,k+1);

                        msh%vols(e+5,1) = ndnum(i,j+1,k);
                        msh%vols(e+5,2) = ndnum(i+1,j+1,k+1);
                        msh%vols(e+5,3) = ndnum(i,j+1,k+1);
                        msh%vols(e+5,4) = ndnum(i+1,j,k+1);
                        e = e + 6

                    end do
                end do
            end do
    end select

    s = 1
    do i = 1, msh%totvols
        msh%surfs(s,1) = msh%vols(i,1)
        msh%surfs(s,2) = msh%vols(i,2)
        msh%surfs(s,3) = msh%vols(i,3)

        msh%surfs(s+1,1) = msh%vols(i,1)
        msh%surfs(s+1,2) = msh%vols(i,2)
        msh%surfs(s+1,3) = msh%vols(i,4)

        msh%surfs(s+2,1) = msh%vols(i,2)
        msh%surfs(s+2,2) = msh%vols(i,3)
        msh%surfs(s+2,3) = msh%vols(i,4)

        msh%surfs(s+3,1) = msh%vols(i,1)
        msh%surfs(s+3,2) = msh%vols(i,3)
        msh%surfs(s+3,3) = msh%vols(i,4)
        s = s + 4
    end do
end subroutine gen_tetgrd

subroutine prnt_vtk(u, ps, fname, t)
!-------------------------------------------------------------------------
! Print pointset in VTK format
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    type(pointset), intent(in)        :: ps
    integer, intent(in)               ::  t
    real(8),dimension(:), intent(in)  ::  u
    character(len=50), intent(in)     :: fname
!-------------------------------------------------------------------------
    integer                           ::  i, j, sum, elemorder
    character(len=50) :: file_name
    character(len=10) :: t_s
    character*29      :: path = 'sim/out/'
!-------------------------------------------------------------------------
        write(t_s, '(i0)') t
        call system('mkdir -p ' // adjustl(trim(path)) // trim(adjustl(fname)) // "_vtk/")
        file_name = trim(adjustl(path)) // trim(adjustl(fname)) // "_vtk/" // trim(adjustl(fname))  // trim(adjustl(t_s)) // ".vtk"
    
        elemorder = 1
        if (ps%surfpnts == 6 .or. ps%surfpnts == 8) elemorder = 2
        sum = ((ps%surfpnts/elemorder)+1)*ps%totsurfs;

        write(*,'(a)') "Printing results to: " // file_name
        open(UNIT=1, FILE=file_name)
        write(1,'(a)') "# vtk DataFile Version 1.0"
        write(1,'(a)') "Cube example"
        write(1,'(a)') "ASCII"
        write(1,*)
        write(1,'(a)') "DATASET POLYDATA"
        write(1,1) "POINTS ",ps%totpnts," double"
        1 format(A,i0,A)
        
        !do i = 1,totpnts
        !  do j = 1,2
        !   write(1,'(2X,F12.6)',advance='no') pnts(i,j)      
        !  end do
        !   write(1,'(2X,F12.6)',advance='no') 0.0
        !   write(1,*)
        !end do
        !write(1,*)
        
        do i = 1,ps%totpnts
            do j = 1,3
                write(1,'(2X,F12.6)',advance='no') ps%pnts(i,j)      
            end do
            write(1,*)
        end do
        write(1,*) 
        if (ps%totsurfs > 0) then
            write(1,2) "POLYGONS ",ps%totsurfs," ",sum
            2 format(A,i0,A,i0)
            do i = 1,ps%totsurfs
                write(1,'(i0)',advance='no') ps%surfpnts/elemorder
                do j = 1,ps%surfpnts/elemorder
                    write(1,'(4X,i0)',advance='no') ps%surfs(i,j) - 1
                end do
                write(1,*)
            end do
            write(1,*)
        end if
        write(1,3) "POINT_DATA ",ps%totpnts
        3 format(A,i0)
        write(1,'(a)') "SCALARS myscalars double"
        write(1,'(a)') "LOOKUP_TABLE custom_table"
        do i=1,ps%totpnts
            write(1,'(F20.6)',advance='yes') u(i)
        end do
        close(unit=1)

end subroutine prnt_vtk

subroutine prnt_unstruct_vtk(u, ps, fname, t)
    use pntst_struct
    type(pointset), intent(in)              :: ps
    integer, intent(IN)                 ::  t
    real(8),dimension(:), intent(IN)       ::  u
    integer                             ::  i, j, sum
    
    character(len=50), intent(IN)::fname
    character(len=50) :: file_name
    character(len=10) :: t_s
    character*29 :: path='sim/out/'
    call system('mkdir -p ' // adjustl(trim(path)) // trim(adjustl(fname)) // "_vtk/")
    !Construct the filename:
    write(t_s, '(i0)') t
    file_name = trim(adjustl(path)) // trim(adjustl(fname)) // "_vtk/" // trim(adjustl(fname))  // trim(adjustl(t_s)) // ".vtk"
   
    sum = (ps%volpnts)*ps%totvols;

    write(*,'(a)') "Printing results to: " // file_name
    open(UNIT=1, FILE=file_name)
    write(1,'(a)') "# vtk DataFile Version 1.0"
    write(1,'(a)') "Cube example"
    write(1,'(a)') "ASCII"
    write(1,*)
    write(1,'(a)') "DATASET UNSTRUCTURED_GRID"
    write(1,1) "POINTS ",ps%totpnts," double"
    1 format(A,i0,A)

    do i = 1,ps%totpnts
      do j = 1,3
       write(1,'(2X,F12.6)',advance='no') ps%pnts(i,j)      
      end do
       write(1,*)
    end do
    write(1,*) 
    write(1,2) "CELLS ",ps%totvols," ", (sum+ps%totvols)
    2 format(A,i0,A,i0)
    do i = 1,ps%totvols
        write(1,'(i0)',advance='no') ps%volpnts
        do j = 1,ps%volpnts
            write(1,'(4X,i0)',advance='no') ps%vols(i,j)-1
        end do
        write(1,*)
    end do
    write(1,*)
    write(1,3) "CELL_TYPES ",ps%totvols
    do i=1,ps%totvols
        write(1,'(i0)') 12
    end do
    write(1,*)
    write(1,3) "POINT_DATA ",ps%totpnts
    3 format(A,i0)
    write(1,'(a)') "SCALARS myscalars double"
    write(1,'(a)') "LOOKUP_TABLE default"
    do i=1,ps%totpnts
        write(1,'(F12.6)',advance='yes') u(i)
    end do
    close(unit=1)
end subroutine prnt_unstruct_vtk

subroutine prnt_vec_vtk(vel, msh, fname, t)
!-------------------------------------------------------------------------
! Print vectors in VTK format
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    type(pointset), intent(in)          :: msh                  
    integer, intent(in)                 :: t
    real(8),dimension(:,:), intent(in)  :: vel
    character(len=50), intent(in)       :: fname
!-------------------------------------------------------------------------
    integer                             :: i, j, sum
    character(len=50)                   :: file_name
    character(len=10)                   :: t_s
    character*29                        :: path = 'sim/out/'
    real(8), dimension(msh%totpnts)     :: v
!-------------------------------------------------------------------------
    write(t_s, '(i0)') t
    call system('mkdir -p ' // adjustl(trim(path)) // trim(adjustl(fname)) // "_vtk/")
    file_name = trim(adjustl(path)) // trim(adjustl(fname)) // "_vtk/" // trim(adjustl(fname))  // trim(adjustl(t_s)) // ".vtk"
   
    sum = (msh%surfpnts+1)*msh%totsurfs

    v = sqrt(vel(:,1) * vel(:,1) + vel(:,2) * vel(:,2) + vel(:,3) * vel(:,3))

    write(*,'(a)') "Printing results to: " // file_name
    open(UNIT=1, FILE=file_name)
    write(1,'(a)') "# vtk DataFile Version 1.0"
    write(1,'(a)') "Cube example"
    write(1,'(a)') "ASCII"
    write(1,*)
    write(1,'(a)') "DATASET POLYDATA"
    write(1,1) "POINTS ",msh%totpnts," double"
    1 format(A,i0,A)

    do i = 1,msh%totpnts
        do j = 1,3
            write(1,'(2X,F12.6)',advance='no') msh%pnts(i,j)      
        end do
       write(1,*)
    end do
    write(1,*) 
    
    write(1,2) "POLYGONS ",msh%totsurfs," ",sum
    2 format(A,i0,A,i0)
    do i = 1,msh%totsurfs
        write(1,'(i0)',advance='no') msh%surfpnts
        do j = 1,msh%surfpnts
            write(1,'(4X,i0)',advance='no') msh%surfs(i,j) - 1
        end do
        write(1,*)
    end do
    write(1,*)

    write(1,3) "POINT_DATA ",msh%totpnts
    3 format(A,i0)
    write(1,'(a)') "SCALARS myscalars double"
    write(1,'(a)') "LOOKUP_TABLE custom_table"
    do i=1,msh%totpnts
        write(1,'(F18.6)',advance='yes') v(i)
    end do
    write(1,*)

    write(1,'(a)') "vec vec double"
    do i = 1, msh%totpnts
        write(1,'(F18.6, a, F18.6, a, F18.6)') vel(i,1), " ", vel(i,2), " ", vel(i,3)
    end do

    close(unit=1)
end subroutine prnt_vec_vtk

subroutine prnt_txt(u, pnts, totpnts, fname, t)
!-------------------------------------------------------------------------
! Print points in TXT format
!-------------------------------------------------------------------------
    integer, intent(in)                 :: totpnts, t
    real(8), dimension(:), intent(in)   :: u
    real(8), dimension(:,:), intent(in) :: pnts
    character(len=50), intent(in)       :: fname
!-------------------------------------------------------------------------
    integer                             :: i, j
    character(len=50)                   :: file_name
    character(len=10)                   :: t_s
    character*29                        :: path = 'sim/out/'
!-------------------------------------------------------------------------
    write(t_s, '(i0)') t
    call system('mkdir -p ' // adjustl(trim(path)) // trim(adjustl(fname)) // "_txt/")
    file_name = trim(adjustl(path)) // trim(adjustl(fname)) // "_txt/" // trim(adjustl(fname))  // trim(adjustl(t_s)) // ".txt"

    write(*,'(a)') "Printing results to: " // file_name
    open(UNIT=1, FILE=file_name)
    do i = 1,totpnts
        do j = 1,3
            write(1,'(2X,F12.6)',advance='no') pnts(i,j)      
        end do
        write(1,'(F12.6)',advance='no') u(i)
        write(1,*)
    end do
    write(1,*) 
    close(unit=1)
end subroutine prnt_txt

subroutine prnt_pnts_vtk(u, pnts, totpnts, fname, t)
!-------------------------------------------------------------------------
! Print pointset in VTK format
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in)                 :: totpnts
    integer, intent(in)                 :: t
    real(8),dimension(:), intent(in)    :: u
!-------------------------------------------------------------------------
    integer                             :: i, j
    character(len=50), intent(in)       :: fname
    character(len=50)                   :: file_name
    character(len=10)                   :: t_s
    character*29                        :: path='sim/out/'
!-------------------------------------------------------------------------
    write(t_s, '(i0)') t
    call system('mkdir -p ' // adjustl(trim(path)) // trim(adjustl(fname)) // "_vtk/")
    file_name = trim(adjustl(path)) // trim(adjustl(fname)) // "_vtk/" // trim(adjustl(fname))  // trim(adjustl(t_s)) // ".vtk"

    write(*,'(a)') "Printing results to: " // file_name
    open(UNIT=1, FILE=file_name)
    write(1,'(a)') "# vtk DataFile Version 1.0"
    write(1,'(a)') "Cube example"
    write(1,'(a)') "ASCII"
    write(1,*)
    write(1,'(a)') "DATASET POLYDATA"
    write(1,1) "POINTS ",totpnts," double"
    1 format(A,i0,A)

    do i = 1,totpnts
        do j = 1,3
            write(1,'(2X,F12.6)',advance='no') pnts(i,j)      
        end do
        write(1,*)
    end do
    write(1,*) 
    write(1,3) "POINT_DATA ",totpnts
    3 format(A,i0)
    write(1,'(a)') "SCALARS myscalars double"
    write(1,'(a)') "LOOKUP_TABLE custom_table"
    do i=1,totpnts
        write(1,'(F20.6)',advance='yes') u(i)
    end do
    close(unit=1)

end subroutine prnt_pnts_vtk

subroutine prnt_pks_vtk(u, pnts, totpnts, fname, t)
!-------------------------------------------------------------------------
! Print u as peak points
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in)    :: pnts
    integer, intent(in)                    :: totpnts, t
    real(8),dimension(:), intent(in)       :: u
    character(len=50), intent(in)          :: fname
!-------------------------------------------------------------------------
    integer                                :: i
    character(len=50)                      :: file_name
    character(len=10)                      :: t_s
    character*29                           :: path='sim/out/'
!-------------------------------------------------------------------------
    write(t_s, '(i0)') t
    call system('mkdir -p ' // adjustl(trim(path)) // trim(adjustl(fname)) // "_vtk/")
    file_name = trim(adjustl(path)) // trim(adjustl(fname)) // "_vtk/" // trim(adjustl(fname))  // trim(adjustl(t_s)) // ".vtk"

    write(*,'(a)') "Printing results to: " // file_name
    open(UNIT=1, FILE=file_name)
    write(1,'(a)') "# vtk DataFile Version 1.0"
    write(1,'(a)') "Cube example"
    write(1,'(a)') "ASCII"
    write(1,*)
    write(1,'(a)') "DATASET POLYDATA"
    write(1,1) "POINTS ",totpnts," double"
    1 format(A,i0,A)

    do i = 1,totpnts
        write(1, '(F12.6, F12.6, F12.6)') pnts(i,1), pnts(i,2), u(i)
        !write(1,*)
    end do
    write(1,*) 
    write(1,3) "POINT_DATA ",totpnts
    3 format(A,i0)
    write(1,'(a)') "SCALARS myscalars double"
    write(1,'(a)') "LOOKUP_TABLE default"
    do i=1,totpnts
        write(1,'(F12.6)',advance='yes') u(i)
    end do
    close(unit=1)
end subroutine prnt_pks_vtk

subroutine read_pnts(pointset_fname, ps)
!-------------------------------------------------------------------------
! Read input points data
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    character(len=50), intent(in)          :: pointset_fname
!-------------------------------------------------------------------------
    type(pointset), intent(out)            :: ps
!-------------------------------------------------------------------------
    real(8), dimension(1,2)                :: read_txt
    real(8)                                :: maxy, maxz
    integer                                :: k
    character(len=50)                      :: file_name
!-------------------------------------------------------------------------

    file_name = trim(adjustl(pointset_fname))

    open(unit=2,file=file_name,access='sequential',status='old',action='read') !read from the existing file
    read(2, *) read_txt(1,:)
    ps%totpnts = int(read_txt(1,1))
    allocate(ps%pnts(ps%totpnts,3))
    ps%totsurfs = 0
    ps%totvols = 0
    maxy = 0.0
    maxz = 0.0
    do k = 1,ps%totpnts
        read(2,*) ps%pnts(k,:)
        if (abs(ps%pnts(k,2)) >= maxy) maxy = ps%pnts(k,2)
        if (abs(ps%pnts(k,3)) >= maxz) maxz= ps%pnts(k,3)
    end do
    write(*,'(a, i0)') "Total Points = ", ps%totpnts
    if (maxy == 0.0 .and. maxz == 0.0) then
        ps%dim = 1
    else if (maxy > 0.0 .and. maxz == 0.0) then
        ps%dim = 2
    else
        ps%dim = 3
    end if
    write(*,'(a, i0)') "Spatial dim = ", ps%dim
    close(unit=2)
end subroutine read_pnts

subroutine read_surfs(pointset_fname, ps)
!-------------------------------------------------------------------------
! Read input surface data
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    character(len=50), intent(in)       :: pointset_fname
!-------------------------------------------------------------------------
    type(pointset), intent(inout)       :: ps
!-------------------------------------------------------------------------
    real(8), dimension(1,4)             :: read_txt
    integer                             :: k
    character(len=50)                   :: file_name
    !character*29                        :: path='sim/in/1_pntst/'
!-------------------------------------------------------------------------

    !file_name = trim(adjustl(path)) // trim(adjustl(pointset_fname))
    file_name = trim(adjustl(pointset_fname))
    
    open(unit=2,file=file_name,access='sequential',status='old',action='read') !read from the existing file
    read(2, *) read_txt(1,:)
    ps%totsurfs = int(read_txt(1,1))
    ps%surfpnts = int(read_txt(1,2))
    ps%surfshape = int(read_txt(1,3))
    ps%order = int(read_txt(1,4))
    ps%totvols = 0
    allocate(ps%surfs(ps%totsurfs,ps%surfpnts))
    do k = 1,ps%totsurfs
        read(2,*) ps%surfs(k,:)
    end do
    write(*,'(a, i0)') "Total Elements = ", ps%totsurfs
    write(*,'(a, i0)') "Element pnts = ", ps%surfpnts
    write(*,'(a, i0)') "Element Shape = ", ps%surfshape
    write(*,'(a, i0)') "Element Order = ", ps%order
    close(unit=2)
end subroutine read_surfs

subroutine read_vols(pointset_fname, ps)
!-------------------------------------------------------------------------
! Read input volumes data
!-------------------------------------------------------------------------
    use pntst_struct
!-------------------------------------------------------------------------
    character(len=50), intent(in)       :: pointset_fname
!-------------------------------------------------------------------------
    type(pointset), intent(inout)       :: ps
!-------------------------------------------------------------------------
    real(8), dimension(1,4)             :: read_txt
    integer                             :: k
    character(len=50)                   :: file_name
    !character*29                        :: path='sim/in/1_pntst/'
!-------------------------------------------------------------------------

    !file_name = trim(adjustl(path)) // trim(adjustl(pointset_fname))
    file_name = trim(adjustl(pointset_fname))
    
    open(unit=2,file=file_name,access='sequential',status='old',action='read') !read from the existing file
    read(2, *) read_txt(1,:)
    ps%totvols = int(read_txt(1,1))
    ps%volpnts = int(read_txt(1,2))
    ps%volshape = int(read_txt(1,3))
    ps%order = int(read_txt(1,4))

    allocate(ps%vols(ps%totvols,ps%volpnts))
    do k = 1,ps%totvols
        read(2,*) ps%vols(k,:)
    end do
    write(*,'(a, i0)') "Total Elements = ", ps%totvols
    write(*,'(a, i0)') "Element pnts = ", ps%volpnts
    write(*,'(a, i0)') "Element Shape = ", ps%volshape
    write(*,'(a, i0)') "Element Order = ", ps%order
    close(unit=2)
end subroutine read_vols

subroutine prnt_pltctrl(pnts, nt, prnt_freq)
!-------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: nt, prnt_freq
!-------------------------------------------------------------------------
    real(8) :: xmin, xmax, ymin, ymax
    character(len=50) :: fname
    character*29 :: path='sim/out/'
!-------------------------------------------------------------------------
    !Construct the filename:
    fname = trim(adjustl(path)) // "pltctrl.txt"
    xmin = minval(pnts(:,1))
    ymin = minval(pnts(:,2))
    xmax = maxval(pnts(:,1))
    ymax = maxval(pnts(:,2))
    
    open(UNIT=1, FILE=fname)
    write(1,'(i0, a, i0, a, f10.3, a, f10.3, a, f10.3, a, f10.3)') nt, " ", prnt_freq, " ", xmin, " ", xmax, " ", ymin, " ", ymax
    close(UNIT=1)
end subroutine prnt_pltctrl

subroutine prnt_elems_txt(elems, totelems)
!-------------------------------------------------------------------------
    integer, dimension(:,:), intent(in) :: elems
    integer, intent(in):: totelems
!-------------------------------------------------------------------------
    integer :: i
    character(len=50) :: fname
    character*29 :: path='sim/out/'
!-------------------------------------------------------------------------
    !Construct the filename:
    fname = trim(adjustl(path)) // "elems.txt"
    
    open(UNIT=1, FILE=fname)
    !nd number is minus 1 to be able to plot it in matplotlib
    do i = 1, totelems
        write(1, '(i0, a, i0, a, i0)') elems(i,1)-1, " ", elems(i,2)-1, " ", elems(i,3)-1
    end do
    close(UNIT=1)
end subroutine prnt_elems_txt
end module pntst


