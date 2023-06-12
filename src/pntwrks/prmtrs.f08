module prmtrs
implicit none
contains

subroutine read_array(fname, data, file_exists)
    use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    character*50, intent(in) :: fname
    integer :: error
    real(8) :: x
    real(8), dimension(:), allocatable, intent(out)      :: data
    logical, intent(out)                                :: file_exists

    allocate(data(0))
    inquire(file = fname, exist = file_exists)
    if (file_exists .eqv. .true.) then
        open (1, file = fname, status = 'old')
        do
            read(1, *, iostat = error) x
            select case(error)
            case (0)
                data = [data , x]
            case(iostat_end)
                exit
            case Default
                write(*,*) 'Error in reading file'
                stop
            end select
        end Do
        close(1)
    end if
end subroutine

subroutine set_rb(sp)
    use slvr_prmtrs_struct
    implicit none

    type(slvr_prmtrs), intent(inout)        :: sp   !solver parameters data structure
    character*50 :: fname = 'sim/in/5_ctrl/rb.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists
    call read_array(fname, data, file_exists)

    if (file_exists .eqv. .true.) then
        sp%t = data(1)      
        sp%thot = data(2)   
        sp%tcold = data(3)  
        sp%pr = int(data(4))
        sp%ra = data(5)     
   else
        write(*,'(a)') "ERROR: sim/in/5_ctrl/krnl.txt does not exist. Can not proceed without this file. Exiting..."
        call exit()
    end if
end subroutine set_rb

subroutine set_krnl(dx, sp)
    use slvr_prmtrs_struct
    implicit none

    real(8), intent(in) :: dx
    type(slvr_prmtrs), intent(inout)        :: sp   !solver parameters data structure
    character*50 :: fname = 'sim/in/5_ctrl/krnl.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists
    call read_array(fname, data, file_exists)

    if (file_exists .eqv. .true.) then
        sp%order = int(data(1))            !intrpolation order 1 or 2
        sp%krnl = int(data(2))             !kernel type: 1 = rbf, 2 = mls, 3 = wls, 4 = sph, 5 = gfd
        sp%h = data(3) * dx        !kernel radius
        sp%rbf = int(data(4))              !rbf weight funct: 1 = mq, 2 = imq, 3 = ge
        sp%rbf_alpha = data(5)     !rbf parameter
        sp%rbf_polyex = .false.
        if (int(data(6)) == 1) sp%rbf_polyex = .true.
        sp%mls = int(data(7))              !mls weight funct: 1 = s3, 2 = s4, 3 = s5, 4 = reg_spln
        sp%wls = int(data(8))              !wls weight funct: 0 = ls, 1 = ge, 2 = s3, 3 = s4, 4 = s5
        sp%sph = int(data(9))              !sph weight funct: 1 = ge, 2 = s3, 3 = w4, 4 = w5, 5 = inv_dis
        sp%gfd = int(data(10))              !gfd weight funct: 1 = s3, 2 = s4, 3 = ge
    else
        write(*,'(a)') "ERROR: sim/in/5_ctrl/krnl.txt does not exist. Can not proceed without this file. Exiting..."
        call exit()
    end if
end subroutine set_krnl

subroutine set_t(sp)

    use slvr_prmtrs_struct
    implicit none

    type(slvr_prmtrs), intent(inout)        :: sp   !solver parameters data structure
    real(8), dimension(:), allocatable      :: data
    character*50 :: fname = 'sim/in/5_ctrl/t.txt'
    logical :: file_exists
    call read_array(fname, data, file_exists)

    if (file_exists .eqv. .true.) then
        sp%dt = data(1)          !time step size
        sp%nt = int(data(2))           !total time steps
        sp%itrs = int(data(3))             !iterations
        sp%tol = data(4)                    !convergence tolerance
        sp%prnt_frq = int(data(5))       !print frequency
        sp%vtk = .false.         !output data format vtk or txt
        if (int(data(6)) == 1) sp%vtk = .true.
    else
        write(*,'(a)') "ERROR: sim/in/5_ctrl/t.txt does not exist. Can not proceed without this file. Exiting..."
        call exit()
    end if
end subroutine set_t

subroutine set_phi(dx, sp)

    use slvr_prmtrs_struct
    implicit none
    real(8), intent(in) :: dx
    type(slvr_prmtrs), intent(inout)        :: sp   !solver parameters data structure
    character*50 :: fname = 'sim/in/5_ctrl/phi.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists
    call read_array(fname, data, file_exists)

    if (file_exists .eqv. .true.) then
        sp%intrf_slvr = int(data(1))          !1 = vof1, 2 = vof2, 3 = ac, 4 = ac2, 5 = ch, 6 = ch2
        sp%w = data(2) * dx        !related to intrf thickness
        sp%m = data(3)             !mobility parameter
        sp%segma = data(4)            !surface tension parameter
        sp%sharp_intrf = .false.   !sharpen the diffuse intrf flag
        if (int(data(5)) == 1) sp%sharp_intrf = .true.
        sp%si_beta = data(6) * dx   !intrf sharpening diffusive parameter
        sp%si_eps = data(7) * dx    !intrf sharpening prmtr
        sp%si_dt = data(8)          !virtual time step size for intrf sharpening
        sp%si_nt = int(data(9))               !virtual time steps for intrf sharpening
        if (size(data) > 9) then
            sp%invert_phi  = int(data(10))           !calculate the inerse of phi
        end if
    else
        write(*,'(a)') "sim/in/5_ctrl/phi.txt does not exist"
    end if
end subroutine set_phi

subroutine set_sld(sp)

    use slvr_prmtrs_struct
    implicit none
    type(slvr_prmtrs), intent(inout)        :: sp   !solver parameters data structure
    character*50 :: fname = 'sim/in/5_ctrl/sld.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists
    call read_array(fname, data, file_exists)

    if (file_exists .eqv. .true.) then
        sp%a_2 = data(1)          !controls the number of primary branches
        sp%delta = data(2)       !the smaller it is the more dendrite splitting
        sp%eps_4 = data(3)        !directly related to the intrf thickness
        sp%theta0 = data(4)        !preferential growth angle
        sp%tau0 = data(5)       !0.0003;
        sp%alpha = data(6)        !smaller it is the smaller the solidification rate
        sp%gamma = data(7)         !the higher it is the more non-constant the intrf thickness. smaller it is the more the surface tension
        sp%Teq = data(8)            !equilibrium temperature
        sp%lhf = data(9)          !latent heat of fusion
    else
        write(*,'(a)') "sim/in/5_ctrl/phi.txt does not exist"
    end if
end subroutine set_sld

subroutine set_v(sp)

    use slvr_prmtrs_struct
    use intrf
    implicit none

    type(slvr_prmtrs), intent(inout)        :: sp   !solver parameters data structure
    character*50 :: fname = 'sim/in/5_ctrl/v.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists
    call read_array(fname, data, file_exists)
   
    if (file_exists .eqv. .true.) then
        sp%lgr = .false.
        if (int(data(1)) == 1) sp%lgr = .true.
        sp%p = data(2)
        sp%av = data(3)
        sp%upwind = .false.
        if (int(data(4)) == 1) sp%upwind = .true.
        sp%shift = .false.       !particle shifting
        if (int(data(5)) == 1) sp%shift = .true.
        sp%shift_prmtr = data(6)
        sp%shift_surf = data(7)
        sp%periodic = .false.
    else
        write(*,'(a)') "sim/in/5_ctrl/v.txt does not exist"
    end if
end subroutine set_v

subroutine set_v0(pnts, phi, totpnts, v)

    use slvr_prmtrs_struct
    use intrf
    implicit none

    real(8), dimension(:,:), intent(in) :: pnts
    real(8), dimension(:), intent(in)   :: phi
    integer, intent(in) :: totpnts
    integer :: i, adv_v
    real(8) :: vx1, vy1, vz1, vx2, vy2, vz2, damp
    real(8), dimension(:,:), allocatable, intent(out) :: v

    character*50 :: fname = 'sim/in/3_init/v0.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists

    call read_array(fname, data, file_exists)

    allocate(v(totpnts, 3))
    v(:,:) = 0.0   

    if (file_exists .eqv. .true.) then
        adv_v = int(data(1))
        damp = data(2)
        vx1 = data(3)
        vy1 = data(4)
        vz1 = data(5)
        !if (size(data) > 5) then
            vx2 = data(6)
            vy2 = data(7)
            vz2 = data(8)
        !end if

        select case (adv_v)
        case (1)
            !if (size(data) > 5) then
                do i = 1, totpnts
                    if (phi(i) <= 0.5) then
                        v(i,1) = vx1
                        v(i,2) = vy1
                        v(i,3) = vz1
                    else
                        v(i,1) = vx2
                        v(i,2) = vy2
                        v(i,3) = vz2
                    end if
                end do
           ! else
           !     v(:,1) = vx1
           !     v(:,2) = vy1
           !     v(:,3) = vz1
           ! end if
        case (2)
            call set_vrtx_vel(v, pnts, totpnts)
        case (3)
            call set_shearing_vel(v, pnts, totpnts)
            v = damp * v
        case (4)
            call set_rot_vel(v, pnts, totpnts, vx1, vy1)
        end select
    else
        write(*,'(a)') "sim/in/3_init/v0.txt does not exist. Setting v(:,:) = 0.0"
    end if
end subroutine set_v0

subroutine set_vn0(pnts, totpnts, vn)

    use slvr_prmtrs_struct
    use intrf
    implicit none

    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: totpnts
    integer :: nrm_v
    real(8) :: damp
    real(8), dimension(:), allocatable, intent(out) :: vn

    character*50 :: fname = 'sim/in/3_init/vn0.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists

    call read_array(fname, data, file_exists)

    allocate(vn(totpnts))
    vn(:) = 0.0
    
    if (file_exists .eqv. .true.) then
        nrm_v = int(data(1))
        damp = data(2)
        
        select case (nrm_v)
        case (1)
            vn(:) = damp
        case (2)
            call set_branch_vn(vn, pnts, totpnts, 4.0d0, damp)
        case (3)
            call set_branch_vn(vn, pnts, totpnts, 6.0d0, damp)
        end select
    else
        write(*,'(a)') "sim/in/3_init/vn0.txt does not exist. Setting vn(:) = 0.0"
    end if

end subroutine set_vn0

subroutine set_u0(phi, totpnts, u)

    use slvr_prmtrs_struct
    use intrf
    use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    real(8), dimension(:), intent(in) :: phi
    integer, intent(in) :: totpnts
    integer                                 :: i
    real(8)                                 :: u1, u2
    real(8), dimension(:), allocatable, intent(out) :: u

    character*50 :: fname = 'sim/in/3_init/u0.txt'
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists

    call read_array(fname, data, file_exists)

    allocate(u(totpnts))
    if (file_exists .eqv. .true.) then
        if (size(data) > 1)  then
            u1 = data(1)
            u2 = data(2)
            u(:) = 0.0
            do i = 1, totpnts
                if (phi(i) <= 0.5) then
                    u(i) = u1
                else 
                    u(i) = u2
                end if
            end do
        else
            u(:) = data(1)
        end if
    else
        write(*,'(a)') "sim/in/3_init/u0.txt does not exist. Setting u(:) = 0.0"
        u(:) = 0.0
    end if
end subroutine set_u0

subroutine set_mat(matname, k)

    use slvr_prmtrs_struct
    use intrf
    implicit none

    character(len=*),  intent(in) :: matname
    character(len=50) :: fname
    real(8), dimension(:), allocatable, intent(out) :: k
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists

    fname = 'sim/in/2_mat/' // trim(adjustl(matname)) // '.txt'
    
    call read_array(fname, data, file_exists)

    if (file_exists .eqv. .true.) then
        if (size(data) > 1) then
            k = data
        else
            allocate(k(2))
            k(:) = data(1)
        end if
    else
        write(*,'(a,a)') trim(adjustl(fname)), " does not exist. Setting it to 0.0"
        allocate(k(2))
        k(:) = 0.0
    end if
end subroutine set_mat

subroutine set_geom(ps)
    use pntst_struct
    use pntst
    implicit none

    integer                    :: i
    real(8), dimension(3)      :: trans, b, l
    integer, dimension(3)      :: n
    real(8)                    :: scale
    type(pointset), intent(out)           :: ps
    character(len=50)                     :: pnts_fname, surf_fname, vol_fname
    character(len=50) :: fname1, fname2, fname3
    real(8), dimension(:), allocatable      :: data
    logical :: pntst_exists, grd_exists, trns_exists
    fname1 = 'sim/in/1_pntst/pnts.txt'
    fname2 = 'sim/in/1_pntst/grd.txt'
    fname3 = 'sim/in/1_pntst/trns.txt'

    inquire(file = fname1, exist = pntst_exists)
    call read_array(fname2, data, grd_exists)
    if (pntst_exists .eqv. .true.) then
        pnts_fname = 'sim/in/1_pntst/pnts.txt'
        surf_fname = 'sim/in/1_pntst/surfs.txt'   
        vol_fname = 'sim/in/1_pntst/vols.txt'           
        call read_pnts(pnts_fname, ps)  
        call read_surfs(surf_fname, ps)
        if (ps%dim == 3) call read_vols(vol_fname, ps)
    else if (grd_exists .eqv. .true.) then
        b = data(1:3)          !base point
        l = data(4:6)          !domain length
        n = int(data(7:9))    !domain discritization
        if (n(2) <= 1 .and. n(3) <= 1) then
            call gen_grdpnts1d(b, l, n, ps)
            write(*,*) ps%totpnts
        else if (n(3) <= 1) then
            call gen_quadgrd(b, l, n, ps)
        else
            call gen_hexgrd(b, l, n, ps)
        end if
    else
        write(*,*) "ERROR: no pntst.txt or grd.txt exist in sim/in/1_pntst/"
        write(*,*) "Solver can not proceed without either file. Exiting..."
        call exit()
    end if

    call read_array(fname3, data, trns_exists)
    if (trns_exists .eqv. .true.) then
        scale = data(1)       !geometry scale factor
        trans = data(2:4)    !geometry translation
        ps%pnts = ps%pnts * scale
        do i = 1, 3
            ps%pnts(:, i) = ps%pnts(:, i) + trans(i)
        end do
        if (size(data) > 4) then
            ps%dx = data(5)    !override average point spacing
        else
            call get_avg_dx(ps%pnts, ps%dx)
         end if
    else
        call get_avg_dx(ps%pnts, ps%dx)
    end if

end subroutine set_geom

subroutine set_geom_arr(ps)
    use pntst_struct
    use pntst
    implicit none

    integer                    :: i
    real(8), dimension(3)      :: trans, b, l
    integer, dimension(3)      :: n
    real(8)                    :: scale
    type(pointset), dimension(2), intent(out)           :: ps
    character(len=50)                     :: pnts_fname, surf_fname, vol_fname
    character(len=50) :: fname1, fname2, fname3
    real(8), dimension(:), allocatable      :: data
    logical :: pntst_exists, grd_exists, trns_exists
    integer :: zz

    fname1 = 'sim/in/1_pntst/pnts_1.txt'
    fname2 = 'sim/in/1_pntst/grd_1.txt'
    fname3 = 'sim/in/1_pntst/trns_1.txt'

    pnts_fname = 'sim/in/1_pntst/pnts_1.txt'
    surf_fname = 'sim/in/1_pntst/surfs_1.txt'   
    vol_fname = 'sim/in/1_pntst/vols_1.txt'    

    do zz = 1, 2
        inquire(file = fname1, exist = pntst_exists)
        call read_array(fname2, data, grd_exists)
        if (pntst_exists .eqv. .true.) then
        
            call read_pnts(pnts_fname, ps(zz))  
            call read_surfs(surf_fname, ps(zz))
            if (ps(zz)%dim == 3) call read_vols(vol_fname, ps(zz))
        else if (grd_exists .eqv. .true.) then
            b = data(1:3)          !base point
            l = data(4:6)          !domain length
            n = int(data(7:9))    !domain discritization
            if (n(2) <= 1 .and. n(3) <= 1) then
                call gen_grdpnts1d(b, l, n, ps(zz))
                write(*,*) ps%totpnts
            else if (n(3) <= 1) then
                call gen_quadgrd(b, l, n, ps(zz))
            else
                call gen_hexgrd(b, l, n, ps(zz))
            end if
        else
            write(*,*) "ERROR: no pntst.txt or grd.txt exist in sim/in/1_pntst/"
            write(*,*) "Solver can not proceed without either file. Exiting..."
            call exit()
        end if

        call read_array(fname3, data, trns_exists)
        if (trns_exists .eqv. .true.) then
            scale = data(1)       !geometry scale factor
            trans = data(2:4)    !geometry translation
            ps(zz)%pnts = ps(zz)%pnts * scale
            do i = 1, 3
                ps(zz)%pnts(:, i) = ps(zz)%pnts(:, i) + trans(i)
            end do
            if (size(data) > 4) then
                ps(zz)%dx = data(5)    !override average point spacing
            else
                call get_avg_dx(ps(zz)%pnts, ps(zz)%dx)
            end if
        else
            call get_avg_dx(ps(zz)%pnts, ps(zz)%dx)
        end if

        fname1 = 'sim/in/1_pntst/pnts_2.txt'
        fname2 = 'sim/in/1_pntst/grd_2.txt'
        fname3 = 'sim/in/1_pntst/trns_2.txt'

        pnts_fname = 'sim/in/1_pntst/pnts_2.txt'
        surf_fname = 'sim/in/1_pntst/surfs_2.txt'   
        vol_fname = 'sim/in/1_pntst/vols_2.txt'
    end do  
end subroutine set_geom_arr

subroutine set_bins(sp)
    use slvr_prmtrs_struct

    use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    type(slvr_prmtrs), intent(inout)        :: sp   !solver parameters data structure
    character(len=50) :: fname
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists
    fname = 'sim/in/4_bc/bins.txt'

    call read_array(fname, data, file_exists)
    if (file_exists .eqv. .true.) then
        sp%bg_nx = int(data(1))    
        sp%bg_b = data(2:4)          !base point
        sp%bg_l = data(5:7)          !domain length
    else
        write(*,*) "sim/in/4_bc/bins.txt does not exist."
    end if
end subroutine set_bins

subroutine set_phi0(pnts, totpnts, dx, phi)
    use intrf

    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: totpnts
    real(8), intent(in) :: dx
    real(8), dimension(:), allocatable, intent(out) :: phi
    real(8), dimension(:,:), allocatable :: qb, sphr, cyl, cosf
    real(8), dimension(:), allocatable :: rand
    real(8) :: r, a0, b0, c0, d0, e0, noise, num !random noise parameters
    real(8), dimension(3) :: c, sz
    real(8) :: sgn, axis, val
    integer :: k, totseeds
    logical :: c1, c2, c3, c4, c5
    character(len=50) :: fname

    allocate(phi(totpnts))
    phi(:) = 1.0

    inquire( file = 'sim/in/3_init/phi0_pln.txt', exist = c1)
    if (c1 .eqv. .true.) then
        open(unit = 1, file = 'sim/in/3_init/phi0_pln.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(1, *) sgn
        read(1, *) axis
        read(1, *) val
        close (1)
        phi = sgn*(pnts(:,int(axis)) + val)
    end if

    inquire( file = 'sim/in/3_init/phi0_sphr.txt', exist = c1)
    if (c1 .eqv. .true.) then
        open(unit = 1, file = 'sim/in/3_init/phi0_sphr.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(1, *) totseeds
        allocate(sphr(totseeds, 5))
        do k = 1, totseeds
            read(1,*) sphr(k,:)
            c = sphr(k, 2:4)
            r = sphr(k, 5)
            if (int(sphr(k,1)) == 1) then
                call add_sphr(c, r, pnts, totpnts, phi)
            else
                call sub_sphr(c, r, pnts, totpnts, phi)
            end if
        end do
        close (1)
    end if

    inquire( file = 'sim/in/3_init/phi0_cyl.txt', exist = c2)
    if (c2 .eqv. .true.) then
        open(unit = 1, file = 'sim/in/3_init/phi0_cyl.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(1, *) totseeds
        allocate(cyl(totseeds, 5))
        do k = 1, totseeds
            read(1,*) cyl(k,:)
            c = cyl(k, 2:4)
            r = cyl(k, 5)
            if (int(cyl(k,1)) == 1) then
                call add_cyl_z(c, r, pnts, totpnts, phi)
            else
                call sub_cyl_z(c, r, pnts, totpnts, phi)
            end if
        end do
        close (1)
    end if

    inquire( file = 'sim/in/3_init/phi0_qb.txt', exist = c3)
    if (c3 .eqv. .true.) then
        open(unit = 2, file = 'sim/in/3_init/phi0_qb.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(2, *) totseeds
        allocate(qb(totseeds, 7))
        do k = 1, totseeds
            read(2,*) qb(k,:)
            c = qb(k, 2:4)
            sz = qb(k, 5:7)
            if (int(qb(k,1)) == 1) then
                call add_cube(c, sz, pnts, totpnts, phi)
            else
                call sub_cube(c, sz, pnts, totpnts, phi)
            end if
        end do
        close (2)
    end if

    inquire( file = 'sim/in/3_init/phi0_cos.txt', exist = c4)
    if (c4 .eqv. .true.) then
        open(unit = 2, file = 'sim/in/3_init/phi0_cos.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(2, *) totseeds
        allocate(cosf(totseeds, 5))
        do k = 1, totseeds
            read(2,*) cosf(k,:)
            a0 = cosf(k, 1)
            b0 = cosf(k, 2)
            c0 = cosf(k, 3)
            d0 = cosf(k, 4)
            e0 = cosf(k, 5)
            if (k == 1) then
                !phi = ((ps%pnts(:,2) - (0.006 * cos(25.0 * ps%pnts(:,1) + 0.1) + 0.0)) - 0.5)
                phi = ((pnts(:,2) - (a0 * cos(b0 * pnts(:,1) + c0) + d0)) - e0)
            else
                !phi = phi - (0.006 * cos(25.0 * ps%pnts(:,3) + 0.1) + 0.0)
                phi = phi - (a0 * cos(b0 * pnts(:,3) + c0) + d0) - e0
            end if
        end do
        close (2)
    end if
    !if (c1 .eqv. .true. .or. c2 .eqv. .true. .or. c3 .eqv. .true. .or. c4 .eqv. .true.) then
        call get_vof(phi, totpnts, dx)
    !end if

    fname = 'sim/in/3_init/phi0_rand.txt'
    call read_array(fname, rand, c5)
    if (c5 .eqv. .true.) then
        c0 = rand(1)   
        noise = rand(2)      
        do k = 1, totpnts
            call random_number(num)
            phi(k) = c0 + noise * (0.5 - num)
        end do
    end if

end subroutine set_phi0

subroutine set_bc(pnts, totpnts, bcs)
    use bc_struct
    use bndry

    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: totpnts
    real(8), dimension(:,:), allocatable :: data
    real(8) :: bcval, bcloc
    integer :: i, totbcs, bcaxis, bctype
    logical :: check
    type(bc), dimension(6), intent(out) :: bcs
    character*100, dimension(:,:), allocatable :: fname
    character*100 :: file_name
    character*50                        :: path='sim/in/4_bc/'

    !bc types: 1 = u, 2 = phi, 3 = vx, 4 = vy, 5 = vz, 6 = p
    do i = 1, 6
        bcs(i)%tot = 0
    end do
    inquire( file = 'sim/in/4_bc/bc.txt', exist = check)
    if (check .eqv. .true.) then
        open(unit = 1, file = 'sim/in/4_bc/bc.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(1, *) totbcs
        allocate(data(totbcs, 4))
        do i = 1, totbcs
            read(1,*) data(i,:)
        end do
        close (1)
        
        !maximum number of bc types possible is 6 (vx, vy, vz, p, u, phi) or (ux, uy, uz, fx, fy, fz)
        do i = 1, 6
            bcs(i)%tot = 0
            allocate(bcs(i)%pnts(0))
            allocate(bcs(i)%vals(0))
        end do

        do i = 1, totbcs
            bctype = int(data(i, 1))
            bcaxis = int(data(i, 2))
            bcloc = data(i, 3)
            bcval = data(i, 4)
            call set_dbc(pnts, totpnts, bcaxis, bcloc, bcval, bcs(bctype))
        end do
        close (1)
    end if
    inquire(file = 'sim/in/4_bc/bc_f.txt', exist = check)
    if (check .eqv. .true.) then
        open(unit = 2, file = 'sim/in/4_bc/bc_f.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(2, *) totbcs                   
        allocate(fname(totbcs, 3))
        do i = 1, totbcs
            read(2,*) fname(i,:)
        end do

        do i = 1, 6 !total possible bc types is 6 (vx, vy, vz, p, u, phi) or (ux, uy, uz, fx, fy, fz)
            bcs(i)%tot = 0
            allocate(bcs(i)%pnts(0))
            allocate(bcs(i)%vals(0))
        end do

        do i = 1, totbcs
            read(fname(i,1),*)  bctype
            read(fname(i,3), '(f10.0)' )  bcval 
            file_name = trim(adjustl(path)) // trim(adjustl(fname(i,2)))
            call read_dbc(file_name, bcval, pnts, totpnts, bcs(bctype))
        end do
    end if

end subroutine set_bc

subroutine set_bc_phi(phi, pnts, totpnts, bcs)
    use bc_struct
    use bndry
    real(8), dimension(:), intent(in) :: phi
    real(8), dimension(:,:), intent(in) :: pnts
    integer, intent(in) :: totpnts
    real(8), dimension(:,:), allocatable :: data
    real(8) :: bcloc
    integer :: i, totbcs, bcaxis, bctype
    logical :: check
    type(bc), dimension(6), intent(inout) :: bcs

    inquire( file = 'sim/in/4_bc/bc_phi.txt', exist = check)
    if (check .eqv. .true.) then
        open(unit = 1, file = 'sim/in/4_bc/bc_phi.txt', access = 'sequential', status='old', action = 'read') !read from the existing file
        read(1, *) totbcs
        allocate(data(totbcs, 3))
        do i = 1, totbcs
            read(1,*) data(i,:)
        end do
        close (1)
        
        !maximum number of bc types possible is 6 (vx, vy, vz, p, u, phi) or (ux, uy, uz, fx, fy, fz)
        do i = 1, totbcs
            bctype = int(data(i, 1))
            bcaxis = int(data(i, 2))
            bcloc = data(i, 3)
            call set_dbc2(pnts, totpnts, bcaxis, bcloc, phi, bcs(bctype))
        end do
        close (1)
    end if
end subroutine set_bc_phi

subroutine set_slvr(slvr)

    implicit none

    integer, intent(out)        :: slvr   !solver parameters data structure
    character(len=50) :: fname
    real(8), dimension(:), allocatable      :: data
    logical :: file_exists
    fname = 'sim/in/5_ctrl/slvr.txt'

    call read_array(fname, data, file_exists)
    if (file_exists .eqv. .true.) then
        slvr = int(data(1))          !time step size
    else
        write(*,'(a)') "ERROR: file sim/in/5_ctrl/slvr.txt does not exist. Cannot proceed. Exiting..."
        call exit()
    end if
  
end subroutine set_slvr

subroutine set_sctr(sctrpnts, totsctr)
    use pntst
    use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    real(8), dimension(3) :: b, l
    real(8), dimension(:,:), allocatable, intent(out) :: sctrpnts
    integer, intent(out) :: totsctr
    character(len=50) :: fname
    real(8), dimension(:), allocatable :: data
    logical :: file_exists
    fname = 'sim/in/1_pntst/sctr.txt'

    call read_array(fname, data, file_exists)
    
    if (file_exists .eqv. .true.) then
        totsctr = int(data(1))    !grid generation flag
        b = data(2:4)          !base point
        l = data(5:7)          !domain length
        call gen_sctr_pnts(b, l, totsctr, sctrpnts)
    else
        write(*,'(a)') "sim/in/1_pntst/sctr.txt does not exist..."
    end if

end subroutine set_sctr

end module prmtrs



