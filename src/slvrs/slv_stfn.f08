module stfn_slvrs
implicit none
contains

subroutine run_stfn ()
!------------------------------------------------------------------------------------------
!   This subroutine solves the sharp-interface stefan equations for dendritic
!   solidification in a pure liquid. Solution is obtained using explicit time 
!   stepping schemes and local strong-form meshfree methods.     
!------------------------------------------------------------------------------------------
    use pntst_struct
    use mat_struct
    use slvr_prmtrs_struct
    use bc_struct
    use pntst
    use intrf
    use prmtrs
    use omp_lib
!------------------------------------------------------------------------------------------
    real(8)                                 :: start, finish !timer parameters
    real(8), dimension(:), allocatable      :: u        !thermal field distribution
    real(8), dimension(:), allocatable      :: phi      !phase-field parameter distribution
    real(8), dimension(:,:), allocatable    :: v        !advective velocity in xyz
    type(pointset)                          :: ps       !pointset data structure
    type(slvr_prmtrs)                       :: sp       !solver parameters data structure
    type(materials)                         :: mat      !phase properties data structures
    type(bc), dimension(6)                  :: bcs      !boundary condition data structure for vx, vy, vz, p
!------------------------------------------------------------------------------------------
    write(*,'(a)') "--------------------------------------------------------"
    write(*,'(a)') "Running Stefan solver..."
    write(*,'(a)') "--------------------------------------------------------"
    !constructing the pointset
    call set_geom(ps)
    !constructing the implicit geometry
    call set_phi(ps%dx, sp)
    !call set_phi0(ps%pnts, ps%totpnts, ps%dx, phi)
    call set_phi0(ps%pnts, ps%totpnts, 0.25*ps%dx, phi)
    !initializing the material properties
    call set_mat('k', mat%d)
    call set_mat('nu', mat%nu)
    call set_mat('ro', mat%ro)
    !initializing the field variables
    call set_v(sp)
    call set_v0(ps%pnts, phi, ps%totpnts, v)
    call set_u0(phi, ps%totpnts, u)
    !assigning the boundary conditions
    call set_bc(ps%pnts, ps%totpnts, bcs)
    !assigning the solver parameters
    call set_sld(sp)
    call set_krnl(ps%dx, sp)
    call set_t(sp)
    !running the solver
    start = omp_get_wtime()
    !call slv_ac_singleseed_2d(ps, sp, mat, bcs, u, phi, v)
    call slv_stfn(ps, sp, mat, bcs, u, phi, v)
    finish = omp_get_wtime()
    print '("Time = ",f8.1," seconds.")',finish-start
end subroutine run_stfn

subroutine slv_stfn(ps, sp, mat, bcs, u, phi, vel0)
!------------------------------------------------------------------------------------------
    use krnl_struct
    use pntst_struct
    use bc_struct
    use slvr_prmtrs_struct
    use mat_struct
    use pntst
    use intrf
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_mls
    use krnl_wls
    use krnl_krg
    use bndry
    use trnsprt
    use ns
    use ns_cbs
!------------------------------------------------------------------------------------------
    type(pointset), intent(in)              :: ps
    type(bc), dimension(5), intent(in)      :: bcs
    type(slvr_prmtrs), intent(in)           :: sp
!------------------------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)    :: phi
    real(8), dimension(:,:), intent(inout)  :: vel0
    type(materials), intent(inout)          :: mat
!------------------------------------------------------------------------------------------
    type (kernel), dimension(ps%totpnts)    :: krnls
    type (kernel), dimension(:), allocatable :: ikrnls
    integer                                 :: c = 1, t
    real(8), dimension(ps%totpnts)          :: vf, u, vf_old, q
    real(8)                                 :: start, finish, delta_ratio, d_nrm, eps
    integer, dimension(:), allocatable      :: pntbins
    integer                                 :: nx, nxy, totbins
    type(kernel), dimension(:), allocatable :: bins
    character(len=50)                       :: phi_fname = 'phi', u_fname = 'u'
    character(len=50)                        :: v_fname = 'v', p_fname = 'p', vn_fname = 'vn'
    !1character(len=50)                        :: ip_fname = 'ip', sol_fname = 'sol', liq_fname = 'liq'
    logical                                 :: fluid_slvr = .false.
    real(8), dimension(:), allocatable      :: v, p, p0, ro, vn, vn_ext, crv
    real(8), dimension(:,:), allocatable    :: vel, f, ip, np_liq, np_sol
    real(8), dimension(:), allocatable :: T_intrf, icrv, ivn
    integer :: totip, i, j, nbr
    logical :: crv_chk = .true., mob_chk = .true.
    real(8) :: d0, a, b, theta0, K_liq, K_sol
    integer :: st_func
    integer, dimension(:), allocatable :: prj_pntnums
    real(8), dimension(:), allocatable :: prj_pntvals
    integer :: totprjpnts
!------------------------------------------------------------------------------------------
    !Initialization of variables
    c = 1
    eps = 1e-6
    !delta_ratio = sp%delta
    !d_nrm = delta_ratio * ps%dx
    !K_sol = mat%d(1)
    !K_liq = mat%d(2)
    !st_func = int(sp%gamma)
    !b = sp%a_2
    !a = sp%alpha
    !d0 = sp%eps_4
    !theta0 = sp%theta0

    delta_ratio = 1.5
    d_nrm = delta_ratio * ps%dx
    K_sol = 0.5
    K_liq = 2.0
    st_func = 1
    b = 4.0
    a = 1
    d0 = 0.001
    theta0 = 0.0

    vf = phi
    vf_old = vf
    phi = -(2.0 * vf - 1.0)
    
    allocate(crv(ps%totpnts))
    allocate(vn_ext(ps%totpnts))
    
    if (mat%nu(1) > 0) then
        fluid_slvr = .true.
        vel0(bcs(1)%pnts, 1) = bcs(1)%vals
        vel0(bcs(2)%pnts, 2) = bcs(2)%vals
        vel0(bcs(3)%pnts ,3) = bcs(3)%vals
        allocate(p(ps%totpnts))
        allocate(ro(ps%totpnts))
        allocate(v(ps%totpnts))
        allocate(f(ps%totpnts, 3))
        p(:) = 0.0
        ro(:) = mat%ro(1)
        f(:,:) = 0.0
        p0 = p
    end if
    vel = vel0
!------------------------------------------------------------------------------------------
    ! Calculate bounds of neighbour search bins
    nx = 30
    call gen_bins(ps%pnts, ps%totpnts, ps%dim, ps%dx, nx, nxy, totbins, bins, pntbins)
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Constructing krnls at fixed points..."
    call cpu_time(start)
    !call get_nbrs_bf(sp%h, ps%pnts, ps%totpnts, ps%pnts, ps%totpnts, krnls)
    call get_nbrs_bg(ps%pnts, pntbins, ps%totpnts, bins, nx, nxy, totbins, ps%dim, sp%h, krnls)
    select case(sp%krnl)
        case (1)
            call get_rbf_krnls(ps%pnts, ps%totpnts, ps%pnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, sp%rbf_polyex, krnls)
        case (2)
            call get_mls_krnls(ps%pnts, ps%totpnts, ps%pnts, sp%mls, sp%order, sp%h, krnls)
        case (3)
            call get_wls_krnls(ps%pnts, ps%totpnts, ps%pnts, sp%wls, sp%order, sp%h, krnls)
        case (4)
            call get_sph_krnls(ps%pnts, ps%totpnts, ps%pnts, sp%sph, ps%dim, sp%h,  krnls)
        case (5)
            call get_krg_krnls(ps%pnts, ps%totpnts, ps%pnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, krnls)
    end select
    call get_intrps_o2(ps%pnts, ps%totpnts, krnls)
    call cpu_time(finish)
    print '("intrpolants construction time = ",f10.1," seconds")',finish-start
!------------------------------------------------------------------------------------------
    if (sp%vtk .eqv. .true.) then
        call prnt_vtk(u, ps, u_fname, 0)
        select case (ps%dim)
            case(3)
                call prnt_unstruct_vtk(phi, ps, phi_fname, 0)
            case default
                call prnt_vtk(vf, ps, phi_fname, 0)
        end select
        if (fluid_slvr .eqv. .true.) then
            v = sqrt(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2))
            call prnt_vtk(v, ps, v_fname, 0)
            call prnt_vtk(p, ps, p_fname, 0)
        end if
    else
        call prnt_pltctrl(ps%pnts, sp%nt, sp%prnt_frq)
        call prnt_elems_txt(ps%surfs, ps%totsurfs)
        call prnt_txt(u, ps%pnts, ps%totpnts, u_fname, 0)
        call prnt_txt(vf, ps%pnts, ps%totpnts, phi_fname, 0)
    end if

    write(*,'(a)') "Computing transient solution..."
    do t = 1, sp%nt

       !generate interface points at the zero-isocontour of phi
       call get_intrf_pnts(ps%pnts, phi, krnls,  ps%totpnts, ps%dx, ip, totip)
       write(*,*) "total interface points generated = ", totip

       allocate(vn(totip))
       allocate(ikrnls(totip))
       allocate(icrv(totip))
       allocate(ivn(totip))

       !Constructing krnls at interface points...
       call get_nbrs_bf(sp%h, ip, totip, ps%pnts, ps%totpnts, ikrnls)
       select case(sp%krnl)
            case (1)
                call get_rbf_krnls(ip, totip, ps%pnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, sp%rbf_polyex, ikrnls)
            case (2)
                call get_mls_krnls(ip, totip, ps%pnts, sp%mls, sp%order, sp%h, ikrnls)
            case (3)
                call get_wls_krnls(ip, totip, ps%pnts, sp%wls, sp%order, sp%h, ikrnls)
            case (4)
                call get_sph_krnls(ip, totip, ps%pnts, sp%sph, ps%dim, sp%h,  ikrnls)
            case (5)
                call get_krg_krnls(ip, totip, ps%pnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, ikrnls)    
        end select
        call get_intrps_o22(ip, totip, ps%pnts, krnls)

	   !Computing interfacial normal points...
       call get_nrm_pnts(ip, totip, d_nrm, phi, ikrnls, np_sol, np_liq)
       !call prnt_txt(uip, np_sol, totip, sol_fname, 0)
       !call prnt_txt(uip, np_liq, totip, liq_fname, 0)

        if (crv_chk .eqv. .true.) then
            !Computing curvature at interface points...
            call get_curv(phi, ps%pnts, krnls, ps%totpnts, crv)
            do i = 1, totip
                icrv(i)  = 0.0
                ivn(i) = 0.0
                do j = 1, ikrnls(i)%totnbrs
                    nbr =  ikrnls(i)%nbrs(j)
                    icrv(i) = icrv(i) + ikrnls(i)%n(j) * crv(nbr)
                    ivn(i) = ivn(i) + ikrnls(i)%n(j) * vn_ext(nbr)
                end do
            end do
        end if

        !Computing thermal field distribution...
        q = sp%lhf * (vf - vf_old) / sp%dt
        call trnsprt_no_upwind(u, krnls, ps%totpnts, mat%d, vf, vel, q, sp%dt)

        !Computing the local interfacial normal velocity...
        call get_gt_pure(sp%Teq, st_func, d0, a, b, theta0, crv_chk, mob_chk, ip, icrv, ivn, totip, T_intrf)
        call get_gt_vn(ip, T_intrf, np_sol, np_liq, totip, sp%lhf, K_sol, K_liq, d_nrm, u, ps%pnts, ikrnls, vn)

        !Projecting the local interfacial velocity onto the fixed domain points close to the interface...
        call get_vn_base(1*ps%dx, ip, vn, totip, ps%pnts, ps%totpnts, ikrnls, prj_pntnums, prj_pntvals, totprjpnts)
        
        !Solving for the extended interfacial velocity at the narrow band...
        call get_vn_ext(krnls, ps%totpnts, prj_pntnums, prj_pntvals, totprjpnts, sp%dt, 500, vn_ext)
        vn_ext = -vn_ext !Flip the sign of velocity when using vof instead of sdf

        vf_old = vf;
        call vof(vf, vel0, vn_ext, krnls, ps%totpnts, sp%dt, 5)
        !if (sp%sharp_intrf .eqv. .true.) then
        !    if(sp%intrf_slvr <= 3 .or. sp%intrf_slvr == 5) then
        !        call shrpn_intrf(phi, krnls, ps%totpnts, .true., sp%si_beta, sp%si_eps, sp%dx, sp%si_dt, sp%si_nt)
        !    else
        !        call shrpn_intrf(phi, krnls, ps%totpnts, .false., sp%si_beta, sp%si_eps, sp%dx, sp%si_dt, sp%si_nt)
        !    end if
        !end if

        do i = 1, ps%totpnts
            if (vf(i) > 1.0) vf(i) = 1.0
            if (vf(i) < 0.0) vf(i) = 0.0
        end do
        phi = -(2.0 * vf - 1.0)

        !Printing solutions to file
        if (t/sp%prnt_frq == c) then
            if (sp%vtk .eqv. .true.) then
                call prnt_vtk(u, ps, u_fname, t)
                call prnt_vtk(vn_ext, ps, vn_fname, t)
                select case (ps%dim)
                    case(3)
                        call prnt_unstruct_vtk(vf, ps, phi_fname, t)
                    case default
                        call prnt_vtk(vf, ps, phi_fname, t)
                end select
            else
                call prnt_txt(u, ps%pnts, ps%totpnts, u_fname, t)
                call prnt_txt(vf, ps%pnts, ps%totpnts, phi_fname, t)
            end if

            c = c + 1
        end if
       !deallocate(uip)
       deallocate(icrv)
       deallocate(ivn)
       deallocate(crv)
       deallocate(vn)
       deallocate(ikrnls)
    end do

end subroutine slv_stfn

subroutine get_nrm_vec(krnls, totip, phi0, dphi)
    use krnl_struct

    real(8), dimension(:), intent(in) :: phi0
    integer, intent(in) :: totip
    type (kernel), dimension(:), intent(in)   :: krnls
    real(8), dimension(:,:), allocatable, intent(out) :: dphi
    real(8) :: dphidx, dphidy, dphidz, absdphi
    integer :: i, mm

    allocate(dphi(totip, 4))
    do i = 1, totip
        ! computing normal vector
        dphidx = 0.0
        dphidy = 0.0
        dphidz = 0.0
        do mm = 1, krnls(i)%totnbrs
            dphidx = dphidx + krnls(i)%nx(mm) * phi0(krnls(i)%nbrs(mm))
            dphidy = dphidy + krnls(i)%ny(mm) * phi0(krnls(i)%nbrs(mm))
            dphidz = dphidz + krnls(i)%nz(mm) * phi0(krnls(i)%nbrs(mm))
        end do
        absdphi = (1e-6)+(dphidx*dphidx + dphidy*dphidy + dphidz*dphidz)**0.5
        dphi(i, 1) = dphidx
        dphi(i, 2) = dphidy
        dphi(i, 3) = dphidz
        dphi(i, 4) = absdphi
    end do
end subroutine get_nrm_vec

subroutine get_nrm_pnts(ip, totip, delta, phi0, krnls, np_sol, np_liq)
    use krnl_struct

    integer, intent(in) :: totip
    real(8), dimension(:,:), intent(in) :: ip
    real(8), dimension(:), intent(in) :: phi0
    real(8), intent(in) :: delta
    type(kernel), dimension(:), intent(in) :: krnls
    real(8), dimension(:,:), allocatable, intent(out) :: np_sol, np_liq
    real(8), dimension(:,:), allocatable :: dphi
    integer :: i
    real(8) :: dphidx, dphidy, dphidz, absdphi

    allocate(dphi(totip, 4))
    allocate(np_sol(totip, 3))
    allocate(np_liq(totip, 3))
    call get_nrm_vec(krnls, totip, phi0, dphi)
    do i = 1, totip
        !generate liquid and solid points normal to the interface point
        dphidx = dphi(i,1)
        dphidy = dphi(i,2)
        dphidz = dphi(i,3)
        absdphi = dphi(i,4)

        np_sol(i,1) = ip(i,1) - (dphidx / absdphi) * delta
        np_sol(i,2) = ip(i,2) - (dphidy / absdphi) * delta
        np_sol(i,3) = ip(i,3) - (dphidz / absdphi) * delta

        !dx = ip(i,1) - np_sol(i,1)
        !dy = ip(i,2) - np_sol(i,2)
        !dz = ip(i,3) - np_sol(i,3)
        !d = (dx*dx + dy*dy + dz*dz)**0.5
        !np_sol(i,4) = -d

        np_liq(i,1) = ip(i,1) + (dphidx / absdphi) * delta
        np_liq(i,2) = ip(i,2) + (dphidy / absdphi) * delta
        np_liq(i,3) = ip(i,3) + (dphidz / absdphi) * delta
        
        !dx = ip(i,1) - np_liq(i,1)
        !dy = ip(i,2) - np_liq(i,2)
        !dz = ip(i,3) - np_liq(i,3)
        !d = (dx*dx + dy*dy + dz*dz)**0.5
        !np_liq(i,4) = d
    end do
end subroutine get_nrm_pnts

subroutine get_gt_pure(T_eq, st_func, d0, a, b, theta0, crv_chk, mob_chk, ip, crv, vn, totip, T_intrf)

    real(8), intent(in) :: T_eq
    integer, intent(in) :: st_func, totip
    real(8), intent(in) :: d0, a, b, theta0
    logical, intent(in) :: crv_chk, mob_chk
    real(8), dimension(:,:), intent(in) :: ip
    real(8), dimension(:), intent(in) :: crv, vn
    real(8), dimension(:), allocatable, intent(out) :: T_intrf
    integer :: i
    real(8) :: theta, ec, trm, curv, ev

    allocate(T_intrf(totip))
    do i = 1, totip 
        theta = atan(ip(i,2)/ip(i,1))
        
        ec = d0
        if (st_func == 1) then
            ec = d0 * (1 - a * cos(b * theta + theta0))
        else if (st_func == 2) then
            trm = 0.5 * b * (theta - theta0)
            trm = sin(trm)**4
            ec = d0 * (a * trm)
        end if

        curv = 1.0
        if (crv_chk .eqv. .true.) then
            curv = crv(i)
        end if
        
        ev = 0
        if (mob_chk .eqv. .true.) then
            ev = ec
        end if
        T_intrf(i) = T_eq - ec * curv - ev * vn(i)
    end do
end subroutine get_gt_pure

subroutine get_gt_vn(ip, T_intrf, np_sol, np_liq, totip, lhf, K_sol, K_liq, dx, T, pnts0, ikrnls, vn)
    use krnl_struct
    
    real(8), dimension(:,:), intent(in) :: ip, np_sol, np_liq, pnts0
    real(8), dimension(:), intent(in) :: T_intrf, T
    real(8), intent(in) :: lhf, dx
    integer, intent(in) :: totip
    real(8), intent(in) :: K_sol, K_liq
    real(8), dimension(:), allocatable, intent(out) :: vn
    type(kernel), dimension(:), intent(in) :: ikrnls
    real(8) :: T_sol_nrm, T_liq_nrm
    integer :: i

    allocate(vn(totip))
    !Computing local normal interface velocity
    do i = 1, totip
        call get_nrm_var(np_sol(i,:), ip, totip, dx, T_intrf, T, pnts0, &
                ikrnls(i)%nbrs, ikrnls(i)%totnbrs, T_sol_nrm)
        call get_nrm_var(np_liq(i,:), ip, totip, dx, T_intrf, T, pnts0, &
                ikrnls(i)%nbrs, ikrnls(i)%totnbrs, T_liq_nrm)
        vn(i) = (1.0 / lhf) * ((K_liq * (T_intrf(i) - T_liq_nrm) / dx) &
                - (K_sol * (T_sol_nrm - T_intrf(i)) / dx))
    end do
end subroutine get_gt_vn

subroutine get_nrm_var(pnt, ip, totip, dxx, C_eq, C, pnts0, nbrs, totnbrs, C_intrp)
    real(8), dimension(:), intent(in) :: pnt
    real(8), dimension(:,:), intent(in) :: ip, pnts0
    integer, intent(in) :: totip, totnbrs
    real(8), intent(in) :: dxx
    real(8), dimension(:), intent(in) :: C, C_eq
    integer, dimension(:), intent(in) :: nbrs
    real(8), intent(out) :: C_intrp
    real(8) :: sum_w, sum_C, w, dx, dy, dz, d
    integer :: mm, nbr

    sum_w = 0
    sum_C = 0
    do mm = 1, totnbrs
        nbr = nbrs(mm)
        dx = pnt(1) - pnts0(nbr, 1)
        dy = pnt(2) - pnts0(nbr, 2)
        dz = pnt(3) - pnts0(nbr, 3)
        d = 1e-6 + (dx*dx + dy*dy + dz*dz)**0.5
        w = 1.0 / d
        sum_w = sum_w + w;
        sum_C = sum_C + C(nbr) * w
    end do

    do mm = 1, totip
        dx = pnt(1) - ip(mm, 1)
        dy = pnt(2) - ip(mm, 2)
        dz = pnt(3) - ip(mm, 3)
        d = 1e-6 + (dx*dx + dy*dy + dz*dz)**0.5
        if (d <= 3 * dxx) then
            w = 1 / d
            sum_w = sum_w + w
            sum_C = sum_C + C_eq(mm) * w
        end if
    end do
    C_intrp = sum_C / sum_w
end subroutine get_nrm_var

subroutine get_vn_base(vn_ext_rad, ip, vn, totip, pnts0, totpnts0, ikrnls, prj_pntnums, prj_pntvals, totprjpnts)
    use krnl_struct

    integer, dimension(:), allocatable :: marker
    type(kernel), dimension(:), intent(in) :: ikrnls
    real(8), intent(in) :: vn_ext_rad
    real(8), dimension(:,:), intent(in) :: ip, pnts0
    integer, intent(in):: totip, totpnts0
    real(8), dimension(:), intent(in) :: vn
    integer, dimension(:), allocatable, intent(out) :: prj_pntnums
    real(8), dimension(:), allocatable, intent(out) :: prj_pntvals
    integer, intent(out) :: totprjpnts
    integer :: i, j, nbr, pntnum
    real(8) :: dx, dy, dz, d, w, sum_w, sum_vn

    allocate(prj_pntnums(0))
    totprjpnts = 0
		
    allocate(marker(totpnts0))
    marker(:) = 0

    do i = 1, totip
        do j = 1, ikrnls(i)%totnbrs
            nbr = ikrnls(i)%nbrs(j)
            dx = ip(i,1) - pnts0(nbr,1)
            dy = ip(i,2) - pnts0(nbr,2)
            dz = ip(i,3) - pnts0(nbr,3)
            d = (dx * dx + dy * dy + dz * dz)**0.5
            if (d <= vn_ext_rad) then
                marker(nbr) = 1
            end if
        end do
    end do
    do i = 1, totpnts0
        if(marker(i) == 1) then
            totprjpnts = totprjpnts + 1
            prj_pntnums = [prj_pntnums, i]
        end if
    end do
		
    allocate(prj_pntvals(totprjpnts))
    do i = 1, totprjpnts
        sum_vn = 0
        sum_w = 0
        pntnum = prj_pntnums(i)
        do j = 1, totip
            dx = ip(j,1) - pnts0(pntnum,1)
            dy = ip(j,2) - pnts0(pntnum,2)
            dz = ip(j,3) - pnts0(pntnum,3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d >= 1e-6 .and. d <= vn_ext_rad) then
                w = 1 / d
                sum_vn = sum_vn + w * vn(j)
                sum_w = sum_w + w
            end if
        end do
        prj_pntvals(i) = sum_vn / sum_w
    end do
end subroutine

subroutine get_vn_ext(krnls, totpnts, prj_pntnums, prj_pntvals, totprjpnts, dt, vnt, vn_ext)
    use krnl_struct

    integer, intent(in) :: totpnts
    type(kernel), dimension(:), intent(in) :: krnls
    integer, dimension(:), intent(in) :: prj_pntnums
    real(8), dimension(:), intent(in) :: prj_pntvals
    integer, intent(in) :: totprjpnts
    real(8), intent(in) :: dt
    integer, intent(in) :: vnt
    integer :: vt, i, j, nbr 
    real(8) :: diff = 0.0
    real(8), dimension(:), allocatable :: vn_fut
    real(8), dimension(:), allocatable, intent(out) :: vn_ext
    
    allocate(vn_ext(totpnts))
    allocate(vn_fut(totpnts))
    vn_ext(:) = 0.0
    vn_fut = vn_ext;
    do vt = 1, vnt
        do i = 1, totprjpnts
            vn_ext(prj_pntnums(i)) = prj_pntvals(i)
        end do

        do i = 1, totpnts
            diff = 0;
            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                diff = diff + (vn_ext(nbr) - vn_ext(i))* krnls(i)%nabla2(j)
            end do
            vn_fut(i) = vn_ext(i) + dt * diff
        end do
        vn_ext = vn_fut
    end do
end subroutine get_vn_ext

end module stfn_slvrs
