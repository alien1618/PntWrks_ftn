module ac_slvrs
implicit none
contains

subroutine run_ac ()
!------------------------------------------------------------------------------------------
!   This subroutine solves the allen-cahn phase-field equations for dendritic
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
    write(*,'(a)') "Running Allen-Cahn solver..."
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
    call slv_ac_singleseed_2d_prv(ps, sp, mat, bcs, u, phi, v)
    finish = omp_get_wtime()
    print '("Time = ",f8.1," seconds.")',finish-start
end subroutine run_ac

subroutine slv_ac_singleseed_2d(ps, sp, mat, bcs, u, phi, vel0)
!-------------------------------------------------------------------------
    use krnl_struct
    use slvr_prmtrs_struct
    use pntst_struct
    use mat_struct
    use bc_struct
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_wls
    use krnl_mls
    use krnl_krg
    use intrf
    use pntst
    use bndry
    use trnsprt
    use ns
!------------------------------------------------------------------------------------------ 
    type(pointset), intent(in)              ::  ps
    type(bc), dimension(4), intent(in)      ::  bcs
    type(slvr_prmtrs), intent(in)           ::  sp 
    type(materials), intent(in)             ::  mat
!------------------------------------------------------------------------------------------ 
    real(8), dimension(:), intent(inout)    ::  u, phi
    real(8), dimension(:,:), intent(inout)  ::  vel0
!------------------------------------------------------------------------------------------
    real(8), dimension(ps%totpnts)          :: phi_old, q
    real(8), dimension(:,:), allocatable    :: vel, f
    integer                                 :: i, c = 1, t
    logical                                 :: fluid_slvr = .false.
    real(8)                                 :: start, finish
    real(8), dimension(:), allocatable      :: v, p, p0, ro
    integer, dimension(:), allocatable      :: pntbins
    integer                                 :: nx, nxy, totbins
    type (kernel), dimension(:), allocatable :: bins
    type (kernel), dimension(ps%totpnts)     :: krnls
    type (kernel), dimension(:), allocatable :: krnls2
    character(len=50)                        :: u_fname = 'u', phi_fname = 'phi'
    character(len=50)                        :: v_fname = 'v', p_fname = 'p'
!------------------------------------------------------------------------------------------  
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
    write(*,'(a)') "Constructing krnls..."
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

    if (fluid_slvr .eqv. .true.) then
        allocate(krnls2(ps%totpnts))
        call get_nbrs_bf(sp%h, ps%pnts, ps%totpnts, ps%pnts, ps%totpnts, krnls2) !find a better way instead of doing this twice
        call get_sph_krnls(ps%pnts, ps%totpnts, ps%pnts, sp%sph, ps%dim, sp%h,  krnls2)
        call get_intrps_o2(ps%pnts, ps%totpnts, krnls2)
    end if
    print '("interpolants construction time = ",f10.1," seconds")',finish-start
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Computing transient solution..."
    if (sp%vtk .eqv. .true.) then
        call prnt_vtk(u, ps, u_fname, 0)
        select case (ps%dim)
            case(3)
                call prnt_unstruct_vtk(phi, ps, phi_fname, 0)
            case default
                call prnt_vtk(phi, ps, phi_fname, 0)
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
        call prnt_txt(phi, ps%pnts, ps%totpnts, phi_fname, 0)
    end if

    do t = 1,sp%nt 
        if (fluid_slvr .eqv. .true.) then
            call slv_ns(vel, v, p, ro, ps%pnts, ps%totpnts, krnls2, vel0, p0, &
            F, phi, mat%nu, mat%ro, mat%total, sp%segma, bcs, &
            sp%dt, sp%av, sp%p, sp%lgr)
            vel(bcs(1)%pnts,1) = bcs(1)%vals
            vel(bcs(2)%pnts,2) = bcs(2)%vals
            vel(bcs(3)%pnts,3) = bcs(3)%vals
            do i = 1, ps%totpnts
                if (phi(i) >= 0.5) vel(i,:) = 0
            end do
            vel0 = vel
            p0 = p
        end if 

        phi_old = phi
        call ac_phstrans(phi, u, krnls, ps%totpnts, sp, sp%dt)
        q = sp%lhf * (phi - phi_old) / sp%dt
        call trnsprt_no_upwind(u, krnls, ps%totpnts, mat%d, phi, vel, q, sp%dt)
        do i = 1, ps%totpnts
            if (phi(i) > 1) phi(i) = 1
            if (phi(i) < 0) phi(i) = 0
        end do
        if (t/sp%prnt_frq == c) then
            if (sp%vtk .eqv. .true.) then
                call prnt_vtk(u, ps, u_fname, t)
                select case (ps%dim)
                    case(3)
                        call prnt_unstruct_vtk(phi, ps, phi_fname, t)
                    case default
                        call prnt_vtk(phi, ps, phi_fname, t)
                end select
                if (fluid_slvr .eqv. .true.) then
                    v = sqrt(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2))
                    call prnt_vtk(v, ps, v_fname, t)
                    call prnt_vtk(p, ps, p_fname, t)
                end if
            else
                call prnt_txt(u, ps%pnts, ps%totpnts, u_fname, t)
                call prnt_txt(phi, ps%pnts, ps%totpnts, phi_fname, t)
            end if

            c = c + 1
        end if
    end do

end subroutine slv_ac_singleseed_2d

subroutine ac_phstrans(phi, T1, krnls, totpnts, pf, dt)
!-----------------------------------------------------------
! solve allen-cahn equation for solidification/melting in a 
! pure liquid/solid
!-----------------------------------------------------------
    use krnl_struct
    use slvr_prmtrs_struct
    use omp_lib
!-----------------------------------------------------------
    real(8), dimension(:), intent(in)          :: T1
    real(8), intent(in)                        :: dt
    integer, intent(in)                        :: totpnts
    type(kernel), dimension(:), intent(in)     :: krnls
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)       :: phi
!-----------------------------------------------------------
    real(8), dimension(totpnts)                :: eps0, eps_trmx, eps_trmy, phi_fut
    real(8)                                    :: eps_deriv, eps_trmx_dy, eps_trmy_dx
    real(8)                                    :: dphidx, dphidy, nabla2_phi, Q, m, theta
    integer                                    :: i, j, nbr
    type(slvr_prmtrs)                          :: pf
!-----------------------------------------------------------
   
    !$omp parallel do private(i)
    do i = 1, totpnts
        dphidx = 0
        dphidy = 0
        do j = 1,krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            dphidx = dphidx + krnls(i)%nx(j)*(phi(nbr)-phi(i))
            dphidy = dphidy + krnls(i)%ny(j)*(phi(nbr)-phi(i))
        end do
        theta = atan2(dphidy, dphidx)
        eps0(i) = pf%eps_4 * (1.0 + pf%delta * cos(pf%a_2 * (theta - pf%theta0)))
        eps_deriv = - pf%eps_4 * pf%a_2 * pf%delta * sin(pf%a_2 * (theta - pf%theta0))
        eps_trmx(i) = eps0(i) * eps_deriv * dphidx
        eps_trmy(i) = eps0(i) * eps_deriv * dphidy
    end do
    !$omp end parallel do

    !$omp parallel do private(i)
    do i = 1, totpnts
        eps_trmx_dy = 0
        eps_trmy_dx = 0
        nabla2_phi = 0
        do j = 1,krnls(i)%totnbrs
            nbr =  krnls(i)%nbrs(j)
            eps_trmx_dy = eps_trmx_dy + krnls(i)%ny(j)*(eps_trmx(nbr)-eps_trmx(i))
            eps_trmy_dx = eps_trmy_dx + krnls(i)%nx(j)*(eps_trmy(nbr)-eps_trmy(i))
            nabla2_phi = nabla2_phi + (phi(nbr)-phi(i))*krnls(i)%nabla2(j)
        end do
        m = pf%alpha / 3.14 * atan(pf%gamma * (pf%Teq - T1(i)))
        Q = eps_trmx_dy - eps_trmy_dx + phi(i) * (1.0 - phi(i)) * (phi(i) - 0.5 + m)
        phi_fut(i) = phi(i) + (dt / pf%tau0) * (eps0(i)**2) * nabla2_phi + (dt / pf%tau0) * Q

        if(phi_fut(i) > 1) phi_fut(i) = 1
        if(phi_fut(i) < 0) phi_fut(i) = 0
    end do
    !$omp end parallel do
    phi = phi_fut
end subroutine ac_phstrans

subroutine slv_ac_singleseed_2d_prv(ps, sp, mat, bcs, u, phi, vel0)
!-------------------------------------------------------------------------
    use krnl_struct
    use slvr_prmtrs_struct
    use pntst_struct
    use mat_struct
    use bc_struct
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_wls
    use krnl_mls
    use krnl_krg
    use intrf
    use pntst
    use bndry
    use trnsprt
    use ns
!------------------------------------------------------------------------------------------ 
    type(pointset), intent(in)              ::  ps
    type(bc), dimension(4), intent(in)      ::  bcs
    type(slvr_prmtrs), intent(in)           ::  sp 
    type(materials), intent(in)             ::  mat
!------------------------------------------------------------------------------------------ 
    real(8), dimension(:), intent(inout)    ::  u, phi
    real(8), dimension(:,:), intent(inout)  ::  vel0
!------------------------------------------------------------------------------------------
    real(8), dimension(ps%totpnts)          :: phi_old, q
    real(8), dimension(:,:), allocatable    :: vel, f
    integer                                 :: i, c = 1, t
    logical                                 :: fluid_slvr = .false.
    real(8)                                 :: start, finish
    real(8), dimension(:), allocatable      :: v, p, p0, ro
    integer, dimension(:), allocatable      :: pntbins
    integer                                 :: nx, nxy, totbins
    type (kernel), dimension(:), allocatable :: bins
    type (kernel), dimension(ps%totpnts)     :: krnls
    type (kernel), dimension(:), allocatable :: krnls2
    character(len=50)                        :: u_fname = 'u', phi_fname = 'phi'
    character(len=50)                        :: v_fname = 'v', p_fname = 'p'
!------------------------------------------------------------------------------------------  
    phi = (2.0 * phi - 1.0)
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
    write(*,'(a)') "Constructing krnls..."
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

    if (fluid_slvr .eqv. .true.) then
        allocate(krnls2(ps%totpnts))
        call get_nbrs_bf(sp%h, ps%pnts, ps%totpnts, ps%pnts, ps%totpnts, krnls2) !find a better way instead of doing this twice
        call get_sph_krnls(ps%pnts, ps%totpnts, ps%pnts, sp%sph, ps%dim, sp%h,  krnls2)
        call get_intrps_o2(ps%pnts, ps%totpnts, krnls2)
    end if
    print '("interpolants construction time = ",f10.1," seconds")',finish-start
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Computing transient solution..."
    if (sp%vtk .eqv. .true.) then
        call prnt_vtk(u, ps, u_fname, 0)
        select case (ps%dim)
            case(3)
                call prnt_unstruct_vtk(phi, ps, phi_fname, 0)
            case default
                call prnt_vtk(phi, ps, phi_fname, 0)
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
        call prnt_txt(phi, ps%pnts, ps%totpnts, phi_fname, 0)
    end if

    do t = 1,sp%nt 
        if (fluid_slvr .eqv. .true.) then
            call slv_ns(vel, v, p, ro, ps%pnts, ps%totpnts, krnls2, vel0, p0, &
            F, 0.5*(phi+1.0), mat%nu, mat%ro, mat%total, sp%segma, bcs, &
            sp%dt, sp%av, sp%p, sp%lgr)
            vel(bcs(1)%pnts,1) = bcs(1)%vals
            vel(bcs(2)%pnts,2) = bcs(2)%vals
            vel(bcs(3)%pnts,3) = bcs(3)%vals
            do i = 1, ps%totpnts
                if (0.5*phi(i)+1.0 >= 0.5) vel(i,:) = 0
            end do
            vel0 = vel
            p0 = p
        end if 

        phi_old = phi
        call ac_phstrans_prv(phi, u, krnls, ps%totpnts, sp, sp%dt)
        !q = sp%lhf * (0.5*(phi+1.0) - 0.5*(phi_old+1.0)) / sp%dt
        !q = sp%lhf * (phi - phi_old) / sp%dt
        !q = (phi - phi_old) / sp%dt
        q = 0.5 * (phi - phi_old)
        call trnsprt_no_upwind(u, krnls, ps%totpnts, 0.6267*mat%D, 0.5*(phi+1.0), vel, q, sp%dt)
        if (t/sp%prnt_frq == c) then
            if (sp%vtk .eqv. .true.) then
                call prnt_vtk(u, ps, u_fname, t)
                select case (ps%dim)
                    case(3)
                        call prnt_unstruct_vtk(phi, ps, phi_fname, t)
                    case default
                        call prnt_vtk(phi, ps, phi_fname, t)
                end select
                if (fluid_slvr .eqv. .true.) then
                    v = sqrt(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2))
                    call prnt_vtk(v, ps, v_fname, t)
                    call prnt_vtk(p, ps, p_fname, t)
                end if
            else
                call prnt_txt(u, ps%pnts, ps%totpnts, u_fname, t)
                call prnt_txt(phi, ps%pnts, ps%totpnts, phi_fname, t)
            end if

            c = c + 1
        end if
    end do

end subroutine slv_ac_singleseed_2d_prv

subroutine ac_phstrans_prv(phi, T1, krnls, totpnts, pf, dt)
!-----------------------------------------------------------
! solve allen-cahn equation for solidification/melting in a 
! pure liquid/solid
!-----------------------------------------------------------
    use krnl_struct
    use slvr_prmtrs_struct
    use omp_lib
!-----------------------------------------------------------
    real(8), dimension(:), intent(in)          :: T1
    real(8), intent(in)                        :: dt
    integer, intent(in)                        :: totpnts
    type(kernel), dimension(:), intent(in)     :: krnls
!-----------------------------------------------------------
    real(8), dimension(:), intent(inout)       :: phi
!-----------------------------------------------------------
    real(8), dimension(totpnts)                :: phi_fut
    real(8)                                    :: dphidx, dphidy, dphidz, absdphi
    real(8)                                    :: a, b, an
    real(8)                                    :: wn_phix, wn_phiy, wn_phiz
    real(8), dimension(totpnts)                :: tn, trmx, trmy, trmz, wn2_phix, wn2_phiy, wn2_phiz
    integer                                    :: i, j, nbr
    type(slvr_prmtrs)                          :: pf
    real(8)             :: dtrmx_dx, dtrmy_dy, dtrmz_dz, phi_xx, phi_yy, phi_zz
    real(8)             :: wn, u, Q, trm, diff, eps, trm0!, trm1 
!-----------------------------------------------------------
    eps = 1e-6

!    !$omp parallel do private(i)
    do i = 1, totpnts
        dphidx = 0
        dphidy = 0
        dphidz = 0
        do j = 1,krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            dphidx = dphidx + krnls(i)%nx(j)*(phi(nbr)-phi(i))
            dphidy = dphidy + krnls(i)%ny(j)*(phi(nbr)-phi(i))
            dphidz = dphidz + krnls(i)%nz(j)*(phi(nbr)-phi(i))
        end do
        absdphi = 1.0e-6 + (dphidx*dphidx + dphidy*dphidy + dphidz*dphidz)**0.5
        !---------------------------------------------------------------
        a = (1.0 - 3.0*pf%eps_4)
        b = (4.0*pf%eps_4)/a
        trm0 = (dphidx**4 + dphidy**4 + dphidz**4)/(absdphi**4)
        an = a * (1.0 + b * trm0)
        wn_phix = pf%w * a * b * ((4.0 * dphidx)/absdphi**2) * ((dphidx*dphidx/(absdphi**2)) - trm0)
        wn_phiy = pf%w * a * b * ((4.0 * dphidy)/absdphi**2) * ((dphidy*dphidy/(absdphi**2)) - trm0)
        wn_phiz = pf%w * a * b * ((4.0 * dphidz)/absdphi**2) * ((dphidz*dphidz/(absdphi**2)) - trm0)
        !---------------------------------------------------------------
        !trm0 = (dphidx**4 + dphidy**4 + dphidz**4)
        !trm1 = (dphidx**2 + dphidy**2 + dphidz**2)
        !an = 1.0 - 3.0 * pf%eps_4 - 4.0 * pf%eps_4 * trm0 / (eps + trm1**2) 
        !wn_phix = pf%w * 16 * pf%eps_4 * dphidx * ((trm0/(eps + trm1**3)) - (dphidx * dphidx / (eps + trm1**2)))
        !wn_phiy = pf%w * 16 * pf%eps_4 * dphidy * ((trm0/(eps + trm1**3)) - (dphidy * dphidy / (eps + trm1**2)))
        !wn_phiz = pf%w * 16 * pf%eps_4 * dphidz * ((trm0/(eps + trm1**3)) - (dphidz * dphidz / (eps + trm1**2)))
        !---------------------------------------------------------------
       
        wn = pf%w * an
        trmx(i) = (absdphi**2)*wn*wn_phix
        trmy(i) = (absdphi**2)*wn*wn_phiy
        trmz(i) = (absdphi**2)*wn*wn_phiz

        wn2_phix(i) = dphidx * wn**2 
        wn2_phiy(i) = dphidy * wn**2
        wn2_phiz(i) = dphidz * wn**2

        tn(i) = an * an + eps
    end do
!    !$omp end parallel do

!    !$omp parallel do private(i)
    do i = 1, totpnts
        dtrmx_dx = 0.0
        dtrmy_dy = 0.0
        dtrmz_dz = 0.0

        phi_xx = 0.0
        phi_yy = 0.0 
        phi_zz = 0.0
        do j = 1,krnls(i)%totnbrs
            nbr =  krnls(i)%nbrs(j)
            dtrmx_dx = dtrmx_dx + krnls(i)%nx(j)*(trmx(nbr)-trmx(i))
            dtrmy_dy = dtrmy_dy + krnls(i)%ny(j)*(trmy(nbr)-trmy(i))
            dtrmz_dz = dtrmz_dz + krnls(i)%nz(j)*(trmz(nbr)-trmz(i))

            phi_xx = phi_xx + krnls(i)%nx(j)*(wn2_phix(nbr)-wn2_phix(i))
            phi_yy = phi_yy + krnls(i)%ny(j)*(wn2_phiy(nbr)-wn2_phiy(i))
            phi_zz = phi_zz + krnls(i)%nz(j)*(wn2_phiz(nbr)-wn2_phiz(i))
        end do

        u = (T1(i) - pf%Teq)

        !Q = (phi(i) - pf%gamma *  u * (1 - phi(i)**2)) * (1-phi(i)**2)
        Q = phi(i) - (phi(i)**3) - pf%gamma * u * ((1.0-(phi(i)**2))**2)
        
        trm = (dtrmx_dx + dtrmy_dy + dtrmz_dz)
        
        diff = phi_xx + phi_yy + phi_zz

        phi_fut(i) = phi(i) + (dt / (pf%tau0*tn(i))) * (diff + trm + Q)

        if(phi_fut(i) > 1.0) phi_fut(i) = 1.0
        if(phi_fut(i) < -1.0) phi_fut(i) = -1.0
    end do
!    !$omp end parallel do
    phi = phi_fut
end subroutine ac_phstrans_prv

end module ac_slvrs
