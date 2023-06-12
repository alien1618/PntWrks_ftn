module ns_2phs_slvrs
implicit none
contains

subroutine run_ns_2phs()
!------------------------------------------------------------------------------------------
! This subroutine solves the eulerian two-phase navier-stokes equations
! local strong-form meshfree methods and iterative explicit time stepping schemes     
!------------------------------------------------------------------------------------------
    use pntst_struct
    use bc_struct
    use slvr_prmtrs_struct
    use mat_struct
    use pntst
    use bndry
    use intrf
    use ns
    use prmtrs
    implicit none

    real(8), dimension(:,:), allocatable       :: v        !velocity
    real(8), dimension(:), allocatable         :: p        !pressure
    real(8), dimension(:,:), allocatable       :: f        !force
    real(8), dimension(:), allocatable         :: phi      !phase order parameter
    type(pointset)                             :: ps       !pointset data structure
    type(bc), dimension(6)                     :: bcs      !boundary conditions data structure for vx, vy, vz, p
    type(slvr_prmtrs)                          :: sp       !solver parameters data structure
    type(materials)                            :: mat      !phase properties
    real*8                                     :: start, finish !timer parameters   
!------------------------------------------------------------------------------------------
    !constructing the pointset
    call set_geom(ps)   
!------------------------------------------------------------------------------------------
    !constructing the implicit geometry
    call set_phi(ps%dx, sp)
    call set_phi0(ps%pnts, ps%totpnts, 0.5*ps%dx, phi)
    if (sp%invert_phi == 1) phi = 1 - phi  
!------------------------------------------------------------------------------------------
    !assigning the material properties to each phase
    mat%total = 2
    call set_mat('k', mat%d)
    call set_mat('nu', mat%nu)
    call set_mat('ro', mat%ro)
    call set_mat('g', sp%g)
!------------------------------------------------------------------------------------------
    !initializing the field variables
    !call set_v(sp)
    call set_v0(ps%pnts, phi, ps%totpnts, v)                            
    allocate(p(ps%totpnts))
    allocate(f(ps%totpnts,3))
    p(:) = 0
    f(:,1) = sp%g(1)
    f(:,2) = sp%g(2)
    f(:,3) = sp%g(3)
!------------------------------------------------------------------------------------------
    !assigning the boundary conditions
    call set_bc(ps%pnts, ps%totpnts, bcs)
    call set_bc_phi(phi, ps%pnts, ps%totpnts, bcs)
!------------------------------------------------------------------------------------------
    !assigning the solver control parameters
    call set_krnl(ps%dx, sp)
    call set_t(sp)
    sp%p = 1.0 / sp%dt        !pressure prmtr
    sp%lgr = .false.            !0 = eulerian, 1 = lagrangian
    sp%av = 0.1                 !artificial viscosity parameter
    write(*,'(a)') "Running solver..."
    call cpu_time(start)
    call slv_ns_2phs(ps, v, p, f, phi, mat, bcs, sp)
    call cpu_time(finish)
    print '("Computation time = ",f8.1," seconds.")',finish-start
end subroutine run_ns_2phs

subroutine slv_ns_2phs(ps, vel1, P1, F, pntsphi, mat, bcs, sp)

    use krnl_struct
    use pntst_struct
    use bc_struct
    use slvr_prmtrs_struct
    use mat_struct
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_mls
    use krnl_wls
    use krnl_krg
    use pntst
    use bndry
    use ns
    use intrf

    implicit none

    type (kernel),dimension(:), allocatable :: krnls
    integer                                 :: t, c = 1, cc = 1
    real(8), dimension(:), allocatable      :: v, p, phi, phi_old, ro, vn
    real(8), dimension(:,:), allocatable    :: vel
    real(8), dimension(:,:), intent(in)     :: F
    real(8), dimension(:), intent(inout)    :: p1, pntsphi
    real(8), dimension(:,:), intent(inout)  :: vel1
    real(8)                                 :: mobility = 1, m, cfl
    character(len=50)                       :: v_fname, p_fname, u_fname, ro_fname, phi_fname
    integer                                 :: i, v_nt = 1000
    real(8)                                 :: sum_vof0, w, dt
    type(pointset), intent(in)              :: ps
    type(bc), dimension(5), intent(in)      :: bcs
    type(slvr_prmtrs), intent(in)           :: sp
    type(materials), intent(in)             :: mat
    real(8)                                 ::  start, finish
!------------------------------------------------------------------------------------------
    v_fname = 'v'
    p_fname = 'p'
    u_fname = 'u'
    ro_fname = 'ro'
    phi_fname = 'phi'

    write(*,'(a)') "Initializing field variable distribution..."
    allocate(phi(ps%totpnts))
    allocate(phi_old(ps%totpnts))
    allocate(v(ps%totpnts))
    allocate(ro(ps%totpnts))
    allocate(vel(ps%totpnts,3))
    allocate(P(ps%totpnts))

    W = ps%dx
    M = mobility*W*W
    dt = sp%dt
    vel = vel1
    call set_var(vel(:,1), bcs(1))
    call set_var(vel(:,2), bcs(2))
    call set_var(vel(:,3), bcs(3))

    allocate(vn(ps%totpnts))
    vn(:) = 0

    sum_vof0 = 0
    do i = 1, ps%totpnts
        sum_vof0 = sum_vof0+phi(i)
    end do

    ! intrf solver method: 1 = vof, 2 = vof2, 3 = AC, 4 = AC2, 5 = CH, 6 = CH2
    if (sp%intrf_slvr == 4 .or. sp%intrf_slvr == 6) then
        !pointset.points.phi is given as vof (0,1), convert to heaveside function for intrf slvrs
        phi = 2.0*pntsphi-1.0
    else
        phi = pntsphi
    end if

    do i = 1,ps%totpnts
        ro(i) = pntsphi(i)*mat%ro(2) + (1-pntsphi(i))*mat%ro(1)
    end do

    p(:) = 0
    v = sqrt(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2)+vel(:,3)*vel(:,3))

!------------------------------------------------------------------------------------------
    write(*,'(a)') "Constructing krnls..."
    call cpu_time(start)
    allocate(krnls(ps%totpnts))
    call get_nbrs_bf(sp%h, ps%pnts, ps%totpnts, ps%pnts, ps%totpnts, krnls)
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
    write(*,'(a)') "Computing transient solution..."
    if (sp%vtk .eqv. .true.) then
        select case (ps%dim)
        case (3)
            call prnt_unstruct_vtk(V, ps, v_fname, 0)
            call prnt_unstruct_vtk(p, ps, p_fname, 0)
            call prnt_unstruct_vtk(ro, ps, ro_fname, 0)
            call prnt_unstruct_vtk(pntsphi, ps, phi_fname, 0)
        case default
            call prnt_vtk(V, ps, v_fname, 0)
            call prnt_vtk(p, ps, p_fname, 0)
            call prnt_vtk(ro, ps, ro_fname, 0)
            call prnt_vtk(pntsphi, ps, phi_fname, 0)
        end select
    else
        call prnt_pltctrl(ps%pnts, sp%nt, sp%prnt_frq)
        call prnt_elems_txt(ps%surfs, ps%totsurfs)
        call prnt_txt(v, ps%pnts, ps%totpnts, v_fname, 0)
        call prnt_txt(p, ps%pnts, ps%totpnts, p_fname, 0)
        call prnt_txt(ro, ps%pnts, ps%totpnts, ro_fname, 0)
        call prnt_txt(phi, ps%pnts, ps%totpnts, phi_fname, 0)
    end if

    do t = 1, sp%nt
        ! cout << "Computing CFL condition..." << endl;
        v = sqrt(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2)+vel(:,3)*vel(:,3))
        CFL = (maxval(v)*dt)/ps%dx
        if (CFL >= 0.001) then
            dt = dt*0.5
        end if

        !"Solving the Navier-Stokes equations..." << endl;
        call slv_ns(vel,V,P,ro, ps%pnts, ps%totpnts, krnls, vel1, P1, &
        F, pntsphi, mat%nu, mat%ro, mat%total, sp%segma, bcs, &
         dt, sp%av, sp%p, sp%lgr)

        call set_var(vel(:,1), bcs(1))
        call set_var(vel(:,2), bcs(2))
        call set_var(vel(:,3), bcs(3))

        vel1 = vel
        P1 = P

        if(t/V_Nt == cc) then
            !cout << "Solving the phase-field/volume of fluid equations..." << endl;
            phi_old = phi
            !cout << "method " << method << endl;
            select case (sp%intrf_slvr)
                case (1)
                    call vof(phi, vel, vn, krnls, ps%totpnts, V_Nt*dt, sp%itrs)
                case (2)
                    call vof2(phi, vel, vn, krnls, ps%totpnts, V_Nt*dt, sp%itrs)
                case (3)
                    call ac(phi, vel, vn, krnls, ps%totpnts, w, m, V_Nt*dt, sp%itrs)
                case (4)
                    call ac2(phi, vel, vn, krnls, ps%totpnts, m, V_Nt*dt, sp%itrs)
                case (5)
                    call ch(phi, vel, vn, krnls, ps%totpnts, w, m, V_Nt*dt, sp%itrs)
                case (6)    
                    call ch2(phi, vel, vn, krnls, ps%totpnts, w, m, V_Nt*dt, sp%itrs)
            end select
            if (sp%sharp_intrf .eqv. .true.) then
                if(sp%intrf_slvr <= 3 .or. sp%intrf_slvr == 5) then
                    call shrpn_intrf(phi, krnls, ps%totpnts, .true., sp%si_beta, sp%si_eps, ps%dx, sp%si_dt, sp%si_nt)
                else
                    call shrpn_intrf(phi, krnls, ps%totpnts, .false., sp%si_beta, sp%si_eps, ps%dx, sp%si_dt, sp%si_nt)
                end if
           end if

            call set_var(phi, bcs(5))

            if (sp%intrf_slvr == 4 .or. sp%intrf_slvr == 6) then
                do i = 1, ps%totpnts
                    if (phi(i) < -1.0) then
                        phi(i) = -1
                    end if
                    if (phi(i) > 1.0) then
                        phi(i) = 1
                    end if
                end do
            else
                do i = 1, ps%totpnts
                    if (phi(i) < 0) then
                        phi(i) = 0
                    end if
                    if (phi(i) > 1.0) then
                        phi(i) = 1.0
                    end if
                end do
            end if

            if (sp%intrf_slvr == 4 .or. sp%intrf_slvr == 6) then
                !convert heavside function (-1,1) to VOF (0,1) for the NS solver
                pntsphi = 0.5*(phi+1.0)
            else
                pntsphi = phi
            end if
            cc = cc + 1
        end if

        if (t/sp%prnt_frq == c) then
            c = c + 1
            write(*,'(a, i0, a, i0, a)') "Computation of time step ", t, " of ", sp%nt, " complete..."
            if (sp%vtk .eqv. .true.) then
                select case (ps%dim)
                case (3)
                    call prnt_unstruct_vtk(V, ps, v_fname, t)
                    call prnt_unstruct_vtk(p, ps, p_fname, t)
                    call prnt_unstruct_vtk(ro, ps, ro_fname, t)
                    call prnt_unstruct_vtk(phi, ps, phi_fname, t)
                case default
                    call prnt_vtk(V, ps, v_fname, t)
                    call prnt_vtk(p, ps, p_fname, t)
                    call prnt_vtk(ro, ps, ro_fname, t)
                    call prnt_vtk(phi, ps, phi_fname, t)
                end select
            else
                call prnt_txt(v, ps%pnts, ps%totpnts, v_fname, t)
                call prnt_txt(p, ps%pnts, ps%totpnts, p_fname, t)
                call prnt_txt(ro, ps%pnts, ps%totpnts, ro_fname, t)
                call prnt_txt(phi, ps%pnts, ps%totpnts, phi_fname, t)
            end if
        end if
    end do
end subroutine slv_ns_2phs

end module ns_2phs_slvrs
