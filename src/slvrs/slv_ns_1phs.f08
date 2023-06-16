module ns_1phs_slvrs
implicit none
contains

subroutine run_ns_1phs() 
!------------------------------------------------------------------------------------------
! This subroutine solves the eulerian single-phase navier-stokes equations
! local strong-form meshfree methods and iterative explicit time stepping schemes    
!------------------------------------------------------------------------------------------
    use pntst_struct
    use bc_struct
    use slvr_prmtrs_struct
    use mat_struct
    use pntst
    use bndry
    use prmtrs
!------------------------------------------------------------------------------------------
    real(8)                                    :: start, finish !timer parameters
    real(8), dimension(:,:), allocatable       :: v        !velocity
    real(8), dimension(:), allocatable         :: p        !pressure
    real(8), dimension(:), allocatable         :: phi      !implicit boundary indicator function
    real(8), dimension(:,:), allocatable       :: f        !force
    type(pointset)                             :: ps       !point set data structure
    type(bc), dimension(6)                     :: bcs      !boundary condition data structure for vx, vy, vz, p
    type(slvr_prmtrs)                          :: sp       !solver parameters data structures
    type(materials)                            :: mat      !phase properties data structures
!------------------------------------------------------------------------------------------
    write(*,'(a)') "--------------------------------------------------------"
    write(*,'(a)') "Running one-phase navier-stokes test..."
    write(*,'(a)') "--------------------------------------------------------"
    !constructing the pointset
    call set_geom(ps)
    !constructing the implicit geometry
    allocate(phi(ps%totpnts))
    phi(:) = 1.0
    !assigning the material properties to each phase
    call set_mat('k', mat%d)
    call set_mat('nu', mat%nu)
    call set_mat('ro', mat%ro)
    call set_mat('g', sp%g)
    mat%total = 1
    if (mat%d(1) > 0.0) then !rayleigh-benard convection
        call set_rb(sp)
        call set_phi0(ps%pnts, ps%totpnts, 0.5*ps%dx, phi)
    end if
    !initializing the field variables
    allocate(p(ps%totpnts))
    allocate(f(ps%totpnts, 3))
    p(:) = 0.0
    f(:,1) = 0.0 !sp%g(1)
    f(:,2) = 0.0 !sp%g(2)
    f(:,3) = 0.0 !sp%g(3)
    call set_v(sp)
    call set_v0(ps%pnts, phi, ps%totpnts, v)
    !assigning the boundary conditions
    call set_bc(ps%pnts, ps%totpnts, bcs)
    !assigning the solver control parameters
    call set_krnl(ps%dx, sp)
    call set_t(sp)
    sp%segma = 0.0
    !running the solver
    call cpu_time(start)
    call slv_ns_1phs(ps, v, p, f, phi, mat, bcs, sp)
    call cpu_time(finish)
    print '("Computation time = ",f8.1," seconds.")',finish-start
end subroutine run_ns_1phs

subroutine slv_ns_1phs(ps, vel1, P1, f, phi, mat, bcs, sp)
!------------------------------------------------------------------------------------------
    use krnl_struct
    use pntst_struct
    use bc_struct
    use slvr_prmtrs_struct
    use mat_struct
    use pntst
    use slvr_cmn
    use bndry
    use trnsprt
    use ns
    use ns_cbs
!------------------------------------------------------------------------------------------
    type(pointset), intent(in)              :: ps
    type(bc), dimension(5), intent(in)      :: bcs
    type(slvr_prmtrs), intent(in)           :: sp
!------------------------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)    :: p1, phi
    real(8), dimension(:,:), intent(inout)  :: vel1, f
    type(materials), intent(inout)          :: mat
!------------------------------------------------------------------------------------------
    type (kernel), dimension(ps%totpnts)    :: krnls
    integer                                 :: c = 1, t
    real(8), dimension(ps%totpnts)          :: vof, q, u
    real(8), dimension(:,:), allocatable    :: vel
    real(8), dimension(:), allocatable      :: v,p, ro
    integer                                 :: i
    real(8)                                 :: nu, k, alpha, ro0 = 0, t0 = 0
    character(len=50)                       :: v_fname = 'v', p_fname = 'p', u_fname = 'u'
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Initializing field variable distribution..."
    allocate(vel(ps%totpnts,3))
    allocate(v(ps%totpnts))
    allocate(P(ps%totpnts))
    allocate(ro(ps%totpnts))
    vof(:) = 1.0
    call set_var(vel(:,1), bcs(1))
    call set_var(vel(:,2), bcs(2))
    call set_var(vel(:,3), bcs(3))
    v = sqrt(vel(:,1) * vel(:,1) + vel(:,2) * vel(:,2) + vel(:,3) * vel(:,3))
!------------------------------------------------------------------------------------------
    if (mat%d(1) > 0) then
        !rayleigh-benard convection
        u(:) = 0.0
        q(:) = 0.0
        u(bcs(5)%pnts) = bcs(5)%vals
        do i = 1, ps%totpnts
            if (phi(i) > 0.5) then
                u(i) = sp%Thot + sp%Thot / 5.0
            end if
        end do
        nu = sqrt(sp%pr / sp%ra) * sp%dt / (ps%dx * ps%dx)
        k = sqrt(1.0 / (sp%pr * sp%ra)) * sp%dt / (ps%dx * ps%dx)
        mat%d = [k, k]
        mat%nu = [nu, nu]
        alpha = 0.01
        ro0 = mat%ro(1)
        T0 = 0.5*(sp%Thot - sp%Tcold)
    end if
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Constructing krnls..."
    call get_krnls(ps, sp, krnls)
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Computing transient solution..."
    if (sp%vtk .eqv. .true.) then
        call prnt_vtk(V, ps, v_fname, 0)
        call prnt_vtk(p, ps, p_fname, 0)
        if (mat%d(1) > 0.0) call prnt_vtk(u, ps, u_fname, 0)
    else
        call prnt_pltctrl(ps%pnts, sp%nt, sp%prnt_frq)
        call prnt_elems_txt(ps%surfs, ps%totsurfs)
        call prnt_txt(v, ps%pnts, ps%totpnts, v_fname, 0)
        call prnt_txt(p, ps%pnts, ps%totpnts, p_fname, 0)
        if (mat%d(1) > 0.0) call prnt_txt(u, ps%pnts, ps%totpnts, u_fname, 0)  
    end if
    do t = 1,sp%nt
        if (mat%d(1) > 0.0) then
            !solve the energy equation
            do i = 1, ps%totpnts
                ro(i) = ro0 - ro0 * alpha * (U(i) - T0)
            end do

            call trnsprt_no_upwind(u, krnls, ps%totpnts, mat%d, vof, vel, q, sp%dt)
            u(bcs(5)%pnts) = bcs(5)%vals

            do i = 1, ps%totpnts
                if (phi(i) > 0.5) then
                    u(i) = sp%thot + sp%thot / 5.0
                end if
            end do

            do i = 1, ps%totpnts
                f(i,1) = sp%g(1) * sp%t * (U(i) - T0) / (sp%thot - sp%tcold)
                f(i,2) = sp%g(2) * sp%t * (U(i) - T0) / (sp%thot - sp%tcold)
                f(i,3) = sp%g(3) * sp%t * (U(i) - T0) / (sp%thot - sp%tcold)
            end do
        end if

        !solve for velcity and pressure
        call slv_ns_cbs(vel, v, p, ro, ps%pnts, ps%totpnts, krnls, vel1, p1, &
        F, vof, mat%nu, mat%ro, mat%total, sp%segma, bcs, &
         sp%dt, sp%p, sp%lgr)

        call set_var(vel(:,1), bcs(1))
        call set_var(vel(:,2), bcs(2))
        call set_var(vel(:,3), bcs(3))
        
        !zero velocity at implicit boundaries
        do i = 1, ps%totpnts
            if (phi(i) < 0.0) vel(i,:) = 0
        end do
        
        vel1 = vel
        p1 = p

        !print solutions
        if (t/sp%prnt_frq == c) then
            if (sp%vtk .eqv. .true.) then 
                call prnt_vtk(V, ps, v_fname, t)
                call prnt_vtk(p, ps, p_fname, t)
                if (mat%d(1) > 0.0) call prnt_vtk(u, ps, u_fname, t)
            else
                call prnt_txt(v, ps%pnts, ps%totpnts, v_fname, t)
                call prnt_txt(p, ps%pnts, ps%totpnts, p_fname, t)
                if (mat%d(1) > 0.0) call prnt_txt(u, ps%pnts, ps%totpnts, u_fname, t)  
            end if
            c = c + 1
        end if
    end do
end subroutine slv_ns_1phs

end module ns_1phs_slvrs
