module ns_lgr_slvrs
implicit none
contains

subroutine run_ns_1phs_lgr() 
!------------------------------------------------------------------------------------------
! This subroutine solves the lagrangian one-phase navier-stokes equations
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
    real(8), dimension(:,:), allocatable    :: v             !velocity
    real(8), dimension(:), allocatable      :: phi           !implicit boundary indicator function
    type(pointset)                          :: ps, walls     !point set data structure
    type(bc), dimension(6)                  :: bcs      !boundary condition data structure for vx, vy, vz, p
    type(slvr_prmtrs)                       :: sp       !solver parameters data structures
    type(materials)                         :: mat      !phase properties data structures
    real(8)                                 :: start, finish !timer parameters
    character(len=50)                       :: walls_fname    !mesh file names
!------------------------------------------------------------------------------------------
    write(*,'(a)') "--------------------------------------------------------"
    write(*,'(a)') "Running navier-stokes lagrangian test..."
    write(*,'(a)') "--------------------------------------------------------"
    !constructing the pointset    
    call set_geom(ps)
    if (ps%dim == 2) then
        walls_fname = 'sim/in/4_bc/walls.txt'  
        call read_pnts(walls_fname, walls) 
    end if
    !constructing the implicit geometry
    allocate(phi(ps%totpnts))
    phi(:) = 1.0
    !assigning the material properties
    mat%total = 1
    call set_mat('nu', mat%nu)
    call set_mat('ro', mat%ro)
    call set_mat('g', sp%g)
    !initializing field variables
    call set_v(sp)
    call set_v0(ps%pnts, phi, ps%totpnts, v)
    !assigning the boundary conditions
    call set_bc(ps%pnts, ps%totpnts, bcs)
    call set_bins(sp)
    !settings the solver parameters
    call set_krnl(ps%dx, sp)
    call set_t(sp)
    sp%segma = 0.0
    !running the solver
    call cpu_time(start)
    call slv_ns_lgr(ps, walls, mat, bcs, sp, v, phi)
    call cpu_time(finish)
    print '("Computation time = ",f8.1," seconds.")',finish-start
end subroutine run_ns_1phs_lgr

subroutine run_ns_2phs_lgr() 
!------------------------------------------------------------------------------------------
! This subroutine solves the lagrangian one-phase navier-stokes equations
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
    real(8), dimension(:,:), allocatable    :: v             !velocity
    real(8), dimension(:), allocatable      :: phi           !implicit boundary indicator function
    type(pointset)                          :: ps, walls     !point set data structure
    type(pointset), dimension(2)            :: ps_arr
    type(bc), dimension(6)                  :: bcs      !boundary condition data structure for vx, vy, vz, p
    type(slvr_prmtrs)                       :: sp       !solver parameters data structures
    type(materials)                         :: mat      !phase properties data structures
    real(8)                                 :: start, finish !timer parameters
    character(len=50)                       :: walls_fname    !mesh file names
    integer :: i, s
!------------------------------------------------------------------------------------------
    write(*,'(a)') "--------------------------------------------------------"
    write(*,'(a)') "Running navier-stokes lagrangian test..."
    write(*,'(a)') "--------------------------------------------------------"
    !constructing the pointset    
    call set_geom_arr(ps_arr)

    ps%totpnts = ps_arr(1)%totpnts + ps_arr(2)%totpnts
    allocate(ps%pnts(ps%totpnts, 3))
    s = 1
    do i = 1, ps_arr(1)%totpnts
        ps%pnts(s, :) = ps_arr(1)%pnts(i, :)
        s = s + 1
    end do
    do i = 1, ps_arr(2)%totpnts
        ps%pnts(s, :) = ps_arr(2)%pnts(i, :)
        s = s + 1
    end do
    ps%dim = ps_arr(1)%dim
    ps%dx = ps_arr(1)%dx

    allocate(phi(ps%totpnts))
    phi(1:ps_arr(1)%totpnts) = 0.0
    phi(ps_arr(1)%totpnts+1:ps_arr(1)%totpnts+ps_arr(2)%totpnts) = 1.0

    if (ps%dim == 2) then
        walls_fname = 'sim/in/4_bc/walls.txt'  
        call read_pnts(walls_fname, walls) 
    end if
    !constructing the implicit geometry

    !assigning the material properties
    mat%total = 1
    call set_mat('nu', mat%nu)
    call set_mat('ro', mat%ro)
    call set_mat('g', sp%g)
    mat%total = 2
    !initializing field variables
    call set_v(sp)
    call set_v0(ps%pnts, phi, ps%totpnts, v)
    !assigning the boundary conditions
    call set_bc(ps%pnts, ps%totpnts, bcs)
    call set_bins(sp)
    !settings the solver parameters
    call set_krnl(ps%dx, sp)
    call set_t(sp)
    sp%segma = 0.0
    !running the solver
    call cpu_time(start)
    call slv_ns_lgr(ps, walls, mat, bcs, sp, v, phi)
    call cpu_time(finish)
    print '("Computation time = ",f8.1," seconds.")',finish-start
end subroutine run_ns_2phs_lgr

subroutine slv_ns_lgr(ps, walls, mat, bcs, sp, vel0, phi0)
    use krnl_struct
    use pntst_struct
    use bc_struct
    use slvr_prmtrs_struct
    use mat_struct
    use pntst
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_mls
    use krnl_wls
    use krnl_krg
    use bndry
    use trnsprt
    use ns
!------------------------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in)     :: vel0
    real(8), dimension(:), intent(in)       :: phi0
    type(pointset), intent(in)              :: ps, walls
    type(bc), dimension(5), intent(in)      :: bcs
    type(slvr_prmtrs), intent(in)           :: sp
!------------------------------------------------------------------------------------------
    type(materials), intent(inout)          :: mat
!------------------------------------------------------------------------------------------
    type (kernel), dimension(:), allocatable :: krnls
    integer                                 :: c = 1, t, totpnts
    real(8), dimension(:), allocatable      :: v, p, p1, phi, ro, ro0
    real(8), dimension(:,:), allocatable    :: vel, vel1, f, pnts, v_shift
    character(len=50)                       :: v_fname = 'v', bg_fname = 'bg'
    character(len=50)                       :: p_fname = 'p', ro_fname = 'ro'
    integer                                 :: i
    integer, dimension(:), allocatable      :: pntbins
    real(8)                                 :: tol
    integer                                 :: nxy, totbins
    type(kernel), dimension(:), allocatable :: bins
    real(8)                                 :: tx, ty, tz
    integer, dimension(3)                   :: n3
    real(8) , dimension(3)                  :: bg_b, bg_l, bg_d
    type(pointset)                          :: bg 
!------------------------------------------------------------------------------------------
!append wall points to fluid points
    totpnts = ps%totpnts
    if (ps%dim == 2) then
        totpnts = ps%totpnts + walls%totpnts
    end if
    allocate(pnts(totpnts, 3))
    pnts(1:ps%totpnts,:) = ps%pnts(1:ps%totpnts,:)
    if (ps%dim == 2) then
        pnts(ps%totpnts+1:totpnts,:) = walls%pnts(:,:)
    end if
!------------------------------------------------------------------------------------------
    !allocate memory
    allocate(vel(totpnts, 3))
    allocate(v(totpnts))
    allocate(p(totpnts))
    allocate(ro(totpnts))
    allocate(f(totpnts,3))
    allocate(v_shift(totpnts, 3))
    allocate(phi(totpnts))
    allocate(ro0(totpnts))
    allocate(vel1(totpnts, 3))
    allocate(p1(totpnts))
    allocate(krnls(totpnts))
!------------------------------------------------------------------------------------------ 
    !initialize variables
    v_shift(:,:) = 0.0d0
    phi(:) = 1.0d0
    phi(1:ps%totpnts) = phi0
    ro = phi*mat%ro(2) + (1.0d0 - phi) * mat%ro(1)
    p(:) = 0.0d0
    vel(:,:) = 0.0d0
    vel(1:ps%totpnts,:) = vel0
    v = sqrt(vel(:,1) * vel(:,1) + vel(:,2) * vel(:,2) + vel(:,3) * vel(:,3))
    f(:,1) = sp%g(1)
    f(:,2) = sp%g(2)
    f(:,3) = sp%g(3)
    vel1 = vel
    p1 = p
!------------------------------------------------------------------------------------------
    ! Calculate bounds of neighbour search bins
    select case (ps%dim)
    case (2)
        call get_bin_bnds(pnts, ps%dx, ps%dim, sp%bg_nx, nxy, totbins, bg_b, bg_l, bg_d)
    case (3)
        !in 3D the bounding box is always assumed from (0,0,0) to (1,1,1)
        bg_b = sp%bg_b
        bg_l = sp%bg_l
        bg_d = [sp%bg_nx, sp%bg_nx, sp%bg_nx]
        nxy = sp%bg_nx**2
        totbins = sp%bg_nx**3
    end select   
    allocate(bins(totbins))
    allocate(pntbins(totpnts)) 
!------------------------------------------------------------------------------------------
    tol = 0.05
    n3 = [2, 2, 2]
    call gen_hexgrd(bg_b, bg_l, n3, bg)
    if (sp%vtk .eqv. .true.) then
        call prnt_pnts_vtk(v, pnts, totpnts, v_fname, 0)
        call prnt_pnts_vtk(p, pnts, totpnts, p_fname, 0)
        call prnt_pnts_vtk(ro, pnts, totpnts, ro_fname, 0)
    else
        call prnt_pltctrl(pnts, sp%nt, sp%prnt_frq)
        call prnt_txt(v, pnts, totpnts, v_fname, 0)
        call prnt_txt(p, pnts, totpnts, p_fname, 0)
        call prnt_txt(ro, pnts, totpnts, ro_fname, 0)
        call prnt_txt(bg%pnts(:,1), bg%pnts, bg%totpnts, bg_fname, 0)
    end if
!------------------------------------------------------------------------------------------
    ! Main loop for transient time evolution
    do t = 1, sp%nt
        ! computing kernel approximation functions
        call get_bins(sp%bg_nx, nxy, totbins, bg_b, bg_d, pnts, totpnts, ps%dim, bins, pntbins)
        call get_nbrs_bg(pnts, pntbins, totpnts, bins, sp%bg_nx, nxy, totbins, ps%dim, sp%h, krnls)
        !call get_nbrs_bf(sp%h, pnts, totpnts, pnts, totpnts, krnls)!source of significant slowdowns

        select case(sp%krnl)
            case (1)
                call get_rbf_krnls(pnts, totpnts, pnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, sp%rbf_polyex, krnls)
            case (2)
                call get_mls_krnls(pnts, totpnts, pnts, sp%mls, sp%order, sp%h, krnls)
            case (3)
                call get_wls_krnls(pnts, totpnts, pnts, sp%wls, sp%order, sp%h, krnls)
            case (4)
                call get_sph_krnls(pnts, totpnts, pnts, sp%sph, ps%dim, sp%h,  krnls)
            case (5)
                call get_krg_krnls(pnts, totpnts, pnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, krnls)
        end select
        call get_intrps_o2(pnts, totpnts, krnls)
        !call get_intrps_o2_v2(pnts, totpnts, krnls)

        call slv_ns(vel, v, p, ro, pnts, totpnts, krnls, vel1, p1, &
        f, phi, mat%nu, mat%ro, mat%total, sp%segma, bcs, &
         sp%dt, sp%av, sp%p, .true.)

        ! computing particle shifting velocity
        if (sp%shift .eqv. .true.) then
            call get_v_shft(sp%shift_prmtr, sp%shift_surf, sp%h, pnts, totpnts, krnls, v_shift)
            !call get_v_shft2(sp%shift_prmtr, sp%shift_surf, sp%h, pnts, totpnts, krnls, v_shift)
            v_shift(ps%totpnts+1:totpnts,:) = 0.0d0       
        else
            v_shift(:,:) = 0.0d8
        end if

        !3d box wall bounce back
        if (ps%dim == 3) then
            do i = 1, totpnts
                tx = pnts(i,1) + sp%dt * (vel(i,1) + v_shift(i,1))
                ty = pnts(i,2) + sp%dt * (vel(i,2) + v_shift(i,2))
                tz = pnts(i,3) + sp%dt * (vel(i,3) + v_shift(i,3))
                
                !in 3D the bounding box is always (0,0,0) to (1,1,1)
                if (tx < bg_b(1) + tol .or. tx > bg_b(1) + bg_l(1) - tol) then
                    vel(i,1) = -vel(i,1)
                    !v_shift(i,1) = -v_shift(i,1)
                end if
                if (ty < bg_b(2) + tol .or. ty > bg_b(2) + bg_l(2) - tol) then
                    vel(i,2) = -vel(i,2)
                    !v_shift(i,2) = -v_shift(i,2)
                end if
                if (tz < bg_b(3) + tol .or. tz > bg_b(3) + bg_l(3) - tol) then
                    vel(i,3) = -vel(i,3)
                    !v_shift(i,3) = -v_shift(i,3)
                end if
            end do
        end if
        !fix wall particles
        if (ps%dim == 2) then
            vel(ps%totpnts+1:totpnts, :) = 0.0d0
        end if
        vel1 = vel
        p1 = p

        ! advecting particles in space
        pnts = pnts + sp%dt * (vel + v_shift)

        !apply periodic boundary conditions
        if (sp%periodic .eqv. .true.) then
            do i = 1, totpnts
                pnts(i,1) = mod(pnts(i,1), bg_l(1))
                if (pnts(i,1) < 0) then
                    pnts(i,1) = bg_l(1) + pnts(i,1)
                endif
                pnts(i,2) = mod(pnts(i,2), bg_l(2))
                if (pnts(i,2) < 0) then
                    pnts(i,2) = bg_l(2) + pnts(i,2)
                endif
            end do
        end if
        v = sqrt(vel(:,1) * vel(:,1) + vel(:,2) * vel(:,2) + vel(:,3) * vel(:,3))
        
        ! print data to file
        if (t/sp%prnt_frq == c) then
            if (sp%vtk .eqv. .true.) then
                call prnt_pnts_vtk(v, pnts, totpnts, v_fname, t)
                call prnt_pnts_vtk(p, pnts, totpnts, p_fname, t)
                call prnt_pnts_vtk(ro, pnts, totpnts, ro_fname, t)
            else
                call prnt_txt(v, pnts, totpnts, v_fname, t)
                call prnt_txt(p, pnts, totpnts, p_fname, t)
                call prnt_txt(ro, pnts, totpnts, ro_fname, t)
            end if
            c = c + 1
        end if

        do i = 1, totbins
            deallocate(bins(i)%nbrs)
            bins(i)%totnbrs = 0
        end do
        do i = 1, totpnts
            deallocate(krnls(i)%nbrs)
            deallocate(krnls(i)%n)
            deallocate(krnls(i)%nx)
            deallocate(krnls(i)%ny)
            deallocate(krnls(i)%nz)
            !deallocate(krnls(i)%nxx)
            !deallocate(krnls(i)%nyy)
            !deallocate(krnls(i)%nzz)
            deallocate(krnls(i)%nabla2)
            !deallocate(krnls(i)%dx)
            !deallocate(krnls(i)%dy)
            !deallocate(krnls(i)%dz)
            !deallocate(krnls(i)%d)
        end do
        
    end do
    write(*,'(a)') "Solver complete..."

end subroutine slv_ns_lgr

pure subroutine get_v_shft(prmtr, srf_mrkr, h, pnts, totpnts, krnls, v_shft)
!------------------------------------------------------------------------------------------ 
    use krnl_struct
!------------------------------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: pnts
    type(kernel), dimension(:), intent(in) :: krnls
    real(8), intent(in) :: prmtr, h, srf_mrkr
    integer, intent(in) :: totpnts
!------------------------------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: v_shft
!------------------------------------------------------------------------------------------ 
    real(8), dimension(:), allocatable :: conc
    real(8) :: dx, dy, dz, d, w, cx, cy
    integer :: i, j, s
!------------------------------------------------------------------------------------------ 
    allocate(v_shft(totpnts, 3))
    allocate(conc(totpnts))
    v_shft(:,:) = 0.0d0

    !calculate particle density
    do i = 1, totpnts
        conc(i) = 0
        do s = 1, krnls(i)%totnbrs
            j = krnls(i)%nbrs(s)
            dx = pnts(i,1) - pnts(j,1)
            dy = pnts(i,2) - pnts(j,2)
            dz = pnts(i,3) - pnts(j,3)
            d = sqrt(dx * dx + dy * dy + dz * dz) + 1e-6
            conc(i) = conc(i) + (((h / d) - 1)**2)
            !conc(i) = conc(i) + ((h / d) - 1)
        end do
    end do
    
    !calculate shifting velocity
    do i = 1, totpnts
        cx = 0.0
        cy = 0.0
        if (conc(i) >= srf_mrkr) then
            do s = 1, krnls(i)%totnbrs
                j = krnls(i)%nbrs(s)
                dx = pnts(i,1) - pnts(j,1)
                dy = pnts(i,2) - pnts(j,2)
                dz = pnts(i,3) - pnts(j,3)
                d = sqrt(dx * dx + dy * dy + dz * dz) + 1e-6
                !double w = ((h/d) - 1)*((h/d) - 1);
                !dC_x = dC_x + (prtcl_conc[j])*((2*h*w)/(d*d))*(points[j].x-px[i]);
                !dC_y = dC_y + (prtcl_conc[j])*((2*h*w)/(d*d))*(points[j].y-py[i]);
                w = ((h/d) - 1)
                cx = cx + conc(j)*(w / (d * d)) * (pnts(j,1) - pnts(i,1))
                cy = cy + conc(j)*(w / (d * d)) * (pnts(j,2) - pnts(i,2))
            end do
        end if
        v_shft(i,1) = - prmtr * h * cx
        v_shft(i,2) = - prmtr * h * cy
    end do
end subroutine get_v_shft

pure subroutine get_v_shft2(prmtr, srf_mrkr, h, pnts, totpnts, krnls, v_shft)
!------------------------------------------------------------------------------------------ 
    use krnl_struct
!------------------------------------------------------------------------------------------ 
    real(8), intent(in) :: prmtr, h, srf_mrkr
    integer, intent(in) :: totpnts
    real(8), dimension(:,:), intent(in) :: pnts
    type(kernel), dimension(:), intent(in) :: krnls
!------------------------------------------------------------------------------------------ 
    real(8), dimension(:,:), allocatable, intent(out) :: v_shft
!------------------------------------------------------------------------------------------ 
    real(8), dimension(:), allocatable :: conc
    real(8) :: dx, dy, dz, d, sumx, sumy, sumz, sumk
    integer :: i, j, nbr
!------------------------------------------------------------------------------------------ 
    allocate(v_shft(totpnts, 3))
    allocate(conc(totpnts))
    v_shft(:,:) = 0.0d0
    
    !calculate particle density
    do i = 1, totpnts
        conc(i) = 0
        do j = 1, krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            dx = pnts(i,1) - pnts(nbr,1)
            dy = pnts(i,2) - pnts(nbr,2)
            dz = pnts(i,3) - pnts(nbr,3)
            d = sqrt(dx * dx + dy * dy + dz * dz) + 1e-6
            conc(i) = conc(i) + (((h / d) - 1)**2)
            !conc(i) = conc(i) + ((h / d) - 1)
        end do
    end do
    
    !calculate shifting velocity
    do i = 1, totpnts
        sumx = 0.0
        sumy = 0.0
        sumz = 0.0
        sumk = 0.0
        if (conc(i) >= srf_mrkr) then
            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                dx = pnts(nbr,1) - pnts(i,1)
                dy = pnts(nbr,2) - pnts(i,2)
                dz = pnts(nbr,3) - pnts(i,3)
                d = sqrt(dx * dx + dy * dy + dz * dz)
                !sumx = sumx + prmtr*((h/d)-1)*dx
                !sumy = sumy + prmtr*((h/d)-1)*dy
                !sumz = sumz + prmtr*((h/d)-1)*dz
                !sumk = sumk + prmtr

                sumx = sumx + ((h / d) - 1) * dx
                sumy = sumy + ((h / d) - 1) * dy
                sumz = sumz + ((h / d) - 1) * dz
            end do
            !v_shft(i,1) = -(sumx/sumk)/dt
            !v_shft(i,2) = -(sumy/sumk)/dt
            !v_shft(i,3) = -(sumz/sumk)/dt

            v_shft(i,1) = -prmtr*sumx
            v_shft(i,2) = -prmtr*sumy
            v_shft(i,3) = -prmtr*sumz
        end if
    end do
end subroutine get_v_shft2

end module ns_lgr_slvrs
