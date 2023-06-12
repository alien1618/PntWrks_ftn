module intrp_slvrs
implicit none
contains

subroutine run_sctr_intrp()
!------------------------------------------------------------------------------------------
! This subroutine performs scattered point interpolation to approximate the peaks 
! function and its first and second order derivatives defined on randomly scattered
! points.                                                                                                                               
!------------------------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_mls
    use krnl_wls
    use krnl_krg
    use pntst_struct
    use slvr_prmtrs_struct
    use prmtrs
    use pntst
!------------------------------------------------------------------------------------------
    integer                                  :: i, s, t, nbr, totsctr
    real(8)                                  :: start, finish !timer start and finish
    real(8)                                  :: df  
    real(8), dimension(:,:), allocatable     :: f, fa       !exact and approximation function
    real(8), dimension(:,:), allocatable     :: sctrpnts    !scatter points
    type(pointset)                           :: ps          !point set data structure
    type(slvr_prmtrs)                        :: sp          !solver parameters data structure
    type (kernel), dimension(:), allocatable :: krnls       !kernels data structure
    character(len=50)                        :: f_str, fx_str, fy_str, fxx_str, fyy_str
    character(len=50)                        :: fa_str, fxa_str, fya_str, fxxa_str, fyya_str
    character(len=50)                        :: nabla2_f_str, nabla2_fa_str
!------------------------------------------------------------------------------------------
    f_str = 'f'; fx_str = 'fx'; fy_str = 'fy'; fxx_str = 'fxx'; fyy_str = 'fyy';
    fa_str = 'fa'; fxa_str = 'fxa'; fya_str = 'fya'; fxxa_str = 'fxxa'; fyya_str = 'fyya';
    nabla2_fa_str = 'nabla2_fa'; nabla2_f_str = 'nabla2_f'
!------------------------------------------------------------------------------------------
    call set_geom(ps)
    call set_sctr(sctrpnts, totsctr) 
    call set_krnl(ps%dx, sp)
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Calculating exact solution of the peaks function and its derivatives at scattered points..."
    allocate(fa(ps%totpnts, 6))
    allocate(f(totsctr,6))
    call get_pks(sctrpnts, totsctr, f)
    !call get_f2(sctrpnts, totsctr, f)
    call prnt_pks_vtk(f(:,1), sctrpnts, totsctr, f_str, 0)
    call prnt_pks_vtk(f(:,2), sctrpnts, totsctr, fx_str, 0)
    call prnt_pks_vtk(f(:,3), sctrpnts, totsctr, fy_str, 0)
    call prnt_pks_vtk(f(:,4), sctrpnts, totsctr, fxx_str, 0)
    call prnt_pks_vtk(f(:,5), sctrpnts, totsctr, fyy_str, 0)
    call prnt_pks_vtk(f(:,6), sctrpnts, totsctr, nabla2_f_str, 0)
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Constructing interpolation kernels..."
    call cpu_time(start)
    allocate(krnls(ps%totpnts))
    call get_nbrs_bf(sp%h, ps%pnts, ps%totpnts, sctrpnts, totsctr, krnls)
    do t = 1, 20
        select case(t)
            case (1)
                sp%rbf = 1 !mq
                sp%rbf_alpha = 0.25
                call get_rbf_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, sp%rbf_polyex, krnls)
            case (2) 
                sp%rbf = 2 !imq
                sp%rbf_alpha = 0.25
                call get_rbf_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, sp%rbf_polyex, krnls)
            case (3)
                sp%rbf = 3 !ge
                sp%rbf_alpha = 10
                call get_rbf_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, sp%rbf_polyex, krnls)
            case (4)
                sp%mls = 1 !s3
                call get_mls_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%mls, sp%order, sp%h, krnls)
            case (5)
                sp%mls = 2 !s4
                call get_mls_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%mls, sp%order, sp%h, krnls)
            case (6)
                sp%mls = 3 !s5
                call get_mls_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%mls, sp%order, sp%h, krnls)
            case (7)
                sp%mls = 4 !reg_spln
                call get_mls_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%mls, sp%order, sp%h, krnls)
            case (8)
                sp%wls = 0 !ls
                call get_wls_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%wls, sp%order, sp%h, krnls)
            case (9)
                sp%wls = 1 !ge
                call get_wls_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%wls, sp%order, sp%h, krnls)
            case (10)
                sp%wls = 2 !s3
                call get_wls_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%wls, sp%order, sp%h, krnls)
            case (11)
                sp%wls = 3 !s4
                call get_wls_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%wls, sp%order, sp%h, krnls)
            case (12)
                sp%wls = 4 !s5
                call get_wls_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%wls, sp%order, sp%h, krnls)
            case (13)
                sp%sph = 1 !ge
                call get_sph_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%sph, ps%dim, sp%h,  krnls)
            case (14)
                sp%sph = 2 !s3
                call get_sph_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%sph, ps%dim, sp%h,  krnls)               
            case (15)
                sp%sph = 3 !w4
                call get_sph_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%sph, ps%dim, sp%h,  krnls)
            case (16)
                sp%sph = 4 !w5
                call get_sph_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%sph, ps%dim, sp%h,  krnls)
            case (17)
                sp%sph = 5 !inv_dis
                call get_sph_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%sph, ps%dim, sp%h,  krnls)
            case (18)
                sp%rbf = 1 !mq
                sp%rbf_alpha = 0.25
                call get_krg_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, krnls)
            case (19) 
                sp%rbf = 2 !imq
                sp%rbf_alpha = 0.25
                call get_krg_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, krnls)
            case (20)
                sp%rbf = 3 !ge
                sp%rbf_alpha = 10
                call get_krg_krnls(ps%pnts, ps%totpnts, sctrpnts, sp%rbf, sp%order, sp%h, sp%rbf_alpha, krnls)
            case default
        end select
        call get_intrps_o2(ps%pnts, ps%totpnts, krnls)
        call cpu_time(finish)
        print '("intrpolants construction time = ",f10.1," seconds")',finish-start
!------------------------------------------------------------------------------------------
        write(*,'(a)') "Approximating solution of peaks function derivatives..."
        fa(:,:) = 0
        do i = 1, ps%totpnts
            do s = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(s)
                fa(i,1) = fa(i,1) + krnls(i)%n(s) * f(nbr,1)
                !df = f(nbr,1)-f(i,1)
                df = f(nbr,1)
                fa(i,2) = fa(i,2) + krnls(i)%nx(s) * df
                fa(i,3) = fa(i,3) + krnls(i)%ny(s) * df
                !fa(i,4) = fa(i,4) + krnls(i)%nxx(s) * df
                !fa(i,5) = fa(i,5) + krnls(i)%nyy(s) * df
            end do
            fa(i,6) = fa(i,4) + fa(i,5)
        end do
!------------------------------------------------------------------------------------------      
        write(*,'(a)') "Printing solutions to file..."
        call prnt_pks_vtk(fa(:,1), ps%pnts, ps%totpnts, fa_str, t)
        call prnt_pks_vtk(fa(:,2), ps%pnts, ps%totpnts, fxa_str, t)
        call prnt_pks_vtk(fa(:,3), ps%pnts, ps%totpnts, fya_str, t)
        !call prnt_pks_vtk(fa(:,4), ps%pnts, ps%totpnts, fxxa_str, t)
        !call prnt_pks_vtk(fa(:,5), ps%pnts, ps%totpnts, fyya_str, t)
        call prnt_pks_vtk(fa(:,6), ps%pnts, ps%totpnts, nabla2_fa_str, t)   
!------------------------------------------------------------------------------------------   
        write(*,'(a)') "Deallocating memory..."    
        do i = 1, ps%totpnts
            deallocate(krnls(i)%n)
            deallocate(krnls(i)%nx)
            deallocate(krnls(i)%ny)
            deallocate(krnls(i)%nz)
            !deallocate(krnls(i)%nxx)
            !deallocate(krnls(i)%nyy)
            !deallocate(krnls(i)%nzz)
            !if (t < 12) then
            !    deallocate(krnls(i)%nxy)
            !    deallocate(krnls(i)%nxz)
            !    deallocate(krnls(i)%nyz)
            !end if
            if (sp%order == 1) then
                deallocate(krnls(i)%nabla2)
            end if
        end do
    end do
end subroutine run_sctr_intrp

subroutine get_pks(pnts, totpnts, f)
!-------------------------------------------------------------------------------------------------
! Subroutine computes the exact solution of the peaks function and its derivatives
! f(:,1) = f, f(:,2) = dfdx, f(:,3)  = dfdy, f(:,4) = dfdxx, f(:,5) = dfdyy, f(:,6) = dfdxx+dfdyy
!-------------------------------------------------------------------------------------------------
    integer, intent(in)              :: totpnts
    real(8), dimension(:,:), intent(in) :: pnts
!-------------------------------------------------------------------------------------------------
    real(8), dimension(totpnts, 6), intent(out) :: f
!-------------------------------------------------------------------------------------------------
    real(8) :: x, y, t
    integer :: i
!-------------------------------------------------------------------------------------------------
    do i = 1, totpnts
        x = pnts(i,1)
        y = pnts(i,2)
        f(i,1) = 3.0*((1.0-x)**2)*exp(-(x**2) - ((y+1)**2))-10.0*(x/5.0 - (x**3) &
                - (y**5))*exp(-(x**2)-(y**2))- 1.0/3.0*exp(-((x+1)**2) - (y**2))
        t = exp(-((y+1)**2)-x*x)
        f(i,2) = -6.0*((1.0-x)**2)*x*t-6.0*(1.0-x)*t +(2.0/3.0)*(x+1.0)*exp(-y*y-((x+1)**2)) &
                +20.0*x*(-(y**5)-(x**3)+x/5.0)*exp(-y*y-x*x)-10.0*(0.2-3.0*x*x)*exp(-y*y-x*x)
        f(i,3) = -6.0*((1.0-x)**2)*(y+1.0)*t+(2.0/3.0)*y*exp(-y*y-((x+1.0)**2)) &
                +20.0*y*(-(y**5)-(x**3)+x/5.0)*exp(-y*y-x*x)+50.0*(y**4)*exp(-y*y-x*x)
        f(i,4) = 12.0*((1.0-x)**2)*x*x*t+24.0*(1.0-x)*x*t-6*((1-x)**2)*t &
                +6.0*t-(4.0/3.0)*((x+1)**2)*exp(-y*y-((x+1)**2)) &
                +(2.0/3.0)*exp(-y*y-((x+1)**2))-40.0*x*x*(-(y**5)-(x**3) &
                +(x/5))*exp(-y*y-x*x)+20*(-(y**5)-(x**3)+(x/5))*exp(-y*y-x*x) &
                +40.0*x*(0.2-3*x*x)*exp(-y*y-x*x)+60*x*exp(-y*y-x*x)
        f(i,5) = 12.0*((1-x)**2)*((y+1)**2)*t-6.0*((1.0-x)**2)*t-(4.0/3.0)*y*y*exp(-y*y-((x+1.0)**2)) &
                +(2.0/3.0)*exp(-y*y-((x+1.0)**2))-200*(y**5)*exp(-y*y-x*x) & 
                -40*y*y*(-(y**5)-(x**3)+(x/5))*exp(-y*y-x*x) &
                +20*(-(y**5)-(x**3)+(x/5))*exp(-y*y-x*x)+200*(y**3)*exp(-y*y-x*x)
        f(i,6) = f(i,4) + f(i,5)
    end do
end subroutine get_pks

subroutine get_f1(pnts, totpnts, f)
!-------------------------------------------------------------------------------------------------
! Subroutine computes the exact solution of a function and its derivatives
! f(:,1) = f, f(:,2) = dfdx, f(:,3)  = dfdy, f(:,4) = dfdxx, f(:,5) = dfdyy, f(:,6) = dfdxx+dfdyy
!-------------------------------------------------------------------------------------------------
    integer, intent(in)                 :: totpnts
    real(8), dimension(:,:), intent(in) :: pnts
!-------------------------------------------------------------------------------------------------
    real(8), dimension(totpnts, 6), intent(out) :: f
!-------------------------------------------------------------------------------------------------
    real(8) :: x, y
    integer :: i
!-------------------------------------------------------------------------------------------------    
    do i = 1, totpnts
        x = pnts(i, 1)
        y = pnts(i, 2)
        f(i, 1) = x * (1.0 - x) * y * (1.0 - y)
        f(i, 2) = y * (1.0 - x) * (1.0 - y) - x * y * (1.0 - y)
        f(i, 3) = x * (1.0 - x) * (1.0 - y) - x * y * (1.0 - x)
        f(i, 4) = - y * (1.0 - y) - y * (1.0 - y)
        f(i, 5) = - x * (1.0 - x) - x * (1.0 - x)
        f(i, 6) = f(i, 4) + f(i, 5)
    end do
end subroutine get_f1

subroutine get_f2(pnts, totpnts, f)
!-------------------------------------------------------------------------------------------------
! Subroutine computes the exact solution of a function and its derivatives
! f(:,1) = f, f(:,2) = dfdx, f(:,3)  = dfdy, f(:,4) = dfdxx, f(:,5) = dfdyy, f(:,6) = dfdxx+dfdyy
!-------------------------------------------------------------------------------------------------
    integer, intent(in)                 :: totpnts
    real(8), dimension(:,:), intent(in) :: pnts
!-------------------------------------------------------------------------------------------------
    real(8), dimension(totpnts, 6), intent(out) :: f
!-------------------------------------------------------------------------------------------------
    real(8) :: x, y
    integer :: i    
!-------------------------------------------------------------------------------------------------
    do i = 1, totpnts
        x = pnts(i, 1)
        y = pnts(i, 2)
        f(i, 1) = x * x * (1.0 - x) * y * ((1.0 - y)**2)
        f(i, 2) = 2 * x * (1.0 - x) * y * ((1.0 - y)**2) - x * x * y * ((1.0 - y)**2)
        f(i, 3) = x * x * (1.0 - x) * ((1.0 - y)**2) - x * x * (1.0 - x) * y * 2 * (1.0 - y)
        f(i, 4) = 2 * (1.0 - x) * y * ((1.0 - y)**2) - 2 * x * y * ((1.0 - y)**2) - 2 * x * y * ((1.0 - y)**2)
        f(i, 5) = - x * x * (1.0 - x) * 2 * (1.0 - y) - x * x * (1.0 - x) * 2 * (1.0 - y) + x * x * (1.0 - x) * y * 2
        f(i, 6) = f(i, 4) + f(i, 5)
    end do
end subroutine get_f2


subroutine run_sctr_reg() 
!------------------------------------------------------------------------------------------
! This subroutine performs solves the parabolic and level set equations on randomly
! scattered points in the domain. 
!------------------------------------------------------------------------------------------
    use pntst_struct
    !use bc_struct
    use slvr_prmtrs_struct
    !use mat_struct
    use pntst
    !use bndry
    !use intrf
    use ns_lgr_slvrs
    !use ls_slvrs
    !use intrf
    use krnl_struct
    use krnl_cmn
    use prmtrs

    implicit none

    real(8), dimension(:,:), allocatable    :: vel        !velocity
    real(8), dimension(:), allocatable      :: v
    type(pointset)                          :: bndr, ps   !point set data structure
    type(slvr_prmtrs)                       :: sp         !solver parameters data structures
    real(8)                                 :: h
    character(len=50)                       :: bndr_fname, v_fname !output file names
    real(8), dimension(:,:), allocatable    :: sctr
    integer                                 :: totsctr, i, t, c
    type (kernel),dimension(:), allocatable :: krnls
    v_fname = 'v'
!------------------------------------------------------------------------------------------
    bndr_fname = 'sim/in/1_pntst/edg.txt'
    call read_pnts(bndr_fname, bndr) 
    call set_sctr(sctr, totsctr) 
    ps%totpnts = totsctr + bndr%totpnts
    allocate(ps%pnts(ps%totpnts, 3))
    ps%pnts(1:totsctr,:) = sctr(1:totsctr,:)
    ps%pnts(totsctr+1:ps%totpnts,:) = bndr%pnts(:,:)
    ps%totsurfs = 0
    ps%dim = 2
    call get_avg_dx(ps%pnts, ps%dx)
    call set_t(sp)
!------------------------------------------------------------------------------------------
    h = 0.05
    c = 1
    allocate(vel(ps%totpnts, 3))
    allocate(v(ps%totpnts))
    v = sqrt(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2)+vel(:,3)*vel(:,3))
    allocate(krnls(ps%totpnts))
    call prnt_vtk(v, ps, v_fname, 0)
    do t = 1, sp%nt
        call get_nbrs_bf(h, ps%pnts, ps%totpnts, ps%pnts, ps%totpnts, krnls)
        !call get_v_shft(0.1d0, 0.0d0, h, ps%pnts, ps%totpnts, krnls, vel)
        call get_v_shft2(1000d0, 0.0d0, h, ps%pnts, ps%totpnts, krnls, vel)
        ps%pnts(1:totsctr,:) = ps%pnts(1:totsctr,:) + sp%dt*vel(1:totsctr,:)
        v = sqrt(vel(:,1)*vel(:,1)+vel(:,2)*vel(:,2)+vel(:,3)*vel(:,3))
        v(totsctr+1:ps%totpnts) = 0.0
        if (t/sp%prnt_frq == c) then
            call prnt_vtk(v, ps, v_fname, t)
            c = c + 1
        end if
        do i = 1, ps%totpnts
            krnls(i)%totnbrs = 0
            deallocate(krnls(i)%nbrs)
        end do
    end do
end subroutine run_sctr_reg

end module intrp_slvrs
