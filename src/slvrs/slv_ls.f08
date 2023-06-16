module ls_slvrs
implicit none
contains

subroutine run_ls()
!------------------------------------------------------------------------------------------
! This subroutine solves the level-set, allen-cahn, cahn-hilliard equations
! for diffuse-interface propagation under a pre-determined velocity field. 
! Solution is obtained using local strong-form meshfree methods and iterative 
! explicit time stepping schemes and 
!------------------------------------------------------------------------------------------
    use pntst_struct
    use slvr_prmtrs_struct
    use prmtrs
    use pntst
    use slvr_cmn
    use omp_lib
!------------------------------------------------------------------------------------------
    real(8),dimension(:),allocatable        ::  phi      !level set function distribution
    real(8), dimension(:,:), allocatable    ::  v        !advective velocity in xyz
    real(8), dimension(:), allocatable      ::  vn       !normal advective velocity
    type(pointset)                          ::  ps       !pointset data structure
    type(slvr_prmtrs)                       ::  sp       !solver parameters data structure
    real(8)                                 ::  start, finish !timer parameters
!------------------------------------------------------------------------------------------
    write(*,'(a)') "--------------------------------------------------------"
    write(*,'(a)') "Running level-set solver..."
    write(*,'(a)') "--------------------------------------------------------"
    !construct pointset
    call set_geom(ps)  
    !construct implicit geometry
    call set_phi(ps%dx, sp)     
    call set_phi0(ps%pnts, ps%totpnts, ps%dx, phi)
    !assign advective velocity
    call set_v(sp)
    call set_v0(ps%pnts, phi, ps%totpnts, v)
    call set_vn0(ps%pnts, ps%totpnts, vn)
    !assign solver parameters
    call set_krnl(ps%dx, sp)    !interpolation kernel parameters
    call set_t(sp)              !time parameters
    !run solver
    !call cpu_time(start)
    start = omp_get_wtime()
    call slv_ls(ps, sp, phi, v, vn)
    finish = omp_get_wtime()
    !call cpu_time(finish)
    print '("Computation time = ",f8.1," seconds.")',finish-start
end subroutine run_ls

subroutine slv_ls(ps, sp, phi, vel, vn)
!------------------------------------------------------------------------------------------
    use krnl_struct
    use pntst_struct
    use slvr_prmtrs_struct
    use slvr_cmn
    use intrf
    use pntst
    use bndry
!------------------------------------------------------------------------------------------
    type(pointset), intent(in)              ::  ps
    type(slvr_prmtrs), intent(in)           ::  sp
    real(8), dimension(:,:), intent(in)     ::  vel
!------------------------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)    ::  vn
    real(8), dimension(:), intent(inout)    ::  phi
!------------------------------------------------------------------------------------------
    type (kernel), dimension(ps%totpnts)    ::  krnls
    integer                                 ::  i, c = 1, t
    real(8), dimension(:), allocatable ::  curv
    character(len=50)                       ::  fname = 'phi'
    real(8)                                 ::  w, m, mobility = 1, vol0, vol
!------------------------------------------------------------------------------------------
    vol0 = sum(phi)
    w = ps%dx
    m = mobility*w*w
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Constructing krnls..."
    call get_krnls(ps, sp, krnls)
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Computing transient solution..."
    if (sp%vtk .eqv. .true.) then
        select case (ps%dim)
        case (3)
            call prnt_unstruct_vtk(phi, ps, fname, 0)
        case default
            call prnt_vtk(phi, ps, fname, 0)
        end select 
    else
        call prnt_pltctrl(ps%pnts, sp%nt, sp%prnt_frq)
        call prnt_elems_txt(ps%surfs, ps%totsurfs)
        call prnt_txt(phi, ps%pnts, ps%totpnts, fname, 0)
    end if

    do t = 1,sp%nt    
        if (abs(sp%segma) > 0.0) then
            select case(sp%krnl)
                case (1:5)
                    call get_curv(phi, ps%pnts, krnls, ps%totpnts, curv)
                case (6)
                    call get_curv_gfd(phi, krnls, ps%totpnts, curv)
            end select
            vn = sp%segma * curv
        end if
        select case(sp%krnl)
            case (1:5)
                select case (sp%intrf_slvr)
                    case (1)
                        call vof(phi, vel, vn, krnls, ps%totpnts, sp%dt, sp%itrs)
                    case (2)
                        call vof2(phi, vel, vn, krnls, ps%totpnts, sp%dt, sp%itrs)
                    case (3)
                        call ac(phi, vel, vn, krnls, ps%totpnts, w, m, sp%dt, sp%itrs)
                    case (4)
                        call ac2(phi, vel, vn, krnls, ps%totpnts, m, sp%dt, sp%itrs)
                    case (5)
                        call ch(phi, vel, vn, krnls, ps%totpnts, w, m, sp%dt, sp%itrs)
                    case (6)    
                        call ch2(phi, vel, vn, krnls, ps%totpnts, w, m, sp%dt, sp%itrs)
                end select
            case (6)
                call vof_gfd(phi, vel, vn, krnls, ps%totpnts, sp%dt, sp%itrs)
        end select
        if (sp%sharp_intrf .eqv. .true.) then
            if(sp%intrf_slvr <= 3 .or. sp%intrf_slvr == 5) then
                call shrpn_intrf(phi, krnls, ps%totpnts, .true., sp%si_beta, sp%si_eps, ps%dx, sp%si_dt, sp%si_nt)
            else
                call shrpn_intrf(phi, krnls, ps%totpnts, .false., sp%si_beta, sp%si_eps, ps%dx, sp%si_dt, sp%si_nt)
            end if
        end if
        do i = 1, ps%totpnts
            if (phi(i) > 1) phi(i) = 1
            if (phi(i) < 0) phi(i) = 0
        end do
        if (t/sp%prnt_frq == c) then
            vol = 100 * (sum(phi) - vol0) / vol0
            write(*,'(a,f10.3)') "Volume loss % = ", vol
            if (sp%vtk .eqv. .true.) then
                select case (ps%dim)
                case (3)
                    call prnt_unstruct_vtk(phi, ps, fname, t)
                case default
                    call prnt_vtk(phi, ps, fname, t)
                end select
            else
                call prnt_txt(phi, ps%pnts, ps%totpnts, fname, t)
            end if
            c = c + 1
        end if
    end do
end subroutine slv_ls

end module ls_slvrs
