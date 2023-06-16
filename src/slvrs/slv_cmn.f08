module slvr_cmn
implicit none
contains

subroutine get_krnls(ps, sp, krnls)
    use pntst_struct
    use slvr_prmtrs_struct
    use krnl_struct
    use krnl_cmn
    use krnl_sph
    use krnl_wls
    use krnl_rbf
    use krnl_mls
    use krnl_krg
    use krnl_gfd

    type(pointset), intent(in)              ::  ps
    type(slvr_prmtrs), intent(in)           ::  sp
    type (kernel), dimension(ps%totpnts), intent(out)    ::  krnls

    integer, dimension(:), allocatable      :: pntbins
    integer                                 :: nx, nxy, totbins
    type(kernel), dimension(:), allocatable :: bins
    real(8)                                 ::  start, finish

    nx = 20
    call cpu_time(start)

    call gen_bins(ps%pnts, ps%totpnts, ps%dim, ps%dx, nx, nxy, totbins, bins, pntbins)
!------------------------------------------------------------------------------------------
    !Collecting neighbouring points
    !call get_nbrs_bf(sp%h, ps%pnts, ps%totpnts, ps%pnts, ps%totpnts, krnls)
    call get_nbrs_bg(ps%pnts, pntbins, ps%totpnts, bins, nx, nxy, totbins, ps%dim, sp%h, krnls)
    
    !Computing interpolation functions and their first derivatives
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
        case (6)
            call get_gfd_krnls(ps%pnts, ps%totpnts, ps%pnts, sp%gfd, sp%h, krnls)
    end select

    !Approximating second order derivatives of interpolants
    if (sp%krnl < 6) then
        call get_intrps_o2(ps%pnts, ps%totpnts, krnls)
    end if
    call cpu_time(finish)
    print '("intrpolants construction time = ",f10.1," seconds")',finish-start
end subroutine get_krnls

end module slvr_cmn
