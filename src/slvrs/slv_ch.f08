module ch_slvrs
implicit none
contains

subroutine run_ch()
!------------------------------------------------------------------------------------------
!   This subroutine solves the cahn-hilliard phase separation equations using       
!   explicit time stepping schemes and local strong-form meshfree methods                                                            
!------------------------------------------------------------------------------------------
    use pntst_struct
    use slvr_prmtrs_struct
    use prmtrs
    use pntst
    use omp_lib
!------------------------------------------------------------------------------------------
    real(8)                                 :: start, finish !timer parameters
    real(8), dimension(:), allocatable      :: phi  !phase-field order parameter
    real(8), dimension(:), allocatable      :: vn   !normal velocity
    real(8), dimension(:,:), allocatable    :: v    !advective velocity in xyz
    type(pointset)                          :: ps   !pointset
    type(slvr_prmtrs)                       :: sp   !solver parameters data structure
!------------------------------------------------------------------------------------------                                                                                                    
    write(*,'(a)') "--------------------------------------------------------"
    write(*,'(a)') "Running the Cahn-Hilliard solver..."
    write(*,'(a)') "--------------------------------------------------------"
    !construct pointset
    call set_geom(ps)
    !construct implicit geometry
    call set_phi(ps%dx, sp) 
    call set_phi0(ps%pnts, ps%totpnts, ps%dx, phi)
    !assign advective velocity
    allocate(vn(ps%totpnts))
    vn(:) = 0.0
    call set_v0(ps%pnts, phi, ps%totpnts, v)
    !assign solver parameters
    call set_krnl(ps%dx, sp)
    call set_t(sp)
    !run solver
    start = omp_get_wtime()
    call slv_ch(ps, sp, phi, v, vn)
    finish = omp_get_wtime()
    print '("Computation time = ",f8.1," seconds.")',finish-start
end subroutine run_ch

subroutine slv_ch(ps, sp, phi, vel, vn)
!------------------------------------------------------------------------------------------  
    use pntst_struct
    use krnl_struct
    use slvr_prmtrs_struct
    use slvr_cmn
    use intrf
    use pntst
    use bndry
!------------------------------------------------------------------------------------------
    type(pointset), intent(in)              :: ps
    type(slvr_prmtrs), intent(in)           :: sp 
    real(8),dimension(:),intent(in)         :: vn
    real(8), dimension(:,:), intent(in)     :: vel
!------------------------------------------------------------------------------------------  
    real(8), dimension(:), intent(inout)    :: phi
!------------------------------------------------------------------------------------------  
    type (kernel), dimension(ps%totpnts)    :: krnls
    integer                                 :: c = 1, t, i
    character(len=50)                       :: fname = 'phi'
!------------------------------------------------------------------------------------------ 
    write(*,'(a)') "Constructing krnls..."
    call get_krnls(ps, sp, krnls)
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Computing transient solution..."
    if (sp%vtk .eqv. .true.) then
        call prnt_vtk(phi, ps, fname, 0)
    else
        call prnt_pltctrl(ps%pnts, sp%nt, sp%prnt_frq)
        call prnt_elems_txt(ps%surfs, ps%totsurfs)
        call prnt_txt(phi, ps%pnts, ps%totpnts, fname, 0)
    end if
    do t = 1,sp%nt    
        select case (sp%krnl)
            case (1:5)
                call ch(phi, vel, vn, krnls, ps%totpnts, sp%w, sp%m, sp%dt, sp%itrs)
            case (6)
                call ch_gfd(phi, vel, vn, krnls, ps%totpnts, sp%w, sp%m, sp%dt, sp%itrs)
        end select
        do i = 1, ps%totpnts
            if (phi(i) > 1) phi(i) = 1
            if (phi(i) < 0) phi(i) = 0
        end do
        if (t/sp%prnt_frq == c) then
            if (sp%vtk .eqv. .true.) then
                call prnt_vtk(phi, ps, fname, t)
            else
                call prnt_txt(phi, ps%pnts, ps%totpnts, fname, t)
            end if
            c = c + 1
        end if
    end do
end subroutine slv_ch

end module ch_slvrs
