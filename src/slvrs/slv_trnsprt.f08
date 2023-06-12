
module trnsprt_slvrs
implicit none
contains

subroutine  run_trnsprt()
!------------------------------------------------------------------------------------------
! This subroutine solves special forms of the general tranport equation using       
! local strong-form meshfree methods and iterative explicit time stepping schemes                                                                                                        
!------------------------------------------------------------------------------------------
    use pntst_struct
    use bc_struct
    use slvr_prmtrs_struct
    use pntst
    use bndry
    use intrf
    use prmtrs
    use omp_lib
!------------------------------------------------------------------------------------------
    real(8), dimension(:), allocatable      :: k         !diffusivity for 2 phases
    real(8), dimension(:), allocatable      :: q         !generation trm
    real(8), dimension(:), allocatable      :: u     !field variable distribution
    real(8), dimension(:,:), allocatable    :: v         !advective velocity in xyz
    real(8), dimension(:), allocatable      :: phi       !volume of fluid distribution (0,1)
    type(pointset)                          :: ps        !point set data structure
    type(bc), dimension(6)                  :: bcs       !boundary conditions data structure
    type(slvr_prmtrs)                       :: sp        !solver parameters data structure
    real(8)                                 :: start, finish !timer parameters
!------------------------------------------------------------------------------------------
    write(*,'(a)') "--------------------------------------------------------"
    write(*,'(a)') "Running transient transport solver..."
    write(*,'(a)') "--------------------------------------------------------"
    !constructing the pointset
    call set_geom(ps)
    !constructing the implicit geometry
    call set_phi0(ps%pnts, ps%totpnts, ps%dx, phi)
    !assigning the material properties
    call set_mat('k', k)
    !initializing field variables
    call set_v0(ps%pnts, phi, ps%totpnts, v)
    call set_u0(phi, ps%totpnts, u)
    allocate(q(ps%totpnts))
    q(:) = 0.0
    !assigning the boundary conditions
    call set_bc(ps%pnts, ps%totpnts, bcs)
    !settings the solver parameters
    call set_krnl(ps%dx, sp)
    call set_t(sp)
    !running the solver
    !call cpu_time(start)
    start = omp_get_wtime()
    call slv_trnsprt(ps, bcs(1), sp, k, u, phi, v, q)
    !call cpu_time(finish)
    finish = omp_get_wtime()
    print '("Computation time = ",f8.1," seconds.")',finish-start
end subroutine run_trnsprt

subroutine  slv_trnsprt(ps, bcs, sp, d, u, vof, vel, q)
!------------------------------------------------------------------------------------------
    use pntst_struct
    use krnl_struct
    use bc_struct
    use slvr_prmtrs_struct
    use pntst
    use krnl_cmn
    use krnl_sph
    use krnl_rbf
    use krnl_mls
    use krnl_wls
    use krnl_krg
    use krnl_gfd
    use bndry
    use trnsprt
!------------------------------------------------------------------------------------------
    real(8), dimension(2), intent(in)          :: d
    real(8), dimension(:), intent(in)          :: q, vof
    real(8), dimension(:,:), intent(in)        :: vel
    type(slvr_prmtrs), intent(in)              :: sp
    type(pointset), intent(in)                 :: ps
    type(bc), intent(in)                       :: bcs
!------------------------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)       :: u
!------------------------------------------------------------------------------------------
    type(kernel), dimension(ps%totpnts)        :: krnls
    character(len=50)                          :: fname = 'u'
    integer                                    :: c = 1, t
    real(8), dimension(ps%totpnts)             :: rhs
    real(8)                                    :: start, finish
!------------------------------------------------------------------------------------------
    rhs(:) = 0.0d0
!------------------------------------------------------------------------------------------
    ! Calculate bounds of neighbour search bins
    !nx = 30
    !call gen_bins(ps%pnts, ps%totpnts, ps%dim, ps%dx, nx, nxy, totbins, bins, pntbins)
    !allocate(krnls(ps%totpnts))
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Constructing krnls..."
    call cpu_time(start)
    call get_nbrs_bf(sp%h, ps%pnts, ps%totpnts, ps%pnts, ps%totpnts, krnls)
    
    !call get_nbrs_sp(sp%h, vof, ps%pnts, ps%totpnts, krnls) !be careful, it causes errors for advection-diffusion case.
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
            call get_gfd_krnls(ps%pnts, ps%totpnts, ps%pnts, sp%gfd, sp%h,  krnls)
    end select
    if (sp%krnl < 6) then
        call get_intrps_o2(ps%pnts, ps%totpnts, krnls)
    end if
    call cpu_time(finish)
    print '("Intrpolants construction time = ",f10.1," seconds")',finish-start
!---------------------------------------------------------------------------------------
    write(*,'(a)') "Computing transient solution..."
    call set_var(u, bcs)
    if (sp%vtk .eqv. .true.) then
        call prnt_vtk(u, ps, fname, 0)
    else
        call prnt_pltctrl(ps%pnts, sp%nt, sp%prnt_frq)
        call prnt_elems_txt(ps%surfs, ps%totsurfs)
        call prnt_txt(u, ps%pnts, ps%totpnts, fname, 0)
    end if
    select case (sp%upwind)
    case(.false.)
        do t = 1,sp%nt
            select case(sp%krnl)
            case (1:5)
                call trnsprt_no_upwind(u, krnls, ps%totpnts, d, vof, vel, q, sp%dt)
                !call trnsprt_adam_bashford(u, rhs, krnls, ps%totpnts, d, vof, vel, q, sp%dt, t)
            case (6)   
                call trnsprt_gfd(u, krnls, ps%totpnts, d, vof, vel, q, sp%dt)
            end select
            call set_var(u, bcs)
            if (t/sp%prnt_frq == c) then
                if (sp%vtk .eqv. .true.) then
                    call prnt_vtk(u, ps, fname, t)
                else
                    call prnt_txt(u, ps%pnts, ps%totpnts, fname, t)
                end if
                c = c + 1
            end if
            if (abs(maxval(u)) >= 1000) then
                write(*,*) "ERROR: garbage output"
                call exit(0)
            end if
        end do
    case(.true.)
        do t = 1,sp%nt    
            call trnsprt_upwind(u, ps%pnts, krnls, ps%totpnts, d, vof, vel, q, sp%dt, ps%dim, 0.75d0)   
            call set_var(u, bcs)
            if (t/sp%prnt_frq == c) then
                if (sp%vtk .eqv. .true.) then
                    call prnt_vtk(u, ps, fname, t)
                else
                    call prnt_txt(u, ps%pnts, ps%totpnts, fname, t)
                end if
                c = c + 1
            end if
        end do
    end select

end subroutine slv_trnsprt

end module trnsprt_slvrs
