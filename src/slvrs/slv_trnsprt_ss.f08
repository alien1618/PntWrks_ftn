
module trnsprt_ss_slvrs
implicit none
contains

subroutine  run_trnsprt_ss()
!------------------------------------------------------------------------------------------
! This subroutine solves special forms of the general tranport equation using       
! local strong-form meshfree methods and iterative explicit time stepping schemes                                                                                                    
!------------------------------------------------------------------------------------------
    use pntst_struct
    use bc_struct
    use slvr_prmtrs_struct
    use pntst
    use bndry
    use prmtrs
    use omp_lib
!------------------------------------------------------------------------------------------
    real(8), dimension(:), allocatable      :: k         !diffusivity for 2 phases
    real(8), dimension(:), allocatable      :: q         !generation trm
    real(8), dimension(:), allocatable      :: u         !field variable distribution
    real(8), dimension(:,:), allocatable    :: v         !advective velocity in xyz
    real(8), dimension(:), allocatable      :: phi       !volume of fluid distribution (0,1)
    type(pointset)                          :: ps        !point set data structure
    type(bc), dimension(6)                  :: bcs       !boundary conditions data structure
    type(slvr_prmtrs)                       :: sp        !solver parameters data structure
    real(8)                                 :: start, finish !timer parameters
!------------------------------------------------------------------------------------------
    write(*,'(a)') "--------------------------------------------------------"
    write(*,'(a)') "Running steady-state transport solver..."
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
    call slv_trnsprt_ss(ps, bcs(1), sp, k, u, phi, v, q)
    !call cpu_time(finish)
    finish = omp_get_wtime()
    print '("Computation time = ",f8.1," seconds.")',finish-start
end subroutine run_trnsprt_ss

subroutine  slv_trnsprt_ss(ps, bcs, sp, d, u, vof, vel, q)
!------------------------------------------------------------------------------------------
    use pntst_struct
    use krnl_struct
    use bc_struct
    use slvr_prmtrs_struct
    use slvr_cmn
    use pntst
    use bndry
    use trnsprt
!------------------------------------------------------------------------------------------
    real(8),dimension(:), intent(in)           ::  q, vof, d
    real(8), dimension(:,:), intent(in)        ::  vel
    type(pointset), intent(in)                 ::  ps
    type(bc), intent(in)                       ::  bcs
    type(slvr_prmtrs), intent(in)              ::  sp
!------------------------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)       ::  u
!------------------------------------------------------------------------------------------
    type(kernel), dimension(ps%totpnts)        ::  krnls
    character(len=50)                          ::  fname = 'u'
    integer                                    ::  c = 1, t
    real(8)                                    ::  tol = 1.0e-3, res = 0.0, u_sum, u_sum_old
!------------------------------------------------------------------------------------------
    u_sum = sum(u)
    u_sum_old = u_sum
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Constructing krnls..."
    call get_krnls(ps, sp, krnls)
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Computing transient solution..."
    call prnt_vtk(u, ps, fname, 0)
    do t = 1,sp%nt    
        select case(sp%krnl)
            case (1:5)
                !call get_trnsprt_gauss_seidel(u, krnls, ps%totpnts, d, vof, vel, q, sp%dt)
                call trnsprt_no_upwind(u, krnls, ps%totpnts, d, vof, vel, q, sp%dt)
            case (6)   
                call trnsprt_gfd(u, krnls, ps%totpnts, d, vof, vel, q, sp%dt)
        end select

        call set_var(u, bcs)
        u_sum_old = u_sum
        u_sum = sum(u);
        res = abs(u_sum-u_sum_old)
        if (res <= tol) exit
        if (t/sp%prnt_frq == c) then
            write(*,*) "res = ", res
            call prnt_vtk(u, ps, fname, t)
            c = c + 1
        end if
    end do
end subroutine slv_trnsprt_ss

end module trnsprt_ss_slvrs
