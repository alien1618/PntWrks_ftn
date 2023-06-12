module sdf_slvrs
implicit none
contains

subroutine run_sdf()
!------------------------------------------------------------------------------------------
!   This subroutine solves the allen-cahn phase-field equations for dendritic
!   solidification in a pure liquid. Solution is obtained using explicit time 
!   stepping schemes and local strong-form meshfree methods.     
!------------------------------------------------------------------------------------------
    use pntst_struct
    use slvr_prmtrs_struct
    use mat_struct
    use bc_struct
    use pntst
    use bndry
    use intrf
    use trnsprt
    use sdf
    use omp_lib
    use prmtrs
!------------------------------------------------------------------------------------------
    type(pointset)                          :: ps, ip       !pointset data structure
    real(8), dimension(:), allocatable      :: phi      !phase-field parameter distribution
    type(slvr_prmtrs)                       :: sp       !solver parameters data structures
    real(8)                                 :: start, finish !timer parameters
    character(len=50)                       :: ip_fname, phi_fname = 'phi'
!------------------------------------------------------------------------------------------
    !read interface points
    ip_fname = 'sim/in/1_pntst/intrf.txt'
    call read_pnts(ip_fname, ip)  
    
    !construct pointset
    call set_geom(ps)  
    
    call cpu_time(start)
    !call get_df(ip%pnts, ip%totpnts, ps%pnts, ps%totpnts, phi)
    call set_krnl(ps%dx, sp)
    call get_sdf_vec(ip%pnts, ip%totpnts, ps%pnts, ps%totpnts, ps%dim, ps%dx, sp, phi)
    !call get_sdf_norm(ip%pnts, ip%totpnts, ps%pnts, ps%totpnts, ps%dim, ps%dx, sp, phi)
    call cpu_time(finish)
    !finish = omp_get_wtime()
    print '("Time = ",f8.1," seconds.")',finish-start
    !call prnt_txt(phi, ps%pnts, ps%totpnts, phi_fname, 0)
    call prnt_vtk(phi, ps, phi_fname, 0)
end subroutine run_sdf

end module sdf_slvrs
