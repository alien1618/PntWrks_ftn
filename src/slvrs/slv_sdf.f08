module sdf_slvrs
implicit none
contains

subroutine run_sdf()
!------------------------------------------------------------------------------------------
! For a pointset ps, this subroutine computes a signed distance function whose zero-isocontour
! correspond to a set of interface points ip.
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
    type(pointset)                          :: ps, ip   !pointset data structure
    real(8), dimension(:), allocatable      :: phi      !signed distance function
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
