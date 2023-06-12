program pntwrks
use intrp_slvrs
use trnsprt_slvrs
use trnsprt_ss_slvrs
use ls_slvrs
use ch_slvrs
use ac_slvrs
use ns_1phs_slvrs
use ns_2phs_slvrs
use ns_lgr_slvrs
use elst_slvrs
use sdf_slvrs 
use stfn_slvrs
use prmtrs

implicit none
integer :: n

call set_slvr(n) 
write(*,'(a,i0)') "Solver number selected is: ", n
select case (n)
    case (1)
        call run_sctr_intrp()
    case (2)
        call run_trnsprt_ss()
    case (3)
        call run_trnsprt()
    case (4)
        call run_ls()
    case (5)
        call run_ch()
    case (6)
        call run_ac()
    case (7)
        call run_ns_1phs()
    case (8)  
        call run_ns_2phs()
    case (9)
        call run_ns_1phs_lgr()
    case (10)
        call run_ns_2phs_lgr() 
    case (11)
        call run_elst()
    case (12)
        call run_sctr_reg()
    case (13)
        call run_sdf()
    case (14)
        call run_stfn()
    case default
        write(*,'(a)') "Invalid tst case number"
end select
   
end program pntwrks
