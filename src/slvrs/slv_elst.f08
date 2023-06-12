
module elst_slvrs
implicit none
contains

subroutine run_elst()
!------------------------------------------------------------------------------------------
! This subroutine solves the transient elasticity equations using
! local strong-form meshfree methods and iterative explicit time stepping schemes                                                                                                                              
!------------------------------------------------------------------------------------------
    use pntst_struct
    use slvr_prmtrs_struct
    use mat_struct
    use bc_struct
    use pntst
    use bndry
    use prmtrs
    use omp_lib
!------------------------------------------------------------------------------------------
    type(pointset)                          :: ps        !pointset
    type(slvr_prmtrs)                       :: sp        !solver parameters data structure
    type(materials)                         :: mat
    type(bc), dimension(6)                  :: bcs       !boundary condition data structure for vx, vy, vz, p
    real(8)                                 :: start, finish !timer parameters
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Running elsticity solver..."
    !constructing the pointset
    call set_geom(ps)
    !assigning the material properties
    call set_mat('moe', mat%moe)
    call set_mat('nu', mat%nu)
    call set_mat('ro', mat%ro)
    !assigning the boundary conditions
    call set_bc(ps%pnts, ps%totpnts, bcs)
    !assigning the solver controls
    call set_krnl(ps%dx, sp)
    call set_t(sp)
    !running the solver
    start = omp_get_wtime()
    call slv_elst(ps, mat, bcs, sp)
    finish = omp_get_wtime()
    print '("Computation time = ",f8.1," seconds.")',finish-start
end subroutine run_elst

subroutine  slv_elst(ps, mat, bcs, sp)
!------------------------------------------------------------------------------------------
    use pntst_struct
    use krnl_struct
    use mat_struct
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
    use prmtrs
    use ns_lgr_slvrs 
!------------------------------------------------------------------------------------------
    real(8)                                    :: moe, nu, ro
    type(materials), intent(in)                :: mat
    real(8), dimension(:,:), allocatable       :: f
    type(bc), dimension(6), intent(in)         :: bcs
    type(slvr_prmtrs), intent(in)              :: sp
!------------------------------------------------------------------------------------------
    type(pointset), intent(inout)              :: ps
!------------------------------------------------------------------------------------------
    type(kernel), dimension(:), allocatable    :: krnls
    character(len=50)                          :: fname = 'u', vms_fname = 'vms'
    integer                                    :: c = 1, i, j, t, nbr
    real(8), dimension(ps%totpnts, 3)          :: vel, vel_new
    real(8), dimension(ps%totpnts)             :: v, vms
    real(8)                                    :: eps = 1.0e-06, t1, t2, start, finish
    real(8)                                    :: dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz
    real(8)                                    :: a_pressure_x, a_pressure_y, a_pressure_z
    real(8)                                    :: art_visc_x, art_visc_y, art_visc_z
    real(8)                                    :: strain_xx, strain_yy, strain_zz, strain_xy, strain_xz, strain_yz
    real(8)                                    :: w_xx, w_yy, w_zz, w_xy, w_xz, w_yz
    real(8)                                    :: w_yx, w_zx, w_zy, adv_x, adv_y, adv_z, dx, dy, dz, d
    real(8), dimension(ps%totpnts)             :: stress_xx, stress_yy, stress_zz
    real(8), dimension(ps%totpnts)             :: stress_xy, stress_xz, stress_yz, conc, conc0
    real(8), dimension(ps%totpnts)             :: s_xx, s_yy, s_zz, s_xy, s_xz, s_yz, p
    real(8), dimension(ps%totpnts)             :: s_xx_new, s_yy_new, s_zz_new, s_xy_new, s_xz_new, s_yz_new
    real(8), dimension(:,:), allocatable       :: d_mat
    real(8), dimension(6)                      :: strain
    !real(8) :: xij, yij, zij, vxij, vyij, vzij, uij_rij, alpha, d0, h0, trm0, trm_f, pii0
    real(8) :: nabla2_vx, nabla2_vy, nabla2_vz
!------------------------------------------------------------------------------------------
    vel(:,:) = 0.0d0
    allocate(f(ps%totpnts, 3))
    f(:,:) = 0.0
    f(bcs(4)%pnts, 1) = bcs(4)%vals
    f(bcs(5)%pnts, 2) = bcs(5)%vals
    f(bcs(6)%pnts, 3) = bcs(6)%vals
    moe = mat%moe(1)
    ro = mat%ro(1)
    nu = mat%nu(1)
!------------------------------------------------------------------------------------------
    if (sp%vtk .eqv. .true.) then
        call prnt_vtk(v, ps, fname, 0)
    else
        call prnt_pltctrl(ps%pnts, sp%nt, sp%prnt_frq)
        call prnt_elems_txt(ps%surfs, ps%totsurfs)
        call prnt_txt(v, ps%pnts, ps%totpnts, fname, 0)
    end if
    
    vel(bcs(1)%pnts, 1) = bcs(1)%vals
    vel(bcs(2)%pnts, 2) = bcs(2)%vals
    vel(bcs(3)%pnts, 3) = bcs(3)%vals

    t1 = moe / (1 + nu)
    t2 = nu / (1 - 2 * nu)

    allocate(d_mat(6,6))
    call get_d3d(moe, nu, d_mat)

    !allocate(d_mat(3,3))
    !call get_d2d(1, moe, nu, 0.01d0, d_mat)
    
    s_xx(:) = 0.0
    s_yy(:) = 0.0
    s_zz(:) = 0.0
    s_xy(:) = 0.0
    s_xz(:) = 0.0
    s_yz(:) = 0.0
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Constructing krnls..."
    call cpu_time(start)
    allocate(krnls(ps%totpnts))
    call get_nbrs_bf(sp%h, ps%pnts, ps%totpnts, ps%pnts, ps%totpnts, krnls)
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
    print '("intrpolants construction time = ",f10.1," seconds")',finish-start
!------------------------------------------------------------------------------------------
    do i = 1, ps%totpnts
        conc(i) = 0
        do j = 1, krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            dx = ps%pnts(i,1) - ps%pnts(nbr,1)
            dy = ps%pnts(i,2) - ps%pnts(nbr,2)
            dz = ps%pnts(i,3) - ps%pnts(nbr,3)
            d = sqrt(dx * dx + dy * dy + dz * dz) + 1e-6
            conc0(i) = conc0(i) + (((sp%h/d) - 1))
        end do
    end do
!------------------------------------------------------------------------------------------
    write(*,'(a)') "Computing transient solution..."
    do t = 1,sp%nt    
            do i = 1, ps%totpnts
            conc(i) = 0
            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                dx = ps%pnts(i,1) - ps%pnts(nbr,1)
                dy = ps%pnts(i,2) - ps%pnts(nbr,2)
                dz = ps%pnts(i,3) - ps%pnts(nbr,3)
                d = sqrt(dx * dx + dy * dy + dz * dz) + 1e-6
                conc(i) = conc(i) + (((sp%h / d) - 1))
            end do
            p(i) = 0.001 * (conc(i) - conc0(i))
        end do
        
        do i = 1, ps%totpnts
            dvxdx = 0.0
            dvxdy = 0.0
            dvxdz = 0.0
            dvydx = 0.0
            dvydy = 0.0
            dvydz = 0.0
            dvzdx = 0.0
            dvzdy = 0.0
            dvzdz = 0.0
            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                dvxdx = dvxdx + krnls(i)%nx(j) * (vel(nbr,1) - vel(i,1))
                dvxdy = dvxdy + krnls(i)%ny(j) * (vel(nbr,1) - vel(i,1))
                dvxdz = dvxdz + krnls(i)%nz(j) * (vel(nbr,1) - vel(i,1))

                dvydx = dvydx + krnls(i)%nx(j) * (vel(nbr,2) - vel(i,2))
                dvydy = dvydy + krnls(i)%ny(j) * (vel(nbr,2) - vel(i,2))
                dvydz = dvydz + krnls(i)%nz(j) * (vel(nbr,2) - vel(i,2))

                dvzdx = dvzdx + krnls(i)%nx(j) * (vel(nbr,3) - vel(i,3))
                dvzdy = dvzdy + krnls(i)%ny(j) * (vel(nbr,3) - vel(i,3))
                dvzdz = dvzdz + krnls(i)%nz(j) * (vel(nbr,3) - vel(i,3))
            end do
            strain_xx = dvxdx
            strain_yy = dvydy
            strain_zz = dvzdz
            strain_xy = 0.5*(dvxdy + dvydx)
            strain_xz = 0.5*(dvxdz + dvzdx)
            strain_yz = 0.5*(dvydz + dvzdy)

            w_xx = 0.0
            w_yy = 0.0
            w_zz = 0.0
            w_xy = 0.5*(dvxdy - dvydx)
            w_xz = 0.5*(dvxdz - dvzdx)
            w_yz = 0.5*(dvydz - dvzdy)
            w_yx = 0.5*(dvydx - dvxdy)
            w_zx = 0.5*(dvzdx - dvxdz)
            w_zy = 0.5*(dvzdy - dvydz)

            strain = [strain_xx, strain_yy, strain_zz, strain_xy, strain_xz, strain_yz]
            s_xx_new(i) = s_xx(i) + sp%dt * (sum(d_mat(:,1) * strain))
            s_yy_new(i) = s_yy(i) + sp%dt * (sum(d_mat(:,2) * strain))
            s_zz_new(i) = s_zz(i) + sp%dt * (sum(d_mat(:,3) * strain))
            s_xy_new(i) = s_xy(i) + sp%dt * (strain_xy * d_mat(4,4) + w_yz * s_xz(i) + w_xz * s_yz(i))
            s_xz_new(i) = s_xz(i) + sp%dt * (strain_xz * d_mat(5,5) + w_zy * s_xy(i) + w_xy * s_yz(i))
            s_yz_new(i) = s_yz(i) + sp%dt * (strain_yz * d_mat(6,6) + w_zx * s_xy(i) + w_yx * s_xz(i))

            !s_xx_new(i) = s_xx(i) + sp%dt * (t1 * (t2 * strain_xx + strain_xx))
            !s_yy_new(i) = s_yy(i) + sp%dt * (t1 * (t2 * strain_yy + strain_yy))
            !s_zz_new(i) = s_zz(i) + sp%dt * (t1 * (t2 * strain_zz + strain_zz))
            !s_xy_new(i) = s_xy(i) + sp%dt * (t1 * strain_xy + w_yz * s_xz(i) + w_xz * s_yz(i))
            !s_xz_new(i) = s_xz(i) + sp%dt * (t1 * strain_xz + w_zy * s_xy(i) + w_xy * s_yz(i))
            !s_yz_new(i) = s_yz(i) + sp%dt * (t1 * strain_yz + w_zx * s_xy(i) + w_yx * s_xz(i))

        end do
        s_xx = s_xx_new
        s_yy = s_yy_new
        s_zz = s_zz_new
        s_xy = s_xy_new
        s_xz = s_xz_new
        s_yz = s_yz_new
        stress_xx = s_xx - p
        stress_yy = s_yy - p
        stress_zz = s_zz - p
        stress_xy = s_xy 
        stress_xz = s_xz 
        stress_yz = s_yz 

        vms = sqrt(((stress_xx - stress_yy)**2) + ((stress_yy - stress_zz)**2) + ((stress_xx - stress_zz)**2) &
                + 6 * (stress_xy * stress_xy + stress_yz * stress_yz + stress_xz * stress_xz)) / sqrt(2.0)

        do i = 1, ps%totpnts
            A_pressure_x = 0
            A_pressure_y = 0
            A_pressure_z = 0
            art_visc_x = 0
            art_visc_y = 0
            art_visc_z = 0

            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)

                !------------------------------------------------------------------------------------------
                ! computing acceleration due to pressure
                !------------------------------------------------------------------------------------------
                A_pressure_x = A_pressure_x + (stress_xx(i) + stress_xx(nbr)) * krnls(i)%nx(j)    !do not touch you idiot
                A_pressure_x = A_pressure_x + (stress_xy(i) + stress_xy(nbr)) * krnls(i)%ny(j)    !do not touch you idiot
                A_pressure_x = A_pressure_x + (stress_xz(i) + stress_xz(nbr)) * krnls(i)%nz(j)    !do not touch you idiot
                
                A_pressure_y = A_pressure_y + (stress_yy(i) + stress_yy(nbr)) * krnls(i)%ny(j)    !do not touch you idiot
                A_pressure_y = A_pressure_y + (stress_xy(i) + stress_xy(nbr)) * krnls(i)%nx(j)    !do not touch you idiot
                A_pressure_y = A_pressure_y + (stress_yz(i) + stress_yz(nbr)) * krnls(i)%nz(j)    !do not touch you idiot
                
                A_pressure_z = A_pressure_z + (stress_zz(i) + stress_zz(nbr)) * krnls(i)%nz(j)    !do not touch you idiot
                A_pressure_z = A_pressure_z + (stress_xz(i) + stress_xz(nbr)) * krnls(i)%nx(j)    !do not touch you idiot
                A_pressure_z = A_pressure_z + (stress_yz(i) + stress_yz(nbr)) * krnls(i)%ny(j)    !do not touch you idiot

                !------------------------------------------------------------------------------------------
                ! computing acceleration due to artificial viscosity
                !------------------------------------------------------------------------------------------
                !xij = ps%pnts(i,1) - ps%pnts(nbr,1)
                !yij = ps%pnts(i,2) - ps%pnts(nbr,2)
                !zij = ps%pnts(i,3) - ps%pnts(nbr,3)
                !vxij = vel(i,1) - vel(nbr,1)
                !vyij = vel(i,2) - vel(nbr,2)
                !vzij = vel(i,3) - vel(nbr,3)
                !uij_rij = xij * vxij + yij * vyij + zij * vzij
  
                !alpha = 1
                !d0 = sqrt(xij * xij + yij * yij + zij * zij);
                !h0 = 2 * d0;
                !trm0 = ((alpha * c * h0) / ro) / (d0 * d0 + h0 * 1e-5)
                !if (lagrangian .eqv. .true.) then
                !    trm0 = ((alpha * c * h0))/(d0 * d0 + h0 * 1e-5)
                !end if
                !trm_f = -uij_rij * trm0
                !pii0 = trm_f
                !if (trm_f < 0) pii0 = 0
                !art_visc_x = art_visc_x + pii0 * krnls(i)%nx(j)
                !art_visc_y = art_visc_y + pii0 * krnls(i)%ny(j)
                !art_visc_z = art_visc_z + pii0 * krnls(i)%nz(j)
            end do

            dvxdx = 0.0
            dvxdy = 0.0
            dvxdz = 0.0
            dvydx = 0.0
            dvydy = 0.0
            dvydz = 0.0
            dvzdx = 0.0
            dvzdy = 0.0
            dvzdz = 0.0
            nabla2_vx = 0.0
            nabla2_vy = 0.0
            nabla2_vz = 0.0
            do j = 1, krnls(i)%totnbrs
                nbr = krnls(i)%nbrs(j)
                dvxdx = dvxdx + krnls(i)%nx(j) * (vel(nbr,1) - vel(i,1))
                dvxdy = dvxdy + krnls(i)%ny(j) * (vel(nbr,1) - vel(i,1))
                dvxdz = dvxdz + krnls(i)%nz(j) * (vel(nbr,1) - vel(i,1))

                dvydx = dvydx + krnls(i)%nx(j) * (vel(nbr,2) - vel(i,2))
                dvydy = dvydy + krnls(i)%ny(j) * (vel(nbr,2) - vel(i,2))
                dvydz = dvydz + krnls(i)%nz(j) * (vel(nbr,2) - vel(i,2))

                dvzdx = dvzdx + krnls(i)%nx(j) * (vel(nbr,3) - vel(i,3))
                dvzdy = dvzdy + krnls(i)%ny(j) * (vel(nbr,3) - vel(i,3))
                dvzdz = dvzdz + krnls(i)%nz(j) * (vel(nbr,3) - vel(i,3))

                !nabla2_vx = nabla2_vx + krnls(i)%nabla2(j) * (vel(nbr,1) - vel(i,1))
                !nabla2_vy = nabla2_vy + krnls(i)%nabla2(j) * (vel(nbr,2) - vel(i,2))
                !nabla2_vz = nabla2_vz + krnls(i)%nabla2(j) * (vel(nbr,3) - vel(i,3))
            end do
            adv_x = 0.0
            adv_y = 0.0
            adv_z = 0.0
            
            !adv_x = vel(i,1) * dvxdx + vel(i,2) * dvxdy + vel(i,3) * dvxdz
            !adv_y = vel(i,1) * dvydx + vel(i,2) * dvydy + vel(i,3) * dvydz
            !adv_z = vel(i,1) * dvzdx + vel(i,2) * dvzdy + vel(i,3) * dvzdz

            vel_new(i,1) = vel(i,1) + sp%dt * ((1 / (ro * ro + eps)) * A_pressure_x &
                - art_visc_x - adv_x + f(i,1) + 0.1 * nabla2_vx)
            vel_new(i,2) = vel(i,2) + sp%dt * ((1 / (ro * ro + eps)) * A_pressure_y &
                - art_visc_y -  adv_y + f(i,2) + 0.1 * nabla2_vy)
            vel_new(i,3) = vel(i,3) + sp%dt * ((1 / (ro * ro + eps)) * A_pressure_z &
                - art_visc_z - adv_z + f(i,3) + 0.1 * nabla2_vz)
        end do
        
        !------------------------------------------------------------------------------------------
        ! computing final velocity
		!------------------------------------------------------------------------------------------
        vel_new(bcs(1)%pnts, 1) = bcs(1)%vals
        vel_new(bcs(2)%pnts, 2) = bcs(2)%vals
        vel_new(bcs(3)%pnts, 3) = bcs(3)%vals
        vel = vel_new
        v = sqrt(vel(:,1) * vel(:,1) + vel(:,2) * vel(:,2) + vel(:,3) * vel(:,3))

        do i = 1, ps%totpnts
            ps%pnts(i,:) = ps%pnts(i,:) + (vel(i,:)) * sp%dt
        end do

        if (t/sp%prnt_frq == c) then
            if (sp%vtk .eqv. .true.) then
                call prnt_vtk(v, ps, fname, t)
                call prnt_vtk(vms, ps, vms_fname, t)
            else
                call prnt_txt(v, ps%pnts, ps%totpnts, fname, t)
                call prnt_vtk(vms, ps, vms_fname, t)
            end if
            c = c + 1
        end if
    end do
end subroutine slv_elst

subroutine get_d2d(load_type, moe, nu, t, dmat)
!------------------------------------------------------------------------------------------
    real(8), intent(IN)                   :: moe, nu, t
    integer, intent(IN)                   :: load_type
!------------------------------------------------------------------------------------------
    real(8), dimension(3,3), intent(OUT)  :: dmat
!------------------------------------------------------------------------------------------
    real(8)                               :: d1 = 0.0, d2 = 0.0, d3 = 0.0
!------------------------------------------------------------------------------------------
    select case (load_type)
        case (1) !plane stress
            d1 = moe * t / (1 - (nu * nu))
            d2 = nu * D1
            d3 = D1 * (1 - nu) / 2
        case (2) !plane strain
            d1 = (1 - nu) * moe * t / ((1 + nu) * (1 - 2 * nu))
            d2 = nu * moe * t / ((1 + nu) * (1 - 2 * nu))
            d3 = (1 - 2 * nu) * moe * t / 2 * ((1 + nu) * (1 - 2 * nu))
        case default
            write(*,'(a)') "Error: load type is invalid"
    end select
    dmat(1,1) = d1; dmat(2,1) = d2; dmat(3,1) = 0;
    dmat(1,2) = d2; dmat(2,2) = d1; dmat(3,2) = 0;
    dmat(1,3) = 0; dmat(2,3) = 0; dmat(3,3) = d3;
end subroutine get_d2d

subroutine get_d3d(moe, nu, dmat)
!------------------------------------------------------------------------------------------
    real(8), intent(IN)                   :: moe, nu
!------------------------------------------------------------------------------------------
    real(8), dimension(6,6), intent(OUT)  :: dmat
!------------------------------------------------------------------------------------------
    real(8)                               :: t1, d1, d2, d3
!------------------------------------------------------------------------------------------
    t1 = (moe / ((1 + nu) * (1 - 2 * nu)))
    D1 = (1 - nu) * t1
    D2 = nu * t1
    D3 = 0.5 * (1 - 2 * nu) * t1
    dmat(1,1) = D1; dmat(2,1) = D2; dmat(3,1) = D2;
    dmat(1,2) = D2; dmat(2,2) = D1; dmat(3,2) = D2;
    dmat(1,3) = D2; dmat(2,3) = D2; dmat(3,3) = D1;
    dmat(4,4) = D3; dmat(5,5) = D3; dmat(6,6) = D3;
end subroutine get_d3d
end module elst_slvrs
