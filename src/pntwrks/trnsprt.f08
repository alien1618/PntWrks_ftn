module trnsprt
implicit none
contains

subroutine trnsprt_no_upwind(u, krnls, totpnts, d, vof, vel, q, dt)
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use omp_lib
!-------------------------------------------------------------------------
    integer, intent(in)                        ::  totpnts
    real(8), intent(in)                        ::  dt
    real(8), dimension(:), intent(in)          ::  d, vof
    type(kernel), dimension(:), intent(in)     ::  krnls
    real(8), dimension(:), intent(in)          ::  q
    real(8), dimension(:,:), intent(in)        ::  vel
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)       ::  u
!-------------------------------------------------------------------------
    real(8), dimension(totpnts)                ::  u_new
    integer                                    ::  i
    real(8)                                    ::  k, ux, uy, uz, adv, nabla2_u
!-------------------------------------------------------------------------
     
    u_new = u  
   ! !$OMP PARALLEL DO private(i)
    do i = 1,totpnts
        call get_grad(u(i), u(krnls(i)%nbrs), krnls(i), ux, uy, uz, nabla2_u)          
        adv = ux * vel(i,1) + uy * vel(i,2) + uz * vel(i,3)
        k = vof(i) * D(1) + (1 - vof(i)) * D(2)
        u_new(i) = u(i) + dt * (k * nabla2_u + q(i) - adv)
    end do
   ! !$OMP END PARALLEL DO
    u = u_new
end subroutine trnsprt_no_upwind

subroutine trnsprt_upwind(u, pnts, krnls, totpnts, d, vof, vel, q, dt, dim, upwind_ratio)
!-------------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!-------------------------------------------------------------------------
    integer, intent(in)                     ::  totpnts, dim
    real(8), intent(in)                        ::  dt, upwind_ratio
    real(8), dimension(:), intent(in)          ::  d, vof
    real(8), dimension(:,:), intent(in)         :: pnts
    type(kernel), dimension(:), intent(in)  ::  krnls
    real(8), dimension(:), intent(in)          ::  q
    real(8), dimension(:,:), intent(in)        ::  vel
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)       ::  u
!-------------------------------------------------------------------------
    real(8), dimension(totpnts)            ::  u_new
    integer                                ::  i, j, nbr
    real(8)                                ::  k, check, du, adv, nabla2_u
    real(8)                                ::  ux, uy, uz, ux_up, uy_up, uz_up, ux_dn, uy_dn, uz_dn
    real(8)                                ::  xj_xi, yj_yi, zj_zi 
!-------------------------------------------------------------------------
    u_new = u
    !$OMP PARALLEL DO private(i)
    do i = 1,totpnts
        ux_up = 0;
        uy_up = 0;
        uz_up = 0;
        ux_dn = 0;
        uy_dn = 0;
        uz_dn = 0;
        nabla2_u = 0
        do j = 1,krnls(i)%totnbrs
            nbr = krnls(i)%nbrs(j)
            xj_xi = pnts(nbr,1) - pnts(i,1)
            yj_yi = pnts(nbr,2) - pnts(i,2)         
            zj_zi = pnts(nbr,3) - pnts(i,3)

            check = vel(i,1) * xj_xi + vel(i,2) * yj_yi + vel(i,3) * zj_zi
            du = u(nbr) - u(i)
            if (check < 0) then !upwind
                ux_up = ux_up + du * krnls(i)%nx(j)
                uy_up = uy_up + du * krnls(i)%ny(j)
                uz_up = uz_up + du * krnls(i)%nz(j)
            end if
            if (check >= 0) then ! downwind
                ux_dn = ux_dn + du * krnls(i)%nx(j)
                uy_dn = uy_dn + du * krnls(i)%ny(j)
                uz_dn = uz_dn + du * krnls(i)%nz(j)
            end if
            nabla2_u = nabla2_u + du * krnls(i)%nabla2(j)
        end do
        ux = dim * (upwind_ratio * ux_up + (1 - upwind_ratio) * ux_dn)
        uy = dim * (upwind_ratio * uy_up + (1 - upwind_ratio) * uy_dn)
        uz = dim * (upwind_ratio * uz_up + (1 - upwind_ratio) * uz_dn)
        adv = ux * vel(i,1) + uy * vel(i,2) + uz * vel(i,3)
        k = vof(i) * D(1) + (1 - vof(i)) * D(2)
        u_new(i) = u(i) + dt * (k * nabla2_u + q(i) - adv)
    end do
    !$OMP END PARALLEL DO
    
    u = u_new
end subroutine trnsprt_upwind

subroutine trnsprt_adam_bashford(u, rhs, krnls, totpnts, d, vof, vel, q, dt, t)
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use omp_lib
!-------------------------------------------------------------------------
    integer, intent(in)                     ::  totpnts, t
    real(8), intent(in)                     ::  dt
    real(8), dimension(:), intent(in)       ::  d, vof
    type(kernel), dimension(:), intent(in)  ::  krnls
    real(8), dimension(:), intent(in)       ::  q
    real(8), dimension(:,:), intent(in)     ::  vel
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)    ::  u
    real(8), dimension(:), intent(inout)    :: rhs
!-------------------------------------------------------------------------
    real(8), dimension(totpnts)         ::  rhs_old
    integer                             ::  i
    real(8)                             ::  k, adv, ux, uy, uz, nabla2_u
!-------------------------------------------------------------------------
    
    rhs_old = rhs
    ! !$OMP PARALLEL DO private(i)
    do i = 1,totpnts
        call get_grad(u(i), u(krnls(i)%nbrs), krnls(i), ux, uy, uz, nabla2_u)          
        adv = ux * vel(i,1) + uy * vel(i,2) + uz * vel(i,3)
        k = vof(i) * D(1) + (1 - vof(i)) * D(2)
        rhs(i) = (k * nabla2_u + q(i) - adv)
    end do
    ! !$OMP END PARALLEL DO   
    
    if (t < 2) then
        u = u + dt*rhs
    else
        u = u + (dt/2.0d0) * (3.0d0 * rhs - rhs_old)
    end if
end subroutine trnsprt_adam_bashford

subroutine get_trnsprt_gauss_seidel(u, krnls, totpnts, d, vof, vel, q, dt)
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
!-------------------------------------------------------------------------
    integer, intent(in)                     ::  totpnts
    real(8), intent(in)                     ::  dt
    type(kernel), dimension(:), intent(in)  ::  krnls
    real(8), dimension(:), intent(in)       ::  q, vof
    real(8), dimension(:,:), intent(in)     ::  vel
    real(8), dimension(:), intent(in)       ::  D
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)    ::  u
!-------------------------------------------------------------------------
    integer                                 ::  i
    real(8)                                 ::  rhs, k, ux, uy, uz, nabla2_u, beta = 3
!-------------------------------------------------------------------------

    !$OMP PARALLEL DO private(i)   
    do i = 1, totpnts
        call get_grad(u(i), u(krnls(i)%nbrs), krnls(i), ux, uy, uz, nabla2_u)          
        k = vof(i) * D(1) + (1 - vof(i)) * D(2)
        RHS = (K * nabla2_u - (vel(i,1) * ux + vel(i,2) * uy + vel(i,3) * uz)) + q(i)
        u(i) = beta*(u(i) + dt * RHS) + (1 - beta) * u(i)
    end do
    !$OMP END PARALLEL DO
end subroutine get_trnsprt_gauss_seidel

subroutine trnsprt_implicit(u, krnls, totpnts, d, vof, vel, q, dt, itrs)
!-------------------------------------------------------------------------
    use krnl_struct
    use krnl_cmn
    use omp_lib
!-------------------------------------------------------------------------
    integer, intent(in)                     ::  totpnts, itrs
    real(8), intent(in)                     ::  dt
    real(8), dimension(:), intent(in)       ::  d, vof
    type(kernel), dimension(:), intent(in)  ::  krnls
    real(8), dimension(:), intent(in)       ::  q
    real(8), dimension(:,:), intent(in)     ::  vel
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)    ::  u
!-------------------------------------------------------------------------
    real(8), dimension(totpnts)  ::  u_itr
    integer                      ::  i, m
    real(8)                      ::  k, ux, uy, uz, adv, nabla2_u
!-------------------------------------------------------------------------
    u_itr = u
    do m = 1, itrs
        ! !$OMP PARALLEL DO private(i)       
        do i = 1,totpnts
            call get_grad(u(i), u(krnls(i)%nbrs), krnls(i), ux, uy, uz, nabla2_u)          
            adv = ux * vel(i,1) + uy * vel(i,2) + uz * vel(i,3)
            k = vof(i) * D(1) + (1 - vof(i)) * D(2)
            u_itr(i) = u(i) + dt * (k * nabla2_u + q(i) - adv)
        end do
        ! !$OMP END PARALLEL DO
    end do
    u = u_itr
end subroutine trnsprt_implicit

pure subroutine trnsprt_gfd(u, krnls, totpnts, d, vof, vel, q, dt)
!-------------------------------------------------------------------------
    use krnl_struct
    use omp_lib
!-------------------------------------------------------------------------
    integer, intent(in)                     ::  totpnts
    real(8), intent(in)                     ::  dt
    real(8), dimension(:), intent(in)       ::  d, vof
    type(kernel), dimension(:), intent(in)  ::  krnls
    real(8), dimension(:), intent(in)       ::  q
    real(8), dimension(:,:), intent(in)     ::  vel
!-------------------------------------------------------------------------
    real(8), dimension(:), intent(inout)    ::  u
!-------------------------------------------------------------------------
    real(8), dimension(totpnts)  ::  u_new
    integer                      ::  i, s1, s2
    real(8)                      ::  k, ux, uy, uz, adv, nabla2_u
    real(8), dimension(9)        :: b, grad
!-------------------------------------------------------------------------

    u_new = u   
    do i = 1,totpnts
        do s1 = 1, 9
            !b(s1) = -u(i) * sum(krnls(i)%gfd_trm(s1,:)) + sum(u(krnls(i)%nbrs)*krnls(i)%gfd_trm(s1,:))
            b(s1) = 0.0
            do s2 = 1, krnls(i)%totnbrs
                b(s1) = b(s1) + krnls(i)%gfd_trm(s1, s2) * (u(krnls(i)%nbrs(s2)) - u(i))
            end do
        end do
        !b = matmul(krnls(i)%gfd_trm,(u(krnls(i)%nbrs)-u(i)))
        
        grad = matmul(krnls(i)%gfd_inv_a,b)
        ux = grad(1)
        uy = grad(2)
        uz = grad(3)
        nabla2_u = sum(grad(4:6))
        adv = ux * vel(i,1) + uy * vel(i,2) + uz * vel(i,3)
        k = vof(i) * D(1) + (1 - vof(i)) * D(2)
        u_new(i) = u(i) + dt * (k * nabla2_u + q(i) - adv)
    end do
    u = u_new
end subroutine trnsprt_gfd

end module trnsprt
