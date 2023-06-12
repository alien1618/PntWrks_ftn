module bndry
implicit none
contains

pure subroutine set_dbc(pnts, totpnts, edge, loc, value, bc_out)
!-----------------------------------------------------------------------
!  assign points with a constant boundary condition value
!-----------------------------------------------------------------------
    use bc_struct
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) ::  pnts
    integer, intent(in)                 ::  totpnts, edge
    real(8), intent(in)                 ::  loc, value
!-----------------------------------------------------------------------
    type(bc), intent(inout)             ::  bc_out
!-----------------------------------------------------------------------
    integer                             ::  i
    logical                             ::  cond1, cond2
!-----------------------------------------------------------------------
    select case (edge)
    case (1)    !x axis
        do i = 1, totpnts
            cond1 = pnts(i, 1) >= loc - 0.01 * abs(loc)
            cond2 = pnts(i, 1) <= loc + 0.01 * abs(loc)
            if (cond1 .and. cond2) then
                bc_out%tot = bc_out%tot + 1
                bc_out%pnts = [bc_out%pnts, i]
                bc_out%vals = [bc_out%vals, value]
            end if
        end do
    case (2)    !y axis
        do i = 1, totpnts
            cond1 = pnts(i, 2) >= loc - 0.01 * abs(loc)
            cond2 = pnts(i, 2) <= loc + 0.01 * abs(loc)
            if (cond1 .and. cond2) then
                bc_out%tot = bc_out%tot + 1
                bc_out%pnts = [bc_out%pnts, i]
                bc_out%vals = [bc_out%vals, value]
            end if
        end do
    case (3)    !z axis
        do i = 1, totpnts
            cond1 = pnts(i, 3) >= loc - 0.01 * abs(loc)
            cond2 = pnts(i, 3) <= loc + 0.01 * abs(loc)
            if (cond1 .and. cond2) then
                bc_out%tot = bc_out%tot + 1
                bc_out%pnts = [bc_out%pnts, i]
                bc_out%vals = [bc_out%vals, value]
            end if
        end do
    end select
end subroutine set_dbc

pure subroutine set_dbc2(pnts, totpnts, edge, loc, values, bc_out)
!-----------------------------------------------------------------------
! assign points with variable boundary condition values
!-----------------------------------------------------------------------
    use bc_struct
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) ::  pnts
    integer, intent(in)                 ::  totpnts, edge
    real(8), intent(in)                 ::  loc
    real(8), dimension(:), intent(in)   ::  values
!-----------------------------------------------------------------------
    type(bc), intent(inout)             ::  bc_out
!-----------------------------------------------------------------------
    integer                             ::  i
    logical                             ::  cond1, cond2
!-----------------------------------------------------------------------
    select case (edge)
    case (1)    !x axis
        do i = 1, totpnts
            cond1 = pnts(i, 1) >= loc - 0.01 * abs(loc)
            cond2 = pnts(i, 1) <= loc + 0.01 * abs(loc)
            if (cond1 .and. cond2) then
                bc_out%tot = bc_out%tot + 1
                bc_out%pnts = [bc_out%pnts, i]
                bc_out%vals = [bc_out%vals, values(i)]
            end if
        end do
    case (2)    !y axis
        do i = 1, totpnts
            cond1 = pnts(i, 2) >= loc - 0.01 * abs(loc)
            cond2 = pnts(i, 2) <= loc + 0.01 * abs(loc)
            if (cond1 .and. cond2) then
                bc_out%tot = bc_out%tot + 1
                bc_out%pnts = [bc_out%pnts, i]
                bc_out%vals = [bc_out%vals, values(i)]
            end if
        end do
    case (3)    !z axis
        do i = 1, totpnts
            cond1 = pnts(i, 3) >= loc - 0.01 * abs(loc)
            cond2 = pnts(i, 3) <= loc + 0.01 * abs(loc)
            if (cond1 .and. cond2) then
                bc_out%tot = bc_out%tot + 1
                bc_out%pnts = [bc_out%pnts, i]
                bc_out%vals = [bc_out%vals, values(i)]
            end if
        end do
    end select
end subroutine set_dbc2

subroutine read_dbc(bc_fname, value, mshpnts, mshtotpnts, bc_out)
!-----------------------------------------------------------------------
! read boundary points and assign them a fixed value
!-----------------------------------------------------------------------
    use pntst
    use pntst_struct
    use bc_struct
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: mshpnts
    integer, intent(in)                 :: mshtotpnts
    real(8), intent(in)                 :: value
    character(len=50), intent(in)       :: bc_fname
!-----------------------------------------------------------------------
    type(bc), intent(inout)             :: bc_out
!-----------------------------------------------------------------------
    integer                             :: i, j
    real(8)                             :: dx, dy, dz, d, tol = 1e-4
    type(pointset)                      :: bcs
!-----------------------------------------------------------------------
    call read_pnts(bc_fname, bcs)
    do i = 1, bcs%totpnts
        do j = 1, mshtotpnts
            dx = bcs%pnts(i, 1) - mshpnts(j, 1)
            dy = bcs%pnts(i, 2) - mshpnts(j, 2)
            dz = bcs%pnts(i, 3) - mshpnts(j, 3)
            d = sqrt(dx * dx + dy * dy + dz * dz)
            if (d <= tol) then
                bc_out%tot = bc_out%tot + 1
                bc_out%pnts = [bc_out%pnts, j]
                bc_out%vals = [bc_out%vals, value]
            end if
        end do
    end do
end subroutine

pure subroutine set_var(u, bcs)
    !-----------------------------------------------------------
    ! set boundary condition values to the field variable u
    !-----------------------------------------------------------
    use bc_struct
    real(8), dimension(:), intent(INOUT)   ::  u
    type(bc), intent(in)                   :: bcs
    
    if (bcs%tot > 0) u(bcs%pnts(:)) = bcs%vals(:)
end subroutine set_var

end module bndry
