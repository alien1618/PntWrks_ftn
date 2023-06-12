module eqslvrs
implicit none
contains

pure subroutine inv(a0, n, b)
!-----------------------------------------------------------------------
! Subroutine calculates the inverse of a square matrix [a0]
! to obtain [b] using gauss elimination
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in) :: a0
    integer, intent(in) :: n
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: b
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable :: a
    integer :: i, j, k, l
    real(8) :: s, t
!-----------------------------------------------------------------------
    allocate(a(n,n))
    allocate(b(n,n))
    a = a0
    b(:, :) = 0.0

    do i = 1, n
      b(i, i) = 1.0
    end do

    !The following code actually performs the matrix inversion
    do j = 1, n
      do i = j, n
        if (a(i, j) /= 0.0) then
          do k = 1, n
            s = a(j, k)
            a(j, k) = a(i, k)
            a(i, k) = s
            s = b(j, k)
            b(j, k) = b(i, k)
            b(i, k) = s
          end do
          t = 1.0 / a(j, j)
          a(j, :) = t * a(j, :)
          b(j, :) = t * b(j, :)
          do L = 1, n
            if (L /= j) then
              t = - a(L, j)
              a(L, :) = a(L, :) + t * a(j, :)
              b(L, :) = b(L, :) + t * b(j, :)
            end if
          end do
        end if
        exit
      end do
    end do
end subroutine inv

subroutine gauss_jordan_inv(a0, n, inv_a)
!-----------------------------------------------------------------------
! Subroutine calculates the inverse of a square matrix
! using gauss-jordan method
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
    integer, intent(in) :: n
    real(8), dimension(:,:), intent(in) :: a0
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable, intent(out) :: inv_a
!-----------------------------------------------------------------------
    real(8), dimension(:,:), allocatable :: a
    integer :: i, j, k
    real(8) :: ratio
!-----------------------------------------------------------------------
    allocate(a(n, 2 * n))
    allocate(inv_a(n, n))
    a(1 : n, 1 : n) = a0
    
    !Augmenting Identity Matrix of Order n
    do i = 1, n
        do j = 1, n
            if(i == j) then
                a(i, j + n) = 1
            else
                a(i, j + n) = 0
            end if
        end do
    end do

    !Applying Gauss Jordan Elimination
    do i = 1, n
        if (a(i, i) == 0.0) then
            write(*,'(a)') "Mathematical Error!"
            !call exit
        end if
        do j = 1, n
            if (i .ne. j) then
                ratio = a(j, i) / a(i, i)
                do k = 1, 2 * n
                    a(j, k) = a(j, k) - ratio * a(i, k)
                end do
            end if
        end do
    end do

    !Row Operation to Make Principal Diagonal to 1
    do i = 1, n
        do j = n + 1, 2 * n
            a(i, j) = a(i, j) / a(i, i)
        end do
    end do

    inv_a = a(1 : n, n + 1 : 2 * n)

end subroutine gauss_jordan_inv

subroutine slv_ge(a0, b, n, x)
!-----------------------------------------------------------
! Subroutine solves a system of equations [a0] {x} = {b}
! to obtain {x} using gaussian elimination
!-----------------------------------------------------------
    implicit none
!-----------------------------------------------------------
    integer, intent(in)                 :: n
    real(8), dimension(:,:), intent(in) :: a0
    real(8), dimension(:), intent(in)   :: b
!-----------------------------------------------------------
    real(8), dimension(n), intent(out) :: x
!-----------------------------------------------------------
    integer :: i, j, k
    real(8), dimension(:,:), allocatable :: a
    real(8) :: ratio
!-----------------------------------------------------------
    allocate(a(n, n+1))
    a(:, 1:n) = a0
    a(:, n+1) = b

    do i = 1, n - 1
        if(a(i, i) == 0.0) then
            write(*,*) "Mathematical Error!"
            call exit
        end if
        do j = i + 1, n
            ratio = a(j, i) / a(i, i)
            do k = 1, n + 1
                a(j, k) = a(j, k) - ratio * a(i, k)
            end do
        end do
    end do
    
    !Obtaining Solution by Back Substitution Method
    x(n) = a(n, n+1) / a(n, n)
    do i = n - 1 , 1, -1
        x(i) = a(i, n + 1)
        do j = i + 1, n
            x(i) = x(i) - a(i, j) * x(j)
        end do
        x(i) = x(i) / a(i, i)
    end do
end subroutine slv_ge

subroutine slv_pcg(A, b, mat_size, cgitr, tol, x)
!-----------------------------------------------------------------------
! Subroutine solves system of equations [A] {x} = {b} with a square
! and sparse  A matrix using pre-conditioned conjugate gradient method.
! Serial code with no OpenMP
!-----------------------------------------------------------------------
    real(8), dimension(:,:), intent(in)   :: A
    real(8), dimension(:), intent(in)     :: b
    integer, intent(in)                   :: mat_size, cgitr
    real(8), intent(in)                   :: tol
!-----------------------------------------------------------------------   
    real(8), dimension(:), intent(inout)  :: x
!-----------------------------------------------------------------------
    real(8), dimension(mat_size)          :: x_old, A_p, r, r_new, p, precond
    integer                               :: i, t
    real(8)                               :: alpha, betta, p_A_p, r_new_r_new, r_r, res
!-----------------------------------------------------------------------

    ! initial guess of x
    x = b

    ! calculate r = b-A*x;
    r = b - matmul(A, x)
    
    do i = 1, mat_size
        precond(i) = 1.0 / A(i,i)
    end do       

    p = precond * r

    do t = 1, cgitr

        x_old = x;

        !calculate r'*r
        r_r = sum(r * r * precond)

        ! calculate p'*A*p
        A_p = matmul(A, p)
        p_A_p = sum(p * A_p)

        alpha=0
        if (p_A_p == 0) then
            p_A_p = 1;
            write(*,'(a)') "WARNING: MATRIX SINGULARITY ISSUE"
        else
            alpha = r_r / p_A_p
        end if

        ! calculate x = x + alpha*p;
        x = x + alpha * p
        r_new = r - alpha * A_p
        !if (r_new[1] <= 0.00001)
        !{ break;}

        ! calculate r_new'*r_new
        r_new_r_new = sum(r_new * r_new * precond)

        betta = 0
        if (r_r <= 0.0) then
            r_r = 1
            write(*,'(a)') "ERROR: MATRIX SINGULARITY ISSUE"
        else
            betta = r_new_r_new / r_r
        end if
        
        p = precond * r_new + betta * p
        r = r_new
        res = sum(abs(x_old - x))
        if (res <= tol) then
            exit
        end if
    end do
    if (t < cgitr) then
        write(*, '(a, i0, a)') "Converged matrix solution obtained after ", t, " itrations..."
    else
        write(*,'(a)') "ERROR: No convergence within the specified tolerance was obtained for the PCG slvr"
        call exit(0)
    end if
end subroutine slv_pcg

end module eqslvrs
