program main
    use matrix
    use lu
    implicit none
    integer :: N
    real(8), allocatable :: A(:, :), B(:, :), L(:, :), U(:, :)
    allocate(A(N, N), B(N, N), L(N, N), U(N, N))
    
    N = 3
    call make_random(A, N)
    call copy_matrix(A, B, N)
    call decompose_v1(L, U, A, N)
    call multiply(L, U, B, N)
    call print_matrix(A, N)
    call print_matrix(B, N)

    deallocate(A)
    deallocate(B)
    deallocate(L)
    deallocate(U)


end program main

! example.f90
! program main
    
!     implicit none (type, external)
!     external :: sscal

!     integer, parameter :: N = 3

!     real :: x(N)
!     real :: a

!     x = [ 5., 6., 7. ]
!     a = 5.

!     print '("a = ", f0.1)', a
!     print '("X = [ ", 3(f0.1, " "), "]")', x

!     call sscal(N, a, x, 1)

!     print '(/, "X = a * X")'
!     print '("X = [ ", 3(f0.1, " "), "]")', x

! end program main