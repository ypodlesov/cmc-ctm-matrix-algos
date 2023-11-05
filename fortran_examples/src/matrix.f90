module matrix
    implicit none
  
    private  ! All entities are now module-private by default
    public multiply  ! Explicitly export public entities
    public make_random  ! Explicitly export public entities
    public print_matrix
    public copy_matrix

    real, parameter :: public_var = 2
    integer :: private_var
  
contains
    subroutine multiply(A, B, C, N)
        integer, intent(in) :: N
        real(8), intent(in) :: A(N, N), B(N, N)
        real(8), intent(inout) :: C(N, N)

        integer :: i, j, k
        do i = 1, N
            do j = 1, N
                C(i, j) = 0
                do k = 1, N
                    C(i, j) = C(i, j) + A(i, k) * B(k, j)
                end do
            end do
        end do

    end subroutine multiply

    subroutine make_random(A, N)
        integer :: N
        real(8), intent(inout) :: A(N, N)

        ! call random_seed(size = N)
        
        integer :: i, j

        do i = 1, N
            do j = 1, N
                call random_number(A(i, j))
            end do
        end do
    end subroutine make_random

    subroutine print_matrix(A, N)
        
        implicit none
        integer :: N
        real(8), dimension(N, N) :: A
        integer i, j
        do i = 1, N
            do j = 1, N
                write (*, "(dc,f12.3)", advance="no") A(i,j)
            enddo
            print *
        enddo
        print *
    
    end subroutine print_matrix

    subroutine copy_matrix(from, to, N)

        implicit none
        integer, intent(in) :: N
        real(8), intent(inout) :: from(N, N), to(N, N)

        integer :: i, j
        do i = 1, N
            do j = 1, N
                to(i, j) = from(i, j)
            end do
        end do

    end subroutine copy_matrix
    

end module matrix
