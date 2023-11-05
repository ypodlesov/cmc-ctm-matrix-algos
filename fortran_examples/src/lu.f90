module lu

    public decompose
    public decompose_v1

    contains
    subroutine decompose(L, U, A, N)

        implicit none
        integer :: N
        real(8), dimension(N, N) :: L, U, A
        integer i, j, k

        do i = 1, N
            L(i, i) = 1
        end do

        do i = 1, N
            do j = i, N
                U(i, j) = A(i,j)
                do k = 1, i - 1
                    U(i, j) = U(i, j) - L(i, k) * U(k, j)
                end do
            end do
            do j = i + 1, N
                L(j, i) = A(j, i)
                do k = 1, i - 1
                    L(j, i) = L(j, i) - L(j, k) * U(k, i)
                end do
                L(j, i) = L(j, i) / U(i, i)
            end do
        end do

    end subroutine decompose

    subroutine decompose_v1(L, U, A, N)

        implicit none
        real(8), external :: ddot
        integer :: N
        real(8), dimension(N, N) :: L, U, A
        integer i, j

        do i = 1, N
            L(i, i) = 1
        end do

        do i = 1, N
            do j = i, N
                U(i, j) = A(i,j) - ddot(i - 1, L(i, 1 : i-1), 1, U(1 : i-1, j), 1) 
            end do
            do j = i + 1, N
                L(j, i) = A(j, i) - ddot(i - 1, L(j, 1 : i-1), 1, U(1 : i-1, i), 1)
                L(j, i) = L(j, i) / U(i, i)
            end do
        end do

    end subroutine decompose_v1

end module lu