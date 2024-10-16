module m_cuda_kernels_fieldops
  use cudafor

  use m_common, only: dp

contains

  attributes(global) subroutine axpby(n, alpha, x, beta, y)
    implicit none

    integer, value, intent(in) :: n
    real(dp), value, intent(in) :: alpha, beta
    real(dp), device, intent(in), dimension(:, :, :) :: x
    real(dp), device, intent(inout), dimension(:, :, :) :: y

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      y(i, j, b) = alpha*x(i, j, b) + beta*y(i, j, b)
    end do

  end subroutine axpby

  attributes(global) subroutine buffer_copy(u_send_s, u_send_e, u, n, n_halo)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: u_send_s, u_send_e
    real(dp), device, intent(in), dimension(:, :, :) :: u
    integer, value, intent(in) :: n, n_halo

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n_halo
      u_send_s(i, j, b) = u(i, j, b)
      u_send_e(i, j, b) = u(i, n - n_halo + j, b)
    end do

  end subroutine buffer_copy

  attributes(global) subroutine scalar_product(s, x, y, n)
    implicit none

    real(dp), device, intent(inout) :: s
    real(dp), device, intent(in), dimension(:, :, :) :: x, y
    integer, value, intent(in) :: n

    real(dp) :: s_pncl !! pencil sum
    integer :: i, j, b, ierr

    i = threadIdx%x
    b = blockIdx%x

    s_pncl = 0._dp
    do j = 1, n
      s_pncl = s_pncl + x(i, j, b)*y(i, j, b)
    end do
    ierr = atomicadd(s, s_pncl)

  end subroutine scalar_product

end module m_cuda_kernels_fieldops
