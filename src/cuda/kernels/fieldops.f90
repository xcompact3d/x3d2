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

  attributes(global) subroutine field_scale(f, alpha, n)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: f
    real(dp), value, intent(in) :: alpha
    integer, value, intent(in) :: n

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      f(i, j, b) = alpha*f(i, j, b)
    end do

  end subroutine field_scale

  attributes(global) subroutine field_shift(f, const, n)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: f
    real(dp), value, intent(in) :: const
    integer, value, intent(in) :: n

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      f(i, j, b) = f(i, j, b) + const
    end do

  end subroutine field_shift

  attributes(global) subroutine scalar_product(s, x, y, n, n_i_pad, n_j)
    implicit none

    real(dp), device, intent(inout) :: s
    real(dp), device, intent(in), dimension(:, :, :) :: x, y
    integer, value, intent(in) :: n, n_i_pad, n_j

    real(dp) :: s_pncl !! pencil sum
    integer :: i, j, b, b_i, b_j, ierr

    i = threadIdx%x
    b_i = blockIdx%x
    b_j = blockIdx%y

    b = b_i + (b_j - 1)*n_i_pad
    s_pncl = 0._dp
    if (i + (b_j - 1)*blockDim%x <= n_j) then
      do j = 1, n
        s_pncl = s_pncl + x(i, j, b)*y(i, j, b)
      end do
    end if
    ierr = atomicadd(s, s_pncl)

  end subroutine scalar_product

  attributes(global) subroutine field_max_sum(max_f, sum_f, f, n, n_i_pad, n_j)
    implicit none

    real(dp), device, intent(inout) :: max_f, sum_f
    real(dp), device, intent(in), dimension(:, :, :) :: f
    integer, value, intent(in) :: n, n_i_pad, n_j

    real(dp) :: max_pncl, sum_pncl, val
    integer :: i, j, b, b_i, b_j, ierr

    i = threadIdx%x
    b_i = blockIdx%x
    b_j = blockIdx%y

    b = b_i + (b_j - 1)*n_i_pad
    max_pncl = 0._dp
    sum_pncl = 0._dp
    if (i + (b_j - 1)*blockDim%x <= n_j) then
      do j = 1, n
        val = abs(f(i, j, b))
        sum_pncl = sum_pncl + val
        max_pncl = max(max_pncl, val)
      end do
    end if
    ierr = atomicadd(sum_f, sum_pncl)
    ierr = atomicmax(max_f, max_pncl)

  end subroutine field_max_sum

end module m_cuda_kernels_fieldops
