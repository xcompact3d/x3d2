module m_cuda_kernels_fieldops
  !! CUDA kernels for field operations (copy, scale, vector arithmetic, reductions).
  !!
  !! Provides GPU kernels for basic field manipulation: copying, scaling, shifting,
  !! linear combinations (AXPBY), pointwise multiplication, scalar products, and
  !! reductions (max, sum, volume integral). All kernels use thread-per-pencil-point
  !! parallelisation with [[m_cuda_common(module):SZ(variable)]] threads per block.
  use cudafor

  use m_common, only: dp
  use m_cuda_common, only: SZ

contains

  attributes(global) subroutine copy(n, dst, src)
    !! Copy field data: dst = src.
    implicit none

    integer, value, intent(in) :: n  !! Pencil length
    real(dp), device, intent(out), dimension(:, :, :) :: dst  !! Destination array
    real(dp), device, intent(in), dimension(:, :, :) :: src  !! Source array

    integer :: i  !! Thread index (pencil point)
    integer :: j  !! Pencil coordinate
    integer :: b  !! Block index (pencil number)

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      dst(i, j, b) = src(i, j, b)
    end do

  end subroutine copy

  attributes(global) subroutine axpby(n, alpha, x, beta, y)
    !! Compute linear combination: y = alpha*x + beta*y.
    implicit none

    integer, value, intent(in) :: n  !! Pencil length
    real(dp), value, intent(in) :: alpha, beta  !! Scalar coefficients
    real(dp), device, intent(in), dimension(:, :, :) :: x  !! Input array
    real(dp), device, intent(inout), dimension(:, :, :) :: y  !! Input/Output array

    integer :: i  !! Thread index (pencil point)
    integer :: j  !! Pencil coordinate
    integer :: b  !! Block index (pencil number)

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      y(i, j, b) = alpha*x(i, j, b) + beta*y(i, j, b)
    end do

  end subroutine axpby

  attributes(global) subroutine pwmul(y, x, n)
    !! Pointwise multiplication: y = y * x.
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: y  !! Input/Output array
    real(dp), device, intent(in), dimension(:, :, :) :: x  !! Multiplier array
    integer, value, intent(in) :: n  !! Pencil length

    integer :: i  !! Thread index (pencil point)
    integer :: j  !! Pencil coordinate
    integer :: b  !! Block index (pencil number)

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      y(i, j, b) = y(i, j, b)*x(i, j, b)
    end do

  end subroutine pwmul

  attributes(global) subroutine buffer_copy(u_send_s, u_send_e, u, n, n_halo)
    !! Copy halo regions into send buffers.
    !!
    !! Extracts first and last n_halo planes into separate buffers for MPI communication.
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: u_send_s  !! Start buffer
    real(dp), device, intent(inout), dimension(:, :, :) :: u_send_e  !! End buffer
    real(dp), device, intent(in), dimension(:, :, :) :: u  !! Source field
    integer, value, intent(in) :: n  !! Pencil length
    integer, value, intent(in) :: n_halo  !! Halo width

    integer :: i  !! Thread index (pencil point)
    integer :: j  !! Halo plane index
    integer :: b  !! Block index (pencil number)

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n_halo
      u_send_s(i, j, b) = u(i, j, b)
      u_send_e(i, j, b) = u(i, n - n_halo + j, b)
    end do

  end subroutine buffer_copy

  attributes(global) subroutine field_scale(f, alpha, n)
    !! Scale field by constant: f = alpha * f.
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: f  !! Field to scale
    real(dp), value, intent(in) :: alpha  !! Scaling factor
    integer, value, intent(in) :: n  !! Pencil length

    integer :: i  !! Thread index (pencil point)
    integer :: j  !! Pencil coordinate
    integer :: b  !! Block index (pencil number)

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      f(i, j, b) = alpha*f(i, j, b)
    end do

  end subroutine field_scale

  attributes(global) subroutine field_shift(f, const, n)
    !! Shift field by constant: f = f + const.
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: f  !! Field to shift
    real(dp), value, intent(in) :: const  !! Shift constant
    integer, value, intent(in) :: n  !! Pencil length

    integer :: i  !! Thread index (pencil point)
    integer :: j  !! Pencil coordinate
    integer :: b  !! Block index (pencil number)

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      f(i, j, b) = f(i, j, b) + const
    end do

  end subroutine field_shift

  attributes(global) subroutine scalar_product(s, x, y, n, n_i_pad, n_j)
    !! Compute scalar product with atomic reduction: s += sum(x * y).
    !!
    !! Uses atomic addition to accumulate partial sums from each pencil.
    implicit none

    real(dp), device, intent(inout) :: s  !! Accumulated scalar product
    real(dp), device, intent(in), dimension(:, :, :) :: x  !! First field
    real(dp), device, intent(in), dimension(:, :, :) :: y  !! Second field
    integer, value, intent(in) :: n  !! Pencil length
    integer, value, intent(in) :: n_i_pad  !! Padded dimension for indexing
    integer, value, intent(in) :: n_j  !! Active pencil count

    real(dp) :: s_pncl  !! Pencil sum
    integer :: i  !! Thread index
    integer :: j  !! Pencil coordinate
    integer :: b  !! Block index (pencil number)
    integer :: b_i, b_j  !! 2D block indices
    integer :: ierr  !! Atomic operation status

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
    !! Compute field maximum and sum with atomic reductions.
    !!
    !! Uses atomic max and add operations to accumulate pencil-wise results.
    implicit none

    real(dp), device, intent(inout) :: max_f  !! Accumulated maximum
    real(dp), device, intent(inout) :: sum_f  !! Accumulated sum
    real(dp), device, intent(in), dimension(:, :, :) :: f  !! Input field
    integer, value, intent(in) :: n  !! Pencil length
    integer, value, intent(in) :: n_i_pad  !! Padded dimension for indexing
    integer, value, intent(in) :: n_j  !! Active pencil count

    real(dp) :: max_pncl  !! Pencil maximum
    real(dp) :: sum_pncl  !! Pencil sum
    real(dp) :: val  !! Absolute value
    integer :: i  !! Thread index
    integer :: j  !! Pencil coordinate
    integer :: b  !! Block index (pencil number)
    integer :: b_i, b_j  !! 2D block indices
    integer :: ierr  !! Atomic operation status

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

  attributes(global) subroutine field_set_y_face(f, c_start, c_end, nx, ny, nz)
    !! Set Y-face boundary values to constants.
    !!
    !! Sets bottom face (y=0) to c_start and top face (y=L) to c_end.
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: f  !! Field to modify
    real(dp), value, intent(in) :: c_start  !! Bottom boundary value
    real(dp), value, intent(in) :: c_end  !! Top boundary value
    integer, value, intent(in) :: nx, ny, nz  !! Grid dimensions

    integer :: i  !! Thread index
    integer :: j  !! X-coordinate
    integer :: b  !! Z-coordinate block
    integer :: n_mod  !! Modulo for top boundary indexing
    integer :: b_end  !! Top boundary block index

    j = threadIdx%x + (blockIdx%x - 1)*blockDim%x ! from 1 to nx
    b = blockIdx%y ! from 1 to nz

    n_mod = mod(ny - 1, SZ) + 1
    b_end = b + (ny - 1)/SZ*nz

    if (j <= nx) then
      f(1, j, b) = c_start
      f(n_mod, j, b_end) = c_end
    end if

  end subroutine field_set_y_face

  attributes(global) subroutine volume_integral(s, f, n, n_i_pad, n_j)
    !! Compute volume integral with atomic reduction: s += sum(f).
    !!
    !! Uses atomic addition to accumulate partial sums from each pencil.
    implicit none

    real(dp), device, intent(inout) :: s  !! Accumulated integral
    real(dp), device, intent(in), dimension(:, :, :) :: f  !! Input field
    integer, value, intent(in) :: n  !! Pencil length
    integer, value, intent(in) :: n_i_pad  !! Padded dimension for indexing
    integer, value, intent(in) :: n_j  !! Active pencil count

    real(dp) :: s_pncl  !! Pencil sum
    integer :: i  !! Thread index
    integer :: j  !! Pencil coordinate
    integer :: b  !! Block index (pencil number)
    integer :: b_i, b_j  !! 2D block indices
    integer :: ierr  !! Atomic operation status

    i = threadIdx%x
    b_i = blockIdx%x
    b_j = blockIdx%y

    b = b_i + (b_j - 1)*n_i_pad
    s_pncl = 0._dp
    if (i + (b_j - 1)*blockDim%x <= n_j) then
      do j = 1, n
        s_pncl = s_pncl + f(i, j, b)
      end do
    end if
    ierr = atomicadd(s, s_pncl)

  end subroutine volume_integral

end module m_cuda_kernels_fieldops
