module m_cuda_kernels_reorder
  !! CUDA kernels for pencil reordering and accumulation between X/Y/Z orientations.
  !!
  !! Provides GPU kernels for rearranging field data between different pencil decompositions
  !! (X-pencils, Y-pencils, Z-pencils, and Cartesian). Most kernels use shared memory tiles
  !! for coalesced memory access. Thread blocks use [[m_cuda_common(module):SZ(variable)]]
  !! configuration (32x1 or 32x32 depending on operation).
  use cudafor

  use m_common, only: dp
  use m_cuda_common, only: SZ

contains

  attributes(global) subroutine reorder_c2x(u_x, u_c, nz)
    !! Reorder from Cartesian to X-pencil orientation.
    !!
    !! Uses shared memory transpose for efficient reordering.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: u_x  !! Output: X-pencil data
    real(dp), device, intent(in), dimension(:, :, :) :: u_c  !! Input: Cartesian data
    integer, value, intent(in) :: nz  !! Z-dimension size

    real(dp), shared :: tile(SZ, SZ)  !! Shared memory for transpose
    integer :: i, j  !! Thread indices
    integer :: b_i, b_j, b_k  !! Block indices

    i = threadIdx%x; j = threadIdx%y; 
    b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

    ! copy into shared
    tile(i, j) = u_c(i + (b_i - 1)*SZ, j + (b_j - 1)*SZ, b_k)
    if (SZ == 64) then
      tile(i + 32, j) = u_c(i + 32 + (b_i - 1)*SZ, j + (b_j - 1)*SZ, b_k)
      tile(i, j + 32) = u_c(i + (b_i - 1)*SZ, j + 32 + (b_j - 1)*SZ, b_k)
      tile(i + 32, j + 32) = &
        u_c(i + 32 + (b_i - 1)*SZ, j + 32 + (b_j - 1)*SZ, b_k)
    end if

    call syncthreads()

    ! copy into output array from shared
    u_x(i, j + (b_i - 1)*SZ, b_k + (b_j - 1)*nz) = tile(j, i)
    if (SZ == 64) then
      u_x(i + 32, j + (b_i - 1)*SZ, b_k + (b_j - 1)*nz) = tile(j, i + 32)
      u_x(i, j + 32 + (b_i - 1)*SZ, b_k + (b_j - 1)*nz) = tile(j + 32, i)
      u_x(i + 32, j + 32 + (b_i - 1)*SZ, b_k + (b_j - 1)*nz) = &
        tile(j + 32, i + 32)
    end if

  end subroutine reorder_c2x

  attributes(global) subroutine reorder_x2c(u_c, u_x, nz)
    !! Reorder from X-pencil to Cartesian orientation.
    !!
    !! Inverse of reorder_c2x. Uses shared memory transpose.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: u_c  !! Output: Cartesian data
    real(dp), device, intent(in), dimension(:, :, :) :: u_x  !! Input: X-pencil data
    integer, value, intent(in) :: nz  !! Z-dimension size

    real(dp), shared :: tile(SZ, SZ)  !! Shared memory for transpose
    integer :: i, j  !! Thread indices
    integer :: b_i, b_j, b_k  !! Block indices

    i = threadIdx%x; j = threadIdx%y; 
    b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

    ! copy into shared
    tile(i, j) = u_x(i, j + (b_i - 1)*SZ, b_k + (b_j - 1)*nz)
    if (SZ == 64) then
      tile(i + 32, j) = u_x(i + 32, j + (b_i - 1)*SZ, b_k + (b_j - 1)*nz)
      tile(i, j + 32) = u_x(i, j + 32 + (b_i - 1)*SZ, b_k + (b_j - 1)*nz)
      tile(i + 32, j + 32) = &
        u_x(i + 32, j + 32 + (b_i - 1)*SZ, b_k + (b_j - 1)*nz)
    end if

    call syncthreads()

    ! copy into output array from shared
    u_c(i + (b_i - 1)*SZ, j + (b_j - 1)*SZ, b_k) = tile(j, i)
    if (SZ == 64) then
      u_c(i + 32 + (b_i - 1)*SZ, j + (b_j - 1)*SZ, b_k) = tile(j, i + 32)
      u_c(i + (b_i - 1)*SZ, j + 32 + (b_j - 1)*SZ, b_k) = tile(j + 32, i)
      u_c(i + 32 + (b_i - 1)*SZ, j + 32 + (b_j - 1)*SZ, b_k) = &
        tile(j + 32, i + 32)
    end if

  end subroutine reorder_x2c

  attributes(global) subroutine reorder_x2y(u_y, u_x, nz)
    !! Reorder from X-pencil to Y-pencil orientation.
    !!
    !! Uses shared memory transpose for efficient reordering.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: u_y  !! Output: Y-pencil data
    real(dp), device, intent(in), dimension(:, :, :) :: u_x  !! Input: X-pencil data
    integer, value, intent(in) :: nz  !! Z-dimension size

    real(dp), shared :: tile(SZ, SZ)  !! Shared memory for transpose
    integer :: i, j  !! Thread indices
    integer :: b_i, b_j, b_k  !! Block indices

    i = threadIdx%x; j = threadIdx%y; 
    b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

    ! copy into shared
    tile(i, j) = u_x(i, j + (b_i - 1)*SZ, b_j + (b_k - 1)*nz)
    if (SZ == 64) then
      tile(i + 32, j) = u_x(i + 32, j + (b_i - 1)*SZ, b_j + (b_k - 1)*nz)
      tile(i, j + 32) = u_x(i, j + 32 + (b_i - 1)*SZ, b_j + (b_k - 1)*nz)
      tile(i + 32, j + 32) = &
        u_x(i + 32, j + 32 + (b_i - 1)*SZ, b_j + (b_k - 1)*nz)
    end if

    call syncthreads()

    ! copy into output array from shared
    u_y(i, j + (b_k - 1)*SZ, b_j + (b_i - 1)*nz) = tile(j, i)
    if (SZ == 64) then
      u_y(i + 32, j + (b_k - 1)*SZ, b_j + (b_i - 1)*nz) = tile(j, i + 32)
      u_y(i, j + 32 + (b_k - 1)*SZ, b_j + (b_i - 1)*nz) = tile(j + 32, i)
      u_y(i + 32, j + 32 + (b_k - 1)*SZ, b_j + (b_i - 1)*nz) = &
        tile(j + 32, i + 32)
    end if

  end subroutine reorder_x2y

  attributes(global) subroutine reorder_x2z(u_z, u_x, nz)
    !! Reorder from X-pencil to Z-pencil orientation.
    !!
    !! No shared memory needed - memory access pattern is already favourable.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: u_z  !! Output: Z-pencil data
    real(dp), device, intent(in), dimension(:, :, :) :: u_x  !! Input: X-pencil data
    integer, value, intent(in) :: nz  !! Z-dimension size

    integer :: i  !! Thread index
    integer :: j  !! Loop index
    integer :: b_i, b_j  !! Block indices
    integer :: nx  !! Grid X-dimension

    i = threadIdx%x; b_i = blockIdx%x; b_j = blockIdx%y
    nx = gridDim%x

    ! Data access pattern for reordering between x and z is quite nice
    ! thus we don't need to use shared memory for this operation.
    do j = 1, nz
      u_z(i, j, b_i + (b_j - 1)*nx) = u_x(i, b_i, j + (b_j - 1)*nz)
    end do

  end subroutine reorder_x2z

  attributes(global) subroutine reorder_y2x(u_x, u_y, nz)
    !! Reorder from Y-pencil to X-pencil orientation.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: u_x  !! Output: X-pencil data
    real(dp), device, intent(in), dimension(:, :, :) :: u_y  !! Input: Y-pencil data
    integer, value, intent(in) :: nz  !! Z-dimension size

    real(dp), shared :: tile(SZ, SZ)  !! Shared memory for transpose
    integer :: i, j, b_i, b_j, b_k  !! Thread and block indices

    i = threadIdx%x; j = threadIdx%y; 
    b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

    ! copy into shared
    tile(i, j) = u_y(i, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k)
    if (SZ == 64) then
      tile(i + 32, j) = u_y(i + 32, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k)
      tile(i, j + 32) = u_y(i, (b_j - 1)*SZ + j + 32, (b_i - 1)*nz + b_k)
      tile(i + 32, j + 32) = &
        u_y(i + 32, (b_j - 1)*SZ + j + 32, (b_i - 1)*nz + b_k)
    end if

    call syncthreads()

    ! copy into output array from shared
    u_x(i, (b_i - 1)*SZ + j, (b_j - 1)*nz + b_k) = tile(j, i)
    if (SZ == 64) then
      u_x(i + 32, (b_i - 1)*SZ + j, (b_j - 1)*nz + b_k) = tile(j, i + 32)
      u_x(i, (b_i - 1)*SZ + j + 32, (b_j - 1)*nz + b_k) = tile(j + 32, i)
      u_x(i + 32, (b_i - 1)*SZ + j + 32, (b_j - 1)*nz + b_k) = &
        tile(j + 32, i + 32)
    end if

  end subroutine reorder_y2x

  attributes(global) subroutine reorder_y2z(u_z, u_y, nx, nz)
    !! Reorder from Y-pencil to Z-pencil orientation.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: u_z  !! Output: Z-pencil data
    real(dp), device, intent(in), dimension(:, :, :) :: u_y  !! Input: Y-pencil data
    integer, value, intent(in) :: nx, nz  !! Grid dimensions

    real(dp), shared :: tile(SZ, SZ)  !! Shared memory for transpose
    integer :: i, j, b_i, b_j, b_k  !! Thread and block indices

    i = threadIdx%x; j = threadIdx%y; 
    b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

    ! copy into shared
    tile(i, j) = u_y(i, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k)
    if (SZ == 64) then
      tile(i + 32, j) = u_y(i + 32, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k)
      tile(i, j + 32) = u_y(i, (b_j - 1)*SZ + j + 32, (b_i - 1)*nz + b_k)
      tile(i + 32, j + 32) = &
        u_y(i + 32, (b_j - 1)*SZ + j + 32, (b_i - 1)*nz + b_k)
    end if

    call syncthreads()

    ! copy into output array from shared
    u_z(i, b_k, (b_i - 1)*SZ + j + (b_j - 1)*nx) = tile(j, i)
    if (SZ == 64) then
      u_z(i + 32, b_k, (b_i - 1)*SZ + j + (b_j - 1)*nx) = tile(j, i + 32)
      u_z(i, b_k, (b_i - 1)*SZ + j + 32 + (b_j - 1)*nx) = tile(j + 32, i)
      u_z(i + 32, b_k, (b_i - 1)*SZ + j + 32 + (b_j - 1)*nx) = &
        tile(j + 32, i + 32)
    end if

  end subroutine reorder_y2z

  attributes(global) subroutine reorder_z2x(u_x, u_z, nz)
    !! Reorder from Z-pencil to X-pencil orientation.
    !!
    !! No shared memory needed - favourable memory access pattern.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: u_x  !! Output: X-pencil data
    real(dp), device, intent(in), dimension(:, :, :) :: u_z  !! Input: Z-pencil data
    integer, value, intent(in) :: nz  !! Z-dimension size

    integer :: i, j, b_i, b_j, nx  !! Thread, loop, block indices and grid size

    i = threadIdx%x; b_i = blockIdx%x; b_j = blockIdx%y
    nx = gridDim%x

    do j = 1, nz
      u_x(i, b_i, j + (b_j - 1)*nz) = u_z(i, j, b_i + (b_j - 1)*nx)
    end do

  end subroutine reorder_z2x

  attributes(global) subroutine reorder_z2y(u_y, u_z, nx, nz)
    !! Reorder from Z-pencil to Y-pencil orientation.
    !!
    !! Uses shared memory tile for coalesced access pattern.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: u_y  !! Output: Y-pencil data
    real(dp), device, intent(in), dimension(:, :, :) :: u_z  !! Input: Z-pencil data
    integer, value, intent(in) :: nx, nz  !! X and Z dimension sizes

    real(dp), shared :: tile(SZ, SZ)  !! Shared memory tile for transpose
    integer :: i, j, b_i, b_j, b_k  !! Thread, block indices

    i = threadIdx%x; j = threadIdx%y; 
    b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

    ! copy into shared
    tile(i, j) = u_z(i, b_k, (b_i - 1)*SZ + j + (b_j - 1)*nx)
    if (SZ == 64) then
      tile(i + 32, j) = u_z(i + 32, b_k, (b_i - 1)*SZ + j + (b_j - 1)*nx)
      tile(i, j + 32) = u_z(i, b_k, (b_i - 1)*SZ + j + 32 + (b_j - 1)*nx)
      tile(i + 32, j + 32) = &
        u_z(i + 32, b_k, (b_i - 1)*SZ + j + 32 + (b_j - 1)*nx)
    end if

    call syncthreads()

    ! copy into output array from shared
    u_y(i, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k) = tile(j, i)
    if (SZ == 64) then
      u_y(i + 32, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k) = tile(j, i + 32)
      u_y(i, (b_j - 1)*SZ + j + 32, (b_i - 1)*nz + b_k) = tile(j + 32, i)
      u_y(i + 32, (b_j - 1)*SZ + j + 32, (b_i - 1)*nz + b_k) = &
        tile(j + 32, i + 32)
    end if

  end subroutine reorder_z2y

  attributes(global) subroutine sum_yintox(u_x, u_y, nz)
    !! Accumulate Y-pencil contributions into X-pencil data.
    !!
    !! Performs u_x += u_y with reordering. Uses shared memory tile
    !! for efficient transpose and coalesced memory access.
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: u_x  !! In/out: X-pencil data to accumulate into
    real(dp), device, intent(in), dimension(:, :, :) :: u_y  !! Input: Y-pencil data to add
    integer, value, intent(in) :: nz  !! Z-dimension size

    real(dp), shared :: tile(SZ, SZ)  !! Shared memory tile for transpose
    integer :: i, j, b_i, b_j, b_k  !! Thread, block indices

    i = threadIdx%x; j = threadIdx%y; 
    b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

    ! copy into shared
    tile(i, j) = u_y(i, (b_j - 1)*SZ + j, (b_k) + nz*(b_i - 1))
    if (SZ == 64) then
      tile(i + 32, j) = u_y(i + 32, (b_j - 1)*SZ + j, (b_k) + nz*(b_i - 1))
      tile(i, j + 32) = u_y(i, (b_j - 1)*SZ + j + 32, (b_k) + nz*(b_i - 1))
      tile(i + 32, j + 32) = &
        u_y(i + 32, (b_j - 1)*SZ + j + 32, (b_k) + nz*(b_i - 1))
    end if

    call syncthreads()

    ! copy into output array from shared
    u_x(i, (b_i - 1)*SZ + j, (b_j - 1)*nz + (b_k)) = &
      u_x(i, (b_i - 1)*SZ + j, (b_j - 1)*nz + (b_k)) + tile(j, i)
    if (SZ == 64) then
      u_x(i + 32, (b_i - 1)*SZ + j, (b_j - 1)*nz + (b_k)) = &
        u_x(i + 32, (b_i - 1)*SZ + j, (b_j - 1)*nz + (b_k)) + tile(j, i + 32)
      u_x(i, (b_i - 1)*SZ + j + 32, (b_j - 1)*nz + (b_k)) = &
        u_x(i, (b_i - 1)*SZ + j + 32, (b_j - 1)*nz + (b_k)) + tile(j + 32, i)
      u_x(i + 32, (b_i - 1)*SZ + j + 32, (b_j - 1)*nz + (b_k)) = &
        u_x(i + 32, (b_i - 1)*SZ + j + 32, (b_j - 1)*nz + (b_k)) + &
        tile(j + 32, i + 32)
    end if

  end subroutine sum_yintox

  attributes(global) subroutine sum_zintox(u_x, u_z, nz)
    !! Accumulate Z-pencil contributions into X-pencil data.
    !!
    !! Performs u_x += u_z with reordering. No shared memory needed
    !! due to favourable memory access pattern.
    implicit none

    ! Arguments
    real(dp), device, intent(inout), dimension(:, :, :) :: u_x  !! In/out: X-pencil data to accumulate into
    real(dp), device, intent(in), dimension(:, :, :) :: u_z  !! Input: Z-pencil data to add
    integer, value, intent(in) :: nz  !! Z-dimension size

    integer :: i, j, b_i, b_j, nx  !! Thread, loop, block indices and grid size

    i = threadIdx%x; b_i = blockIdx%x; b_j = blockIdx%y
    nx = gridDim%x

    do j = 1, nz
      u_x(i, b_i, j + (b_j - 1)*nz) = u_x(i, b_i, j + (b_j - 1)*nz) &
                                      + u_z(i, j, b_i + (b_j - 1)*nx)
    end do

  end subroutine sum_zintox

end module m_cuda_kernels_reorder
