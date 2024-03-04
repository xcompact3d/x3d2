module m_cuda_kernels_reorder
   use cudafor

   use m_common, only: dp
   use m_cuda_common, only: SZ

contains

   attributes(global) subroutine reorder_x2y(u_y, u_x, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_y
      real(dp), device, intent(in), dimension(:, :, :) :: u_x
      integer, value, intent(in) :: nz

      real(dp), shared :: tile(SZ, SZ)
      integer :: i, j, b_i, b_j, b_k

      i = threadIdx%x; j = threadIdx%y
      b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

      ! copy into shared
      tile(i, j) = u_x(i, j + (b_i - 1)*SZ, b_j + (b_k - 1)*nz)

      call syncthreads()

      ! copy into output array from shared
      u_y(i, j + (b_k - 1)*SZ, b_j + (b_i - 1)*nz) = tile(j, i)

   end subroutine reorder_x2y

   attributes(global) subroutine reorder_x2z(u_z, u_x, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_z
      real(dp), device, intent(in), dimension(:, :, :) :: u_x
      integer, value, intent(in) :: nz

      integer :: i, j, b_i, b_j, nx

      i = threadIdx%x; b_i = blockIdx%x; b_j = blockIdx%y
      nx = gridDim%x

      ! Data access pattern for reordering between x and z is quite nice
      ! thus we don't need to use shared memory for this operation.
      do j = 1, nz
         u_z(i, j, b_i + (b_j - 1)*nx) = u_x(i, b_i, j + (b_j - 1)*nz)
      end do

   end subroutine reorder_x2z

   attributes(global) subroutine reorder_y2x(u_x, u_y, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_x
      real(dp), device, intent(in), dimension(:, :, :) :: u_y
      integer, value, intent(in) :: nz

      real(dp), shared :: tile(SZ, SZ)
      integer :: i, j, b_i, b_j, b_k

      i = threadIdx%x; j = threadIdx%y
      b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

      ! copy into shared
      tile(i, j) = u_y(i, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k)

      call syncthreads()

      ! copy into output array from shared
      u_x(i, (b_i - 1)*SZ + j, (b_j - 1)*nz + b_k) = tile(j, i)

   end subroutine reorder_y2x

   attributes(global) subroutine reorder_y2z(u_z, u_y, nx, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_z
      real(dp), device, intent(in), dimension(:, :, :) :: u_y
      integer, value, intent(in) :: nx, nz

      real(dp), shared :: tile(SZ, SZ)
      integer :: i, j, b_i, b_j, b_k

      i = threadIdx%x; j = threadIdx%y
      b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

      ! copy into shared
      tile(i, j) = u_y(i, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k)

      call syncthreads()

      ! copy into output array from shared
      u_z(i, b_k, (b_i - 1)*SZ + j + (b_j - 1)*nx) = tile(j, i)

   end subroutine reorder_y2z

   attributes(global) subroutine reorder_z2x(u_x, u_z, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_x
      real(dp), device, intent(in), dimension(:, :, :) :: u_z
      integer, value, intent(in) :: nz

      integer :: i, j, b_i, b_j, nx

      i = threadIdx%x; b_i = blockIdx%x; b_j = blockIdx%y
      nx = gridDim%x

      do j = 1, nz
         u_x(i, b_i, j + (b_j - 1)*nz) = u_z(i, j, b_i + (b_j - 1)*nx)
      end do

   end subroutine reorder_z2x

   attributes(global) subroutine reorder_z2y(u_y, u_z, nx, nz)
      implicit none

      real(dp), device, intent(out), dimension(:, :, :) :: u_y
      real(dp), device, intent(in), dimension(:, :, :) :: u_z
      integer, value, intent(in) :: nx, nz

      real(dp), shared :: tile(SZ, SZ)
      integer :: i, j, b_i, b_j, b_k

      i = threadIdx%x; j = threadIdx%y
      b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

      ! copy into shared
      tile(i, j) = u_z(i, b_k, (b_i - 1)*SZ + j + (b_j - 1)*nx)

      call syncthreads()

      ! copy into output array from shared
      u_y(i, (b_j - 1)*SZ + j, (b_i - 1)*nz + b_k) = tile(j, i)

   end subroutine reorder_z2y

   attributes(global) subroutine sum_yintox(u_x, u_y, nz)
      implicit none

      real(dp), device, intent(inout), dimension(:, :, :) :: u_x
      real(dp), device, intent(in), dimension(:, :, :) :: u_y
      integer, value, intent(in) :: nz

      real(dp), shared :: tile(SZ,SZ)
      integer :: i, j, b_i, b_j, b_k

      i = threadIdx%x; j = threadIdx%y
      b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

      ! copy into shared
      tile(i, j) = u_y(i, (b_j-1)*SZ+j, (b_k)+nz*(b_i-1))

      call syncthreads()

      ! copy into output array from shared
      u_x(i, (b_i-1)*SZ+j, (b_j-1)*nz+(b_k)) = &
         u_x(i, (b_i-1)*SZ+j, (b_j-1)*nz+(b_k)) + tile(j, i)

   end subroutine sum_yintox

   attributes(global) subroutine sum_zintox(u_x, u_z, nz)
      implicit none

      ! Arguments
      real(dp), device, intent(inout), dimension(:, :, :) :: u_x
      real(dp), device, intent(in), dimension(:, :, :) :: u_z
      integer, value, intent(in) :: nz

      integer :: i, j, b_i, b_j, nx

      i = threadIdx%x; b_i = blockIdx%x; b_j = blockIdx%y
      nx = gridDim%x

      do j = 1, nz
         u_x(i, b_i, j+(b_j-1)*nz) = u_x(i, b_i, j+(b_j-1)*nz) &
                                     + u_z(i, j, b_i+(b_j-1)*nx)
      end do

   end subroutine sum_zintox

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

end module m_cuda_kernels_reorder
