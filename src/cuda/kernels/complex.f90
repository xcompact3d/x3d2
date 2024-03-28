module m_cuda_complex
  use cudafor

  use m_common, only: dp
  use m_cuda_common, only: SZ

  implicit none

contains

  attributes(global) subroutine process_spectral_div_u( &
    div, waves, nx, ny, nz, ax, bx, ay, by, az, bz &
    )
    implicit none

    ! Arguments
    complex(dp), device, intent(inout), dimension(:, :, :) :: div
    complex(dp), device, intent(in), dimension(:, :, :) :: waves
    real(dp), device, intent(in), dimension(:) :: ax, bx, ay, by, az, bz
    integer, value, intent(in) :: nx, ny, nz

    ! Local variables
    integer :: i, j, k
    real(dp) :: tmp_r, tmp_c, div_r, div_c

    i = threadIdx%x
    k = blockIdx%x

    do j = 1, ny
      ! normalisation
      div_r = real(div(k, j, i), kind=dp)/(nx*ny*nz)
      div_c = aimag(div(k, j, i))/(nx*ny*nz)

      ! post-process forward
      ! post-process in z
      tmp_r = div_r
      tmp_c = div_c
      div_r = tmp_r*bz(k) + tmp_c*az(k)
      div_c = tmp_c*bz(k) - tmp_r*az(k)

      ! post-process in y
      tmp_r = div_r
      tmp_c = div_c
      div_r = tmp_r*by(j) + tmp_c*ay(j)
      div_c = tmp_c*by(j) - tmp_r*ay(j)
      if (j > ny/2 + 1) div_r = -div_r
      if (j > ny/2 + 1) div_c = -div_c

      ! post-process in x
      tmp_r = div_r
      tmp_c = div_c
      div_r = tmp_r*bx(i) + tmp_c*ax(i)
      div_c = tmp_c*bx(i) - tmp_r*ax(i)
      if (i > nx/2 + 1) div_r = -div_r
      if (i > nx/2 + 1) div_c = -div_c

      ! Solve Poisson
      tmp_r = real(waves(k, j, i), kind=dp)
      tmp_c = aimag(waves(k, j, i))
      if ((tmp_r < 1.e-16_dp) .or. (tmp_c < 1.e-16_dp)) then
        div_r = 0._dp; div_c = 0._dp
      else
        div_r = -div_r/tmp_r
        div_c = -div_c/tmp_c
      end if

      ! post-process backward
      ! post-process in z
      tmp_r = div_r
      tmp_c = div_c
      div_r = tmp_r*bz(k) - tmp_c*az(k)
      div_c = -tmp_c*bz(k) - tmp_r*az(k)

      ! post-process in y
      tmp_r = div_r
      tmp_c = div_c
      div_r = tmp_r*by(j) + tmp_c*ay(j)
      div_c = tmp_c*by(j) - tmp_r*ay(j)
      if (j > ny/2 + 1) div_r = -div_r
      if (j > ny/2 + 1) div_c = -div_c

      ! post-process in x
      tmp_r = div_r
      tmp_c = div_c
      div_r = tmp_r*bx(i) + tmp_c*ax(i)
      div_c = -tmp_c*bx(i) + tmp_r*ax(i)
      if (i > nx/2 + 1) div_r = -div_r
      if (i > nx/2 + 1) div_c = -div_c

      ! update the entry
      div(k, j, i) = cmplx(div_r, div_c, kind=dp)
    end do

  end subroutine process_spectral_div_u

  attributes(global) subroutine reorder_cmplx_x2y_T(u_y, u_x, nz)
    implicit none

    complex(dp), device, intent(out), dimension(:, :, :) :: u_y
    complex(dp), device, intent(in), dimension(:, :, :) :: u_x
    integer, value, intent(in) :: nz

    complex(dp), shared :: tile(SZ, SZ)

    integer :: i, j, b_i, b_j, b_k

    i = threadIdx%x; j = threadIdx%y
    b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

    ! copy into shared
    tile(i, j) = u_x((b_i - 1)*SZ + j, i, b_k + nz*(b_j - 1))

    call syncthreads()

    ! copy into output array from shared
    u_y((b_j - 1)*SZ + j, i, (b_i - 1)*nz + b_k) = tile(j, i)

  end subroutine reorder_cmplx_x2y_T

  attributes(global) subroutine reorder_cmplx_y2x_T(u_x, u_y, nz)
    implicit none

    complex(dp), device, intent(out), dimension(:, :, :) :: u_x
    complex(dp), device, intent(in), dimension(:, :, :) :: u_y
    integer, value, intent(in) :: nz

    complex(dp), shared :: tile(SZ, SZ)

    integer :: i, j, b_i, b_j, b_k

    i = threadIdx%x; j = threadIdx%y
    b_i = blockIdx%x; b_j = blockIdx%y; b_k = blockIdx%z

    ! copy into shared
    tile(i, j) = u_y((b_j - 1)*SZ + j, i, b_k + nz*(b_i - 1))

    call syncthreads()

    ! copy into output array from shared
    u_x((b_i - 1)*SZ + j, i, (b_j - 1)*nz + b_k) = tile(j, i)

  end subroutine reorder_cmplx_y2x_T

  attributes(global) subroutine reorder_cmplx_y2z_T(u_z, u_y, nx, nz)
    implicit none

    complex(dp), device, intent(out), dimension(:, :, :) :: u_z
    complex(dp), device, intent(in), dimension(:, :, :) :: u_y
    integer, value, intent(in) :: nx, nz

    complex(dp), shared :: tile(SZ, SZ)

    integer :: i, j, k, b_i, b_j, b_k, b_x, b_y, b_z

    i = threadIdx%x
    j = threadIdx%y
    k = threadIdx%z

    b_x = blockIdx%z
    b_y = blockIdx%y
    b_z = blockIdx%x

    ! copy into shared
    if (j + (b_z - 1)*SZ <= nz) &
      tile(i, j) = u_y(i + (b_y - 1)*SZ, mod(b_x - 1, SZ) + 1, &
                       j + (b_z - 1)*SZ + ((b_x - 1)/SZ)*nz)

    call syncthreads()

    ! copy into output array from shared
    if (i + (b_z - 1)*SZ <= nz) &
      u_z(i + (b_z - 1)*SZ, j, b_x + (b_y - 1)*nx) = tile(j, i)

  end subroutine reorder_cmplx_y2z_T

  attributes(global) subroutine reorder_cmplx_z2y_T(u_y, u_z, nx, nz)
    implicit none

    complex(dp), device, intent(out), dimension(:, :, :) :: u_y
    complex(dp), device, intent(in), dimension(:, :, :) :: u_z
    integer, value, intent(in) :: nx, nz

    complex(dp), shared :: tile(SZ, SZ)

    integer :: i, j, k, b_x, b_y, b_z

    i = threadIdx%x
    j = threadIdx%y
    k = threadIdx%z

    b_x = blockIdx%z
    b_y = blockIdx%y
    b_z = blockIdx%x

    ! copy into shared
    if (i + (b_z - 1)*SZ <= nz) &
      tile(i, j) = u_z(i + (b_z - 1)*SZ, j, b_x + (b_y - 1)*nx)

    call syncthreads()

    ! copy into output array from shared
    if (j + (b_z - 1)*SZ <= nz) &
      u_y(i + (b_y - 1)*SZ, mod(b_x - 1, SZ) + 1, &
          j + (b_z - 1)*SZ + ((b_x - 1)/SZ)*nz) = tile(j, i)

  end subroutine reorder_cmplx_z2y_T

  attributes(global) subroutine reshapeDSF(uout, uin)
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: uout
    real(dp), device, intent(in), dimension(:, :, :) :: uin

    real(dp), shared :: tile(SZ + 1, SZ)

    integer :: i, j, b_i, b

    i = threadIdx%x; j = threadIdx%y
    b_i = blockIdx%x; b = blockIdx%y

    tile(i, j) = uin(i, j + (b_i - 1)*SZ, b)

    call syncthreads()

    uout(i + (b_i - 1)*SZ, j, b) = tile(j, i)

  end subroutine reshapeDSF

  attributes(global) subroutine reshapeDSB(uout, uin)
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: uout
    real(dp), device, intent(in), dimension(:, :, :) :: uin

    real(dp), shared :: tile(SZ + 1, SZ)

    integer :: i, j, b_i, b

    i = threadIdx%x; j = threadIdx%y
    b_i = blockIdx%x; b = blockIdx%y

    tile(i, j) = uin(i + (b_i - 1)*SZ, j, b)

    call syncthreads()

    uout(i, j + (b_i - 1)*SZ, b) = tile(j, i)

  end subroutine reshapeDSB

  attributes(global) subroutine reshapeCDSF(uout, uin)
    implicit none

    complex(dp), device, intent(out), dimension(:, :, :) :: uout
    complex(dp), device, intent(in), dimension(:, :, :) :: uin

    complex(dp), shared :: tile(SZ + 1, SZ)

    integer :: i, j, b_i, b

    i = threadIdx%x; j = threadIdx%y
    b_i = blockIdx%x; b = blockIdx%y

    tile(i, j) = uin(i, j + (b_i - 1)*SZ, b)

    call syncthreads()

    uout(i + (b_i - 1)*SZ, j, b) = tile(j, i)

  end subroutine reshapeCDSF

  attributes(global) subroutine reshapeCDSB(uout, uin)
    implicit none

    complex(dp), device, intent(out), dimension(:, :, :) :: uout
    complex(dp), device, intent(in), dimension(:, :, :) :: uin

    complex(dp), shared :: tile(SZ + 1, SZ)

    integer :: i, j, b_i, b

    i = threadIdx%x; j = threadIdx%y
    b_i = blockIdx%x; b = blockIdx%y

    tile(i, j) = uin(i + (b_i - 1)*SZ, j, b)

    call syncthreads()

    uout(i, j + (b_i - 1)*SZ, b) = tile(j, i)

  end subroutine reshapeCDSB

end module m_cuda_complex
