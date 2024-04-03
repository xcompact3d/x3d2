module m_cuda_complex
  use cudafor

  use m_common, only: dp

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
      div_r = real(div(i, j, k), kind=dp)/(nx*ny*nz)
      div_c = aimag(div(i, j, k))/(nx*ny*nz)

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
      tmp_r = real(waves(i, j, k), kind=dp)
      tmp_c = aimag(waves(i, j, k))
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
      div(i, j, k) = cmplx(div_r, div_c, kind=dp)
    end do

  end subroutine process_spectral_div_u

end module m_cuda_complex
