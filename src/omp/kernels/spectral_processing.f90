module m_omp_spectral
  use m_common, only: dp
  implicit none

contains
  subroutine process_spectral_div_u( &
    div_u, waves, nx_spec, ny_spec, nz_spec, x_sp_st, y_sp_st, z_sp_st, &
    nx, ny, nz, ax, bx, ay, by, az, bz &
    )

    complex(dp), intent(inout), dimension(:, :, :) :: div_u
    complex(dp), intent(in), dimension(:, :, :) :: waves
    real(dp), intent(in), dimension(:) :: ax, bx, ay, by, az, bz
    integer, intent(in) :: nx_spec, ny_spec, nz_spec
    integer, intent(in) :: x_sp_st, y_sp_st, z_sp_st
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k, ix, iy, iz
    real(dp) :: tmp_r, tmp_c, div_r, div_c

    !$omp parallel do private(div_r, div_c, ix, iy, iz, tmp_r, tmp_c) collapse(3)
    do k = 1, nz_spec
      do j = 1, ny_spec
        do i = 1, nx_spec
          ! normalisation
          div_r = real(div_u(i, j, k), kind=dp)/nx/ny/nz
          div_c = aimag(div_u(i, j, k))/nx/ny/nz

          ix = i + x_sp_st
          iy = j + y_sp_st
          iz = k + z_sp_st

          ! post-process forward
          ! post-process in z
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(iz) + tmp_c*az(iz)
          div_c = tmp_c*bz(iz) - tmp_r*az(iz)
          if (iz > nz/2 + 1) div_r = -div_r
          if (iz > nz/2 + 1) div_c = -div_c

          ! post-process in y
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*by(iy) + tmp_c*ay(iy)
          div_c = tmp_c*by(iy) - tmp_r*ay(iy)
          if (iy > ny/2 + 1) div_r = -div_r
          if (iy > ny/2 + 1) div_c = -div_c

          ! post-process in x
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bx(ix) + tmp_c*ax(ix)
          div_c = tmp_c*bx(ix) - tmp_r*ax(ix)

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
          div_r = tmp_r*bz(iz) - tmp_c*az(iz)
          div_c = -tmp_c*bz(iz) - tmp_r*az(iz)
          if (iz > nz/2 + 1) div_r = -div_r
          if (iz > nz/2 + 1) div_c = -div_c

          ! post-process in y
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*by(iy) + tmp_c*ay(iy)
          div_c = tmp_c*by(iy) - tmp_r*ay(iy)
          if (iy > ny/2 + 1) div_r = -div_r
          if (iy > ny/2 + 1) div_c = -div_c

          ! post-process in x
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bx(ix) + tmp_c*ax(ix)
          div_c = -tmp_c*bx(ix) + tmp_r*ax(ix)

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

end module m_omp_spectral
