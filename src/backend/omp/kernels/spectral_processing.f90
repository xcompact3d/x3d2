module m_omp_spectral
  use m_common, only: dp
  implicit none

contains

  subroutine process_spectral_000( &
    div_u, waves, nx_spec, ny_spec, nz_spec, x_sp_st, y_sp_st, z_sp_st, &
    nx, ny, nz, ax, bx, ay, by, az, bz &
    )
    !! Post-process div U* in spectral space for all periodic BCs.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    !> Divergence of velocity in spectral space
    complex(dp), intent(inout), dimension(:, :, :) :: div_u
    !> Spectral equivalence constants
    complex(dp), intent(in), dimension(:, :, :) :: waves
    real(dp), intent(in), dimension(:) :: ax, bx, ay, by, az, bz
    !> Grid size in spectral space
    integer, intent(in) :: nx_spec, ny_spec, nz_spec
    !> Offsets in the permuted pencils in spectral space
    integer, intent(in) :: x_sp_st, y_sp_st, z_sp_st
    !> Global cell size
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

  end subroutine process_spectral_000

  subroutine process_spectral_010( &
    div_u, waves, nx_spec, ny_spec, nz_spec, x_sp_st, y_sp_st, z_sp_st, &
    nx, ny, nz, ax, bx, ay, by, az, bz &
    )
    !! Post-process div U* in spectral space, for non-periodic BC in y-dir.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    !> Divergence of velocity in spectral space
    complex(dp), intent(inout), dimension(:, :, :) :: div_u
    !> Spectral equivalence constants
    complex(dp), intent(in), dimension(:, :, :) :: waves
    real(dp), intent(in), dimension(:) :: ax, bx, ay, by, az, bz
    !> Grid size in spectral space
    integer, intent(in) :: nx_spec, ny_spec, nz_spec
    !> Offsets in the permuted pencils in spectral space
    integer, intent(in) :: x_sp_st, y_sp_st, z_sp_st
    !> Global cell size
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k, ix, iy, iz, iy_r
    real(dp) :: tmp_r, tmp_c, div_r, div_c, l_r, l_c, r_r, r_c

    !$omp parallel do private(div_r, div_c, ix, iz, tmp_r, tmp_c) collapse(3)
    do k = 1, nz_spec
      do j = 1, ny_spec
        do i = 1, nx_spec
          ix = i + x_sp_st
          iz = k + z_sp_st

          ! normalisation
          div_r = real(div_u(i, j, k), kind=dp)/nx/ny/nz
          div_c = aimag(div_u(i, j, k))/nx/ny/nz

          ! postprocess in z
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(iz) + tmp_c*az(iz)
          div_c = tmp_c*bz(iz) - tmp_r*az(iz)
          if (iz > nz/2 + 1) div_r = -div_r
          if (iz > nz/2 + 1) div_c = -div_c

          ! postprocess in x
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bx(ix) + tmp_c*ax(ix)
          div_c = tmp_c*bx(ix) - tmp_r*ax(ix)
          if (ix > nx/2 + 1) div_r = -div_r
          if (ix > nx/2 + 1) div_c = -div_c

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(div_r, div_c, iy, iy_r, l_r, l_c, r_r, r_c) collapse(3)
    do k = 1, nz_spec
      do j = 2, ny_spec/2 + 1
        do i = 1, nx_spec
          iy = j + y_sp_st
          iy_r = ny_spec - j + 2 + y_sp_st

          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, ny_spec - j + 2, k), kind=dp)
          r_c = aimag(div_u(i, ny_spec - j + 2, k))

          ! update the entry
          div_u(i, j, k) = 0.5_dp*cmplx( & !&
            l_r*by(iy) + l_c*ay(iy) + r_r*by(iy) - r_c*ay(iy), &
            -l_r*ay(iy) + l_c*by(iy) + r_r*ay(iy) + r_c*by(iy), kind=dp &
            )
          div_u(i, ny_spec - j + 2, k) = 0.5_dp*cmplx( & !&
            r_r*by(iy_r) + r_c*ay(iy_r) + l_r*by(iy_r) - l_c*ay(iy_r), &
            -r_r*ay(iy_r) + r_c*by(iy_r) + l_r*ay(iy_r) + l_c*by(iy_r), &
           kind=dp &
           )
        end do
      end do
    end do
    !$omp end parallel do

    ! Solve Poisson
    !$omp parallel do private(div_r, div_c, tmp_r, tmp_c) collapse(3)
    do k = 1, nz_spec
      do j = 1, ny_spec
        do i = 1, nx_spec
          div_r = real(div_u(i, j, k), kind=dp)
          div_c = aimag(div_u(i, j, k))

          tmp_r = real(waves(i, j, k), kind=dp)
          tmp_c = aimag(waves(i, j, k))
          if (abs(tmp_r) < 1.e-16_dp) then
            div_r = 0._dp
          else
            div_r = -div_r/tmp_r
          end if
          if (abs(tmp_c) < 1.e-16_dp) then
            div_c = 0._dp
          else
            div_c = -div_c/tmp_c
          end if

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
          if (i == nx/2 + 1 .and. k == nz/2 + 1) div_u(i, j, k) = 0._dp
        end do
      end do
    end do
    !$omp end parallel do

    ! post-process backward
    !$omp parallel do private(div_r, div_c, iy, iy_r, l_r, l_c, r_r, r_c) collapse(3)
    do k = 1, nz_spec
      do j = 2, ny_spec/2 + 1
        do i = 1, nx_spec
          iy = j + y_sp_st
          iy_r = ny_spec - j + 2 + y_sp_st

          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, ny_spec - j + 2, k), kind=dp)
          r_c = aimag(div_u(i, ny_spec - j + 2, k))

          ! update the entry
          div_u(i, j, k) = cmplx( & !&
            l_r*by(iy) - l_c*ay(iy) + r_r*ay(iy) + r_c*by(iy), &
            l_r*ay(iy) + l_c*by(iy) - r_r*by(iy) + r_c*ay(iy), kind=dp &
            )
          div_u(i, ny_spec - j + 2, k) = cmplx( & !&
            r_r*by(iy_r) - r_c*ay(iy_r) + l_r*ay(iy_r) + l_c*by(iy_r), &
            r_r*ay(iy_r) + r_c*by(iy_r) - l_r*by(iy_r) + l_c*ay(iy_r), &
            kind=dp &
            )
        end do
      end do
    end do
    !$omp end parallel do

    !$omp parallel do private(div_r, div_c, ix, iz, tmp_r, tmp_c) collapse(3)
    do k = 1, nz_spec
      do j = 1, ny_spec
        do i = 1, nx_spec
          ix = i + x_sp_st
          iz = k + z_sp_st

          div_r = real(div_u(i, j, k), kind=dp)
          div_c = aimag(div_u(i, j, k))

          ! post-process in z
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(iz) - tmp_c*az(iz)
          div_c = tmp_c*bz(iz) + tmp_r*az(iz)
          if (iz > nz/2 + 1) div_r = -div_r
          if (iz > nz/2 + 1) div_c = -div_c

          ! post-process in x
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bx(ix) - tmp_c*ax(ix)
          div_c = tmp_c*bx(ix) + tmp_r*ax(ix)
          if (ix > nx/2 + 1) div_r = -div_r
          if (ix > nx/2 + 1) div_c = -div_c

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_010

! ===========================================================================
  ! 110 spectral kernels (OMP, natural (nx, ny, nz/2+1) layout).
  ! Seven phases mirroring the CUDA driver, but each is one OMP parallel-do
  ! block. Between phases we rely on the implicit join at the end of each
  ! parallel-do; non-local pair dependencies between phases need that
  ! barrier and cannot be fused into a single loop.
  ! ===========================================================================

  subroutine process_spectral_110_norm_z( &
    div_u, nx, ny, nz_h, nz, az, bz)
    !! Step 1: normalise + Z periodic post-process (forward).
    !! Element-local: each (i, j, k) updates independently.
    implicit none
    complex(dp), intent(inout) :: div_u(:, :, :)   ! (nx, ny, nz/2+1)
    real(dp), intent(in) :: az(:), bz(:)
    integer, intent(in) :: nx, ny, nz_h, nz

    integer :: i, j, k
    real(dp) :: tmp_r, tmp_c, div_r, div_c

    !$omp parallel do collapse(2) private(i, k, div_r, div_c, tmp_r, tmp_c)
    do j = 1, ny
      do i = 1, nx
        do k = 1, nz_h
          div_r = real(div_u(i, j, k), kind=dp)/(nx*ny*nz)
          div_c = aimag(div_u(i, j, k))/(nx*ny*nz)
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(k) + tmp_c*az(k)
          div_c = tmp_c*bz(k) - tmp_r*az(k)
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_norm_z

  subroutine process_spectral_110_x_pair_fw( &
    div_u, nx, ny, nz_h, x_sp_st, ax, bx)
    !! Step 2 (forward): X paired even/odd split. Only i in [2, nx/2+1] runs;
    !! pairs i with ix_pair = nx - i + 2. At i == nx/2+1 (even nx) the pair
    !! self-coincides and both writes produce the same value -- harmless.
    implicit none
    complex(dp), intent(inout) :: div_u(:, :, :)
    real(dp), intent(in) :: ax(:), bx(:)
    integer, intent(in) :: nx, ny, nz_h, x_sp_st

    integer :: i, j, k, ix, ix_pair
    real(dp) :: l_r, l_c, r_r, r_c

    !$omp parallel do collapse(2) private(i, ix, ix_pair, l_r, l_c, r_r, r_c)
    do k = 1, nz_h
      do j = 1, ny
        do i = 2, nx/2 + 1
          ix = i + x_sp_st
          ix_pair = nx - i + 2
          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(ix_pair, j, k), kind=dp)
          r_c = aimag(div_u(ix_pair, j, k))
          div_u(i, j, k) = 0.5_dp*cmplx( &
            l_r*bx(ix) + l_c*ax(ix) + r_r*bx(ix) - r_c*ax(ix), &
            -l_r*ax(ix) + l_c*bx(ix) + r_r*ax(ix) + r_c*bx(ix), kind=dp)
          div_u(ix_pair, j, k) = 0.5_dp*cmplx( &
            r_r*bx(ix_pair + x_sp_st) + r_c*ax(ix_pair + x_sp_st) &
              + l_r*bx(ix_pair + x_sp_st) - l_c*ax(ix_pair + x_sp_st), &
            -r_r*ax(ix_pair + x_sp_st) + r_c*bx(ix_pair + x_sp_st) &
              + l_r*ax(ix_pair + x_sp_st) + l_c*bx(ix_pair + x_sp_st), &
            kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_x_pair_fw

  subroutine process_spectral_110_y_pair_fw( &
    div_u, nx, ny, nz_h, y_sp_st, ay, by)
    !! Step 3 (forward): Y paired even/odd split.
    implicit none
    complex(dp), intent(inout) :: div_u(:, :, :)
    real(dp), intent(in) :: ay(:), by(:)
    integer, intent(in) :: nx, ny, nz_h, y_sp_st

    integer :: i, j, k, iy, iy_pair
    real(dp) :: l_r, l_c, r_r, r_c

    !$omp parallel do collapse(2) private(j, iy, iy_pair, l_r, l_c, r_r, r_c)
    do k = 1, nz_h
      do i = 1, nx
        do j = 2, ny/2 + 1
          iy = j + y_sp_st
          iy_pair = ny - j + 2
          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, iy_pair, k), kind=dp)
          r_c = aimag(div_u(i, iy_pair, k))
          div_u(i, j, k) = 0.5_dp*cmplx( &
            l_r*by(iy) + l_c*ay(iy) + r_r*by(iy) - r_c*ay(iy), &
            -l_r*ay(iy) + l_c*by(iy) + r_r*ay(iy) + r_c*by(iy), kind=dp)
          div_u(i, iy_pair, k) = 0.5_dp*cmplx( &
            r_r*by(iy_pair + y_sp_st) + r_c*ay(iy_pair + y_sp_st) &
              + l_r*by(iy_pair + y_sp_st) - l_c*ay(iy_pair + y_sp_st), &
            -r_r*ay(iy_pair + y_sp_st) + r_c*by(iy_pair + y_sp_st) &
              + l_r*ay(iy_pair + y_sp_st) + l_c*by(iy_pair + y_sp_st), &
            kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_y_pair_fw

  subroutine process_spectral_110_poisson( &
    div_u, waves, nx, ny, nz_h, nz, x_sp_st)
    !! Step 4: divide div_u by waves componentwise (real and imag parts
    !! separately, with an epsilon guard). Zero the Nyquist self-pair.
    implicit none
    complex(dp), intent(inout) :: div_u(:, :, :)
    complex(dp), intent(in) :: waves(:, :, :)
    integer, intent(in) :: nx, ny, nz_h, nz, x_sp_st

    integer :: i, j, k, ix
    real(dp) :: div_r, div_c, tmp_r, tmp_c

    !$omp parallel do collapse(2) private(i, k, ix, div_r, div_c, tmp_r, tmp_c)
    do j = 1, ny
      do i = 1, nx
        ix = i + x_sp_st
        do k = 1, nz_h
          div_r = real(div_u(i, j, k), kind=dp)
          div_c = aimag(div_u(i, j, k))
          tmp_r = real(waves(i, j, k), kind=dp)
          tmp_c = aimag(waves(i, j, k))
          if (abs(tmp_r) < 1.e-16_dp) then
            div_r = 0._dp
          else
            div_r = -div_r/tmp_r
          end if
          if (abs(tmp_c) < 1.e-16_dp) then
            div_c = 0._dp
          else
            div_c = -div_c/tmp_c
          end if
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
          if (ix == nx/2 + 1 .and. k == nz/2 + 1) div_u(i, j, k) = 0._dp
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_poisson

  subroutine process_spectral_110_y_pair_bw( &
    div_u, nx, ny, nz_h, y_sp_st, ay, by)
    !! Step 5 (backward): Y paired even/odd recombine.
    implicit none
    complex(dp), intent(inout) :: div_u(:, :, :)
    real(dp), intent(in) :: ay(:), by(:)
    integer, intent(in) :: nx, ny, nz_h, y_sp_st

    integer :: i, j, k, iy, iy_pair
    real(dp) :: l_r, l_c, r_r, r_c

    !$omp parallel do collapse(2) private(j, iy, iy_pair, l_r, l_c, r_r, r_c)
    do k = 1, nz_h
      do i = 1, nx
        do j = 2, ny/2 + 1
          iy = j + y_sp_st
          iy_pair = ny - j + 2
          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, iy_pair, k), kind=dp)
          r_c = aimag(div_u(i, iy_pair, k))
          div_u(i, j, k) = cmplx( &
            l_r*by(iy) - l_c*ay(iy) + r_r*ay(iy) + r_c*by(iy), &
            l_r*ay(iy) + l_c*by(iy) - r_r*by(iy) + r_c*ay(iy), kind=dp)
          div_u(i, iy_pair, k) = cmplx( &
            r_r*by(iy_pair + y_sp_st) - r_c*ay(iy_pair + y_sp_st) &
              + l_r*ay(iy_pair + y_sp_st) + l_c*by(iy_pair + y_sp_st), &
            r_r*ay(iy_pair + y_sp_st) + r_c*by(iy_pair + y_sp_st) &
              - l_r*by(iy_pair + y_sp_st) + l_c*ay(iy_pair + y_sp_st), &
            kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_y_pair_bw

  subroutine process_spectral_110_x_pair_bw( &
    div_u, nx, ny, nz_h, x_sp_st, ax, bx)
    !! Step 6 (backward): X paired even/odd recombine.
    implicit none
    complex(dp), intent(inout) :: div_u(:, :, :)
    real(dp), intent(in) :: ax(:), bx(:)
    integer, intent(in) :: nx, ny, nz_h, x_sp_st

    integer :: i, j, k, ix, ix_pair
    real(dp) :: l_r, l_c, r_r, r_c

    !$omp parallel do collapse(2) private(i, ix, ix_pair, l_r, l_c, r_r, r_c)
    do k = 1, nz_h
      do j = 1, ny
        do i = 2, nx/2 + 1
          ix = i + x_sp_st
          ix_pair = nx - i + 2
          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(ix_pair, j, k), kind=dp)
          r_c = aimag(div_u(ix_pair, j, k))
          div_u(i, j, k) = cmplx( &
            l_r*bx(ix) - l_c*ax(ix) + r_r*ax(ix) + r_c*bx(ix), &
            l_r*ax(ix) + l_c*bx(ix) - r_r*bx(ix) + r_c*ax(ix), kind=dp)
          div_u(ix_pair, j, k) = cmplx( &
            r_r*bx(ix_pair + x_sp_st) - r_c*ax(ix_pair + x_sp_st) &
              + l_r*ax(ix_pair + x_sp_st) + l_c*bx(ix_pair + x_sp_st), &
            r_r*ax(ix_pair + x_sp_st) + r_c*bx(ix_pair + x_sp_st) &
              - l_r*bx(ix_pair + x_sp_st) + l_c*ax(ix_pair + x_sp_st), &
            kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_x_pair_bw

  subroutine process_spectral_110_z_bw( &
    div_u, nx, ny, nz_h, az, bz)
    !! Step 7 (backward): Z periodic undo. Element-local.
    implicit none
    complex(dp), intent(inout) :: div_u(:, :, :)
    real(dp), intent(in) :: az(:), bz(:)
    integer, intent(in) :: nx, ny, nz_h

    integer :: i, j, k
    real(dp) :: tmp_r, tmp_c, div_r, div_c

    !$omp parallel do collapse(2) private(i, k, div_r, div_c, tmp_r, tmp_c)
    do j = 1, ny
      do i = 1, nx
        do k = 1, nz_h
          div_r = real(div_u(i, j, k), kind=dp)
          div_c = aimag(div_u(i, j, k))
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(k) - tmp_c*az(k)
          div_c = tmp_c*bz(k) + tmp_r*az(k)
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_z_bw

end module m_omp_spectral
