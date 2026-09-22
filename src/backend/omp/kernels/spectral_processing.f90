module m_omp_spectral
  use m_common, only: dp
  implicit none

  ! Threshold below which a spectral coefficient is treated as zero. Must
  ! scale with the working precision: single-precision roundoff leaves
  ! ~1e-13 noise in wave numbers that should vanish (e.g. Nyquist modes),
  ! and dividing by that noise corrupts the pressure field.
  real(dp), parameter :: eps_wave = epsilon(1._dp)

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
          if ((tmp_r < eps_wave) .or. (tmp_c < eps_wave)) then
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
          if (abs(tmp_r) < eps_wave) then
            div_r = 0._dp
          else
            div_r = -div_r/tmp_r
          end if
          if (abs(tmp_c) < eps_wave) then
            div_c = 0._dp
          else
            div_c = -div_c/tmp_c
          end if

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
          if (i + x_sp_st == nx/2 + 1 .and. k + z_sp_st == nz/2 + 1) &
            div_u(i, j, k) = 0._dp
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

  subroutine process_spectral_100( &
    div_u, waves, nx_spec, ny_spec, nz_spec, x_sp_st, y_sp_st, z_sp_st, &
    nx, ny, nz, ax, bx, ay, by, az, bz &
    )
    !! Post-process div U* in spectral space, for non-periodic BC in x-dir.
    !!
    !! This is the process_spectral_010 arrangement with x and y swapped:
    !! the caller passes the transposed buffer (dim1 = y modes, dim2 = x
    !! modes, dim3 = z), so this wrapper maps the global sizes nx <-> ny
    !! and the coefficient pairs ax, bx <-> ay, by onto process_spectral_010's
    !! non-periodic-y convention before delegating to it.
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

    call process_spectral_010(div_u=div_u, waves=waves, nx_spec=nx_spec, &
                              ny_spec=ny_spec, nz_spec=nz_spec, &
                              x_sp_st=x_sp_st, y_sp_st=y_sp_st, &
                              z_sp_st=z_sp_st, nx=ny, ny=nx, nz=nz, &
                              ax=ay, bx=by, ay=ax, by=bx, az=az, bz=bz)

  end subroutine process_spectral_100

  subroutine process_spectral_110_fw( &
    div_u, n1, n2, n3, st1, st2, st3, nx, ny, nz, ax, bx, az, bz &
    )
    !! Forward third of the spectral post-process for the 110 case
    !! (non-periodic x and y, periodic z, R2C Z-transpose layout):
    !! normalise, Z periodic post-process forward, then the X paired
    !! even/odd split.
    !!
    !! Layout: div_u(n1, n2, n3) with dim1 = Z r2c modes (az, bz), dim2 =
    !! X modes (ax, bx), dim3 = Y modes. Dim2 must span the whole global
    !! x range (nx), which the caller guarantees, since the X pairing
    !! writes both a mode and its mirror in the same pass.
    implicit none

    !> Divergence of velocity in spectral space, (nz/2+1, nx, ny) layout
    complex(dp), intent(inout), dimension(:, :, :) :: div_u
    real(dp), intent(in), dimension(:) :: ax, bx, az, bz
    !> Local extents of div_u
    integer, intent(in) :: n1, n2, n3
    !> Offsets of the local block within the global spectral array
    integer, intent(in) :: st1, st2, st3
    !> Global cell size
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k, iz, ix, ix_pair
    real(dp) :: tmp_r, tmp_c, div_r, div_c, l_r, l_c, r_r, r_c

    ! normalise + Z periodic post-process (forward). No sign flip: dim1
    ! never exceeds nz/2+1.
    !$omp parallel do private(iz, div_r, div_c, tmp_r, tmp_c) collapse(3)
    do k = 1, n3
      do j = 1, n2
        do i = 1, n1
          iz = i + st1

          div_r = real(div_u(i, j, k), kind=dp)/nx/ny/nz
          div_c = aimag(div_u(i, j, k))/nx/ny/nz

          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(iz) + tmp_c*az(iz)
          div_c = tmp_c*bz(iz) - tmp_r*az(iz)

          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

    ! X paired even/odd split (forward)
    !$omp parallel do private(ix, ix_pair, l_r, l_c, r_r, r_c) collapse(3)
    do k = 1, n3
      do j = 2, nx/2 + 1
        do i = 1, n1
          ix = j + st2
          ix_pair = nx - j + 2

          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, ix_pair, k), kind=dp)
          r_c = aimag(div_u(i, ix_pair, k))

          div_u(i, j, k) = 0.5_dp*cmplx( & !&
            l_r*bx(ix) + l_c*ax(ix) + r_r*bx(ix) - r_c*ax(ix), &
            -l_r*ax(ix) + l_c*bx(ix) + r_r*ax(ix) + r_c*bx(ix), kind=dp &
            )
          div_u(i, ix_pair, k) = 0.5_dp*cmplx( & !&
            r_r*bx(ix_pair + st2) + r_c*ax(ix_pair + st2) &
              + l_r*bx(ix_pair + st2) - l_c*ax(ix_pair + st2), &
            -r_r*ax(ix_pair + st2) + r_c*bx(ix_pair + st2) &
              + l_r*ax(ix_pair + st2) + l_c*bx(ix_pair + st2), &
            kind=dp &
            )
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_fw

  subroutine process_spectral_110_solve( &
    div_u, waves, n1, n2, n3, st1, st2, st3, nx, ny, nz, ay, by &
    )
    !! Middle third of the spectral post-process for the 110 case: the Y
    !! paired even/odd split, the Poisson divide, and the Y paired
    !! even/odd recombine.
    !!
    !! Layout: div_u(n1, n2, n3) with dim1 = Z r2c modes, dim2 = X modes,
    !! dim3 = Y modes (ay, by). Dim3 must span the whole global y range
    !! (ny), which the caller guarantees, since the Y pairing writes both
    !! a mode and its mirror in the same pass.
    implicit none

    !> Divergence of velocity in spectral space, (nz/2+1, nx, ny) layout
    complex(dp), intent(inout), dimension(:, :, :) :: div_u
    !> Spectral equivalence constants, (nz/2+1, nx, ny) layout
    complex(dp), intent(in), dimension(:, :, :) :: waves
    real(dp), intent(in), dimension(:) :: ay, by
    !> Local extents of div_u
    integer, intent(in) :: n1, n2, n3
    !> Offsets of the local block within the global spectral array
    integer, intent(in) :: st1, st2, st3
    !> Global cell size
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k, iy, iy_pair
    real(dp) :: div_r, div_c, tmp_r, tmp_c, l_r, l_c, r_r, r_c

    ! Y paired even/odd split (forward)
    !$omp parallel do private(iy, iy_pair, l_r, l_c, r_r, r_c) collapse(3)
    do k = 2, ny/2 + 1
      do j = 1, n2
        do i = 1, n1
          iy = k + st3
          iy_pair = ny - k + 2

          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, j, iy_pair), kind=dp)
          r_c = aimag(div_u(i, j, iy_pair))

          div_u(i, j, k) = 0.5_dp*cmplx( & !&
            l_r*by(iy) + l_c*ay(iy) + r_r*by(iy) - r_c*ay(iy), &
            -l_r*ay(iy) + l_c*by(iy) + r_r*ay(iy) + r_c*by(iy), kind=dp &
            )
          div_u(i, j, iy_pair) = 0.5_dp*cmplx( & !&
            r_r*by(iy_pair + st3) + r_c*ay(iy_pair + st3) &
              + l_r*by(iy_pair + st3) - l_c*ay(iy_pair + st3), &
            -r_r*ay(iy_pair + st3) + r_c*by(iy_pair + st3) &
              + l_r*ay(iy_pair + st3) + l_c*by(iy_pair + st3), &
            kind=dp &
            )
        end do
      end do
    end do
    !$omp end parallel do

    ! Solve Poisson
    !$omp parallel do private(div_r, div_c, tmp_r, tmp_c) collapse(3)
    do k = 1, n3
      do j = 1, n2
        do i = 1, n1
          div_r = real(div_u(i, j, k), kind=dp)
          div_c = aimag(div_u(i, j, k))

          tmp_r = real(waves(i, j, k), kind=dp)
          tmp_c = aimag(waves(i, j, k))
          if (abs(tmp_r) < eps_wave) then
            div_r = 0._dp
          else
            div_r = -div_r/tmp_r
          end if
          if (abs(tmp_c) < eps_wave) then
            div_c = 0._dp
          else
            div_c = -div_c/tmp_c
          end if

          ! update the entry
          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
          ! Zero Nyquist modes
          if (j + st2 == nx/2 + 1 .and. i + st1 == nz/2 + 1) &
            div_u(i, j, k) = 0._dp
        end do
      end do
    end do
    !$omp end parallel do

    ! Y paired even/odd recombine (backward)
    !$omp parallel do private(iy, iy_pair, l_r, l_c, r_r, r_c) collapse(3)
    do k = 2, ny/2 + 1
      do j = 1, n2
        do i = 1, n1
          iy = k + st3
          iy_pair = ny - k + 2

          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, j, iy_pair), kind=dp)
          r_c = aimag(div_u(i, j, iy_pair))

          div_u(i, j, k) = cmplx( & !&
            l_r*by(iy) - l_c*ay(iy) + r_r*ay(iy) + r_c*by(iy), &
            l_r*ay(iy) + l_c*by(iy) - r_r*by(iy) + r_c*ay(iy), kind=dp &
            )
          div_u(i, j, iy_pair) = cmplx( & !&
            r_r*by(iy_pair + st3) - r_c*ay(iy_pair + st3) &
              + l_r*ay(iy_pair + st3) + l_c*by(iy_pair + st3), &
            r_r*ay(iy_pair + st3) + r_c*by(iy_pair + st3) &
              - l_r*by(iy_pair + st3) + l_c*ay(iy_pair + st3), &
            kind=dp &
            )
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_solve

  subroutine process_spectral_110_bw( &
    div_u, n1, n2, n3, st1, st2, st3, nx, ny, nz, ax, bx, az, bz &
    )
    !! Backward third of the spectral post-process for the 110 case: the
    !! X paired even/odd recombine, then the Z periodic post-process
    !! backward.
    !!
    !! Layout: div_u(n1, n2, n3) with dim1 = Z r2c modes (az, bz), dim2 =
    !! X modes (ax, bx), dim3 = Y modes. Dim2 must span the whole global
    !! x range (nx), which the caller guarantees, since the X pairing
    !! writes both a mode and its mirror in the same pass.
    implicit none

    !> Divergence of velocity in spectral space, (nz/2+1, nx, ny) layout
    complex(dp), intent(inout), dimension(:, :, :) :: div_u
    real(dp), intent(in), dimension(:) :: ax, bx, az, bz
    !> Local extents of div_u
    integer, intent(in) :: n1, n2, n3
    !> Offsets of the local block within the global spectral array
    integer, intent(in) :: st1, st2, st3
    !> Global cell size
    integer, intent(in) :: nx, ny, nz

    integer :: i, j, k, iz, ix, ix_pair
    real(dp) :: tmp_r, tmp_c, div_r, div_c, l_r, l_c, r_r, r_c

    ! X paired even/odd recombine (backward)
    !$omp parallel do private(ix, ix_pair, l_r, l_c, r_r, r_c) collapse(3)
    do k = 1, n3
      do j = 2, nx/2 + 1
        do i = 1, n1
          ix = j + st2
          ix_pair = nx - j + 2

          l_r = real(div_u(i, j, k), kind=dp)
          l_c = aimag(div_u(i, j, k))
          r_r = real(div_u(i, ix_pair, k), kind=dp)
          r_c = aimag(div_u(i, ix_pair, k))

          div_u(i, j, k) = cmplx( & !&
            l_r*bx(ix) - l_c*ax(ix) + r_r*ax(ix) + r_c*bx(ix), &
            l_r*ax(ix) + l_c*bx(ix) - r_r*bx(ix) + r_c*ax(ix), kind=dp &
            )
          div_u(i, ix_pair, k) = cmplx( & !&
            r_r*bx(ix_pair + st2) - r_c*ax(ix_pair + st2) &
              + l_r*ax(ix_pair + st2) + l_c*bx(ix_pair + st2), &
            r_r*ax(ix_pair + st2) + r_c*bx(ix_pair + st2) &
              - l_r*bx(ix_pair + st2) + l_c*ax(ix_pair + st2), &
            kind=dp &
            )
        end do
      end do
    end do
    !$omp end parallel do

    ! Z periodic post-process (backward). No sign flip: dim1 never
    ! exceeds nz/2+1.
    !$omp parallel do private(iz, div_r, div_c, tmp_r, tmp_c) collapse(3)
    do k = 1, n3
      do j = 1, n2
        do i = 1, n1
          iz = i + st1

          div_r = real(div_u(i, j, k), kind=dp)
          div_c = aimag(div_u(i, j, k))

          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*bz(iz) - tmp_c*az(iz)
          div_c = tmp_c*bz(iz) + tmp_r*az(iz)

          div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine process_spectral_110_bw

  subroutine process_spectral_110( &
    div_u, waves, n1, n2, n3, st1, st2, st3, nx, ny, nz, &
    ax, bx, ay, by, az, bz &
    )
    !! Post-process div U* in spectral space, for non-periodic BC in both
    !! x- and y-directions (R2C Z-transpose layout: dim1 = Z r2c modes,
    !! dim2 = X modes, dim3 = Y modes).
    implicit none

    !> Divergence of velocity in spectral space, (nz/2+1, nx, ny) layout
    complex(dp), intent(inout), dimension(:, :, :) :: div_u
    !> Spectral equivalence constants, (nz/2+1, nx, ny) layout
    complex(dp), intent(in), dimension(:, :, :) :: waves
    real(dp), intent(in), dimension(:) :: ax, bx, ay, by, az, bz
    !> Local extents of div_u
    integer, intent(in) :: n1, n2, n3
    !> Offsets of the local block within the global spectral array
    integer, intent(in) :: st1, st2, st3
    !> Global cell size
    integer, intent(in) :: nx, ny, nz

    call process_spectral_110_fw(div_u, n1, n2, n3, st1, st2, st3, &
                                 nx, ny, nz, ax, bx, az, bz)
    call process_spectral_110_solve(div_u, waves, n1, n2, n3, st1, st2, st3, &
                                    nx, ny, nz, ay, by)
    call process_spectral_110_bw(div_u, n1, n2, n3, st1, st2, st3, &
                                 nx, ny, nz, ax, bx, az, bz)

  end subroutine process_spectral_110

end module m_omp_spectral
