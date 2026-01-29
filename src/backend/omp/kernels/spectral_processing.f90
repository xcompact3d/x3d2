module m_omp_spectral
  !! Spectral space processing for FFT-based Poisson solver.
  !!
  !! Provides kernels for solving Poisson equation in Fourier space with
  !! spectral equivalence transformations. Handles different boundary
  !! condition combinations: fully periodic (000) and Dirichlet in Y (010).
  !!
  !! **Spectral equivalence:** Modified wavenumbers for finite-difference
  !! grid (Lele 1992). Ensures spectral solver matches compact FD schemes.
  !!
  !! **Reference:** JCP 228 (2009), 5989-6015, Section 4
  !!
  !! **Processing steps:**
  !! 1. Forward spectral equivalence transform (physical $\rightarrow$ modified wavenumbers)
  !! 2. Solve: $\hat{\phi}_k = -\hat{f}_k / k^2$
  !! 3. Backward spectral equivalence transform (modified wavenumbers $\rightarrow$ physical)
  use m_common, only: dp
  implicit none

contains

  subroutine process_spectral_000( &
    div_u, waves, nx_spec, ny_spec, nz_spec, x_sp_st, y_sp_st, z_sp_st, &
    nx, ny, nz, ax, bx, ay, by, az, bz &
    )
    !! Solve Poisson in spectral space for (0,0,0) boundary conditions.
    !!
    !! Processes fully periodic case. Applies spectral equivalence transforms
    !! in all three directions, divides by squared wavenumber, then applies
    !! inverse transforms.
    !!
    !! **Algorithm:**
    !! 1. Normalise by grid size (FFT convention)
    !! 2. Forward spectral equivalence: physical $\rightarrow$ modified waves (Z, Y, X order)
    !! 3. Solve: $\phi_k = -f_k / k^2$ (handle zero mode specially)
    !! 4. Backward spectral equivalence: modified waves $\rightarrow$ physical
    !!
    !! **Special case:** Zero wavenumber (k=0) set to zero to remove constant mode.
    !!
    !! **Ref.** JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    !> Divergence of velocity in spectral space
    complex(dp), intent(inout), dimension(:, :, :) :: div_u  !! In: RHS, Out: Solution
    !> Spectral equivalence constants
    complex(dp), intent(in), dimension(:, :, :) :: waves  !! Modified wavenumbers squared
    real(dp), intent(in), dimension(:) :: ax, bx, ay, by, az, bz  !! Spectral equivalence coefficients
    !> Grid size in spectral space
    integer, intent(in) :: nx_spec, ny_spec, nz_spec  !! Local spectral dimensions
    !> Offsets in the permuted pencils in spectral space
    integer, intent(in) :: x_sp_st, y_sp_st, z_sp_st  !! Global offsets
    !> Global cell size
    integer, intent(in) :: nx, ny, nz  !! Global grid dimensions

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
    !! Solve Poisson in spectral space for (0,1,0) boundary conditions.
    !!
    !! Processes Dirichlet in Y, periodic in X and Z. Uses sine series
    !! in Y-direction (symmetry/antisymmetry transform) combined with
    !! Fourier in X and Z.
    !!
    !! **Algorithm:**
    !! 1. Normalise by grid size
    !! 2. Forward spectral equivalence in Z and X (not Y, handled separately)
    !! 3. Apply Y symmetry transform (combine left/right halves)
    !! 4. Solve: $\phi_k = -f_k / k^2$
    !! 5. Inverse Y symmetry transform
    !! 6. Backward spectral equivalence in X and Z
    !!
    !! **Y-direction:** Sine series requires special symmetric processing
    !! to maintain real-valued solution with Dirichlet BCs.
    !!
    !! **Ref.** JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    !> Divergence of velocity in spectral space
    complex(dp), intent(inout), dimension(:, :, :) :: div_u  !! In: RHS, Out: Solution
    !> Spectral equivalence constants
    complex(dp), intent(in), dimension(:, :, :) :: waves  !! Modified wavenumbers squared
    real(dp), intent(in), dimension(:) :: ax, bx, ay, by, az, bz  !! Spectral equivalence coefficients
    !> Grid size in spectral space
    integer, intent(in) :: nx_spec, ny_spec, nz_spec  !! Local spectral dimensions
    !> Offsets in the permuted pencils in spectral space
    integer, intent(in) :: x_sp_st, y_sp_st, z_sp_st  !! Global offsets
    !> Global cell size
    integer, intent(in) :: nx, ny, nz  !! Global grid dimensions

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

end module m_omp_spectral
