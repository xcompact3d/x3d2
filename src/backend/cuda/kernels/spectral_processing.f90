module m_cuda_spectral
  !! CUDA kernels for spectral space processing and FFT post-processing.
  !!
  !! This module contains kernels for:
  !!
  !! - Post-processing spectral transforms (forward/backward)
  !! - Solving Poisson equations in spectral space
  !! - Enforcing and undoing periodicity in Y-direction
  !!
  !! Implements spectral equivalence method from JCP 228 (2009), 5989-6015, Sec 4.
  !! Handles both periodic (000) and non-periodic (010) boundary conditions.
  use cudafor

  use m_common, only: dp

  implicit none

contains

  attributes(global) subroutine memcpy3D(dst, src, nx, ny, nz)
    !! Copy data between x3d2 padded arrays and cuFFTMp descriptors.
    !!
    !! Each thread handles one Y-Z plane position, looping over X.
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: dst  !! Output: Destination array
    real(dp), device, intent(in), dimension(:, :, :) :: src  !! Input: Source array
    integer, value, intent(in) :: nx, ny, nz  !! Grid dimensions

    integer :: i, j, k  !! Loop and thread indices

    j = threadIdx%x + (blockIdx%x - 1)*blockDim%x !ny
    k = blockIdx%y !nz

    if (j <= ny) then
      do i = 1, nx
        dst(i, j, k) = src(i, j, k)
      end do
    end if
  end subroutine memcpy3D

  attributes(global) subroutine process_spectral_000( &
    div_u, waves, nx_spec, ny_spec, y_sp_st, nx, ny, nz, &
    ax, bx, ay, by, az, bz &
    )
    !! Post-processes the divergence of velocity in spectral space, including
    !! scaling w.r.t. grid size.
    !!
    !! Performs forward post-processing, Poisson solve, and backward post-processing
    !! using spectral equivalence method. Ref: JCP 228 (2009), 5989-6015, Sec 4.
    implicit none

    complex(dp), device, intent(inout), dimension(:, :, :) :: div_u  !! In/out: Divergence of velocity in spectral space
    complex(dp), device, intent(in), dimension(:, :, :) :: waves  !! Input: Spectral wavenumbers for Poisson solve
    real(dp), device, intent(in), dimension(:) :: ax, bx, ay, by, az, bz  !! Input: Spectral equivalence constants
    integer, value, intent(in) :: nx_spec, ny_spec  !! Spectral space grid size
    integer, value, intent(in) :: y_sp_st  !! Y-direction offset in the permuted slabs in spectral space
    integer, value, intent(in) :: nx, ny, nz  !! Physical space grid size

    integer :: i, j, k, ix, iy, iz  !! Loop and spectral mode indices
    real(dp) :: tmp_r, tmp_c, div_r, div_c  !! Temporary real/imaginary components

    j = threadIdx%x + (blockIdx%x - 1)*blockDim%x
    k = blockIdx%y ! nz_spec

    if (j <= ny_spec) then
      do i = 1, nx_spec
        ! normalisation
        div_r = real(div_u(i, j, k), kind=dp)/nx/ny/nz
        div_c = aimag(div_u(i, j, k))/nx/ny/nz

        ix = i; iy = j + y_sp_st; iz = k

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
    end if

  end subroutine process_spectral_000

  attributes(global) subroutine process_spectral_010( &
    div_u, waves, nx_spec, ny_spec, y_sp_st, nx, ny, nz, &
    ax, bx, ay, by, az, bz &
    )
    !! Post-process divergence field and solve Poisson equation in spectral space
    !! for non-periodic boundary conditions in Y-direction (010).
    !!
    !! Performs forward post-processing with odd/even mode handling, Poisson solve,
    !! and backward post-processing. Ref: JCP 228 (2009), 5989-6015, Sec 4.
    implicit none

    complex(dp), device, intent(inout), dimension(:, :, :) :: div_u  !! In/out: Divergence field / pressure solution
    complex(dp), device, intent(in), dimension(:, :, :) :: waves  !! Input: Spectral wavenumbers for Poisson solve
    real(dp), device, intent(in), dimension(:) :: ax, bx, ay, by, az, bz  !! Input: Spectral equivalence constants
    integer, value, intent(in) :: nx_spec, ny_spec  !! Spectral space grid size
    integer, value, intent(in) :: y_sp_st  !! Y-direction offset in spectral slabs
    integer, value, intent(in) :: nx, ny, nz  !! Physical space grid size

    integer :: i, j, k, ix, iy, iz, iy_rev  !! Loop, spectral, and reversed mode indices
    real(dp) :: tmp_r, tmp_c, div_r, div_c, l_r, l_c, r_r, r_c  !! Temporary components for left/right modes

    i = threadIdx%x + (blockIdx%x - 1)*blockDim%x
    k = blockIdx%y ! nz_spec

    if (i <= nx_spec) then
      do j = 1, ny_spec
        ix = i; iy = j + y_sp_st; iz = k

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
    end if

    if (i <= nx_spec) then
      do j = 2, ny_spec/2 + 1
        ix = i; iy = j + y_sp_st; iz = k
        iy_rev = ny_spec - j + 2 + y_sp_st

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
         r_r*by(iy_rev) + r_c*ay(iy_rev) + l_r*by(iy_rev) - l_c*ay(iy_rev), &
         -r_r*ay(iy_rev) + r_c*by(iy_rev) + l_r*ay(iy_rev) + l_c*by(iy_rev), &
         kind=dp &
         )
      end do
    end if

    ! Solve Poisson
    if (i <= nx_spec) then
      do j = 1, ny_spec
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
    end if

    ! post-process backward
    if (i <= nx_spec) then
      do j = 2, ny_spec/2 + 1
        ix = i; iy = j + y_sp_st; iz = k
        iy_rev = ny_spec - j + 2 + y_sp_st

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
          r_r*by(iy_rev) - r_c*ay(iy_rev) + l_r*ay(iy_rev) + l_c*by(iy_rev), &
          r_r*ay(iy_rev) + r_c*by(iy_rev) - l_r*by(iy_rev) + l_c*ay(iy_rev), &
          kind=dp &
          )
      end do
    end if

    if (i <= nx_spec) then
      do j = 1, ny_spec
        ix = i; iy = j + y_sp_st; iz = k

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
    end if

  end subroutine process_spectral_010

  attributes(global) subroutine process_spectral_010_fw( &
    div_u, nx_spec, ny_spec, y_sp_st, nx, ny, nz, ax, bx, ay, by, az, bz &
    )
    !! Forward post-processing only for non-periodic Y-direction (010).
    !!
    !! Performs normalisation, post-processing in X and Z, and odd/even mode handling
    !! in Y. Used when Poisson solve and backward processing are separate steps.
    implicit none

    complex(dp), device, intent(inout), dimension(:, :, :) :: div_u  !! In/out: Divergence field to post-process
    real(dp), device, intent(in), dimension(:) :: ax, bx, ay, by, az, bz  !! Input: Spectral equivalence constants
    !> Grid size in spectral space
    integer, value, intent(in) :: nx_spec, ny_spec  !! Spectral space grid size
    !> Offset in y direction in the permuted slabs in spectral space
    integer, value, intent(in) :: y_sp_st  !! Y-direction offset in spectral slabs
    !> Grid size
    integer, value, intent(in) :: nx, ny, nz  !! Physical space grid size

    integer :: i, j, k, ix, iy, iz, iy_rev  !! Loop, spectral, and reversed mode indices
    real(dp) :: tmp_r, tmp_c, div_r, div_c, l_r, l_c, r_r, r_c  !! Temporary real/imaginary components

    i = threadIdx%x + (blockIdx%x - 1)*blockDim%x
    k = blockIdx%y ! nz_spec

    if (i <= nx_spec) then
      do j = 1, ny_spec
        ix = i; iy = j + y_sp_st; iz = k

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
    end if

    if (i <= nx_spec) then
      do j = 2, ny_spec/2 + 1
        ix = i; iy = j + y_sp_st; iz = k
        iy_rev = ny_spec - j + 2 + y_sp_st

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
         r_r*by(iy_rev) + r_c*ay(iy_rev) + l_r*by(iy_rev) - l_c*ay(iy_rev), &
         -r_r*ay(iy_rev) + r_c*by(iy_rev) + l_r*ay(iy_rev) + l_c*by(iy_rev), &
         kind=dp &
         )
      end do
    end if

  end subroutine process_spectral_010_fw

  attributes(global) subroutine process_spectral_010_poisson( &
    div_u, a_re, a_im, off, inc, nx_spec, n, nx, ny, nz &
    )
    !! Solve Poisson equation for non-periodic Y-direction using pentadiagonal solver.
    !!
    !! Handles odd/even mode separation using offset and increment parameters.
    !! Modifies pentadiagonal coefficients in-place during forward/backward passes.
    implicit none

    complex(dp), device, intent(inout), dimension(:, :, :) :: div_u  !! In/out: RHS / Solution
    real(dp), device, intent(inout), dimension(:, :, :, :) :: a_re, a_im  !! In/out: Pentadiagonal coefficients (real/imag)
    integer, value, intent(in) :: off, inc  !! Offset and increment for odd/even modes
    integer, value, intent(in) :: nx_spec, n, nx, ny, nz  !! Grid dimensions

    integer :: i, j, k, jm, nm  !! Loop indices and mapped indices
    real(dp) :: tmp_r, tmp_c, div_r, div_c, epsilon  !! Temporary variables and tolerance

    i = threadIdx%x + (blockIdx%x - 1)*blockDim%x
    k = blockIdx%y ! nz_spec

    epsilon = 1.e-16_dp

    ! Solve Poisson
    if (i <= nx_spec) then
      ! Forward pass for the pentadiagonal matrix
      do j = 1, n - 2
        ! j mapping based on odd/even
        ! inc=2, off=0 ---> j => 2j - 1
        ! inc=2, off=1 ---> j => 2j
        ! inc=1, off=0 ---> j => j
        jm = inc*j + off - inc/2
        ! eliminate diag-1
        tmp_r = 0._dp
        if (abs(a_re(i, j, k, 3)) > epsilon) then
          tmp_r = a_re(i, j + 1, k, 2)/a_re(i, j, k, 3)
        end if
        tmp_c = 0._dp
        if (abs(a_im(i, j, k, 3)) > epsilon) then
          tmp_c = a_im(i, j + 1, k, 2)/a_im(i, j, k, 3)
        end if
        div_r = real(div_u(i, jm + inc, k) - tmp_r*div_u(i, jm, k), kind=dp)
        div_c = aimag(div_u(i, jm + inc, k) - tmp_c*div_u(i, jm, k))
        div_u(i, jm + inc, k) = cmplx(div_r, div_c, kind=dp)
        ! modify pentadiagonal coefficients in-place
        a_re(i, j + 1, k, 3) = a_re(i, j + 1, k, 3) - tmp_r*a_re(i, j, k, 4)
        a_im(i, j + 1, k, 3) = a_im(i, j + 1, k, 3) - tmp_c*a_im(i, j, k, 4)
        a_re(i, j + 1, k, 4) = a_re(i, j + 1, k, 4) - tmp_r*a_re(i, j, k, 5)
        a_im(i, j + 1, k, 4) = a_im(i, j + 1, k, 4) - tmp_c*a_im(i, j, k, 5)

        ! eliminate diag-2
        tmp_r = 0._dp
        if (abs(a_re(i, j, k, 3)) > epsilon) then
          tmp_r = a_re(i, j + 2, k, 1)/a_re(i, j, k, 3)
        end if
        tmp_c = 0._dp
        if (abs(a_im(i, j, k, 3)) > epsilon) then
          tmp_c = a_im(i, j + 2, k, 1)/a_im(i, j, k, 3)
        end if
        div_r = real(div_u(i, jm + 2*inc, k) - tmp_r*div_u(i, jm, k), kind=dp)
        div_c = aimag(div_u(i, jm + 2*inc, k) - tmp_c*div_u(i, jm, k))
        div_u(i, jm + 2*inc, k) = cmplx(div_r, div_c, kind=dp)
        ! modify pentadiagonal coefficients in-place
        a_re(i, j + 2, k, 2) = a_re(i, j + 2, k, 2) - tmp_r*a_re(i, j, k, 4)
        a_im(i, j + 2, k, 2) = a_im(i, j + 2, k, 2) - tmp_c*a_im(i, j, k, 4)
        a_re(i, j + 2, k, 3) = a_re(i, j + 2, k, 3) - tmp_r*a_re(i, j, k, 5)
        a_im(i, j + 2, k, 3) = a_im(i, j + 2, k, 3) - tmp_c*a_im(i, j, k, 5)
      end do

      ! handle the last row
      if (abs(a_re(i, n - 1, k, 3)) > epsilon) then
        tmp_r = a_re(i, n, k, 2)/a_re(i, n - 1, k, 3)
      else
        tmp_r = 0._dp
      end if
      if (abs(a_im(i, n - 1, k, 3)) > epsilon) then
        tmp_c = a_im(i, n, k, 2)/a_im(i, n - 1, k, 3)
      else
        tmp_c = 0._dp
      end if
      div_r = a_re(i, n, k, 3) - tmp_r*a_re(i, n - 1, k, 4)
      div_c = a_im(i, n, k, 3) - tmp_c*a_im(i, n - 1, k, 4)

      ! j mapping based on odd/even for last point j=n
      nm = inc*n + off - inc/2
      if (abs(div_r) > epsilon) then
        tmp_r = tmp_r/div_r
        div_r = real(div_u(i, nm, k), kind=dp)/div_r &
                - tmp_r*real(div_u(i, nm - inc, k), kind=dp)
      else
        tmp_r = 0._dp
        div_r = 0._dp
      end if
      if (abs(div_c) > epsilon) then
        tmp_c = tmp_c/div_c
        div_c = aimag(div_u(i, nm, k))/div_c &
                - tmp_c*aimag(div_u(i, nm - inc, k))
      else
        tmp_c = 0._dp
        div_c = 0._dp
      end if
      div_u(i, nm, k) = cmplx(div_r, div_c, kind=dp)

      if (abs(a_re(i, n - 1, k, 3)) > epsilon) then
        tmp_r = 1._dp/a_re(i, n - 1, k, 3)
      else
        tmp_r = 0._dp
      end if
      if (abs(a_im(i, n - 1, k, 3)) > epsilon) then
        tmp_c = 1._dp/a_im(i, n - 1, k, 3)
      else
        tmp_c = 0._dp
      end if
      div_r = a_re(i, n - 1, k, 4)*tmp_r
      div_c = a_im(i, n - 1, k, 4)*tmp_c
      div_u(i, nm - inc, k) = cmplx( & !&
        real(div_u(i, nm - inc, k), kind=dp)*tmp_r &
        - real(div_u(i, nm, k), kind=dp)*div_r, &
        aimag(div_u(i, nm - inc, k))*tmp_c &
        - aimag(div_u(i, nm, k))*div_c, &
        kind=dp &
        )

      if (i == nx/2 + 1 .and. k == nz/2 + 1) then
        div_u(i, nm, k) = 0._dp
        div_u(i, nm - inc, k) = 0._dp
      end if

      ! backward pass
      do j = n - 2, 1, -1
        ! j mapping based on odd/even
        jm = inc*j + off - inc/2
        if (abs(a_re(i, j, k, 3)) > epsilon) then
          tmp_r = 1._dp/a_re(i, j, k, 3)
        else
          tmp_r = 0._dp
        end if
        if (abs(a_im(i, j, k, 3)) > epsilon) then
          tmp_c = 1._dp/a_im(i, j, k, 3)
        else
          tmp_c = 0._dp
        end if
        div_u(i, jm, k) = cmplx( & !&
          tmp_r*(real(div_u(i, jm, k), kind=dp) &
                 - a_re(i, j, k, 4)*real(div_u(i, jm + inc, k), kind=dp) &
                 - a_re(i, j, k, 5)*real(div_u(i, jm + 2*inc, k), kind=dp)), &
          tmp_c*(aimag(div_u(i, jm, k)) &
                 - a_im(i, j, k, 4)*aimag(div_u(i, jm + inc, k)) &
                 - a_im(i, j, k, 5)*aimag(div_u(i, jm + 2*inc, k))), &
          kind=dp &
          )
        if (i == nx/2 + 1 .and. k == nz/2 + 1) div_u(i, jm, k) = 0._dp
      end do
    end if

  end subroutine process_spectral_010_poisson

  attributes(global) subroutine process_spectral_010_bw( &
    div_u, nx_spec, ny_spec, y_sp_st, nx, ny, nz, ax, bx, ay, by, az, bz &
    )
    !! Backward post-processing only for non-periodic Y-direction (010).
    !!
    !! Performs odd/even mode recombination and post-processing in X and Z directions.
    !! Completes the spectral-to-physical transformation after Poisson solve.
    implicit none

    complex(dp), device, intent(inout), dimension(:, :, :) :: div_u  !! In/out: Solution field to post-process
    real(dp), device, intent(in), dimension(:) :: ax, bx, ay, by, az, bz  !! Input: Spectral equivalence constants
    !> Grid size in spectral space
    integer, value, intent(in) :: nx_spec, ny_spec  !! Spectral space grid size
    !> Offset in y direction in the permuted slabs in spectral space
    integer, value, intent(in) :: y_sp_st  !! Y-direction offset in spectral slabs
    !> Grid size
    integer, value, intent(in) :: nx, ny, nz  !! Physical space grid size

    integer :: i, j, k, ix, iy, iz, iy_rev  !! Loop, spectral, and reversed mode indices
    real(dp) :: tmp_r, tmp_c, div_r, div_c, l_r, l_c, r_r, r_c  !! Temporary real/imaginary components

    i = threadIdx%x + (blockIdx%x - 1)*blockDim%x
    k = blockIdx%y ! nz_spec

    ! post-process backward
    if (i <= nx_spec) then
      do j = 2, ny_spec/2 + 1
        ix = i; iy = j + y_sp_st; iz = k
        iy_rev = ny_spec - j + 2 + y_sp_st

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
          r_r*by(iy_rev) - r_c*ay(iy_rev) + l_r*ay(iy_rev) + l_c*by(iy_rev), &
          r_r*ay(iy_rev) + r_c*by(iy_rev) - l_r*by(iy_rev) + l_c*ay(iy_rev), &
          kind=dp &
          )
      end do
    end if

    if (i <= nx_spec) then
      do j = 1, ny_spec
        ix = i; iy = j + y_sp_st; iz = k

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
    end if

  end subroutine process_spectral_010_bw

  attributes(global) subroutine enforce_periodicity_y(f_out, f_in, ny)
    !! Enforce Y-direction periodicity by reordering data for non-periodic transforms.
    !!
    !! Maps full domain [1:ny] to symmetric layout required by non-periodic FFT.
    !! First half: odd points, second half: even points in reverse order.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: f_out  !! Output: Reordered field
    real(dp), device, intent(in), dimension(:, :, :) :: f_in  !! Input: Original field
    integer, value, intent(in) :: ny  !! Y-dimension size

    integer :: i, j, k  !! Thread and loop indices

    i = threadIdx%x
    k = blockIdx%x

    do j = 1, ny/2
      f_out(i, j, k) = f_in(i, 2*j - 1, k)
    end do
    do j = ny/2 + 1, ny
      f_out(i, j, k) = f_in(i, 2*ny - 2*j + 2, k)
    end do

  end subroutine enforce_periodicity_y

  attributes(global) subroutine undo_periodicity_y(f_out, f_in, ny)
    !! Undo Y-direction periodicity reordering after non-periodic transforms.
    !!
    !! Inverse of enforce_periodicity_y: reconstructs original domain layout
    !! from symmetric FFT ordering. Restores odd/even point positions.
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: f_out  !! Output: Restored field
    real(dp), device, intent(in), dimension(:, :, :) :: f_in  !! Input: Reordered field
    integer, value, intent(in) :: ny  !! Y-dimension size

    integer :: i, j, k  !! Thread and loop indices

    i = threadIdx%x
    k = blockIdx%x

    do j = 1, ny/2
      f_out(i, 2*j - 1, k) = f_in(i, j, k)
      f_out(i, 2*j, k) = f_in(i, ny - j + 1, k)
    end do

  end subroutine undo_periodicity_y

end module m_cuda_spectral
