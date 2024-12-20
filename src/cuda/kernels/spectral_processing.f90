module m_cuda_spectral
  use cudafor

  use m_common, only: dp

  implicit none

contains

  attributes(global) subroutine process_spectral_000( &
    div_u, waves, nx_spec, ny_spec, y_sp_st, nx, ny, nz, &
    ax, bx, ay, by, az, bz &
    )
    !! Post-processes the divergence of velocity in spectral space, including
    !! scaling w.r.t. grid size.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    !> Divergence of velocity in spectral space
    complex(dp), device, intent(inout), dimension(:, :, :) :: div_u
    !> Spectral equivalence constants
    complex(dp), device, intent(in), dimension(:, :, :) :: waves
    real(dp), device, intent(in), dimension(:) :: ax, bx, ay, by, az, bz
    !> Grid size in spectral space
    integer, value, intent(in) :: nx_spec, ny_spec
    !> Offset in y direction in the permuted slabs in spectral space
    integer, value, intent(in) :: y_sp_st
    !> Grid size
    integer, value, intent(in) :: nx, ny, nz

    integer :: i, j, k, ix, iy, iz
    real(dp) :: tmp_r, tmp_c, div_r, div_c

    j = threadIdx%x + (blockIdx%x - 1)*blockDim%x
    k = blockIdx%y ! nz_spec

    if (j <= ny_spec) then
      do i = 1, nx_spec
        ! normalisation
        div_r = real(div_u(i, j, k), kind=dp)/(nx*ny*nz)
        div_c = aimag(div_u(i, j, k))/(nx*ny*nz)

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
    !! Post-processes the divergence of velocity in spectral space, including
    !! scaling w.r.t. grid size.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    !> Divergence of velocity in spectral space
    complex(dp), device, intent(inout), dimension(:, :, :) :: div_u
    !> Spectral equivalence constants
    complex(dp), device, intent(in), dimension(:, :, :) :: waves
    real(dp), device, intent(in), dimension(:) :: ax, bx, ay, by, az, bz
    !> Grid size in spectral space
    integer, value, intent(in) :: nx_spec, ny_spec
    !> Offset in y direction in the permuted slabs in spectral space
    integer, value, intent(in) :: y_sp_st
    !> Grid size
    integer, value, intent(in) :: nx, ny, nz

    integer :: i, j, k, ix, iy, iz
    real(dp) :: tmp_r, tmp_c, div_r, div_c

    i = threadIdx%x + (blockIdx%x - 1)*blockDim%x
    k = blockIdx%y ! nz_spec

    if (i <= nx_spec) then
      do j = 1, ny_spec
        ! normalisation
        div_r = real(div_u(i, j, k), kind=dp)/(nx*ny*nz)
        div_c = aimag(div_u(i, j, k))/(nx*ny*nz)

        ix = i; iy = j + y_sp_st; iz = k

        ! post-process forward
        ! post-process in y
        if (j /= 1) then
          tmp_r = div_r
          tmp_c = div_c
          div_r = (tmp_r*by(iy) + tmp_c*ay(iy))/2
          div_c = (tmp_c*by(iy) - tmp_r*ay(iy))/2

          tmp_r = real(div_u(i, ny_spec - j + 2, k), kind=dp)/(nx*ny*nz)
          tmp_c = aimag(div_u(i, ny_spec - j + 2, k))/(nx*ny*nz)
          div_r = div_r + (tmp_r*by(iy) - tmp_c*ay(iy))/2
          div_c = div_c + (tmp_c*by(iy) + tmp_r*ay(iy))/2
        end if

        ! post-process in z
        tmp_r = div_r
        tmp_c = div_c
        div_r = tmp_r*bz(iz) + tmp_c*az(iz)
        div_c = tmp_c*bz(iz) - tmp_r*az(iz)
        if (iz > nz/2 + 1) div_r = -div_r
        if (iz > nz/2 + 1) div_c = -div_c

        ! post-process in x
        tmp_r = div_r
        tmp_c = div_c
        div_r = tmp_r*bx(ix) + tmp_c*ax(ix)
        div_c = tmp_c*bx(ix) - tmp_r*ax(ix)

        ! Solve Poisson
        tmp_r = real(waves(i, j, k), kind=dp)
        tmp_c = aimag(waves(i, j, k))
        if (tmp_r < 1.e-16_dp) then
          div_r = 0._dp
        else
          div_r = -div_r/tmp_r
        end if
        if (tmp_c < 1.e-16_dp) then
          div_c = 0._dp
        else
          div_c = -div_c/tmp_c
        end if

        ! update the entry
        div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
      end do
    end if

    ! post-process backward
    if (i <= nx_spec) then
      do j = 1, ny_spec
        div_r = real(div_u(i, j, k), kind=dp)
        div_c = aimag(div_u(i, j, k))

        ix = i; iy = j + y_sp_st; iz = k

        ! post-process in y
        if (j /= 1) then
          tmp_r = div_r
          tmp_c = div_c
          div_r = tmp_r*by(iy) - tmp_c*ay(iy)
          div_c = tmp_c*by(iy) + tmp_r*ay(iy)

          tmp_r = real(div_u(i, ny_spec - j + 2, k), kind=dp)
          tmp_c = aimag(div_u(i, ny_spec - j + 2, k))
          div_r = div_r + tmp_r*ay(iy) + tmp_c*by(iy)
          div_c = div_c - tmp_c*ay(iy) + tmp_r*by(iy)
        end if

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

        ! update the entry
        div_u(i, j, k) = cmplx(div_r, div_c, kind=dp)
      end do
    end if

  end subroutine process_spectral_010

  attributes(global) subroutine enforce_periodicity_y(f_out, f_in, ny)
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: f_out
    real(dp), device, intent(in), dimension(:, :, :) :: f_in
    integer, value, intent(in) :: ny

    integer :: i, j, k

    i = threadIdx%x
    k = blockIdx%x

    do j = 1, ny/2
      f_out(i, j, k) = f_in(i, 2*j - 1, k)
      f_out(i, j + ny/2, k) = f_in(i, ny - 2*j + 2, k)
    end do

  end subroutine enforce_periodicity_y

  attributes(global) subroutine undo_periodicity_y(f_out, f_in, ny)
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: f_out
    real(dp), device, intent(in), dimension(:, :, :) :: f_in
    integer, value, intent(in) :: ny

    integer :: i, j, k

    i = threadIdx%x
    k = blockIdx%x

    do j = 1, ny/2
      f_out(i, 2*j - 1, k) = f_in(i, j, k)
      f_out(i, 2*j, k) = f_in(i, ny - j + 1, k)
    end do

  end subroutine undo_periodicity_y

end module m_cuda_spectral
