module m_test_memory_sizing
  !! Byte-exact sizes for the cuFFT/cuFFTMp spectral-array terms that change
  !! with GPU count, mirroring the allocate statements in
  !! src/backend/cuda/poisson_fft.f90. Deliberately narrow: only the terms
  !! that vary with ng and are needed to correct the memory-estimate
  !! probe's per-GPU table are covered here. Terms that do not vary with ng
  !! (the y-stretching arrays, the 110 case's extra real workspace) are
  !! already exact in the probe's ng=1 measurement and cancel out of the
  !! ng>1 delta, so they are deliberately not modelled here - see
  !! docs/memory_model_journey.md for the full methodology.
  use m_common, only: i8, nbytes
  implicit none

  private
  public :: cell_dims, spectral_slab_bytes, mirror_buffer_bytes_100

contains

  pure function cell_dims(vert_dims, periodic) result(c)
    !! Global CELL dims from global VERT dims. Mirrors src/mesh.f90:90-100:
    !! a periodic direction has as many cells as vertices, a non-periodic
    !! direction has one fewer.
    integer, intent(in) :: vert_dims(3)
    logical, intent(in) :: periodic(3)
    integer :: c(3)
    integer :: dir

    do dir = 1, 3
      if (periodic(dir)) then
        c(dir) = vert_dims(dir)
      else
        c(dir) = vert_dims(dir) - 1
      end if
    end do
  end function cell_dims

  pure function spectral_slab_bytes(bc_is_100, bc_is_110, cdims, ng) &
    result(nbytes8)
    !! Bytes of ONE complex spectral array (waves_dev, and c_dev on the
    !! plain-cuFFT fallback path), shape per
    !! src/backend/cuda/poisson_fft.f90:259-291:
    !!   100: (cy/2+1, cx/ng, cz)  - decomposition axis maps to dim2
    !!   110: (cz/2+1, cx, cy)    - never runs at ng>1 in this solver
    !!   000/010: (cx/2+1, cy/ng, cz) - decomposition axis maps to dim2
    !! ng is the number of ranks the z-direction is split over
    !! (mesh%par%nproc_dir(3)); note this differs from the physical
    !! decomposition axis (z), which is why ng divides cx or cy here, not
    !! cz.
    logical, intent(in) :: bc_is_100, bc_is_110
    integer, intent(in) :: cdims(3)
    integer, intent(in) :: ng
    integer(i8) :: nbytes8
    integer(i8) :: nx_spec, ny_spec, nz_spec

    if (bc_is_100) then
      nx_spec = cdims(2)/2 + 1
      ny_spec = cdims(1)/ng
      nz_spec = cdims(3)
    else if (bc_is_110) then
      nx_spec = cdims(3)/2 + 1
      ny_spec = cdims(1)
      nz_spec = cdims(2)
    else
      ! 000, 010
      nx_spec = cdims(1)/2 + 1
      ny_spec = cdims(2)/ng
      nz_spec = cdims(3)
    end if
    nbytes8 = nx_spec*ny_spec*nz_spec*2_i8*int(nbytes, i8)
  end function spectral_slab_bytes

  pure function mirror_buffer_bytes_100(cdims, ng) result(nbytes8)
    !! 100-case multi-rank-only mirror/exchange buffers, zero at ng<=1 -
    !! matches the nproc>1 guard at
    !! src/backend/cuda/poisson_fft.f90:332-349:
    !!   c_mirror_dev + c_slab_send_dev: (nx_spec, ny_spec, nz_spec) x2
    !!   c_plane_send_dev + c_plane_recv_dev: (nx_spec, nz_spec) x2
    !! all complex(dp).
    integer, intent(in) :: cdims(3)
    integer, intent(in) :: ng
    integer(i8) :: nbytes8
    integer(i8) :: nx_spec, ny_spec, nz_spec

    if (ng <= 1) then
      nbytes8 = 0_i8
      return
    end if
    nx_spec = cdims(2)/2 + 1
    ny_spec = cdims(1)/ng
    nz_spec = cdims(3)
    nbytes8 = 2_i8*nx_spec*ny_spec*nz_spec*2_i8*int(nbytes, i8) + &
              2_i8*nx_spec*nz_spec*2_i8*int(nbytes, i8)
  end function mirror_buffer_bytes_100

end module m_test_memory_sizing
