module m_omp_poisson_fft

  use decomp_2d_constants, only: PHYSICAL_IN_X, PHYSICAL_IN_Z
  use decomp_2d_fft, only: decomp_2d_fft_init, decomp_2d_fft_3d, &
                           decomp_2d_fft_get_size

  use m_common, only: dp, CELL
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t

  use m_omp_spectral, only: process_spectral_000, process_spectral_010

  implicit none

  type, extends(poisson_fft_t) :: omp_poisson_fft_t
      !! FFT based Poisson solver
    complex(dp), allocatable, dimension(:, :, :) :: c_x, c_y, c_z
    ! Contiguous real workspace for the 110 FFT path. The fields handed
    ! to the Poisson solver come from the OMP allocator with SZ-aware
    ! padding (e.g. logical 256x128, stride 272x144), but 2decomp's
    ! decomp_2d_fft_3d expects a contiguous logical-sized array. The
    ! 110 wrappers copy through r_x. 000/010 happen to work with the
    ! padded buffer directly when the case dims don't trigger padding;
    ! that is a pre-existing latent issue not addressed here.
    real(dp), allocatable, dimension(:, :, :) :: r_x
  contains
    procedure :: fft_forward => fft_forward_omp
    procedure :: fft_forward_010 => fft_forward_omp
    procedure :: fft_forward_100 => fft_forward_100_omp
    procedure :: fft_forward_110 => fft_forward_110_omp
    procedure :: fft_backward => fft_backward_omp
    procedure :: fft_backward_010 => fft_backward_omp
    procedure :: fft_backward_100 => fft_backward_100_omp
    procedure :: fft_backward_110 => fft_backward_110_omp
    procedure :: fft_postprocess_000 => fft_postprocess_000_omp
    procedure :: fft_postprocess_010 => fft_postprocess_010_omp
    procedure :: fft_postprocess_100 => fft_postprocess_100_omp
    procedure :: fft_postprocess_110 => fft_postprocess_110_omp
    procedure :: enforce_periodicity_x => enforce_periodicity_x_omp
    procedure :: undo_periodicity_x => undo_periodicity_x_omp
    procedure :: enforce_periodicity_y => enforce_periodicity_y_omp
    procedure :: undo_periodicity_y => undo_periodicity_y_omp
    procedure :: enforce_periodicity_xy => enforce_periodicity_xy_omp
    procedure :: undo_periodicity_xy => undo_periodicity_xy_omp
  end type omp_poisson_fft_t

  interface omp_poisson_fft_t
    module procedure init
  end interface omp_poisson_fft_t

  private :: init

contains

  function init(mesh, xdirps, ydirps, zdirps, lowmem) result(poisson_fft)
    implicit none

    type(mesh_t), intent(in) :: mesh
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, optional, intent(in) :: lowmem
    integer, dimension(3) :: istart, iend, isize
    integer :: dims(3)

    type(omp_poisson_fft_t) :: poisson_fft

    if (mesh%par%is_root()) then
      print *, "Initialising 2decomp&fft"
    end if

    if (present(lowmem)) then
      print *, 'lowmem_fft has no impact in the OpenMP backend.'
    end if

    ! Get global cell dims
    dims = mesh%get_global_dims(CELL)

    ! Pick the physical-input axis for the R2C transform.
    ! 2decomp folds the first physical axis into nx/2+1 (or nz/2+1) modes,
    ! so the fold has to land on a periodic axis. The 110 case (non-periodic
    ! X, non-periodic Y, periodic Z) is the only configuration in this
    ! backend that needs PHYSICAL_IN_Z; 000 and 010 are fine with the
    ! default PHYSICAL_IN_X because X is periodic in both.
    if ((.not. mesh%grid%periodic_BC(1)) .and. &
        (.not. mesh%grid%periodic_BC(2)) .and. &
        mesh%grid%periodic_BC(3)) then
      if (mesh%par%is_root()) then
        print *, "2decomp FFT engine: PHYSICAL_IN_Z (110 case)"
      end if
      call decomp_2d_fft_init(PHYSICAL_IN_Z, dims(1), dims(2), dims(3))
    else
      call decomp_2d_fft_init(PHYSICAL_IN_X, dims(1), dims(2), dims(3))
    end if

    call decomp_2d_fft_get_size(istart, iend, isize)
    ! Converts a start position into an offset
    istart(:) = istart(:) - 1

    call poisson_fft%base_init(mesh, xdirps, ydirps, zdirps, isize, istart)

    if (mesh%geo%stretched(2)) then
      error stop 'OpenMP backends FFT based Poisson solver does not support&
                  & stretching in y-direction yet!'
    end if

    allocate (poisson_fft%c_x(poisson_fft%nx_spec, poisson_fft%ny_spec, &
                              poisson_fft%nz_spec))
    allocate (poisson_fft%r_x(poisson_fft%nx_loc, poisson_fft%ny_loc, &
                              poisson_fft%nz_loc))

  end function init
subroutine fft_roundtrip_check(self, f_in)
    !! TODO: remove once 110 is verified end-to-end.
    !! Confirms that a forward + manual /N + backward FFT round-trips
    !! a real input to itself, within round-off. Run once at startup
    !! via a save-guarded call from fft_forward_omp. Rank-local only;
    !! non-periodic BC paths run single-rank by upstream construction.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    real(dp), allocatable :: r_tmp(:, :, :)
    complex(dp), allocatable :: c_tmp(:, :, :)
    real(dp) :: max_diff, norm, val
    integer :: i, j, k
    integer :: nx_r, ny_r, nz_r

    nx_r = size(f_in%data, 1)
    ny_r = size(f_in%data, 2)
    nz_r = size(f_in%data, 3)

    ! Diagnostics: how does the padded array shape compare to the
    ! logical grid? Are there large values sitting in the padded region?
    if (self%mesh%par%is_root()) then
      print '(A, 3I5)', ' f_in shape (padded): ', nx_r, ny_r, nz_r
      print '(A, 3I5)', ' f_in logical (glob): ', &
        self%nx_glob, self%ny_glob, self%nz_glob
      print '(A, ES12.5)', ' max |f_in| over full padded array: ', &
        maxval(abs(f_in%data))
      print '(A, ES12.5)', ' max |f_in| over logical region:    ', &
        maxval(abs(f_in%data(1:self%nx_glob, &
                             1:self%ny_glob, &
                             1:self%nz_glob)))
    end if

    allocate (r_tmp(nx_r, ny_r, nz_r))
    allocate (c_tmp(self%nx_spec, self%ny_spec, self%nz_spec))

    if (self%mesh%par%is_root()) then
      print '(A, 3I5)', ' 2decomp pencil dims (xstart..xend over CELL): ', &
        self%nx_loc, self%ny_loc, self%nz_loc
    end if

    if (self%mesh%par%is_root()) then
      print '(A, I0)', ' f_in%dir = ', f_in%dir
      print '(A, 3I5)', ' actual stride in memory: ', &
        size(f_in%data, 1), size(f_in%data, 2), size(f_in%data, 3)
    end if
    ! Forward: real -> complex
    call decomp_2d_fft_3d(f_in%data, c_tmp)

    ! Manual 2decomp normalisation (forward + backward both unscaled)
    norm = 1._dp/real(self%nx_glob, dp)/real(self%ny_glob, dp) &
           /real(self%nz_glob, dp)
    c_tmp = c_tmp*norm

    ! Backward: complex -> real
    call decomp_2d_fft_3d(c_tmp, r_tmp)

    ! Diff only over the logical region; the padded entries are
    ! uninitialised on entry and may not survive the FFT round-trip.
    max_diff = 0._dp
    do k = 1, self%nz_glob
      do j = 1, self%ny_glob
        do i = 1, self%nx_glob
          val = abs(r_tmp(i, j, k) - f_in%data(i, j, k))
          if (val > max_diff) max_diff = val
        end do
      end do
    end do

    if (self%mesh%par%is_root()) then
      print '(A, ES12.5)', ' FFT round-trip max abs diff: ', max_diff
    end if

    deallocate (r_tmp, c_tmp)

  end subroutine fft_roundtrip_check
  subroutine fft_forward_omp(self, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    ! TODO: remove once the 110 path is committed and verified end-to-end.
    logical, save :: did_check = .false.

    if (.not. did_check) then
      call fft_roundtrip_check(self, f_in)
      did_check = .true.
    end if

    call decomp_2d_fft_3d(f_in%data, self%c_x)

  end subroutine fft_forward_omp

  subroutine fft_forward_010_omp(self, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    error stop 'OpenMP backend does not support fft_forward_010 yet!'

  end subroutine fft_forward_010_omp

  subroutine fft_forward_100_omp(self, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    error stop 'OpenMP backend does not support fft_forward_100 yet!'

  end subroutine fft_forward_100_omp

  subroutine fft_forward_110_omp(self, f_in)
    !! Forward FFT for the 110 case (non-periodic X, non-periodic Y,
    !! periodic Z). The poisson_fft engine was initialised with
    !! PHYSICAL_IN_Z so 2decomp folds along Z; the spectral output
    !! lands in self%c_x with shape (nx, ny, nz/2+1).
    !!
    !! Mirrors fft_forward_110_cuda in shape: dedicated wrapper, not an
    !! alias of the general fft_forward. The copy through self%r_x is
    !! the OMP-side equivalent of the CUDA Z-transpose -- a layout
    !! adjustment that gets the data into the form 2decomp expects.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    ! TODO: remove once 110 is verified end-to-end.
    logical, save :: did_check = .false.

    if (.not. did_check) then
      call fft_roundtrip_check_110(self, f_in)
      did_check = .true.
    end if

    ! Copy the logical region out of the allocator-padded buffer into a
    ! contiguous workspace before handing to 2decomp.
    call copy_logical_to_contig(self%r_x, f_in%data, &
                                self%nx_loc, self%ny_loc, self%nz_loc)
    call decomp_2d_fft_3d(self%r_x, self%c_x)

  end subroutine fft_forward_110_omp

  subroutine fft_backward_110_omp(self, f_out)
    !! Backward FFT for the 110 case. Mirrors fft_backward_110_cuda:
    !! dedicated wrapper, copies back through self%r_x.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out

    call decomp_2d_fft_3d(self%c_x, self%r_x)
    ! Copy back from the contiguous workspace into the logical region of
    ! the allocator-padded destination buffer. Padded entries of
    ! f_out%data are left untouched.
    call copy_contig_to_logical(f_out%data, self%r_x, &
                                self%nx_loc, self%ny_loc, self%nz_loc)

  end subroutine fft_backward_110_omp

  subroutine copy_logical_to_contig(dst, src, nx, ny, nz)
    !! Copy the (nx, ny, nz) logical region out of an allocator-padded
    !! source array into a contiguous destination of the same logical size.
    implicit none
    real(dp), intent(out) :: dst(:, :, :)
    real(dp), intent(in) :: src(:, :, :)
    integer, intent(in) :: nx, ny, nz
    integer :: i, j, k

    !$omp parallel do collapse(2) private(i)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          dst(i, j, k) = src(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine copy_logical_to_contig

  subroutine copy_contig_to_logical(dst, src, nx, ny, nz)
    !! Copy the (nx, ny, nz) logical region from a contiguous source array
    !! back into the (possibly padded) destination buffer. Padded entries
    !! of dst are left untouched.
    implicit none
    real(dp), intent(inout) :: dst(:, :, :)
    real(dp), intent(in) :: src(:, :, :)
    integer, intent(in) :: nx, ny, nz
    integer :: i, j, k

    !$omp parallel do collapse(2) private(i)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          dst(i, j, k) = src(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine copy_contig_to_logical

  subroutine fft_roundtrip_check_110(self, f_in)
    !! TODO: remove once 110 is verified end-to-end.
    !! Confirms that a forward + manual /N + backward FFT round-trips
    !! a real input to itself, within round-off. Run once via a
    !! save-guarded call from fft_forward_110_omp. Rank-local only;
    !! non-periodic BC paths run single-rank by upstream construction.
    !! Exercises the same copy-through-r_x code path the production
    !! 110 wrappers take, so a clean result here means the production
    !! wrappers are also handling the layout correctly.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    real(dp), allocatable :: r_in(:, :, :), r_out(:, :, :)
    complex(dp), allocatable :: c_tmp(:, :, :)
    real(dp) :: max_diff, norm, val
    integer :: i, j, k

    allocate (r_in(self%nx_loc, self%ny_loc, self%nz_loc))
    allocate (r_out(self%nx_loc, self%ny_loc, self%nz_loc))
    allocate (c_tmp(self%nx_spec, self%ny_spec, self%nz_spec))

    call copy_logical_to_contig(r_in, f_in%data, &
                                self%nx_loc, self%ny_loc, self%nz_loc)

    call decomp_2d_fft_3d(r_in, c_tmp)
    norm = 1._dp/real(self%nx_glob, dp)/real(self%ny_glob, dp) &
           /real(self%nz_glob, dp)
    c_tmp = c_tmp*norm
    call decomp_2d_fft_3d(c_tmp, r_out)

    max_diff = 0._dp
    do k = 1, self%nz_loc
      do j = 1, self%ny_loc
        do i = 1, self%nx_loc
          val = abs(r_out(i, j, k) - r_in(i, j, k))
          if (val > max_diff) max_diff = val
        end do
      end do
    end do

    if (self%mesh%par%is_root()) then
      print '(A, ES12.5)', ' FFT round-trip max abs diff (110): ', max_diff
    end if

    deallocate (r_in, r_out, c_tmp)

  end subroutine fft_roundtrip_check_110

  subroutine fft_backward_omp(self, f_out)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out

    call decomp_2d_fft_3d(self%c_x, f_out%data)

  end subroutine fft_backward_omp

  subroutine fft_backward_010_omp(self, f_out)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out

    error stop 'OpenMP backend does not support fft_backward_010 yet!'

  end subroutine fft_backward_010_omp

  subroutine fft_backward_100_omp(self, f_out)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out

    error stop 'OpenMP backend does not support fft_backward_100 yet!'

  end subroutine fft_backward_100_omp


  subroutine fft_postprocess_000_omp(self)
    implicit none

    class(omp_poisson_fft_t) :: self

    call process_spectral_000( &
      self%c_x, self%waves, self%nx_spec, self%ny_spec, self%nz_spec, &
      self%sp_st(1), self%sp_st(2), self%sp_st(3), &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax, self%bx, self%ay, self%by, self%az, self%bz &
      )

  end subroutine fft_postprocess_000_omp

  subroutine fft_postprocess_010_omp(self)
    implicit none

    class(omp_poisson_fft_t) :: self

    call process_spectral_010( &
      self%c_x, self%waves, self%nx_spec, self%ny_spec, self%nz_spec, &
      self%sp_st(1), self%sp_st(2), self%sp_st(3), &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax, self%bx, self%ay, self%by, self%az, self%bz &
      )

  end subroutine fft_postprocess_010_omp

  subroutine fft_postprocess_100_omp(self)
    implicit none

    class(omp_poisson_fft_t) :: self

    error stop 'OpenMP backend does not support fft_postprocess_100 yet!'

  end subroutine fft_postprocess_100_omp

  subroutine fft_postprocess_110_omp(self)
    implicit none

    class(omp_poisson_fft_t) :: self

    error stop 'OpenMP backend does not support fft_postprocess_110 yet!'

  end subroutine fft_postprocess_110_omp

  subroutine enforce_periodicity_x_omp(self, f_out, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    error stop 'OpenMP backend does not support enforce_periodicity_x yet!'

  end subroutine enforce_periodicity_x_omp

  subroutine undo_periodicity_x_omp(self, f_out, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    error stop 'OpenMP backend does not support undo_periodicity_x yet!'

  end subroutine undo_periodicity_x_omp

  subroutine enforce_periodicity_y_omp(self, f_out, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k

    !$omp parallel do
    do k = 1, self%nz_loc
      do j = 1, self%ny_glob/2
        do i = 1, self%nx_loc
          f_out%data(i, j, k) = f_in%data(i, 2*(j - 1) + 1, k)
        end do
      end do
      do j = self%ny_glob/2 + 1, self%ny_glob
        do i = 1, self%nx_loc
          f_out%data(i, j, k) = f_in%data(i, 2*self%ny_glob - 2*j + 2, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine enforce_periodicity_y_omp

  subroutine undo_periodicity_y_omp(self, f_out, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k

    !$omp parallel do
    do k = 1, self%nz_loc
      do i = 1, self%nx_loc
        do j = 1, self%ny_glob/2
          f_out%data(i, 2*j - 1, k) = f_in%data(i, j, k)
        end do
        do j = 1, self%ny_glob/2
          f_out%data(i, 2*j, k) = f_in%data(i, self%ny_glob - j + 1, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine undo_periodicity_y_omp

subroutine enforce_periodicity_xy_omp(self, f_out, f_in)
    !! Combined X and Y periodicity enforcement (interleave shuffle).
    !! Port of the CUDA kernel: pure index reordering, no FFT.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k, n2x, n2y, src_i, src_j
    integer :: nx, ny

    nx = self%nx_glob
    ny = self%ny_glob
    n2x = nx/2
    n2y = ny/2

    !$omp parallel do collapse(2) private(i, j, src_i, src_j)
    do k = 1, self%nz_loc
      do j = 1, ny
        if (j <= n2y) then
          src_j = 2*j - 1
        else if (mod(ny, 2) == 1 .and. j == n2y + 1) then
          src_j = ny
        else
          src_j = 2*ny - 2*j + 2
        end if
        do i = 1, nx
          if (i <= n2x) then
            src_i = 2*i - 1
          else if (mod(nx, 2) == 1 .and. i == n2x + 1) then
            src_i = nx
          else
            src_i = 2*nx - 2*i + 2
          end if
          f_out%data(i, j, k) = f_in%data(src_i, src_j, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine enforce_periodicity_xy_omp

  subroutine undo_periodicity_xy_omp(self, f_out, f_in)
    !! Combined X and Y periodicity undo (reverse interleave shuffle).
    !! Port of the CUDA kernel: pure index reordering, no FFT.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k, n2x, n2y, src_i, src_j
    integer :: nx, ny

    nx = self%nx_glob
    ny = self%ny_glob
    n2x = nx/2
    n2y = ny/2

    !$omp parallel do collapse(2) private(i, j, src_i, src_j)
    do k = 1, self%nz_loc
      do j = 1, ny
        if (mod(ny, 2) == 1 .and. j == ny) then
          src_j = n2y + 1
        else if (mod(j, 2) == 1) then
          src_j = (j + 1)/2
        else
          src_j = ny - j/2 + 1
        end if
        do i = 1, nx
          if (mod(nx, 2) == 1 .and. i == nx) then
            src_i = n2x + 1
          else if (mod(i, 2) == 1) then
            src_i = (i + 1)/2
          else
            src_i = nx - i/2 + 1
          end if
          f_out%data(i, j, k) = f_in%data(src_i, src_j, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine undo_periodicity_xy_omp
end module m_omp_poisson_fft
