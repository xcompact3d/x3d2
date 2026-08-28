module m_omp_poisson_fft

  use decomp_2d_constants, only: PHYSICAL_IN_X
  use decomp_2d, only: decomp_info, get_decomp_dims, &
                       transpose_x_to_y, transpose_y_to_x, &
                       transpose_z_to_y, transpose_y_to_z
  use decomp_2d_fft, only: decomp_2d_fft_init, decomp_2d_fft_3d, &
                           decomp_2d_fft_get_ph, decomp_2d_fft_get_sp

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
    !> Non-periodic in x, periodic in y and z
    logical :: is_100_case = .false.
    !> 2decomp processor grid. p_row splits y in the x-pencil, p_col splits z.
    integer :: p_row = 1, p_col = 1
    !> Physical space decomposition owned by the FFT. The 100 case builds
    !! the transform on the swapped cell grid (ny, nx, nz), so ph%xsz is
    !! exactly the layout the transform reads, and ph%ysz is what a local
    !! dim1/dim2 swap of the untransposed x-pencil lands on. decomp_main
    !! cannot serve here because it is built on vertex dims.
    type(decomp_info), pointer :: ph => null()
    !> Y-pencil staging buffer of ph, only needed when y is decomposed.
    real(dp), allocatable, dimension(:, :, :) :: r_yp
    !> Spectral decomposition, and its y-pencil staging buffers. The FFT
    !! leaves the spectral slab in a z-pencil, where dim2 (the x modes) is
    !! split across p_col. The paired split in the postprocess couples a
    !! mode with its mirror in that same direction, so it needs dim2 whole.
    !! The y-pencil of the very same decomposition has dim2 complete for any
    !! processor grid, so 2decomp's own transpose delivers exactly the
    !! layout the pairing wants. waves is redistributed once, at init.
    !! sp is associated for every case, since init also reads the spectral
    !! slab extents from it, but c_pair and waves_pair exist only when the
    !! pairing actually needs the hop.
    type(decomp_info), pointer :: sp => null()
    complex(dp), allocatable, dimension(:, :, :) :: c_pair, waves_pair
    !> Exact sized, x-y transposed real work buffer for the 100 case.
    !! The FFT is initialised with the x and y dimensions swapped so that
    !! the r2c transform runs along the periodic y direction and the
    !! non-periodic x direction keeps its full mode range in dim2, which
    !! is what the paired split in the spectral postprocess needs. The
    !! buffer is also exactly (ny, nx, nz) rather than padded up to a
    !! multiple of SZ, because 2decomp plans its FFTW transforms over the
    !! flat array and a padded one would be read with the wrong stride.
    real(dp), allocatable, dimension(:, :, :) :: r_tr
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

    type(mesh_t), target, intent(in) :: mesh
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, optional, intent(in) :: lowmem
    integer :: dims(3), grid_dims(2)

    type(omp_poisson_fft_t) :: poisson_fft

    if (mesh%par%is_root()) then
      print *, "Initialising 2decomp&fft"
    end if

    if (present(lowmem)) then
      print *, 'lowmem_fft has no impact in the OpenMP backend.'
    end if

    ! Get global cell dims
    dims = mesh%get_global_dims(CELL)

    poisson_fft%is_100_case = (.not. mesh%grid%periodic_BC(1)) &
                              .and. mesh%grid%periodic_BC(2) &
                              .and. mesh%grid%periodic_BC(3)

    ! Work out the spectral dimensions in the permuted state.
    ! The 100 case initialises the transform with x and y swapped, so the
    ! spectral slab comes back as (ny/2 + 1, nx, nz). That is the layout
    ! waves_set builds for this case, so the base class needs no changes.
    if (poisson_fft%is_100_case) then
      ! get_decomp_dims returns the 2decomp processor grid (p_row, p_col).
      grid_dims = get_decomp_dims()
      poisson_fft%p_row = grid_dims(1)
      poisson_fft%p_col = grid_dims(2)
      call decomp_2d_fft_init(PHYSICAL_IN_X, dims(2), dims(1), dims(3))
    else
      call decomp_2d_fft_init(PHYSICAL_IN_X, dims(1), dims(2), dims(3))
    end if
    ! For PHYSICAL_IN_X the spectral slab is a z-pencil of sp, which is all
    ! decomp_2d_fft_get_size hands back. Read sp directly, subtracting one
    ! to convert the start position into an offset.
    poisson_fft%sp => decomp_2d_fft_get_sp()
    call poisson_fft%base_init(mesh, xdirps, ydirps, zdirps, &
                               poisson_fft%sp%zsz, poisson_fft%sp%zst - 1)

    if (mesh%geo%stretched(2)) then
      error stop 'OpenMP backends FFT based Poisson solver does not support&
                  & stretching in y-direction yet!'
    end if

    allocate (poisson_fft%c_x(poisson_fft%nx_spec, poisson_fft%ny_spec, &
                              poisson_fft%nz_spec))

    if (poisson_fft%is_100_case) then
      poisson_fft%ph => decomp_2d_fft_get_ph()
      if (poisson_fft%p_row > 1) then
        ! Reaching the transposed x-pencil takes a real redistribution once
        ! y is split: local dim1/dim2 swap -> y-pencil of ph -> x-pencil.
        allocate (poisson_fft%r_yp(poisson_fft%ph%ysz(1), &
                                   poisson_fft%ph%ysz(2), &
                                   poisson_fft%ph%ysz(3)))
        allocate (poisson_fft%r_tr(poisson_fft%ph%xsz(1), &
                                   poisson_fft%ph%xsz(2), &
                                   poisson_fft%ph%xsz(3)))
      else
        ! p_row == 1 makes the y-pencil identical to the x-pencil, so the
        ! swap alone reaches the transposed x-pencil with no communication.
        allocate (poisson_fft%r_tr(poisson_fft%ny_loc, poisson_fft%nx_loc, &
                                   poisson_fft%nz_loc))
      end if

      if (poisson_fft%p_col > 1) then
        allocate (poisson_fft%c_pair(poisson_fft%sp%ysz(1), &
                                  poisson_fft%sp%ysz(2), &
                                  poisson_fft%sp%ysz(3)))
        allocate (poisson_fft%waves_pair(poisson_fft%sp%ysz(1), &
                                      poisson_fft%sp%ysz(2), &
                                      poisson_fft%sp%ysz(3)))
        call transpose_z_to_y(poisson_fft%waves, poisson_fft%waves_pair, &
                              poisson_fft%sp)
      end if
    end if

  end function init

  subroutine fft_forward_omp(self, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    call decomp_2d_fft_3d(f_in%data, self%c_x)

  end subroutine fft_forward_omp

  subroutine fft_forward_010_omp(self, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    error stop 'OpenMP backend does not support fft_forward_010 yet!'

  end subroutine fft_forward_010_omp

  subroutine fft_forward_100_omp(self, f_in)
    !! Forward FFT for the non-periodic x case.
    !! Transposes x and y so that the non-periodic direction sits in dim2,
    !! where it keeps its full mode range, then transforms the (ny, nx, nz)
    !! buffer with the swapped plan set up in init.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    integer :: i, j, k, n1, n2, n3

    if (self%p_row > 1) then
      ! Strip the SZ padding and swap dim1 and dim2 in one pass, which
      ! lands on the y-pencil of ph, then redistribute to its x-pencil.
      n1 = self%ph%ysz(1); n2 = self%ph%ysz(2)
      n3 = self%ph%ysz(3)
      !$omp parallel do collapse(2)
      do k = 1, n3
        do i = 1, n2
          do j = 1, n1
            self%r_yp(j, i, k) = f_in%data(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do

      call transpose_y_to_x(self%r_yp, self%r_tr, self%ph)
    else
      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do i = 1, self%nx_loc
          do j = 1, self%ny_loc
            self%r_tr(j, i, k) = f_in%data(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do
    end if

    call decomp_2d_fft_3d(self%r_tr, self%c_x)

  end subroutine fft_forward_100_omp

  subroutine fft_forward_110_omp(self, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    error stop 'OpenMP backend does not support fft_forward_110 yet!'

  end subroutine fft_forward_110_omp

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
    !! Backward FFT for the non-periodic x case, undoing the x-y transpose
    !! applied by fft_forward_100_omp. Only the unpadded extents are
    !! written, which is all undo_periodicity_x_omp goes on to read.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out

    integer :: i, j, k, n1, n2, n3

    call decomp_2d_fft_3d(self%c_x, self%r_tr)

    if (self%p_row > 1) then
      call transpose_x_to_y(self%r_tr, self%r_yp, self%ph)

      n1 = self%ph%ysz(1); n2 = self%ph%ysz(2)
      n3 = self%ph%ysz(3)
      !$omp parallel do collapse(2)
      do k = 1, n3
        do j = 1, n1
          do i = 1, n2
            f_out%data(i, j, k) = self%r_yp(j, i, k)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do j = 1, self%ny_loc
          do i = 1, self%nx_loc
            f_out%data(i, j, k) = self%r_tr(j, i, k)
          end do
        end do
      end do
      !$omp end parallel do
    end if

  end subroutine fft_backward_100_omp

  subroutine fft_backward_110_omp(self, f_out)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out

    error stop 'OpenMP backend does not support fft_backward_110 yet!'

  end subroutine fft_backward_110_omp

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
    !! Post-process div U* in spectral space for the non-periodic x case.
    !!
    !! After the x-y transpose applied by fft_forward_100_omp, dim1 holds
    !! the periodic y r2c modes and dim2 holds the non-periodic x modes.
    !! That is the same arrangement process_spectral_010 already handles,
    !! with x and y swapped, so this reuses that routine rather than
    !! adding a second copy of the same maths. The swap is:
    !!   global sizes nx <-> ny
    !!   dim1 coefficients ax, bx -> ay, by
    !!   dim2 coefficients ay, by -> ax, bx  (the paired direction)
    !!
    !! The CUDA backend stages the equivalent work across three kernels so
    !! that a mode can fetch its pairing partner from another GPU. On a
    !! single rank every partner is local, so the staging collapses back
    !! into this one fused pass.
    implicit none

    class(omp_poisson_fft_t) :: self

    if (self%p_col > 1) then
      ! Hop to the y-pencil, where dim2 is whole, pair there, and hop back.
      call transpose_z_to_y(self%c_x, self%c_pair, self%sp)
      call process_spectral_010( &
        self%c_pair, self%waves_pair, self%sp%ysz(1), self%sp%ysz(2), &
        self%sp%ysz(3), self%sp%yst(1) - 1, self%sp%yst(2) - 1, &
        self%sp%yst(3) - 1, &
        self%ny_glob, self%nx_glob, self%nz_glob, &
        self%ay, self%by, self%ax, self%bx, self%az, self%bz &
        )
      call transpose_y_to_z(self%c_pair, self%c_x, self%sp)
    else
      call process_spectral_010( &
        self%c_x, self%waves, self%nx_spec, self%ny_spec, self%nz_spec, &
        self%sp_st(1), self%sp_st(2), self%sp_st(3), &
        self%ny_glob, self%nx_glob, self%nz_glob, &
        self%ay, self%by, self%ax, self%bx, self%az, self%bz &
        )
    end if

  end subroutine fft_postprocess_100_omp

  subroutine fft_postprocess_110_omp(self)
    implicit none

    class(omp_poisson_fft_t) :: self

    error stop 'OpenMP backend does not support fft_postprocess_110 yet!'

  end subroutine fft_postprocess_110_omp

  subroutine enforce_periodicity_x_omp(self, f_out, f_in)
    !! Gathers the non-periodic x line into a periodic one so that a plain
    !! FFT can stand in for the cosine transform. Unlike the y version
    !! below this handles an odd nx, which the 100 case requires.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k, n2

    n2 = self%nx_glob/2

    !$omp parallel do collapse(2)
    do k = 1, self%nz_loc
      do j = 1, self%ny_loc
        do i = 1, n2
          f_out%data(i, j, k) = f_in%data(2*i - 1, j, k)
        end do
        if (mod(self%nx_glob, 2) == 1) then
          ! odd-size centre entry
          f_out%data(n2 + 1, j, k) = f_in%data(self%nx_glob, j, k)
          do i = n2 + 2, self%nx_glob
            f_out%data(i, j, k) = f_in%data(2*self%nx_glob - 2*i + 2, j, k)
          end do
        else
          do i = n2 + 1, self%nx_glob
            f_out%data(i, j, k) = f_in%data(2*self%nx_glob - 2*i + 2, j, k)
          end do
        end if
      end do
    end do
    !$omp end parallel do

  end subroutine enforce_periodicity_x_omp

  subroutine undo_periodicity_x_omp(self, f_out, f_in)
    !! Scatters the gathered x line back to its original ordering.
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k, n2

    n2 = self%nx_glob/2

    ! The gathered ordering only covers 1..nx_glob, but the block is padded
    ! up to a multiple of SZ and reorder copies the whole padded extent back
    ! out to the pressure field, so the untouched entries have to be cleared
    ! rather than left holding whatever the recycled block came with.
    f_out%data = 0._dp

    !$omp parallel do collapse(2)
    do k = 1, self%nz_loc
      do j = 1, self%ny_loc
        do i = 1, n2
          f_out%data(2*i - 1, j, k) = f_in%data(i, j, k)
          f_out%data(2*i, j, k) = f_in%data(self%nx_glob - i + 1, j, k)
        end do
        if (mod(self%nx_glob, 2) == 1) then
          ! odd-size centre entry
          f_out%data(self%nx_glob, j, k) = f_in%data(n2 + 1, j, k)
        end if
      end do
    end do
    !$omp end parallel do

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
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    error stop 'OpenMP backend does not support enforce_periodicity_xy yet!'

  end subroutine enforce_periodicity_xy_omp

  subroutine undo_periodicity_xy_omp(self, f_out, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    error stop 'OpenMP backend does not support undo_periodicity_xy yet!'

  end subroutine undo_periodicity_xy_omp
end module m_omp_poisson_fft