module m_omp_poisson_fft

  use decomp_2d_constants, only: PHYSICAL_IN_X
  use decomp_2d_fft, only: decomp_2d_fft_init, decomp_2d_fft_3d, &
                           decomp_2d_fft_get_size

  use m_common, only: dp
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t

  use m_omp_spectral, only: process_spectral_000, process_spectral_010

  implicit none

  type, extends(poisson_fft_t) :: omp_poisson_fft_t
      !! FFT based Poisson solver
    complex(dp), allocatable, dimension(:, :, :) :: c_x, c_y, c_z
  contains
    procedure :: fft_forward => fft_forward_omp
    procedure :: fft_backward => fft_backward_omp
    procedure :: fft_postprocess_000 => fft_postprocess_000_omp
    procedure :: fft_postprocess_010 => fft_postprocess_010_omp
    procedure :: enforce_periodicity_y => enforce_periodicity_y_omp
    procedure :: undo_periodicity_y => undo_periodicity_y_omp
  end type omp_poisson_fft_t

  interface omp_poisson_fft_t
    module procedure init
  end interface omp_poisson_fft_t

  private :: init

contains

  function init(mesh, xdirps, ydirps, zdirps) result(poisson_fft)
    implicit none

    type(mesh_t), intent(in) :: mesh
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    integer, dimension(3) :: istart, iend, isize

    type(omp_poisson_fft_t) :: poisson_fft

    if (mesh%par%is_root()) then
      print *, "Initialising 2decomp&fft"
    end if

    ! Work out the spectral dimensions in the permuted state
    call decomp_2d_fft_init(PHYSICAL_IN_X)
    call decomp_2d_fft_get_size(istart, iend, isize)
    ! Converts a start position into an offset
    istart(:) = istart(:) - 1

    call poisson_fft%base_init(mesh, xdirps, ydirps, zdirps, isize, istart)

    allocate (poisson_fft%c_x(poisson_fft%nx_spec, poisson_fft%ny_spec, &
                              poisson_fft%nz_spec))

  end function init

  subroutine fft_forward_omp(self, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    call decomp_2d_fft_3d(f_in%data, self%c_x)

  end subroutine fft_forward_omp

  subroutine fft_backward_omp(self, f_out)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out

    call decomp_2d_fft_3d(self%c_x, f_out%data)

  end subroutine fft_backward_omp

  subroutine fft_postprocess_000_omp(self)
    implicit none

    class(omp_poisson_fft_t) :: self

    call process_spectral_000( &
      self%c_x, self%waves, self%nx_spec, self%ny_spec, self%nz_spec, &
      self%x_sp_st, self%y_sp_st, self%z_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax, self%bx, self%ay, self%by, self%az, self%bz &
      )

  end subroutine fft_postprocess_000_omp

  subroutine fft_postprocess_010_omp(self)
    implicit none

    class(omp_poisson_fft_t) :: self

    call process_spectral_010( &
      self%c_x, self%waves, self%nx_spec, self%ny_spec, self%nz_spec, &
      self%x_sp_st, self%y_sp_st, self%z_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax, self%bx, self%ay, self%by, self%az, self%bz &
      )

  end subroutine fft_postprocess_010_omp

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

end module m_omp_poisson_fft
