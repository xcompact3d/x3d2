module m_cuda_poisson_fft
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use cufftXt
  use cufft
  use mpi

  use m_common, only: dp, CELL
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t

  use m_cuda_allocator, only: cuda_field_t
  use m_cuda_spectral, only: process_spectral_000, process_spectral_010, &
                             memcpy3D

  implicit none

  type, extends(poisson_fft_t) :: cuda_poisson_fft_t
    !! FFT based Poisson solver

    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), device, allocatable, dimension(:, :, :) :: waves_dev
    !> Wave numbers in x, y, and z
    real(dp), device, allocatable, dimension(:) :: ax_dev, bx_dev, &
                                                   ay_dev, by_dev, &
                                                   az_dev, bz_dev
    !> Forward and backward FFT transform plans
    integer :: plan3D_fw, plan3D_bw

    !> cuFFTMp object manages decomposition and data storage
    type(cudaLibXtDesc), pointer :: xtdesc
  contains
    procedure :: fft_forward => fft_forward_cuda
    procedure :: fft_backward => fft_backward_cuda
    procedure :: fft_postprocess_000 => fft_postprocess_000_cuda
    procedure :: fft_postprocess_010 => fft_postprocess_010_cuda
    procedure :: enforce_periodicity_y => enforce_periodicity_y_cuda
    procedure :: undo_periodicity_y => undo_periodicity_y_cuda
  end type cuda_poisson_fft_t

  interface cuda_poisson_fft_t
    module procedure init
  end interface cuda_poisson_fft_t

  private :: init

contains

  function init(mesh, xdirps, ydirps, zdirps) result(poisson_fft)
    implicit none

    type(mesh_t), intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    type(cuda_poisson_fft_t) :: poisson_fft

    integer :: nx, ny, nz

    integer :: ierr
    integer(int_ptr_kind()) :: worksize

    integer :: dims_glob(3), dims_loc(3), n_spec(3), n_sp_st(3)

    ! 1D decomposition along Z in real domain, and along Y in spectral space
    if (mesh%par%nproc_dir(2) /= 1) print *, 'nproc_dir in y-dir must be 1'

    ! Work out the spectral dimensions in the permuted state
    dims_glob = mesh%get_global_dims(CELL)
    dims_loc = mesh%get_dims(CELL)

    n_spec(1) = dims_loc(1)/2 + 1
    n_spec(2) = dims_loc(2)/mesh%par%nproc_dir(3)
    n_spec(3) = dims_glob(3)

    n_sp_st(1) = 0
    n_sp_st(2) = dims_loc(2)/mesh%par%nproc_dir(3)*mesh%par%nrank_dir(3)
    n_sp_st(3) = 0

    call poisson_fft%base_init(mesh, xdirps, ydirps, zdirps, n_spec, n_sp_st)

    nx = poisson_fft%nx_glob
    ny = poisson_fft%ny_glob
    nz = poisson_fft%nz_glob

    allocate (poisson_fft%waves_dev(poisson_fft%nx_spec, &
                                    poisson_fft%ny_spec, &
                                    poisson_fft%nz_spec))
    poisson_fft%waves_dev = poisson_fft%waves

    allocate (poisson_fft%ax_dev(nx), poisson_fft%bx_dev(nx))
    allocate (poisson_fft%ay_dev(ny), poisson_fft%by_dev(ny))
    allocate (poisson_fft%az_dev(nz), poisson_fft%bz_dev(nz))
    poisson_fft%ax_dev = poisson_fft%ax; poisson_fft%bx_dev = poisson_fft%bx
    poisson_fft%ay_dev = poisson_fft%ay; poisson_fft%by_dev = poisson_fft%by
    poisson_fft%az_dev = poisson_fft%az; poisson_fft%bz_dev = poisson_fft%bz

    ! 3D plans
    ierr = cufftCreate(poisson_fft%plan3D_fw)
    ierr = cufftMpAttachComm(poisson_fft%plan3D_fw, CUFFT_COMM_MPI, &
                             MPI_COMM_WORLD)
    ierr = cufftMakePlan3D(poisson_fft%plan3D_fw, nz, ny, nx, CUFFT_D2Z, &
                           worksize)
    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Forward 3D FFT plan generation failed'
    end if

    ierr = cufftCreate(poisson_fft%plan3D_bw)
    ierr = cufftMpAttachComm(poisson_fft%plan3D_bw, CUFFT_COMM_MPI, &
                             MPI_COMM_WORLD)
    ierr = cufftMakePlan3D(poisson_fft%plan3D_bw, nz, ny, nx, CUFFT_Z2D, &
                           worksize)
    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Backward 3D FFT plan generation failed'
    end if

    ! allocate storage for cuFFTMp
    ierr = cufftXtMalloc(poisson_fft%plan3D_fw, poisson_fft%xtdesc, &
                         CUFFT_XT_FORMAT_INPLACE)
    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'cufftXtMalloc failed'
    end if

  end function init

  subroutine fft_forward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(in) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)

    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), d_dev, &
                     [self%nx_loc + 2, self%ny_loc, self%nz_loc])

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%ny_loc - 1)/tsize + 1, self%nz_loc, 1)
    threads = dim3(tsize, 1, 1)

    call memcpy3D<<<blocks, threads>>>(d_dev, padded_dev, & !&
                                       self%nx_loc, self%ny_loc, self%nz_loc)

    ierr = cufftXtExecDescriptor(self%plan3D_fw, self%xtdesc, self%xtdesc, &
                                 CUFFT_FORWARD)

    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Forward 3D FFT execution failed'
    end if

  end subroutine fft_forward_cuda

  subroutine fft_backward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)

    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    ierr = cufftXtExecDescriptor(self%plan3D_bw, self%xtdesc, self%xtdesc, &
                                 CUFFT_INVERSE)
    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Backward 3D FFT execution failed'
    end if

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), d_dev, &
                     [self%nx_loc + 2, self%ny_loc, self%nz_loc])

    tsize = 16
    blocks = dim3((self%ny_loc - 1)/tsize + 1, self%nz_loc, 1)
    threads = dim3(tsize, 1, 1)
    call memcpy3D<<<blocks, threads>>>(padded_dev, d_dev, & !&
                                       self%nx_loc, self%ny_loc, self%nz_loc)

  end subroutine fft_backward_cuda

  subroutine fft_postprocess_000_cuda(self)
    implicit none

    class(cuda_poisson_fft_t) :: self

    type(cudaXtDesc), pointer :: descriptor

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize

    ! obtain a pointer to descriptor so that we can carry out postprocessing
    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), c_dev, &
                     [self%nx_spec, self%ny_spec, self%nz_spec])

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%ny_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)

    ! Postprocess div_u in spectral space
    call process_spectral_000<<<blocks, threads>>>( & !&
      c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )

  end subroutine fft_postprocess_000_cuda

  subroutine fft_postprocess_010_cuda(self)
    implicit none

    class(cuda_poisson_fft_t) :: self

    type(cudaXtDesc), pointer :: descriptor

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize

    ! obtain a pointer to descriptor so that we can carry out postprocessing
    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), c_dev, &
                     [self%nx_spec, self%ny_spec, self%nz_spec])

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%nx_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)

    ! Postprocess div_u in spectral space
    call process_spectral_010<<<blocks, threads>>>( & !&
      c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )

  end subroutine fft_postprocess_010_cuda

  subroutine enforce_periodicity_y_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

  end subroutine enforce_periodicity_y_cuda

  subroutine undo_periodicity_y_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

  end subroutine undo_periodicity_y_cuda

end module m_cuda_poisson_fft
