module m_cuda_poisson_fft
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use cufftXt
  use cufft
  use mpi

  use m_allocator, only: field_t
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, CELL
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t

  use m_cuda_allocator, only: cuda_field_t
  use m_cuda_spectral, only: process_spectral_div_u

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
    procedure :: fft_postprocess => fft_postprocess_cuda
  end type cuda_poisson_fft_t

  interface cuda_poisson_fft_t
    module procedure init
  end interface cuda_poisson_fft_t

  private :: init

contains

  function init(mesh, xdirps, ydirps, zdirps) result(poisson_fft)
    implicit none

    class(mesh_t), intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    type(cuda_poisson_fft_t) :: poisson_fft

    integer :: nx, ny, nz

    integer :: ierr
    integer(int_ptr_kind()) :: worksize

    call poisson_fft%base_init(mesh, xdirps, ydirps, zdirps)

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

    real(dp), device, pointer :: flat_dev(:, :), d_dev(:, :, :)

    type(cudaXtDesc), pointer :: descriptor

    integer :: ierr

    select type (f)
    type is (cuda_field_t)
      flat_dev(1:self%nx_loc, 1:self%ny_loc*self%nz_loc) => f%data_d
    end select

    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), d_dev, &
                     [self%nx_loc + 2, self%ny_loc*self%nz_loc])
    ierr = cudaMemcpy2D(d_dev, self%nx_loc + 2, flat_dev, self%nx_loc, &
                        self%nx_loc, self%ny_loc*self%nz_loc)
    if (ierr /= 0) then
      print *, 'cudaMemcpy2D error code: ', ierr
      error stop 'cudaMemcpy2D failed'
    end if

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

    real(dp), device, pointer :: flat_dev(:, :), d_dev(:, :, :)

    type(cudaXtDesc), pointer :: descriptor

    integer :: ierr

    ierr = cufftXtExecDescriptor(self%plan3D_bw, self%xtdesc, self%xtdesc, &
                                 CUFFT_INVERSE)
    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Backward 3D FFT execution failed'
    end if

    select type (f)
    type is (cuda_field_t)
      flat_dev(1:self%nx_loc, 1:self%ny_loc*self%nz_loc) => f%data_d
    end select

    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), d_dev, &
                     [self%nx_loc + 2, self%ny_loc*self%nz_loc])
    ierr = cudaMemcpy2D(flat_dev, self%nx_loc, d_dev, self%nx_loc + 2, &
                        self%nx_loc, self%ny_loc*self%nz_loc)
    if (ierr /= 0) then
      print *, 'cudaMemcpy2D error code: ', ierr
      error stop 'cudaMemcpy2D failed'
    end if

  end subroutine fft_backward_cuda

  subroutine fft_postprocess_cuda(self)
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
    call process_spectral_div_u<<<blocks, threads>>>( & !&
      c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )

  end subroutine fft_postprocess_cuda

end module m_cuda_poisson_fft
