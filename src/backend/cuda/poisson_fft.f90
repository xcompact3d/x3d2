module m_cuda_poisson_fft
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer, c_int, c_float, &
                           c_double_complex, c_float_complex
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use cufftXt
  use cufft
  use mpi

  use m_common, only: dp, CELL, is_sp
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t

  use m_cuda_allocator, only: cuda_field_t
  use m_cuda_spectral, only: memcpy3D, &
                             process_spectral_000, process_spectral_010, &
                             enforce_periodicity_y, undo_periodicity_y, &
                             process_spectral_010_fw, &
                             process_spectral_010_poisson, &
                             process_spectral_010_bw

  implicit none

  type, extends(poisson_fft_t) :: cuda_poisson_fft_t
    !! FFT based Poisson solver

    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), device, allocatable, dimension(:, :, :) :: waves_dev
    !> Wave numbers in x, y, and z
    real(dp), device, allocatable, dimension(:) :: ax_dev, bx_dev, &
                                                   ay_dev, by_dev, &
                                                   az_dev, bz_dev
    !> Stretching operator matrices stores
    real(dp), device, allocatable, dimension(:, :, :, :) :: &
      store_a_odd_re_dev, store_a_odd_im_dev, &
      store_a_even_re_dev, store_a_even_im_dev, &
      store_a_re_dev, store_a_im_dev
    !> Stretching operator matrices
    real(dp), device, allocatable, dimension(:, :, :, :) :: &
      a_odd_re_dev, a_odd_im_dev, a_even_re_dev, a_even_im_dev, &
      a_re_dev, a_im_dev
    !> Forward and backward FFT transform plans
    integer :: plan3D_fw, plan3D_bw

    !> Flag to indicate whether cuFFTMp (multi-GPU) is used
    logical :: use_cufftmp = .true.

    !> cuFFTMp object manages decomposition and data storage (multi-GPU)
    type(cudaLibXtDesc), pointer :: xtdesc

    !> Standard cuFFT storage for single-GPU mode
    complex(dp), device, allocatable, dimension(:, :, :) :: c_dev
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

  ! Explicit C interfaces for cuFFT functions that nvfortran has trouble with
  interface
    integer(c_int) function cufftExecR2C_C(plan, idata, odata) &
      bind(C, name='cufftExecR2C')
      use iso_c_binding
      integer(c_int), value :: plan
      type(c_ptr), value :: idata
      type(c_ptr), value :: odata
    end function cufftExecR2C_C

    integer(c_int) function cufftExecC2R_C(plan, idata, odata) &
      bind(C, name='cufftExecC2R')
      use iso_c_binding
      integer(c_int), value :: plan
      type(c_ptr), value :: idata
      type(c_ptr), value :: odata
    end function cufftExecC2R_C
  end interface

  private :: init, create_fft_plan

contains

  subroutine create_fft_plan(plan, use_cufftmp, nx, ny, nz, &
                             plan_type, is_root, plan_name)
    !! Helper subroutine to create FFT plan with automatic cuFFTMp fallback
    implicit none

    integer, intent(inout) :: plan
    logical, intent(inout) :: use_cufftmp
    integer, intent(in) :: nx, ny, nz, plan_type
    logical, intent(in) :: is_root
    character(*), intent(in) :: plan_name

    integer :: ierr
    integer(int_ptr_kind()) :: worksize
    logical :: cufftmp_failed

    cufftmp_failed = .false.
    ierr = cufftCreate(plan)

    if (use_cufftmp) then
      ! Try to attach MPI communicator for cuFFTMp
      ierr = cufftMpAttachComm(plan, CUFFT_COMM_MPI, MPI_COMM_WORLD)
      if (ierr /= 0) then
        if (is_root) then
          print *, 'cuFFTMp: MPI attach failed (error ', ierr, ')'
        end if
        cufftmp_failed = .true.
      else
        ! MPI attach succeeded, create the plan
        ierr = cufftMakePlan3D(plan, nz, ny, nx, plan_type, worksize)
        if (ierr /= 0) then
          if (is_root) then
            print *, 'cuFFTMp: Plan creation failed (error ', ierr, ')'
          end if
          cufftmp_failed = .true.
        end if
      end if

      ! if cuFFTMp failed at any stage, fall back to cuFFT
      if (cufftmp_failed) then
        if (is_root) then
          print *, 'Falling back to single-GPU cuFFT'
        end if
        use_cufftmp = .false.
        ierr = cufftDestroy(plan)
        ierr = cufftCreate(plan)
      end if
    end if

    ! create plan with cuFFT
    if (.not. use_cufftmp) then
      ierr = cufftMakePlan3D(plan, nz, ny, nx, plan_type, worksize)
      if (ierr /= 0) then
        write (stderr, *), 'cuFFT Error Code: ', ierr
        error stop trim(plan_name)//' 3D FFT plan generation failed'
      end if
    end if

  end subroutine create_fft_plan

  function init(mesh, xdirps, ydirps, zdirps, lowmem) &
    result(poisson_fft)
    implicit none

    type(mesh_t), intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, optional, intent(in) :: lowmem

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

    ! will store the a matrix coefficients in GPU memory if (.not. lowmem)
    ! and do a device-to-device copy at each iter. Otherwise copy from host.
    ! lowmem is .false. by default
    if (present(lowmem)) poisson_fft%lowmem = lowmem

    ! Try cuFFTMp (multi-GPU) first, with automatic fallback to cuFFT if not supported
    poisson_fft%use_cufftmp = .true.

    ! if stretching in y is 'centred' or 'top-bottom'
    if (poisson_fft%stretched_y .and. poisson_fft%stretched_y_sym) then
      poisson_fft%a_odd_re_dev = poisson_fft%a_odd_re
      poisson_fft%a_odd_im_dev = poisson_fft%a_odd_im
      poisson_fft%a_even_re_dev = poisson_fft%a_even_re
      poisson_fft%a_even_im_dev = poisson_fft%a_even_im
      if (.not. poisson_fft%lowmem) then
        poisson_fft%store_a_odd_re_dev = poisson_fft%a_odd_re
        poisson_fft%store_a_odd_im_dev = poisson_fft%a_odd_im
        poisson_fft%store_a_even_re_dev = poisson_fft%a_even_re
        poisson_fft%store_a_even_im_dev = poisson_fft%a_even_im
      end if
    !! if stretching in y is 'bottom'
    else if (poisson_fft%stretched_y .and. &
             (.not. poisson_fft%stretched_y_sym)) then
      poisson_fft%a_re_dev = poisson_fft%a_re
      poisson_fft%a_im_dev = poisson_fft%a_im
      if (.not. poisson_fft%lowmem) then
        poisson_fft%store_a_re_dev = poisson_fft%a_re
        poisson_fft%store_a_im_dev = poisson_fft%a_im
      end if
    end if

    ! Create forward FFT plan with automatic cuFFTMp detection/fallback
    if (is_sp) then
      call create_fft_plan(poisson_fft%plan3D_fw, &
                           poisson_fft%use_cufftmp, &
                           nx, ny, nz, CUFFT_R2C, &
                           mesh%par%is_root(), 'Forward')
    else
      call create_fft_plan(poisson_fft%plan3D_fw, &
                           poisson_fft%use_cufftmp, &
                           nx, ny, nz, CUFFT_D2Z, &
                           mesh%par%is_root(), 'Forward')
    end if

    ! Create backward FFT plan with automatic cuFFTMp detection/fallback
    if (is_sp) then
      call create_fft_plan(poisson_fft%plan3D_bw, &
                           poisson_fft%use_cufftmp, &
                           nx, ny, nz, CUFFT_C2R, &
                           mesh%par%is_root(), 'Backward')
    else
      call create_fft_plan(poisson_fft%plan3D_bw, &
                           poisson_fft%use_cufftmp, &
                           nx, ny, nz, CUFFT_Z2D, &
                           mesh%par%is_root(), 'Backward')
    end if

    ! Allocate storage - cuFFTMp uses xtdesc, single-GPU uses c_dev
    if (poisson_fft%use_cufftmp) then
      ! allocate storage for cuFFTMp
      ierr = cufftXtMalloc(poisson_fft%plan3D_fw, poisson_fft%xtdesc, &
                           CUFFT_XT_FORMAT_INPLACE)
      if (ierr /= 0) then
        write (stderr, *), 'cuFFT Error Code: ', ierr
        error stop 'cufftXtMalloc failed'
      end if
    else
      ! allocate storage for single-GPU cuFFT
      allocate (poisson_fft%c_dev(poisson_fft%nx_spec, &
                                  poisson_fft%ny_spec, &
                                  poisson_fft%nz_spec))
    end if

    ! Print final status
    if (mesh%par%is_root()) then
      if (poisson_fft%use_cufftmp) then
        print *, 'Using cuFFTMp for multi-GPU FFT'
      else
        print *, 'Using cuFFT for single-GPU FFT'
      end if
    end if

  end function init

  subroutine fft_forward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(in) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)
    real(dp), device, pointer :: f_ptr
    type(c_ptr) :: f_c_ptr

    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    if (self%use_cufftmp) then
      ! Multi-GPU path using cuFFTMp
      ! tsize is different than SZ, because here we work on a 3D Cartesian
      ! data structure, and free to specify any suitable thread/block size.
      tsize = 16
      blocks = dim3((self%ny_loc - 1)/tsize + 1, self%nz_loc, 1)
      threads = dim3(tsize, 1, 1)

      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), d_dev, &
                       [self%nx_loc + 2, self%ny_loc, self%nz_loc])

      call memcpy3D<<<blocks, threads>>>(d_dev, padded_dev, & !&
                                         self%nx_loc, self%ny_loc, self%nz_loc)

      ierr = cufftXtExecDescriptor(self%plan3D_fw, self%xtdesc, self%xtdesc, &
                                   CUFFT_FORWARD)
    else
      ! Single-GPU path using standard cuFFT
      ! Using padded_dev directly causes segfault, use pointer workaround
      f_c_ptr = c_loc(padded_dev)
      call c_f_pointer(f_c_ptr, f_ptr)

      ! Use explicit C interface for single precision, Fortran interface for double
#ifdef SINGLE_PREC
      ierr = cufftExecR2C_C(self%plan3D_fw, c_loc(f_ptr), c_loc(self%c_dev))
#else
      ierr = cufftExecD2Z(self%plan3D_fw, f_ptr, self%c_dev)
#endif
    end if

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
    real(dp), device, pointer :: f_ptr
    type(c_ptr) :: f_c_ptr

    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    if (self%use_cufftmp) then
      ! Multi-GPU path using cuFFTMp
      ierr = cufftXtExecDescriptor(self%plan3D_bw, self%xtdesc, self%xtdesc, &
                                   CUFFT_INVERSE)
    else
      ! Single-GPU path using standard cuFFT
      ! Using padded_dev directly causes segfault, use pointer workaround
      f_c_ptr = c_loc(padded_dev)
      call c_f_pointer(f_c_ptr, f_ptr)

      ! Use explicit C interface for single precision, Fortran interface for double
#ifdef SINGLE_PREC
      ierr = cufftExecC2R_C(self%plan3D_bw, c_loc(self%c_dev), c_loc(f_ptr))
#else
      ierr = cufftexecZ2D(self%plan3D_bw, self%c_dev, f_ptr)
#endif
    end if

    if (ierr /= 0) then
      write (stderr, *), 'cuFFT Error Code: ', ierr
      error stop 'Backward 3D FFT execution failed'
    end if

    if (self%use_cufftmp) then
      ! Multi-GPU path: copy from cuFFTMp storage
      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), d_dev, &
                       [self%nx_loc + 2, self%ny_loc, self%nz_loc])

      tsize = 16
      blocks = dim3((self%ny_loc - 1)/tsize + 1, self%nz_loc, 1)
      threads = dim3(tsize, 1, 1)
      call memcpy3D<<<blocks, threads>>>(padded_dev, d_dev, & !&
                                         self%nx_loc, self%ny_loc, self%nz_loc)
    end if
    ! Single-GPU path: data already in padded_dev from cufftExec

  end subroutine fft_backward_cuda

  subroutine fft_postprocess_000_cuda(self)
    implicit none

    class(cuda_poisson_fft_t) :: self

    type(cudaXtDesc), pointer :: descriptor

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%ny_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)

    ! Get pointer to the appropriate FFT data storage
    if (self%use_cufftmp) then
      ! Multi-GPU path: get pointer from xtdesc
      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    else
      ! Single-GPU path: get pointer to self%c_dev
      call c_f_pointer(c_loc(self%c_dev), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    end if

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
    integer :: tsize, off, inc

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%nx_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)

    ! Get pointer to the appropriate FFT data storage
    if (self%use_cufftmp) then
      ! Multi-GPU path: get pointer from xtdesc
      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    else
      ! Single-GPU path: get pointer to self%c_dev
      call c_f_pointer(c_loc(self%c_dev), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    end if

    ! Postprocess div_u in spectral space
    if (.not. self%stretched_y) then
      call process_spectral_010<<<blocks, threads>>>( & !&
        c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
        self%nx_glob, self%ny_glob, self%nz_glob, &
        self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
        self%az_dev, self%bz_dev &
        )
    else
      call process_spectral_010_fw<<<blocks, threads>>>( & !&
        c_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
        self%nx_glob, self%ny_glob, self%nz_glob, &
        self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
        self%az_dev, self%bz_dev &
        )

      ! if stretching in y is 'centred' or 'top-bottom'
      if (self%stretched_y_sym) then
        ! copy from host to device if lowmem else from device stores
        if (self%lowmem) then
          self%a_odd_re_dev = self%a_odd_re
          self%a_odd_im_dev = self%a_odd_im
          self%a_even_re_dev = self%a_even_re
          self%a_even_im_dev = self%a_even_im
        else
          self%a_odd_re_dev = self%store_a_odd_re_dev
          self%a_odd_im_dev = self%store_a_odd_im_dev
          self%a_even_re_dev = self%store_a_even_re_dev
          self%a_even_im_dev = self%store_a_even_im_dev
        end if
        ! start from the first odd entry and continue with odd ones
        off = 0
        inc = 2
        call process_spectral_010_poisson<<<blocks, threads>>>( & !&
          c_dev, self%a_odd_re_dev, self%a_odd_im_dev, off, inc, &
          self%nx_spec, self%ny_spec/2, &
          self%nx_glob, self%ny_glob, self%nz_glob &
          )
        ! start from the first even entry and continue with even ones
        off = 1
        call process_spectral_010_poisson<<<blocks, threads>>>( & !&
          c_dev, self%a_even_re_dev, self%a_even_im_dev, off, inc, &
          self%nx_spec, self%ny_spec/2, &
          self%nx_glob, self%ny_glob, self%nz_glob &
          )
      !! if stretching in y is 'bottom'
      else
        ! copy from host to device if lowmem else from device stores
        if (self%lowmem) then
          self%a_re_dev = self%a_re
          self%a_im_dev = self%a_im
        else
          self%a_re_dev = self%store_a_re_dev
          self%a_im_dev = self%store_a_im_dev
        end if
        ! start from the first entry and increment 1
        off = 0
        inc = 1
        call process_spectral_010_poisson<<<blocks, threads>>>( & !&
          c_dev, self%a_re_dev, self%a_im_dev, off, inc, &
          self%nx_spec, self%ny_spec, &
          self%nx_glob, self%ny_glob, self%nz_glob &
          )
      end if

      call process_spectral_010_bw<<<blocks, threads>>>( & !&
        c_dev, self%nx_spec, self%ny_spec, self%y_sp_st, &
        self%nx_glob, self%ny_glob, self%nz_glob, &
        self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
        self%az_dev, self%bz_dev &
        )
    end if

  end subroutine fft_postprocess_010_cuda

  subroutine enforce_periodicity_y_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    blocks = dim3(self%nz_loc, 1, 1)
    threads = dim3(self%nx_loc, 1, 1)
    call enforce_periodicity_y<<<blocks, threads>>>( & !&
      f_out_dev, f_in_dev, self%ny_glob &
      )

  end subroutine enforce_periodicity_y_cuda

  subroutine undo_periodicity_y_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    blocks = dim3(self%nz_loc, 1, 1)
    threads = dim3(self%nx_loc, 1, 1)
    call undo_periodicity_y<<<blocks, threads>>>( & !&
      f_out_dev, f_in_dev, self%ny_glob &
      )

  end subroutine undo_periodicity_y_cuda

end module m_cuda_poisson_fft
