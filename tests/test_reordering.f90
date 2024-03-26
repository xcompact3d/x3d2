  program test_reorder
    use iso_fortran_env, only: stderr => error_unit
    use mpi

    use m_allocator, only: allocator_t, field_t
    use m_base_backend, only: base_backend_t
    use m_tdsops, only: dirps_t

    use m_common, only: dp, pi, globs_t, set_pprev_pnext, &
      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y, &
      DIR_X, DIR_Y, DIR_Z

    use m_ordering, only: get_index_dir, get_index_ijk

#ifdef CUDA
    use cudafor

    use m_cuda_allocator, only: cuda_allocator_t, cuda_field_t
    use m_cuda_backend, only: cuda_backend_t
    use m_cuda_common, only: SZ
#else
    use m_omp_common, only: SZ
    use m_omp_backend, only: omp_backend_t
#endif

    implicit none

    logical :: allpass = .true.
    class(field_t), pointer :: u_x, u_y, u_z
    class(field_t), pointer :: u_x_original

    real(dp), allocatable, dimension(:, :, :) :: u_array, temp_1, temp_2

    integer :: dims(3)

    integer :: nrank, nproc
    integer :: ierr, i, j, k

    real(dp) :: dx, dx_per

    type(globs_t) :: globs
    class(base_backend_t), pointer :: backend
    class(allocator_t), pointer :: allocator
    type(dirps_t), target :: xdirps, ydirps, zdirps
    logical :: pass_X, pass_Y, pass_Z

#ifdef CUDA
    type(cuda_backend_t), target :: cuda_backend
    type(cuda_allocator_t), target :: cuda_allocator
    integer :: ndevs, devnum
#else
    type(omp_backend_t), target :: omp_backend
    type(allocator_t), target :: omp_allocator
#endif

    ! Initialise variables and arrays
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

#ifdef CUDA
    ierr = cudaGetDeviceCount(ndevs)
    ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
    ierr = cudaGetDevice(devnum)
#endif

    globs%nx = 32
    globs%ny = 64
    globs%nz = 96

    globs%nx_loc = globs%nx/nproc
    globs%ny_loc = globs%ny/nproc
    globs%nz_loc = globs%nz/nproc

    globs%n_groups_x = globs%ny_loc*globs%nz_loc/SZ
    globs%n_groups_y = globs%nx_loc*globs%nz_loc/SZ
    globs%n_groups_z = globs%nx_loc*globs%ny_loc/SZ

    xdirps%n = globs%nx_loc
    ydirps%n = globs%ny_loc
    zdirps%n = globs%nz_loc

    xdirps%n_blocks = globs%n_groups_x
    ydirps%n_blocks = globs%n_groups_y
    zdirps%n_blocks = globs%n_groups_z


#ifdef CUDA
    cuda_allocator = cuda_allocator_t(globs%nx_loc, globs%ny_loc, globs%nz_loc, SZ)
    allocator => cuda_allocator
    print*, 'CUDA allocator instantiated'

    cuda_backend = cuda_backend_t(globs, allocator)
    backend => cuda_backend
    print*, 'CUDA backend instantiated'
#else
    omp_allocator = allocator_t(globs%nx_loc, globs%ny_loc, globs%nz_loc, SZ)
    allocator => omp_allocator
    print*, 'OpenMP allocator instantiated'

    omp_backend = omp_backend_t(globs, allocator)
    backend => omp_backend
    print*, 'OpenMP backend instantiated'
#endif

    backend%xdirps => xdirps
    backend%ydirps => ydirps
    backend%zdirps => zdirps

    if (nrank == 0) print*, 'Parallel run with', nproc, 'ranks'
    pass_X = .true.
    pass_Y = .true.
    pass_Z = .true.

    ! Test indexing only
    do k=1, zdirps%n
      do j=1, ydirps%n
        do i=1, xdirps%n
          call test_index_reversing(pass_X, i, j, k, DIR_X, SZ, xdirps%n, ydirps%n, zdirps%n)
          call test_index_reversing(pass_Y, i, j, k, DIR_Y, SZ, xdirps%n, ydirps%n, zdirps%n)
          call test_index_reversing(pass_Z, i, j, k, DIR_Z, SZ, xdirps%n, ydirps%n, zdirps%n)
        end do
      end do
    end do
    if (.not. pass_X) print *, "Error in X direction for index reversing"
    if (.not. pass_Y) print *, "Error in Y direction for index reversing"
    if (.not. pass_Z) print *, "Error in Z direction for index reversing"

    allpass = (pass_X .and. pass_Y .and. pass_Z)


    ! Test reordering
    u_x => allocator%get_block(DIR_X)
    u_y => allocator%get_block(DIR_Y)
    u_z => allocator%get_block(DIR_Z)
    u_x_original => allocator%get_block(DIR_X)

    dims(:) = allocator%xdims_padded(:)
    allocate (u_array(dims(1), dims(2), dims(3)))

    call random_number(u_array)

#ifdef CUDA
    allocate (temp_1(dims(1), dims(2), dims(3)))
    allocate (temp_2(dims(1), dims(2), dims(3)))

    select type (u_x_original)
     type is (cuda_field_t)
      u_x_original%data_d = u_array
    end select
#else
    select type (u_x_original)
     type is (field_t)
      u_x_original%data = u_array
    end select
#endif

    call backend%reorder(u_y, u_x_original, RDR_X2Y)
    call backend%reorder(u_x, u_y, RDR_Y2X)
    call check_reorder(allpass, u_x, u_x_original, "testing X2Y and Y2X failed")

    call backend%reorder(u_z, u_x, RDR_X2Z)
    call backend%reorder(u_x, u_z, RDR_Z2X)
    call check_reorder(allpass, u_x, u_x_original, "testing X2Z and Z2X failed")

    call backend%reorder(u_z, u_y, RDR_Y2Z)
    call backend%reorder(u_x, u_z, RDR_Z2X)
    call check_reorder(allpass, u_x, u_x_original, "testing Y2Z and Z2X failed")

    call backend%reorder(u_y, u_z, RDR_Z2Y)
    call backend%reorder(u_x, u_y, RDR_Y2X)
    call check_reorder(allpass, u_x, u_x_original, "testing Z2Y and Y2X failed")

    if (allpass) then
      if (nrank == 0) write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
    else
      error stop 'SOME TESTS FAILED.'
    end if

    call MPI_Finalize(ierr)

  contains

    subroutine test_index_reversing(pass, i, j, k, dir, SZ, nx, ny, nz)
      logical, intent(inout) :: pass
      integer, intent(in) :: i, j, k    ! original indices in the cartesian space
      integer, intent(in) :: dir
      integer, intent(in) :: SZ, nx, ny, nz
      integer :: dir_i, dir_j, dir_k    ! indices in the applicatin storage direction
      integer :: cart_i, cart_j, cart_k ! newly computed indices in the cartesian space

      call get_index_dir(dir_i, dir_j, dir_k, i, j, k, dir, SZ, nx, ny, nz)
      call get_index_ijk(cart_i, cart_j, cart_k, dir_i, dir_j, dir_k, dir, SZ, nx, ny, nz)

      if (i /= cart_i .or. j /= cart_j .or. k /= cart_k) then
        pass = .false.
      end if

    end subroutine

    subroutine check_reorder(allpass, a, b, message)
      logical, intent(inout) :: allpass
      class(field_t), intent(in) :: a, b
      character(len=*), intent(in) :: message
      real(dp) :: tol = 1d-8

#ifdef CUDA
      select type (a); type is (cuda_field_t); temp_1 = a%data_d; end select
      select type (b); type is (cuda_field_t); temp_2 = b%data_d; end select
      if (norm2(temp_1 - temp_2) > tol) then
        allpass = .false.
        write(stderr, '(a)') message
      end if
#else
      if (norm2(a%data - b%data) > tol) then
        allpass = .false.
        write(stderr, '(a)') message
      end if
#endif

    end subroutine


  end program test_reorder

