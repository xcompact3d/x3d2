program test_reorder
    use iso_fortran_env, only: stderr => error_unit
    use mpi

    use m_allocator, only: allocator_t, field_t
    use m_base_backend, only: base_backend_t
    use m_solver, only: allocate_tdsops
    use m_tdsops, only: dirps_t, tdsops_t

    use m_common, only: dp, pi, globs_t, set_pprev_pnext, &
                       RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y, dir_X, dir_Y, dir_Z

    use m_ordering, only: get_index_dir, get_index_ijk

#ifdef CUDA
   use m_cuda_common, only: SZ
   use m_cuda_tdsops, only: cuda_tdsops_t
#else
   use m_omp_common, only: SZ
   use m_omp_backend, only: omp_backend_t
#endif

    implicit none

    logical :: allpass = .true.
    class(field_t), pointer :: u_x, u_y, u_z, u_x_original

    integer :: nrank, nproc
    integer :: ierr, i, j, k, i_, j_, k_, i__, j__, k__

    real(dp) :: dx, dx_per

    type(globs_t) :: globs
    class(base_backend_t), pointer :: backend
    class(allocator_t), pointer :: allocator
    type(dirps_t), target :: xdirps, ydirps, zdirps

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
 
    globs%nx = 96
    globs%ny = 96
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
    cuda_allocator = cuda_allocator_t([SZ, globs%nx_loc, globs%n_groups_x])
    allocator => cuda_allocator
    print*, 'CUDA allocator instantiated'

    cuda_backend = cuda_backend_t(globs, allocator)
    backend => cuda_backend
    print*, 'CUDA backend instantiated'
#else
    omp_allocator = allocator_t([SZ, globs%nx_loc, globs%n_groups_x])
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

   ! Test indexing only
    do k=1, xdirps%n
        do j=1, xdirps%n
            do i=1, xdirps%n
                call get_index_dir(i_, j_, k_, i, j, k, dir_X, SZ, xdirps%n, xdirps%n, xdirps%n)
                call get_index_ijk(i__, j__, k__, i_, j_, k_, dir_X, SZ, xdirps%n, xdirps%n, xdirps%n)
                if (i /= i__ .or. j /= j__ .or. k/= k__) then
                    print *, "Error in dir X"
                    allpass = .false.
                end if
                call get_index_dir(i_, j_, k_, i, j, k, dir_Y, SZ, xdirps%n, xdirps%n, xdirps%n)
                call get_index_ijk(i__, j__, k__, i_, j_, k_, dir_Y, SZ, xdirps%n, xdirps%n, xdirps%n)
                if (i /= i__ .or. j /= j__ .or. k/= k__) then
                    print *, "Error in dir Y"
                    allpass = .false.
                end if
                call get_index_dir(i_, j_, k_, i, j, k, dir_Z, SZ, xdirps%n, xdirps%n, xdirps%n)
                call get_index_ijk(i__, j__, k__, i_, j_, k_, dir_Z, SZ, xdirps%n, xdirps%n, xdirps%n)
                if (i /= i__ .or. j /= j__ .or. k/= k__) then
                    print *, "Error in dir Z"
                    allpass = .false.
                end if
            end do
        end do
    end do
 

    ! Test reordering
    u_x => allocator%get_block()
    u_y => allocator%get_block()
    u_z => allocator%get_block()
    u_x_original => allocator%get_block()

    call random_number(u_x_original%data)

    u_x%data(:, :, :) = u_x_original%data(:, :, :)

    call backend%reorder(u_y, u_x, RDR_X2Y)
    call backend%reorder(u_x, u_y, RDR_Y2X)
    call check_reorder(allpass, u_x, u_x_original, "X2Y, Y2X")

    call backend%reorder(u_z, u_x, RDR_X2Z)
    call backend%reorder(u_x, u_z, RDR_Z2X)
    call check_reorder(allpass, u_x, u_x_original, "X2Z, Z2X")

    call backend%reorder(u_z, u_y, RDR_Y2Z)
    call backend%reorder(u_x, u_z, RDR_Z2X)
    call check_reorder(allpass, u_x, u_x_original, "Y2Z, Z2X")

    call backend%reorder(u_y, u_z, RDR_Z2Y)
    call backend%reorder(u_x, u_y, RDR_Y2X)
    call check_reorder(allpass, u_x, u_x_original, "Z2Y, Y2X")

    if (allpass) then
        if (nrank == 0) write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
    else
        error stop 'SOME TESTS FAILED.'
    end if

    call MPI_Finalize(ierr)

    contains

    subroutine check_reorder(allpass, a, b, message)
        logical, intent(inout) :: allpass
        class(field_t), intent(in) :: a, b
        character(len=*), intent(in) :: message
        real(dp) :: tol = 1d-8

        if (norm2(a%data - b%data) > tol) then
            allpass = .false.
            write(stderr, '(a)') message
        end if

    end subroutine

end program test_reorder

