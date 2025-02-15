program test_adios2
    use mpi
    use m_base_adios2
    use m_adios2_writer
    use m_adios2_reader
    use m_common, only: dp
    use iso_fortran_env, only: int64, real64, stderr => error_unit
    implicit none

    ! ADIOS2 handlers
    type(adios2_writer_t) :: writer
    type(adios2_reader_t) :: reader

    ! MPI variables
    integer :: ierr, irank, isize
    integer(kind=int64), dimension(2) :: shape_dims, start_dims, count_dims
    integer(kind=int64), dimension(2) :: sel_start, sel_count
    real(kind=real64), dimension(:,:), allocatable :: data_write, data_read

    ! application variables
    integer :: i, j, rank_id, inx = 3, iny = 4
    logical :: allpass = .true.
    real(kind=real64) :: expected

    ! launch MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

    ! data initialization
    allocate(data_write(inx, iny))

    ! initialize data (unique per rank)
    do j = 1, iny
        do i = 1, inx
            data_write(i,j) = irank*inx*iny + (j-1)*inx + (i-1)
        end do
    end do

    ! global shape and local offsets for parallel I/O
    shape_dims = [isize*inx, iny]
    start_dims = [irank*inx, 0]
    count_dims = [inx, iny]

    ! write data
    call writer%init_adios2(MPI_COMM_WORLD, "test_io_write")
    call writer%open_adios2("test_output.bp", adios2_mode_write)
    call writer%begin_step_adios2()
    call writer%write_data("data2D", data_write, shape_dims, start_dims, count_dims)
    call writer%end_step_adios2()
    call writer%close_adios2()

    if (allocated(data_write)) deallocate(data_write)

    ! read data (rank 0 only)
    if (irank == 0) then
        call reader%init_adios2(MPI_COMM_SELF, "test_io_read")
        call reader%open_adios2("test_output.bp", adios2_mode_read)
        call reader%begin_step_adios2()

        sel_start = [0, 0]
        sel_count = [shape_dims(1), shape_dims(2)]

        ! read entire dataset
        allocate(data_read(sel_count(1), sel_count(2)))
        call reader%read_data("data2D", data_read, sel_start, sel_count)

        ! verify data
        do j = 1, sel_count(2)
            do i = 1, sel_count(1)
                ! calculate expected value based on original data pattern
                rank_id = (i-1) / inx ! determine which rank wrote this data
                expected = rank_id*inx*iny + (j-1)*inx + mod(i-1,inx)
                if (abs(data_read(i,j) - expected ) > 1.e-8_dp) then
                    allpass = .false.
                    print *, "Data mismatch at (", i, ",", j, "): ", data_read(i,j)
                end if
            end do
        end do

        call reader%end_step_adios2()
        call reader%close_adios2()

        if (allocated(data_read)) deallocate(data_read)
    end if

    ! Cleanup and finalize
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Finalize(ierr)

    ! Test result
    if (irank == 0) then
        if (allpass) then
            write(stderr,'(a)') 'ADIOS2 TEST PASSED SUCCESSFULLY.'
        else
            error stop 'ADIOS2 TEST FAILED.'
        end if
    end if

end program test_adios2