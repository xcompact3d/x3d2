program test_adios2
   use mpi
   use adios2
   use m_common, only: dp
   use iso_fortran_env, only: int64, stderr => error_unit

   implicit none

   ! adios2 handlers
   type(adios2_adios) :: adios
   type(adios2_io)    :: io_write, io_read
   type(adios2_variable) :: var
   type(adios2_engine) :: engine

   ! MPI variables
   integer :: ierr, irank, isize
   integer(kind=int64), dimension(2) :: shape_dims, start_dims, count_dims
   integer(kind=int64), dimension(2) :: sel_start, sel_count
   real(dp), dimension(:, :), allocatable :: data_write, data_read

   ! application variables
   integer :: i, j, inx = 3, iny = 4
   logical :: allpass = .true.

   ! launch MPI
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, isize, ierr)

   ! data initialisation
   allocate (data_write(inx, iny))

   ! initialize data (unique per rank)
   do j = 1, iny
      do i = 1, inx
         data_write(i, j) = irank*inx*iny + (j - 1)*inx + (i - 1)
      end do
   end do

   ! global shape and local offsets for parallel I/O
   shape_dims = [isize*inx, iny]
   start_dims = [irank*inx, 0]
   count_dims = [inx, iny]

   ! create adios handler passing the communicator and error flag
   call adios2_init(adios, MPI_COMM_WORLD, ierr)

   ! declare IO process configuration inside adios
   call adios2_declare_io(io_write, adios, "test_io_write", ierr)

   ! define adios2 variable to be written in bp format
   call adios2_define_variable(var, io_write, "data2D", adios2_type_real, 2, &
                               shape_dims, start_dims, count_dims, &
                               adios2_constant_dims, ierr)

   ! open in write mode, this launches an engine
   call adios2_open(engine, io_write, "test_output.bp", adios2_mode_write, ierr)

   call adios2_begin_step(engine, ierr)

   ! put data content to bp buffer
   call adios2_put(engine, var, data_write, adios2_mode_sync, ierr)
   call adios2_end_step(engine, ierr)

   ! closes engine1 and deallocates it, becomes unreachable
   call adios2_close(engine, ierr)
   if (allocated(data_write)) deallocate (data_write)

   ! reading back the data (only rank 0)
   if (irank == 0) then
      call adios2_init(adios, MPI_COMM_SELF, ierr)
      call adios2_declare_io(io_read, adios, "test_io_read", ierr)

      ! open in read mode, this launches an engine
      call adios2_open(engine, io_read, "test_output.bp", adios2_mode_read, ierr)

      call adios2_begin_step(engine, ierr)
      call adios2_inquire_variable(var, io_read, "data2D", ierr)

      if (ierr == adios2_found) then

         sel_start = [0, 0]
         sel_count = [shape_dims(1), shape_dims(2)]

         ! select entire dataset
         allocate (data_read(sel_count(1), sel_count(2)))
         call adios2_set_selection(var, 2, sel_start, sel_count, ierr)
         call adios2_get(engine, var, data_read, adios2_mode_sync, ierr)

         ! data verification
         do j = 1, sel_count(2)
            do i = 1, sel_count(1)
               if (abs(data_read(i, j) - real((j - 1)*shape_dims(1) + (i - 1), kind=8)) > 1.e-8_dp) then
                  allpass = .false.
                  print *, "Data mismatch at (", i, ",", j, "): ", data_read(i, j)
               end if
            end do
         end do

         if (allocated(data_read)) deallocate (data_read)
      end if

      call adios2_end_step(engine, ierr)
      call adios2_close(engine, ierr)
   end if

   ! cleanup
   call adios2_finalize(adios, ierr)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   call MPI_Finalize(ierr)

   ! test result
   if (irank == 0) then
      if (allpass) then
         write (stderr, '(a)') 'ADIOS2 TEST PASSED SUCCESSFULLY.'
      else
         error stop 'ADIOS2 TEST FAILED.'
      end if
   end if

end program test_adios2
