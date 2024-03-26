program test_omp_io
   use iso_fortran_env, only: stderr => error_unit
   use adios2
   use mpi

   use m_allocator, only: allocator_t, field_t
   use m_common, only: dp, DIR_C
   use m_omp_common, only: SZ

   implicit none

   logical :: allpass !! flag indicating all tests pass

   class(field_t), pointer :: arr_to_write !! array to save & restore

   !> MPI vars
   integer :: irank, nproc
   !> Error code, used for MPI and ADIOS
   integer :: ierr

   type(adios2_adios) :: adios_ctx

   type(allocator_t), target :: omp_allocator

   integer(kind=8), dimension(3) :: ishape, istart, icount
   integer :: n, n_block, n_glob
   integer, parameter :: nx = 64, ny = 32, nz = 16
   integer, parameter :: n_groups_x = 1
   character(*), parameter :: varname = "TestArr"
   character(*), parameter :: fname = "__FILE__.bp"

   integer :: i, j, k
   character(len=80)::fmt

   allpass = .true.

   call init_mpi_adios()

   ! if (nrank == 0 && nproc > 1) print*, 'Parallel run with', nproc, 'ranks'
   if (nproc > 1) call abort_test("Test does not support multiple MPI processes")

   !================ 
   ! SETUP TEST DATA
   !================

   icount = (/ nx, ny, nz/) ! global size
   istart = (/ 0, 0, 0/) ! local offset
   ishape = (/ nx, ny, nz/) ! local size

   omp_allocator = allocator_t(nx, ny, nz, SZ)
   print*, 'OpenMP allocator instantiated'

   arr_to_write => omp_allocator%get_block(DIR_C)

   ! Initialise with a simple index to verify read later
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            arr_to_write%data(i, j, k) = i + j*nx + k*nx*ny
         end do
      end do
   end do

   call write_test_file()
   call read_test_file()

   if (allpass) then
      if (irank == 0) write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
   else
      error stop 'SOME TESTS FAILED.'
   end if

   call deinit_mpi_adios()

   contains

     subroutine init_mpi_adios()
       call MPI_Init(ierr)
       call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
       call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

       call adios2_init(adios_ctx, MPI_COMM_WORLD, ierr)
       if (.not.adios_ctx%valid) then
          call abort_test("Cannot create ADIOS2 context")
       endif
     end subroutine

     subroutine deinit_mpi_adios()
       if (adios_ctx%valid) then
         call adios2_finalize(adios_ctx, ierr)
       endif

       call MPI_Finalize(ierr)
     end subroutine

     subroutine abort_test(msg)
       character(*) :: msg

       call deinit_mpi_adios()

       error stop msg
     end subroutine

     subroutine write_test_file()
       type(adios2_io) :: io
       type(adios2_engine) :: engine
       type(adios2_variable) :: adios_var

       integer :: ierr

       call adios2_declare_io (io, adios_ctx, 'TestWrite', ierr)
       if (.not.io%valid) then
          call abort_test("Cannot create ADIOS2 IO")
       endif

       call adios2_open(engine, io, fname, adios2_mode_write, ierr)
       if (.not.engine%valid) then
          call abort_test("Cannot create ADIOS2 engine")
       endif

       ! TODO tie double precision variable to a
       call adios2_define_variable(adios_var, io, varname, adios2_type_dp, &
                                   3, ishape, istart, icount, adios2_constant_dims, ierr)

       call adios2_begin_step(engine, adios2_step_mode_append, ierr)
       call adios2_put(engine, adios_var, arr_to_write%data, ierr)
       call adios2_end_step(engine, ierr)

       if (engine%valid) then
          call adios2_close(engine, ierr)
       endif

     end subroutine

     subroutine read_test_file()
       type(adios2_io) :: io
       type(adios2_engine) :: engine
       type(adios2_variable) :: adios_var
       class(field_t), pointer :: arr_to_read !! array to save & restore
       integer :: ierr

       if (irank /= 0) then
         return
       endif

       arr_to_read => omp_allocator%get_block(DIR_C)

       call adios2_declare_io (io, adios_ctx, 'TestRead', ierr)
       if (.not.io%valid) then
          call abort_test("Cannot create ADIOS2 IO")
       endif

       call adios2_open(engine, io, fname, adios2_mode_read, MPI_COMM_SELF, ierr)
       if (.not.engine%valid) then
          call abort_test("Cannot create ADIOS2 engine")
       endif

       call adios2_begin_step(engine, adios2_step_mode_read, ierr)
       call adios2_inquire_variable(adios_var, io, varname, ierr)
       if (.not.adios_var%valid) then
          call abort_test("Cannot fetch ADIOS2 IO")
       endif
       call adios2_set_step_selection(adios_var, 0_8, 1_8, ierr)
       call adios2_get(engine, adios_var, arr_to_read%data, ierr)
       call adios2_end_step(engine, ierr)

       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                if(arr_to_read%data(i, j, k) /= i + j*nx + k*nx*ny) then
                   if (irank == 0) write(stderr, '(a, f8.4, a, f8.4, a, i5, i5, i5)') &
                     'Mismatch between read array(', arr_to_write%data(i, j, k), &
                     ") and expected index (", i + j*nx + k*nx*ny, "at (i,j,k) = ", i, j, k
                   allpass = .false.
                   return ! end test on first issue
                end if
             end do
          end do
       end do
       !
       if (engine%valid) then
          call adios2_close(engine, ierr)
       endif
     end subroutine
end program
