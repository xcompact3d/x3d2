program test_omp_io
   use iso_fortran_env, only: stderr => error_unit
   use adios2
   use mpi

   use m_allocator, only: allocator_t, field_t
   use m_common, only: dp, globs_t
   use m_omp_common, only: SZ

   implicit none

   logical :: allpass !! flag indicating all tests pass

   class(field_t), pointer :: arr !! array to save & restore

   !> MPI vars
   integer :: nrank, nproc
   !> Error code, used for MPI and ADIOS
   integer :: ierr

   !> ADIOS2 context, io and variable
   integer, parameter :: NSTEPS = 16

   integer :: n, n_block, n_glob
   integer :: nx, ny, nz
   integer :: nx_loc, ny_loc, nz_loc

   integer :: i
   character(len=80)::fmt

   allpass = .true.

   ! nx = 96
   ! ny = 96
   ! nz = 96
   !
   ! nx_loc = nx/nproc
   ! ny_loc = ny/nproc
   ! nz_loc = nz/nproc

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

   !================ 
   ! SETUP TEST DATA
   !================

   ! omp_allocator = allocator_t([SZ, globs%nx_loc, globs%n_groups_x])
   ! allocator => omp_allocator
   ! print*, 'OpenMP allocator instantiated'
   !
   ! if (nrank == 0 && nproc > 1) print*, 'Parallel run with', nproc, 'ranks'
   !
   ! arr => omp_allocator%get_block()
   !
   ! n_glob = globs%nx
   ! n = n_glob/nproc
   ! n_block = globs%n_groups_x 
   !
   ! ! Initialise with a simple index (value of index not important)
   ! do k = 1, n_block
   !    do j = 1, n
   !       do i = 1, SZ
   !          arr%data(i, j, k) = i + j*SZ + k*SZ*n
   !       end do
   !    end do
   ! end do

   !================ 
   ! WRITE TEST DATA
   !================

   call write_test_file()

   !================ 
   ! READ TEST DATA
   !===============

   ! call read_test_file()

   if (allpass) then
      if (nrank == 0) write(stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
   else
      error stop 'SOME TESTS FAILED.'
   end if

   call MPI_Finalize(ierr)

   contains

     subroutine write_test_file()
       type(adios2_adios) :: adios_ctx
       type(adios2_io) :: io
       type(adios2_engine) :: engine
       type(adios2_variable) :: adios_var

       integer :: ierr

       call adios2_init(adios_ctx, MPI_COMM_WORLD, ierr)
       if (.not.adios_ctx%valid) then
          write(*,*) "Error when creating ADIOS2 context"
       endif

       call adios2_declare_io (io, adios_ctx, 'TestOutput', ierr)
       if (.not.io%valid) then
          write(*,*) "Error when creating ADIOS2 io"
       endif

       call adios2_open(engine, io, "test_omp_io.bp", adios2_mode_write, ierr)
       if (.not.engine%valid) then
          write(*,*) "Error when creating ADIOS2 engine"
       endif

       ! call adios2_define_variable(adios_var, io, "test", adios2_type_dp, &
       !                             3, [globs%nx, globs%ny, globs%nz], [SZ, globs%nx_loc, globs%n_groups_x], ierr)

       call adios2_define_variable(adios_var, io, "test", adios2_type_integer4, ierr)

       do i=0, NSTEPS-1
          call adios2_begin_step(engine, adios2_step_mode_append, ierr)
          ! call adios2_put(engine, "test", arr%data, ierr)
          call adios2_put(engine, "test", 42+i, ierr)
          call adios2_end_step(engine, ierr)
       enddo

       if (engine%valid) then
          call adios2_close(engine, ierr)
       endif

       if (adios_ctx%valid) then
         call adios2_finalize(adios_ctx, ierr)
       endif
     end subroutine

     subroutine read_test_file()
       type(adios2_adios) :: adios_ctx
       type(adios2_io) :: io
       type(adios2_engine) :: engine
       type(adios2_variable) :: adios_var
       integer*8 :: numsteps
       integer, dimension(:), allocatable :: steps

       integer :: ierr

       call adios2_init(adios_ctx, MPI_COMM_WORLD, ierr)
       if (.not.adios_ctx%valid) then
          write(*,*) "Error when creating ADIOS2 context"
       endif

       call adios2_declare_io (io, adios_ctx, 'TestInput', ierr)
       if (.not.io%valid) then
          write(*,*) "Error when creating ADIOS2 io"
       endif

       call adios2_open(engine, io, "test_omp_io.bp", adios2_mode_read, ierr)
       if (.not.engine%valid) then
          write(*,*) "Error when creating ADIOS2 engine"
       endif

       call adios2_inquire_variable(adios_var, io, "test", ierr)
       if (.not.adios_var%valid) then
          write(*,*) "ierr = ", ierr
          allpass = .false.
          error stop "Error when fetching adios variable"
       endif

       call adios2_variable_steps(numsteps, adios_var, ierr)
       if(numsteps /= NSTEPS) then
          if (nrank == 0) write(stderr, '(a)') 'numsteps /= NSTEPS.'
          if (nrank == 0) write(stderr, '(a, i5)') 'numsteps =', numsteps
          allpass = .false.
          return
       endif
       ! call adios2_set_step_selection(adios_var, 0_8, numsteps, ierr)
       ! allocate(steps(numsteps))
       ! call adios2_get(engine, adios_var, steps , ierr)
       ! write(*,'(a,i5)') "Number of steps in diag.bp = ", numsteps
       ! write(fmt,'(a,i5,a)') '(a10,',numsteps,'i4,a2)'
       ! !write(*,'(a)') fmt
       ! write(*,fmt) "Steps = [", steps, " ]"
       ! deallocate(steps)
       !
       if (engine%valid) then
          call adios2_close(engine, ierr)
       endif

       if (adios_ctx%valid) then
         call adios2_finalize(adios_ctx, ierr)
       endif
     end subroutine
end program
