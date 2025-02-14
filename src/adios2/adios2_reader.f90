module m_adios2_reader
    use adios2
    use m_base_adios2
    use iso_fortran_env, only: real64, int64
    implicit none

    private
    public :: adios2_reader_t

    !> @brief ADIOS2 reader type
    !> @details Extends `base_adios2_t` and implements reading functionality
    type, extends(base_adios2_t) :: adios2_reader_t
    contains
        procedure, public :: begin_step_adios2 => begin_step_reader  !! Begins a step in the ADIOS2 engine

        ! generic interface for write data
        generic, public :: read_data => read_scalar_real, read_array2d_real

        ! specific interfaces for writing data
        procedure, private :: read_scalar_real      !! Reads scalar single precision real data
        procedure, private :: read_array2d_real     !! Reads 2d array single precision real data
    end type adios2_reader_t

contains

    subroutine begin_step_reader(self)
        class(adios2_reader_t), intent(inout) :: self

        if (self%is_step_active) return

        call adios2_begin_step(self%engine, adios2_step_mode_read, self%ierr)
        call self%handle_error(self%ierr, "Error beginning  ADIOS2 step")
        self%is_step_active = .true.  ! set flag after successful step begin
    end subroutine begin_step_reader

    subroutine read_scalar_real(self, name, data)
        class(adios2_reader_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=real64), intent(out) :: data
        type(adios2_variable) :: var

        ! retrieve a variable hanlder within current io handler
        call adios2_inquire_variable(var, self%io, name, self%ierr)
        call self%handle_error(self%ierr, "Failed to inquire variable"//name)

        if (self%ierr == adios2_found) then
            call adios2_get(self%engine, var, data, adios2_mode_sync, self%ierr)
            call self%handle_error(self%ierr, "Failed to read variable"//name)
        end if
    end subroutine read_scalar_real

    subroutine read_array2d_real(self, name, data, start_dims, count_dims)
        class(adios2_reader_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=real64), dimension(:,:), allocatable, intent(out) :: data
        integer(kind=int64), dimension(ndims_2d), intent(in) :: start_dims, count_dims
        type(adios2_variable) :: var

        ! retrieve a variable hanlder within current io handler
        call adios2_inquire_variable(var, self%io, name, self%ierr)
        call self%handle_error(self%ierr, "Failed to inquire variable"//name)

        if (self%ierr == adios2_found) then
            ! allocate data array
            if (.not. allocated(data)) allocate(data(count_dims(1), count_dims(2)))

            ! set selection for reading portion of data
            call adios2_set_selection(var, ndims_2d, start_dims, count_dims, self%ierr)
            call self%handle_error(self%ierr, "Failed to set selection for variable"//name)

            ! read data
            call adios2_get(self%engine, var, data, adios2_mode_sync, self%ierr)
            call self%handle_error(self%ierr, "Failed to read variable"//name)
        end if
    end subroutine read_array2d_real

end module m_adios2_reader
