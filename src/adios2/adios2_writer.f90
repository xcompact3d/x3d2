module m_adios2_writer
    use adios2
    use m_base_adios2
    use iso_fortran_env, only: real64, int64
    implicit none

    private
    public :: adios2_writer_t

    !> @brief ADIOS2 writer type
    !> @details Extends `base_adios2_t` and implements writing functionality
    type, extends(base_adios2_t) :: adios2_writer_t
    contains
        procedure, public :: begin_step_adios2 => begin_step_writer  !! Begins a step in the ADIOS2 engine

        ! generic interface for write data
        generic, public :: write_data => write_scalar_real, write_array_2d_real

        ! specific interfaces for writing data
        procedure, private :: write_scalar_real    !! Writes scalar single precision real data
        procedure, private :: write_array_2d_real  !! Writes 2d array single precision real data
    end type adios2_writer_t

contains

    subroutine begin_step_writer(self)
        class(adios2_writer_t), intent(inout) :: self

        if (self%is_step_active) return

        call adios2_begin_step(self%engine, adios2_step_mode_append, self%ierr)
        call self%handle_error(self%ierr, "Error beginning  ADIOS2 step")
        self%is_step_active = .true.  ! set flag after successful step begin
    end subroutine begin_step_writer

    subroutine write_scalar_real(self, name, data, shape_dims, start_dims, count_dims)
        class(adios2_writer_t), intent(inout) :: self
        character(len=*), intent(in) :: name  !! unique variable identifier within io
        real(kind=real64), intent(in) :: data
        integer(kind=int64), dimension(:), optional, intent(in) :: shape_dims, start_dims, count_dims
        type(adios2_variable) :: var                   !! hanlder to newly defined variable

        ! define adios2 variable to be written in given file format
        call adios2_define_variable(var, self%io, name, &
                                    adios2_type_dp, self%ierr)
        call self%handle_error(self%ierr, "Error defining ADIOS2 scalar single precision real variable")

        call adios2_put(self%engine, var, data, adios2_mode_deferred, self%ierr)
        call self%handle_error(self%ierr, "Error writing ADIOS2 scalar single precision real data")
    end subroutine write_scalar_real

    subroutine write_array_2d_real(self, name, data, shape_dims, start_dims, count_dims)
        class(adios2_writer_t), intent(inout) :: self
        character(len=*), intent(in) :: name
        real(kind=real64), dimension(:,:), intent(in) :: data
        integer(kind=int64), dimension(ndims_2d), intent(in) :: shape_dims, start_dims, count_dims
        type(adios2_variable) :: var

        print*, "writing ADIOS2 2d array"
        call adios2_define_variable(var, self%io, name, adios2_type_dp, ndims_2d, &
                                    shape_dims, start_dims, count_dims, &
                                    adios2_constant_dims, self%ierr)
        call self%handle_error(self%ierr, "Error defining ADIOS2 2D array single precision real variable")

        call adios2_put(self%engine, var, data, adios2_mode_deferred, self%ierr)
        call self%handle_error(self%ierr, "Error writing ADIOS2 2D array single precision real data")
    end subroutine write_array_2d_real

end module m_adios2_writer
