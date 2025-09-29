module m_io_field_utils
!! Common utilities for field I/O operations
  use m_common, only: dp, i8, DIR_C
  use m_field, only: field_t
  use m_solver, only: solver_t
  use m_io_base, only: io_file_t, io_writer_t

  implicit none

  private
  public :: field_buffer_map_t, field_ptr_t
  public :: stride_data, stride_data_to_buffer, get_output_dimensions
  public :: generate_coordinates

  type :: field_buffer_map_t
    !! Race-free field buffer mapping for async I/O operations.
    !! Each field gets its own dedicated buffer to prevent data races
    !! when multiple async write operations are in flight.
    character(len=32) :: field_name
    real(dp), dimension(:, :, :), allocatable :: buffer
  end type field_buffer_map_t

  type :: field_ptr_t
    class(field_t), pointer :: ptr => null()
  end type field_ptr_t

contains

  function stride_data( &
    input_data, dims, stride, output_dims_out) &
    result(output_data)
    !! stride the input data based on the specified stride
    real(dp), dimension(:, :, :), intent(in) :: input_data
    integer, dimension(3), intent(in) :: dims
    integer, dimension(3), intent(in) :: stride
    integer, dimension(3), intent(out) :: output_dims_out

    real(dp), dimension(:, :, :), allocatable :: output_data
    integer :: i_stride, j_stride, k_stride
    integer :: i_max, j_max, k_max

    if (all(stride == 1)) then
      allocate (output_data(dims(1), dims(2), dims(3)))
      output_data = input_data
      output_dims_out = dims
      return
    end if

    i_stride = stride(1); j_stride = stride(2); k_stride = stride(3)

    i_max = (dims(1) - 1)/i_stride + 1
    j_max = (dims(2) - 1)/j_stride + 1
    k_max = (dims(3) - 1)/k_stride + 1

    output_dims_out = [i_max, j_max, k_max]
    allocate (output_data(i_max, j_max, k_max))

    output_data = input_data(1:dims(1):i_stride, &
                             1:dims(2):j_stride, 1:dims(3):k_stride)
  end function stride_data

  subroutine stride_data_to_buffer( &
    input_data, dims, stride, out_buffer, output_dims_out &
    )
    real(dp), dimension(:, :, :), intent(in) :: input_data
    integer, dimension(3), intent(in) :: dims
    integer, dimension(3), intent(in) :: stride
    real(dp), dimension(:, :, :), allocatable, intent(inout) :: out_buffer
    integer, dimension(3), intent(out) :: output_dims_out

    integer :: i_stride, j_stride, k_stride
    integer :: i_max, j_max, k_max

    if (all(stride == 1)) then
      if (allocated(out_buffer)) then
        if (size(out_buffer, 1) /= dims(1) &
            .or. size(out_buffer, 2) /= dims(2) .or. &
            size(out_buffer, 3) /= dims(3)) then
          deallocate (out_buffer)
          allocate (out_buffer(dims(1), dims(2), dims(3)))
        end if
      else
        allocate (out_buffer(dims(1), dims(2), dims(3)))
      end if
      out_buffer = input_data
      output_dims_out = dims
      return
    end if

    i_stride = stride(1); j_stride = stride(2); k_stride = stride(3)

    i_max = (dims(1) + i_stride - 1)/i_stride
    j_max = (dims(2) + j_stride - 1)/j_stride
    k_max = (dims(3) + k_stride - 1)/k_stride

    output_dims_out = [i_max, j_max, k_max]

    if (allocated(out_buffer)) then
      if (size(out_buffer, 1) /= i_max &
          .or. size(out_buffer, 2) /= j_max .or. &
          size(out_buffer, 3) /= k_max) then
        deallocate (out_buffer)
        allocate (out_buffer(i_max, j_max, k_max))
      end if
    else
      allocate (out_buffer(i_max, j_max, k_max))
    end if

    out_buffer = input_data(1:dims(1):i_stride, &
                            1:dims(2):j_stride, 1:dims(3):k_stride)
  end subroutine stride_data_to_buffer

  subroutine get_output_dimensions( &
    shape_dims, start_dims, count_dims, stride_factors, &
    output_shape, output_start, output_count, output_dims_local, &
    last_shape_dims, last_stride_factors, last_output_shape)

    integer(i8), dimension(3), intent(in) :: shape_dims, start_dims, count_dims
    integer, dimension(3), intent(in) :: stride_factors
    integer(i8), dimension(3), intent(out) :: output_shape, output_start
    integer(i8), dimension(3), intent(out) :: output_count
    integer, dimension(3), intent(out) :: output_dims_local
    ! Cache parameters - optional for performance optimization
    integer(i8), dimension(3), intent(inout), optional :: last_shape_dims
    integer, dimension(3), intent(inout), optional :: last_stride_factors
    integer(i8), dimension(3), intent(inout), optional :: last_output_shape

    if (all(stride_factors == 1)) then
      output_shape = shape_dims
      output_start = start_dims
      output_count = count_dims
      output_dims_local = int(count_dims)
      return
    end if

    ! Use cache if available and valid
    if (present(last_shape_dims) .and. present(last_stride_factors) .and. &
        present(last_output_shape)) then
      if (all(shape_dims == last_shape_dims) .and. &
          all(stride_factors == last_stride_factors) .and. &
          all(last_output_shape > 0)) then
        output_shape = last_output_shape
      else
        output_shape = [(shape_dims(1) + stride_factors(1) - 1_i8) &
                        /int(stride_factors(1), i8), &
                        (shape_dims(2) + stride_factors(2) - 1_i8) &
                        /int(stride_factors(2), i8), &
                        (shape_dims(3) + stride_factors(3) - 1_i8) &
                        /int(stride_factors(3), i8)]

        last_shape_dims = shape_dims
        last_stride_factors = stride_factors
        last_output_shape = output_shape
      end if
    else
      ! No cache available, compute directly
      output_shape = [(shape_dims(1) + stride_factors(1) - 1_i8) &
                      /int(stride_factors(1), i8), &
                      (shape_dims(2) + stride_factors(2) - 1_i8) &
                      /int(stride_factors(2), i8), &
                      (shape_dims(3) + stride_factors(3) - 1_i8) &
                      /int(stride_factors(3), i8)]
    end if

    output_start = [start_dims(1)/int(stride_factors(1), i8), &
                    start_dims(2)/int(stride_factors(2), i8), &
                    start_dims(3)/int(stride_factors(3), i8)]

    output_dims_local = [(int(count_dims(1)) + stride_factors(1) - 1) &
                         /stride_factors(1), &
                         (int(count_dims(2)) + stride_factors(2) - 1) &
                         /stride_factors(2), &
                         (int(count_dims(3)) + stride_factors(3) - 1) &
                         /stride_factors(3)]

    output_count = [int(output_dims_local(1), i8), &
                    int(output_dims_local(2), i8), &
                    int(output_dims_local(3), i8)]
  end subroutine get_output_dimensions

  subroutine generate_coordinates( &
    solver, writer, file, shape_dims, start_dims, count_dims, data_loc, &
    coords_x, coords_y, coords_z &
    )
    class(solver_t), intent(in) :: solver
    class(io_writer_t), intent(inout) :: writer
    class(io_file_t), intent(inout) :: file
    integer(i8), dimension(3), intent(in) :: shape_dims, start_dims, count_dims
    integer, intent(in) :: data_loc
    real(dp), dimension(:, :, :), allocatable, intent(inout) :: coords_x, coords_y, coords_z

    integer :: i, nx, ny, nz
    real(dp), dimension(3) :: coords
    integer(i8), dimension(3) :: x_shape, y_shape, z_shape
    integer(i8), dimension(3) :: x_start, y_start, z_start
    integer(i8), dimension(3) :: x_count, y_count, z_count

    nx = int(count_dims(1))
    ny = int(count_dims(2))
    nz = int(count_dims(3))

    ! coordinates are structured as 3D arrays for ParaView ADIOS2 reader compatibility
    if (.not. allocated(coords_x) .or. size(coords_x) /= nx) then
      if (allocated(coords_x)) deallocate (coords_x)
      allocate (coords_x(nx, 1, 1))
    end if

    if (.not. allocated(coords_y) .or. size(coords_y) /= ny) then
      if (allocated(coords_y)) deallocate (coords_y)
      allocate (coords_y(1, ny, 1))
    end if

    if (.not. allocated(coords_z) .or. size(coords_z) /= nz) then
      if (allocated(coords_z)) deallocate (coords_z)
      allocate (coords_z(1, 1, nz))
    end if

    do i = 1, nx
      coords = solver%mesh%get_coordinates(i, 1, 1, data_loc)
      coords_x(i, 1, 1) = coords(1)
    end do

    do i = 1, ny
      coords = solver%mesh%get_coordinates(1, i, 1, data_loc)
      coords_y(1, i, 1) = coords(2)
    end do

    do i = 1, nz
      coords = solver%mesh%get_coordinates(1, 1, i, data_loc)
      coords_z(1, 1, i) = coords(3)
    end do

    x_shape = [1_i8, 1_i8, shape_dims(1)]
    x_start = [0_i8, 0_i8, start_dims(1)]
    x_count = [1_i8, 1_i8, count_dims(1)]

    y_shape = [1_i8, shape_dims(2), 1_i8]
    y_start = [0_i8, start_dims(2), 0_i8]
    y_count = [1_i8, count_dims(2), 1_i8]

    z_shape = [shape_dims(3), 1_i8, 1_i8]
    z_start = [start_dims(3), 0_i8, 0_i8]
    z_count = [count_dims(3), 1_i8, 1_i8]

    call writer%write_data( &
      "coordinates/x", coords_x, file, x_shape, x_start, x_count &
      )
    call writer%write_data( &
      "coordinates/y", coords_y, file, y_shape, y_start, y_count &
      )
    call writer%write_data( &
      "coordinates/z", coords_z, file, z_shape, z_start, z_count &
      )
  end subroutine generate_coordinates

end module m_io_field_utils
