module m_omptgt_kernels_fieldops
  !! Field operation kernels for the OpenMP target backend, following
  !! cuda/kernels/fieldops.f90 and the sum/reorder kernels of
  !! cuda/kernels/reorder.f90.
  !!
  !! Arrays are explicit-shape or assumed-size dummies named in
  !! `is_device_ptr` (see m_omptgt_common). Pointwise kernels see a field as
  !! a flat array of `n` points, as the CUDA ones sweep the whole block.
  !! Kernels on the pencil layout index a DIR_X/Y/Z field as (SZ, n, *) in
  !! the layout of m_ordering, which is not the one the CUDA kernels use: a
  !! DIR_X group index is b = y_block + (z - 1)*n_yb_pad, y-blocks fastest
  !! with n_yb_pad = ny_padded/SZ of them per z station.
  use m_common, only: dp
  use m_ordering, only: get_index_reordering
  use m_omptgt_common, only: SZ

  implicit none

  private
  public :: copy, axpby, pwmul, field_scale, field_shift
  public :: vorticity_from_gradients, qcriterion_from_gradients, &
            smagorinsky_from_gradients, sgs_stress_from_gradients
  public :: buffer_copy
  public :: scalar_product, scalar_product_cart, vector_norm_squared, &
            field_max_sum, volume_integral, slice_max_sum
  public :: field_set_x_face, field_set_y_face, field_set_y_plane, &
            field_set_abl_wall_stress, field_set_x_face_from_field, &
            field_set_y_face_from_field, field_add_x_face_from_field, &
            field_add_y_face_from_field
  public :: field_plane_partial_sums, field_plane_sums_reduce
  public :: reorder, sum_intox

contains

  ! --------------------------------------------------------------------------
  ! Pointwise kernels
  ! --------------------------------------------------------------------------

  subroutine copy(n, dst, src)
    integer, intent(in) :: n
    real(dp), intent(inout) :: dst(*)
    real(dp), intent(in) :: src(*)

    integer :: i

    !$omp target teams distribute parallel do is_device_ptr(dst, src)
    do i = 1, n
      dst(i) = src(i)
    end do
    !$omp end target teams distribute parallel do

  end subroutine copy

  subroutine axpby(n, alpha, x, beta, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(in) :: x(*)
    real(dp), intent(inout) :: y(*)

    integer :: i

    !$omp target teams distribute parallel do is_device_ptr(x, y)
    do i = 1, n
      y(i) = alpha*x(i) + beta*y(i)
    end do
    !$omp end target teams distribute parallel do

  end subroutine axpby

  subroutine pwmul(y, x, n)
    real(dp), intent(inout) :: y(*)
    real(dp), intent(in) :: x(*)
    integer, intent(in) :: n

    integer :: i

    !$omp target teams distribute parallel do is_device_ptr(x, y)
    do i = 1, n
      y(i) = y(i)*x(i)
    end do
    !$omp end target teams distribute parallel do

  end subroutine pwmul

  subroutine field_scale(f, alpha, n)
    real(dp), intent(inout) :: f(*)
    real(dp), intent(in) :: alpha
    integer, intent(in) :: n

    integer :: i

    !$omp target teams distribute parallel do is_device_ptr(f)
    do i = 1, n
      f(i) = alpha*f(i)
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_scale

  subroutine field_shift(f, const, n)
    real(dp), intent(inout) :: f(*)
    real(dp), intent(in) :: const
    integer, intent(in) :: n

    integer :: i

    !$omp target teams distribute parallel do is_device_ptr(f)
    do i = 1, n
      f(i) = f(i) + const
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_shift

  subroutine vorticity_from_gradients(vort, dudy, dudz, dvdx, dvdz, dwdx, &
                                      dwdy, n)
    real(dp), intent(inout) :: vort(*)
    real(dp), intent(in) :: dudy(*), dudz(*), dvdx(*), dvdz(*), dwdx(*), &
                            dwdy(*)
    integer, intent(in) :: n

    integer :: i

    !$omp target teams distribute parallel do &
    !$omp is_device_ptr(vort, dudy, dudz, dvdx, dvdz, dwdx, dwdy)
    do i = 1, n
      vort(i) = sqrt((dwdy(i) - dvdz(i))**2 + (dudz(i) - dwdx(i))**2 &
                     + (dvdx(i) - dudy(i))**2)
    end do
    !$omp end target teams distribute parallel do

  end subroutine vorticity_from_gradients

  subroutine qcriterion_from_gradients(qcrit, dudx, dudy, dudz, dvdx, dvdy, &
                                       dvdz, dwdx, dwdy, dwdz, n)
    real(dp), intent(inout) :: qcrit(*)
    real(dp), intent(in) :: dudx(*), dudy(*), dudz(*)
    real(dp), intent(in) :: dvdx(*), dvdy(*), dvdz(*)
    real(dp), intent(in) :: dwdx(*), dwdy(*), dwdz(*)
    integer, intent(in) :: n

    integer :: i

    !$omp target teams distribute parallel do &
    !$omp is_device_ptr(qcrit, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, &
    !$omp               dwdy, dwdz)
    do i = 1, n
      qcrit(i) = -0.5_dp*(dudx(i)*dudx(i) + dvdy(i)*dvdy(i) &
                          + dwdz(i)*dwdz(i)) &
                 - dudy(i)*dvdx(i) - dudz(i)*dwdx(i) - dvdz(i)*dwdy(i)
    end do
    !$omp end target teams distribute parallel do

  end subroutine qcriterion_from_gradients

  subroutine smagorinsky_from_gradients(nut, mixing_length_sq, dudx, dudy, &
                                        dudz, dvdx, dvdy, dvdz, dwdx, dwdy, &
                                        dwdz, n)
    real(dp), intent(inout) :: nut(*)
    real(dp), intent(in) :: mixing_length_sq(*)
    real(dp), intent(in) :: dudx(*), dudy(*), dudz(*)
    real(dp), intent(in) :: dvdx(*), dvdy(*), dvdz(*)
    real(dp), intent(in) :: dwdx(*), dwdy(*), dwdz(*)
    integer, intent(in) :: n

    integer :: i

    !$omp target teams distribute parallel do &
    !$omp is_device_ptr(nut, mixing_length_sq, dudx, dudy, dudz, dvdx, dvdy, &
    !$omp               dvdz, dwdx, dwdy, dwdz)
    do i = 1, n
      if (mixing_length_sq(i) > 0._dp) then
        ! nut = l^2*sqrt(2*S_ij*S_ij)
        nut(i) = mixing_length_sq(i) &
                 *sqrt(2._dp*(dudx(i)**2 + dvdy(i)**2 + dwdz(i)**2 &
                              + 0.5_dp*(dudy(i) + dvdx(i))**2 &
                              + 0.5_dp*(dudz(i) + dwdx(i))**2 &
                              + 0.5_dp*(dvdz(i) + dwdy(i))**2))
      else
        nut(i) = 0._dp
      end if
    end do
    !$omp end target teams distribute parallel do

  end subroutine smagorinsky_from_gradients

  subroutine sgs_stress_from_gradients(stress, nut, gradient_a, gradient_b, &
                                       scale_a, scale_b, n)
    real(dp), intent(inout) :: stress(*)
    real(dp), intent(in) :: nut(*), gradient_a(*), gradient_b(*)
    real(dp), intent(in) :: scale_a, scale_b
    integer, intent(in) :: n

    integer :: i

    !$omp target teams distribute parallel do &
    !$omp is_device_ptr(stress, nut, gradient_a, gradient_b)
    do i = 1, n
      if (nut(i) > 0._dp) then
        stress(i) = nut(i)*(scale_a*gradient_a(i) + scale_b*gradient_b(i))
      else
        stress(i) = 0._dp
      end if
    end do
    !$omp end target teams distribute parallel do

  end subroutine sgs_stress_from_gradients

  ! --------------------------------------------------------------------------
  ! Halo buffers
  ! --------------------------------------------------------------------------

  subroutine buffer_copy(u_send_s, u_send_e, u, n_u, n, n_groups)
    !! Copies the first and last four points of every line into the start
    !! and end send buffers. `n` is the line length, `n_u` the padded extent.
    integer, intent(in) :: n_u, n, n_groups
    real(dp), intent(inout) :: u_send_s(SZ, 4, *), u_send_e(SZ, 4, *)
    real(dp), intent(in) :: u(SZ, n_u, *)

    integer :: i, j, b

    !$omp target teams distribute parallel do collapse(3) &
    !$omp is_device_ptr(u_send_s, u_send_e, u)
    do b = 1, n_groups
      do j = 1, 4
        do i = 1, SZ
          u_send_s(i, j, b) = u(i, j, b)
          u_send_e(i, j, b) = u(i, n - 4 + j, b)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine buffer_copy

  ! --------------------------------------------------------------------------
  ! Reductions over the physical points of a pencil field. Lines are stacked
  ! SZ at a time across n_j and repeat across n_i, so the group of line
  ! (b_j, b_i) is b_j + (b_i - 1)*n_groups_pad and only i + (b_j - 1)*SZ <=
  ! n_j are real. n_groups_pad is the padded stacked extent over SZ.
  ! --------------------------------------------------------------------------

  real(dp) function scalar_product(x, y, n_f, n, n_i, n_groups_pad, n_j) &
    result(s)
    integer, intent(in) :: n_f, n, n_i, n_groups_pad, n_j
    real(dp), intent(in) :: x(SZ, n_f, *), y(SZ, n_f, *)

    integer :: i, j, b_i, b_j

    s = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:s) map(tofrom:s) is_device_ptr(x, y)
    do b_j = 1, (n_j - 1)/SZ + 1
      do b_i = 1, n_i
        do i = 1, SZ
          if (i + (b_j - 1)*SZ <= n_j) then
            do j = 1, n
              s = s + x(i, j, b_j + (b_i - 1)*n_groups_pad) &
                  *y(i, j, b_j + (b_i - 1)*n_groups_pad)
            end do
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end function scalar_product

  real(dp) function scalar_product_cart(x, y, n_c, dims) result(s)
    !! Scalar product of two DIR_C fields, whose storage is the padded
    !! Cartesian box n_c with the physical points in 1:dims.
    integer, intent(in) :: n_c(3), dims(3)
    real(dp), intent(in) :: x(*), y(*)

    s = scalar_product_cart_(x, y, n_c(1), n_c(2), dims(1), dims(2), dims(3))

  end function scalar_product_cart

  real(dp) function scalar_product_cart_(x, y, n1, n2, d1, d2, d3) result(s)
    integer, intent(in) :: n1, n2, d1, d2, d3
    real(dp), intent(in) :: x(n1, n2, *), y(n1, n2, *)

    integer :: i, j, k

    s = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:s) map(tofrom:s) is_device_ptr(x, y)
    do k = 1, d3
      do j = 1, d2
        do i = 1, d1
          s = s + x(i, j, k)*y(i, j, k)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end function scalar_product_cart_

  real(dp) function vector_norm_squared(a, b, c, n_f, n, n_i, n_groups_pad, &
                                        n_j) result(s)
    integer, intent(in) :: n_f, n, n_i, n_groups_pad, n_j
    real(dp), intent(in) :: a(SZ, n_f, *), b(SZ, n_f, *), c(SZ, n_f, *)

    integer :: i, j, b_i, b_j

    s = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:s) map(tofrom:s) is_device_ptr(a, b, c)
    do b_j = 1, (n_j - 1)/SZ + 1
      do b_i = 1, n_i
        do i = 1, SZ
          if (i + (b_j - 1)*SZ <= n_j) then
            do j = 1, n
              s = s + a(i, j, b_j + (b_i - 1)*n_groups_pad)**2 &
                  + b(i, j, b_j + (b_i - 1)*n_groups_pad)**2 &
                  + c(i, j, b_j + (b_i - 1)*n_groups_pad)**2
            end do
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end function vector_norm_squared

  subroutine field_max_sum(max_f, sum_f, f, n_f, n, n_i, n_groups_pad, &
                           n_j)
    !! Maximum and sum of |f|.
    real(dp), intent(out) :: max_f, sum_f
    integer, intent(in) :: n_f, n, n_i, n_groups_pad, n_j
    real(dp), intent(in) :: f(SZ, n_f, *)

    integer :: i, j, b_i, b_j

    max_f = 0._dp
    sum_f = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(max:max_f) reduction(+:sum_f) map(tofrom:max_f, sum_f) &
    !$omp is_device_ptr(f)
    do b_j = 1, (n_j - 1)/SZ + 1
      do b_i = 1, n_i
        do i = 1, SZ
          if (i + (b_j - 1)*SZ <= n_j) then
            do j = 1, n
              sum_f = sum_f + abs(f(i, j, b_j + (b_i - 1)*n_groups_pad))
              max_f = max(max_f, &
                          abs(f(i, j, b_j + (b_i - 1)*n_groups_pad)))
            end do
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_max_sum

  real(dp) function volume_integral(f, n_f, n, n_i, n_groups_pad, n_j) &
    result(s)
    integer, intent(in) :: n_f, n, n_i, n_groups_pad, n_j
    real(dp), intent(in) :: f(SZ, n_f, *)

    integer :: i, j, b_i, b_j

    s = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:s) map(tofrom:s) is_device_ptr(f)
    do b_j = 1, (n_j - 1)/SZ + 1
      do b_i = 1, n_i
        do i = 1, SZ
          if (i + (b_j - 1)*SZ <= n_j) then
            do j = 1, n
              s = s + f(i, j, b_j + (b_i - 1)*n_groups_pad)
            end do
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end function volume_integral

  subroutine slice_max_sum(max_f, sum_f, f, n_f, i_slice, n_i, &
                           n_groups_pad, n_j)
    !! Signed maximum and sum of the slice f(:, i_slice, :).
    real(dp), intent(out) :: max_f, sum_f
    integer, intent(in) :: n_f, i_slice, n_i, n_groups_pad, n_j
    real(dp), intent(in) :: f(SZ, n_f, *)

    integer :: i, b_i, b_j

    max_f = -huge(1._dp)
    sum_f = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(max:max_f) reduction(+:sum_f) map(tofrom:max_f, sum_f) &
    !$omp is_device_ptr(f)
    do b_j = 1, (n_j - 1)/SZ + 1
      do b_i = 1, n_i
        do i = 1, SZ
          if (i + (b_j - 1)*SZ <= n_j) then
            sum_f = sum_f + f(i, i_slice, b_j + (b_i - 1)*n_groups_pad)
            max_f = max(max_f, f(i, i_slice, b_j + (b_i - 1)*n_groups_pad))
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine slice_max_sum

  ! --------------------------------------------------------------------------
  ! Faces and planes of DIR_X fields. The y-block holding y is
  ! (y - 1)/SZ + 1, at lane mod(y - 1, SZ) + 1, and its group at z station z
  ! is y_block + (z - 1)*n_yb_pad.
  ! --------------------------------------------------------------------------

  subroutine field_set_x_face(f, c_start, c_end, bc_start, bc_end, &
                              flow_rate_diff, n_f, nx, ny, nz, n_yb_pad)
    !! Sets the x faces: a Dirichlet value c_start at the inlet (pencil index
    !! 1) and the convective outflow du/dt + Uc*du/dx = 0, Uc = c_end, at
    !! the outlet (pencil index nx). BC_NEUMANN (1) leaves a face untouched.
    integer, intent(in) :: n_f, nx, ny, nz, n_yb_pad, bc_start, bc_end
    real(dp), intent(inout) :: f(SZ, n_f, *)
    real(dp), intent(in) :: c_start, c_end, flow_rate_diff

    integer :: i, y_block, z, b

    !$omp target teams distribute parallel do collapse(3) private(b) &
    !$omp is_device_ptr(f)
    do z = 1, nz
      do y_block = 1, (ny - 1)/SZ + 1
        do i = 1, SZ
          ! Lanes past ny in the last y-block are padding
          if (i + (y_block - 1)*SZ <= ny) then
            b = y_block + (z - 1)*n_yb_pad
            if (bc_start == 2) f(i, 1, b) = c_start
            if (bc_end == 2) then
              f(i, nx, b) = f(i, nx, b) &
                            - c_end*(f(i, nx, b) - f(i, nx - 1, b)) &
                            + flow_rate_diff
            end if
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_set_x_face

  subroutine field_set_x_face_from_field(f, f_start, c_end, bc_start, &
                                         bc_end, flow_rate_diff, n_f, nx, &
                                         ny, nz, n_yb_pad)
    !! As field_set_x_face, with the inlet taken from the inlet plane of
    !! f_start rather than a constant.
    integer, intent(in) :: n_f, nx, ny, nz, n_yb_pad, bc_start, bc_end
    real(dp), intent(inout) :: f(SZ, n_f, *)
    real(dp), intent(in) :: f_start(SZ, n_f, *)
    real(dp), intent(in) :: c_end, flow_rate_diff

    integer :: i, y_block, z, b

    !$omp target teams distribute parallel do collapse(3) private(b) &
    !$omp is_device_ptr(f, f_start)
    do z = 1, nz
      do y_block = 1, (ny - 1)/SZ + 1
        do i = 1, SZ
          if (i + (y_block - 1)*SZ <= ny) then
            b = y_block + (z - 1)*n_yb_pad
            if (bc_start == 2) f(i, 1, b) = f_start(i, 1, b)
            if (bc_end == 2) then
              f(i, nx, b) = f(i, nx, b) &
                            - c_end*(f(i, nx, b) - f(i, nx - 1, b)) &
                            + flow_rate_diff
            end if
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_set_x_face_from_field

  subroutine field_add_x_face_from_field(f, g, set_start, set_end, n_f, nx, &
                                         ny, nz, n_yb_pad)
    !! f += g on the x faces flagged by set_start/set_end.
    integer, intent(in) :: n_f, nx, ny, nz, n_yb_pad
    real(dp), intent(inout) :: f(SZ, n_f, *)
    real(dp), intent(in) :: g(SZ, n_f, *)
    logical, intent(in) :: set_start, set_end

    integer :: i, y_block, z, b

    !$omp target teams distribute parallel do collapse(3) private(b) &
    !$omp is_device_ptr(f, g)
    do z = 1, nz
      do y_block = 1, (ny - 1)/SZ + 1
        do i = 1, SZ
          if (i + (y_block - 1)*SZ <= ny) then
            b = y_block + (z - 1)*n_yb_pad
            if (set_start) f(i, 1, b) = f(i, 1, b) + g(i, 1, b)
            if (set_end) f(i, nx, b) = f(i, nx, b) + g(i, nx, b)
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_add_x_face_from_field

  subroutine field_set_y_face(f, c_start, c_end, set_start, set_end, n_f, &
                              nx, ny, nz, n_yb_pad)
    !! Sets the bottom (y = 1) and top (y = ny) faces where flagged, so the
    !! two y faces can differ (ABL: no-slip floor, free-slip lid).
    integer, intent(in) :: n_f, nx, ny, nz, n_yb_pad
    real(dp), intent(inout) :: f(SZ, n_f, *)
    real(dp), intent(in) :: c_start, c_end
    logical, intent(in) :: set_start, set_end

    integer :: j, z, i_top, yb_top

    i_top = mod(ny - 1, SZ) + 1
    yb_top = (ny - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) is_device_ptr(f)
    do z = 1, nz
      do j = 1, nx
        if (set_start) f(1, j, 1 + (z - 1)*n_yb_pad) = c_start
        if (set_end) f(i_top, j, yb_top + (z - 1)*n_yb_pad) = c_end
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_set_y_face

  subroutine field_set_y_face_from_field(f, f_start, n_f, nx, ny, nz, &
                                         n_yb_pad)
    !! Copies the bottom and top y planes of f_start onto those of f.
    integer, intent(in) :: n_f, nx, ny, nz, n_yb_pad
    real(dp), intent(inout) :: f(SZ, n_f, *)
    real(dp), intent(in) :: f_start(SZ, n_f, *)

    integer :: j, z, i_top, yb_top

    i_top = mod(ny - 1, SZ) + 1
    yb_top = (ny - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) &
    !$omp is_device_ptr(f, f_start)
    do z = 1, nz
      do j = 1, nx
        f(1, j, 1 + (z - 1)*n_yb_pad) = f_start(1, j, 1 + (z - 1)*n_yb_pad)
        f(i_top, j, yb_top + (z - 1)*n_yb_pad) = &
          f_start(i_top, j, yb_top + (z - 1)*n_yb_pad)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_set_y_face_from_field

  subroutine field_add_y_face_from_field(f, g, set_start, set_end, n_f, nx, &
                                         ny, nz, n_yb_pad)
    !! f += g on the y faces flagged by set_start/set_end.
    integer, intent(in) :: n_f, nx, ny, nz, n_yb_pad
    real(dp), intent(inout) :: f(SZ, n_f, *)
    real(dp), intent(in) :: g(SZ, n_f, *)
    logical, intent(in) :: set_start, set_end

    integer :: j, z, i_top, yb_top

    i_top = mod(ny - 1, SZ) + 1
    yb_top = (ny - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) is_device_ptr(f, g)
    do z = 1, nz
      do j = 1, nx
        if (set_start) then
          f(1, j, 1 + (z - 1)*n_yb_pad) = f(1, j, 1 + (z - 1)*n_yb_pad) &
                                          + g(1, j, 1 + (z - 1)*n_yb_pad)
        end if
        if (set_end) then
          f(i_top, j, yb_top + (z - 1)*n_yb_pad) = &
            f(i_top, j, yb_top + (z - 1)*n_yb_pad) &
            + g(i_top, j, yb_top + (z - 1)*n_yb_pad)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_add_y_face_from_field

  subroutine field_set_y_plane(f, c, plane, n_f, nx, nz, n_yb_pad)
    !! Sets the y plane `plane` to c.
    integer, intent(in) :: plane, n_f, nx, nz, n_yb_pad
    real(dp), intent(inout) :: f(SZ, n_f, *)
    real(dp), intent(in) :: c

    integer :: j, z, i_p, yb_p

    i_p = mod(plane - 1, SZ) + 1
    yb_p = (plane - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) is_device_ptr(f)
    do z = 1, nz
      do j = 1, nx
        f(i_p, j, yb_p + (z - 1)*n_yb_pad) = c
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_set_y_plane

  subroutine field_set_abl_wall_stress(stress, u, w, drag_coeff, &
                                       sample_plane, stress_plane, n_f, nx, &
                                       nz, n_yb_pad, component)
    !! Drag law tau = drag_coeff*u_c*|u_h| from each column's own sampled
    !! velocity, written into one stress plane (component 1: x, 3: z).
    integer, intent(in) :: sample_plane, stress_plane
    integer, intent(in) :: n_f, nx, nz, n_yb_pad, component
    real(dp), intent(inout) :: stress(SZ, n_f, *)
    real(dp), intent(in) :: u(SZ, n_f, *), w(SZ, n_f, *)
    real(dp), intent(in) :: drag_coeff

    integer :: j, z, i_sa, yb_sa, i_st, yb_st

    i_sa = mod(sample_plane - 1, SZ) + 1
    yb_sa = (sample_plane - 1)/SZ + 1
    i_st = mod(stress_plane - 1, SZ) + 1
    yb_st = (stress_plane - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) &
    !$omp is_device_ptr(stress, u, w)
    do z = 1, nz
      do j = 1, nx
        if (component == 1) then
          stress(i_st, j, yb_st + (z - 1)*n_yb_pad) = &
            drag_coeff*u(i_sa, j, yb_sa + (z - 1)*n_yb_pad) &
            *sqrt(u(i_sa, j, yb_sa + (z - 1)*n_yb_pad)**2 &
                  + w(i_sa, j, yb_sa + (z - 1)*n_yb_pad)**2)
        else
          stress(i_st, j, yb_st + (z - 1)*n_yb_pad) = &
            drag_coeff*w(i_sa, j, yb_sa + (z - 1)*n_yb_pad) &
            *sqrt(u(i_sa, j, yb_sa + (z - 1)*n_yb_pad)**2 &
                  + w(i_sa, j, yb_sa + (z - 1)*n_yb_pad)**2)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_set_abl_wall_stress

  subroutine field_plane_partial_sums(partial, f, n_f, nx, ny, nz, n_yb_pad)
    !! First pass of field_plane_sums: partial(j, z) is the sum over x of y
    !! index j at z, summed in a fixed order.
    integer, intent(in) :: n_f, nx, ny, nz, n_yb_pad
    real(dp), intent(inout) :: partial(ny, *)
    real(dp), intent(in) :: f(SZ, n_f, *)

    integer :: j, z, k

    !$omp target teams distribute parallel do collapse(2) &
    !$omp is_device_ptr(partial, f)
    do z = 1, nz
      do j = 1, ny
        partial(j, z) = 0._dp
        do k = 1, nx
          partial(j, z) = partial(j, z) &
                          + f(mod(j - 1, SZ) + 1, k, &
                              (j - 1)/SZ + 1 + (z - 1)*n_yb_pad)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_plane_partial_sums

  subroutine field_plane_sums_reduce(sums, partial, ny, nz)
    !! Second pass of field_plane_sums: sums(j) is the sum over z of
    !! partial(j, z), in a fixed order.
    integer, intent(in) :: ny, nz
    real(dp), intent(inout) :: sums(*)
    real(dp), intent(in) :: partial(ny, *)

    integer :: j, z

    !$omp target teams distribute parallel do is_device_ptr(sums, partial)
    do j = 1, ny
      sums(j) = 0._dp
      do z = 1, nz
        sums(j) = sums(j) + partial(j, z)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_plane_sums_reduce

  ! --------------------------------------------------------------------------
  ! Reordering. Every point is moved to the slot get_index_reordering gives
  ! it, which covers all twelve layout pairs in one kernel.
  ! --------------------------------------------------------------------------

  subroutine reorder(u_o, n_o, u_i, n_i, dims, dir_from, dir_to, cart_padded)
    !! Reorders u_i, whose layout dir_from has padded extents `dims`, into
    !! u_o in layout dir_to. n_o and n_i are the storage shapes.
    integer, intent(in) :: n_o(3), n_i(3), dims(3)
    integer, intent(in) :: dir_from, dir_to
    integer, intent(in) :: cart_padded(3)
    real(dp), intent(inout) :: u_o(*)
    real(dp), intent(in) :: u_i(*)

    ! Scalar copies enter the region as kernel arguments, where the arrays
    ! would each be mapped with a blocking copy per launch.
    call reorder_(u_o, n_o(1), n_o(2), u_i, n_i(1), n_i(2), &
                  dims(1), dims(2), dims(3), dir_from, dir_to, &
                  cart_padded(1), cart_padded(2), cart_padded(3))

  end subroutine reorder

  subroutine reorder_(u_o, n_o1, n_o2, u_i, n_i1, n_i2, d1, d2, d3, &
                      dir_from, dir_to, cp1, cp2, cp3)
    integer, intent(in) :: n_o1, n_o2, n_i1, n_i2, d1, d2, d3
    integer, intent(in) :: dir_from, dir_to, cp1, cp2, cp3
    real(dp), intent(inout) :: u_o(n_o1, n_o2, *)
    real(dp), intent(in) :: u_i(n_i1, n_i2, *)

    integer :: i, j, k

    !$omp target teams distribute parallel do collapse(3) &
    !$omp is_device_ptr(u_o, u_i)
    do k = 1, d3
      do j = 1, d2
        do i = 1, d1
          call reorder_point(u_o, n_o1, n_o2, u_i, n_i1, n_i2, i, j, k, &
                             dir_from, dir_to, cp1, cp2, cp3)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine reorder_

  subroutine reorder_point(u_o, n_o1, n_o2, u_i, n_i1, n_i2, i, j, k, &
                           dir_from, dir_to, cp1, cp2, cp3)
    !$omp declare target
    !! Moves one point into its reordered slot.
    !!
    !! The reordered indices are locals of this routine rather than scalars
    !! privatised on the loop: NVHPC miscompiles the private-clause spelling,
    !! handing get_index_reordering an invalid address for its intent(out)
    !! arguments.
    integer, intent(in) :: n_o1, n_o2, n_i1, n_i2
    real(dp), intent(inout) :: u_o(n_o1, n_o2, *)
    real(dp), intent(in) :: u_i(n_i1, n_i2, *)
    integer, intent(in) :: i, j, k
    integer, intent(in) :: dir_from, dir_to, cp1, cp2, cp3

    integer :: out_i, out_j, out_k

    call get_index_reordering(out_i, out_j, out_k, i, j, k, &
                              dir_from, dir_to, SZ, [cp1, cp2, cp3])
    u_o(out_i, out_j, out_k) = u_i(i, j, k)

  end subroutine reorder_point

  subroutine sum_intox(u, n_u, u_d, n_d, dims, dir_to, cart_padded)
    !! u += u_d, where u is DIR_X with padded extents `dims` and u_d is the
    !! same field in layout dir_to (DIR_Y or DIR_Z).
    integer, intent(in) :: n_u(3), n_d(3), dims(3)
    integer, intent(in) :: dir_to
    integer, intent(in) :: cart_padded(3)
    real(dp), intent(inout) :: u(*)
    real(dp), intent(in) :: u_d(*)

    call sum_intox_(u, n_u(1), n_u(2), u_d, n_d(1), n_d(2), &
                    dims(1), dims(2), dims(3), dir_to, &
                    cart_padded(1), cart_padded(2), cart_padded(3))

  end subroutine sum_intox

  subroutine sum_intox_(u, n_u1, n_u2, u_d, n_d1, n_d2, d1, d2, d3, &
                        dir_to, cp1, cp2, cp3)
    integer, intent(in) :: n_u1, n_u2, n_d1, n_d2, d1, d2, d3
    integer, intent(in) :: dir_to, cp1, cp2, cp3
    real(dp), intent(inout) :: u(n_u1, n_u2, *)
    real(dp), intent(in) :: u_d(n_d1, n_d2, *)

    integer :: i, j, k

    !$omp target teams distribute parallel do collapse(3) &
    !$omp is_device_ptr(u, u_d)
    do k = 1, d3
      do j = 1, d2
        do i = 1, d1
          call sum_point(u, n_u1, n_u2, u_d, n_d1, n_d2, i, j, k, &
                         dir_to, cp1, cp2, cp3)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine sum_intox_

  subroutine sum_point(u, n_u1, n_u2, u_d, n_d1, n_d2, i, j, k, dir_to, &
                       cp1, cp2, cp3)
    !$omp declare target
    !! Adds to one DIR_X point its counterpart in layout dir_to.
    use m_common, only: DIR_X

    integer, intent(in) :: n_u1, n_u2, n_d1, n_d2
    real(dp), intent(inout) :: u(n_u1, n_u2, *)
    real(dp), intent(in) :: u_d(n_d1, n_d2, *)
    integer, intent(in) :: i, j, k
    integer, intent(in) :: dir_to, cp1, cp2, cp3

    integer :: d_i, d_j, d_k

    call get_index_reordering(d_i, d_j, d_k, i, j, k, DIR_X, dir_to, SZ, &
                              [cp1, cp2, cp3])
    u(i, j, k) = u(i, j, k) + u_d(d_i, d_j, d_k)

  end subroutine sum_point

end module m_omptgt_kernels_fieldops
