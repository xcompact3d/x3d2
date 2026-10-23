module m_omptgt_kernels_fieldops
  !! Field operation kernels for the OpenMP target backend, following
  !! cuda/kernels/fieldops.f90 and the sum/reorder kernels of
  !! cuda/kernels/reorder.f90.
  !!
  !! Arrays are explicit-shape or assumed-size dummies named in
  !! `OMP_DEVICE_ADDR` (see m_omptgt_common). Pointwise kernels see a field as
  !! a flat array of `n` points, as the CUDA ones sweep the whole block.
  !! Kernels on the pencil layout index a DIR_X/Y/Z field as (SZ, n, groups) in
  !! the layout of m_ordering, which is not the one the CUDA kernels use: a
  !! DIR_X group index is b = y_block + (z - 1)*n_yb_pad, y-blocks fastest
  !! with n_yb_pad = ny_padded/SZ of them per z station.
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z
  use m_omptgt_common, only: SZ

  implicit none

  private
  public :: copy, axpby, pwmul, field_scale, field_shift
  public :: sample_points, add_gaussian_sources
  public :: vorticity_from_gradients, qcriterion_from_gradients, &
            smagorinsky_from_gradients, sgs_stress_from_gradients
  public :: buffer_copy
  public :: scalar_product, scalar_product_cart, vector_norm_squared, &
            field_max_sum, volume_integral, slice_max_sum
  public :: labelled_sums, labelled_axpy
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
    real(dp), intent(inout) :: dst(n)
    real(dp), intent(in) :: src(n)

    integer :: i

    !$omp target teams distribute parallel do OMP_DEVICE_ADDR(dst, src)
    do i = 1, n
      dst(i) = src(i)
    end do
    !$omp end target teams distribute parallel do

  end subroutine copy

  subroutine axpby(n, alpha, x, beta, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(in) :: x(n)
    real(dp), intent(inout) :: y(n)

    integer :: i

    !$omp target teams distribute parallel do OMP_DEVICE_ADDR(x, y)
    do i = 1, n
      y(i) = alpha*x(i) + beta*y(i)
    end do
    !$omp end target teams distribute parallel do

  end subroutine axpby

  subroutine pwmul(y, x, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: y(n)
    real(dp), intent(in) :: x(n)

    integer :: i

    !$omp target teams distribute parallel do OMP_DEVICE_ADDR(x, y)
    do i = 1, n
      y(i) = y(i)*x(i)
    end do
    !$omp end target teams distribute parallel do

  end subroutine pwmul

  subroutine field_scale(f, alpha, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: f(n)
    real(dp), intent(in) :: alpha

    integer :: i

    !$omp target teams distribute parallel do OMP_DEVICE_ADDR(f)
    do i = 1, n
      f(i) = alpha*f(i)
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_scale

  subroutine field_shift(f, const, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: f(n)
    real(dp), intent(in) :: const

    integer :: i

    !$omp target teams distribute parallel do OMP_DEVICE_ADDR(f)
    do i = 1, n
      f(i) = f(i) + const
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_shift

  subroutine vorticity_from_gradients(vort, dudy, dudz, dvdx, dvdz, dwdx, &
                                      dwdy, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: vort(n)
    real(dp), intent(in) :: dudy(n), dudz(n), dvdx(n), dvdz(n), dwdx(n), &
                            dwdy(n)

    integer :: i

    !$omp target teams distribute parallel do &
    !$omp OMP_DEVICE_ADDR(vort, dudy, dudz, dvdx, dvdz, dwdx, dwdy)
    do i = 1, n
      vort(i) = sqrt((dwdy(i) - dvdz(i))**2 + (dudz(i) - dwdx(i))**2 &
                     + (dvdx(i) - dudy(i))**2)
    end do
    !$omp end target teams distribute parallel do

  end subroutine vorticity_from_gradients

  subroutine qcriterion_from_gradients(qcrit, dudx, dudy, dudz, dvdx, dvdy, &
                                       dvdz, dwdx, dwdy, dwdz, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: qcrit(n)
    real(dp), intent(in) :: dudx(n), dudy(n), dudz(n)
    real(dp), intent(in) :: dvdx(n), dvdy(n), dvdz(n)
    real(dp), intent(in) :: dwdx(n), dwdy(n), dwdz(n)

    integer :: i

    !$omp target teams distribute parallel do &
    !$omp OMP_DEVICE_ADDR(qcrit, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, &
    !$omp                 dwdy, dwdz)
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
    integer, intent(in) :: n
    real(dp), intent(inout) :: nut(n)
    real(dp), intent(in) :: mixing_length_sq(n)
    real(dp), intent(in) :: dudx(n), dudy(n), dudz(n)
    real(dp), intent(in) :: dvdx(n), dvdy(n), dvdz(n)
    real(dp), intent(in) :: dwdx(n), dwdy(n), dwdz(n)

    integer :: i

    !$omp target teams distribute parallel do &
    !$omp OMP_DEVICE_ADDR(nut, mixing_length_sq, dudx, dudy, dudz, dvdx, &
    !$omp                 dvdy, dvdz, dwdx, dwdy, dwdz)
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
    integer, intent(in) :: n
    real(dp), intent(inout) :: stress(n)
    real(dp), intent(in) :: nut(n), gradient_a(n), gradient_b(n)
    real(dp), intent(in) :: scale_a, scale_b

    integer :: i

    !$omp target teams distribute parallel do &
    !$omp OMP_DEVICE_ADDR(stress, nut, gradient_a, gradient_b)
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
    real(dp), intent(inout) :: u_send_s(SZ, 4, n_groups), &
                               u_send_e(SZ, 4, n_groups)
    real(dp), intent(in) :: u(SZ, n_u, n_groups)

    integer :: i, j, b

    !$omp target teams distribute parallel do collapse(3) &
    !$omp OMP_DEVICE_ADDR(u_send_s, u_send_e, u)
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
    real(dp), intent(in) :: x(SZ, n_f, n_groups_pad*n_i), &
                            y(SZ, n_f, n_groups_pad*n_i)

    integer :: i, j, b_i, b_j

    s = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:s) map(tofrom:s) OMP_DEVICE_ADDR(x, y)
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

    s = scalar_product_cart_(x, y, n_c(1), n_c(2), n_c(3), dims(1), dims(2), &
                             dims(3))

  end function scalar_product_cart

  real(dp) function scalar_product_cart_(x, y, n1, n2, n3, d1, d2, d3) &
    result(s)
    integer, intent(in) :: n1, n2, n3, d1, d2, d3
    real(dp), intent(in) :: x(n1, n2, n3), y(n1, n2, n3)

    integer :: i, j, k

    s = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:s) map(tofrom:s) OMP_DEVICE_ADDR(x, y)
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
    real(dp), intent(in) :: a(SZ, n_f, n_groups_pad*n_i), &
                            b(SZ, n_f, n_groups_pad*n_i), &
                            c(SZ, n_f, n_groups_pad*n_i)

    integer :: i, j, b_i, b_j

    s = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:s) map(tofrom:s) OMP_DEVICE_ADDR(a, b, c)
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

  subroutine labelled_sums(s, g, label, a, b, c, n_label, n_f, n, n_i, &
                           n_groups_pad, n_j)
    !! s(1:3, l) += sums of g*a, g*b and g*c over the points labelled l, as
    !! cuda/kernels/fieldops.f90: each thread accumulates runs of equal labels
    !! along its pencil and flushes a run when the label changes, so the
    !! atomics scale with the label changes rather than the points. Atomics
    !! instead of reduction(+:s), which flang cannot do on the GPU for an
    !! array whose shape is only known at run time.
    integer, intent(in) :: n_label, n_f, n, n_i, n_groups_pad, n_j
    real(dp), intent(inout) :: s(3, n_label)
    real(dp), intent(in) :: g(SZ, n_f, n_groups_pad*n_i), &
                            label(SZ, n_f, n_groups_pad*n_i), &
                            a(SZ, n_f, n_groups_pad*n_i), &
                            b(SZ, n_f, n_groups_pad*n_i), &
                            c(SZ, n_f, n_groups_pad*n_i)

    real(dp) :: sum_a, sum_b, sum_c
    integer :: i, j, k, l, l_run, b_i, b_j

    !$omp target teams distribute parallel do collapse(3) &
    !$omp private(k, l, l_run, sum_a, sum_b, sum_c) &
    !$omp OMP_DEVICE_ADDR(s, g, label, a, b, c)
    do b_j = 1, (n_j - 1)/SZ + 1
      do b_i = 1, n_i
        do i = 1, SZ
          if (i + (b_j - 1)*SZ <= n_j) then
            k = b_j + (b_i - 1)*n_groups_pad
            l_run = 0
            sum_a = 0._dp
            sum_b = 0._dp
            sum_c = 0._dp
            do j = 1, n
              l = nint(label(i, j, k))
              if (l < 1 .or. l > n_label) l = 0
              if (l /= l_run) then
                if (l_run /= 0) then
                  !$omp atomic update
                  s(1, l_run) = s(1, l_run) + sum_a
                  !$omp atomic update
                  s(2, l_run) = s(2, l_run) + sum_b
                  !$omp atomic update
                  s(3, l_run) = s(3, l_run) + sum_c
                end if
                l_run = l
                sum_a = 0._dp
                sum_b = 0._dp
                sum_c = 0._dp
              end if
              if (l /= 0) then
                sum_a = sum_a + g(i, j, k)*a(i, j, k)
                sum_b = sum_b + g(i, j, k)*b(i, j, k)
                sum_c = sum_c + g(i, j, k)*c(i, j, k)
              end if
            end do
            if (l_run /= 0) then
              !$omp atomic update
              s(1, l_run) = s(1, l_run) + sum_a
              !$omp atomic update
              s(2, l_run) = s(2, l_run) + sum_b
              !$omp atomic update
              s(3, l_run) = s(3, l_run) + sum_c
            end if
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine labelled_sums

  subroutine labelled_axpy(coeffs, g, label, a, b, c, n_label, n_f, n, n_i, &
                           n_groups_pad, n_j)
    !! a, b, c += coeffs(1:3, l)*g at the points labelled l.
    integer, intent(in) :: n_label, n_f, n, n_i, n_groups_pad, n_j
    real(dp), intent(in) :: coeffs(3, n_label)
    real(dp), intent(in) :: g(SZ, n_f, n_groups_pad*n_i), &
                            label(SZ, n_f, n_groups_pad*n_i)
    real(dp), intent(inout) :: a(SZ, n_f, n_groups_pad*n_i), &
                               b(SZ, n_f, n_groups_pad*n_i), &
                               c(SZ, n_f, n_groups_pad*n_i)

    integer :: i, j, k, l, b_i, b_j

    !$omp target teams distribute parallel do collapse(3) private(k, l) &
    !$omp OMP_DEVICE_ADDR(coeffs, g, label, a, b, c)
    do b_j = 1, (n_j - 1)/SZ + 1
      do b_i = 1, n_i
        do i = 1, SZ
          if (i + (b_j - 1)*SZ <= n_j) then
            k = b_j + (b_i - 1)*n_groups_pad
            do j = 1, n
              l = nint(label(i, j, k))
              if (l < 1 .or. l > n_label) cycle
              a(i, j, k) = a(i, j, k) + coeffs(1, l)*g(i, j, k)
              b(i, j, k) = b(i, j, k) + coeffs(2, l)*g(i, j, k)
              c(i, j, k) = c(i, j, k) + coeffs(3, l)*g(i, j, k)
            end do
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine labelled_axpy

  subroutine field_max_sum(max_f, sum_f, f, n_f, n, n_i, n_groups_pad, &
                           n_j)
    !! Maximum and sum of |f|.
    real(dp), intent(out) :: max_f, sum_f
    integer, intent(in) :: n_f, n, n_i, n_groups_pad, n_j
    real(dp), intent(in) :: f(SZ, n_f, n_groups_pad*n_i)

    integer :: i, j, b_i, b_j

    max_f = 0._dp
    sum_f = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(max:max_f) reduction(+:sum_f) map(tofrom:max_f, sum_f) &
    !$omp OMP_DEVICE_ADDR(f)
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
    real(dp), intent(in) :: f(SZ, n_f, n_groups_pad*n_i)

    integer :: i, j, b_i, b_j

    s = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(+:s) map(tofrom:s) OMP_DEVICE_ADDR(f)
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
    real(dp), intent(in) :: f(SZ, n_f, n_groups_pad*n_i)

    integer :: i, b_i, b_j

    max_f = -huge(1._dp)
    sum_f = 0._dp
    !$omp target teams distribute parallel do collapse(3) &
    !$omp reduction(max:max_f) reduction(+:sum_f) map(tofrom:max_f, sum_f) &
    !$omp OMP_DEVICE_ADDR(f)
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
    real(dp), intent(inout) :: f(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: c_start, c_end, flow_rate_diff

    integer :: i, y_block, z, b

    !$omp target teams distribute parallel do collapse(3) private(b) &
    !$omp OMP_DEVICE_ADDR(f)
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
    real(dp), intent(inout) :: f(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: f_start(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: c_end, flow_rate_diff

    integer :: i, y_block, z, b

    !$omp target teams distribute parallel do collapse(3) private(b) &
    !$omp OMP_DEVICE_ADDR(f, f_start)
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
    real(dp), intent(inout) :: f(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: g(SZ, n_f, n_yb_pad*nz)
    logical, intent(in) :: set_start, set_end

    integer :: i, y_block, z, b

    !$omp target teams distribute parallel do collapse(3) private(b) &
    !$omp OMP_DEVICE_ADDR(f, g)
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
    real(dp), intent(inout) :: f(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: c_start, c_end
    logical, intent(in) :: set_start, set_end

    integer :: j, z, i_top, yb_top

    i_top = mod(ny - 1, SZ) + 1
    yb_top = (ny - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) OMP_DEVICE_ADDR(f)
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
    real(dp), intent(inout) :: f(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: f_start(SZ, n_f, n_yb_pad*nz)

    integer :: j, z, i_top, yb_top

    i_top = mod(ny - 1, SZ) + 1
    yb_top = (ny - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(f, f_start)
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
    real(dp), intent(inout) :: f(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: g(SZ, n_f, n_yb_pad*nz)
    logical, intent(in) :: set_start, set_end

    integer :: j, z, i_top, yb_top

    i_top = mod(ny - 1, SZ) + 1
    yb_top = (ny - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) OMP_DEVICE_ADDR(f, g)
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
    real(dp), intent(inout) :: f(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: c

    integer :: j, z, i_p, yb_p

    i_p = mod(plane - 1, SZ) + 1
    yb_p = (plane - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) OMP_DEVICE_ADDR(f)
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
    real(dp), intent(inout) :: stress(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: u(SZ, n_f, n_yb_pad*nz), w(SZ, n_f, n_yb_pad*nz)
    real(dp), intent(in) :: drag_coeff

    integer :: j, z, i_sa, yb_sa, i_st, yb_st

    i_sa = mod(sample_plane - 1, SZ) + 1
    yb_sa = (sample_plane - 1)/SZ + 1
    i_st = mod(stress_plane - 1, SZ) + 1
    yb_st = (stress_plane - 1)/SZ + 1

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(stress, u, w)
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
    real(dp), intent(inout) :: partial(ny, nz)
    real(dp), intent(in) :: f(SZ, n_f, n_yb_pad*nz)

    integer :: j, z, k

    !$omp target teams distribute parallel do collapse(2) &
    !$omp OMP_DEVICE_ADDR(partial, f)
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
    real(dp), intent(inout) :: sums(ny)
    real(dp), intent(in) :: partial(ny, nz)

    integer :: j, z

    !$omp target teams distribute parallel do OMP_DEVICE_ADDR(sums, partial)
    do j = 1, ny
      sums(j) = 0._dp
      do z = 1, nz
        sums(j) = sums(j) + partial(j, z)
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine field_plane_sums_reduce

  ! --------------------------------------------------------------------------
  ! Reordering, in the layouts of m_ordering. With nx, ny, nz the padded
  ! Cartesian extents and lane = mod(. - 1, SZ) + 1, block = (. - 1)/SZ + 1:
  !   DIR_X: u(lane(y), x, block(y) + (z - 1)*ny/SZ)
  !   DIR_Y: u(lane(x), y, block(x) + (z - 1)*nx/SZ)
  !   DIR_Z: u(lane(x), z, block(x) + (y - 1)*nx/SZ)
  !   DIR_C: u(x, y, z)
  ! Each kernel walks its output in storage order, so the writes are
  ! coalesced, and gathers the input through `offset`. The gather is strided
  ! only where the fast index changes (to or from DIR_X), and those strides
  ! stay within an SZ x SZ tile that L2 serves; staging that tile in team
  ! shared memory measured slower under NVHPC than this direct form.
  ! --------------------------------------------------------------------------

  pure integer function offset(dir, x, y, z, nx, ny, nz) result(l)
    !$omp declare target
    !! 1-based linear offset of Cartesian point (x, y, z) in layout `dir`.

    integer, intent(in) :: dir, x, y, z, nx, ny, nz

    select case (dir)
    case (DIR_X)
      l = mod(y - 1, SZ) + 1 + SZ*(x - 1) &
          + SZ*nx*((y - 1)/SZ + (z - 1)*(ny/SZ))
    case (DIR_Y)
      l = mod(x - 1, SZ) + 1 + SZ*(y - 1) &
          + SZ*ny*((x - 1)/SZ + (z - 1)*(nx/SZ))
    case (DIR_Z)
      l = mod(x - 1, SZ) + 1 + SZ*(z - 1) &
          + SZ*nz*((x - 1)/SZ + (y - 1)*(nx/SZ))
    case default
      l = x + nx*(y - 1) + nx*ny*(z - 1)
    end select
  end function offset

  subroutine reorder(u_o, u_i, dir_from, dir_to, cart_padded)
    !! Reorders u_i in layout dir_from into u_o in layout dir_to.

    real(dp), intent(inout) :: u_o(*)
    real(dp), intent(in) :: u_i(*)
    integer, intent(in) :: dir_from, dir_to
    integer, intent(in) :: cart_padded(3)

    integer :: nx, ny, nz

    nx = cart_padded(1); ny = cart_padded(2); nz = cart_padded(3)

    select case (dir_to)
    case (DIR_X)
      call gather_x(u_o, u_i, dir_from, .false., nx, ny, nz)
    case (DIR_Y)
      call gather_y(u_o, u_i, dir_from, nx, ny, nz)
    case (DIR_Z)
      call gather_z(u_o, u_i, dir_from, nx, ny, nz)
    case default
      call gather_c(u_o, u_i, dir_from, nx, ny, nz)
    end select

  end subroutine reorder

  subroutine sum_intox(u, u_d, dir_from, cart_padded)
    !! u += u_d, where u is DIR_X and u_d the same field in layout dir_from.
    real(dp), intent(inout) :: u(*)
    real(dp), intent(in) :: u_d(*)
    integer, intent(in) :: dir_from
    integer, intent(in) :: cart_padded(3)

    call gather_x(u, u_d, dir_from, .true., cart_padded(1), &
                  cart_padded(2), cart_padded(3))

  end subroutine sum_intox

  subroutine gather_x(u_o, u_i, dir_from, accumulate, nx, ny, nz)
    !! Writes (or with `accumulate` adds) u_i into the DIR_X field u_o.
    integer, intent(in) :: dir_from, nx, ny, nz
    logical, intent(in) :: accumulate
    real(dp), intent(inout) :: u_o(SZ, nx, ny/SZ*nz)
    real(dp), intent(in) :: u_i(nx*ny*nz)

    integer :: lane, x, yb, z

    !$omp target teams distribute parallel do collapse(4) &
    !$omp OMP_DEVICE_ADDR(u_o, u_i)
    do z = 1, nz
      do yb = 1, ny/SZ
        do x = 1, nx
          do lane = 1, SZ
            if (accumulate) then
              u_o(lane, x, yb + (z - 1)*(ny/SZ)) = &
                u_o(lane, x, yb + (z - 1)*(ny/SZ)) &
                + u_i(offset(dir_from, x, (yb - 1)*SZ + lane, z, nx, ny, nz))
            else
              u_o(lane, x, yb + (z - 1)*(ny/SZ)) = &
                u_i(offset(dir_from, x, (yb - 1)*SZ + lane, z, nx, ny, nz))
            end if
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine gather_x

  subroutine gather_y(u_o, u_i, dir_from, nx, ny, nz)
    integer, intent(in) :: dir_from, nx, ny, nz
    real(dp), intent(inout) :: u_o(SZ, ny, nx/SZ*nz)
    real(dp), intent(in) :: u_i(nx*ny*nz)

    integer :: lane, xb, y, z

    !$omp target teams distribute parallel do collapse(4) &
    !$omp OMP_DEVICE_ADDR(u_o, u_i)
    do z = 1, nz
      do xb = 1, nx/SZ
        do y = 1, ny
          do lane = 1, SZ
            u_o(lane, y, xb + (z - 1)*(nx/SZ)) = &
              u_i(offset(dir_from, (xb - 1)*SZ + lane, y, z, nx, ny, nz))
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine gather_y

  subroutine gather_z(u_o, u_i, dir_from, nx, ny, nz)
    !! The y lane runs just outside the x lane, rather than y outermost, so
    !! the SZ x SZ tile a DIR_X input is gathered from is read by
    !! neighbouring threads while it is still in cache.
    integer, intent(in) :: dir_from, nx, ny, nz
    real(dp), intent(inout) :: u_o(SZ, nz, nx/SZ*ny)
    real(dp), intent(in) :: u_i(nx*ny*nz)

    integer :: lane, xb, y_lane, yb, z

    !$omp target teams distribute parallel do collapse(5) &
    !$omp OMP_DEVICE_ADDR(u_o, u_i)
    do yb = 1, ny/SZ
      do z = 1, nz
        do xb = 1, nx/SZ
          do y_lane = 1, SZ
            do lane = 1, SZ
              u_o(lane, z, xb + ((yb - 1)*SZ + y_lane - 1)*(nx/SZ)) = &
                u_i(offset(dir_from, (xb - 1)*SZ + lane, &
                           (yb - 1)*SZ + y_lane, z, nx, ny, nz))
            end do
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine gather_z

  subroutine gather_c(u_o, u_i, dir_from, nx, ny, nz)
    integer, intent(in) :: dir_from, nx, ny, nz
    real(dp), intent(inout) :: u_o(nx, ny, nz)
    real(dp), intent(in) :: u_i(nx*ny*nz)

    integer :: x, y, z

    !$omp target teams distribute parallel do collapse(3) &
    !$omp OMP_DEVICE_ADDR(u_o, u_i)
    do z = 1, nz
      do y = 1, ny
        do x = 1, nx
          u_o(x, y, z) = u_i(offset(dir_from, x, y, z, nx, ny, nz))
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine gather_c

  subroutine sample_points(vals, points, a, b, c, n_points, n_f, n_k, d, &
                           grid)
    !! vals(1:3, p) = this rank's part of the trilinear interpolation of the
    !! DIR_X VERT fields a, b and c at points(:, p), one thread per point;
    !! see field_sample_points. grid holds the local VERT dims, global VERT
    !! dims, rank offsets and periodic flags (1/0), three entries each, then
    !! the padded y-block count of the DIR_X layout.
    integer, intent(in) :: n_points, n_f, n_k
    real(dp), intent(out) :: vals(3, n_points)
    real(dp), intent(in) :: points(3, n_points), d(3)
    integer, intent(in) :: grid(13)
    real(dp), intent(in) :: a(SZ, n_f, n_k), b(SZ, n_f, n_k), c(SZ, n_f, n_k)

    real(dp) :: fx, fy, fz, w, sa, sb, sc
    integer :: p, o1, o2, o3, x0, y0, z0, lx, ly, lz, i, k

    !$omp target teams distribute parallel do &
    !$omp private(fx, fy, fz, w, sa, sb, sc, o1, o2, o3, x0, y0, z0, &
    !$omp         lx, ly, lz, i, k) &
    !$omp map(to: points, d, grid) map(from: vals) OMP_DEVICE_ADDR(a, b, c)
    do p = 1, n_points
      fx = points(1, p)/d(1)
      fy = points(2, p)/d(2)
      fz = points(3, p)/d(3)
      x0 = floor(fx)
      y0 = floor(fy)
      z0 = floor(fz)
      fx = fx - x0
      fy = fy - y0
      fz = fz - z0
      sa = 0._dp
      sb = 0._dp
      sc = 0._dp
      do o3 = 0, 1
        lz = z0 + o3
        if (grid(12) /= 0) lz = modulo(lz, grid(6))
        if (lz < 0 .or. lz >= grid(6)) cycle
        lz = lz - grid(9) + 1
        if (lz < 1 .or. lz > grid(3)) cycle
        do o2 = 0, 1
          ly = y0 + o2
          if (grid(11) /= 0) ly = modulo(ly, grid(5))
          if (ly < 0 .or. ly >= grid(5)) cycle
          ly = ly - grid(8) + 1
          if (ly < 1 .or. ly > grid(2)) cycle
          do o1 = 0, 1
            lx = x0 + o1
            if (grid(10) /= 0) lx = modulo(lx, grid(4))
            if (lx < 0 .or. lx >= grid(4)) cycle
            lx = lx - grid(7) + 1
            if (lx < 1 .or. lx > grid(1)) cycle
            w = merge(fx, 1._dp - fx, o1 == 1) &
                *merge(fy, 1._dp - fy, o2 == 1) &
                *merge(fz, 1._dp - fz, o3 == 1)
            i = mod(ly - 1, SZ) + 1
            k = (ly - 1)/SZ + 1 + (lz - 1)*grid(13)
            sa = sa + w*a(i, lx, k)
            sb = sb + w*b(i, lx, k)
            sc = sc + w*c(i, lx, k)
          end do
        end do
      end do
      vals(1, p) = sa
      vals(2, p) = sb
      vals(3, p) = sc
    end do
    !$omp end target teams distribute parallel do

  end subroutine sample_points

  subroutine add_gaussian_sources(a, b, c, src, box, n_src, n_box_y, &
                                  n_box_z, inv_eps2, r_cut, n_f, n_k, d, &
                                  grid)
    !! a, b, c += src(4:6, s)*exp(-dist^2/eps^2) within r_cut of the point
    !! src(1:3, s), over the global index box box(1:3, s)..box(4:6, s); see
    !! source_index_boxes. Box indices are wrapped across periodic
    !! boundaries and skipped where this rank does not own them. One thread
    !! per (source, z, y) of the box, each sweeping x; sources overlap, so
    !! the adds are atomic. grid holds the local VERT dims, global VERT dims,
    !! rank offsets and periodic flags (1/0), three entries each, then the
    !! padded y-block count.
    integer, intent(in) :: n_src, n_box_y, n_box_z, n_f, n_k
    real(dp), intent(in) :: src(6, n_src), inv_eps2, r_cut, d(3)
    integer, intent(in) :: box(6, n_src), grid(13)
    real(dp), intent(inout) :: a(SZ, n_f, n_k), b(SZ, n_f, n_k), &
                               c(SZ, n_f, n_k)

    real(dp) :: dz2, dyz2, dist2, kern
    integer :: s, ky, kz, g_x, g_y, g_z, x, y, z, i, k

    !$omp target teams distribute parallel do collapse(3) &
    !$omp private(dz2, dyz2, dist2, kern, g_x, g_y, g_z, x, y, z, i, k) &
    !$omp map(to: src, box, d, grid) OMP_DEVICE_ADDR(a, b, c)
    do s = 1, n_src
      do kz = 0, n_box_z - 1
        do ky = 0, n_box_y - 1
          g_z = box(3, s) + kz
          g_y = box(2, s) + ky
          if (g_z > box(6, s) .or. g_y > box(5, s)) cycle
          ! Wrap periodic indices; skip the ones this rank does not own
          z = g_z
          if (grid(12) /= 0) z = modulo(g_z, grid(6))
          z = z - grid(9) + 1
          if (z < 1 .or. z > grid(3)) cycle
          y = g_y
          if (grid(11) /= 0) y = modulo(g_y, grid(5))
          y = y - grid(8) + 1
          if (y < 1 .or. y > grid(2)) cycle
          dz2 = (src(3, s) - g_z*d(3))**2
          dyz2 = dz2 + (src(2, s) - g_y*d(2))**2
          if (dyz2 > r_cut**2) cycle
          i = mod(y - 1, SZ) + 1
          k = (y - 1)/SZ + 1 + (z - 1)*grid(13)
          do g_x = box(1, s), box(4, s)
            x = g_x
            if (grid(10) /= 0) x = modulo(g_x, grid(4))
            x = x - grid(7) + 1
            if (x < 1 .or. x > grid(1)) cycle
            dist2 = dyz2 + (src(1, s) - g_x*d(1))**2
            if (dist2 > r_cut**2) cycle
            kern = exp(-dist2*inv_eps2)
            !$omp atomic update
            a(i, x, k) = a(i, x, k) + src(4, s)*kern
            !$omp atomic update
            b(i, x, k) = b(i, x, k) + src(5, s)*kern
            !$omp atomic update
            c(i, x, k) = c(i, x, k) + src(6, s)*kern
          end do
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

  end subroutine add_gaussian_sources

end module m_omptgt_kernels_fieldops
