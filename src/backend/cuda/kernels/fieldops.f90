module m_cuda_kernels_fieldops
  use cudafor

  use m_common, only: dp
  use m_cuda_common, only: SZ

contains

  attributes(global) subroutine copy(n, dst, src)
    implicit none

    integer, value, intent(in) :: n
    real(dp), device, intent(out), dimension(:, :, :) :: dst
    real(dp), device, intent(in), dimension(:, :, :) :: src

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      dst(i, j, b) = src(i, j, b)
    end do

  end subroutine copy

  attributes(global) subroutine axpby(n, alpha, x, beta, y)
    implicit none

    integer, value, intent(in) :: n
    real(dp), value, intent(in) :: alpha, beta
    real(dp), device, intent(in), dimension(:, :, :) :: x
    real(dp), device, intent(inout), dimension(:, :, :) :: y

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      y(i, j, b) = alpha*x(i, j, b) + beta*y(i, j, b)
    end do

  end subroutine axpby

  attributes(global) subroutine pwmul(y, x, n)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: y
    real(dp), device, intent(in), dimension(:, :, :) :: x
    integer, value, intent(in) :: n

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      y(i, j, b) = y(i, j, b)*x(i, j, b)
    end do

  end subroutine pwmul

  attributes(global) subroutine buffer_copy(u_send_s, u_send_e, u, n, n_halo)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: u_send_s, u_send_e
    real(dp), device, intent(in), dimension(:, :, :) :: u
    integer, value, intent(in) :: n, n_halo

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n_halo
      u_send_s(i, j, b) = u(i, j, b)
      u_send_e(i, j, b) = u(i, n - n_halo + j, b)
    end do

  end subroutine buffer_copy

  attributes(global) subroutine field_scale(f, alpha, n)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: f
    real(dp), value, intent(in) :: alpha
    integer, value, intent(in) :: n

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      f(i, j, b) = alpha*f(i, j, b)
    end do

  end subroutine field_scale

  attributes(global) subroutine field_shift(f, const, n)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: f
    real(dp), value, intent(in) :: const
    integer, value, intent(in) :: n

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      f(i, j, b) = f(i, j, b) + const
    end do

  end subroutine field_shift

  attributes(global) subroutine vorticity_from_gradients( &
    vort, dudy, dudz, dvdx, dvdz, dwdx, dwdy, n)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: vort
    real(dp), device, intent(in), dimension(:, :, :) :: &
      dudy, dudz, dvdx, dvdz, dwdx, dwdy
    integer, value, intent(in) :: n

    integer :: i, j, b
    real(dp) :: wx, wy, wz

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      wx = dwdy(i, j, b) - dvdz(i, j, b)
      wy = dudz(i, j, b) - dwdx(i, j, b)
      wz = dvdx(i, j, b) - dudy(i, j, b)
      vort(i, j, b) = sqrt(wx*wx + wy*wy + wz*wz)
    end do

  end subroutine vorticity_from_gradients

  attributes(global) subroutine qcriterion_from_gradients( &
    qcrit, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, n)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: qcrit
    real(dp), device, intent(in), dimension(:, :, :) :: dudx, dudy, dudz
    real(dp), device, intent(in), dimension(:, :, :) :: dvdx, dvdy, dvdz
    real(dp), device, intent(in), dimension(:, :, :) :: dwdx, dwdy, dwdz
    integer, value, intent(in) :: n

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      qcrit(i, j, b) = -0.5_dp*(dudx(i, j, b)*dudx(i, j, b) + &
                                dvdy(i, j, b)*dvdy(i, j, b) + &
                                dwdz(i, j, b)*dwdz(i, j, b)) - &
                       dudy(i, j, b)*dvdx(i, j, b) - &
                       dudz(i, j, b)*dwdx(i, j, b) - &
                       dvdz(i, j, b)*dwdy(i, j, b)
    end do

  end subroutine qcriterion_from_gradients

  attributes(global) subroutine smagorinsky_from_gradients( &
    nut, mixing_length_sq, dudx, dudy, dudz, dvdx, dvdy, dvdz, &
    dwdx, dwdy, dwdz, n)
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: nut
    real(dp), device, intent(in), dimension(:, :, :) :: mixing_length_sq
    real(dp), device, intent(in), dimension(:, :, :) :: dudx, dudy, dudz
    real(dp), device, intent(in), dimension(:, :, :) :: dvdx, dvdy, dvdz
    real(dp), device, intent(in), dimension(:, :, :) :: dwdx, dwdy, dwdz
    integer, value, intent(in) :: n

    real(dp) :: sij_sq
    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      if (mixing_length_sq(i, j, b) > 0._dp) then
        sij_sq = dudx(i, j, b)**2 + dvdy(i, j, b)**2 + &
                 dwdz(i, j, b)**2 + &
                 0.5_dp*(dudy(i, j, b) + dvdx(i, j, b))**2 + &
                 0.5_dp*(dudz(i, j, b) + dwdx(i, j, b))**2 + &
                 0.5_dp*(dvdz(i, j, b) + dwdy(i, j, b))**2
        nut(i, j, b) = mixing_length_sq(i, j, b)*sqrt(2._dp*sij_sq)
      else
        nut(i, j, b) = 0._dp
      end if
    end do

  end subroutine smagorinsky_from_gradients

  attributes(global) subroutine sgs_stress_from_gradients( &
    stress, nut, gradient_a, gradient_b, scale_a, scale_b, n)
    implicit none

    real(dp), device, intent(out), dimension(:, :, :) :: stress
    real(dp), device, intent(in), dimension(:, :, :) :: &
      nut, gradient_a, gradient_b
    real(dp), value, intent(in) :: scale_a, scale_b
    integer, value, intent(in) :: n

    integer :: i, j, b

    i = threadIdx%x
    b = blockIdx%x

    do j = 1, n
      if (nut(i, j, b) > 0._dp) then
        stress(i, j, b) = nut(i, j, b)* &
                          (scale_a*gradient_a(i, j, b) &
                           + scale_b*gradient_b(i, j, b))
      else
        stress(i, j, b) = 0._dp
      end if
    end do
  end subroutine sgs_stress_from_gradients

  attributes(global) subroutine abl_wall_boundary_correction( &
    sgs_u, sgs_v, sgs_w, nut, dudy, dvdx, dwdy, dvdz, &
    wall_u_sample, wall_w_sample, drag_coeff, sampling_height, nx, ny)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: sgs_u, sgs_v, sgs_w
    real(dp), device, intent(in), dimension(:, :, :) :: nut
    real(dp), device, intent(in), dimension(:, :, :) :: dudy, dvdx, dwdy, dvdz
    real(dp), value, intent(in) :: wall_u_sample, wall_w_sample
    real(dp), value, intent(in) :: drag_coeff, sampling_height
    integer, value, intent(in) :: nx, ny

    integer :: i, z, b, n_y_blocks
    real(dp) :: speed, tau_x, tau_z
    real(dp) :: sxy_2, sxy_3, syz_2, syz_3

    i = threadIdx%x + (blockIdx%x - 1)*blockDim%x
    z = blockIdx%y
    n_y_blocks = (ny - 1)/SZ + 1
    b = 1 + (z - 1)*n_y_blocks

    if (i <= nx) then
      speed = sqrt(wall_u_sample**2 + wall_w_sample**2)
      tau_x = -drag_coeff*wall_u_sample*speed
      tau_z = -drag_coeff*wall_w_sample*speed

      sxy_2 = 0.5_dp*(dudy(2, i, b) + dvdx(2, i, b))
      sxy_3 = 0.5_dp*(dudy(3, i, b) + dvdx(3, i, b))
      syz_2 = 0.5_dp*(dwdy(2, i, b) + dvdz(2, i, b))
      syz_3 = 0.5_dp*(dwdy(3, i, b) + dvdz(3, i, b))

      sgs_u(1, i, b) = -( &
                       -0.5_dp*(-2._dp*nut(3, i, b)*sxy_3) + &
                       2._dp*(-2._dp*nut(2, i, b)*sxy_2) - &
                       1.5_dp*tau_x)/(2._dp*sampling_height)
      sgs_v(1, i, b) = 0._dp
      sgs_w(1, i, b) = -( &
                       -0.5_dp*(-2._dp*nut(3, i, b)*syz_3) + &
                       2._dp*(-2._dp*nut(2, i, b)*syz_2) - &
                       1.5_dp*tau_z)/(2._dp*sampling_height)
    end if
  end subroutine abl_wall_boundary_correction

  attributes(global) subroutine scalar_product(s, x, y, n, n_i_pad, n_j)
    implicit none

    real(dp), device, intent(inout) :: s
    real(dp), device, intent(in), dimension(:, :, :) :: x, y
    integer, value, intent(in) :: n, n_i_pad, n_j

    real(dp) :: s_pncl !! pencil sum
    integer :: i, j, b, b_i, b_j, ierr

    i = threadIdx%x
    b_i = blockIdx%x
    b_j = blockIdx%y

    b = b_i + (b_j - 1)*n_i_pad
    s_pncl = 0._dp
    if (i + (b_j - 1)*blockDim%x <= n_j) then
      do j = 1, n
        s_pncl = s_pncl + x(i, j, b)*y(i, j, b)
      end do
    end if
    ierr = atomicadd(s, s_pncl)

  end subroutine scalar_product

  attributes(global) subroutine vector_norm_squared(s, a, b, c, &
                                                    n, n_i_pad, n_j)
    implicit none

    real(dp), device, intent(inout) :: s
    real(dp), device, intent(in), dimension(:, :, :) :: a, b, c
    integer, value, intent(in) :: n, n_i_pad, n_j

    real(dp) :: pencil_sum
    integer :: i, j, pencil, block_i, block_j, ierr

    i = threadIdx%x
    block_i = blockIdx%x
    block_j = blockIdx%y
    pencil = block_i + (block_j - 1)*n_i_pad

    pencil_sum = 0._dp
    if (i + (block_j - 1)*blockDim%x <= n_j) then
      do j = 1, n
        pencil_sum = pencil_sum + a(i, j, pencil)**2 + &
                     b(i, j, pencil)**2 + c(i, j, pencil)**2
      end do
    end if
    ierr = atomicadd(s, pencil_sum)
  end subroutine vector_norm_squared

  attributes(global) subroutine field_max_sum(max_f, sum_f, f, n, n_i_pad, n_j)
    implicit none

    real(dp), device, intent(inout) :: max_f, sum_f
    real(dp), device, intent(in), dimension(:, :, :) :: f
    integer, value, intent(in) :: n, n_i_pad, n_j

    real(dp) :: max_pncl, sum_pncl, val
    integer :: i, j, b, b_i, b_j, ierr

    i = threadIdx%x
    b_i = blockIdx%x
    b_j = blockIdx%y

    b = b_i + (b_j - 1)*n_i_pad
    max_pncl = 0._dp
    sum_pncl = 0._dp
    if (i + (b_j - 1)*blockDim%x <= n_j) then
      do j = 1, n
        val = abs(f(i, j, b))
        sum_pncl = sum_pncl + val
        max_pncl = max(max_pncl, val)
      end do
    end if
    ierr = atomicadd(sum_f, sum_pncl)
    ierr = atomicmax(max_f, max_pncl)

  end subroutine field_max_sum

  attributes(global) subroutine field_set_y_face( &
    f, c_start, c_end, flow_rate_diff, nx, ny, nz)
  !! Set domain Y_FACE boundary values.
  !! c_start: Dirichlet value applied at the bottom face (j = 1)
  !! c_end:   Dirichlet value applied at the top face (j = ny)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: f
    real(dp), value, intent(in) :: c_start, c_end, flow_rate_diff
    integer, value, intent(in) :: nx, ny, nz

    integer :: i, j, b, n_mod, b_end

    j = threadIdx%x + (blockIdx%x - 1)*blockDim%x ! from 1 to nx
    b = blockIdx%y ! from 1 to nz

    n_mod = mod(ny - 1, SZ) + 1
    b_end = b + (ny - 1)/SZ*nz

    if (j <= nx) then
      f(1, j, b) = c_start
      f(n_mod, j, b_end) = c_end
    end if

  end subroutine field_set_y_face
  attributes(global) subroutine field_set_y_face_from_field( &
    f, f_start, nx, ny, nz)
!! Set domain Y_FACE boundary values from another field.
!! Both f and f_start are DIR_X VERT pencil-layout. Only the bottom
!! (j = 1) and top (j = ny) y-pencil planes of f_start are read.
!!
!! Launch convention matches field_set_y_face exactly:
!!   threads = dim3(64, 1, 1)
!!   blocks  = dim3((nx-1)/64 + 1, nz, 1)
    implicit none

    real(dp), device, intent(inout), dimension(:, :, :) :: f
    real(dp), device, intent(in), dimension(:, :, :) :: f_start
    integer, value, intent(in) :: nx, ny, nz

    integer :: j, b, n_mod, b_end

    j = threadIdx%x + (blockIdx%x - 1)*blockDim%x ! from 1 to nx
    b = blockIdx%y ! from 1 to nz

    n_mod = mod(ny - 1, SZ) + 1
    b_end = b + (ny - 1)/SZ*nz

    if (j <= nx) then
      f(1, j, b) = f_start(1, j, b)
      f(n_mod, j, b_end) = f_start(n_mod, j, b_end)
    end if

  end subroutine field_set_y_face_from_field
  attributes(global) subroutine field_set_x_face( &
    f, c_start, c_end, bc_start, bc_end, flow_rate_diff, nx, ny, nz)
  !! Set domain X_FACE boundary values.
  !! c_start: Dirichlet value applied at the left face (i = 1)
  !! c_end:   convective velocity Uc = uxmax * gdt / dx,
  !!          used as multiplier in the outflow scheme du/dt + Uc*du/dx = 0
    implicit none
    real(dp), device, intent(inout), dimension(:, :, :) :: f
    real(dp), value, intent(in) :: c_start, c_end, flow_rate_diff
    integer, value, intent(in) :: bc_start, bc_end
    integer, value, intent(in) :: nx, ny, nz
    integer :: i, b, n_mod, n_y_blocks, y_block, i_max

    i = threadIdx%x + (blockIdx%x - 1)*blockDim%x  ! 1..SZ
    b = blockIdx%y                                 ! 1..n_y_blocks*nz

    n_mod = mod(ny - 1, SZ) + 1
    n_y_blocks = (ny - 1)/SZ + 1

    ! Determine if this b is the last y-block (padding present)
    ! With b = (y_block - 1)*nz + z, y_block = (b - 1)/nz + 1
    y_block = (b - 1)/nz + 1
    if (y_block == n_y_blocks) then
      i_max = n_mod    ! last y-block: only i in [1, n_mod] are real
    else
      i_max = SZ       ! interior y-blocks: all i in [1, SZ] are real
    end if

    if (i <= i_max) then
      ! --- Left face (j = 1) ---
      select case (bc_start)
      case (1) ! BC_NEUMANN
        !this can be empty for now future TODO
      case (2) ! BC_DIRICHLET
        f(i, 1, b) = c_start
      end select

      ! --- Right face (j = nx) ---
      select case (bc_end)
      case (1) !BC_NEUMANN
        !this can be empty for now future TODO
      case (2) ! BC_DIRICHLET
        f(i, nx, b) = f(i, nx, b) &
                      - c_end*(f(i, nx, b) - f(i, nx - 1, b)) &
                      + flow_rate_diff
      end select
    end if
  end subroutine field_set_x_face

  attributes(global) subroutine field_set_x_face_from_field( &
    f, f_start, c_end, bc_start, bc_end, flow_rate_diff, nx, ny, nz)
  !! Set domain X_FACE boundary values with a spatially-varying inlet.
  !!
  !! Both `f` and `f_start` are DIR_X VERT fields with shape
  !! (SZ, nx, n_y_blocks*nz). The inlet plane is pencil index 1 of each,
  !! so f(i, 1, b) <- f_start(i, 1, b) for every real (i, b). Only the
  !! i = 1 pencil-index plane of f_start is read; the rest is ignored.
  !!
  !! Launch convention (matches field_set_x_face exactly):
  !!   threads = dim3(64, 1, 1)
  !!   blocks  = dim3((SZ-1)/64 + 1, n_y_blocks*nz, 1)
  !! Threads beyond the real pencil width (SZ in interior y-blocks,
  !! n_mod in the last y-block) are masked out by the i_max guard.
  !!
  !! c_end is the convective velocity Uc = uxmax * gdt / dx used in the
  !! outflow scheme du/dt + Uc*du/dx = 0 at the right face.
    implicit none
    real(dp), device, intent(inout), dimension(:, :, :) :: f
    real(dp), device, intent(in), dimension(:, :, :) :: f_start
    real(dp), value, intent(in) :: c_end, flow_rate_diff
    integer, value, intent(in) :: bc_start, bc_end
    integer, value, intent(in) :: nx, ny, nz

    integer :: i, b, n_mod, n_y_blocks, y_block, i_max

    i = threadIdx%x + (blockIdx%x - 1)*blockDim%x  ! 1..SZ (or n_mod)
    b = blockIdx%y                                 ! 1..n_y_blocks*nz

    n_mod = mod(ny - 1, SZ) + 1
    n_y_blocks = (ny - 1)/SZ + 1
    y_block = (b - 1)/nz + 1

    if (y_block == n_y_blocks) then
      i_max = n_mod    ! last y-block: only i in [1, n_mod] are real
    else
      i_max = SZ       ! interior y-blocks: all i in [1, SZ] are real
    end if

    if (i <= i_max) then
      ! --- Left face (pencil index 1, i.e. global x = 1) ---
      select case (bc_start)
      case (1) ! BC_NEUMANN
        ! future TODO
      case (2) ! BC_DIRICHLET
        f(i, 1, b) = f_start(i, 1, b)
      end select

      ! --- Right face (pencil index nx, i.e. global x = nx) ---
      select case (bc_end)
      case (1) ! BC_NEUMANN
        ! future TODO
      case (2) ! BC_DIRICHLET
        f(i, nx, b) = f(i, nx, b) &
                      - c_end*(f(i, nx, b) - f(i, nx - 1, b)) &
                      + flow_rate_diff
      end select
    end if
  end subroutine field_set_x_face_from_field

  attributes(global) subroutine volume_integral(s, f, n, n_i_pad, n_j)
    implicit none

    real(dp), device, intent(inout) :: s
    real(dp), device, intent(in), dimension(:, :, :) :: f
    integer, value, intent(in) :: n, n_i_pad, n_j

    real(dp) :: s_pncl !! pencil sum
    integer :: i, j, b, b_i, b_j, ierr

    i = threadIdx%x
    b_i = blockIdx%x
    b_j = blockIdx%y

    b = b_i + (b_j - 1)*n_i_pad
    s_pncl = 0._dp
    if (i + (b_j - 1)*blockDim%x <= n_j) then
      do j = 1, n
        s_pncl = s_pncl + f(i, j, b)
      end do
    end if
    ierr = atomicadd(s, s_pncl)

  end subroutine volume_integral

end module m_cuda_kernels_fieldops
