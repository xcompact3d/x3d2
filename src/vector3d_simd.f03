module m_vector3d_simd
  use m_vector3d, only: vector3d
  type, extends(vector3d) :: vector3d_simd
     real, allocatable :: data(:,:,:, 3)
   contains
     procedure, public :: transport
     procedure, public :: div
  end type vector3d_simd

  interface vector3d_simd
     module procedure construct
  end interface vector3d_simd

  real, allocatable :: du(:, :, :)
  real, allocatable :: d2u(:, :, :)
  real, allocatable :: u2(:, :, :)
  real, allocatable :: du2(:, :, :)

contains

  subroutine ensure_work_arrays(dims)
    integer, intent(in) :: dims(3)

    if (.not. allocated(du)) then
       allocate(du(dims(1), dims(2), dims(3)))
       allocate(d2u(dims(1), dims(2), dims(3)))
       allocate(u2(dims(1), dims(2), dims(3)))
       allocate(du2(dims(1), dims(2), dims(3)))
    else
       if ( &
            & dims(1) > dim(du, 1) .or. &
            & dims(2) > dim(du, 2) .or. &
            & dims(3) > dim(du, 3) &
            & ) then
          du = reshape(du, dims)
          d2u = reshape(d2u, dims)
          u2 = reshape(u2, dims)
          du2 = reshape(du2, dims)
       end if
    end if
  end subroutine ensure_work_arrays

  type(vector3d_simd) function construct(name, dims)
    character(*), intent(in) :: name
    integer, intent(in) :: dims(3)

    call ensure_work_arrays(dims)

    call mpi_comm_rank(MPI_COMM_WORLD, construct%rankid, errcode)
    call mpi_comm_size(MPI_COMM_WORLD, construct%nranks, errcode)

    allocate(construct%u(dims(1), dims(2), dims(3)))
    allocate(construct%v(dims(1), dims(2), dims(3)))
    allocate(construct%w(dims(1), dims(2), dims(3)))

    construct%name = name
  end function construct

  subroutine transport(self, rslt)
    class(vector3d_simd), intent(in) :: self
    class(vector3d), intent(inout) :: rslt
    real, pointer :: u(:, :, :), u_dir(:, :, :)

    u_dir => self%data(:, :, :, 1)
    select type (rslt)
    type is (vector3d_simd)
       components: do l = 1, size(self%data, 4)
          u => self%data(:, :, :, l)
          rslt%data(:, :, :, l) = transport_dir(u, u_dir)
       end do components
    class default
       error stop
    end select
  end subroutine transport

  pure function transport_dir(u, u_dir)
    real, intent(in) :: u(:, :, :)
    real, intent(in) :: u_dir(:, :, :)
    real :: transport_dir(size(u))

    layers: do k = 1, size(u, 3)
       call diff(u(:, :, k), du)
       call diff2(u(:, :, k), d2u)
       do j = 1, size(u, 2)
          !$omp simd
          do i = 1, size(u, 1)
             u2(i, j) = u(i, j, k) * u_dir(i, j, k)
          end do
          !$omp end simd
       end do
       call diff(du2, u2)
       reshape(transport_dir, shape(u))
       do j = 1, size(u, 2)
          !$omp simd
          do i = 1, size(u, 1)
             rslt%u(i, j, k) = -0.5 * &
                  (u(i, j, k) * du(i, j) + du2(i, j)) &
                  & + xnu * d2u(i, j)
          !$omp end simd
          end do
       end do
    end do layers
  end function transport_dir

  subroutine diff(self, f, df)
    class(vector3d_simd), intent(in) :: self
    real, intent(in) :: f(:, :)
    real, intent(out) :: df(:, :)
    integer, parameter :: n = size(f, 2)
    type(stencil), pointer :: s
    !$omp simd
      do i = 1, size(f, 1)
         s => self%stencils(1)
         du(i, 1) = dot_product(s%coeffs, u(i, s%nodes + 1))
         s => self%stencils(2)
         du(i, 2) = dot_product(s%coeffs, u(i, s%nodes + 2))
      end do
      !$omp end simd
      s => self%stencils(3)
      do j = 3, n - 2
         !$omp simd
         do i = 1, size(f, 1)
            du(i, j) = dot_product(s%coeffs, u(i, s%nodes + j))
         end do
         !$omp end simd
      end do
      !$omp simd
      do i = 1, size(f, 1)
         s => self%stencils(4)
         du(i, n - 1) = dot_product(s%coeffs, u(i, s%nodes + n - 1))
         s => self%stencils(5)
         du(i, n) = dot_product(s%coeffs, u(i, s%nodes + n))
      end do
      !$omp end simd

      call self%thomas_solver%solve(du, df)
  end subroutine diff

  subroutine div(self, rslt)
    class(vector3d_simd), intent(in) :: self
    class(vector3d), intent(inout) :: rslt

    select type (rslt)
    type is (vector3d_simd)
       rslt%u = self%u + 1.
       rslt%v = self%v + 1.
       rslt%w = self%w + 1.
    class default
       error stop
    end select
  end subroutine div
end module m_vector3d_simd
