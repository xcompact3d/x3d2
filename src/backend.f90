module m_base_backend
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_common, only: dp, DIR_C, get_rdr_from_dirs
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: tdsops_t, dirps_t
  use m_mesh, only: mesh_t

  implicit none

  type, abstract :: base_backend_t
      !! base_backend class defines all the abstract operations that the
      !! solver class requires.
      !!
      !! For example, transport equation in solver class evaluates the
      !! derivatives in x, y, and z directions, and reorders the input
      !! fields as required. Then finally, combines all the directional
      !! derivatives to obtain the divergence of U*.
      !!
      !! All these high level operations solver class executes are
      !! defined here using the abstract interfaces. Every backend
      !! implementation extends the present abstact backend class to
      !! define the specifics of these operations based on the target
      !! architecture.

    real(dp) :: nu
    class(mesh_t), pointer :: mesh
    class(allocator_t), pointer :: allocator
    class(poisson_fft_t), pointer :: poisson_fft
  contains
    procedure(transeq_ders), deferred :: transeq_x
    procedure(transeq_ders), deferred :: transeq_y
    procedure(transeq_ders), deferred :: transeq_z
    procedure(tds_solve), deferred :: tds_solve
    procedure(reorder), deferred :: reorder
    procedure(sum_intox), deferred :: sum_yintox
    procedure(sum_intox), deferred :: sum_zintox
    procedure(vecadd), deferred :: vecadd
    procedure(scalar_product), deferred :: scalar_product
    procedure(copy_data_to_f), deferred :: copy_data_to_f
    procedure(copy_f_to_data), deferred :: copy_f_to_data
    procedure(alloc_tdsops), deferred :: alloc_tdsops
    procedure(init_poisson_fft), deferred :: init_poisson_fft
    procedure :: base_init
    procedure :: get_field_data
    procedure :: set_field_data
  end type base_backend_t

  abstract interface
    subroutine transeq_ders(self, du, dv, dw, u, v, w, dirps)
         !! transeq equation obtains the derivatives direction by
         !! direction, and the exact algorithm used to obtain these
         !! derivatives are decided at runtime. Backend implementations
         !! are responsible from directing calls to transeq_ders into
         !! the correct algorithm.
      import :: base_backend_t
      import :: field_t
      import :: dirps_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
      type(dirps_t), intent(in) :: dirps
    end subroutine transeq_ders
  end interface

  abstract interface
    subroutine tds_solve(self, du, u, dirps, tdsops)
         !! transeq equation obtains the derivatives direction by
         !! direction, and the exact algorithm used to obtain these
         !! derivatives are decided at runtime. Backend implementations
         !! are responsible from directing calls to transeq_ders into
         !! the correct algorithm.
      import :: base_backend_t
      import :: field_t
      import :: dirps_t
      import :: tdsops_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: du
      class(field_t), intent(in) :: u
      type(dirps_t), intent(in) :: dirps
      class(tdsops_t), intent(in) :: tdsops
    end subroutine tds_solve
  end interface

  abstract interface
    subroutine reorder(self, u_, u, direction)
         !! reorder subroutines are straightforward, they rearrange
         !! data into our specialist data structure so that regardless
         !! of the direction tridiagonal systems are solved efficiently
         !! and fast.
      import :: base_backend_t
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: u_
      class(field_t), intent(in) :: u
      integer, intent(in) :: direction
    end subroutine reorder
  end interface

  abstract interface
    subroutine sum_intox(self, u, u_)
         !! sum9into3 subroutine combines all the directional velocity
         !! derivatives into the corresponding x directional fields.
      import :: base_backend_t
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(inout) :: u
      class(field_t), intent(in) :: u_
    end subroutine sum_intox
  end interface

  abstract interface
    subroutine vecadd(self, a, x, b, y)
         !! adds two vectors together: y = a*x + b*y
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      real(dp), intent(in) :: a
      class(field_t), intent(in) :: x
      real(dp), intent(in) :: b
      class(field_t), intent(inout) :: y
    end subroutine vecadd
  end interface

  abstract interface
    real(dp) function scalar_product(self, x, y) result(s)
         !! Calculates the scalar product of two input fields
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t) :: self
      class(field_t), intent(in) :: x, y
    end function scalar_product
  end interface

  abstract interface
    subroutine copy_data_to_f(self, f, data)
         !! Copy the specialist data structure from device or host back
         !! to a regular 3D data array in host memory.
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t), intent(inout) :: self
      class(field_t), intent(inout) :: f
      real(dp), dimension(:, :, :), intent(in) :: data
    end subroutine copy_data_to_f

    subroutine copy_f_to_data(self, data, f)
         !! Copy a regular 3D array in host memory into the specialist
         !! data structure field that lives on device or host
      import :: base_backend_t
      import :: dp
      import :: field_t
      implicit none

      class(base_backend_t), intent(inout) :: self
      real(dp), dimension(:, :, :), intent(out) :: data
      class(field_t), intent(in) :: f
    end subroutine copy_f_to_data
  end interface

  abstract interface
    subroutine alloc_tdsops(self, tdsops, dir, operation, scheme, n_halo, &
                            from_to, bc_start, bc_end, sym, c_nu, nu0_nu)
      import :: base_backend_t
      import :: dp
      import :: tdsops_t
      implicit none

      class(base_backend_t) :: self
      class(tdsops_t), allocatable, intent(inout) :: tdsops
      integer, intent(in) :: dir
      character(*), intent(in) :: operation, scheme
      integer, optional, intent(in) :: n_halo
      character(*), optional, intent(in) :: from_to, bc_start, bc_end
      logical, optional, intent(in) :: sym
      real(dp), optional, intent(in) :: c_nu, nu0_nu
    end subroutine alloc_tdsops
  end interface

  abstract interface
    subroutine init_poisson_fft(self, mesh, xdirps, ydirps, zdirps)
      import :: base_backend_t
      import :: dirps_t
      import :: mesh_t
      implicit none

      class(base_backend_t) :: self
      class(mesh_t), intent(in) :: mesh
      type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    end subroutine init_poisson_fft
  end interface

contains

  subroutine base_init(self)
    implicit none

    class(base_backend_t) :: self

  end subroutine base_init

  subroutine get_field_data(self, data, f, dir)
   !! Extract data from field `f` optionally reordering into `dir` orientation.
   !! To output in same orientation as `f`, use `call ...%get_field_data(data, f, f%dir)`
    implicit none

    class(base_backend_t) :: self
    real(dp), dimension(:, :, :), intent(out) :: data !! Output array
    class(field_t), intent(in) :: f !! Field
    integer, optional, intent(in) :: dir !! Desired orientation of output array (defaults to Cartesian)

    class(field_t), pointer :: f_temp
    integer :: direction, rdr_dir

    if (present(dir)) then
      direction = dir
    else
      direction = DIR_C
    end if

    ! Returns 0 if no reorder required
    rdr_dir = get_rdr_from_dirs(f%dir, direction)

    ! Carry out a reorder if we need, and copy from field to data array
    if (rdr_dir /= 0) then
      f_temp => self%allocator%get_block(direction)
      call self%reorder(f_temp, f, rdr_dir)
      call self%copy_f_to_data(data, f_temp)
      call self%allocator%release_block(f_temp)
    else
      call self%copy_f_to_data(data, f)
    end if

  end subroutine get_field_data

  subroutine set_field_data(self, f, data, dir)
    implicit none

    class(base_backend_t) :: self
    class(field_t), intent(inout) :: f !! Field
    real(dp), dimension(:, :, :), intent(in) :: data !! Input array
    integer, optional, intent(in) :: dir !! Orientation of input array (defaults to Cartesian)

    class(field_t), pointer :: f_temp
    integer :: direction, rdr_dir

    if (present(dir)) then
      direction = dir
    else
      direction = DIR_C
    end if

    ! Returns 0 if no reorder required
    rdr_dir = get_rdr_from_dirs(direction, f%dir)

    ! Carry out a reorder if we need, and copy from data array to field
    if (rdr_dir /= 0) then
      f_temp => self%allocator%get_block(direction)
      call self%copy_data_to_f(f_temp, data)
      call self%reorder(f, f_temp, rdr_dir)
      call self%allocator%release_block(f_temp)
    else
      call self%copy_data_to_f(f, data)
    end if

  end subroutine set_field_data

end module m_base_backend
