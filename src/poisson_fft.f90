module m_poisson_fft
  use m_common, only: dp, pi, CELL
  use m_field, only: field_t
  use m_mesh, only: mesh_t, geo_t
  use m_tdsops, only: dirps_t

  implicit none

  type, abstract :: poisson_fft_t
    !! FFT based Poisson solver
    type(mesh_t), pointer :: mesh => null()
    !> Global dimensions
    integer :: nx_glob, ny_glob, nz_glob
    !> Local dimensions
    integer :: nx_loc, ny_loc, nz_loc
    !> Local dimensions in the permuted slabs
    integer :: nx_perm, ny_perm, nz_perm
    !> Local dimensions in the permuted slabs in spectral space
    integer :: nx_spec, ny_spec, nz_spec
    !> Offset in y and z directions in the permuted slabs in spectral space
    integer :: x_sp_st, y_sp_st, z_sp_st
    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), allocatable, dimension(:, :, :) :: waves
    !> Wave numbers in x, y, and z
    real(dp), allocatable, dimension(:) :: ax, bx, ay, by, az, bz
    !> Wave numbers in x, y, and z
    complex(dp), allocatable, dimension(:) :: kx, ky, kz, exs, eys, ezs, &
                                              k2x, k2y, k2z
    !> Staggared grid transformation
    real(dp), allocatable, dimension(:) :: trans_x_re, trans_x_im, &
                                           trans_y_re, trans_y_im, &
                                           trans_z_re, trans_z_im
    !> Periodicity in x, y, and z
    logical :: periodic_x, periodic_y, periodic_z, &
               stretched_y = .false., stretched_y_sym
    !> Stretching operator matrices
    real(dp), allocatable, dimension(:, :, :, :) :: a_odd_re, a_odd_im, &
                                                    a_even_re, a_even_im, &
                                                    a_re, a_im
    !> lowmem option, only used in CUDA backend
    logical :: lowmem = .false.
    !> Procedure pointer to BC specific poisson solvers
    procedure(poisson_xxx), pointer :: poisson => null()
  contains
    procedure(fft_forward), deferred :: fft_forward_010
    procedure(fft_forward), deferred :: fft_forward_100
    procedure(fft_forward), deferred :: fft_forward_110
    procedure(fft_forward), deferred :: fft_forward
    procedure(fft_backward), deferred :: fft_backward_010
    procedure(fft_backward), deferred :: fft_backward_100
    procedure(fft_backward), deferred :: fft_backward_110
    procedure(fft_backward), deferred :: fft_backward
    procedure(fft_postprocess), deferred :: fft_postprocess_000
    procedure(fft_postprocess), deferred :: fft_postprocess_010
    procedure(fft_postprocess), deferred :: fft_postprocess_100
    procedure(fft_postprocess), deferred :: fft_postprocess_110
    procedure(field_process), deferred :: enforce_periodicity_x
    procedure(field_process), deferred :: undo_periodicity_x
    procedure(field_process), deferred :: enforce_periodicity_y
    procedure(field_process), deferred :: undo_periodicity_y
    procedure :: base_init
    procedure :: solve_poisson
    procedure :: stretching_matrix
    procedure :: waves_set
    procedure :: get_km
    procedure :: get_km_re
    procedure :: get_km_im
  end type poisson_fft_t

  abstract interface
    subroutine fft_forward(self, f_in)
      import :: poisson_fft_t
      import :: field_t
      implicit none

      class(poisson_fft_t) :: self
      class(field_t), intent(in) :: f_in
    end subroutine fft_forward

    subroutine fft_backward(self, f_out)
      import :: poisson_fft_t
      import :: field_t
      implicit none

      class(poisson_fft_t) :: self
      class(field_t), intent(inout) :: f_out
    end subroutine fft_backward

    subroutine fft_postprocess(self)
      import :: poisson_fft_t
      implicit none

      class(poisson_fft_t) :: self
    end subroutine fft_postprocess
  end interface

  abstract interface
    subroutine poisson_xxx(self, f, temp)
      import :: poisson_fft_t
      import :: field_t

      class(poisson_fft_t) :: self
      class(field_t), intent(inout) :: f, temp
    end subroutine poisson_xxx

    subroutine field_process(self, f_out, f_in)
      import :: poisson_fft_t
      import :: field_t

      class(poisson_fft_t) :: self
      class(field_t), intent(inout) :: f_out
      class(field_t), intent(in) :: f_in
    end subroutine field_process
  end interface

contains

  subroutine base_init(self, mesh, xdirps, ydirps, zdirps, n_spec, n_sp_st)
    implicit none
    class(poisson_fft_t) :: self
    type(mesh_t), intent(in), target :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    integer, dimension(3), intent(in) :: n_spec ! Size of the spectral pencil
    integer, dimension(3), intent(in) :: n_sp_st ! Offset of the spectral pencil

    integer :: dims(3)
    self%mesh => mesh
    ! Decomposition is in y- and z-directions
    if (mesh%par%nproc_dir(1) /= 1) print *, 'nproc_dir in x-dir must be 1'

    dims = mesh%get_global_dims(CELL)
    self%nx_glob = dims(1); self%ny_glob = dims(2); self%nz_glob = dims(3)
    dims = mesh%get_dims(CELL)
    self%nx_loc = dims(1); self%ny_loc = dims(2); self%nz_loc = dims(3)

    self%nx_spec = n_spec(1)
    self%ny_spec = n_spec(2)
    self%nz_spec = n_spec(3)

    self%periodic_x = mesh%grid%periodic_BC(1)
    self%periodic_y = mesh%grid%periodic_BC(2)
    self%periodic_z = mesh%grid%periodic_BC(3)

    self%x_sp_st = n_sp_st(1)
    self%y_sp_st = n_sp_st(2)
    self%z_sp_st = n_sp_st(3)

    allocate (self%ax(self%nx_glob), self%bx(self%nx_glob))
    allocate (self%ay(self%ny_glob), self%by(self%ny_glob))
    allocate (self%az(self%nz_glob), self%bz(self%nz_glob))

    allocate (self%kx(self%nx_glob), self%k2x(self%nx_glob))
    allocate (self%ky(self%ny_glob), self%k2y(self%ny_glob))
    allocate (self%kz(self%nz_glob), self%k2z(self%nz_glob))
    allocate (self%exs(self%nx_glob))
    allocate (self%eys(self%ny_glob))
    allocate (self%ezs(self%nz_glob))

    allocate (self%trans_x_re(self%nx_spec), self%trans_x_im(self%nx_spec))
    allocate (self%trans_y_re(self%ny_spec), self%trans_y_im(self%ny_spec))
    allocate (self%trans_z_re(self%nz_spec), self%trans_z_im(self%nz_spec))

    allocate (self%waves(self%nx_spec, self%ny_spec, self%nz_spec))

    ! waves_set requires some of the preprocessed tdsops variables.
    call self%waves_set(mesh%geo, xdirps, ydirps, zdirps)

    if (mesh%geo%stretched(1) .or. mesh%geo%stretched(3)) then
      error stop 'FFT based Poisson solver does not support stretching in x-&
                  & or z-directions!'
    end if

    ! use correct procedure based on BCs
    if (self%periodic_x .and. self%periodic_y .and. self%periodic_z) then
      self%poisson => poisson_000
    else if (self%periodic_x .and. (.not. self%periodic_y) &
             .and. (self%periodic_z)) then
      if (mesh%par%nproc > 1) then
        error stop 'Multiple ranks are not yet supported for non-periodic BCs!'
      end if
      self%poisson => poisson_010
      ! stretching requires some coefficients matrices
      if (mesh%geo%stretched(2)) then
        self%stretched_y = .true.
        call self%stretching_matrix(mesh%geo, xdirps, ydirps, zdirps)
      end if
    else if ((.not. self%periodic_x) .and. (self%periodic_y) &
             .and. (self%periodic_z)) then
      if (mesh%par%nproc > 1) then
        error stop 'Multiple ranks are not yet supported for non-periodic BCs!'
      end if

      self%poisson => poisson_100
    else if ((.not. self%periodic_x) .and. (.not. self%periodic_y) &
             .and. (self%periodic_z)) then
      if (mesh%par%nproc > 1) then
        error stop 'Multiple ranks are not yet supported for non-periodic BCs!'
      end if
      self%poisson => poisson_110

    else
      error stop 'Requested BCs are not supported in FFT-based Poisson solver!'
    end if
  end subroutine base_init

  subroutine solve_poisson(self, f, temp)
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%poisson(f, temp)

  end subroutine solve_poisson

  subroutine poisson_000(self, f, temp)
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%fft_forward(f)
    call self%fft_postprocess_000
    call self%fft_backward(f)

  end subroutine poisson_000

  subroutine poisson_010(self, f, temp)
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%enforce_periodicity_y(temp, f)

    call self%fft_forward_010(temp)
    call self%fft_postprocess_010
    call self%fft_backward_010(temp)

    call self%undo_periodicity_y(f, temp)

  end subroutine poisson_010

  subroutine poisson_100(self, f, temp)
    implicit none
    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%enforce_periodicity_x(temp, f)

    call self%fft_forward_100(temp)
    call self%fft_postprocess_100
    call self%fft_backward_100(temp)

    call self%undo_periodicity_x(f, temp)

  end subroutine poisson_100

  subroutine poisson_110(self, f, temp)
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f, temp

    ! Apply periodicity enforcement for both X and Y
    call self%enforce_periodicity_x(temp, f)
    call self%enforce_periodicity_y(f, temp)

    call self%fft_forward_110(f)
    call self%fft_postprocess_110
    call self%fft_backward_110(f)

    ! Undo periodicity for both X and Y
    call self%undo_periodicity_y(temp, f)
    call self%undo_periodicity_x(f, temp)

  end subroutine poisson_110

  subroutine stretching_matrix(self, geo, xdirps, ydirps, zdirps)
    !! Stretching necessitates a special operation in spectral space.
    !! The coefficients for the operation are stored in matrix form.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 5
    implicit none

    class(poisson_fft_t) :: self
    type(geo_t), intent(in) :: geo
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    real(dp) :: temp, a0, a1, c1_od, c2_od, c1_ev, c2_ev
    complex(dp) :: km_a1, km_a1_od, km_a1_ev
    integer :: i, j, k, ix, iy, iz, iy_od, iy_ev

    do i = 1, self%nx_spec
      temp = real(self%exs(i), kind=dp)*geo%d(1)
      self%trans_x_re(i) = 2*(xdirps%interpl_v2p%a*cos(temp*0.5_dp) &
                              + xdirps%interpl_v2p%b*cos(temp*1.5_dp) &
                              + xdirps%interpl_v2p%c*cos(temp*2.5_dp) &
                              + xdirps%interpl_v2p%d*cos(temp*3.5_dp)) &
                           /(1._dp + 2*xdirps%interpl_v2p%alpha*cos(temp))
      self%trans_x_im(i) = self%trans_x_re(i)
    end do

    do j = 1, self%ny_spec
      temp = real(self%eys(j), kind=dp)*geo%d(2)
      self%trans_y_re(j) = 2*(ydirps%interpl_v2p%a*cos(temp*0.5_dp) &
                              + ydirps%interpl_v2p%b*cos(temp*1.5_dp) &
                              + ydirps%interpl_v2p%c*cos(temp*2.5_dp) &
                              + ydirps%interpl_v2p%d*cos(temp*3.5_dp)) &
                           /(1._dp + 2*ydirps%interpl_v2p%alpha*cos(temp))
      self%trans_y_im(j) = self%trans_y_re(j)
    end do

    do k = 1, self%nz_spec
      temp = real(self%ezs(k), kind=dp)*geo%d(3)
      self%trans_z_re(k) = 2*(zdirps%interpl_v2p%a*cos(temp*0.5_dp) &
                              + zdirps%interpl_v2p%b*cos(temp*1.5_dp) &
                              + zdirps%interpl_v2p%c*cos(temp*2.5_dp) &
                              + zdirps%interpl_v2p%d*cos(temp*3.5_dp)) &
                           /(1._dp + 2*zdirps%interpl_v2p%alpha*cos(temp))
      self%trans_z_im(k) = self%trans_z_re(k)
    end do

    if (trim(geo%stretching(2)) == 'bottom') then
      self%stretched_y_sym = .false.
      allocate (self%a_re(self%nx_spec, self%ny_spec, self%nz_spec, 5))
      allocate (self%a_im(self%nx_spec, self%ny_spec, self%nz_spec, 5))

      a0 = (geo%alpha(2)/pi + 1._dp/(2*pi*geo%beta(2)))*geo%L(2)
      a1 = -1._dp/(4*pi*geo%beta(2))*geo%L(2)

      ! diagonal
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st
            if (iy == 1) then
              km_a1 = self%get_km(ix, 2, iz)
            else if (iy == self%ny_spec) then
              km_a1 = self%get_km(ix, self%ny_spec - 1, iz)
            else
              km_a1 = self%get_km(ix, iy - 1, iz) + self%get_km(ix, iy + 1, iz)
            end if

            self%a_re(i, j, k, 3) = &
              -(get_real(self%kx(ix)) &
                *self%trans_y_re(iy)*self%trans_z_re(iz))**2 &
              - (get_real(self%kz(iz)) &
                 *self%trans_y_re(iy)*self%trans_x_re(ix))**2 &
              - a0**2*self%get_km_re(ix, iy, iz)**2 &
              - a1**2*self%get_km_re(ix, iy, iz)*get_real(km_a1)
            self%a_im(i, j, k, 3) = &
              -(get_imag(self%kx(ix)) &
                *self%trans_y_im(iy)*self%trans_z_im(iz))**2 &
              - (get_imag(self%kz(iz)) &
                 *self%trans_y_im(iy)*self%trans_x_im(ix))**2 &
              - a0**2*self%get_km_im(ix, iy, iz)**2 &
              - a1**2*self%get_km_im(ix, iy, iz)*get_imag(km_a1)
          end do
        end do
      end do

      ! diagonal + 1
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st

            self%a_re(i, j, k, 4) = &
              a0*a1*self%get_km_re(ix, iy + 1, iz) &
              *(self%get_km_re(ix, iy, iz) + self%get_km_re(ix, iy + 1, iz))
            self%a_im(i, j, k, 4) = &
              a0*a1*self%get_km_im(ix, iy + 1, iz) &
              *(self%get_km_im(ix, iy, iz) + self%get_km_im(ix, iy + 1, iz))
          end do
        end do
      end do

      ! diagonal + 2
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec - 2
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st

            self%a_re(i, j, k, 5) = -a1*a1*self%get_km_re(ix, iy + 1, iz) &
                                    *self%get_km_re(ix, iy + 2, iz)
            self%a_im(i, j, k, 5) = -a1*a1*self%get_km_im(ix, iy + 1, iz) &
                                    *self%get_km_im(ix, iy + 2, iz)
          end do
        end do
      end do

      ! diagonal - 1
      do k = 1, self%nz_spec
        do j = 2, self%ny_spec
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st

            self%a_re(i, j, k, 2) = a0*a1*self%get_km_re(ix, iy - 1, iz) &
                                    *(self%get_km_re(ix, iy, iz) &
                                      + self%get_km_re(ix, iy - 1, iz))
            self%a_im(i, j, k, 2) = a0*a1*self%get_km_im(ix, iy - 1, iz) &
                                    *(self%get_km_im(ix, iy, iz) &
                                      + self%get_km_im(ix, iy - 1, iz))
          end do
        end do
      end do

      ! diagonal - 2
      do k = 1, self%nz_spec
        do j = 3, self%ny_spec
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st

            self%a_re(i, j, k, 1) = -a1*a1*self%get_km_re(ix, iy - 1, iz) &
                                    *self%get_km_re(ix, iy - 2, iz)
            self%a_im(i, j, k, 1) = -a1*a1*self%get_km_im(ix, iy - 1, iz) &
                                    *self%get_km_im(ix, iy - 2, iz)
          end do
        end do
      end do

      ! tweak the matrix to make it not singular
      self%a_re(1, 1, 1, 3) = 1._dp; self%a_im(1, 1, 1, 3) = 1._dp
      self%a_re(1, 1, 1, 4) = 0; self%a_im(1, 1, 1, 4) = 0
      self%a_re(1, 1, 1, 5) = 0; self%a_im(1, 1, 1, 5) = 0
    else
      self%stretched_y_sym = .true.
      allocate (self%a_odd_re(self%nx_spec, self%ny_spec/2, self%nz_spec, 5))
      allocate (self%a_odd_im(self%nx_spec, self%ny_spec/2, self%nz_spec, 5))
      allocate (self%a_even_re(self%nx_spec, self%ny_spec/2, self%nz_spec, 5))
      allocate (self%a_even_im(self%nx_spec, self%ny_spec/2, self%nz_spec, 5))

      self%a_odd_re = 0._dp
      self%a_odd_im = 0._dp
      self%a_even_re = 0._dp
      self%a_even_im = 0._dp

      a0 = (geo%alpha(2)/pi + 1._dp/(2*pi*geo%beta(2)))*geo%L(2)
      select case (trim(geo%stretching(2)))
      case ('centred')
        a1 = 1._dp/(4*pi*geo%beta(2))*geo%L(2)
      case ('top-bottom')
        a1 = -1._dp/(4*pi*geo%beta(2))*geo%L(2)
      case default
        a1 = 0._dp
      end select

      ! diagonal
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec/2
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st
            iy_od = 2*iy - 1
            iy_ev = 2*iy
            c1_od = a0*a0
            c2_od = a1*a1
            c1_ev = a0*a0
            c2_ev = a1*a1
            if (iy == 1) then
              c1_ev = a0*a0 - a1*a1
              km_a1_od = self%get_km(ix, 3, iz)
              km_a1_ev = self%get_km(ix, 4, iz)
            else if (iy == self%ny_spec/2) then
              c1_ev = (a0 + a1)*(a0 + a1)
              km_a1_od = self%get_km(ix, iy_od - 2, iz)
              km_a1_ev = self%get_km(ix, iy_ev - 2, iz)
            else
              km_a1_od = self%get_km(ix, iy_od - 2, iz) &
                         + self%get_km(ix, iy_od + 2, iz)
              km_a1_ev = self%get_km(ix, iy_ev - 2, iz) &
                         + self%get_km(ix, iy_ev + 2, iz)
            end if

            self%a_odd_re(i, j, k, 3) = &
              -(get_real(self%kx(ix)) &
                *self%trans_y_re(iy_od)*self%trans_z_re(iz))**2 &
              - (get_real(self%kz(iz)) &
                 *self%trans_y_re(iy_od)*self%trans_x_re(ix))**2 &
              - c1_od*self%get_km_re(ix, iy_od, iz)**2 &
              - c2_od*self%get_km_re(ix, iy_od, iz)*get_real(km_a1_od)
            self%a_odd_im(i, j, k, 3) = &
              -(get_imag(self%kx(ix)) &
                *self%trans_y_im(iy_od)*self%trans_z_im(iz))**2 &
              - (get_imag(self%kz(iz)) &
                 *self%trans_y_im(iy_od)*self%trans_x_im(ix))**2 &
              - c1_od*self%get_km_im(ix, iy_od, iz)**2 &
              - c2_od*self%get_km_im(ix, iy_od, iz)*get_imag(km_a1_od)
            self%a_even_re(i, j, k, 3) = &
              -(get_real(self%kx(ix)) &
                *self%trans_y_re(iy_ev)*self%trans_z_re(iz))**2 &
              - (get_real(self%kz(iz)) &
                 *self%trans_y_re(iy_ev)*self%trans_x_re(ix))**2 &
              - c1_ev*self%get_km_re(ix, iy_ev, iz)**2 &
              - c2_ev*self%get_km_re(ix, iy_ev, iz)*get_real(km_a1_ev)
            self%a_even_im(i, j, k, 3) = &
              -(get_imag(self%kx(ix)) &
                *self%trans_y_im(iy_ev)*self%trans_z_im(iz))**2 &
              - (get_imag(self%kz(iz)) &
                 *self%trans_y_im(iy_ev)*self%trans_x_im(ix))**2 &
              - c1_ev*self%get_km_im(ix, iy_ev, iz)**2 &
              - c2_ev*self%get_km_im(ix, iy_ev, iz)*get_imag(km_a1_ev)
          end do
        end do
      end do

      ! diagonal + 1
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec/2
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st
            iy_od = 2*iy - 1
            iy_ev = 2*iy
            c1_od = a0*a1
            c2_od = a0*a1
            c1_ev = a0*a1
            c2_ev = a0*a1
            if (iy == 1) then
              c1_od = 2*a0*a1
              c2_od = 2*a0*a1
              c1_ev = a0*a1 - a1*a1
              c2_ev = a0*a1
            else if (iy == self%ny_spec/2 - 1) then
              c1_ev = a0*a1
              c2_ev = (a0 + a1)*a1
            else if (iy == self%ny_spec/2) then
              c1_ev = 0; c2_ev = 0
            end if
            self%a_odd_re(i, j, k, 4) = &
              c1_od*(self%get_km_re(ix, iy_od, iz) &
                     *self%get_km_re(ix, iy_od + 2, iz)) &
              + c2_od*self%get_km_re(ix, iy_od + 2, iz)**2
            self%a_odd_im(i, j, k, 4) = &
              c1_od*(self%get_km_im(ix, iy_od, iz) &
                     *self%get_km_im(ix, iy_od + 2, iz)) &
              + c2_od*self%get_km_im(ix, iy_od + 2, iz)**2
            self%a_even_re(i, j, k, 4) = &
              c1_ev*(self%get_km_re(ix, iy_ev, iz) &
                     *self%get_km_re(ix, iy_ev + 2, iz)) &
              + c2_ev*self%get_km_re(ix, iy_ev + 2, iz)**2
            self%a_even_im(i, j, k, 4) = &
              c1_ev*(self%get_km_im(ix, iy_ev, iz) &
                     *self%get_km_im(ix, iy_ev + 2, iz)) &
              + c2_ev*self%get_km_im(ix, iy_ev + 2, iz)**2
          end do
        end do
      end do

      ! diagonal + 2
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec/2 - 2
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st
            iy_od = 2*iy - 1
            iy_ev = 2*iy
            c1_od = a1*a1
            c1_ev = a1*a1
            if (iy == 1) then
              c1_od = 2*a1*a1
            end if
            self%a_odd_re(i, j, k, 5) = &
              -(c1_od*self%get_km_re(ix, iy_od + 2, iz) &
                *self%get_km_re(ix, iy_od + 4, iz))
            self%a_odd_im(i, j, k, 5) = &
              -(c1_od*self%get_km_im(ix, iy_od + 2, iz) &
                *self%get_km_im(ix, iy_od + 4, iz))
            self%a_even_re(i, j, k, 5) = &
              -(c1_ev*self%get_km_re(ix, iy_ev + 2, iz) &
                *self%get_km_re(ix, iy_ev + 4, iz))
            self%a_even_im(i, j, k, 5) = &
              -(c1_ev*self%get_km_im(ix, iy_ev + 2, iz) &
                *self%get_km_im(ix, iy_ev + 4, iz))
          end do
        end do
      end do

      ! diagonal - 1
      do k = 1, self%nz_spec
        do j = 2, self%ny_spec/2
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st
            iy_od = 2*iy - 1
            iy_ev = 2*iy
            c1_od = a0*a1
            c2_od = a0*a1
            c1_ev = a0*a1
            c2_ev = a0*a1
            if (iy == 1) then
              c1_od = 0; c2_od = 0; c1_ev = 0; c2_ev = 0
            else if (iy == 2) then
              c1_ev = a0*a1
              c2_ev = (a0 + a1)*a1
            else if (iy == self%ny_spec/2) then
              c1_ev = (a0 + a1)*a1
              c2_ev = a0*a1
            end if
            self%a_odd_re(i, j, k, 2) = &
              c1_od*(self%get_km_re(ix, iy_od, iz) &
                     *self%get_km_re(ix, iy_od - 2, iz)) &
              + c2_od*self%get_km_re(ix, iy_od - 2, iz)**2
            self%a_odd_im(i, j, k, 2) = &
              c1_od*(self%get_km_im(ix, iy_od, iz) &
                     *self%get_km_im(ix, iy_od - 2, iz)) &
              + c2_od*self%get_km_im(ix, iy_od - 2, iz)**2
            self%a_even_re(i, j, k, 2) = &
              c1_ev*(self%get_km_re(ix, iy_ev, iz) &
                     *self%get_km_re(ix, iy_ev - 2, iz)) &
              + c2_ev*self%get_km_re(ix, iy_ev - 2, iz)**2
            self%a_even_im(i, j, k, 2) = &
              c1_ev*(self%get_km_im(ix, iy_ev, iz) &
                     *self%get_km_im(ix, iy_ev - 2, iz)) &
              + c2_ev*self%get_km_im(ix, iy_ev - 2, iz)**2
          end do
        end do
      end do

      ! diagonal - 2
      do k = 1, self%nz_spec
        do j = 3, self%ny_spec/2
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st
            iy_od = 2*iy - 1
            iy_ev = 2*iy
            self%a_odd_re(i, j, k, 1) = &
              -(a1*a1*self%get_km_re(ix, iy_od - 2, iz) &
                *self%get_km_re(ix, iy_od - 4, iz))
            self%a_odd_im(i, j, k, 1) = &
              -(a1*a1*self%get_km_im(ix, iy_od - 2, iz) &
                *self%get_km_im(ix, iy_od - 4, iz))
            self%a_even_re(i, j, k, 1) = &
              -(a1*a1*self%get_km_re(ix, iy_ev - 2, iz) &
                *self%get_km_re(ix, iy_ev - 4, iz))
            self%a_even_im(i, j, k, 1) = &
              -(a1*a1*self%get_km_im(ix, iy_ev - 2, iz) &
                *self%get_km_im(ix, iy_ev - 4, iz))
          end do
        end do
      end do

      do k = 1, self%nz_spec
        do i = 1, self%nx_spec
          ix = i + self%x_sp_st; iz = k + self%z_sp_st
          if (get_real(self%k2x(ix)) < 1e-15 &
              .and. get_real(self%k2z(iz)) < 1e-15) then
            self%a_odd_re(i, 1, k, 3) = 1._dp
            self%a_odd_im(i, 1, k, 3) = 1._dp
            self%a_odd_re(i, 1, k, 4) = 0
            self%a_odd_im(i, 1, k, 4) = 0
            self%a_odd_re(i, 1, k, 5) = 0
            self%a_odd_im(i, 1, k, 5) = 0
          end if
        end do
      end do
    end if

  end subroutine stretching_matrix

  subroutine waves_set(self, geo, xdirps, ydirps, zdirps)
  !! Spectral equivalence constants
  !!
  !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    class(poisson_fft_t) :: self
    type(geo_t), intent(in) :: geo
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    integer :: i, j, k, ix, iy, iz
    real(dp) :: rlexs, rleys, rlezs, xtt, ytt, ztt, xt1, yt1, zt1
    complex(dp) :: xt2, yt2, zt2, xyzk

    call wave_numbers( &
      self%ax, self%bx, self%kx, self%exs, self%k2x, &
      self%nx_glob, geo%L(1), geo%d(1), self%periodic_x, &
      xdirps%stagder_v2p%a, xdirps%stagder_v2p%b, xdirps%stagder_v2p%alpha &
      )

    call wave_numbers( &
      self%ay, self%by, self%ky, self%eys, self%k2y, &
      self%ny_glob, geo%L(2), geo%d(2), self%periodic_y, &
      ydirps%stagder_v2p%a, ydirps%stagder_v2p%b, ydirps%stagder_v2p%alpha &
      )

    call wave_numbers( &
      self%az, self%bz, self%kz, self%ezs, self%k2z, &
      self%nz_glob, geo%L(3), geo%d(3), self%periodic_z, &
      zdirps%stagder_v2p%a, zdirps%stagder_v2p%b, zdirps%stagder_v2p%alpha &
      )

    ! Determine which case we're in and compute waves accordingly
    if ((.not. self%periodic_x) .and. self%periodic_y .and. &
        self%periodic_z) then
      ! =========================================================================
      ! 100 case: Non-periodic X, Periodic Y, Periodic Z
      ! Uses TRANSPOSED indexing because data is transposed before FFT
      ! =========================================================================
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec  ! This iterates over X (Dirichlet) after transpose
          do i = 1, self%nx_spec  ! This iterates over Y (periodic, R2C) after transpose
            ! After transpose: array is (ny, nx, nz), R2C gives (ny/2+1, nx, nz)
            ! So i indexes into Y direction, j indexes into X direction
            iy = i + self%y_sp_st  ! Use for ky (first dim after transpose)
            ix = j + self%x_sp_st  ! Use for kx (second dim after transpose)
            iz = k + self%z_sp_st

            rlexs = real(self%exs(ix), kind=dp)*geo%d(1)
            rleys = real(self%eys(iy), kind=dp)*geo%d(2)
            rlezs = real(self%ezs(iz), kind=dp)*geo%d(3)

            xtt = 2*(xdirps%interpl_v2p%a*cos(rlexs*0.5_dp) &
                     + xdirps%interpl_v2p%b*cos(rlexs*1.5_dp) &
                     + xdirps%interpl_v2p%c*cos(rlexs*2.5_dp) &
                     + xdirps%interpl_v2p%d*cos(rlexs*3.5_dp))
            ytt = 2*(ydirps%interpl_v2p%a*cos(rleys*0.5_dp) &
                     + ydirps%interpl_v2p%b*cos(rleys*1.5_dp) &
                     + ydirps%interpl_v2p%c*cos(rleys*2.5_dp) &
                     + ydirps%interpl_v2p%d*cos(rleys*3.5_dp))
            ztt = 2*(zdirps%interpl_v2p%a*cos(rlezs*0.5_dp) &
                     + zdirps%interpl_v2p%b*cos(rlezs*1.5_dp) &
                     + zdirps%interpl_v2p%c*cos(rlezs*2.5_dp) &
                     + zdirps%interpl_v2p%d*cos(rlezs*3.5_dp))

            xt1 = 1._dp + 2*xdirps%interpl_v2p%alpha*cos(rlexs)
            yt1 = 1._dp + 2*ydirps%interpl_v2p%alpha*cos(rleys)
            zt1 = 1._dp + 2*zdirps%interpl_v2p%alpha*cos(rlezs)

            xt2 = self%k2x(ix)*((ytt/yt1)*(ztt/zt1))**2
            yt2 = self%k2y(iy)*((xtt/xt1)*(ztt/zt1))**2
            zt2 = self%k2z(iz)*((xtt/xt1)*(ytt/yt1))**2

            xyzk = xt2 + yt2 + zt2
            self%waves(i, j, k) = xyzk
          end do
        end do
      end do

    else if (self%periodic_z) then
      ! =========================================================================
      ! 000, 010, 110 cases: Periodic Z (standard indexing, no transpose)
      ! 000: Periodic X, Periodic Y, Periodic Z
      ! 010: Periodic X, Non-Periodic Y, Periodic Z
      ! 110: Non-Periodic X, Non-Periodic Y, Periodic Z
      ! =========================================================================
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st
            iy = j + self%y_sp_st
            iz = k + self%z_sp_st

            rlexs = real(self%exs(ix), kind=dp)*geo%d(1)
            rleys = real(self%eys(iy), kind=dp)*geo%d(2)
            rlezs = real(self%ezs(iz), kind=dp)*geo%d(3)

            xtt = 2*(xdirps%interpl_v2p%a*cos(rlexs*0.5_dp) &
                     + xdirps%interpl_v2p%b*cos(rlexs*1.5_dp) &
                     + xdirps%interpl_v2p%c*cos(rlexs*2.5_dp) &
                     + xdirps%interpl_v2p%d*cos(rlexs*3.5_dp))
            ytt = 2*(ydirps%interpl_v2p%a*cos(rleys*0.5_dp) &
                     + ydirps%interpl_v2p%b*cos(rleys*1.5_dp) &
                     + ydirps%interpl_v2p%c*cos(rleys*2.5_dp) &
                     + ydirps%interpl_v2p%d*cos(rleys*3.5_dp))
            ztt = 2*(zdirps%interpl_v2p%a*cos(rlezs*0.5_dp) &
                     + zdirps%interpl_v2p%b*cos(rlezs*1.5_dp) &
                     + zdirps%interpl_v2p%c*cos(rlezs*2.5_dp) &
                     + zdirps%interpl_v2p%d*cos(rlezs*3.5_dp))

            xt1 = 1._dp + 2*xdirps%interpl_v2p%alpha*cos(rlexs)
            yt1 = 1._dp + 2*ydirps%interpl_v2p%alpha*cos(rleys)
            zt1 = 1._dp + 2*zdirps%interpl_v2p%alpha*cos(rlezs)

            xt2 = self%k2x(ix)*((ytt/yt1)*(ztt/zt1))**2
            yt2 = self%k2y(iy)*((xtt/xt1)*(ztt/zt1))**2
            zt2 = self%k2z(iz)*((xtt/xt1)*(ytt/yt1))**2

            xyzk = xt2 + yt2 + zt2
            self%waves(i, j, k) = xyzk
          end do
        end do
      end do

    else if (.not. (self%periodic_x .and. self%periodic_y &
                    .and. self%periodic_z)) then
      ! poisson 111
      error stop 'No support for all non-periodic BCs yet!'

    else
      ! poisson 001, 011, 101
      error stop 'FFT Poisson solver does not support specified BCs!'
    end if

  end subroutine waves_set

  subroutine wave_numbers(a, b, k, e, k2, n, L, d, periodic, c_a, c_b, c_alpha)
    implicit none

    real(dp), dimension(:), intent(out) :: a, b
    complex(dp), dimension(:), intent(out) :: k, e, k2
    integer, intent(in) :: n
    real(dp), intent(in) :: c_a, c_b, c_alpha, L, d
    logical, intent(in) :: periodic

    real(dp) :: w, wp
    integer :: i

    do i = 1, n
      if (periodic) then
        a(i) = sin((i - 1)*pi/n)
        b(i) = cos((i - 1)*pi/n)
      else
        a(i) = sin((i - 1)*pi/2/n)
        b(i) = cos((i - 1)*pi/2/n)
      end if
    end do

    if (periodic) then
      do i = 1, n/2 + 1
        w = 2*pi*(i - 1)/n
        wp = c_a*2*d*sin(0.5_dp*w) + c_b*2*d*sin(1.5_dp*w)
        wp = wp/(1._dp + 2*c_alpha*cos(w))

        k(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*wp/L)
        e(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*w/L)
        k2(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*wp/L)**2
      end do
      do i = n/2 + 2, n
        k(i) = k(n - i + 2)
        e(i) = e(n - i + 2)
        k2(i) = k2(n - i + 2)
      end do
    else
      do i = 1, n
        w = pi*(i - 1)/n
        wp = c_a*2*d*sin(0.5_dp*w) + c_b*2*d*sin(1.5_dp*w)
        wp = wp/(1._dp + 2*c_alpha*cos(w))

        k(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*wp/L)
        e(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*w/L)
        k2(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*wp/L)**2
      end do
    end if

  end subroutine wave_numbers

  real(dp) function get_km_re(self, i, j, k) result(re)
    implicit none

    class(poisson_fft_t) :: self
    integer, intent(in) :: i, j, k

    re = get_real(self%get_km(i, j, k))
  end function get_km_re

  real(dp) function get_km_im(self, i, j, k) result(re)
    implicit none

    class(poisson_fft_t) :: self
    integer, intent(in) :: i, j, k

    re = get_imag(self%get_km(i, j, k))
  end function get_km_im

  complex(dp) function get_km(self, i, j, k) result(km)
    implicit none

    class(poisson_fft_t) :: self
    integer, intent(in) :: i, j, k

    km = cmplx(self%trans_x_re(i)*get_real(self%ky(j))*self%trans_z_re(k), &
               self%trans_x_im(i)*get_imag(self%ky(j))*self%trans_z_im(k), &
               kind=dp)
  end function get_km

  real(dp) function get_real(complx) result(re)
    implicit none

    complex(dp), intent(in) :: complx

    re = real(complx, kind=dp)
  end function get_real

  real(dp) function get_imag(complx) result(im)
    implicit none

    complex(dp), intent(in) :: complx

    im = aimag(complx)
  end function get_imag

end module m_poisson_fft
