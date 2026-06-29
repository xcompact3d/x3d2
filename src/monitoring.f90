module m_monitoring
  !! Computes and logs global scalar monitoring quantities to a file.
  !!
  !! Tracks the following quantities:
  !!
  !! - Enstrophy: \( \mathcal{E} = \frac{1}{2N} \sum |\nabla \times \mathbf{u}|^2 \)
  !! - Divergence: \( \max |\nabla \cdot \mathbf{u}| \) and mean (divergence-free check)

  use m_common, only: dp, DIR_X, DIR_Z, VERT
  use m_field, only: field_t
  use m_scalar_series, only: scalar_series_t
  use m_solver, only: solver_t

  implicit none

  type :: monitoring_t
    logical, private :: is_root = .false.
    type(scalar_series_t), private :: series
  contains
    procedure :: init
    procedure :: write_step
    procedure :: finalise
  end type monitoring_t

contains

  subroutine init(self, solver, append)
    class(monitoring_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    logical, intent(in), optional :: append

    logical :: append_output
    character(len=16), parameter :: columns(3) = &
                                    ['enstrophy      ', &
                                     'div_u_max      ', &
                                     'div_u_mean     ']

    self%is_root = solver%mesh%par%is_root()
    append_output = .false.
    if (present(append)) append_output = append
    call self%series%init('monitoring.csv', &
                          columns, self%is_root, append_output)

  end subroutine init

  subroutine write_step(self, solver, t, u, v, w)
    !! Computes and reports monitoring quantities, writing to both stdout
    !! and the monitoring file.
    class(monitoring_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    real(dp), intent(in) :: t
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: du, dv, dw, div_u
    real(dp) :: enstrophy, div_u_max, div_u_mean

    !! Enstrophy: \( \mathcal{E} = \frac{1}{2N} \sum |\nabla \times \mathbf{u}|^2 \)
    du => solver%backend%allocator%get_block(DIR_X, VERT)
    dv => solver%backend%allocator%get_block(DIR_X, VERT)
    dw => solver%backend%allocator%get_block(DIR_X, VERT)

    call solver%curl(du, dv, dw, u, v, w)

    enstrophy = 0.5_dp*(solver%backend%scalar_product(du, du) &
                        + solver%backend%scalar_product(dv, dv) &
                        + solver%backend%scalar_product(dw, dw)) &
                /solver%ngrid

    call solver%backend%allocator%release_block(du)
    call solver%backend%allocator%release_block(dv)
    call solver%backend%allocator%release_block(dw)

    !! Divergence: \( \nabla \cdot \mathbf{u} \) max and mean (should be ~0)
    div_u => solver%backend%allocator%get_block(DIR_Z)

    call solver%divergence_v2p(div_u, u, v, w)

    call solver%backend%field_max_mean(div_u_max, div_u_mean, div_u)

    call solver%backend%allocator%release_block(div_u)

    ! Print to stdout and write to file (root only)
    if (self%is_root) then
      print *, 'enstrophy:', enstrophy
      print *, 'div u max mean:', div_u_max, div_u_mean

      call self%series%write_step(t, [enstrophy, div_u_max, div_u_mean])
    end if

  end subroutine write_step

  subroutine finalise(self)
    class(monitoring_t), intent(inout) :: self

    call self%series%finalise()

  end subroutine finalise

end module m_monitoring
