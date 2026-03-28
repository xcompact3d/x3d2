module m_monitoring
  !! Computes and logs global scalar monitoring quantities to a file.
  !!
  !! Tracks the following quantities:
  !!
  !! - Kinetic energy: \( E_k = \frac{1}{2N} \sum (u^2 + v^2 + w^2) \)
  !! - Enstrophy: \( \mathcal{E} = \frac{1}{2N} \sum |\nabla \times \mathbf{u}|^2 \)
  !! - Dissipation rate: \( \varepsilon = -\frac{\nu}{N} \sum \mathbf{u} \cdot \nabla^2 \mathbf{u} \)
  !! - Divergence: \( \max |\nabla \cdot \mathbf{u}| \) and mean (divergence-free check)

  use m_common, only: dp, DIR_X, DIR_Z, VERT
  use m_field, only: field_t
  use m_solver, only: solver_t

  implicit none

  type :: monitoring_t
    integer, private :: file_unit = -1
    logical, private :: is_root = .false.
  contains
    procedure :: init
    procedure :: write_step
    procedure :: finalise
  end type monitoring_t

contains

  subroutine init(self, solver)
    class(monitoring_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver

    self%is_root = solver%mesh%par%is_root()

    if (self%is_root) then
      open (newunit=self%file_unit, file='monitoring.csv', &
            status='replace', action='write')
      write (self%file_unit, '(A)') &
        '# time, enstrophy, ke, eps, div_u_max, div_u_mean'
    end if

  end subroutine init

  subroutine write_step(self, solver, t, u, v, w)
    !! Computes and reports monitoring quantities, writing to both stdout
    !! and the monitoring file.
    class(monitoring_t), intent(inout) :: self
    class(solver_t), intent(inout) :: solver
    real(dp), intent(in) :: t
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: du, dv, dw, div_u, lapl
    real(dp) :: enstrophy, ke, eps, div_u_max, div_u_mean

    !! Kinetic energy: \( E_k = \frac{1}{2N} \sum (u^2 + v^2 + w^2) \)
    ke = 0.5_dp*(solver%backend%scalar_product(u, u) &
                 + solver%backend%scalar_product(v, v) &
                 + solver%backend%scalar_product(w, w)) &
         /solver%ngrid

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

    !! Dissipation rate (2nd derivative form):
    !! \( \varepsilon = -\frac{\nu}{N} \sum \mathbf{u} \cdot \nabla^2 \mathbf{u} \)
    !! Equivalent to the strain-rate form for periodic flows.
    lapl => solver%backend%allocator%get_block(DIR_X, VERT)

    call solver%vector_calculus%laplacian( &
      lapl, u, solver%xdirps%der2nd, solver%ydirps%der2nd, &
      solver%zdirps%der2nd &
      )
    eps = solver%backend%scalar_product(u, lapl)

    call solver%vector_calculus%laplacian( &
      lapl, v, solver%xdirps%der2nd, solver%ydirps%der2nd, &
      solver%zdirps%der2nd &
      )
    eps = eps + solver%backend%scalar_product(v, lapl)

    call solver%vector_calculus%laplacian( &
      lapl, w, solver%xdirps%der2nd, solver%ydirps%der2nd, &
      solver%zdirps%der2nd &
      )
    eps = eps + solver%backend%scalar_product(w, lapl)

    call solver%backend%allocator%release_block(lapl)

    eps = -solver%nu*eps/solver%ngrid

    !! Divergence: \( \nabla \cdot \mathbf{u} \) max and mean (should be ~0)
    div_u => solver%backend%allocator%get_block(DIR_Z)

    call solver%divergence_v2p(div_u, u, v, w)

    call solver%backend%field_max_mean(div_u_max, div_u_mean, div_u)

    call solver%backend%allocator%release_block(div_u)

    ! Print to stdout and write to file (root only)
    if (self%is_root) then
      print *, 'enstrophy:', enstrophy
      print *, 'div u max mean:', div_u_max, div_u_mean

      write (self%file_unit, '(ES20.12,5(",",ES20.12))') &
        t, enstrophy, ke, eps, div_u_max, div_u_mean
      flush (self%file_unit)
    end if

  end subroutine write_step

  subroutine finalise(self)
    class(monitoring_t), intent(inout) :: self

    if (self%is_root .and. self%file_unit /= -1) then
      close (self%file_unit)
      self%file_unit = -1
    end if

  end subroutine finalise

end module m_monitoring
