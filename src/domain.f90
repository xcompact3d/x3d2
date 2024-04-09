module m_domain
  use mpi

  use m_common, only: dp
  use m_tdsops, only: dirps_t

  implicit none

contains

  subroutine domain_decomposition(xdirps, ydirps, zdirps, nrank, nproc)
    !! Supports 1D, 2D, and 3D domain decomposition.
    !!
    !! Current implementation allows only constant sub-domain size across a
    !! given direction.
    type(dirps_t), intent(inout) :: xdirps, ydirps, zdirps
    integer, intent(in) :: nrank, nproc

    integer, allocatable, dimension(:, :, :) :: global_ranks
    integer :: subd_pos(3), i, iprev, inext

    ! A 3D array corresponding to each region in the global domain
    allocate (global_ranks(xdirps%nproc, ydirps%nproc, zdirps%nproc))

    ! set the corresponding global rank for each sub-domain
    global_ranks = reshape([(i, i=0, nproc - 1)], &
                           shape=[xdirps%nproc, ydirps%nproc, zdirps%nproc])

    ! subdomain position in the global domain
    subd_pos = findloc(global_ranks, nrank)

    ! local/directional position of the subdomain
    xdirps%nrank = subd_pos(1) - 1
    ydirps%nrank = subd_pos(2) - 1
    zdirps%nrank = subd_pos(3) - 1

    xdirps%n_offset = xdirps%n*xdirps%nrank
    ydirps%n_offset = ydirps%n*ydirps%nrank
    zdirps%n_offset = zdirps%n*zdirps%nrank

    iprev = modulo(subd_pos(1) - 2, xdirps%nproc) + 1
    inext = modulo(subd_pos(1) - xdirps%nproc, xdirps%nproc) + 1
    xdirps%pprev = global_ranks(iprev, subd_pos(2), subd_pos(3))
    xdirps%pnext = global_ranks(inext, subd_pos(2), subd_pos(3))

    iprev = modulo(subd_pos(2) - 2, ydirps%nproc) + 1
    inext = modulo(subd_pos(2) - ydirps%nproc, ydirps%nproc) + 1
    ydirps%pprev = global_ranks(subd_pos(1), iprev, subd_pos(3))
    ydirps%pnext = global_ranks(subd_pos(1), inext, subd_pos(3))

    iprev = modulo(subd_pos(3) - 2, zdirps%nproc) + 1
    inext = modulo(subd_pos(3) - zdirps%nproc, zdirps%nproc) + 1
    zdirps%pprev = global_ranks(subd_pos(1), subd_pos(2), iprev)
    zdirps%pnext = global_ranks(subd_pos(1), subd_pos(2), inext)

  end subroutine domain_decomposition

end module m_domain
