module m_memcheck_scratch
  !! Scratch directory handling of the x3d2-memcheck --build path.
  use iso_c_binding, only: c_char, c_int, c_size_t, c_ptr, c_null_char
  use m_memcheck_context, only: memcheck_ctx_t
  implicit none

  interface
    function c_chdir(path) bind(C, name='chdir') result(rc)
      import :: c_char, c_int
      character(kind=c_char), dimension(*), intent(in) :: path
      integer(c_int) :: rc
    end function c_chdir
    function c_getcwd(buf, size) bind(C, name='getcwd') result(ptr)
      import :: c_char, c_size_t, c_ptr
      character(kind=c_char), dimension(*), intent(inout) :: buf
      integer(c_size_t), value :: size
      type(c_ptr) :: ptr
    end function c_getcwd
    function c_getpid() bind(C, name='getpid') result(pid)
      import :: c_int
      integer(c_int) :: pid
    end function c_getpid
  end interface

contains

  function current_dir() result(path)
    !! Working directory via libc getcwd, trimmed at the first embedded NUL
    !! (the C string terminator - Fortran character variables are not
    !! themselves NUL-terminated, so the raw buffer must be cut there before
    !! use). Used by validate_build_inputs to record the invoking directory.
    character(len=4096) :: path
    character(kind=c_char, len=4096) :: buf
    type(c_ptr) :: ptr
    integer :: nul_pos

    ptr = c_getcwd(buf, 4096_c_size_t)
    nul_pos = index(buf, c_null_char)
    if (nul_pos > 0) then
      path = buf(1:nul_pos - 1)
    else
      path = buf
    end if
  end function current_dir

  pure function has_dotdot_component(path) result(found)
    !! True if any '/'-separated component of path is exactly '..' - used
    !! by make_scratch to reject a relative input path it cannot
    !! safely mirror by symlink.
    character(len=*), intent(in) :: path
    logical :: found
    integer :: start, slash_pos

    found = .false.
    start = 1
    do
      slash_pos = index(path(start:), '/')
      if (slash_pos == 0) then
        if (path(start:) == '..') found = .true.
        exit
      end if
      if (path(start:start + slash_pos - 2) == '..') found = .true.
      start = start + slash_pos
    end do
  end function has_dotdot_component

  subroutine run_sh(cmd, status)
    !! Run cmd through the shell; status (optional) receives its exit
    !! status, leave it out for a best effort cleanup command.
    character(len=*), intent(in) :: cmd
    integer, intent(out), optional :: status

    call execute_command_line(cmd, exitstat=status)
  end subroutine run_sh

  function scratch_name(parent) result(path)
    !! <parent>/x3d2-memcheck-build.<pid>, the scratch directory of this
    !! process, shared by make_scratch and leave_build_scratch.
    character(len=*), intent(in) :: parent
    character(len=4096) :: path
    integer(c_int) :: pid
    character(len=32) :: pid_str

    pid = c_getpid()
    write (pid_str, '(i0)') pid
    path = trim(parent)//'/x3d2-memcheck-build.'//trim(pid_str)
  end function scratch_name

  subroutine skip_build(ctx, msg, cleanup, ok)
    !! Tail of every "Real build skipped" exit: print msg, remove the
    !! scratch directory if cleanup is set, and flag ok=.false.
    type(memcheck_ctx_t), intent(in) :: ctx
    character(len=*), intent(in) :: msg
    logical, intent(in) :: cleanup
    logical, intent(out) :: ok

    print '(a)', msg
    if (cleanup) call run_sh("rm -rf '"//trim(ctx%scratch_dir)//"'")
    ok = .false.
  end subroutine skip_build

end module m_memcheck_scratch
