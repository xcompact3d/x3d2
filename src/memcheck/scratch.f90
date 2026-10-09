module m_memcheck_scratch
  !! Scratch directory handling of the x3d2-memcheck --build path.
  use iso_c_binding, only: c_char, c_int, c_size_t, c_ptr, c_null_char
  use m_memcheck_context, only: memcheck_ctx_t
  implicit none

  private
  public :: ibm_mask_filename, validate_build_inputs, make_scratch, &
            leave_build_scratch

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

  pure function ibm_mask_filename(px, py, pz) result(fname)
    !! Build the ibm_<BC-suffix>.bp mask filename from the three periodic_BC
    !! flags, e.g. mesh%grid%periodic_BC(1:3) in build_and_measure or the
    !! module periodic_x/y/z classify_bc sets - both carry the same
    !! information since both derive from domain_cfg%BC_x/y/z via
    !! periodic_dir. Shared by build_and_measure's ibm_on pre-check and
    !! make_scratch's mask mirroring so the suffix logic can never
    !! diverge between the two (src/module/ibm.f90:70-75 is the third,
    !! authoritative copy this mirrors).
    logical, intent(in) :: px, py, pz
    character(len=16) :: fname
    character(len=3) :: bc_suffix

    bc_suffix(1:1) = '0'
    if (.not. px) bc_suffix(1:1) = '1'
    bc_suffix(2:2) = '0'
    if (.not. py) bc_suffix(2:2) = '1'
    bc_suffix(3:3) = '0'
    if (.not. pz) bc_suffix(3:3) = '1'
    fname = "ibm_"//bc_suffix//".bp"
  end function ibm_mask_filename

  subroutine validate_build_inputs(ctx, ok)
    !! Validates the known failure causes of a real build up front, before
    !! touching the filesystem: the input file must exist,
    !! domain_cfg%flow_case_name must be one this tool (and
    !! build_and_measure's own select case) can dispatch, and (ibm_on=T) the
    !! matching ibm_<BC-suffix>.bp mask file must already be present in the
    !! invoking directory. Records the invoking directory in orig_dir for
    !! make_scratch. Sets ok=.false. (and the estimate above stands, like
    !! the other run_tier3 skips) on any failure.
    type(memcheck_ctx_t), intent(inout) :: ctx
    logical, intent(out) :: ok

    character(len=16) :: ibm_file
    logical :: ibm_file_exists, input_exists, flow_case_supported

    ok = .true.
    ctx%orig_dir = current_dir()

    inquire (file=trim(ctx%input_path), exist=input_exists)
    if (.not. input_exists) then
      call skip_build(ctx, 'Real build skipped: input file not found: '// &
                      trim(ctx%input_path), .false., ok)
      return
    end if

    select case (trim(ctx%domain_cfg%flow_case_name))
    case ('tgv', 'generic', 'channel', 'cylinder')
      flow_case_supported = .true.
    case default
      flow_case_supported = .false.
    end select
    if (.not. flow_case_supported) then
      call skip_build(ctx, "Real build skipped: flow case '"// &
                      trim(ctx%domain_cfg%flow_case_name)// &
                      "' has no dispatch in x3d2-memcheck", .false., ok)
      return
    end if

    if (ctx%solver_cfg%ibm_on) then
      ibm_file = ibm_mask_filename(ctx%periodic_x, ctx%periodic_y, &
                                   ctx%periodic_z)
      inquire (file=trim(ctx%orig_dir)//'/'//trim(ibm_file), &
               exist=ibm_file_exists)
      if (.not. ibm_file_exists) then
        call skip_build(ctx, 'Real build skipped: ibm_on=T but the matching &
                        &ibm_<BC-suffix>.bp mask file was not found in the &
                        &working directory; the estimate above stands.', &
                        .false., ok)
        return
      end if
    end if
  end subroutine validate_build_inputs

  subroutine make_scratch(ctx, ok)
    !! Real build runs inside a throwaway x3d2-memcheck-build.<pid>
    !! subdirectory of the invoking directory, because a case build
    !! initialises monitoring (writes monitoring.csv) and calls
    !! postprocess(0), which clobbered run directories on 2026-09-15. A
    !! relative input path is mirrored into the scratch directory by
    !! symlink; a relative path containing a '..' component (or a single
    !! quote, which the shell quoting below cannot handle) cannot be
    !! mirrored this way and is rejected - pass an absolute path instead.
    !! Sets ok=.false. (and the estimate above stands, like the other
    !! run_tier3 skips) on any failure; leave_build_scratch removes the
    !! scratch directory once the build finishes. validate_build_inputs
    !! must have run first. If an unexpected error stop happens after the
    !! chdir below regardless - a failure mode the pre-validation does not
    !! cover - the scratch directory is left behind under the invoking
    !! directory; the next run with the same pid in the same directory
    !! removes it as a stale leftover before creating its own.
    type(memcheck_ctx_t), intent(inout) :: ctx
    logical, intent(out) :: ok

    integer(c_int) :: rc
    integer :: st, slash_pos
    character(len=16) :: ibm_file
    logical :: input_exists

    ok = .true.
    ctx%scratch_dir = scratch_name(ctx%orig_dir)

    ! Stale leftover (see the docstring) would make a plain mkdir fail.
    call run_sh("test -d '"//trim(ctx%scratch_dir)//"'", st)
    if (st == 0) then
      call run_sh("rm -rf '"//trim(ctx%scratch_dir)//"'", st)
      if (st /= 0) then
        call skip_build(ctx, 'Real build skipped: could not remove stale &
                        &scratch directory '//trim(ctx%scratch_dir), &
                        .false., ok)
        return
      end if
      print '(a,a)', 'Removed stale scratch directory ', trim(ctx%scratch_dir)
    end if

    call run_sh("mkdir '"//trim(ctx%scratch_dir)//"'", st)
    if (st /= 0) then
      call skip_build(ctx, 'Real build skipped: could not create scratch &
                      &directory '//trim(ctx%scratch_dir), .false., ok)
      return
    end if

    if (ctx%input_path(1:1) /= '/') then
      if (has_dotdot_component(trim(ctx%input_path)) .or. &
          index(trim(ctx%input_path), "'") > 0) then
        call skip_build(ctx, "Real build skipped: relative input path with &
                        &'..' cannot be mirrored; pass an absolute path", &
                        .true., ok)
        return
      end if
      slash_pos = index(trim(ctx%input_path), '/', back=.true.)
      if (slash_pos > 0) then
        call run_sh("mkdir -p '"//trim(ctx%scratch_dir)//'/'// &
                    trim(ctx%input_path(1:slash_pos - 1))//"'", st)
        if (st /= 0) then
          call skip_build(ctx, 'Real build skipped: could not prepare scratch &
                          &directory (mkdir of the input''s parent failed)', &
                          .true., ok)
          return
        end if
      end if
      call run_sh("ln -s '"//trim(ctx%orig_dir)//'/'// &
                  trim(ctx%input_path)//"' '"//trim(ctx%scratch_dir)// &
                  '/'//trim(ctx%input_path)//"'", st)
      if (st /= 0) then
        call skip_build(ctx, 'Real build skipped: could not prepare scratch &
                        &directory (input symlink failed)', .true., ok)
        return
      end if
      inquire (file=trim(ctx%scratch_dir)//'/'//trim(ctx%input_path), &
               exist=input_exists)
      if (.not. input_exists) then
        call skip_build(ctx, 'Real build skipped: could not prepare scratch &
                        &directory (input symlink failed)', .true., ok)
        return
      end if
    end if

    if (ctx%solver_cfg%ibm_on) then
      ibm_file = ibm_mask_filename(ctx%periodic_x, ctx%periodic_y, &
                                   ctx%periodic_z)
      call run_sh("ln -s '"//trim(ctx%orig_dir)//'/'// &
                  trim(ibm_file)//"' '"//trim(ctx%scratch_dir)// &
                  '/'//trim(ibm_file)//"'", st)
      if (st /= 0) then
        call skip_build(ctx, 'Real build skipped: could not prepare scratch &
                        &directory (ibm mask symlink failed)', .true., ok)
        return
      end if
    end if

    rc = c_chdir(trim(ctx%scratch_dir)//c_null_char)
    if (rc /= 0) then
      call skip_build(ctx, 'Real build skipped: could not chdir into scratch &
                      &directory '//trim(ctx%scratch_dir), .true., ok)
      return
    end if

    print '(a,a,a)', 'Real build scratch directory: ', trim(ctx%scratch_dir), &
      ' (removed after the build)'
  end subroutine make_scratch

  subroutine leave_build_scratch(ctx)
    !! Restore the invoking directory and remove the scratch directory
    !! make_scratch created. Called after build_and_measure returns,
    !! on both the normal path and the ibm_missing early return.
    type(memcheck_ctx_t), intent(inout) :: ctx
    integer(c_int) :: rc
    character(len=4096) :: expected_scratch_dir

    rc = c_chdir(trim(ctx%orig_dir)//c_null_char)
    if (rc /= 0) &
      error stop 'x3d2-memcheck: could not chdir back to the invoking &
        &directory after the real build; state is unknown, not removing &
        &the scratch directory.'

    ! Never remove anything other than the exact scratch directory
    ! make_scratch created and chdir'd into.
    expected_scratch_dir = scratch_name(ctx%orig_dir)
    if (trim(ctx%scratch_dir) == trim(expected_scratch_dir)) &
      call run_sh("rm -rf '"//trim(ctx%scratch_dir)//"'")

    ctx%orig_dir = ''
    ctx%scratch_dir = ''
  end subroutine leave_build_scratch

end module m_memcheck_scratch
