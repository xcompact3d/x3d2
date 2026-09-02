module m_checkpoint_state
  !! Interface for case- or model-specific state stored in solver checkpoints.
  use m_io_session, only: reader_session_t, writer_session_t

  implicit none

  type, abstract :: checkpoint_state_t
  contains
    procedure(write_checkpoint_iface), deferred :: write_checkpoint
    procedure(read_checkpoint_iface), deferred :: read_checkpoint
  end type checkpoint_state_t

  abstract interface
    subroutine write_checkpoint_iface(self, writer)
      import :: checkpoint_state_t, writer_session_t
      class(checkpoint_state_t), intent(inout) :: self
      type(writer_session_t), intent(inout) :: writer
    end subroutine write_checkpoint_iface

    subroutine read_checkpoint_iface(self, reader)
      import :: checkpoint_state_t, reader_session_t
      class(checkpoint_state_t), intent(inout) :: self
      type(reader_session_t), intent(inout) :: reader
    end subroutine read_checkpoint_iface
  end interface

end module m_checkpoint_state
