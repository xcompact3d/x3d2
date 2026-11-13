module m_cuda_errors
  !! CUDA and cuFFT counterparts of m_logger's check_mpi: turn a failed
  !! return code into a report naming the operation, the rank and the error,
  !! followed by an abort of the whole run.
  !!
  !! A CUDA failure is local to the rank that sees it, while its partners
  !! carry on into MPI and would wait for it forever, so these always abort
  !! the job rather than stopping the one rank.
  use cudafor, only: cudaSuccess, cudaGetErrorString

  use m_logger, only: abort_run, str

  implicit none

  private
  public :: check_cuda, check_cufft

contains

  subroutine check_cuda(ierr, operation)
    !! Abort the run if a CUDA runtime call did not return cudaSuccess.
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: operation

    if (ierr == cudaSuccess) return

    call abort_run(operation//' failed: '//trim(cudaGetErrorString(ierr)) &
                   //' (CUDA error code '//str(ierr)//')')
  end subroutine check_cuda

  subroutine check_cufft(ierr, operation)
    !! Abort the run if a cuFFT call did not return CUFFT_SUCCESS (0).
    !! cuFFT has no error-string function, so only the code is reported.
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: operation

    if (ierr == 0) return

    call abort_run(operation//' failed (cuFFT error code '//str(ierr)//')')
  end subroutine check_cufft

end module m_cuda_errors
