module M_Equil_Specia
  use M_Kinds
  use M_Trace,only: iDebug,fTrc
  implicit none
  !
  private
  !
  public:: Equil_Specia
  !
contains

subroutine Equil_Specia(iErr)
!--
!-- calculate equilibrium speciation of the aqueous phase
!-- (i.e. no consideration of saturation with minerals or gases)
!--
  use M_Equil_Solve,only: Equil_Solve
  !
  integer,intent(out):: iErr
  !
  integer :: nIts
  !! real(dp):: CpuBegin,CpuEnd
  
  if(idebug>1) write(fTrc,'(/,A)') "< Equil_Specia >"
  !! if(idebug>1) call CPU_TIME(CpuBegin)
  !
  call Equil_Solve(nIts,iErr)
  !
  !! if(idebug>1) call CPU_TIME(CpuEnd)
  !! if(idebug>1) write(fTrc,'(/,A,G15.6,/)') "CPU-Time",1.D3*(CpuEnd-CpuBegin)
  
  if(idebug>1) write(fTrc,'(A,/)') "</ Equil_Specia >"
  
  return
end subroutine Equil_Specia

end module M_Equil_Specia


