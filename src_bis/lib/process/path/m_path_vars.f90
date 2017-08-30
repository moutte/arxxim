module M_Path_Vars !variables for path calculations
  use M_Kinds
  use M_T_Tpcond,only: T_TPCond
  !
  implicit none
  !
  public
  !
  logical, allocatable:: vLPath(:),vFasFound(:)
  integer, allocatable:: vPhasBegin(:),vPhasFinal(:) !index of phases in vMixFas
  real(dp),allocatable:: tPathData(:,:)
  real(dp),allocatable:: vPathLogK(:)
  !
  integer:: DimPath
  integer:: iLogK
  !
  integer,parameter:: TotalMixStep= 100
  !
  real(dp),allocatable:: tGrt(:,:)
  real(dp),allocatable:: tPathResults(:,:)
  !
  !character(len=15):: PathMode
  !
  type(T_TPCond),allocatable:: vTPpath(:)
  !
contains

subroutine Path_Vars_Clean
  if(allocated(vLPath))      deallocate(vLPath)
  if(allocated(vFasFound))   deallocate(vFasFound)
  if(allocated(vPhasBegin))  deallocate(vPhasBegin)
  if(allocated(vPhasFinal))  deallocate(vPhasFinal)
  if(allocated(tPathData))   deallocate(tPathData)
  if(allocated(vTPpath))     deallocate(vTPpath)
  if(allocated(vPathLogK))   deallocate(vPathLogK)
end subroutine Path_Vars_Clean 

end module M_Path_Vars

