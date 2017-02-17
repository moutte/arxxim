module M_System_Vars
!--
!-- variables for system definition
!--
  use M_Kinds
  !
  use M_T_System
  use M_T_Component,only: T_Component
  !
  implicit none
  !
  private
  !
  public:: System_Zero
  public:: System_Clean
  !
  real(dp),public:: TdgK, Pbar
  !
  character(len=7),public:: System_Type="AQUEOUS" !AQUEOUS, SIMPLEX, MOMAS, GLOBAL
  !
  logical,public:: BufferIsExtern= .false.  !.true.  !
  logical,public:: CpnIsSpc= .true.
  !
  type(T_Component),allocatable,public:: vCpn   (:) !-> default system
  type(T_Component),allocatable,public:: vCpnTot(:) !-> "global" system
  !
  ! stoikiometry table of the species in terms of the components
  real(dp),allocatable,public:: tStoikioCpn(:,:)
  !
  type(T_System):: SysDefault,SysTotal,SysInject,SysBox
  !
  !-- vars for GEM computations
  !-- ? should be in a specific M_GEM_Vars module ?
  !! type(T_Component),allocatable,public:: vCpnGEM(:)
  !! real(dp),         allocatable,public:: tStoikioGEM(:,:)
  
contains

subroutine System_Zero
  use M_Dtb_Const,only: Tref,Pref
  
  call System_Clean
  
  allocate(vCpn(0))
  allocate(vCpnTot(0))
  !allocate(vCpnInj(0))
  !allocate(vCpnBox(0))
  
  TdgK= Tref
  Pbar= Pref
  
end subroutine System_Zero

subroutine System_Clean
  if(allocated(vCpn))     deallocate(vCpn)
  if(allocated(vCpnTot))  deallocate(vCpnTot)
  if(allocated(tStoikioCpn)) deallocate(tStoikioCpn)
end subroutine System_Clean

end module M_System_Vars



