module M_Dynam
  use M_Kinds
  use M_Trace,      only: iDebug,Pause_,T_
  use M_T_Component,only: T_Component
  !
  implicit none
  !
  private
  !
  public:: Dynam_Box
  public:: Dynam_Initialize
  !
  ! public:: TUnit, Tfinal, Dtime, Time
  ! public:: PhiF, Fout, vQTotInj
  ! public:: CoupledCoores
  
contains

subroutine Dynam_Box
  use M_Dynam_Tools !Dynam_Zero_X,Dynam_Close
  use M_Dynam_Read  !Dynam_Read
  use M_Dynam_Tools !Dynam_Init
  use M_Dynam_Files !Dynam_Files_Init,Dynam_Files_Close
  use M_Dynam_Calc  !Dynam_Calc
  use M_Equil_Write, only: Equil_Write_Detail
  !
  use M_System_Vars, only: TdgK,Pbar,vCpn
  use M_Dynam_Vars,  only: DynTime, DynBox, DynBoxUser, Dynam_nStep
  !
  use M_Dynam_Vars, only: & !
  !-----------------------import from M_Dynam_Vars to publish in M_Dynam
  & TUnit, Tfinal, Dtime, & !
  & Time,                 & !
  & PhiF, Fout, vQTotInj, & !
  & CoupledCoores
  !---------------------------------------------------------------------
  logical::Ok
  !
  !---------------------------------------------------------Dynamic calc
  call Dynam_Zero_Numeric
  !
  call Dynam_Zero_Time(DynTime) !-> default values
  !
  call Dynam_Zero_Box(DynBox)
  !
  call Dynam_Read(DynTime,DynBox,DynBoxUser,Ok)
  !
  if(.not. Ok) then
    ! write warning some where ... !!!
    return !------------------------------------------------------return
  end if
  !
  call Dynam_Alloc
  !
  call Dynam_Init(TdgK,Pbar,DynTime,DynBox)
  !
  !----------------------------------write speciation of initial fluid--
  call Equil_Write_Detail("INI",TdgK,Pbar,vCpn) !----------------------/
  !
  call Dynam_ReStart_Time(DynTime)
  !
  call Dynam_ReStart_Box(DynBox,DynTime%TimeFactor)
  !
  call Dynam_Files_Init !(bCell=.false.)
  !
  !---------------------------------------------------------Dynam Calc--
  call Dynam_Calc_Init
  call Dynam_Calc_Main
  call Dynam_Calc_Close
  !--------------------------------------------------------------------/
  call Dynam_Save
  !
  !------------------------------------write speciation of final fluid--
  call Equil_Write_Detail("DYN",TdgK,Pbar,vCpn) !----------------------/
  !
  call Dynam_Files_Close
  !
  call Dynam_Save_ToFile
  !
  call Dynam_Clean
  Dynam_nStep = 0
  !
  !--------------------------------------------------------/Dynamic calc
  !
  return
end subroutine Dynam_Box

subroutine Dynam_Initialize
  use M_Dynam_Tools
  use M_System_Tools
  use M_Basis,only: Basis_Change
  use M_Equil
  !
  use M_Global_Vars,only: vEle,vSpc
  use M_System_Vars,only: vCpn,vCpnTot
  use M_Dynam_Vars, only: vCpnInj,vCpnBox
  !
  logical:: OkInj,OkBox
  integer:: iCp,J,N
  !
  call System_Build !-> build 'master' vCpn
  N= size(vCpn)
  !
  if(allocated(vCpnTot)) deallocate(vCpnTot)
  allocate(vCpnTot(1:N)); vCpnTot= vCpn
  !
  if(allocated(vCpnInj)) deallocate(vCpnInj)
  allocate(vCpnInj(1:N))
  call Dynam_Init_FluidInject(OkInj)
  !
  if(allocated(vCpnBox)) deallocate(vCpnBox)
  allocate(vCpnBox(1:N))
  call Dynam_Init_FluidBox(OkBox)
  !
  if(OkInj) then
    !permute vCpnInj -> its element order must be same as vCpnBox
    vCpn= vCpnInj !vCpn used as permutation buffer
    do iCp=1,N
      do J=1,N
        if(vCpn(J)%iEle==vCpnBox(iCp)%iEle) exit
      end do
      vCpnInj(iCp)= vCpn(J)
    end do
    vCpn= vCpnBox !-> vCpnBox will be the 'default' vCpn
  else
    vCpnInj= vCpnBox
  end if
  !
  ! print *, "Dynam_Initialize"
  ! print '(20G15.7)',(vCpnInj(iCp)%Mole,iCp=1,size(vCpnInj))   ;   pause
  !when no SYSTEM.INJECT is given, vCpnInj is the vCpn read from SYSTEM
  !-> as default,
  !  fluid in box is main system,
  !  fluid injected is also fluid in box
  !
  if(iDebug>2) then
    print '(A)',"vCpnBox,vCpnInj"
    do iCp=1,N
      print '(2(2A,G15.6,1X))', &
      & "CpnBox ", vEle(vCpnBox(iCp)%iEle)%NamEl,vCpnBox(iCp)%Mole, &
      & "CpnInj ", vEle(vCpnInj(iCp)%iEle)%NamEl,vCpnInj(iCp)%Mole
    end do
    call Pause_
  end if
  !
end subroutine Dynam_Initialize

subroutine Dynam_Init_FluidInject(Ok)
!--
!-- read constraints on injected fluid,
!-- then calc' speciation, then make all cpn inert,
!-- CAVEAT:
!-- Dynam_Init_FluidInject must be called BEFORE Dynam_Init_FluidBox,
!-- otherwise the speciation of fluid in box
!.  (calc'ted in Dynam_Init_FluidBox and used in dynamic) 
!-- would be overwritten by the call to Equil_Calc in Dynam_InitFluidInject
!--
  use M_Basis,       only: Basis_Change
  use M_System_Tools,only: System_Build_Custom
  use M_Equil
  !
  use M_Global_Vars, only: vEle,vSpc,vMixFas
  use M_System_Vars, only: TdgK,Pbar,vCpn,vCpnTot
  use M_Dynam_Vars,  only: vCpnInj
  !
  logical,intent(out):: Ok
  
  integer:: i
  
  call System_Build_Custom( &
  & "SYSTEM.INJECT", &
  & vEle,vSpc,vMixFas,vCpnTot, &
  & TdgK,Pbar, &
  & vCpnInj,Ok) !,SysType)
  !
  if(Ok) then
    vCpn= vCpnInj
    !
    call Equil_Calc("INJ") != calc' speciation for inject system
    !
    call Basis_Change("DYN",vCpn) != make all components inert
    !
    vCpnInj= vCpn
  end if
  
  ! print *,"Dynam_Init_FluidInject"
  ! print '(20G15.7)',(vCpn(i)%Mole,i=1,size(vCpn))   ;   pause
  
  return
end subroutine Dynam_Init_FluidInject

subroutine Dynam_Init_FluidBox(Ok)
  use M_System_Tools,only: System_Build_Custom
  use M_Basis,       only: Basis_Change
  use M_Equil
  !
  use M_Global_Vars, only: vEle,vSpc,vMixFas
  use M_System_Vars, only: TdgK,Pbar,vCpn,vCpnTot
  use M_Dynam_Vars,  only: vCpnBox
  !
  logical,intent(out):: Ok
  
  call System_Build_Custom( &
  & "SYSTEM.BOX", &
  & vEle,vSpc,vMixFas,vCpnTot, &
  & TdgK,Pbar, &
  & vCpnBox,Ok) !,SysType)
  
  if(Ok) then  ;  vCpn= vCpnBox
  else         ;  vCpn= vCpnTot !work with current vCpn !!!
  end if
  
  call Equil_Calc("BOX") != calc' speciation for fluid in box
  call Basis_Change("DYN",vCpn) != make all components inert
  
  vCpnBox= vCpn
  
  return
end subroutine Dynam_Init_FluidBox

end module M_Dynam

