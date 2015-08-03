MODULE M_Dynam
  USE M_Kinds
  USE M_Trace,      ONLY: iDebug,Pause_,T_
  USE M_T_Component,ONLY: T_Component
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Dynam_Box
  PUBLIC:: Dynam_Initialize
  !
  ! PUBLIC:: TUnit, Tfinal, Dtime, Time
  ! PUBLIC:: PhiF, Fout, vQTotInj
  ! PUBLIC:: CoupledCoores
  
CONTAINS

SUBROUTINE Dynam_Box
  USE M_Dynam_Tools !Dynam_Zero_X,Dynam_Close
  USE M_Dynam_Read  !Dynam_Read
  USE M_Dynam_Tools !Dynam_Init
  USE M_Dynam_Files !Dynam_Files_Init,Dynam_Files_Close
  USE M_Dynam_Calc  !Dynam_Calc
  USE M_Equil_Write, ONLY: Equil_Write_Detail
  !
  USE M_System_Vars, ONLY: TdgK,Pbar,vCpn
  USE M_Dynam_Vars,  ONLY: DynTime, DynBox, DynBoxUser, Dynam_nStep
  !
  USE M_Dynam_Vars, ONLY: & !
  !-----------------------import from M_Dynam_Vars to publish in M_Dynam
  & TUnit, Tfinal, Dtime, & !
  & Time,                 & !
  & PhiF, Fout, vQTotInj, & !
  & CoupledCoores
  !---------------------------------------------------------------------
  LOGICAL::Ok
  !
  !---------------------------------------------------------Dynamic calc
  CALL Dynam_Zero_Numeric
  !
  CALL Dynam_Zero_Time(DynTime) !-> default values
  !
  CALL Dynam_Zero_Box(DynBox)
  !
  CALL Dynam_Read(DynTime,DynBox,DynBoxUser,Ok)
  !
  IF(.NOT. Ok) THEN
    ! write warning some where ... !!!
    RETURN !------------------------------------------------------RETURN
  ENDIF
  !
  CALL Dynam_Alloc
  !
  CALL Dynam_Init(TdgK,Pbar,DynTime,DynBox)
  !
  !----------------------------------write speciation of initial fluid--
  CALL Equil_Write_Detail("INI",TdgK,Pbar,vCpn) !----------------------/
  !
  CALL Dynam_ReStart_Time(DynTime)
  !
  CALL Dynam_ReStart_Box(DynBox,DynTime%TimeFactor)
  !
  CALL Dynam_Files_Init !(bCell=.FALSE.)
  !
  !---------------------------------------------------------Dynam Calc--
  CALL Dynam_Calc_Init
  CALL Dynam_Calc_Main
  CALL Dynam_Calc_Close
  !--------------------------------------------------------------------/
  CALL Dynam_Save
  !
  !------------------------------------write speciation of final fluid--
  CALL Equil_Write_Detail("DYN",TdgK,Pbar,vCpn) !----------------------/
  !
  CALL Dynam_Files_Close
  !
  CALL Dynam_Save_ToFile
  !
  CALL Dynam_Clean
  Dynam_nStep = 0
  !
  !--------------------------------------------------------/Dynamic calc
  !
  RETURN
ENDSUBROUTINE Dynam_Box

SUBROUTINE Dynam_Initialize
  USE M_Dynam_Tools
  USE M_System_Tools
  USE M_Basis,ONLY: Basis_Change
  USE M_Equil
  !
  USE M_Global_Vars,ONLY: vEle,vSpc
  USE M_System_Vars,ONLY: vCpn,vCpnTot
  USE M_Dynam_Vars, ONLY: vCpnInj,vCpnBox
  !
  LOGICAL:: OkInj,OkBox
  INTEGER:: iCp,J,N
  !
  CALL System_Build !-> build 'master' vCpn
  N= SIZE(vCpn)
  !
  IF(ALLOCATED(vCpnTot)) DEALLOCATE(vCpnTot)
  ALLOCATE(vCpnTot(1:N)); vCpnTot= vCpn
  !
  IF(ALLOCATED(vCpnInj)) DEALLOCATE(vCpnInj)
  ALLOCATE(vCpnInj(1:N))
  CALL Dynam_Init_FluidInject(OkInj)
  !
  IF(ALLOCATED(vCpnBox)) DEALLOCATE(vCpnBox)
  ALLOCATE(vCpnBox(1:N))
  CALL Dynam_Init_FluidBox(OkBox)
  !
  IF(OkInj) THEN
    !permute vCpnInj -> its element order must be same as vCpnBox
    vCpn= vCpnInj !vCpn USEd as permutation buffer
    DO iCp=1,N
      DO J=1,N
        IF(vCpn(J)%iEle==vCpnBox(iCp)%iEle) EXIT
      ENDDO
      vCpnInj(iCp)= vCpn(J)
    ENDDO
    vCpn= vCpnBox !-> vCpnBox will be the 'default' vCpn
  ELSE
    vCpnInj= vCpnBox
  ENDIF
  !
  ! print *, "Dynam_Initialize"
  ! print '(20G15.7)',(vCpnInj(iCp)%Mole,iCp=1,size(vCpnInj))   ;   pause
  !when no SYSTEM.INJECT is given, vCpnInj is the vCpn READ from SYSTEM
  !-> as default,
  !  fluid in box is main system,
  !  fluid injected is also fluid in box
  !
  IF(iDebug>2) THEN
    PRINT '(A)',"vCpnBox,vCpnInj"
    DO iCp=1,N
      PRINT '(2(2A,G15.6,1X))', &
      & "CpnBox ", vEle(vCpnBox(iCp)%iEle)%NamEl,vCpnBox(iCp)%Mole, &
      & "CpnInj ", vEle(vCpnInj(iCp)%iEle)%NamEl,vCpnInj(iCp)%Mole
    ENDDO
    CALL Pause_
  ENDIF
  !
ENDSUBROUTINE Dynam_Initialize

SUBROUTINE Dynam_Init_FluidInject(Ok)
!--
!-- read constraints on injected fluid,
!-- then calc' speciation, then make all cpn inert,
!-- CAVEAT:
!-- Dynam_Init_FluidInject must be called BEFORE Dynam_Init_FluidBox,
!-- otherwise the speciation of fluid in box
!.  (calc'ted in Dynam_Init_FluidBox and used in dynamic) 
!-- would be overwritten by the call to Equil_Calc in Dynam_InitFluidInject
!--
  USE M_Basis,       ONLY: Basis_Change
  USE M_System_Tools,ONLY: System_Build_Custom
  USE M_Equil
  !
  USE M_Global_Vars, ONLY: vEle,vSpc,vMixFas
  USE M_System_Vars, ONLY: TdgK,Pbar,vCpn,vCpnTot
  USE M_Dynam_Vars,  ONLY: vCpnInj
  !
  LOGICAL,INTENT(OUT):: Ok
  
  integer:: i
  
  CALL System_Build_Custom( &
  & "SYSTEM.INJECT", &
  & vEle,vSpc,vMixFas,vCpnTot, &
  & TdgK,Pbar, &
  & vCpnInj,Ok) !,SysType)
  !
  IF(Ok) THEN
    vCpn= vCpnInj
    !
    CALL Equil_Calc("INJ") != calc' speciation for inject system
    !
    CALL Basis_Change("DYN",vCpn) != make all components inert
    !
    vCpnInj= vCpn
  ENDIF
  
  ! print *,"Dynam_Init_FluidInject"
  ! print '(20G15.7)',(vCpn(i)%Mole,i=1,size(vCpn))   ;   pause
  
  RETURN
ENDSUBROUTINE Dynam_Init_FluidInject

SUBROUTINE Dynam_Init_FluidBox(Ok)
  USE M_System_Tools,ONLY: System_Build_Custom
  USE M_Basis,       ONLY: Basis_Change
  USE M_Equil
  !
  USE M_Global_Vars, ONLY: vEle,vSpc,vMixFas
  USE M_System_Vars, ONLY: TdgK,Pbar,vCpn,vCpnTot
  USE M_Dynam_Vars,  ONLY: vCpnBox
  !
  LOGICAL,INTENT(OUT):: Ok
  
  CALL System_Build_Custom( &
  & "SYSTEM.BOX", &
  & vEle,vSpc,vMixFas,vCpnTot, &
  & TdgK,Pbar, &
  & vCpnBox,Ok) !,SysType)
  
  IF(Ok) THEN  ;  vCpn= vCpnBox
  ELSE         ;  vCpn= vCpnTot !work with current vCpn !!!
  ENDIF
  
  CALL Equil_Calc("BOX") != calc' speciation for fluid in box
  CALL Basis_Change("DYN",vCpn) != make all components inert
  
  vCpnBox= vCpn
  
  RETURN
ENDSUBROUTINE Dynam_Init_FluidBox

ENDMODULE M_Dynam

