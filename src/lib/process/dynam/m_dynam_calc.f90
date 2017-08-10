module M_Dynam_Calc
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_
  use M_Numeric_Const, only: Ln10
  use M_IOTools
  use M_Trace, only: DebugCoores
  !
  !! use M_Stockvar_Kinxim,only: LSTOCK,SET_STOCKVAR
  !
  use M_Global_Vars,only: vSpc,vKinFas
  use M_System_Vars,only: TdgK,Pbar
  use M_Basis_Vars, only: nCi,nAs,vPrmBk
  use M_Basis_Vars, only: tAlfPr,tAlfAs
  !
  use M_Dynam_Files,only: Dynam_Files_OpenLogs,Dynam_Files_WriteFasAff
  use M_Dynam_Vars, only: Time,TFinal,dTime,dTmin,dTMax,TUnit,dTSav,TimeFactor
  use M_Dynam_Vars, only: TdgK0,Pbar0
  use M_Dynam_Vars, only: UDarcy,VBox,dX,RhoF
  use M_Dynam_Vars, only: FOut,PhiF0,PhiF,PhiInert
  use M_Dynam_Vars, only: vMolF,vTotF,vTotInj
  use M_Dynam_Vars, only: vLnAct,vLnBuf,vLnGam,pH_
  use M_Dynam_Vars, only: vMolarVol,vKinQsk,vStatusK
  use M_Dynam_Vars, only: tAlfKin,vLKinActiv,vLEquActiv,vKinPrm
  use M_Dynam_Vars, only: vMolK,vMolK0,vSurfK,vSurfK0
  use M_Dynam_Vars, only: vVmQsK,vVmAct,UpdateMassFluid
  use M_Dynam_Vars, only: &
  & Dynam_nStep,     &
  & Dynam_nNewtIter, &
  & Dynam_nTotalNewtIter, &
  & nEvalFunc
  use M_Dynam_Vars, only: Extrapole
  use M_Dynam_Vars, only: SteadyState_Stop
  !
  use M_Dynam_Tools,only: Dynam_TP_Update
  use M_Dynam_OneStep_Tools
  use M_Dynam_OneStep_Solve
  !
  use M_Dynam_Solve_Vars,only: Dynam_Solve_Vars_Alloc, Dynam_Solve_Vars_Clean
  !
  implicit none

  private

  public:: Dynam_Calc
  public:: Dynam_Calc_Init
  public:: Dynam_Calc_Main
  public:: Dynam_Calc_Close

  !---------------------------------------------------------------vars--
  real(dp),dimension(:),allocatable:: & !
  & vMolF_0, & !species mole nrs in fluid, prev.step -> to compute variations
  & vMolK_0, & !mole nrs of minerals, prev.step -> to compute variations
  & vTotF_0, & !component mole nrs in fluid, prev.step -> to compute variations
  & vTotK,   & !component mole nrs in minerals
  & vTotK_0, & !component mole nrs in minerals, prev.step -> to compute variations
  & vMolM_Slope
  !
  logical,public:: DTimeTooSmall,PorosTooSmall
  !
  real(dp):: PhiF_0,X !,UDarcyOut
  real(dp):: VarMax,dTimeOut
  real(dp):: MaxAqu=Zero, MaxMin=Zero
  !
  real(dp):: ZeroVarMax= 1.0D-15
  !
  integer :: nCp,nMk,nAq
  ! integer :: nMkA0,nEqA0
  integer :: f1,f2,f3,f4
  integer :: iStep,Newt_iDo,NZeroVar,Newt_iErr
  integer :: NZeroMin= 12
  !
  real(dp):: TSave
  real(dp):: TAgain
  ! TAgain= duration of an additional simulation starting at time TFinal
  !
  real(dp):: CpuBegin,CpuEnd
  !--------------------------------------------------------------------/
  
contains

subroutine Dynam_Calc
  call Dynam_Calc_Init
  call Dynam_Calc_Main
  call Dynam_Calc_Close
end subroutine Dynam_Calc

subroutine Dynam_Calc_Init
  !
  if(idebug>1) write(fTrc,'(/,A)') "< CalcDynamic"
  !
  SteadyState_Stop= (iDebug>2)
  if(DebugCoores) SteadyState_Stop= .false.
  !
  nCp= size(vTotF)
  nMk= size(vKinFas)
  nAq= count(vSpc%Typ=="AQU")
  !
  call Dynam_Solve_Vars_Alloc(nMk,nCi)
  !
  allocate(vMolF_0(1:nAq))
  allocate(vMolK_0(1:nMk))
  allocate(vTotF_0(1:nCp))
  !
  allocate(vTotK(1:nCp))
  allocate(vTotK_0(1:nCp))
  !
  allocate(vMolM_Slope(1:nMk))
  vMolM_Slope(:)= Zero
  !
  f1= 0  ;  f2= 0
  f3= 0  ;  f4= 0
  if(iDebug>2) call Dynam_Files_OpenLogs(f1,f2,f3,f4)
  !
  TAgain= TFinal
  TSave=Zero !------------------for output at definite time span dTSav--
  !-------------------------------(provisional, needs improvement ...)--
  !
  DTimeTooSmall= .false.
  PorosTooSmall= .false.
  !
  NZeroVar= 0
  Newt_iDo= 0
  VarMax=   Zero
  !
  !! !---------- message output for Coores --------------
  !! if (DebugCoores) &
  !! & print '(1A6, 1I6, 1A, 1F14.6, 1A)', &
  !! & '[ STEP',iStep,'] - [T = ',time,', ------ START    ]'
  !! !--- store variable and output flux for Coores -----
  !! if (LSTOCK) &
  !! & call set_stockvar( vMolF (vPrmBk(1:nAq)) , Fout/VBox, Time, dTime, VarMax , 'D' )
  !! !---------------------------------------------------
  !
  !------------------------------------------initialize kinetic phases--
  X= maxval(vMolF(1:nAq))
  call Dynam_OneStep_KinFas_Update( & !
  & log(vMolF(1:nAq)/X),   & !in
  & vMolK,                 & !in
  & vKinQsk,vStatusK,      & !out
  & vVmQsK,vVmAct,         & !out
  & vLKinActiv,vLEquActiv, & !out
  & vKinPrm,PhiInert)        !out
  !
  call Dynam_OneStep_KinFasSurf_Update( & !
  & vLKinActiv,  & ! IN
  & vStatusK,    & ! IN
  & vMolK0,vMolK,& ! IN
  & PhiF0,PhiF,  & ! IN
  & vSurfK0,     & ! IN
  & vSurfK)        ! OUT
  !
  !! nMkA0= count(vLKinActiv)
  !! nEqA0= count(vLEquActiv)
  !! !! !
  !! if(LogForAqu) then
  !!   call Dynam_OneStep_Alloc(nAq +nMkA0 +nEqA0)
  !! else
  !!   call Dynam_OneStep_Alloc(nAq +nMkA0 +nEqA0)
  !! end if
  !
  !----------------------------------------/initialize kinetic phases--
  !
  call Dynam_OneStep_TotK( & !
  & vMolK, & !in
  & vTotK)   !out
  !
  call Dynam_OneStep_TotF( & !
  & vMolF, & !in
  & vTotF)   !out
  !
  if(idebug>1) write(fTrc,'(/,A)') "</CalcDynamic"
  !
end subroutine Dynam_Calc_Init

subroutine Dynam_Calc_Main
  use M_Dynam_Vars,only: NewtTolF,NewtTolF_Equil,NewtTolX,NewtTolMin
  use M_Dynam_Vars,only: NewtMaxIts,NewtIterMax,NewtIterMin
  use M_Dynam_Vars,only: Time_Decrease,Time_Increase,iCtrlTime
  use M_Dynam_Vars,only: fSavRate,fSavTime,DebNewt
  !----------------------------------------------for test basis change--
  use M_System_Vars,only: vCpn
  use M_Basis,only: Basis_Change_Wrk
  !--------------------------------------------------------------------/
  
  if(iDebug>2) call CPU_TIME(CpuBegin)
  if(iDebug>2) nEvalFunc= 0
  !
  if(iDebug>1) print '(/,A,/)', &
  & "OUTPUT= iStep,Time,dTime,PhiF,nEqA,nMkA,iDo2,Newt_iDo,Newt_iErr, &
  & iMaxDelta,VarMax,NewtErrF,NewtErrX"
  !
  iStep= 0
  !-------------------------------------------------------------DoStep--
  DoStep: do !while (Time<tFinal)
    !
    iStep= iStep+1
    DebNewt= (iStep<4) !write jacobian for the first 4 steps ...
    !if(Time<TFinal) call CtrlDtime(iCtrlTime,Newt_iDo,MaxAqu,MaxMin,dTime)
    !
    !--------adjust mole numbers to available volume, based on density--
    if(UpdateMassFluid) call Dynam_OneStep_Adjust( & !
    & RhoF,VBox,PhiF,  & !in
    & vMolF,           & !inout
    & vTotF)             !inout
    !
    !----------------save current value for calculations of variations--
    vMolF_0(:)= vMolF(:)
    vTotF_0(:)= vTotF(:)
    PhiF_0=     PhiF
    vMolK_0(:)= vMolK(:)
    vTotK_0(:)= vTotK(:)
    !
    if(TP_Changed(TdgK,Pbar)) call Dynam_TP_Update(TdgK,Pbar,TimeFactor)
    !
    !---------------------------------------------------solve one step--
    call Dynam_OneStep_Solve(       & !
    & F2,iStep,                     & !IN
    & vLKinActiv,vLEquActiv,        & !IN
    & vMolM_Slope,                  & !IN
    & NewtTolF,NewtTolF_Equil,      & !IN
    & NewtTolX,NewtTolMin,          & !IN
    & Time_Decrease,NewtMaxIts,     & !IN
    & Time,dTmin,dTime,             & !IN
    & dTimeOut,Newt_iDo,Newt_iErr,  & !OUT
    & DTimeTooSmall,VarMax)           !OUT
    !------------------------------------------------------------------/
    !
    dTime= dTimeOut
    !
    Time= Time +dTime !------------------------------update time elapsed
    !
    !-------------------------------------------------------basis change
    if(iDebug>2) call Basis_Change_Wrk(vCpn)
    !------------------------------------------------------/basis change
    !
    if(iDebug>2) &
    & call Dynam_OneStep_WriteElements(iStep,Time,dTime,vMolF,vMolK)
    !
    call Dynam_OneStep_TotF( & !
    & vMolF, & !in
    & vTotF)   !out
    !
    call Dynam_OneStep_TotK( & !
    & vMolK, & !in
    & vTotK)   !out
    !
    PhiF= One - sum(vMolK(:)*vMolarVol(:)) /Vbox
    !
    if(Extrapole) &
    & vMolM_Slope(:)= (vMolK(:) -vMolK_0(:)) /dTime
    !vMolM_Slope(:)= Zero
    !
    if(iStep>3 .and. F1>0) &
    & call Dynam_OneStep_Write_Bilan( & !
    & F1,iStep,Time,dTimeOut,               &  !IN
    & PhiF,UDarcy,VBox,dX,                  &  !IN
    & vTotInj,vTotF_0,vTotF,vTotK_0,vTotK)     !IN
    !
    !----------------------------------------------update kin'phase data
    X= maxval(vMolF(1:nAq))
    !
    call Dynam_OneStep_KinFas_Update( & !
    & log(vMolF(1:nAq)/X),   & !in
    & vMolK,                 & !in
    & vKinQsk,vStatusK,      & !out
    & vVmQsK,vVmAct,         & !out
    & vLKinActiv,vLEquActiv, & !out
    & vKinPrm,PhiInert)        !out
    !
    call Dynam_OneStep_KinFasSurf_Update( & !
    & vLKinActiv,  & ! IN
    & vStatusK,    & ! IN
    & vMolK0,vMolK,& ! IN
    & PhiF0,PhiF,  & ! IN
    & vSurfK0,     & ! IN
    & vSurfK)        ! OUT
    !
    !--------------------------------------------------change dimensions
    !! if(OnlyActiv) then
    !!   if(count(vLKinActiv)/=nMkA0 .or. count(vLEquActiv)/=nEqA0) then
    !!     !
    !!     nMkA0= count(vLKinActiv)
    !!     nEqA0= count(vLEquActiv)
    !!     !
    !!     call Dynam_OneStep_Clean
    !!     !
    !!     if(LogForAqu) then
    !!       call Dynam_OneStep_Alloc(nAq +nMkA0 +nEqA0)
    !!     else
    !!       call Dynam_OneStep_Alloc(nAq -nAs +nMkA0 +nEqA0)
    !!     end if
    !!     !
    !!   end if
    !! end if
    !-------------------------------------------------/change dimensions
    !
    !---------------------------------------------/update kin'phase data
    !
    if(dTimeTooSmall) exit DoStep !--------------------------exit DoStep
    !
    !! !---------- message output for COORES ----------------------!
    !! if (DebugCoores) &
    !! & print '(1A6, 1I6, 1A, 1F14.6, 1A, 1G14.6, 1A, 1I3, 1A, 1G14.6, 1A)', &
    !! & '[ STEP', iStep,'] - [T = ',time,',DT =', dtime,', NBit=', Newt_iDo,' ][',VarMAx,']'
    !! !--- store variable and output flux for Coores -----
    !! if (LSTOCK) &
    !! & call set_stockvar( vMolF (vPrmBk(1:nAq)) , Fout/VBox, Time, dTime, VarMax , 'C' )
    !! !-----------------------------------------------------------!
    !
    !------------------------save control parameters for the time step--
    Dynam_nStep=          iStep                                        !
    Dynam_nNewtIter=      Newt_iDo                                     !
    Dynam_nTotalNewtIter= Dynam_nTotalNewtIter  + Newt_iDo             !
    !------------------------------------------------------------------/
    !
    call Dynam_OneStep_Variations( &
    & F3,iStep,Time,dTime, &
    & vMolF_0,vMolK_0,vMolF,vMolK, &
    & MaxAqu,MaxMin)
    !
    call Dynam_OneStep_Write( &
    & iStep,Time,dTime,       &
    & fSavRate,fSavTime,      &
    & PhiF,RhoF,VBox,dX,      &
    & vTotF,vMolF,vLnAct,vLnGam, &
    & vTotK,vMolK,vMolarVol,vKinQsk,vSurfK,vSurfK0,vVmQsK,vVmAct, &
    & dTSav,TSave)
    !
    !--save current vSurfK, vMolK, PhiF for surface calculations
    !vSurfK0(:)=  vSurfK(:)  !update vSurfK0
    !PhiF0=       PhiF       !update PhiF0
    !vMolK0(:)=   vMolK(:)   !update vMolK0
    !
    if(iDebug>2) call Dynam_Files_WriteFasAff(iStep,Time,vLnAct,vLnBuf)
    !
    !---------------------------------------------check for steady state
    if(SteadyState_Stop) then
      !
      if(VarMax<ZeroVarMax) NZeroVar= NZeroVar + 1
      if(VarMax>ZeroVarMax) NZeroVar= 0 !counter reset
      !ensure ZeroVarMax is true on several consecutive time steps
      !
      if(NZeroVar>NZeroMin .or. Time>=TFinal) then
        !-------------------------------------------ask for prolongation
        if(YesAnswer("Continue")) then
          NZeroVar= 0
          TFinal= TFinal + TAgain
          !-> extend end of calc'
        else
          Time= Tfinal
          exit DoStep
        end if
      end if
      !
    end if
    !--------------------------------------------/check for steady state
    !
    !-----------------------------------------------Time>=TFinal -> exit
    if(Time>=TFinal) then
      Time= Tfinal
      exit DoStep
    end if
    !-----------------------------------------------------------------//
    !
    if(UpdateMassFluid .and. (PhiF<1.0D-6)) then
      PorosTooSmall=.true.
      exit DoStep !--------------------------PorosTooSmall > exit DoStep
    end if
    !
    !----------------------------------------calc. new time step-> dTime
    if(Time<TFinal) call Dynam_CtrlDtime( &
    & NewtIterMax,NewtIterMin, &
    & Time_Decrease,Time_Increase, &
    & iCtrlTime,Newt_iDo, &
    & Time,TFinal,dTMax,MaxAqu,MaxMin, &
    & dTime)
    !
    if(dTime<dTmin) then
      DTimeTooSmall=.true.
      exit DoStep
    end if
    !
  end do DoStep !->Time>=tFinal, end of run
  !--------------------------------------------------------------/DoStep
  !
  if(iDebug>2) call CPU_TIME(CpuEnd)
  if(iDebug>2) print *,"CPU_TIME=",CpuEnd-CpuBegin
  if(iDebug>2) print *,"nEvalFunc=", nEvalFunc
  !
  !------------------------------------------- message output for COORES
  if(DebugCoores) &
  & print '(1A6, 1I6, 1A, 1F14.6, 1A)', &
  & '[ STEP',iStep,'] - [T = ',time,', ------ end      ]'
  !-------------------------------------------------------------------//
  !
  if(PorosTooSmall) then
    print '(A)',"Porosity Too Low ..."
    if(idebug>1) write(fTrc,'(A)') "Porosity Too Low !!!"
  end if
  if(DTimeTooSmall) then
    print '(A,G12.3)',"Time Step= ", dTmin
    print '(A)',"Time Step Too Small !!!"
    if(idebug>1) write(fTrc,'(A)') "Time Step Too Small !!!"
  end if
  
  return
end subroutine Dynam_Calc_Main

subroutine Dynam_Calc_Close
  !
  deallocate(vMolF_0)
  deallocate(vMolK_0)
  deallocate(vTotF_0)
  deallocate(vTotK)
  deallocate(vMolM_Slope)
  !
  ! call Dynam_OneStep_Clean
  call Dynam_Solve_Vars_Clean
  !
  !cpu! if(iDebug>2) print      '(A,G15.3)', "CPU-Time (sec)=",CpuEnd-CpuBegin
  !cpu! if(idebug>1) write(fTrc,'(A,G15.3)') "CPU-Time (sec)=",CpuEnd-CpuBegin
  !
  if(f1>0) close(f1)
  if(f2>0) close(f2)
  if(f3>0) close(f3)
  if(f4>0) close(f4)
  !
  if(idebug>1) write(fTrc,'(A,/)') "</ CalcDynamic"
  
  return
end subroutine Dynam_Calc_Close

logical function YesAnswer(S)
  use M_VarStr,only: T_VarStr,VarStr_Get,VarStr_Char
  
  character(len=*),intent(in):: S
  !
  type(T_VarStr)::Str
  
  write(*,'(A)',advance="no") trim(S)//" (Y/N) ? "
  call VarStr_GET(STRING=Str)
  YesAnswer= (VarStr_Char(Str)=="Y" .or. VarStr_Char(Str)=="y")
  
  return
end function YesAnswer

logical function TP_Changed(TdgK,Pbar)
  use M_Dynam_Vars,only: TdgK0,Pbar0,TolTdgK,TolPbar
  real(dp),intent(in):: TdgK,Pbar
  !
  TP_Changed= ( abs(TdgK - TdgK0)>TolTdgK .or. abs(Pbar - Pbar0)>TolPbar)
  if(TP_Changed) then
    TdgK0= TdgK
    Pbar0= Pbar
  end if
  
end function TP_Changed

!wrk! character function CharFromKbd(S)
!wrk!   use M_VarStr,only: T_VarStr,VarStr_Get,VarStr_Char
!wrk!   !
!wrk!   character(len=*),intent(in):: S
!wrk!   !
!wrk!   type(T_VarStr)  :: Str
!wrk!   character(len=1):: CC
!wrk!   !
!wrk!   write(*,'(A)',advance="no") trim(S) !!//" (Y/N) ? "
!wrk!   call VarStr_GET(STRING=Str)
!wrk!   CC=VarStr_Char(Str)
!wrk!   CharFromKbd= CC(1:1)
!wrk!   return
!wrk! end function CharFromKbd

end module M_Dynam_Calc
