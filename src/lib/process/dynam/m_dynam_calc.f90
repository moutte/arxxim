MODULE M_Dynam_Calc
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_
  USE M_Numeric_Const, ONLY: Ln10
  USE M_IOTools
  USE M_Trace, ONLY: DebugCoores
  !
  !! USE M_Stockvar_Kinxim,ONLY: LSTOCK,SET_STOCKVAR
  !
  USE M_Global_Vars,ONLY: vSpc,vKinFas
  USE M_System_Vars,ONLY: TdgK,Pbar
  USE M_Basis_Vars, ONLY: nCi,nAs,vPrmBk
  USE M_Basis_Vars, ONLY: tAlfPr,tAlfAs
  !
  USE M_Dynam_Files,ONLY: Dynam_Files_OpenLogs,Dynam_Files_WriteFasAff
  USE M_Dynam_Vars, ONLY: Time,TFinal,dTime,dTmin,dTMax,TUnit,dTSav,TimeFactor
  USE M_Dynam_Vars, ONLY: TdgK0,Pbar0
  USE M_Dynam_Vars, ONLY: UDarcy,VBox,dX,RhoF
  USE M_Dynam_Vars, ONLY: FOut,PhiF0,PhiF,PhiInert
  USE M_Dynam_Vars, ONLY: vMolF,vTotF,vTotInj
  USE M_Dynam_Vars, ONLY: vLnAct,vLnBuf,vLnGam,pH_
  USE M_Dynam_Vars, ONLY: vMolarVol,vKinQsk,vStatusK
  USE M_Dynam_Vars, ONLY: tAlfKin,vLKinActiv,vLEquActiv,vKinPrm
  USE M_Dynam_Vars, ONLY: vMolK,vMolK0,vSurfK,vSurfK0
  USE M_Dynam_Vars, ONLY: vVmQsK,vVmAct,UpdateMassFluid
  USE M_Dynam_Vars, ONLY: &
  & Dynam_nStep,     &
  & Dynam_nNewtIter, &
  & Dynam_nTotalNewtIter, &
  & nEvalFunc
  USE M_Dynam_Vars, ONLY: Extrapole
  USE M_Dynam_Vars, ONLY: SteadyState_Stop
  !
  USE M_Dynam_Tools,ONLY: Dynam_TP_Update
  USE M_Dynam_OneStep_Tools
  USE M_Dynam_OneStep_Solve
  !
  USE M_Dynam_Solve_Vars,ONLY: Dynam_Solve_Vars_Alloc, Dynam_Solve_Vars_Clean
  !
  IMPLICIT NONE

  PRIVATE

  PUBLIC:: Dynam_Calc
  PUBLIC:: Dynam_Calc_Init
  PUBLIC:: Dynam_Calc_Main
  PUBLIC:: Dynam_Calc_Close

  !---------------------------------------------------------------vars--
  REAL(dp),DIMENSION(:),ALLOCATABLE:: & !
  & vMolF_0, & !species mole nrs in fluid, prev.step -> to compute variations
  & vMolK_0, & !mole nrs of minerals, prev.step -> to compute variations
  & vTotF_0, & !component mole nrs in fluid, prev.step -> to compute variations
  & vTotK,   & !component mole nrs in minerals
  & vTotK_0, & !component mole nrs in minerals, prev.step -> to compute variations
  & vMolM_Slope
  !
  LOGICAL,PUBLIC:: DTimeTooSmall,PorosTooSmall
  !
  REAL(dp):: PhiF_0,X !,UDarcyOut
  REAL(dp):: VarMax,dTimeOut
  REAL(dp):: MaxAqu=Zero, MaxMin=Zero
  !
  REAL(dp):: ZeroVarMax= 1.0D-15
  !
  INTEGER :: nCp,nMk,nAq
  ! INTEGER :: nMkA0,nEqA0
  INTEGER :: f1,f2,f3,f4
  INTEGER :: iStep,Newt_iDo,NZeroVar,Newt_iErr
  INTEGER :: NZeroMin= 12
  !
  REAL(dp):: TSave
  REAL(dp):: TAgain
  ! TAgain= duration of an additional simulation starting at time TFinal
  !
  REAL(dp):: CpuBegin,CpuEnd
  !--------------------------------------------------------------------/
  
CONTAINS

SUBROUTINE Dynam_Calc
  CALL Dynam_Calc_Init
  CALL Dynam_Calc_Main
  CALL Dynam_Calc_Close
END SUBROUTINE Dynam_Calc

SUBROUTINE Dynam_Calc_Init
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< CalcDynamic"
  !
  SteadyState_Stop= (iDebug>2)
  IF(DebugCoores) SteadyState_Stop= .FALSE.
  !
  nCp= SIZE(vTotF)
  nMk= SIZE(vKinFas)
  nAq= COUNT(vSpc%Typ=="AQU")
  !
  CALL Dynam_Solve_Vars_Alloc(nMk,nCi)
  !
  ALLOCATE(vMolF_0(1:nAq))
  ALLOCATE(vMolK_0(1:nMk))
  ALLOCATE(vTotF_0(1:nCp))
  !
  ALLOCATE(vTotK(1:nCp))
  ALLOCATE(vTotK_0(1:nCp))
  !
  ALLOCATE(vMolM_Slope(1:nMk))
  vMolM_Slope(:)= Zero
  !
  f1= 0  ;  f2= 0
  f3= 0  ;  f4= 0
  IF(iDebug>2) CALL Dynam_Files_OpenLogs(f1,f2,f3,f4)
  !
  TAgain= TFinal
  TSave=Zero !------------------for output at definite time span dTSav--
  !-------------------------------(provisional, needs improvement ...)--
  !
  DTimeTooSmall= .FALSE.
  PorosTooSmall= .FALSE.
  !
  NZeroVar= 0
  Newt_iDo= 0
  VarMax=   Zero
  !
  !! !---------- message output for Coores --------------
  !! IF (DebugCoores) &
  !! & PRINT '(1A6, 1I6, 1A, 1F14.6, 1A)', &
  !! & '[ STEP',iStep,'] - [T = ',time,', ------ START    ]'
  !! !--- store variable and output flux for Coores -----
  !! IF (LSTOCK) &
  !! & CALL set_stockvar( vMolF (vPrmBk(1:nAq)) , Fout/VBox, Time, dTime, VarMax , 'D' )
  !! !---------------------------------------------------
  !
  !------------------------------------------initialize kinetic phases--
  X= MAXVAL(vMolF(1:nAq))
  CALL Dynam_OneStep_KinFas_Update( & !
  & LOG(vMolF(1:nAq)/X),   & !in
  & vMolK,                 & !in
  & vKinQsk,vStatusK,      & !out
  & vVmQsK,vVmAct,         & !out
  & vLKinActiv,vLEquActiv, & !out
  & vKinPrm,PhiInert)        !out
  !
  CALL Dynam_OneStep_KinFasSurf_Update( & !
  & vLKinActiv,  & ! IN
  & vStatusK,    & ! IN
  & vMolK0,vMolK,& ! IN
  & PhiF0,PhiF,  & ! IN
  & vSurfK0,     & ! IN
  & vSurfK)        ! OUT
  !
  !!vVm0(:)= vSurfK(:) *vVmAct(:) *vVmQsK(:)
  !
  !! nMkA0= COUNT(vLKinActiv)
  !! nEqA0= COUNT(vLEquActiv)
  !! !~ !
  !! IF(LogForAqu) THEN
  !!   CALL Dynam_OneStep_Alloc(nAq +nMkA0 +nEqA0)
  !! ELSE
  !!   CALL Dynam_OneStep_Alloc(nAq +nMkA0 +nEqA0)
  !! ENDIF
  !
  !----------------------------------------/initialize kinetic phases--
  !
  CALL Dynam_OneStep_TotK( & !
  & vMolK, & !in
  & vTotK)   !out
  !
  CALL Dynam_OneStep_TotF( & !
  & vMolF, & !in
  & vTotF)   !out
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "</CalcDynamic"
  !
END SUBROUTINE Dynam_Calc_Init

SUBROUTINE Dynam_Calc_Main
  USE M_Dynam_Vars,ONLY: NewtTolF,NewtTolF_Equil,NewtTolX,NewtTolMin
  USE M_Dynam_Vars,ONLY: NewtMaxIts,NewtIterMax,NewtIterMin
  USE M_Dynam_Vars,ONLY: Time_Decrease,Time_Increase,iCtrlTime
  USE M_Dynam_Vars,ONLY: fSavRate,fSavTime,DebNewt
  !----------------------------------------------for test basis change--
  USE M_System_Vars,ONLY: vCpn
  USE M_Basis,ONLY: Basis_Change_Wrk
  USE M_Basis_Vars,ONLY: isW
  !--------------------------------------------------------------------/
  
  IF(iDebug>2) CALL CPU_TIME(CpuBegin)
  IF(iDebug>2) nEvalFunc= 0
  !
  IF(iDebug>1) PRINT '(/,A,/)', &
  & "OUTPUT= iStep,Time,dTime,PhiF,nEqA,nMkA,iDo2,Newt_iDo,Newt_iErr, &
  & iMaxDelta,VarMax,NewtErrF,NewtErrX"
  !
  iStep= 0
  !-------------------------------------------------------------DoStep--
  DoStep: DO !WHILE (Time<tFinal)
    !
    iStep= iStep+1
    DebNewt= (iStep<4) !WRITE jacobian for the first 4 steps ...
    !IF(Time<TFinal) CALL CtrlDtime(iCtrlTime,Newt_iDo,MaxAqu,MaxMin,dTime)
    !
    !--------adjust mole numbers to available volume, based on density--
    IF(UpdateMassFluid) CALL Dynam_OneStep_Adjust( & !
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
    IF(TP_Changed(TdgK,Pbar)) CALL Dynam_TP_Update(TdgK,Pbar,TimeFactor)
    !
    !---------------------------------------------------solve one step--
    CALL Dynam_OneStep_Solve(       & !
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
    IF(iDebug>2) CALL Basis_Change_Wrk(isW,vCpn)
    !------------------------------------------------------/basis change
    !
    IF(iDebug>2) &
    & CALL Dynam_OneStep_WriteElements(iStep,Time,dTime,vMolF,vMolK)
    !
    CALL Dynam_OneStep_TotF( & !
    & vMolF, & !in
    & vTotF)   !out
    !
    CALL Dynam_OneStep_TotK( & !
    & vMolK, & !in
    & vTotK)   !out
    !
    PhiF= One - SUM(vMolK(:)*vMolarVol(:)) /Vbox
    !
    IF(Extrapole) &
    & vMolM_Slope(:)= (vMolK(:) -vMolK_0(:)) /dTime
    !vMolM_Slope(:)= Zero
    !
    IF(iStep>3 .AND. F1>0) &
    & CALL Dynam_OneStep_Write_Bilan( & !
    & F1,iStep,Time,dTimeOut,               &  !IN
    & PhiF,UDarcy,VBox,dX,                  &  !IN
    & vTotInj,vTotF_0,vTotF,vTotK_0,vTotK)     !IN
    !
    !----------------------------------------------update kin'phase data
    X= MAXVAL(vMolF(1:nAq))
    !
    CALL Dynam_OneStep_KinFas_Update( & !
    & LOG(vMolF(1:nAq)/X),   & !in
    & vMolK,                 & !in
    & vKinQsk,vStatusK,      & !out
    & vVmQsK,vVmAct,         & !out
    & vLKinActiv,vLEquActiv, & !out
    & vKinPrm,PhiInert)        !out
    !
    CALL Dynam_OneStep_KinFasSurf_Update( & !
    & vLKinActiv,  & ! IN
    & vStatusK,    & ! IN
    & vMolK0,vMolK,& ! IN
    & PhiF0,PhiF,  & ! IN
    & vSurfK0,     & ! IN
    & vSurfK)        ! OUT
    !
    !--------------------------------------------------change dimensions
    !! IF(OnlyActiv) THEN
    !!   IF(COUNT(vLKinActiv)/=nMkA0 .OR. COUNT(vLEquActiv)/=nEqA0) THEN
    !!     !
    !!     nMkA0= COUNT(vLKinActiv)
    !!     nEqA0= COUNT(vLEquActiv)
    !!     !
    !!     CALL Dynam_OneStep_Clean
    !!     !
    !!     IF(LogForAqu) THEN
    !!       CALL Dynam_OneStep_Alloc(nAq +nMkA0 +nEqA0)
    !!     ELSE
    !!       CALL Dynam_OneStep_Alloc(nAq -nAs +nMkA0 +nEqA0)
    !!     ENDIF
    !!     !
    !!   ENDIF
    !! ENDIF
    !-------------------------------------------------/change dimensions
    !
    !---------------------------------------------/update kin'phase data
    !
    !! vVm0(:)= vSurfK(:) *vVmAct(:) *vVmQsK(:)
    !
    IF(dTimeTooSmall) EXIT DoStep !--------------------------EXIT DoStep
    !
    !! !---------- message output for COORES ----------------------!
    !! IF (DebugCoores) &
    !! & PRINT '(1A6, 1I6, 1A, 1F14.6, 1A, 1G14.6, 1A, 1I3, 1A, 1G14.6, 1A)', &
    !! & '[ STEP', iStep,'] - [T = ',time,',DT =', dtime,', NBit=', Newt_iDo,' ][',VarMAx,']'
    !! !--- store variable and output flux for Coores -----
    !! IF (LSTOCK) &
    !! & CALL set_stockvar( vMolF (vPrmBk(1:nAq)) , Fout/VBox, Time, dTime, VarMax , 'C' )
    !! !-----------------------------------------------------------!
    !
    !------------------------save control parameters for the time step--
    Dynam_nStep=          iStep                                        !
    Dynam_nNewtIter=      Newt_iDo                                     !
    Dynam_nTotalNewtIter= Dynam_nTotalNewtIter  + Newt_iDo             !
    !------------------------------------------------------------------/
    !
    CALL Dynam_OneStep_Variations( &
    & F3,iStep,Time,dTime, &
    & vMolF_0,vMolK_0,vMolF,vMolK, &
    & MaxAqu,MaxMin)
    !
    CALL Dynam_OneStep_Write( &
    & iStep,Time,dTime,       &
    & fSavRate,fSavTime,      &
    & pH_,PhiF,RhoF,VBox,dX,  &
    & vTotF,vMolF,vLnAct,vLnGam, &
    & vTotK,vMolK,vMolarVol,vKinQsk,vSurfK,vSurfK0,vVmQsK,vVmAct, &
    & dTSav,TSave)
    !
    !--save current vSurfK, vMolK, PhiF for surface calculations
    !vSurfK0(:)=  vSurfK(:)  !update vSurfK0
    !PhiF0=       PhiF       !update PhiF0
    !vMolK0(:)=   vMolK(:)   !update vMolK0
    !
    IF(iDebug>2) CALL Dynam_Files_WriteFasAff(iStep,Time,vLnAct,vLnBuf)
    !
    !---------------------------------------------check for steady state
    IF(SteadyState_Stop) THEN
      !
      IF(VarMax<ZeroVarMax) NZeroVar= NZeroVar + 1
      IF(VarMax>ZeroVarMax) NZeroVar= 0 !counter reset
      !ensure ZeroVarMax is true on several consecutive time steps
      !
      IF(NZeroVar>NZeroMin .OR. Time>=TFinal) THEN
        !-------------------------------------------ask for prolongation
        IF(YesAnswer("Continue")) THEN
          NZeroVar= 0
          TFinal= TFinal + TAgain
          !-> extend end of calc'
        ELSE
          Time= Tfinal
          EXIT DoStep
        ENDIF
      ENDIF
      !
    ENDIF
    !--------------------------------------------/check for steady state
    !
    !-----------------------------------------------Time>=TFinal -> exit
    IF(Time>=TFinal) THEN
      Time= Tfinal
      EXIT DoStep
    ENDIF
    !-------------------------------------------------------------------
    !
    IF(UpdateMassFluid .AND. (PhiF<1.0D-6)) THEN
      PorosTooSmall=.TRUE.
      EXIT DoStep !--------------------------PorosTooSmall > EXIT DoStep
    ENDIF
    !
    !----------------------------------------calc. new time step-> dTime
    IF(Time<TFinal) CALL Dynam_CtrlDtime( &
    & NewtIterMax,NewtIterMin, &
    & Time_Decrease,Time_Increase, &
    & iCtrlTime,Newt_iDo, &
    & Time,TFinal,dTMax,MaxAqu,MaxMin, &
    & dTime)
    !
    IF(dTime<dTmin) THEN
      DTimeTooSmall=.TRUE.
      EXIT DoStep
    ENDIF
    !
  ENDDO DoStep !->Time>=tFinal, end of run
  !--------------------------------------------------------------/DoStep
  !
  IF(iDebug>2) CALL CPU_TIME(CpuEnd)
  IF(iDebug>2) PRINT *,"CPU_TIME=",CpuEnd-CpuBegin
  IF(iDebug>2) PRINT *,"nEvalFunc=", nEvalFunc
  !
  !------------------------------------------- message output for COORES
  IF(DebugCoores) &
  & PRINT '(1A6, 1I6, 1A, 1F14.6, 1A)', &
  & '[ STEP',iStep,'] - [T = ',time,', ------ END      ]'
  !--------------------------------------------------------------------/
  !
  IF(PorosTooSmall) THEN
    PRINT '(A)',"Porosity Too Low ..."
    IF(iDebug>0) WRITE(fTrc,'(A)') "Porosity Too Low !!!"
  ENDIF
  IF(DTimeTooSmall) THEN
    PRINT '(A,G12.3)',"Time Step= ", dTmin
    PRINT '(A)',"Time Step Too Small !!!"
    IF(iDebug>0) WRITE(fTrc,'(A)') "Time Step Too Small !!!"
  ENDIF
  
  RETURN
END SUBROUTINE Dynam_Calc_Main

SUBROUTINE Dynam_Calc_Close
  !
  DEALLOCATE(vMolF_0)
  DEALLOCATE(vMolK_0)
  DEALLOCATE(vTotF_0)
  DEALLOCATE(vTotK)
  DEALLOCATE(vMolM_Slope)
  !
  ! CALL Dynam_OneStep_Clean
  CALL Dynam_Solve_Vars_Clean
  !
  !cpu! IF(iDebug>2) PRINT      '(A,G15.3)', "CPU-Time (sec)=",CpuEnd-CpuBegin
  !cpu! IF(iDebug>0) WRITE(fTrc,'(A,G15.3)') "CPU-Time (sec)=",CpuEnd-CpuBegin
  !
  IF(f1>0) CLOSE(f1)
  IF(f2>0) CLOSE(f2)
  IF(f3>0) CLOSE(f3)
  IF(f4>0) CLOSE(f4)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ CalcDynamic"
  
  RETURN
END SUBROUTINE Dynam_Calc_Close

LOGICAL FUNCTION YesAnswer(S)
  USE M_VarStr,ONLY: T_VarStr,VarStr_Get,VarStr_Char
  
  CHARACTER(LEN=*),INTENT(IN):: S
  !
  TYPE(T_VarStr)::Str
  
  WRITE(*,'(A)',ADVANCE='NO') TRIM(S)//" (Y/N) ? "
  CALL VarStr_GET(STRING=Str)
  YesAnswer= (VarStr_Char(Str)=="Y" .OR. VarStr_Char(Str)=="y")
  
  RETURN
ENDFUNCTION YesAnswer

LOGICAL FUNCTION TP_Changed(TdgK,Pbar)
  USE M_Dynam_Vars,ONLY: TdgK0,Pbar0,TolTdgK,TolPbar
  REAL(dp),INTENT(IN):: TdgK,Pbar
  !
  TP_Changed= ( ABS(TdgK - TdgK0)>TolTdgK .OR. ABS(Pbar - Pbar0)>TolPbar)
  IF(TP_Changed) THEN
    TdgK0= TdgK
    Pbar0= Pbar
  ENDIF
  
ENDFUNCTION TP_Changed

!wrk! CHARACTER FUNCTION CharFromKbd(S)
!wrk!   USE M_VarStr,ONLY: T_VarStr,VarStr_Get,VarStr_Char
!wrk!   !
!wrk!   CHARACTER(LEN=*),INTENT(IN):: S
!wrk!   !
!wrk!   TYPE(T_VarStr)  :: Str
!wrk!   CHARACTER(LEN=1):: CC
!wrk!   !
!wrk!   WRITE(*,'(A)',ADVANCE='NO') TRIM(S) !!//" (Y/N) ? "
!wrk!   CALL VarStr_GET(STRING=Str)
!wrk!   CC=VarStr_Char(Str)
!wrk!   CharFromKbd= CC(1:1)
!wrk!   RETURN
!wrk! ENDFUNCTION CharFromKbd

ENDMODULE M_Dynam_Calc
