MODULE M_Dynam_OneStep_Tools
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC:: Dynam_OneStep_Write
  PUBLIC:: Dynam_OneStep_WriteElements
  PUBLIC:: Dynam_OneStep_TotF
  PUBLIC:: Dynam_OneStep_TotK
  PUBLIC:: Dynam_OneStep_Adjust
  PUBLIC:: Dynam_OneStep_KinFas_Update
  PUBLIC:: Dynam_OneStep_KinFasSurf_Update
  PUBLIC:: Dynam_OneStep_SecondSpecies
  PUBLIC:: Dynam_OneStep_GammaUpd
  PUBLIC:: Dynam_OneStep_DLGammaUpd
  PUBLIC:: Dynam_OneStep_Variations
  PUBLIC:: Dynam_OneStep_Write_Bilan
  PUBLIC:: Dynam_CtrlDtime

CONTAINS

SUBROUTINE Dynam_OneStep_SecondSpecies(N,vX,vXsec)
  USE M_Basis_Vars,ONLY: nAs,nAx,nMx,tNuAs,nCi,MWSv
  USE M_Dynam_Vars,ONLY: vDG_As,vDLGam_As,vLnGam,vLnAct,LogForAqu,vLnBuf
  !
  INTEGER, INTENT(IN) :: N
  REAL(dp),INTENT(IN) :: vX(:)
  REAL(dp),INTENT(OUT):: vXsec(:)
  !
  INTEGER :: iAs,iCi,J
  REAL(dp):: X
  !
  IF(LogForAqu) THEN
    !
    DO iAs=1,nAs
      !
      X=   vLnAct(1) *tNuAs(iAs,1) &
      &  + (One - SUM(tNuAs(iAs,2:nCi)))*(vX(1) +LOG(MWSv)) &
      &  - vDG_As(iAs) -vDLGam_As(iAs)
      DO iCi=1,N
        IF(iCi/=1) X= X + vX(iCi) *tNuAs(iAs,iCi)
      ENDDO
      DO J=1,nAx+nMx
        X=  X + vLnBuf(J) *tNuAs(iAs,nCi+J)
      ENDDO
      !
      vXsec(iAs)= EXP(X)
      != (mole nr species iAs)
      !
    END DO
    !
  ELSE
    !
    DO iAs=1,nAs
      !
      X=  vLnAct(1) *tNuAs(iAs,1) &
      & - vDG_As(iAs) -vDLGam_As(iAs)
      DO J=1,nAx+nMx
        X= X + vLnBuf(J) *tNuAs(iAs,nCi+J)
      ENDDO
      X= EXP(X)
      !
      DO iCi=1,N
        IF(iCi/=1 .AND. tNuAs(iAs,iCi) /= Zero) &
        & X= X *(vX(iCi))**tNuAs(iAs,iCi)
      END DO
      vXsec(iAs)= X *(vX(1)*MWSv)**(One-SUM(tNuAs(iAs,2:nCi)))
      != (mole nr species iAs)
      !
    END DO
    !
  ENDIF
    !
  !---------------------------------------------------------- trace --
  if(iDebug>2) then
    if(LogForAqu) then
      write(71,'(A)') "====================SOLVE=="
      do j=1,N
        write(71,'(A,G15.6)') "resPrim=",EXP(vX(j))
      enddo
      do j=1,size(vXsec)
        write(71,'(A,G15.6)') "resSec==",EXP(vXsec(j))
      enddo
    else
      write(71,'(A)') "====================SOLVE=="
      do j=1,N
        write(71,'(A,G15.6)') "resPrim=",vX(j)
      enddo
      do j=1,size(vXsec)
        write(71,'(A,G15.6)') "resSec==",vXsec(j)
      enddo
    endif
  endif
  !---------------------------------------------------------/ trace --
  !pause
  !
END SUBROUTINE Dynam_OneStep_SecondSpecies

SUBROUTINE Dynam_OneStep_Write( &
& iStep,Time,dTime,       &
& fSavRate,fSavTime,      &
& pH_,PhiF,RhoF,VBox,dX,  &
& vTotF,vMolF,vLnAct,vLnGam, &
& vTotK,vMolK,vKinMVol,vKinQsk,vSurfK,vSurfK0,vVmQsK,vVmAct, &
& dTSav,TSave)
  !------
  USE M_Numeric_Const,ONLY: Ln10
  USE M_IoTools,      ONLY: OutStrVec
  USE M_Safe_Functions
  !
  USE M_Dynam_Files,ONLY: fDynEle,fDynMol,fDynAct,fDynGam,fDynMnK
  !
  USE M_Global_Vars,ONLY: vKinFas
  USE M_Basis_Vars, ONLY: isW,MWsv
  USE M_Dynam_Vars, ONLY: vWeitCp,TimeScale
  !------
  INTEGER,              INTENT(IN):: iStep
  REAL(dp),             INTENT(IN):: Time,dTime
  INTEGER,              INTENT(IN):: fSavRate,fSavTime
  REAL(dp),             INTENT(IN):: pH_,PhiF,RhoF,VBox,dX
  REAL(dp),DIMENSION(:),INTENT(IN):: vTotF,vMolF,vLnAct,vLnGam
  REAL(dp),DIMENSION(:),INTENT(IN):: vTotK,vMolK,vKinMVol,vKinQsk
  REAL(dp),DIMENSION(:),INTENT(IN):: vSurfK,vSurfK0,vVmQsK,vVmAct
  REAL(dp),             INTENT(IN):: dTSav
  !
  REAL(dp),             INTENT(INOUT):: TSave
  !
  INTEGER :: nAq,I
  REAL(dp):: R,UDelta,VolF
  REAL(dp):: vX(SIZE(vLnGam))
  !------
  nAq= SIZE(vMolF)
  !
  VolF= DOT_PRODUCT(vWeitCp(:),vTotF(:))  /RhoF  !volume of fluid in box
  UDelta= (VolF -PhiF) /dTime /(VBox/dX)
  R= vTotF(isW)*MWsv
  !
  !-------------- write mole nrs of elements in fluid and in minerals --
  IF(fDynEle>0) THEN
  WRITE(fDynEle,'(I7,A1,2(G15.6,A1),$)') iStep,T_,Time/TimeScale,T_,UDelta,T_
  CALL OutStrVec(fDynEle, vTotF(:),      Opt_C="G", CR=.FALSE.) !->mole nrs in fluid
  CALL OutStrVec(fDynEle, vTotF(:)/R,    Opt_C="G", CR=.FALSE.) !->molality
  CALL OutStrVec(fDynEle, vTotK(:),      Opt_C="G")
  END IF
  !
  !-- nr'mole component injected in time step:
  !-- - vTotInj(iPr) *UDarcy /dX *dTime
  !
  !--------------------------------------- write mole nrs aqu'species --
  IF(fDynMol>0) &
  & CALL OutStrVec(fDynMol, vMolF(1:nAq), Opt_C="G", Opt_I=iStep, Opt_R=Time/TimeScale)
  !
  !------------------------------------ write activities aqu'species ---
  IF(fDynAct>0) &
  & CALL OutStrVec(fDynAct, vLnAct(1:nAq)/Ln10,Opt_C="F", Opt_I=iStep, Opt_R=Time/TimeScale)
  !
  !----------------------------------- write activ'coeffs aqu'species --
  DO I=1,nAq
    vX(i)= FSafe_Exp(vLnGam(i))
  ENDDO
  IF(fDynGam>0) &
  & CALL OutStrVec(fDynGam, vX(1:nAq),Opt_C="F", Opt_I=iStep, Opt_R=Time/TimeScale)
  !
  !------------------ write fluid and min'species fractions, QsK, etc --
  IF(fDynMnK>0) &
  & CALL ShoMin(fDynMnK,iStep,VBox,Time/TimeScale,pH_,PhiF,vMolK,vKinMVol,vKinQsk)
  !
  IF(fSavRate>0) &
  & CALL ShoMinRate( &
  & fSavRate,iStep,VBox,Time/TimeScale,pH_,PhiF, &
  & vMolK,vKinMVol,vKinQsk, &
  & vSurfK,vSurfK0,vVmQsK,vVmAct, &
  & vKinFas)
  !
  IF(dTSav>Zero .AND. Time>=TSave) THEN
    TSave= TSave +dTSav
    CALL ShoMin(fSavTime,iStep,VBox,Time/TimeScale,pH_,PhiF,vMolK,vKinMVol,vKinQsk)
  ENDIF
  !
  RETURN
ENDSUBROUTINE Dynam_OneStep_Write

SUBROUTINE Dynam_OneStep_WriteElements( &
& iStep,Time,dTime, &
& vMolF,vMolK)
  USE M_Dynam_Files,ONLY: fDynElement
  USE M_Dynam_Vars, ONLY: tStoikioKin,tStoikioAqu,TimeScale
  !
  INTEGER, INTENT(IN):: iStep
  REAL(dp),INTENT(IN):: Time,dTime
  REAL(dp),INTENT(IN):: vMolF(:)
  REAL(dp),INTENT(IN):: vMolK(:)
  !
  REAL(dp),ALLOCATABLE:: vX(:)
  INTEGER :: iEl,iSp,N
  REAL(dp):: X
  !
  N= SIZE(tStoikioAqu,2) != number of elements
  ALLOCATE(vX(N))
  !
  WRITE(fDynElement,'(I7,A1,G15.6,A1,$)') &
  & iStep,T_,Time/TimeScale,T_
  !
  DO iEl=1,N
    X= Zero
    DO iSp=1,SIZE(vMolF)
      X= X + tStoikioAqu(iSp,iEl)*vMolF(iSp)
    ENDDO
    vX(iEl)= X
    WRITE(fDynElement,'(G15.6,A1)',ADVANCE="NO") X,T_
  ENDDO
  !
  DO iEl=1,N
    X= Zero
    DO iSp=1,SIZE(vMolK)
      X= X + tStoikioKin(iSp,iEl)*vMolK(iSp)
    ENDDO
    vX(iEl)= vX(iEl)+X
    WRITE(fDynElement,'(G15.6,A1)',ADVANCE="NO") X,T_
  ENDDO
  !
  DO iEl=1,N
    WRITE(fDynElement,'(G15.6,A1)',ADVANCE="NO") vX(iEl),T_
  ENDDO
  !
  WRITE(fDynElement,*)
  !
  DEALLOCATE(vX)
  
  RETURN
ENDSUBROUTINE Dynam_OneStep_WriteElements

SUBROUTINE Dynam_OneStep_Write_Bilan( &
& F,iStep,Time,dTime,  &
& PhiF,UDarcy,VBox,dX, &
& vTotInj,vTotF_0,vTotF,vTotK_0,vTotK)
!--
!-- write balances on components
!--
  USE M_IoTools
  !
  INTEGER, INTENT(IN) :: F,iStep
  REAL(dp),INTENT(IN) :: Time,dTime,PhiF,UDarcy,VBox,dX
  REAL(dp),INTENT(IN) :: vTotInj(:),vTotF_0(:),vTotF(:)
  !                      !mole nrs in fluid, prev. & current step
  REAL(dp),INTENT(IN) :: vTotK_0(:),vTotK(:)
  !                      !mole nrs in minerals, prev. & current step
  !
  REAL(dp),ALLOCATABLE:: vX(:)
  !
  ALLOCATE(vX(SIZE(vTotF)))
  !
  WRITE(F,'(I7,A1,G15.6,A1,$)') iStep,T_,Time,T_
  !
  !
  vX(:)= vTotInj(:) *UDarcy /dX
  CALL OutStrVec(F, vX(:), Opt_C="G", CR=.FALSE.)
  !
  vX(:)=  vTotF(:)   + vTotK(:)   &
  &    - (vTotF_0(:) + vTotK_0(:))
  vX(:)= vX(:) /dTime
  CALL OutStrVec(F, vX(:), Opt_C="G", CR=.FALSE.)
  !
  vX(:)= vX(:) -vTotInj(:) *UDarcy /dX
  CALL OutStrVec(F, vX(:), Opt_C="G", CR=.FALSE.)
  !vX(:)=  vTotF(:) -vTotF_0(:) &
  !&    +  vTotK(:) -vTotK_0(:) &
  !&    + (vTotF(:) /PhiF  - vTotInj(:) *dTime *UDarcy /dX)
  !
  WRITE(F,*)
  !
  DEALLOCATE(vX)
  !
  RETURN
ENDSUBROUTINE Dynam_OneStep_Write_Bilan

SUBROUTINE Dynam_OneStep_Variations( &
& F,iStep,Time,dTime, &
& vMolF_0,vMolK_0,vMolF,vMolK,&
& MaxAqu,MaxMin)
!--
!-- compute max.variations during time step
!--
  USE M_Numeric_Tools,ONLY: iMaxLoc_R
  !
  USE M_Global_Vars,ONLY: vSpc,vKinFas
  USE M_Basis_Vars, ONLY: vOrdAq
  !
  INTEGER, INTENT(IN) :: F,iStep
  REAL(dp),INTENT(IN) :: Time,dTime
  REAL(dp),INTENT(IN) :: vMolF_0(:),vMolK_0(:)
  REAL(dp),INTENT(IN) :: vMolF(:),vMolK(:)
  REAL(dp),INTENT(OUT):: MaxAqu,MaxMin
  !
  INTEGER :: iMaxAqu,iMaxMin
  REAL(dp),ALLOCATABLE:: vDeltaF(:),vDeltaK(:)
  !
  !do i=1,size(vMolF)
  !  print '(3(G15.6,1X))',vMolF_0(i),vMolF(i),ABS(vMolF(i) -vMolF_0(i))
  !enddo
  !pause
  !
  !do i=1,size(vMolK)
  !  print '(3(G15.6,1X))',vMolK_0(i),vMolK(i),ABS(vMolK(i) -vMolK_0(i))
  !enddo
  !pause
  !
  !print *,"MAXVAL",  MAXVAL(ABS(vMolF(:) -vMolF_0(:)))     ;   pause
  !
  ALLOCATE(vDeltaF(SIZE(vMolF)))
  IF(SIZE(vKinFas)>0) ALLOCATE(vDeltaK(SIZE(vMolK)))
  !
  vDeltaF(:)= ABS(vMolF(:) -vMolF_0(:))
  IF(SIZE(vKinFas)>0) vDeltaK(:)= ABS(vMolK(:) -vMolK_0(:))
  !
  iMaxAqu= iMaxLoc_R(vDeltaF(:))
  !-> which aqu.species varies most
  IF(SIZE(vKinFas)>0) iMaxMin= iMaxLoc_R(vDeltaK(:))
  !-> which min.species varies most
  MaxAqu= MAXVAL(vDeltaF(:))
  IF(SIZE(vKinFas)>0) MaxMin=  MAXVAL(vDeltaK(:))
  !
  MaxAqu= MaxAqu /vMolF(iMaxAqu)
  IF(SIZE(vKinFas)>0) MaxMin= MaxMin /vMolK(iMaxMin)
  !
  IF(F>0) THEN
  
    IF(SIZE(vKinFas)>0) THEN
      
      WRITE(F,'(2(A,A1),I6,A1,2(G15.6,A1),4(G15.6,A1))') &
      & TRIM(vSpc(vOrdAq(iMaxAqu))%NamSp),T_, &
      & TRIM(vKinFas(iMaxMin)%NamKF),     T_, &
      & iStep,T_,           &
      & Time,T_,  dTime,T_, &
      & MaxAqu,T_,MaxMin,T_,&
      & MaxAqu /vMolF(iMaxAqu),T_,MaxMin /vMolK(iMaxMin),T_
      !& MAXVAL(vDeltaMolF),T_,MAXVAL(vDeltaMolM),T_
      
    ELSE
      
      WRITE(F,'(A,A1,I3,A1,2(G15.6,A1),2(G15.6,A1))') &
      & TRIM(vSpc(vOrdAq(iMaxAqu))%NamSp),T_, &
      & iStep,T_, &
      & Time,T_,  dTime,T_, &
      & MaxAqu,T_, &
      & MaxAqu /vMolF(iMaxAqu),T_
      !& MAXVAL(vDeltaMolF),T_,MAXVAL(vDeltaMolM),T_
    
    ENDIF
  
  ENDIF
  !
  DEALLOCATE(vDeltaF)
  IF(SIZE(vKinFas)>0) DEALLOCATE(vDeltaK)
  !
ENDSUBROUTINE Dynam_OneStep_Variations

SUBROUTINE Dynam_OneStep_TotF( & !
& vMolF,  & !in
& vTotF)    !out
!--
!-- update mole nrs COMPONENTS IN FLUID at time T
!-- will give vTotF, result of previous step used as RHS in FuncPsi(1:nCi)
!--
  USE M_Basis_Vars,ONLY: nCi,tAlfPr,tAlfAs,nAs
  !
  REAL(dp),INTENT(IN) :: vMolF(:)
  REAL(dp),INTENT(OUT):: vTotF(:)
  !
  INTEGER :: nCp
  !
  nCp= SIZE(vTotF)
  !
  !update mole nrs ALL COMPONENTS in Fluid at time T (= result of previous step)
  vTotF(1:nCp)= MATMUL(tAlfPr(1:nCp,1:nCi),vMolF(1:    nCi))     &
  &           + MATMUL(tAlfAs(1:nCp,1:nAs),vMolF(nCi+1:nCi+nAs))
  !
ENDSUBROUTINE Dynam_OneStep_TotF

SUBROUTINE Dynam_OneStep_TotK( & !
& vMolK,  & !in
& vTotK)    !out
!--
!-- update mole nrs COMPONENTS IN MINERALS at time T
!-- vTotK is used for balance computations
!--
  USE M_Dynam_Vars, ONLY: tAlfKin
  !------
  REAL(dp),INTENT(IN) :: vMolK(:)
  REAL(dp),INTENT(OUT):: vTotK(:)
  !------
  INTEGER :: I
  !
  DO I=1,SIZE(vTotK)
    vTotK(I)= DOT_PRODUCT(tAlfKin(I,:),vMolK(:))
  ENDDO
  !
ENDSUBROUTINE Dynam_OneStep_TotK

SUBROUTINE Dynam_OneStep_Adjust( & !
& RhoF,VBox,PhiF,  & !in
& vMolF,           & !inout (IF UpdateMassFluid THEN adjust)
& vTotF)             !inout
!--
!-- IF volume is "not free",
!-- THEN adjust mole numbers to available volume, based on density
!-- will give vTotF, result of previous step used as RHS in FuncPsi(1:nCi)
!--
  USE M_Global_Vars,ONLY: nAq !,nMk
  USE M_Dynam_Vars, ONLY: vWeitSp
  !------
  REAL(dp),INTENT(IN)::    RhoF,VBox,PhiF
  REAL(dp),INTENT(INOUT):: vMolF(:)
  REAL(dp),INTENT(INOUT):: vTotF(:)
  !
  REAL(dp):: R
  !------
  ! vWeitSp(:)=vSpc(vOrdAq(:))%WeitKg
  ! -> Aqu'Species' weitkg sorted according to local Species' order
  !
  ! MFluid= mass of fluid at time T
  !       = DOT_PRODUCT(vWeitSp(1:nAq),vMolF(1:nAq))
  ! MFluidBox= mass of fluid of density RhoF that fit in VBox*PhiF
  !          = RhoF*VBox*PhiF !=mass fluid in box of volume VBox
  !
  ! -> conversion factor
  ! R= MFluidBox / MFluid
  !  = RhoF*VBox*PhiF /DOT_PRODUCT(vWeitSp(1:nAq),vMolF(1:nAq))
  !
  R= RhoF*VBox*PhiF /DOT_PRODUCT(vWeitSp(:),vMolF(:))
  vTotF(:)= vTotF(:) *R !nr moles COMPONENTS in fluid in box
  vMolF(:)= vMolF(:) *R !nr moles SPECIES in fluid in box
  !
ENDSUBROUTINE Dynam_OneStep_Adjust

SUBROUTINE Dynam_OneStep_KinFas_Update( & !
& vLnX,                   & !in
& vMolK,                  & !in
& vKinQsk,vStatusK,       & !out
& vVmQsK,vVmAct,          & !out
& vLKinActiv,vLEquActiv,  & !out
& vKinPrm,PhiInert)         !out
!--
!-- compute mineral Qsk, rate
!--
  USE M_Numeric_Const,ONLY: Ln10
  USE M_KinRate
  USE M_T_KinFas
  !
  USE M_Global_Vars,ONLY: vFas,vKinFas,nAq
  USE M_Basis_Vars, ONLY: nCi,nMx,nAx,isW,MWSv
  !
  USE M_Dynam_Vars, ONLY: &
  & sModelSurf,Implicit_Surface, &
  & VBox,vDG_Kin,tNu_Kin, &
  & vMolarVol,vKinMod,    &
  & vKinMinim,            &
  & vLnAct,vLnGam,vLnBuf, &
  & vMolK0,vSurfK0,PhiF0, &
  & QsK_Iota
  !------
  REAL(dp),INTENT(IN):: vLnX(:)
  REAL(dp),INTENT(IN):: vMolK(:)
  !
  REAL(dp), INTENT(OUT):: vKinQsk(:)
  REAL(dp), INTENT(OUT):: vVmQsK(:)
  REAL(dp), INTENT(OUT):: vVmAct(:)
  LOGICAL,  INTENT(OUT):: vLKinActiv(:)
  LOGICAL,  INTENT(OUT):: vLEquActiv(:)
  INTEGER,  INTENT(OUT):: vKinPrm(:)
  REAL(dp), INTENT(OUT):: PhiInert
  CHARACTER,INTENT(OUT):: vStatusK(:)
  !
  REAL(dp):: dQsKdLnXi(1:nCi)
  REAL(dp):: Dum(1:nAq)
  REAL(dp):: QsK, QskIota
  INTEGER :: iMk,iActiv
  !------
  QskIota=    QsK_Iota
  iActiv=     0
  !
  vLKinActiv(:)= .FALSE.
  vLEquActiv(:)= .FALSE.
  !
  vKinQsk(:)= Zero
  vVmQsK(:)=  Zero
  vVmAct(:)=  Zero
  vKinPrm(:)= 0
  !
  != compute solute activities =================================================
  vLnAct(2:nAq)= vLnX(2:nAq) +vLnGam(2:nAq) -vLnX(isW) -LOG(MWSv)
  !---------------------------------------------------------------------
  !
  DO iMk=1,SIZE(vKinFas)
    !
    !-------------------------------------------------- calc. vVmQsK ---
    CALL KinRate_CalcQsK( & !
    & nCi,nAx+nMx,    & !
    & vDG_Kin(iMk),   & !
    & tNu_Kin(iMk,:), & !
    & vLnX,           & !IN:  Ln(aqu.species mole nrs) -> used for CalcQsK
    & vLnGam,         & !IN:  Ln(ActCoeff,aqu.species)
    & vLnBuf,         & !IN:  Ln(ActivitySpecies), used for buffered species
    & vLnAct(isW),    & !IN:  Ln(ActivitySolvent)
    & QsK,            & !OUT, mineral saturation
    & dQsKdLnXi)        !OUT
    !
    !!------------------------------------- cut extreme values of QsK ---
    !IF(iDebug>0 .AND. (QsK>Qsk_Max .OR. QsK<-Qsk_Min)) &
    !& WRITE(fTrc,'(A,G15.6,A)') vKinFas(iMk)%Name, QsK, "-> off Limits -> cut-off"
    !QsK= MIN(QsK,QsK_Max)
    !QsK= MAX(QsK,QsK_Min)
    !
    vKinQsK(iMk)= QsK
    !
    !------------------------------------- kinetic phases (%iKin/-0) ---
    IF(vKinFas(iMk)%iKin>0) THEN
      !------------------------------------- from QsK deduce cSatur, ---
      !------------------------------ for branching to dissol/precip ---
      CALL KinRate_SatState( &
      & vKinFas(iMk),   & !IN: M= mineral: for cSatur, QsKseuil, cMode
      & vMolK(iMk),     & !IN: nMol= Nr Moles min', to check nMol<=MolMinim
      & vKinMinim(iMk), & !IN: nMolMinim
      & QsK,            & !IN: QsK
      & QskIota,        & !IN
      & vStatusK(iMk))    !OUT: saturation state
      !---/
      !
      !---------------------------------------- from cSatur and QsK, ---
      !---------------------- calculate the QsK part of the rate law ---
      IF (vStatusK(iMk)=="D" .OR. vStatusK(iMk)=="P") &
      & CALL KinRate_CalcQsKFactor( & 
      & vStatusK(iMk),              & !IN
      & vKinMod(vKinFas(iMk)%iKin),& !IN
      & QsK,                        & !IN
      & dQsKdLnXi,                  & !IN
      & vVmQsK(iMk),                & !OUT: QsK function in the rate law
      & dum)
      !---/
      !
      vLKinActiv(iMk)= &
      & vStatusK(iMk)=="D" .OR. &  !"DISSOLU"
      & vStatusK(iMk)=="P"         !"PRECIPI"
      !
      IF(vLKinActiv(iMk)) THEN
        iActiv= iActiv + 1
        vKinPrm(iActiv)= iMk
      ENDIF
      !------------------------------------------------ calc. vVmAct ---
      IF(vLKinActiv(iMk)) THEN
        CALL KinRate_CalcActivFactor( & !
        & vStatusK(iMk),              & !IN
        & vKinMod(vKinFas(iMk)%iKin), & !IN: kinetic parameters
        & vLnAct,                     & !IN: activate / inhibate
        & vVmAct(iMk))                  !OUT
        !! & Dum) !dVmAdLnXf(iMk,:))       !OUT
      ENDIF
      !
      vKinFas(iMk)%Dat%cSat= vStatusK(iMk)
    ELSE
    !---------------------------------- non kinetic phases (%iKin-0) ---
      !-------------------------------------------------- vLEquActiv ---
      IF(vKinFas(iMk)%cMode=="E") THEN
        vLEquActiv(iMk)= & !
        & (vMolK(iMk)> vKinMinim(iMk) *Two) .OR. & ! minereal is present
        & (vKinQsK(iMk)> QsKIota +One)             ! or it has become saturated
      ENDIF
      !
      !! if(iDebug>2) write(51,'(2(G15.6,1X))',advance="no") vMolK(iMk),vKinMinim(iMk)
      !
    ENDIF !!IF(vKinFas(iMk)%iKin>0)
    !
  ENDDO
  !
  !! if(iDebug>2) write(51,*)
  !---------------------------------------------------------------------
  !
  !-------- compute fraction of box occupied by "non-active" minerals --
  PhiInert= Zero
  DO iMk=1,SIZE(vKinFas)
    IF(.NOT. ( vLKinActiv(iMk) .OR. vLEquActiv(iMk) )) &
    & PhiInert= PhiInert + vMolK(iMk) * vMolarVol(iMk) /VBox
  ENDDO
  !--------/
  !
  RETURN
ENDSUBROUTINE Dynam_OneStep_KinFas_Update

SUBROUTINE Dynam_OneStep_KinFasSurf_Update( & !
& vLKinActiv,  & ! IN
& vStatusK,    & ! IN
& vMolK0,vMolK,& ! IN
& PhiF0,PhiF,  & ! IN
& vSurfK0,     & ! IN
& vSurfK)        ! OUT
!--
!-- update mineral surface
!--
  !
  USE M_KinFas_Surf,ONLY: KinFas_Surf_Calc,KinFas_Surf_Crunch,KinFas_Surf_Update
  USE M_KinRate
  !
  USE M_Global_Vars,ONLY: vFas,vKinFas
  USE M_Dynam_Vars, ONLY: sModelSurf,Implicit_Surface
  !---------------------------------------------------------------------
  LOGICAL,  INTENT(IN) :: vLKinActiv(:)
  CHARACTER,INTENT(IN) :: vStatusK(:)
  REAL(dp), INTENT(IN) :: vMolK0(:)
  REAL(dp), INTENT(IN) :: vMolK(:)
  REAL(dp), INTENT(IN) :: PhiF0,PhiF
  REAL(dp), INTENT(IN) :: vSurfK0(:)
  REAL(dp), INTENT(OUT):: vSurfK(:)
  !
  REAL(dp):: x1, x2
  INTEGER :: iMk
  !---------------------------------------------------------------------
  
  !----------------------------------------------------- calc.surface --
  !----------------(IF Implicit_Surface, it's already done in FuncPsi)--
  DO iMk=1,SIZE(vKinFas)
    !
    IF(vLKinActiv(iMk)) THEN
      !
      SELECT CASE(sModelSurf)
      !
      CASE("CRUNCH")
        CALL KinFas_Surf_Crunch( &
        !."CRUNCH mode" currently implemented ONLY in explicit surface mode
        & vStatusK(iMk),            & !IN
        & vMolK(iMk),               & !IN: Nr Moles of mineral
        & vMolK0(iMk),              & !IN: Nr Moles at given stage 0
        & vSurfK0(iMk),             & !IN: Surf at stage 0
        & PhiF,PhiF0,               & !IN: PhiF current /PhiF at stage 0
        & vSurfK(iMk))                !OUT: Surface
      !
      CASE("SPHERE")
        CALL KinFas_Surf_Calc( & !
        & .false.,            & !IN: bImplicit
        & vKinFas(iMk)%cMode, & !IN: kinetic phase
        & vMolK(iMk),         & !IN: Nr Moles of mineral
        & vMolK0(iMk),        & !IN: Nr Moles of phase at given stage
        & vSurfK0(iMk),       & !IN: Surface of phase at that stage
        & vSurfK(iMk),        & !OUT: Surface
        & x1, x2) !dumm
      !
      END SELECT
      !
      vKinFas(iMk)%Dat%Surf= vSurfK(iMk)
      !
    ENDIF !!IF(vLKinActiv(iMk)) 
    !
  ENDDO
  !----------------------------------------------------/ calc.surface --
  
  DO iMk=1,SIZE(vKinFas)
    IF(vLKinActiv(iMk)) THEN
      != update specific surface and equivalent radius
      CALL KinFas_Surf_Update( & !
      & vFas,vMolK(iMk),vSurfK(iMk), &  !IN
      & vKinFas(iMk)) !INOUT
    ENDIF
  ENDDO
  
  RETURN
ENDSUBROUTINE Dynam_OneStep_KinFasSurf_Update

SUBROUTINE Dynam_OneStep_GammaUpd( &
& TdgK,Pbar, &
& vMolF,vLnAct,vLnGam, &
& Z_Plus,Z_Minus, &
& pH_)

  USE M_Numeric_Const,ONLY: Ln10
  USE M_Dtb_Const,   ONLY: T_CK 
  USE M_Global_Vars, ONLY: vSpc,SolModel
  USE M_SolModel_Calc,ONLY: Solmodel_CalcGamma
  USE M_Basis_Vars,  ONLY: isH_,isW,vLAx,vOrdAq
  !------
  REAL(dp),INTENT(IN)   :: TdgK,Pbar
  REAL(dp),INTENT(INOUT):: vMolF(:),vLnAct(:),vLnGam(:)
  REAL(dp),INTENT(OUT)  :: Z_Plus,Z_Minus
  REAL(dp),INTENT(OUT)  :: pH_
  !
  REAL(dp),DIMENSION(SIZE(vMolF)):: vMole,vLAct,vLGam
  LOGICAL, DIMENSION(SIZE(vMolF)):: vTooLow
  REAL(dp):: OsmoSv
  INTEGER :: N
  !------
  
  !--- permute aqu'species data to vSpc order
  vMole(vOrdAq(:))= vMolF(:)
  vLAct(vOrdAq(:))= vLnAct(:)
  vLGam(vOrdAq(:))= vLnGam(:)
  !---/
  !
  CALL Solmodel_CalcGamma( &
  & TdgK,Pbar,   &
  & SolModel,    & !IN
  & isW,vSpc,    & !IN
  & vLAx,        & !IN
  & vMole,       & !IN
  & vLAct,vLGam, & !INOUT
  & vTooLow,OsmoSv)
  !
  N= SIZE(vMole)
  Z_Plus= SUM(vMole(1:N)*vSpc(1:N)%Z, MASK=(vSpc(1:N)%Z >0)) /vMole(isW)
  Z_Minus=SUM(vMole(1:N)*vSpc(1:N)%Z, MASK=(vSpc(1:N)%Z <0)) /vMole(isW)
  !WRITE(16,'(2G15.6)') Z_Plus,Z_Minus
  !
  pH_=-vLAct(isH_)/Ln10
  !
  !--- back-permute aqu'species data
  vMolF(:)=  vMole(vOrdAq(:))
  vLnAct(:)= vLAct(vOrdAq(:))
  vLnGam(:)= vLGam(vOrdAq(:))
  !---/
  !
ENDSUBROUTINE Dynam_OneStep_GammaUpd

SUBROUTINE Dynam_OneStep_DLGammaUpd( &
& vLnGam, &
& vDLGam_As)
  USE M_Basis_Vars,ONLY:  nCi,nAs,tNuAs
  !
  REAL(dp),INTENT(IN) :: vLnGam(:)
  REAL(dp),INTENT(OUT):: vDLGam_As(:)
  !
  INTEGER:: iAs

  DO iAs=1,nAs
    vDLGam_As(iAs)= vLnGam(nCi+iAs) &
    &             - DOT_PRODUCT(tNuAs(iAs,2:nCi), vLnGam(2:nCi))
  ENDDO

  RETURN
ENDSUBROUTINE Dynam_OneStep_DLGammaUpd

!-----------------------------------------------------------------------
!---adjust time step----------------------------------------------------
!-----------------------------------------------------------------------
SUBROUTINE Dynam_CtrlDtime( &
& NewtIterMax,NewtIterMin, &
& Time_Decrease,Time_Increase, &
& iCtrl,iter, &
& Time,TFinal,dTMax,&
& VarAqu_,VarMin_,&
& dTime_)
  INTEGER, INTENT(IN)   :: NewtIterMax,NewtIterMin
  REAL(dp),INTENT(IN)   :: Time_Decrease,Time_Increase
  INTEGER, INTENT(IN)   :: iCtrl,iter
  REAL(dp),INTENT(IN)   :: Time,TFinal,dTMax,VarAqu_,VarMin_
  REAL(dp),INTENT(INOUT):: dTime_ 
  !
  REAL(dp)::VarMax_,dT0
  !
  dT0=dTime_
  !
  IF(dT0<dTMax) THEN
    SELECT CASE(iCtrl)
    !
    CASE(1)
      !----------------- when convergence is slow, decrease time step --
      IF(iter>NewtIterMax) dT0= dT0 *Time_Decrease
      !
      !---------------- when convergence is quick, increase time step --
      IF(iter<NewtIterMin) dT0= dT0 *Time_Increase
      !
    CASE(2)
      VarMax_=MAX(VarAqu_,VarMin_)
      IF(VarMax_<1.0E-3) dT0=dT0*2.0D0 !LOG(1.030D0)
      IF(VarMax_>1.0E-2) dT0=dT0/2.0D0 !LOG(1.080D0)
      !
    END SELECT
    !
  ENDIF
  !
  IF(dT0>dTMax) dT0=dTMax
  !
  dTime_=dT0
  IF (Time +dTime_ > tFinal) THEN !!!??? dTime ???
    dTime_=tFinal -Time ! go directly to the solution now
  ELSE IF(Time +2*dTime_ > tFinal) THEN
    dTime_= 0.5*(tFinal -Time) ! slow DOwn two steps before the wall
  END IF
  !
ENDSUBROUTINE Dynam_CtrlDtime

SUBROUTINE ShoMin(f,iStep_,VBox,TimeScaled,pH_,PhiF,vMolK,vKinMVol,vKinQsk)
  USE M_Numeric_Const, ONLY: Ln10, TinyDP
  !
  INTEGER, INTENT(IN):: f,iStep_
  REAL(dp),INTENT(IN):: VBox,TimeScaled,pH_,PhiF
  REAL(dp),DIMENSION(:),INTENT(IN):: vMolK,vKinMVol,vKinQsk
  !
  INTEGER::I,N
  !
  N=SIZE(vMolK)
  WRITE(f,'(I7,A1,3(G15.8,A1),$)') iStep_,T_,TimeScaled,T_,pH_,T_,PhiF,T_
  DO i=1,N; WRITE(f,'(G15.8,A1,$)') vMolK(i)*vKinMVol(i)/VBox,T_ ; ENDDO
  DO i=1,N; WRITE(f,'(G15.8,A1,$)') LOG(vKinQsK(i))/Ln10,T_      ; ENDDO
  DO i=1,N; WRITE(f,'(G15.8,A1,$)') vMolK(i),T_                  ; ENDDO
  WRITE(f,*)
  !
ENDSUBROUTINE ShoMin

SUBROUTINE ShoMinRate( & !
!--
!-- write details of mineral rates (surface, vVmQsK, vVmAct, ...)
!--
& f,iStep_,VBox,TimeScaled,pH_,PhiF, & !
& vMolK,vKinMVol,vKinQsk,vSurfK,vSurfK0,vVmQsK,vVmAct, & !
& vKinFas)

  USE M_Numeric_Const,ONLY: Ln10 !, TinyDP
  USE M_T_KinFas,ONLY: T_KinFas
  !------
  INTEGER,                    INTENT(IN):: f,iStep_
  REAL(dp),                   INTENT(IN):: VBox,TimeScaled,pH_,PhiF
  REAL(dp),      DIMENSION(:),INTENT(IN):: vMolK,vKinMVol,vKinQsk,vSurfK,vSurfK0
  REAL(dp),      DIMENSION(:),INTENT(IN):: vVmQsK,vVmAct
  TYPE(T_KinFas),DIMENSION(:),INTENT(IN):: vKinFas
  !
  INTEGER::I,N
  !------
  
  N=SIZE(vMolK)
  WRITE(f,'(I7,A1,3(G15.8,A1),$)') iStep_,T_,TimeScaled,T_,pH_,T_,PhiF,T_
  !
  DO i=1,n; WRITE(f,'(G15.8,A1,$)') vMolK(i)*vKinMVol(i)/VBox,  T_ ; ENDDO
  DO i=1,n; WRITE(f,'(G15.8,A1,$)') vSurfK(i)/vSurfK0(i),       T_ ; ENDDO
  DO i=1,n; WRITE(f,'(G15.8,A1,$)') vKinFas(i)%Dat%SurfKg/1.D3, T_ ; ENDDO
  DO i=1,n; WRITE(f,'(G15.8,A1,$)') LOG(vKinQsK(i))/Ln10,       T_ ; ENDDO
  DO i=1,n; WRITE(f,'(G15.8,A1,$)') vVmQsK(i),                  T_; ENDDO
  DO i=1,n; WRITE(f,'(G15.8,A1,$)') vVmAct(i),                  T_; ENDDO
  DO i=1,n; WRITE(f,'(G15.8,A1,$)') LOG(vKinFas(i)%Dat%Radius)/Ln10,T_; ENDDO
  !
  WRITE(f,*)
  
  RETURN
ENDSUBROUTINE ShoMinRate

ENDMODULE M_Dynam_OneStep_Tools

