module M_Dynam_OneStep_Tools
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_
  implicit none

  private
  
  public:: Dynam_OneStep_Write
  public:: Dynam_OneStep_WriteElements
  public:: Dynam_OneStep_TotF
  public:: Dynam_OneStep_TotK
  public:: Dynam_OneStep_Adjust
  public:: Dynam_OneStep_KinFas_Update
  public:: Dynam_OneStep_KinFasSurf_Update
  public:: Dynam_OneStep_SecondSpecies
  public:: Dynam_OneStep_GammaUpd
  public:: Dynam_OneStep_DLGammaUpd
  public:: Dynam_OneStep_Variations
  public:: Dynam_OneStep_Write_Bilan
  public:: Dynam_CtrlDtime

contains

subroutine Dynam_OneStep_SecondSpecies(N,vX,vXsec)
  use M_Global_Vars,only: SolModel
  use M_Basis_Vars, only: nAs,nAx,nMx,tNuAs,nCi
  use M_Dynam_Vars, only: vDG_As,vDLGam_As,vLnGam,vLnAct,LogForAqu,vLnBuf
  !
  integer, intent(in) :: N
  real(dp),intent(in) :: vX(:)
  real(dp),intent(out):: vXsec(:)
  !
  integer :: iAs,iCi,J
  real(dp):: X,MWSv
  !
  MWSv= SolModel%MolWeitSv
  !
  if(LogForAqu) then
    !
    do iAs=1,nAs
      !
      X=   vLnAct(1) *tNuAs(iAs,1) &
      &  + (One - sum(tNuAs(iAs,2:nCi)))*(vX(1) +log(MWSv)) &
      &  - vDG_As(iAs) -vDLGam_As(iAs)
      do iCi=1,N
        if(iCi/=1) X= X + vX(iCi) *tNuAs(iAs,iCi)
      end do
      do J=1,nAx+nMx
        X=  X + vLnBuf(J) *tNuAs(iAs,nCi+J)
      end do
      !
      vXsec(iAs)= exp(X)
      != (mole nr species iAs)
      !
    end do
    !
  else
    !
    do iAs=1,nAs
      !
      X=  vLnAct(1) *tNuAs(iAs,1) &
      & - vDG_As(iAs) -vDLGam_As(iAs)
      do J=1,nAx+nMx
        X= X + vLnBuf(J) *tNuAs(iAs,nCi+J)
      end do
      X= exp(X)
      !
      do iCi=1,N
        if(iCi/=1 .and. tNuAs(iAs,iCi) /= Zero) &
        & X= X *(vX(iCi))**tNuAs(iAs,iCi)
      end do
      vXsec(iAs)= X *(vX(1)*MWSv)**(One-sum(tNuAs(iAs,2:nCi)))
      != (mole nr species iAs)
      !
    end do
    !
  end if
    !
  !---------------------------------------------------------- trace --
  if(iDebug>2) then
    if(LogForAqu) then
      write(71,'(A)') "====================SOLVE=="
      do j=1,N
        write(71,'(A,G15.6)') "resPrim=",exp(vX(j))
      end do
      do j=1,size(vXsec)
        write(71,'(A,G15.6)') "resSec==",exp(vXsec(j))
      end do
    else
      write(71,'(A)') "====================SOLVE=="
      do j=1,N
        write(71,'(A,G15.6)') "resPrim=",vX(j)
      end do
      do j=1,size(vXsec)
        write(71,'(A,G15.6)') "resSec==",vXsec(j)
      end do
    end if
  end if
  !---------------------------------------------------------/ trace --
  !pause
  !
end subroutine Dynam_OneStep_SecondSpecies

subroutine Dynam_OneStep_Write( &
& iStep,Time,dTime,       &
& fSavRate,fSavTime,      &
& PhiF,RhoF,VBox,dX,      &
& vTotF,vMolF,vLnAct,vLnGam, &
& vTotK,vMolK,vKinMVol,vKinQsk,vSurfK,vSurfK0,vVmQsK,vVmAct, &
& dTSav,TSave)
  !------
  use M_Numeric_Const,only: Ln10
  use M_IoTools,      only: OutStrVec
  use M_T_Species,    only: Species_Index
  use M_Safe_Functions
  !
  use M_Dynam_Files,only: fDynEle,fDynMol,fDynAct,fDynGam,fDynMnK
  !
  use M_Global_Vars,only: vKinFas,SolModel,vSpc
  use M_Dynam_Vars, only: vWeitCp,TimeScale
  !------
  integer,              intent(in):: iStep
  real(dp),             intent(in):: Time,dTime
  integer,              intent(in):: fSavRate,fSavTime
  real(dp),             intent(in):: PhiF,RhoF,VBox,dX
  real(dp),dimension(:),intent(in):: vTotF,vMolF,vLnAct,vLnGam
  real(dp),dimension(:),intent(in):: vTotK,vMolK,vKinMVol,vKinQsk
  real(dp),dimension(:),intent(in):: vSurfK,vSurfK0,vVmQsK,vVmAct
  real(dp),             intent(in):: dTSav
  !
  real(dp),             intent(inout):: TSave
  !
  integer :: isW,nAq,I,isH_
  real(dp):: R,UDelta,VolF,MWsv,pH_
  real(dp):: vX(size(vLnGam))
  !------
  isW=  SolModel%iSolvent
  MWSv= SolModel%MolWeitSv
  isH_= SolModel%isH_
  nAq=  size(vMolF)
  !
  VolF= dot_product(vWeitCp(:),vTotF(:))  /RhoF  !volume of fluid in box
  UDelta= (VolF -PhiF) /dTime /(VBox/dX)
  R= vTotF(isW)*MWsv
  !
  !-------------- write mole nrs of elements in fluid and in minerals --
  if(fDynEle>0) then
  write(fDynEle,'(I7,A1,2(G15.6,A1),$)') iStep,T_,Time/TimeScale,T_,UDelta,T_
  call OutStrVec(fDynEle, vTotF(:),      Opt_C="G", CR=.false.) !->mole nrs in fluid
  call OutStrVec(fDynEle, vTotF(:)/R,    Opt_C="G", CR=.false.) !->molality
  call OutStrVec(fDynEle, vTotK(:),      Opt_C="G")
  end if
  !
  !-- nr'mole component injected in time step:
  !-- - vTotInj(iPr) *UDarcy /dX *dTime
  !
  !--------------------------------------- write mole nrs aqu'species --
  if(fDynMol>0) &
  & call OutStrVec(fDynMol, vMolF(1:nAq), Opt_C="G", Opt_I=iStep, Opt_R=Time/TimeScale)
  !
  !------------------------------------ write activities aqu'species ---
  if(fDynAct>0) &
  & call OutStrVec(fDynAct, vLnAct(1:nAq)/Ln10,Opt_C="F", Opt_I=iStep, Opt_R=Time/TimeScale)
  !
  !----------------------------------- write activ'coeffs aqu'species --
  do I=1,nAq
    vX(i)= FSafe_Exp(vLnGam(i))
  end do
  if(fDynGam>0) &
  & call OutStrVec(fDynGam, vX(1:nAq),Opt_C="F", Opt_I=iStep, Opt_R=Time/TimeScale)
  !
  !------------------ write fluid and min'species fractions, QsK, etc --
  pH_= -vLnAct(isH_)/Ln10
  if(fDynMnK>0) &
  & call ShoMin(fDynMnK,iStep,VBox,Time/TimeScale,pH_,PhiF,vMolK,vKinMVol,vKinQsk)
  !
  if(fSavRate>0) &
  & call ShoMinRate( &
  & fSavRate,iStep,VBox,Time/TimeScale,pH_,PhiF, &
  & vMolK,vKinMVol,vKinQsk, &
  & vSurfK,vSurfK0,vVmQsK,vVmAct, &
  & vKinFas)
  !
  if(dTSav>Zero .and. Time>=TSave) then
    TSave= TSave +dTSav
    call ShoMin(fSavTime,iStep,VBox,Time/TimeScale,pH_,PhiF,vMolK,vKinMVol,vKinQsk)
  end if
  !
  return
end subroutine Dynam_OneStep_Write

subroutine Dynam_OneStep_WriteElements( &
& iStep,Time,dTime, &
& vMolF,vMolK)
  use M_Dynam_Files,only: fDynElement
  use M_Dynam_Vars, only: tStoikioKin,tStoikioAqu,TimeScale
  !
  integer, intent(in):: iStep
  real(dp),intent(in):: Time,dTime
  real(dp),intent(in):: vMolF(:)
  real(dp),intent(in):: vMolK(:)
  !
  real(dp),allocatable:: vX(:)
  integer :: iEl,iSp,N
  real(dp):: X
  !
  N= size(tStoikioAqu,2) != number of elements
  allocate(vX(N))
  !
  write(fDynElement,'(I7,A1,G15.6,A1,$)') &
  & iStep,T_,Time/TimeScale,T_
  !
  do iEl=1,N
    X= Zero
    do iSp=1,size(vMolF)
      X= X + tStoikioAqu(iSp,iEl)*vMolF(iSp)
    end do
    vX(iEl)= X
    write(fDynElement,'(G15.6,A1)',advance="NO") X,T_
  end do
  !
  do iEl=1,N
    X= Zero
    do iSp=1,size(vMolK)
      X= X + tStoikioKin(iSp,iEl)*vMolK(iSp)
    end do
    vX(iEl)= vX(iEl)+X
    write(fDynElement,'(G15.6,A1)',advance="NO") X,T_
  end do
  !
  do iEl=1,N
    write(fDynElement,'(G15.6,A1)',advance="NO") vX(iEl),T_
  end do
  !
  write(fDynElement,*)
  !
  deallocate(vX)
  
  return
end subroutine Dynam_OneStep_WriteElements

subroutine Dynam_OneStep_Write_Bilan( &
& F,iStep,Time,dTime,  &
& PhiF,UDarcy,VBox,dX, &
& vTotInj,vTotF_0,vTotF,vTotK_0,vTotK)
!--
!-- write balances on components
!--
  use M_IoTools
  !
  integer, intent(in) :: F,iStep
  real(dp),intent(in) :: Time,dTime,PhiF,UDarcy,VBox,dX
  real(dp),intent(in) :: vTotInj(:),vTotF_0(:),vTotF(:)
  !                      !mole nrs in fluid, prev. & current step
  real(dp),intent(in) :: vTotK_0(:),vTotK(:)
  !                      !mole nrs in minerals, prev. & current step
  !
  real(dp),allocatable:: vX(:)
  !
  allocate(vX(size(vTotF)))
  !
  write(F,'(I7,A1,G15.6,A1,$)') iStep,T_,Time,T_
  !
  !
  vX(:)= vTotInj(:) *UDarcy /dX
  call OutStrVec(F, vX(:), Opt_C="G", CR=.false.)
  !
  vX(:)=  vTotF(:)   + vTotK(:)   &
  &    - (vTotF_0(:) + vTotK_0(:))
  vX(:)= vX(:) /dTime
  call OutStrVec(F, vX(:), Opt_C="G", CR=.false.)
  !
  vX(:)= vX(:) -vTotInj(:) *UDarcy /dX
  call OutStrVec(F, vX(:), Opt_C="G", CR=.false.)
  !vX(:)=  vTotF(:) -vTotF_0(:) &
  !&    +  vTotK(:) -vTotK_0(:) &
  !&    + (vTotF(:) /PhiF  - vTotInj(:) *dTime *UDarcy /dX)
  !
  write(F,*)
  !
  deallocate(vX)
  !
  return
end subroutine Dynam_OneStep_Write_Bilan

subroutine Dynam_OneStep_Variations( &
& F,iStep,Time,dTime, &
& vMolF_0,vMolK_0,vMolF,vMolK,&
& MaxAqu,MaxMin)
!--
!-- compute max.variations during time step
!--
  use M_Numeric_Tools,only: iMaxLoc_R
  !
  use M_Global_Vars,only: vSpc,vKinFas
  use M_Basis_Vars, only: vOrdAq
  !
  integer, intent(in) :: F,iStep
  real(dp),intent(in) :: Time,dTime
  real(dp),intent(in) :: vMolF_0(:),vMolK_0(:)
  real(dp),intent(in) :: vMolF(:),vMolK(:)
  real(dp),intent(out):: MaxAqu,MaxMin
  !
  integer :: iMaxAqu,iMaxMin
  real(dp),allocatable:: vDeltaF(:),vDeltaK(:)
  !
  !do i=1,size(vMolF)
  !  print '(3(G15.6,1X))',vMolF_0(i),vMolF(i),abs(vMolF(i) -vMolF_0(i))
  !end do
  !pause
  !
  !do i=1,size(vMolK)
  !  print '(3(G15.6,1X))',vMolK_0(i),vMolK(i),abs(vMolK(i) -vMolK_0(i))
  !end do
  !pause
  !
  !print *,"MAXVAL",  maxval(abs(vMolF(:) -vMolF_0(:)))     ;   pause
  !
  allocate(vDeltaF(size(vMolF)))
  if(size(vKinFas)>0) allocate(vDeltaK(size(vMolK)))
  !
  vDeltaF(:)= abs(vMolF(:) -vMolF_0(:))
  if(size(vKinFas)>0) vDeltaK(:)= abs(vMolK(:) -vMolK_0(:))
  !
  iMaxAqu= iMaxLoc_R(vDeltaF(:))
  !-> which aqu.species varies most
  if(size(vKinFas)>0) iMaxMin= iMaxLoc_R(vDeltaK(:))
  !-> which min.species varies most
  MaxAqu= maxval(vDeltaF(:))
  if(size(vKinFas)>0) MaxMin=  maxval(vDeltaK(:))
  !
  MaxAqu= MaxAqu /vMolF(iMaxAqu)
  if(size(vKinFas)>0) MaxMin= MaxMin /vMolK(iMaxMin)
  !
  if(F>0) then
  
    if(size(vKinFas)>0) then
      
      write(F,'(2(A,A1),I6,A1,2(G15.6,A1),4(G15.6,A1))') &
      & trim(vSpc(vOrdAq(iMaxAqu))%NamSp),T_, &
      & trim(vKinFas(iMaxMin)%NamKF),     T_, &
      & iStep,T_,           &
      & Time,T_,  dTime,T_, &
      & MaxAqu,T_,MaxMin,T_,&
      & MaxAqu /vMolF(iMaxAqu),T_,MaxMin /vMolK(iMaxMin),T_
      !& maxval(vDeltaMolF),T_,maxval(vDeltaMolM),T_
      
    else
      
      write(F,'(A,A1,I3,A1,2(G15.6,A1),2(G15.6,A1))') &
      & trim(vSpc(vOrdAq(iMaxAqu))%NamSp),T_, &
      & iStep,T_, &
      & Time,T_,  dTime,T_, &
      & MaxAqu,T_, &
      & MaxAqu /vMolF(iMaxAqu),T_
      !& maxval(vDeltaMolF),T_,maxval(vDeltaMolM),T_
    
    end if
  
  end if
  !
  deallocate(vDeltaF)
  if(size(vKinFas)>0) deallocate(vDeltaK)
  !
end subroutine Dynam_OneStep_Variations

subroutine Dynam_OneStep_TotF( & !
& vMolF,  & !in
& vTotF)    !out
!--
!-- update mole nrs COMPONENTS IN FLUID at time T
!-- will give vTotF, result of previous step used as RHS in FuncPsi(1:nCi)
!--
  use M_Basis_Vars,only: nCi,tAlfPr,tAlfAs,nAs
  !
  real(dp),intent(in) :: vMolF(:)
  real(dp),intent(out):: vTotF(:)
  !
  integer :: nCp
  !
  nCp= size(vTotF)
  !
  !update mole nrs ALL COMPONENTS in Fluid at time T (= result of previous step)
  vTotF(1:nCp)= matmul(tAlfPr(1:nCp,1:nCi),vMolF(1:    nCi))     &
  &           + matmul(tAlfAs(1:nCp,1:nAs),vMolF(nCi+1:nCi+nAs))
  !
end subroutine Dynam_OneStep_TotF

subroutine Dynam_OneStep_TotK( & !
& vMolK,  & !in
& vTotK)    !out
!--
!-- update mole nrs COMPONENTS IN MINERALS at time T
!-- vTotK is used for balance computations
!--
  use M_Dynam_Vars, only: tAlfKin
  !------
  real(dp),intent(in) :: vMolK(:)
  real(dp),intent(out):: vTotK(:)
  !------
  integer :: I
  !
  do I=1,size(vTotK)
    vTotK(I)= dot_product(tAlfKin(I,:),vMolK(:))
  end do
  !
end subroutine Dynam_OneStep_TotK

subroutine Dynam_OneStep_Adjust( & !
& RhoF,VBox,PhiF,  & !in
& vMolF,           & !inout (if UpdateMassFluid then adjust)
& vTotF)             !inout
!--
!-- if volume is "not free",
!-- then adjust mole numbers to available volume, based on density
!-- will give vTotF, result of previous step used as RHS in FuncPsi(1:nCi)
!--
  use M_Global_Vars,only: nAq !,nMk
  use M_Dynam_Vars, only: vWeitSp
  !------
  real(dp),intent(in)::    RhoF,VBox,PhiF
  real(dp),intent(inout):: vMolF(:)
  real(dp),intent(inout):: vTotF(:)
  !
  real(dp):: R
  !------
  ! vWeitSp(:)=vSpc(vOrdAq(:))%WeitKg
  ! -> Aqu'Species' weitkg sorted according to local Species' order
  !
  ! MFluid= mass of fluid at time T
  !       = dot_product(vWeitSp(1:nAq),vMolF(1:nAq))
  ! MFluidBox= mass of fluid of density RhoF that fit in VBox*PhiF
  !          = RhoF*VBox*PhiF !=mass fluid in box of volume VBox
  !
  ! -> conversion factor
  ! R= MFluidBox / MFluid
  !  = RhoF*VBox*PhiF /dot_product(vWeitSp(1:nAq),vMolF(1:nAq))
  !
  R= RhoF*VBox*PhiF /dot_product(vWeitSp(:),vMolF(:))
  vTotF(:)= vTotF(:) *R !nr moles COMPONENTS in fluid in box
  vMolF(:)= vMolF(:) *R !nr moles SPECIES in fluid in box
  !
end subroutine Dynam_OneStep_Adjust

subroutine Dynam_OneStep_KinFas_Update( & !
& vLnX,                   & !in
& vMolK,                  & !in
& vKinQsk,vStatusK,       & !out
& vVmQsK,vVmAct,          & !out
& vLKinActiv,vLEquActiv,  & !out
& vKinPrm,PhiInert)         !out
!--
!-- compute mineral Qsk, rate
!--
  use M_Numeric_Const,only: Ln10
  use M_KinRate
  use M_T_KinFas
  !
  use M_Global_Vars,only: SolModel,vFas,vKinFas,nAq
  use M_Basis_Vars, only: nCi,nMx,nAx
  !
  use M_Dynam_Vars, only: &
  & sModelSurf,Implicit_Surface, &
  & VBox,vDG_Kin,tNu_Kin, &
  & vMolarVol,vKinMod,    &
  & vKinMinim,            &
  & vLnAct,vLnGam,vLnBuf, &
  & vMolK0,vSurfK0,PhiF0, &
  & QsK_Iota
  !------
  real(dp),intent(in):: vLnX(:)
  real(dp),intent(in):: vMolK(:)
  !
  real(dp), intent(out):: vKinQsk(:)
  real(dp), intent(out):: vVmQsK(:)
  real(dp), intent(out):: vVmAct(:)
  logical,  intent(out):: vLKinActiv(:)
  logical,  intent(out):: vLEquActiv(:)
  integer,  intent(out):: vKinPrm(:)
  real(dp), intent(out):: PhiInert
  character,intent(out):: vStatusK(:)
  !
  real(dp):: dQsKdLnXi(1:nCi)
  real(dp):: Dum(1:nAq)
  real(dp):: QsK, QskIota, MWSv
  integer :: iMk,iActiv,isW
  !------
  isW=  SolModel%iSolvent
  MWSv= SolModel%MolWeitSv
  !
  QskIota=    QsK_Iota
  iActiv=     0
  !
  vLKinActiv(:)= .false.
  vLEquActiv(:)= .false.
  !
  vKinQsk(:)= Zero
  vVmQsK(:)=  Zero
  vVmAct(:)=  Zero
  vKinPrm(:)= 0
  !
  != compute solute activities =================================================
  vLnAct(2:nAq)= vLnX(2:nAq) +vLnGam(2:nAq) -vLnX(isW) -log(MWSv)
  !---------------------------------------------------------------------
  !
  do iMk=1,size(vKinFas)
    !
    !-------------------------------------------------- calc. vVmQsK ---
    call KinRate_CalcQsK( & !
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
    !if(idebug>1 .and. (QsK>Qsk_Max .or. QsK<-Qsk_Min)) &
    !& write(fTrc,'(A,G15.6,A)') vKinFas(iMk)%Name, QsK, "-> off Limits -> cut-off"
    !QsK= min(QsK,QsK_Max)
    !QsK= max(QsK,QsK_Min)
    !
    vKinQsK(iMk)= QsK
    !
    !------------------------------------- kinetic phases (%iKin/-0) ---
    if(vKinFas(iMk)%iKin>0) then
      !------------------------------------- from QsK deduce cSatur, ---
      !------------------------------ for branching to dissol/precip ---
      call KinRate_SatState( &
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
      if (vStatusK(iMk)=="D" .or. vStatusK(iMk)=="P") &
      & call KinRate_CalcQsKFactor( & 
      & vStatusK(iMk),              & !IN
      & vKinMod(vKinFas(iMk)%iKin),& !IN
      & QsK,                        & !IN
      & dQsKdLnXi,                  & !IN
      & vVmQsK(iMk),                & !OUT: QsK function in the rate law
      & dum)
      !---/
      !
      vLKinActiv(iMk)= &
      & vStatusK(iMk)=="D" .or. &  !"DISSOLU"
      & vStatusK(iMk)=="P"         !"PRECIPI"
      !
      if(vLKinActiv(iMk)) then
        iActiv= iActiv + 1
        vKinPrm(iActiv)= iMk
      end if
      !------------------------------------------------ calc. vVmAct ---
      if(vLKinActiv(iMk)) then
        call KinRate_CalcActivFactor( & !
        & vStatusK(iMk),              & !IN
        & vKinMod(vKinFas(iMk)%iKin), & !IN: kinetic parameters
        & vLnAct,                     & !IN: activate / inhibate
        & vVmAct(iMk))                  !OUT
        !! & Dum) !dVmAdLnXf(iMk,:))       !OUT
      end if
      !
      vKinFas(iMk)%Dat%cSat= vStatusK(iMk)
    else
    !---------------------------------- non kinetic phases (%iKin-0) ---
      !-------------------------------------------------- vLEquActiv ---
      if(vKinFas(iMk)%cMode=="E") then
        vLEquActiv(iMk)= & !
        & (vMolK(iMk)> vKinMinim(iMk) *Two) .or. & ! minereal is present
        & (vKinQsK(iMk)> QsKIota +One)             ! or it has become saturated
      end if
      !
      !! if(iDebug>2) write(51,'(2(G15.6,1X))',advance="no") vMolK(iMk),vKinMinim(iMk)
      !
    end if !!if(vKinFas(iMk)%iKin>0)
    !
  end do
  !
  !! if(iDebug>2) write(51,*)
  !---------------------------------------------------------------------
  !
  !-------- compute fraction of box occupied by "non-active" minerals --
  PhiInert= Zero
  do iMk=1,size(vKinFas)
    if(.not. ( vLKinActiv(iMk) .or. vLEquActiv(iMk) )) &
    & PhiInert= PhiInert + vMolK(iMk) * vMolarVol(iMk) /VBox
  end do
  !--------/
  !
  return
end subroutine Dynam_OneStep_KinFas_Update

subroutine Dynam_OneStep_KinFasSurf_Update( & !
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
  use M_KinFas_Surf,only: KinFas_Surf_Calc,KinFas_Surf_Crunch,KinFas_Surf_Update
  use M_KinRate
  !
  use M_Global_Vars,only: vFas,vKinFas
  use M_Dynam_Vars, only: sModelSurf
  !---------------------------------------------------------------------
  logical,  intent(in) :: vLKinActiv(:)
  character,intent(in) :: vStatusK(:)
  real(dp), intent(in) :: vMolK0(:)
  real(dp), intent(in) :: vMolK(:)
  real(dp), intent(in) :: PhiF0,PhiF
  real(dp), intent(in) :: vSurfK0(:)
  real(dp), intent(out):: vSurfK(:)
  !
  real(dp):: x1, x2
  integer :: iMk
  !---------------------------------------------------------------------
  
  !----------------------------------------------------- calc.surface --
  !----------------(if Implicit_Surface, it's already done in FuncPsi)--
  do iMk=1,size(vKinFas)
    !
    if(vLKinActiv(iMk)) then
      !
      select case(sModelSurf)
      !
      case("CRUNCH")
        call KinFas_Surf_Crunch( &
        !."CRUNCH mode" currently implemented only in explicit surface mode
        & vStatusK(iMk),            & !IN
        & vMolK(iMk),               & !IN: Nr Moles of mineral
        & vMolK0(iMk),              & !IN: Nr Moles at given stage 0
        & vSurfK0(iMk),             & !IN: Surf at stage 0
        & PhiF,PhiF0,               & !IN: PhiF current /PhiF at stage 0
        & vSurfK(iMk))                !OUT: Surface
      !
      case("SPHERE")
        call KinFas_Surf_Calc( & !
        & .false.,            & !IN: bImplicit
        & vKinFas(iMk)%cMode, & !IN: kinetic phase
        & vMolK(iMk),         & !IN: Nr Moles of mineral
        & vMolK0(iMk),        & !IN: Nr Moles of phase at given stage
        & vSurfK0(iMk),       & !IN: Surface of phase at that stage
        & vSurfK(iMk),        & !OUT: Surface
        & x1, x2) !dumm
      !
      end select
      !
      vKinFas(iMk)%Dat%Surf= vSurfK(iMk)
      !
    end if !!if(vLKinActiv(iMk)) 
    !
  end do
  !----------------------------------------------------/ calc.surface --
  
  do iMk=1,size(vKinFas)
    if(vLKinActiv(iMk)) then
      != update specific surface and equivalent radius
      call KinFas_Surf_Update( & !
      & vFas,vMolK(iMk),vSurfK(iMk), &  !IN
      & vKinFas(iMk)) !INOUT
    end if
  end do
  
  return
end subroutine Dynam_OneStep_KinFasSurf_Update

subroutine Dynam_OneStep_GammaUpd( &
& TdgK,Pbar, &
& vMolF,vLnAct,vLnGam, &
& Z_Plus,Z_Minus)

  use M_Numeric_Const,only: Ln10
  use M_Dtb_Const,    only: T_CK 
  use M_Global_Vars,  only: vSpc,SolModel
  use M_SolModel_Calc,only: Solmodel_CalcGamma
  use M_Basis_Vars,   only: vLAx,vOrdAq
  !------
  real(dp),intent(in)   :: TdgK,Pbar
  real(dp),intent(inout):: vMolF(:),vLnAct(:),vLnGam(:)
  real(dp),intent(out)  :: Z_Plus,Z_Minus
  !
  real(dp),dimension(size(vMolF)):: vMole,vLAct,vLGam
  logical, dimension(size(vMolF)):: vTooLow
  real(dp):: OsmoSv
  integer :: N
  integer :: isW
  !------
  isW= SolModel%iSolvent
  !-------------------------------permute aqu'species data to vSpc order
  vMole(vOrdAq(:))= vMolF(:)
  vLAct(vOrdAq(:))= vLnAct(:)
  vLGam(vOrdAq(:))= vLnGam(:)
  !--------------------------------------------------------------------/
  !
  call Solmodel_CalcGamma( &
  & TdgK,Pbar,   &
  & SolModel,    & !IN
  & vSpc,        & !IN
  & vLAx,        & !IN
  & vMole,       & !IN
  & vLAct,vLGam, & !INOUT
  & vTooLow,OsmoSv)
  !
  N= size(vMole)
  Z_Plus= sum(vMole(1:N)*vSpc(1:N)%Z, mask=(vSpc(1:N)%Z >0)) /vMole(isW)
  Z_Minus=sum(vMole(1:N)*vSpc(1:N)%Z, mask=(vSpc(1:N)%Z <0)) /vMole(isW)
  !write(16,'(2G15.6)') Z_Plus,Z_Minus
  !
  ! pH_=-vLAct(isH_)/Ln10
  !
  !----------------------------------------back-permute aqu'species data
  vMolF(:)=  vMole(vOrdAq(:))
  vLnAct(:)= vLAct(vOrdAq(:))
  vLnGam(:)= vLGam(vOrdAq(:))
  !--------------------------------------------------------------------/
  !
end subroutine Dynam_OneStep_GammaUpd

subroutine Dynam_OneStep_DLGammaUpd( &
& vLnGam, &
& vDLGam_As)
  use M_Basis_Vars,only:  nCi,nAs,tNuAs
  !
  real(dp),intent(in) :: vLnGam(:)
  real(dp),intent(out):: vDLGam_As(:)
  !
  integer:: iAs

  do iAs=1,nAs
    vDLGam_As(iAs)= vLnGam(nCi+iAs) &
    &             - dot_product(tNuAs(iAs,2:nCi), vLnGam(2:nCi))
  end do

  return
end subroutine Dynam_OneStep_DLGammaUpd

!-----------------------------------------------------------------------
!---adjust time step----------------------------------------------------
!-----------------------------------------------------------------------
subroutine Dynam_CtrlDtime( &
& NewtIterMax,NewtIterMin, &
& Time_Decrease,Time_Increase, &
& iCtrl,iter, &
& Time,TFinal,dTMax,&
& VarAqu_,VarMin_,&
& dTime_)
  integer, intent(in)   :: NewtIterMax,NewtIterMin
  real(dp),intent(in)   :: Time_Decrease,Time_Increase
  integer, intent(in)   :: iCtrl,iter
  real(dp),intent(in)   :: Time,TFinal,dTMax,VarAqu_,VarMin_
  real(dp),intent(inout):: dTime_ 
  !
  real(dp)::VarMax_,dT0
  !
  dT0=dTime_
  !
  if(dT0<dTMax) then
    select case(iCtrl)
    !
    case(1)
      !----------------- when convergence is slow, decrease time step --
      if(iter>NewtIterMax) dT0= dT0 *Time_Decrease
      !
      !---------------- when convergence is quick, increase time step --
      if(iter<NewtIterMin) dT0= dT0 *Time_Increase
      !
    case(2)
      VarMax_=MAX(VarAqu_,VarMin_)
      if(VarMax_<1.0E-3) dT0=dT0*2.0D0 !log(1.030D0)
      if(VarMax_>1.0E-2) dT0=dT0/2.0D0 !log(1.080D0)
      !
    end select
    !
  end if
  !
  if(dT0>dTMax) dT0=dTMax
  !
  dTime_=dT0
  if (Time +dTime_ > tFinal) then !!!??? dTime ???
    dTime_=tFinal -Time ! go directly to the solution now
  else if(Time +2*dTime_ > tFinal) then
    dTime_= 0.5*(tFinal -Time) ! slow down two steps before the wall
  end if
  !
end subroutine Dynam_CtrlDtime

subroutine ShoMin(f,iStep_,VBox,TimeScaled,pH_,PhiF,vMolK,vKinMVol,vKinQsk)
  use M_Numeric_Const, only: Ln10, TinyDP
  !
  integer, intent(in):: f,iStep_
  real(dp),intent(in):: VBox,TimeScaled,pH_,PhiF
  real(dp),dimension(:),intent(in):: vMolK,vKinMVol,vKinQsk
  !
  integer::I,N
  !
  N=size(vMolK)
  write(f,'(I7,A1,3(G15.8,A1),$)') iStep_,T_,TimeScaled,T_,pH_,T_,PhiF,T_
  do i=1,N; write(f,'(G15.8,A1,$)') vMolK(i)*vKinMVol(i)/VBox,T_ ; end do
  do i=1,N; write(f,'(G15.8,A1,$)') log(vKinQsK(i))/Ln10,T_      ; end do
  do i=1,N; write(f,'(G15.8,A1,$)') vMolK(i),T_                  ; end do
  write(f,*)
  !
end subroutine ShoMin

subroutine ShoMinRate( & !
!--
!-- write details of mineral rates (surface, vVmQsK, vVmAct, ...)
!--
& f,iStep_,VBox,TimeScaled,pH_,PhiF, & !
& vMolK,vKinMVol,vKinQsk,vSurfK,vSurfK0,vVmQsK,vVmAct, & !
& vKinFas)
  use M_Numeric_Const,only: Ln10 !, TinyDP
  use M_T_KinFas,only: T_KinFas
  !------
  integer,                    intent(in):: f,iStep_
  real(dp),                   intent(in):: VBox,TimeScaled,pH_,PhiF
  real(dp),      dimension(:),intent(in):: vMolK,vKinMVol,vKinQsk,vSurfK,vSurfK0
  real(dp),      dimension(:),intent(in):: vVmQsK,vVmAct
  type(T_KinFas),dimension(:),intent(in):: vKinFas
  !
  integer::I,N
  !------
  
  N=size(vMolK)
  write(f,'(I7,A1,3(G15.8,A1),$)') iStep_,T_,TimeScaled,T_,pH_,T_,PhiF,T_
  !
  do i=1,n; write(f,'(G15.8,A1,$)') vMolK(i)*vKinMVol(i)/VBox,  T_ ; end do
  do i=1,n; write(f,'(G15.8,A1,$)') vSurfK(i)/vSurfK0(i),       T_ ; end do
  do i=1,n; write(f,'(G15.8,A1,$)') vKinFas(i)%Dat%SurfKg/1.D3, T_ ; end do
  do i=1,n; write(f,'(G15.8,A1,$)') log(vKinQsK(i))/Ln10,       T_ ; end do
  do i=1,n; write(f,'(G15.8,A1,$)') vVmQsK(i),                  T_; end do
  do i=1,n; write(f,'(G15.8,A1,$)') vVmAct(i),                  T_; end do
  do i=1,n; write(f,'(G15.8,A1,$)') log(vKinFas(i)%Dat%Radius)/Ln10,T_; end do
  !
  write(f,*)
  
  return
end subroutine ShoMinRate

end module M_Dynam_OneStep_Tools

