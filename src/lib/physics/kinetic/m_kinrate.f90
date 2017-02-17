module M_KinRate
!--
!-- rate calculations on kinetic species
!--

  use M_Kinds
  use M_Trace,only: T_,fHtm,Pause_
  
  implicit none
  
  private
  
  public:: KinRate_CalcActivFactor
  public:: KinRate_CalcQsK
  public:: KinRate_SatState
  public:: KinRate_CalcQsKFactor
  !
  public:: KinRate_ActivTest
  
  real(dp),parameter:: QsK_Max=  1.D+12 !
  logical, parameter:: LimitQsK= .true.
  
contains

subroutine KinRate_CalcActivFactor(&
!-----------------------------------------------------------------------
!-- compute VmAct, the factor depending on fluid chemistry only
!-----------------------------------------------------------------------
& cSatur, &  !IN:  "D" / "P" ! "DISSOLU"/"PRECIPI"
& M,      &  !IN:  kinetic model of mineral
& vLnAct, &  !IN:  fluid species activities
& VmAct)     !OUT: factor depending on fluid chemistry only
!! & dVmAdLnX_M) !OUT: its derivative vs ln(X)
  !
  use M_T_KinFas,  only: T_KinFas
  use M_T_Kinmodel,only: T_KinModel
  !
  character,       intent(in) :: cSatur
  type(T_KinModel),intent(in) :: M
  real(dp),        intent(in) :: vLnAct(:)
  !
  real(dp), intent(out):: VmAct
  ! real(dp), intent(out):: dVmAdLnX_M(:) !(1:nAq)
  !
  integer ::I
  real(dp)::X
  !
  ! dVmAdLnX_M=Zero
  !provisionally, dVmAdLnX_M is not computed in this new version  !!!!!!!!!!!!!
  !-> we assume that this factor will not be implicited           !!!!!!!!!!!!!
  !
  !! write(51,*) -LnActH_/Ln10, char(9), -LnActOH/Ln10
  VmAct=Zero
  select case(cSatur)

  case("D") !"DISSOLU")
    do I=1,M%NTermD
      X= exp( M%N_d(I) *vLnAct(M%ISpcDiss(I)) )
      if(M%JSpcDiss(I)>0) X= X * exp( M%NJ_d(I)*vLnAct(M%JSpcDiss(I)) )
      VmAct= VmAct + M%Kd(I) * X
    end do
    !if(M%Special/=0) then !-> additional terms for carbonates, sulfides ...
    !  select case(M%Special)
    !    case(1) !CALCITE
    !    case(2) !etcetera
    !    .......
    !  end select
    !end if

  case("P") !"PRECIPI")
    do I=1,M%NTermP
      X= exp( M%N_p(I) *vLnAct(M%ISpcPrec(I)) )
      if(M%JSpcPrec(I)>0) X= X * exp( M%NJ_p(I) *vLnAct(M%JSpcPrec(I)) )
      VmAct= VmAct + M%Kp(I) * X
    end do

  end select
  !
  return
end subroutine KinRate_CalcActivFactor

!-----------------------------------------------------------------------
!-- calc. Q/K and d(Q/K)/dLnX of a phase of properties DG,vNu
!-- for a fluid composition vLnX
!-- (composition input in log(mole numbers))
!-----------------------------------------------------------------------
subroutine KinRate_CalcQsK( &
& nCi,       & !IN
& nCx,       & !IN
& DG0rt,     & !IN
& vNu,       & !IN
& vLnX,      & !IN
& vLnGam,    & !IN
& vLnActBuf, & !IN
& LnActW,    & !IN
& QsK,       & !OUT
& dQsKdLnXi)   !OUT

  use M_Numeric_Const,     only: MinExpDP,MaxExpDP
  use M_Basis_Vars,only: isW,MWSv
  !
  integer, intent(in):: nCi          !nr inert comp'nt
  integer, intent(in):: nCx          !nr mobile comp'nt
  real(dp),intent(in):: DG0rt        !deltaG/RT of formation reaction
  real(dp),intent(in):: vNu(:)       !stoikio of formation reaction
  real(dp),intent(in):: vLnX(:)      !log(Mole Nrs. Prim.Species)  !1:nCi
  real(dp),intent(in):: vLnGam(:)    !log(Act.Coeff. Prim.Species) !1:nCi
  real(dp),intent(in):: vLnActBuf(:) !log(ActivityBufferSpecies)   !1:nCx
  real(dp),intent(in):: LnActW       !log(Activity Solvent)
  !
  real(dp),intent(out):: QsK
  real(dp),intent(out):: dQsKdLnXi(:) !1:nCi
  !
  real(dp)::X
  !
  !lnAct= lnMolal + lnGamma=  lnMol + lnGam - lnMolH2O - lnMWH2O
  !QsK=exp(dot_product(tNuMk(iMk,1:nCp),vLnAct(1:nCp)) - vDG_Mk(iMk))
  !or,directly from G0,
  !QsK=exp(dot_product(tNuMk(iMk,1:nCp),vLnAct(1:nCp) +vSpc(1:nCp)%G0) - vKinFas(iMk)%G0)
  !
  !caveat:
  !Solvent should follow mole fraction scale, not molality ...
  !-> here, we assume explicited Solvent activity:
  X=            vNu(isW)  *LnActW &
  + dot_product(vNu(2:nCi),vLnX(2:nCi)) & !log(mole numbers)
  + dot_product(vNu(2:nCi),vLnGam(2:nCi)) & !log(gammas)
  -         SUM(vNu(2:nCi))*(vLnX(isW)+log(MWSv)) & !log(mass Solvent)
  - DG0rt !equiv -logK
  !
  if(nCx>0) & != mobile species
  & X= X &
  &  + dot_product(vNu(nCi+1:nCi+nCx),vLnActBuf(1:nCx))
  !
  if(X>MinExpDP .and. X<MaxExpDP) then
    QsK=exp(X)
  else
    if(X<=MinExpDP) QsK=exp(MinExpDP)
    if(X>=MaxExpDP) QsK=exp(MaxExpDP)
  end if
  !
  !<new 200911>
  if(LimitQsK) QsK= MIN(QsK,QsK_Max)
  !</new>
  !
  dQsKdLnXi(:)=     Zero
  dQsKdLnXi(isW)= - SUM(vNu(2:nCi))*QsK
  dQsKdLnXi(2:nCi)=     vNu(2:nCi) *QsK
  !
  return
end subroutine KinRate_CalcQsK

!-----------------------------------------------------------------------
!= from QsK deduce cSatur, for branching to dissol/precip
!-----------------------------------------------------------------------
subroutine KinRate_SatState( &
& M,         & !IN
& nMol,      & !IN
& nMolMinim, & !IN
& QsK,       & !IN
& Iota,      & !IN !mod 10/06/2008 17:02 added
& cSatur)      !OUT
  use M_T_KinFas,only: T_KinFas
  !
  type(T_KinFas), intent(in) :: M      !IN
  real(dp),       intent(in) :: nMol   !IN, Nr Moles of mineral M
  real(dp),       intent(in) :: nMolMinim !IN
  real(dp),       intent(in) :: QsK    !IN
  real(dp),       intent(in) :: Iota   !IN
  character,      intent(out):: cSatur !OUT
  !
  cSatur=M%Dat%cSat != current value, cSatur will be the new value
  !
  cSatur="I" ! INERT, normal assumption by default
  if(QsK < One - Iota) then
    if(nMol<=nMolMinim) then  ;  cSatur="M" ! "MINIMAL"
    else                      ;   cSatur="D"  ! "DISSOLU"
    end if
  elseif(QsK >= M%QsKSeuil + Iota) then
    cSatur="P" ! PRECIPI
  end if
  !
end subroutine KinRate_SatState
  !
  ! if ( QsK < One - Iota ) then
  ! if ( nMol > nMolMinim ) then
  ! cSatur="DISSOLU"
  ! else
  ! cSatur="MINIMAL"
  ! end if
  ! elseif ( QsK > QsKseuil + Iota ) then
  ! cSatur="PRECIPI"
  ! else
  ! cSatur="INERT"
  ! end if


  !one possible sequence:
  !!!  if (QsK<One .and. nMol>nMolMinim .and. cSatur/="MINIMAL") then
  !!!    cSatur="DISSOLU"
  !!!  elseif (QsK<=One.and.nMol<=nMolMinim) then
  !!!    cSatur="MINIMAL" !QsK<1, but amount Too low to dissolve anymore
  !!!  elseif (QsK>M%Dat%QsKSeuil .or. (QsK>One .and. M%Dat%cSat=="PRECIPI")) then
  !!!    cSatur="PRECIPI" !supersaturation or min. already precipitating
  !!!  end if
  !another possible sequence:
  !
  !!  !! if(ABS(Qsk - One) < Iota) then
  !!  !!   cSatur="INERT"
  !!  !! else!
  !!  !!   select case(cSatur)
  !!  !!     case("PRECIPI"); if(QsK<One)   cSatur="DISSOLU" !else it continues as "PRECIPI"
  !!  !!     case("DISSOLU"); if(QsK>One)   cSatur="PRECIPI" !else it continues as "DISSOLU"
  !!  !!     case("MINIMAL")
  !!  !!       if(QsK > One + Iota)        cSatur="PRECIPI"
  !!  !!     case("INERT")
  !!  !!       if(QsK > M%Dat%QsKSeuil + Iota) cSatur="PRECIPI" !else it continues as "INERT_"
  !!  !!       if(QsK < One - Iota)        cSatur="DISSOLU" !
  !!  !!     !!case("PRIMARY")
  !!  !!     !!  if(QsK > One + Iota)        cSatur="PRECIPI" !
  !!  !!     !!  if(QsK < One - Iota)        cSatur="DISSOLU" !
  !!  !!     !!case("SECONDA")
  !!  !!     !!  if(QsK > M%Dat%QsKSeuil + Iota) cSatur="PRECIPI"
  !!  !!   end select
  !!  !! end if
  !!  !! if(cSatur=="DISSOLU" .and. nMol<nMolMinim) cSatur="MINIMAL"
  !
  !select case(cSatur)
  !case("PRECIPI"); if(QsK<One)   cSatur="DISSOLU" !else it continues as "PRECIPI"
  !case("DISSOLU"); if(QsK>=One)  cSatur="PRECIPI" !else it continues as "DISSOLU"
  !case("MINIMAL")
  !  if(QsK >= One + Iota)        cSatur="PRECIPI"
  !case("INERT")
    !if(QsK >= M%QsKSeuil + Iota) cSatur="PRECIPI" !else it continues as "INERT_"
    !if(QsK < One - Iota)         cSatur="DISSOLU" !
  !case("PRIMARY")
  !  if(QsK >= One + Iota)        cSatur="PRECIPI" !
  !  if(QsK < One - Iota)        cSatur="DISSOLU" !
  !case("SECONDA")
  !  if(QsK >= M%QsKSeuil + Iota) cSatur="PRECIPI"
  !end select
  !! if(cSatur=="DISSOLU" .and. nMol<=nMolMinim) cSatur="MINIMAL"
  !! if(cSatur=="DISSOLU" .and. nMol<=1.E-6) cSatur="MINIMAL"

subroutine KinRate_CalcQsKFactor(&
& cSatur,        & !IN
& M,             & !IN: Kinetic model
& QsK,           & !IN
& dQsKdLnXi,     & !IN
& VmQsK,         & !OUT
& dVmQdLnX_M)      !OUT
  !
  use M_T_Kinmodel,only: T_KinModel
  !
  character,       intent(in):: cSatur
  type(T_KinModel),intent(in):: M
  real(dp),        intent(in):: QsK
  real(dp),        intent(in):: dQsKdLnXi(:) !(1:nCi)
  !
  real(dp), intent(out):: VmQsK
  real(dp), intent(out):: dVmQdLnX_M(:) !(1:nAq)
  !
  real(dp):: dVmQdQsK
  integer :: nCi_
  !
  nCi_= size(dQsKdLnXi)
  VmQsK=     Zero !-> values for cSatur=="MINIMAL"
  dVmQdLnX_M=Zero !
  dVmQdQsK=  Zero
  !
  select case(cSatur)

  case("D") !("DISSOLU")
  ! QsK<1 -> (One - QsK**M%AlfaD)>0 -> VmQsK<0
  ! QsK<1 -> -log(QsK)>0 -> VmQsK<0
    
    if(M%BetaD /=Zero) then
      VmQsK=    - (One - QsK**M%AlfaD)**M%BetaD
      dVmQdQsK= M%BetaD &
      &         *(One - QsK**M%AlfaD)**(M%BetaD-One) &
      &         *M%AlfaD *QsK**(M%AlfaD-One)
    else
    !cf Soler_Lasaga
      VmQsK=    (-log(QsK))**M%AlfaD
      dVmQdQsK= M%AlfaD *VmQsK /QsK /log(QsK)
    end if

  case("P") !("PRECIPI")
  ! QsK>1 -> (QsK**M%AlfaD - One)>0 -> VmQsK>0
  ! QsK>1 -> log(QsK)>0 -> VmQsK>0
    
    if(M%BetaP /=Zero) then
      VmQsK=(QsK**M%AlfaP - One)**M%BetaP
      !if(VmF<0.0) then
      !  cSatur="INERT" !"I"nert mineral
      !  VmF=0.0
      !  !test!!!dVSdQsK=0.0
      !  return
      !end if
      dVmQdQsK= M%BetaP &
      &       *(QsK**M%AlfaP - One)**(M%BetaP-One) &
      &       *M%AlfaP *QsK**(M%AlfaP-One)
    else
    !cf Soler_Lasaga
      VmQsK=    (log(QsK))**M%AlfaD
      dVmQdQsK= M%AlfaD *VmQsK /QsK /log(QsK)
    end if

  case("M") !("MINIMAL")
    VmQsK=    Zero
    dVmQdQsK= Zero
    
  end select

  dVmQdLnX_M(1:nCi_)= dVmQdQsK*dQsKdLnXi(1:nCi_)

  !! write(51,'(A,3A1,3(G15.9,A1))') &
  !! & trim(M%Name),T_,cSatur,T_,log(QsK)/Ln10,T_,VmQsK,T_,QsK**M%AlfaD,T_
  return
end subroutine KinRate_CalcQsKFactor

subroutine KinRate_ActivTest(fMnK,iCount) !,LnActH_,LnActOH,LnActCO2)
!-----------------------------------------------------------------------
!= test the "validity" of KinRate_CalcActivFactor
!= of module M_T_KinFas for rate calculations
!= (= the factor depending on species'activities)
!-----------------------------------------------------------------------
  use M_IoTools,  only: GetUnit
  use M_Files,    only: DirOut,cTitle,Files_Index_Write !,bOpenMnk
  !
  use M_Numeric_Const,     only: Ln10,TinyDP
  use M_Dtb_Const, only: T_CK
  !
  use M_T_Species,  only: Species_Index
  use M_Basis_Vars, only: nMx, isH_
  use M_T_Kinmodel, only: T_KinModel,KinModel_CalcCoeffs,KinModel_PrmIndex
  !
  use M_System_Vars,only: TdgK
  use M_Global_Vars,only: nAq,vSpc,vKinModel
  
  ! use M_KinRate,    only: KinRate_CalcActivFactor
  
  integer,intent(inout):: fMnK
  integer,intent(in)   :: iCount
  !
  type(T_KinModel),allocatable::vKinMod(:)
  type(T_KinModel)::M
  integer :: I,iKm,nKm
  real(dp):: VmAct,X
  ! real(dp):: dum(1:nAq)
  real(dp):: vLnAct(1:nAq) !for dissolution rates
  integer :: vPrm(1:nAq)
  real(dp):: pH
  
  !! print *,'KinRate_ActivTest'  ;  call Pause_
  
  nKm=size(vKinModel)
  !
  if(nKm==0) return
  !
  do iKm=1,nKm
    call KinModel_CalcCoeffs(vKinModel(iKm),One,TdgK)
  end do
  !
  vLnAct(1:nAq)=vSpc(1:nAq)%Dat%LAct
  do I=1,nAq; vPrm(i)=i; end do
  !
  allocate(vKinMod(1:size(vKinModel))); vKinMod=vKinModel
  !
  do iKm=1,nKm
    call KinModel_PrmIndex(vPrm,vKinModel(iKm),vKinMod(iKm))
  end do
  !
  pH=-vLnAct(isH_)/Ln10
  !
  if(fMnk==0) then
  
    call GetUnit(fMnk)
    
    open(fMnk,file=trim(DirOut)//"_vmact.restab")
    
    call Files_Index_Write(fHtm, &
    & trim(DirOut)//"_vmact.restab", &
    & "dependence of mineral dissolution rates on solution chemistry" )
    write(fMnk,'(3(A,A1))',advance= 'NO') "count",T_,"TdgC",T_,"pH",T_ !,"pOH",T_
    
    do iKm=1,nKm
      write(fMnk,'(A12,A1)',advance= 'NO') vKinMod(iKm)%Name,T_
    end do
    write(fMnk,*)
  
  end if
  !
  write(fMnk,'(I3,A1,2(G15.3,A1))',advance= 'NO') iCount,T_,TdgK-T_CK,T_,pH,T_
  !
  do iKm=1,nKm
    
    M=vKinMod(iKm)
    !iMn=Species_Index(M%Name,vSpc)-nAq-nMx
    
    !if(iMn>0) then !mineral is not "mobile"
      call KinRate_CalcActivFactor(&
      & "D",      & !IN: "D"/"P"
      & M,        & !IN: kinetic model of mineral
      & vLnAct,   & !LnActH_,LnActOH,LnActCO2, & !"activators"
      & VmAct)      !OUT
      !! & dum) !dVmAdLnX_M) !OUT
      !X=ABS(VmAct)
      X=VmAct
      if(X>TinyDP) then; X=log(X)/Ln10
      else             ; X=Zero
      end if
      write(fMnk,'(G15.6,A1)',advance="no") X,T_ !/(QsK-One)
    !end if
    
  end do
  !
  write(fMnk,*)
  !
  deallocate(vKinMod)
  
  return
end subroutine KinRate_ActivTest

end module M_KinRate
