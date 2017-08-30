module M_SolModel_Tools
!--
!--
  use M_Kinds
  use M_Trace,only: Stop_,fTrc,T_,iDebug,pause_
  use M_T_DtbLogKTbl,only: DimLogK_Max
  !
  implicit none
  !
  private
  !
  public:: Solmodel_pHpE
  public:: Solmodel_CalcMolal
  public:: Solmodel_CalcMolal2
  public:: SolModel_TP_Update
  public:: Solmodel_Init_TotO_TotH
  
contains

!subroutine Solmodel_Zero
!  !default SolModel parameters, in case not found in input
!  isW= 1; iH_= 0; iOx= 0
!  isH_=0; isOH=0; isO2=0
!end subroutine Solmodel_Zero

subroutine Solmodel_Init_TotO_TotH( & !
& AdjustMode, &                       !
& iBal,ic_W,ic_H,vEle, &              !IN
& vCpn)                               !INOUT
!--
!-- adjust total molar amounts of O and H according 
!-- to molar amounts of other elements
!-- -> values of vTotF
!--
  use M_T_Element,  only: T_Element
  use M_T_Component,only: T_Component
  !
  integer,          intent(in)   :: AdjustMode
  integer,          intent(in)   :: iBal
  integer,          intent(in)   :: ic_W,ic_H
  type(T_Element),  intent(in)   :: vEle(:)
  type(T_Component),intent(inout):: vCpn(:)
  !
  !
  real(dp):: Zbalance,TotO_
  integer :: I,iBal_

  if(iDebug>0) write(fTrc,'(/,A)') "< Solmodel_Init_TotO_TotH"

  TotO_= vCpn(ic_W)%Mole

  if(iBal==0) then ;  iBal_= ic_H
  else             ;  iBal_= iBal
  end if

  if(iBal_ /= ic_H .and. vCpn(ic_H)%Statut=="INERT") then
    vCpn(ic_W)%Mole= TotO_
    vCpn(ic_H)%Mole= TotO_*Two
    !return
  end if
  
  Zbalance= Zero
  do I=1,size(vCpn)
    if(vCpn(I)%Statut=="INERT" .and. I/=ic_W .and. I/=ic_H) then
      Zbalance= Zbalance + vEle(vCpn(I)%iEle)%Z *vCpn(I)%Mole
    end if
  enddo

  !! select case(AdjustMode)
  !!   !
  !!   case(1)
  !!     !--- if excess cation, Na>Cl, adjust with OH
  !!     if(Zbalance > Zero) then
  !!       vCpn(ic_W)%Mole= TotO_     +Zbalance
  !!       vCpn(ic_H)%Mole= TotO_*Two +Zbalance
  !!     !--- if excess anion, Cl>Na, adjust with H
  !!     else
  !!       vCpn(ic_W)%Mole= TotO_
  !!       vCpn(ic_H)%Mole= TotO_*Two -Zbalance
  !!     end if
  !!   !
  !!   case(2)
  !!     !--- Oxygen number unchanged, equilibrate using hydrogen only
  !!     vCpn(ic_W)%Mole= TotO_
  !!     vCpn(ic_H)%Mole= TotO_*Two -Zbalance
  !!   !
  !! end select

  select case(AdjustMode)
  
  case(1)
    !--- if excess cation, Na>Cl, adjust with OH
    if(Zbalance > Zero) then
      vCpn(ic_W)%Mole= TotO_     +Zbalance
      vCpn(iBal_)%Mole= TotO_*Two +Zbalance
    !--- if excess anion, Cl>Na, adjust with H
    else
      vCpn(ic_W)%Mole=  TotO_
      vCpn(iBal_)%Mole= (TotO_*Two -Zbalance) /real(vEle(iBal_)%Z)
    end if

  case(2)
    !--- Oxygen number unchanged, equilibrate using hydrogen only
    vCpn(ic_W)%Mole= TotO_
    vCpn(iBal_)%Mole= TotO_*Two -Zbalance

  end select
  
  if(iDebug>2) then
    write(fTrc,'(A,G15.6)')  "Zbalance       =", Zbalance
    write(fTrc,'(A,G15.6)')  "vCpn(ic_W)%Mole=", vCpn(ic_W)%Mole
    write(fTrc,'(A,G15.6)')  "vCpn(ic_H)%Mole=", vCpn(ic_H)%Mole
  end if

  if(iDebug>0) write(fTrc,'(A,/)') "</ Solmodel_Init_TotO_TotH"
  
end subroutine Solmodel_Init_TotO_TotH
  !
subroutine Solmodel_pHpE(isW,iOx,isH_,isO2,vSpc,pH_,pE_)
  use M_Numeric_Const,only: LN10
  use M_T_Species,    only: T_Species
  !
  integer,        intent(in) :: isW,iOx,isH_,isO2
  type(T_Species),intent(in) :: vSpc(:)
  real(dp),       intent(out):: pH_, pE_
  
  !pH_= - (log(vMolF(isH_)) + vLnGam(isH_) - log(vMolF(isW)*MolWeitSv))/LN10
  pH_= -vSpc(isH_)%Dat%LAct/LN10
  pE_=  Zero
  
  if(iOx/=0 .and. isO2/=0) &
    pE_= & 
    & (  vSpc(isO2)%G0rt + vSpc(isO2)%Dat%LAct          &
    &  -(vSpc(isW )%G0rt + vSpc(isW )%Dat%LAct) *Two  ) &
    & /LN10 /4.0D0 &
    & - pH_
  
  !2.H2O=  O2 + 4.H+ + 4.e-
  !2.G(H2O)+lnA_H2O=  G(O2)+lnA_O2  +4.lnA_H+ +4.lnA_e-
  !4.(pH + pE)=  [G(O2) +lnA_O2 -2.G(H2O) -2.lnA_H2O]/ln10
  !log10(Act(O2))=4*(pH_+pE_)+ 2[G(H2O) + lnA_H2O] - [G(O2) + lnA_O2]
  !
  !if Redox species is H2aq:
  !2H+ + 2E-=  H2 -> pH + pE=  -( vSpc(iOx)%G0 + vLnAct(iOx)) /2
  
end subroutine Solmodel_pHpE

subroutine Solmodel_CalcMolal( &
& vSpc,isW, &
& vMolal,IonStrTrue)
!-- from S%Mole, S%Z,
!-- calc. Molalities vMolal, and ionic strengths ("true")
!--
!-- we assume that SolModel has index 1 in vSpc 
!-- and solute species are indexed 2:nAq !!!
!--
  use M_T_Species,only: T_Species
  !
  type(T_Species),intent(in) :: vSpc(:)
  integer,        intent(in) :: isW
  real(dp),       intent(out):: vMolal(:)
  real(dp),       intent(out):: IonStrTrue !,IonStrStoik
  !
  real(dp):: MolWeitSv
  integer :: I,J
  !
  IonStrTrue= Zero
  MolWeitSv=  vSpc(isW)%WeitKg
  J= 0
  do I=1,size(vSpc)
    if(vSpc(I)%Typ=="AQU") then
      J= J+1
      vMolal(J)= vSpc(I)%Dat%Mole /vSpc(isW)%Dat%Mole /MolWeitSv
      IonStrTrue= IonStrTrue + vSpc(I)%Z*vSpc(I)%Z*vMolal(I)
    end if
  enddo
  
  IonStrTrue= IonStrTrue /Two
  
end subroutine Solmodel_CalcMolal

subroutine Solmodel_CalcMolal2( &
& vSpc,nAq, &
& vMolal,IonStrTrue) !,IonStrStoik
!-- from S%Mole, S%Z, 
!-- calc. Molalities vMolal, and ionic strengths ("true")
!-- we assume that SolModel has index 1 in vSpc 
!-- and solute species are indexed 2:nAq !!!
  use M_T_Species,only: T_Species,Species_Index
  !
  type(T_Species),intent(in) :: vSpc(:)
  integer,        intent(in) :: nAq
  real(dp),       intent(out):: vMolal(:)
  real(dp),       intent(out):: IonStrTrue !,IonStrStoik
  !
  integer, dimension(1:nAq):: vZ2
  real(dp):: MolWeitSv
  integer :: isW
  !
  isW=  Species_Index("H2O",vSpc)
  MolWeitSv=vSpc(isW)%WeitKg
  
  vMolal(1:nAq)= vSpc(1:nAq)%Dat%Mole /vSpc(isW)%Dat%Mole /MolWeitSv
  
  vZ2(1:nAq)=    vSpc(1:nAq)%Z *vSpc(1:nAq)%Z
  
  IonStrTrue=    dot_product(vZ2(2:nAq),vMolal(2:nAq))/Two
  
end subroutine Solmodel_CalcMolal2

subroutine Solmodel_TP_Update_LogK( &
& TdgK,Pbar, &
& Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl, &
& Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot, &
& SolModel)

  use M_Dtb_Const,only: T_CK
  use M_Solmodel_Vars,only: T_Spline
  use M_T_SolModel,only: T_SolModel
  use M_CMM_Spline
  
  real(dp),        intent(in)   :: TdgK,Pbar
  type(T_Spline),  intent(in)   :: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl
  logical,         intent(in)   :: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  type(T_SolModel),intent(inout):: SolModel
  
  integer:: N
  
  !--
  
  if(Ok_Rho) then
    N= Rho_Spl%Dimm
    call CMM_Spline_Eval(n, &
    & Rho_Spl%vX(1:n), Rho_Spl%vY(1:n), &
    & Rho_Spl%vSplineB(1:n),Rho_Spl%vSplineC(1:n),Rho_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%Rho)
  end if
  if(Ok_Eps) then
    N= Eps_Spl%Dimm
    call CMM_Spline_Eval(n, &
    & Eps_Spl%vX(1:n), Eps_Spl%vY(1:n), &
    & Eps_Spl%vSplineB(1:n),Eps_Spl%vSplineC(1:n),Eps_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%Eps)
  end if
  if(Ok_DHA) then
    N= DHA_Spl%Dimm
    call CMM_Spline_Eval(n, &
    & DHA_Spl%vX(1:n), DHA_Spl%vY(1:n), &
    & DHA_Spl%vSplineB(1:n),DHA_Spl%vSplineC(1:n),DHA_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%DHA)
  end if
  if(Ok_Rho) then
    N= DHB_Spl%Dimm
    call CMM_Spline_Eval(n, &
    & DHB_Spl%vX(1:n), DHB_Spl%vY(1:n), &
    & DHB_Spl%vSplineB(1:n),DHB_Spl%vSplineC(1:n),DHB_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%DHB)
  end if
  if(Ok_BDot) then
    N= BDot_Spl%Dimm
    call CMM_Spline_Eval(n, &
    & BDot_Spl%vX(1:n), BDot_Spl%vY(1:n), &
    & BDot_Spl%vSplineB(1:n),BDot_Spl%vSplineC(1:n),BDot_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%BDot)
  end if
  
  if(iDebug>2) then
    if(Ok_Rho)  print '(A,G15.6)',"Rho= ",SolModel%Dat%Rho
    if(Ok_Eps)  print '(A,G15.6)',"Eps= ",SolModel%Dat%Eps
    if(Ok_DHA)  print '(A,G15.6)',"DHA= ",SolModel%Dat%DHA
    if(Ok_DHB)  print '(A,G15.6)',"DHB= ",SolModel%Dat%DHB
    if(Ok_BDot) print '(A,G15.6)',"BDot=",SolModel%Dat%BDot
  end if
  
  return
end subroutine Solmodel_TP_Update_LogK

subroutine SolModel_TP_Update( &
& TdgK,Pbar, &
& SolModel)
!--
!-- For given (T,P) condition,
!-- compute SolModel properties using DtbH2OHkf_Calc and write in SolModel.
!-- If the database is of logK format,
!-- it may contains tabulated values for either Rho, Eps, etc.
!-- apply special procedure that will overwrite SolModel
!--
  use M_Dtb_Const,  only: T_CK
  use M_T_DtbH2OHkf,only: DtbH2OHkf_Calc,T_DtbH2OHkf
  use M_T_SolModel, only: T_SolModel
  use M_Solmodel_Pitzer_Dtb, only: Solmodel_Pitzer_Dtb_TPUpdate
  
  use M_Dtb_Vars,    only: DtbFormat
  use M_Solmodel_Vars,only: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  use M_Solmodel_Vars,only: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl,T_Spline
    
  real(dp),        intent(in)   :: TdgK,Pbar
  type(T_SolModel),intent(inout):: SolModel
  
  type(T_DtbH2OHkf)::PropsH2O
  
  !-- compute solvent (-H2O) properties at T,P
  call DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O)
  
  SolModel%Dat%Rho=  PropsH2O%Rho
  SolModel%Dat%Eps=  PropsH2O%Eps
  SolModel%Dat%dhA=  PropsH2O%dhA
  SolModel%Dat%dhB=  PropsH2O%dhB
  
  SolModel%Dat%bDot= Solmodel_BDot(TdgK)
  
  ! LOGK databases may contain specific values 
  ! for RHO, EPS, DHA, DHB, BdoT
  
  if(DtbFormat=="LOGKTBL") &
  & call Solmodel_TP_Update_LogK( &
  & TdgK,Pbar, &
  & Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl, &
  & Ok_Rho, Ok_Eps, Ok_DHA, Ok_DHB, Ok_BDot, &
  & SolModel)
  
  if(iDebug>2) then
    write(fTrc,'(/,A,/)') "< SolModel_TP_Update"
    write(fTrc,'(5(A,G15.6,/))') &
    & "Rho= ",SolModel%Dat%Rho, & 
    & "Eps= ",SolModel%Dat%Eps, &
    & "dhA= ",SolModel%Dat%dhA, &
    & "dhB= ",SolModel%Dat%dhB, &
    & "Bdot=",SolModel%Dat%bDot
    write(fTrc,'(A,/)') "</ SolModel_TP_Update"
  end if
  
  !call pause_ !!debug!!
  call flush(fTrc)
  
  if(SolModel%iActModel== 8  .or. &    !! "PITZER"
  &  SolModel%iActModel== 11 ) then    !! "SIT"
    call Solmodel_Pitzer_Dtb_TPUpdate(TdgK,Pbar)
  end if
    
  return
end subroutine SolModel_TP_Update

real(dp) function Solmodel_Bdot(TdgK)
  real(dp),intent(in):: TdgK
  
  Solmodel_Bdot=  &
  & -0.616616D0  &
  & +6.90440D-3  *TdgK &
  & -2.73898D-5  *TdgK**2 &
  & +4.87313D-8  *TdgK**3 &
  & -3.26030D-11 *TdgK**4

end function Solmodel_Bdot

end module M_SolModel_Tools

