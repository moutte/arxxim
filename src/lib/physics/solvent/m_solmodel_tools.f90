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
  
contains

subroutine Solmodel_pHpE(SM,vSpc,vSpcDat,pH_,pE_)
  use M_Numeric_Const,only: LN10
  use M_T_Species,    only: T_Species, T_SpcData
  use M_T_SolModel,   only: T_SolModel
  !
  type(T_SolModel),intent(in) :: SM !SolModel
  type(T_Species), intent(in) :: vSpc(:)
  type(T_SpcData), intent(in) :: vSpcDat(:)
  real(dp),        intent(out):: pH_, pE_
  
  !pH_= - (log(vMolF(isH_)) + vLnGam(isH_) - log(vMolF(isW)*MolWeitSv))/LN10
  pH_= -vSpcDat(SM%isH_)%LAct/LN10
  pE_=  Zero
  
  if(SM%isO2/=0) &
  & pE_= & 
  & (  vSpc(SM%isO2)%G0rt + vSpcDat(SM%isO2)%LAct &
  &  -(vSpc(SM%iSolvent)%G0rt + vSpcDat(SM%iSolvent)%LAct) *Two  ) &
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
& SolModel, &
& vSpc,vSpcDat, &
& vMolal,IonStrTrue)
!-- from S%Mole, S%Z,
!-- calc. Molalities vMolal, and "true" ionic strength
!--
!-- we assume that SolModel has index 1 in vSpc 
!-- and solute species are indexed 2:nAq !!!
!--
  use M_T_SolModel,only: T_SolModel
  use M_T_Species, only: T_Species,T_SpcData
  !---------------------------------------------------------------------
  type(T_SolModel),intent(in) :: SolModel
  type(T_Species), intent(in) :: vSpc(:)
  type(T_SpcData), intent(in) :: vSpcDat(:)
  real(dp),        intent(out):: vMolal(:)
  real(dp),        intent(out):: IonStrTrue !,IonStrStoik
  !---------------------------------------------------------------------
  real(dp):: MolWeitSv
  integer :: isW
  integer :: I,J
  !
  isW= SolModel%iSolvent
  MolWeitSv=  SolModel%MolWeitSv
  !
  IonStrTrue= Zero
  J= 0
  do I=1,size(vSpc)
    if(vSpc(I)%Typ=="AQU") then
      J= J+1
      vMolal(J)= vSpcDat(I)%Mole /vSpcDat(isW)%Mole /MolWeitSv
      IonStrTrue= IonStrTrue + vSpc(I)%Z*vSpc(I)%Z*vMolal(J)
    end if
  end do
  IonStrTrue= IonStrTrue /Two
  !
end subroutine Solmodel_CalcMolal

subroutine Solmodel_CalcMolal2( &
& SolModel, &
& vSpc,vSpcDat,nAq, &
& vMolal,IonStrTrue) !,IonStrStoik
!-- from S%Mole, S%Z, 
!-- calc. Molalities vMolal, and ionic strengths ("true")
!-- we assume that SolModel has index 1 in vSpc 
!-- and solute species are indexed 2:nAq !!!
  use M_T_SolModel,only: T_SolModel
  use M_T_Species, only: T_Species,Species_Index,T_SpcData
  !
  type(T_SolModel),intent(in) :: SolModel
  type(T_Species), intent(in) :: vSpc(:)
  type(T_SpcData), intent(in) :: vSpcDat(:)
  integer,         intent(in) :: nAq
  real(dp),        intent(out):: vMolal(:)
  real(dp),        intent(out):: IonStrTrue !,IonStrStoik
  !
  integer, dimension(1:nAq):: vZ2
  real(dp):: MolWeitSv
  integer :: isW
  !
  isW= SolModel%iSolvent
  MolWeitSv=  SolModel%MolWeitSv
  !
  vMolal(1:nAq)= vSpcDat(1:nAq)%Mole /vSpcDat(isW)%Mole /MolWeitSv
  !
  vZ2(1:nAq)=    vSpc(1:nAq)%Z *vSpc(1:nAq)%Z
  !
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
  
  ! if(iDebug>2) then
  !   if(Ok_Rho)  print '(A,G15.6)',"Rho= ",SolModel%Dat%Rho
  !   if(Ok_Eps)  print '(A,G15.6)',"Eps= ",SolModel%Dat%Eps
  !   if(Ok_DHA)  print '(A,G15.6)',"DHA= ",SolModel%Dat%DHA
  !   if(Ok_DHB)  print '(A,G15.6)',"DHB= ",SolModel%Dat%DHB
  !   if(Ok_BDot) print '(A,G15.6)',"BDot=",SolModel%Dat%BDot
  ! end if
  
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

