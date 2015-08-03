MODULE M_SolModel_Tools
!--
!--
  USE M_Kinds
  USE M_Trace,ONLY: Stop_,fTrc,T_,iDebug
  USE M_T_DtbLogKTbl,ONLY: DimLogK_Max
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Solmodel_pHpE
  PUBLIC:: Solmodel_CalcMolal
  PUBLIC:: Solmodel_CalcMolal2
  PUBLIC:: SolModel_TP_Update
  PUBLIC:: Solmodel_Init_TotO_TotH
  
CONTAINS

!SUBROUTINE Solmodel_Zero
!  !default SolModel parameters, in case not found in input
!  isW= 1; iH_= 0; iOx= 0
!  isH_=0; isOH=0; isO2=0
!ENDSUBROUTINE Solmodel_Zero

SUBROUTINE Solmodel_Init_TotO_TotH( & !
& AdjustMode, &                       !
& iBal,ic_W,ic_H,vEle, &              !IN
& vCpn)                               !INOUT
!--
!-- adjust total molar amounts of O and H according 
!-- to molar amounts of other elements
!-- -> values of vTotF
!--
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Component,ONLY: T_Component
  !
  INTEGER,          INTENT(IN)   :: AdjustMode
  INTEGER,          INTENT(IN)   :: iBal
  INTEGER,          INTENT(IN)   :: ic_W,ic_H
  TYPE(T_Element),  INTENT(IN)   :: vEle(:)
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  !
  REAL(dp):: Zbalance,TotO_
  INTEGER :: I,iBal_

  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Solmodel_Init_TotO_TotH"

  TotO_= vCpn(ic_W)%Mole

  IF(iBal==0) THEN ;  iBal_= ic_H
  ELSE             ;  iBal_= iBal
  ENDIF

  IF(iBal_ /= ic_H .AND. vCpn(ic_H)%Statut=="INERT") THEN
    vCpn(ic_W)%Mole= TotO_
    vCpn(ic_H)%Mole= TotO_*Two
    !RETURN
  ENDIF
  
  Zbalance= Zero
  DO I=1,SIZE(vCpn)
    IF(vCpn(I)%Statut=="INERT" .AND. I/=ic_W .AND. I/=ic_H) THEN
      Zbalance= Zbalance + vEle(vCpn(I)%iEle)%Z *vCpn(I)%Mole
    ENDIF
  ENDDO

  !! SELECT CASE(AdjustMode)
  !!   !
  !!   CASE(1)
  !!     !--- IF excess cation, Na>Cl, adjust with OH
  !!     IF(Zbalance > Zero) THEN
  !!       vCpn(ic_W)%Mole= TotO_     +Zbalance
  !!       vCpn(ic_H)%Mole= TotO_*Two +Zbalance
  !!     !--- IF excess anion, Cl>Na, adjust with H
  !!     ELSE
  !!       vCpn(ic_W)%Mole= TotO_
  !!       vCpn(ic_H)%Mole= TotO_*Two -Zbalance
  !!     ENDIF
  !!   !
  !!   CASE(2)
  !!     !--- Oxygen number unchanged, equilibrate using hydrogen only
  !!     vCpn(ic_W)%Mole= TotO_
  !!     vCpn(ic_H)%Mole= TotO_*Two -Zbalance
  !!   !
  !! ENDSELECT

  SELECT CASE(AdjustMode)
  
  CASE(1)
    !--- IF excess cation, Na>Cl, adjust with OH
    IF(Zbalance > Zero) THEN
      vCpn(ic_W)%Mole= TotO_     +Zbalance
      vCpn(iBal_)%Mole= TotO_*Two +Zbalance
    !--- IF excess anion, Cl>Na, adjust with H
    ELSE
      vCpn(ic_W)%Mole=  TotO_
      vCpn(iBal_)%Mole= (TotO_*Two -Zbalance) /REAL(vEle(iBal_)%Z)
    ENDIF

  CASE(2)
    !--- Oxygen number unchanged, equilibrate using hydrogen only
    vCpn(ic_W)%Mole= TotO_
    vCpn(iBal_)%Mole= TotO_*Two -Zbalance

  ENDSELECT
  
  IF(iDebug>2) THEN
    WRITE(fTrc,'(A,G15.6)')  "Zbalance       =", Zbalance
    WRITE(fTrc,'(A,G15.6)')  "vCpn(ic_W)%Mole=", vCpn(ic_W)%Mole
    WRITE(fTrc,'(A,G15.6)')  "vCpn(ic_H)%Mole=", vCpn(ic_H)%Mole
  ENDIF

  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Solmodel_Init_TotO_TotH"
  
ENDSUBROUTINE Solmodel_Init_TotO_TotH
  !
SUBROUTINE Solmodel_pHpE(isW,iOx,isH_,isO2,vSpc,pH_,pE_)
  USE M_Numeric_Const,ONLY: LN10
  USE M_T_Species,    ONLY: T_Species
  !
  INTEGER,        INTENT(IN) :: isW,iOx,isH_,isO2
  TYPE(T_Species),INTENT(IN) :: vSpc(:)
  REAL(dp),       INTENT(OUT):: pH_, pE_
  
  !pH_= - (LOG(vMolF(isH_)) + vLnGam(isH_) - LOG(vMolF(isW)*MolWeitSv))/LN10
  pH_= -vSpc(isH_)%Dat%LAct/LN10
  pE_=  Zero
  
  IF(iOx/=0 .AND. isO2/=0) &
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
  !IF ReDOx species is H2aq:
  !2H+ + 2E-=  H2 -> pH + pE=  -( vSpc(iOx)%G0 + vLnAct(iOx)) /2
  
ENDSUBROUTINE Solmodel_pHpE

SUBROUTINE Solmodel_CalcMolal( &
& vSpc,isW, &
& vMolal,IonStrTrue)
!-- from S%Mole, S%Z,
!-- calc. Molalities vMolal, and ionic strengths ("true")
!--
!-- we assume that SolModel has index 1 in vSpc 
!-- and solute species are indexed 2:nAq !!!
!--
  USE M_T_Species,ONLY: T_Species
  !
  TYPE(T_Species),INTENT(IN) :: vSpc(:)
  INTEGER,        INTENT(IN) :: isW
  REAL(dp),       INTENT(OUT):: vMolal(:)
  REAL(dp),       INTENT(OUT):: IonStrTrue !,IonStrStoik
  !
  REAL(dp):: MolWeitSv
  INTEGER :: I,J
  !
  IonStrTrue= Zero
  MolWeitSv=  vSpc(isW)%WeitKg
  J= 0
  DO I=1,SIZE(vSpc)
    IF(vSpc(I)%Typ=="AQU") THEN
      J= J+1
      vMolal(J)= vSpc(I)%Dat%Mole /vSpc(isW)%Dat%Mole /MolWeitSv
      IonStrTrue= IonStrTrue + vSpc(I)%Z*vSpc(I)%Z*vMolal(I)
    ENDIF
  ENDDO
  
  IonStrTrue= IonStrTrue /Two
  
ENDSUBROUTINE Solmodel_CalcMolal

SUBROUTINE Solmodel_CalcMolal2( &
& vSpc,nAq, &
& vMolal,IonStrTrue) !,IonStrStoik
!-- from S%Mole, S%Z, 
!-- calc. Molalities vMolal, and ionic strengths ("true")
!-- we assume that SolModel has index 1 in vSpc 
!-- and solute species are indexed 2:nAq !!!
  USE M_T_Species,ONLY: T_Species,Species_Index
  !
  TYPE(T_Species),INTENT(IN) :: vSpc(:)
  INTEGER,        INTENT(IN) :: nAq
  REAL(dp),       INTENT(OUT):: vMolal(:)
  REAL(dp),       INTENT(OUT):: IonStrTrue !,IonStrStoik
  !
  INTEGER, DIMENSION(1:nAq):: vZ2
  REAL(dp):: MolWeitSv
  INTEGER :: isW
  !
  isW=  Species_Index("H2O",vSpc)
  MolWeitSv=vSpc(isW)%WeitKg
  
  vMolal(1:nAq)= vSpc(1:nAq)%Dat%Mole /vSpc(isW)%Dat%Mole /MolWeitSv
  
  vZ2(1:nAq)=    vSpc(1:nAq)%Z *vSpc(1:nAq)%Z
  
  IonStrTrue=    DOT_PRODUCT(vZ2(2:nAq),vMolal(2:nAq))/Two
  
ENDSUBROUTINE Solmodel_CalcMolal2

SUBROUTINE Solmodel_TP_Update_LogK( &
& TdgK,Pbar, &
& Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl, &
& Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot, &
& SolModel)

  USE M_Dtb_Const,ONLY: T_CK
  USE M_Solmodel_Vars,ONLY: T_Spline
  USE M_T_SolModel,ONLY: T_SolModel
  USE M_CMM_Spline
  
  REAL(dp),        INTENT(IN)   :: TdgK,Pbar
  TYPE(T_Spline),  INTENT(IN)   :: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl
  LOGICAL,         INTENT(IN)   :: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  TYPE(T_SolModel),INTENT(INOUT):: SolModel
  
  INTEGER:: N
  
  !--
  
  IF(Ok_Rho) THEN
    N= Rho_Spl%Dimm
    CALL CMM_Spline_Eval(n, &
    & Rho_Spl%vX(1:n), Rho_Spl%vY(1:n), &
    & Rho_Spl%vSplineB(1:n),Rho_Spl%vSplineC(1:n),Rho_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%Rho)
  ENDIF
  IF(Ok_Eps) THEN
    N= Eps_Spl%Dimm
    CALL CMM_Spline_Eval(n, &
    & Eps_Spl%vX(1:n), Eps_Spl%vY(1:n), &
    & Eps_Spl%vSplineB(1:n),Eps_Spl%vSplineC(1:n),Eps_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%Eps)
  ENDIF
  IF(Ok_DHA) THEN
    N= DHA_Spl%Dimm
    CALL CMM_Spline_Eval(n, &
    & DHA_Spl%vX(1:n), DHA_Spl%vY(1:n), &
    & DHA_Spl%vSplineB(1:n),DHA_Spl%vSplineC(1:n),DHA_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%DHA)
  ENDIF
  IF(Ok_Rho) THEN
    N= DHB_Spl%Dimm
    CALL CMM_Spline_Eval(n, &
    & DHB_Spl%vX(1:n), DHB_Spl%vY(1:n), &
    & DHB_Spl%vSplineB(1:n),DHB_Spl%vSplineC(1:n),DHB_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%DHB)
  ENDIF
  IF(Ok_BDot) THEN
    N= BDot_Spl%Dimm
    CALL CMM_Spline_Eval(n, &
    & BDot_Spl%vX(1:n), BDot_Spl%vY(1:n), &
    & BDot_Spl%vSplineB(1:n),BDot_Spl%vSplineC(1:n),BDot_Spl%vSplineD(1:n),&
    & TdgK-T_CK, SolModel%Dat%BDot)
  ENDIF
  
  IF(iDebug>2) THEN
    IF(Ok_Rho)  PRINT '(A,G15.6)',"Rho= ",SolModel%Dat%Rho
    IF(Ok_Eps)  PRINT '(A,G15.6)',"Eps= ",SolModel%Dat%Eps
    IF(Ok_DHA)  PRINT '(A,G15.6)',"DHA= ",SolModel%Dat%DHA
    IF(Ok_DHB)  PRINT '(A,G15.6)',"DHB= ",SolModel%Dat%DHB
    IF(Ok_BDot) PRINT '(A,G15.6)',"BDot=",SolModel%Dat%BDot
  ENDIF
  
  RETURN
END SUBROUTINE Solmodel_TP_Update_LogK

SUBROUTINE SolModel_TP_Update( &
& TdgK,Pbar, &
& SolModel)
!--
!-- For given (T,P) condition,
!-- compute SolModel properties using DtbH2OHkf_Calc and write in SolModel.
!-- If the database is of logK format,
!-- it may contains tabulated values for either Rho, Eps, etc.
!-- apply special procedure that will overwrite SolModel
!--
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_T_DtbH2OHkf,ONLY: DtbH2OHkf_Calc,T_DtbH2OHkf
  USE M_T_SolModel, ONLY: T_SolModel
  USE M_Solmodel_Pitzer_Dtb, ONLY: Solmodel_Pitzer_Dtb_TPUpdate
  
  USE M_Dtb_Vars,    ONLY: DtbFormat
  USE M_Solmodel_Vars,ONLY: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  USE M_Solmodel_Vars,ONLY: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl,T_Spline
    
  REAL(dp),        INTENT(IN)   :: TdgK,Pbar
  TYPE(T_SolModel),INTENT(INOUT):: SolModel
  
  TYPE(T_DtbH2OHkf)::PropsH2O
  
  !-- compute solvent (-H2O) properties at T,P
  CALL DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O)
  
  SolModel%Dat%Rho=  PropsH2O%Rho
  SolModel%Dat%Eps=  PropsH2O%Eps
  SolModel%Dat%dhA=  PropsH2O%dhA
  SolModel%Dat%dhB=  PropsH2O%dhB
  
  SolModel%Dat%bDot= Solmodel_BDot(TdgK)
  
  ! LOGK databases may contain specific values 
  ! for RHO, EPS, DHA, DHB, BDOT
  
  IF(DtbFormat=="LOGKTBL") &
  & CALL Solmodel_TP_Update_LogK( &
  & TdgK,Pbar, &
  & Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl, &
  & Ok_Rho, Ok_Eps, Ok_DHA, Ok_DHB, Ok_BDot, &
  & SolModel)
  
  IF(iDebug>2) THEN
    WRITE(fTrc,'(/,A,/)') "< SolModel_TP_Update"
    WRITE(fTrc,'(5(A,G15.6,/))') &
    & "Rho= ",SolModel%Dat%Rho, & 
    & "Eps= ",SolModel%Dat%Eps, &
    & "dhA= ",SolModel%Dat%dhA, &
    & "dhB= ",SolModel%Dat%dhB, &
    & "Bdot=",SolModel%Dat%bDot
    WRITE(fTrc,'(A,/)') "</ SolModel_TP_Update"
  ENDIF
  
  IF(SolModel%iActModel== 8  .OR. &    !! "PITZER"
  &  SolModel%iActModel== 11 ) THEN    !! "SIT"
    CALL Solmodel_Pitzer_Dtb_TPUpdate(TdgK,Pbar)
  ENDIF
    
  RETURN
ENDSUBROUTINE SolModel_TP_Update

REAL(dp) FUNCTION Solmodel_Bdot(TdgK)
  REAL(dp),INTENT(IN):: TdgK
  
  Solmodel_Bdot=  &
  & -0.616616D0  &
  & +6.90440D-3  *TdgK &
  & -2.73898D-5  *TdgK**2 &
  & +4.87313D-8  *TdgK**3 &
  & -3.26030D-11 *TdgK**4

END FUNCTION Solmodel_Bdot

ENDMODULE M_SolModel_Tools

