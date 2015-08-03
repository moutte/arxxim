MODULE M_Eos_H2O_Rho_Psat
!! was M_CalcPSatRho_H2O

!--------------------------------------------------------------
! Purpose : Compute PSatH2O and RhoH2O for Pure Water
!--------------------------------------------------------------
! PSatH2O : Saturation Pressure of Steam in Water
! RhoH2O  : Pure Water Density  
!
! Reference : deCapitani Theriak [fsol.f90 or gcalc.f90]
!--------------------------------------------------------------
! Arxim Integration : J.Moutte
!---------------------------------------------------------------

  USE M_Kinds
  IMPLICIT NONE
  PRIVATE

  !// Public Functions
  PUBLIC :: Eos_H2O_rho
  PUBLIC :: Eos_H2O_psat

CONTAINS

!---

SUBROUTINE Eos_H2O_psat(Tk,PSatBar) 
!--
!-- ... should give the *source* of these equations !!! ...
!--
  REAL(dp),INTENT(IN) :: Tk
  REAL(dp),INTENT(OUT):: PSatBar
  
  REAL(dp):: W,WSQ,V,FF
  INTEGER :: I
  REAL(dp):: A(1:8)=(/ &
  & -7.8889166D0, 2.5514255D0, -6.716169D0,  33.239495D0, &
  & -105.38479D0, 174.35319D0, -148.39348D0, 48.631602D0 /)
  !
  IF (Tk <= 314.0D0) THEN
  
    PSatBar= EXP( 6.3573118D0 -8858.843D0/Tk +607.56335D0/(Tk**0.6D0) )
    
  ELSE
  
    V=  Tk/647.25D0
    W=  ABS(One-V)
    WSQ=SQRT(W)
    FF= Zero
    DO I=1,8
      FF=FF+A(I)*W
      W=W*WSQ
    ENDDO
    !
    PSatBar= 220.93D0*EXP(FF/V)
  
  ENDIF
  !
  RETURN
ENDSUBROUTINE Eos_H2O_psat

!---

SUBROUTINE Eos_H2O_rho( &
& TdgK,Pbar, & !IN
& RH2O,Alfa,Beta,dAlfdT,dBetdT,dBetdP) !OUT
!-----------------------------------------------------------------------
!--retrieved from deCapitani / TheriakSuite, similar to routine in SupCrt
!--For given Pbar-Pbar,TdgK-TdegK,
!--compute density of water, RH2O, 
!--and Alfa, Beta, dAlfa/dT, dBeta/dT, dAlfa/dP_, dBeta/dP_
!--which are input DATA for :
!--__G_Shok92, to compute Gf, dG/dP_, dG/dT, d2G/dT2
!--__JohnsonNorton91, to compute epsilon, and its derivatives
!-----------------------------------------------------------------------
  USE M_Eos_Utils, ONLY : Cubic
  
  IMPLICIT NONE
  
  REAL(dp),INTENT(IN) ::Pbar,TdgK
  REAL(dp),INTENT(OUT)::RH2O,Alfa,Beta,dAlfdT,dBetdT,dBetdP
  
  INTEGER::KI(40),LI(40)
  REAL(dp)::&
  & TAUI(0:6),ERMI(0:9),&
  & GI(40), &
  & RHOI(37:40),TTTI(37:40),ALPI(37:40),BETI(37:40),&
  & TIT(7),PIT(7),RHIT(7),DPRIT(7)
  INTEGER:: ITER,LOO,I
  REAL(dp):: &
  & RR,R,T,T0,P,RT,PR,PS,S,&
  & ARK,BRK,OFT,   &
  & BUK,CUK,DUK,   &
  & DEL,DELT,DELP, &
  & DP_,DR_,       &
  & B,BB,          &
  & ALY,BETY,      &
  & X1,X2,X2I,X3,  &
  & Y,Y3,F1,F2,    &
  & VOL,RH,RH2,RHN,RHOI2, &
  & DPR,ER,    &
  & TAU,QHEX,Q10,QM
  !
  KI=(/&
    1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5, &
    6,6,6,6, 7,7,7,7, 9,9,9,9, 3,3,1,5, 2,2,2,4  &
      /)
  LI=(/&
    1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, &
    1,2,4,6, 1,2,4,6, 1,2,4,6, 0,3,3,3, 0,2,0,0  &
      /)
  !GI ARE IN (bar cc / g)=  10 * (J / g)
  GI= (/&
  & -.53062968529023D4,  .22744901424408D5, .78779333020687D4, &
  & -.69830527374994D3,  .17863832875422D6,-.39514731563338D6, &
  &  .33803884280753D6, -.13855050202703D6,-.25637436613260D7, &
  &  .48212575981415D7, -.34183016969660D7, .12223156417448D7, &
  &  .11797433655832D8, -.21734810110373D8, .10829952168620D8, &
  & -.25441998064049D7, -.31377774947767D8, .52911910757704D8, &
  & -.13802577177877D8, -.25109914369001D7, .46561826115608D8, &
  & -.72752773275387D8,  .41774246148294D7, .14016358244614D8, &
  & -.31555231392127D8,  .47929666384584D8, .40912664781209D7, &
  & -.13626369388386D8,  .69625220862664D7,-.10834900096447D8, &
  & -.22722827401688D7,  .38365486000660D7, .68833257944332D5, &
  &  .21757245522644D6, -.26627944829770D5,-.70730418082074D6, &
  & -.225D1,            -1.68D1,            .055D1,           -93.0D1 &
  &   /) 
  RHOI=(/0.319D0, 0.310D0, 0.310D0, 1.550D0 /)
  TTTI=(/640.0D0, 640.0D0, 641.6D0, 270.0D0 /)
  ALPI=(/34.0D0,  40.0D0,  30.0D0,  1050.0D0/)
  BETI=(/2.0D4,   2.0D4,   4.0D4,   25.0D0  /)
  !
  R= 4.6152D0
  !4.6152= 10 * 8.31441 / 18.0152 = R in bar.cc/ gramH2O /K
  RR=8.31441D0 !R in J/mole/K
  !
  T=TdgK
  P=Pbar
  RT=R*T
  !
  PS=220.55D0
  !
  IF (T<=647.25D0) CALL Eos_H2O_psat(T,PS) !for subcritic. water
  !
  !SET INITIAL GUESS FOR RHO USING THB.HOLLAND-FIT TO REDLICH-KWONG
  ARK= 1.279186D8 -2.241415D04*T !REDLICH-KWONG constant A
  BRK= 1.428062D1 +6.092237D-4*T !REDLICH-KWONG constant B
  OFT= ARK/(P*SQRT(T))
  BUK=-10D0*RR*T/P
  CUK= OFT -BRK*BRK +BRK*BUK
  DUK=     -BRK*OFT
  !
  CALL CUBIC(BUK,CUK,DUK,X1,X2,X2I,X3)
  !
  IF (X2I /= 0.0D0) THEN
    VOL=X1
  ELSE
    !IF (P<PS) THEN; VOL=MAX(X1,X2,X3)
    !ELSE;           VOL=MIN(X1,X2,X3)
    !END IF
    VOL=DMIN1(X1,X2,X3);
  END IF
  IF (VOL<=0.0D0) THEN; RHN=1.9D0
  ELSE                ; RHN=1D0/VOL*18.0152D0
  ENDIF
  !
  DELT=0.001D0; DELP=0.01D0
  TIT(1)=TdgK;             PIT(1)=Pbar
  TIT(2)=TdgK -DELT;       PIT(2)=Pbar
  TIT(3)=TdgK +DELT;       PIT(3)=Pbar
  TIT(4)=TdgK -DELT*2.0D0; PIT(4)=Pbar
  TIT(5)=TdgK +DELT*2.0D0; PIT(5)=Pbar
  TIT(6)=TdgK;             PIT(6)=Pbar-DELP
  TIT(7)=TdgK;             PIT(7)=Pbar+DELP
  T0=    647.073D0
  !
  DO100: DO ITER=1,7
  
    IF (ITER /= 1) RHN=RHIT(1)
    
    T=TIT(ITER)
    P=PIT(ITER)
    RT=R*T
    
    !The values (T/T0)**i are stored in the array TAUI(i)
    TAUI(0)=1.D0
    TAUI(1)=T/T0
    DO I=2,6; TAUI(I)=TAUI(I-1)*TAUI(1); ENDDO
    B = -0.3540782D0 *LOG(TAUI(1)) &
    & + 0.7478629D0 &
    & + 0.007159876D0/TAUI(3) &
    & - 0.003528426D0/TAUI(5)
    BB= 1.1278334D0  &
    & - 0.5944001D0 /TAUI(1) &
    & - 5.010996D0  /TAUI(2)  &
    & + 0.63684256D0/TAUI(4)
    
    !FIND THE TRUE(?) RH(T,P)
    !NOTE: PR= PRESSURE CORRESPONDING TO GUESSED RH
    !      DPR= (dP_ / dRH)
    !      the values (1-EXP(-RH))**i are stored in the array ERMI(i)
    DO20: DO LOO=1,100
    
      RH=RHN
      IF (RH<=1D-8) RH=1D-8
      IF (RH> 1.9D0) RH=1.9D0
      RH2= RH*RH
      Y=   RH*B/4.D0
      ER=  EXP(-RH)
      Y3=  (1D0-Y)**3
      ALY= 11.D0*Y
      BETY=44.33333333333333D0*Y*Y
      F1= (1.D0+ALY+BETY)/Y3
      F2= 4.D0*Y*(BB/B-3.5D0)
      ERMI(0)=1D0
      ERMI(1)=1D0-ER
      DO I=2,9; ERMI(I)=ERMI(I-1)*ERMI(1); ENDDO
      PR=0.0D0
      DPR=0.0D0
      DO I=1,36
        S=  GI(I)/TAUI(LI(I))*ERMI(KI(I)-1)
        PR= PR+S
        DPR=DPR+(2D0+RH*(KI(I)*ER-1D0)/ERMI(1))*S
      END DO
      DO I=37,40
        DEL=  RH/RHOI(I)-One
        RHOI2=RHOI(I)*RHOI(I)
        TAU=  T/TTTI(I)-One
        QHEX=(-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
        IF (QHEX>-150.0D0) THEN
          Q10= GI(I) *DEL**LI(I) &
          &          *EXP( -ALPI(I)*DEL**LI(I) -BETI(I)*TAU*TAU )
        ELSE
          Q10=0.0D0
        END IF
        QM= LI(I)/DEL-KI(I)*ALPI(I)*DEL**(KI(I)-1)
        S= Q10*QM*RH2/RHOI(I)
        PR= PR+S
        DPR= DPR &
        &  + S*(2.0D0/RH+QM/RHOI(I)) &
        &  - RH2 /RHOI2 *Q10 &
        &    *( LI(I)/DEL/DEL +KI(I)*(KI(I)-1)*ALPI(I)*DEL**(KI(I)-2) )
      END DO
      PR=RH*(RH*ER*PR+RT*(F1+F2))
      DPR= RH*ER*DPR &
      &  + RT * ((1.D0 +2.D0*ALY +3D0*BETY)/Y3 &
      &         + 3.D0*Y*F1/(1D0-Y) &
      &         + 2.D0*F2)
      !
      IF (DPR<=0.0D0) THEN
        IF (P<=PS) THEN; RHN=RHN*0.95D0
        ELSE           ; RHN=RHN*1.05D0
        ENDIF
      ELSE
        IF(DPR<0.01D0) DPR=0.01D0
        S=(P-PR)/DPR
        IF(ABS(S)>0.1D0) S=0.1D0*S/ABS(S)
        RHN=RH+S
      END IF
      DP_=ABS(One-PR/P)
      DR_=ABS(One-RHN/RH)
      IF (DP_<1D-5 .AND. DR_<1D-12) EXIT DO20 !GOTO 30
    
    ENDDO DO20 !20 CONTINUE
    
    RHIT(ITER)=RHN !30
    DPRIT(ITER)=DPR
  
  ENDDO DO100
  !
  RH2O=RHIT(1)
  !
  Beta=  One/DPRIT(1)/RHIT(1)
  F1=    One/DPRIT(6)/RHIT(6)
  F2=    One/DPRIT(7)/RHIT(7)
  !
  dBetdP=(F2-F1)/DELP/2.0D0
  !
  F1=    One/DPRIT(2)/RHIT(2)
  F2=    One/DPRIT(3)/RHIT(3)
  !
  dBetdT=(F2-F1)/DELT/2.0D0
  !
  Alfa= -(RHIT(3)-RHIT(2))/DELT/2.0D0/RHIT(1)
  !
  F1=   -(RHIT(1)-RHIT(4))/DELT/2.0D0/RHIT(2)
  F2=   -(RHIT(5)-RHIT(1))/DELT/2.0D0/RHIT(3)
  !
  dAlfdT=(F2-F1)/DELT/2.0D0
  
  RETURN
ENDSUBROUTINE Eos_H2O_rho

ENDMODULE M_Eos_H2O_Rho_Psat
