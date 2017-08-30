module M_Eos_H2O_Rho_Psat
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

  use M_Kinds
  implicit none
  private

  !// Public Functions
  public :: Eos_H2O_rho
  public :: Eos_H2O_psat

contains

!---

subroutine Eos_H2O_psat(Tk,PSatBar) 
!--
!-- ... should give the *source* of these equations !!! ...
!--
  real(dp),intent(in) :: Tk
  real(dp),intent(out):: PSatBar
  
  real(dp):: W,WSQ,V,FF
  integer :: I
  real(dp):: A(1:8)=(/ &
  & -7.8889166D0, 2.5514255D0, -6.716169D0,  33.239495D0, &
  & -105.38479D0, 174.35319D0, -148.39348D0, 48.631602D0 /)
  !
  if (Tk <= 314.0D0) then
  
    PSatBar= exp( 6.3573118D0 -8858.843D0/Tk +607.56335D0/(Tk**0.6D0) )
    
  else
  
    V=  Tk/647.25D0
    W=  ABS(One-V)
    WSQ=SQRT(W)
    FF= Zero
    do I=1,8
      FF=FF+A(I)*W
      W=W*WSQ
    enddo
    !
    PSatBar= 220.93D0*exp(FF/V)
  
  end if
  !
  return
end subroutine Eos_H2O_psat

!---

subroutine Eos_H2O_rho( &
& TdgK,Pbar, & !IN
& RH2O,Alfa,Beta,dAlfdT,dBetdT,dBetdP) !OUT
!-----------------------------------------------------------------------
!--retrieved from deCapitani / TheriakSuite, similar to routine in SupCrt
!--For given Pbar-Pbar,TdgK-TdegK,
!--compute density of water, RH2O, 
!--and Alfa, Beta, dAlfa/dT, dBeta/dT, dAlfa/dP_, dBeta/dP_
!--which are input data for :
!--__G_Shok92, to compute Gf, dG/dP_, dG/dT, d2G/dT2
!--__JohnsonNorton91, to compute epsilon, and its derivatives
!-----------------------------------------------------------------------
  use M_Eos_Utils, only : Cubic
  
  implicit none
  
  real(dp),intent(in) ::Pbar,TdgK
  real(dp),intent(out)::RH2O,Alfa,Beta,dAlfdT,dBetdT,dBetdP
  
  integer::KI(40),LI(40)
  real(dp)::&
  & TAUI(0:6),ERMI(0:9),&
  & GI(40), &
  & RHOI(37:40),TTTI(37:40),ALPI(37:40),BETI(37:40),&
  & TIT(7),PIT(7),RHIT(7),DPRIT(7)
  integer:: ITER,LOO,I
  real(dp):: &
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
  if (T<=647.25D0) call Eos_H2O_psat(T,PS) !for subcritic. water
  !
  !SET INITIAL GUESS FOR RHO USING THB.HOLLAND-FIT TO REDLICH-KWONG
  ARK= 1.279186D8 -2.241415D04*T !REDLICH-KWONG constant A
  BRK= 1.428062D1 +6.092237D-4*T !REDLICH-KWONG constant B
  OFT= ARK/(P*SQRT(T))
  BUK=-10D0*RR*T/P
  CUK= OFT -BRK*BRK +BRK*BUK
  DUK=     -BRK*OFT
  !
  call CUBIC(BUK,CUK,DUK,X1,X2,X2I,X3)
  !
  if (X2I /= 0.0D0) then
    VOL=X1
  else
    !if (P<PS) then; VOL=MAX(X1,X2,X3)
    !else;           VOL=MIN(X1,X2,X3)
    !end if
    VOL=DMIN1(X1,X2,X3);
  end if
  if (VOL<=0.0D0) then; RHN=1.9D0
  else                ; RHN=1D0/VOL*18.0152D0
  end if
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
  do100: do ITER=1,7
  
    if (ITER /= 1) RHN=RHIT(1)
    
    T=TIT(ITER)
    P=PIT(ITER)
    RT=R*T
    
    !The values (T/T0)**i are stored in the array TAUI(i)
    TAUI(0)=1.D0
    TAUI(1)=T/T0
    do I=2,6; TAUI(I)=TAUI(I-1)*TAUI(1); enddo
    B = -0.3540782D0 *log(TAUI(1)) &
    & + 0.7478629D0 &
    & + 0.007159876D0/TAUI(3) &
    & - 0.003528426D0/TAUI(5)
    BB= 1.1278334D0  &
    & - 0.5944001D0 /TAUI(1) &
    & - 5.010996D0  /TAUI(2)  &
    & + 0.63684256D0/TAUI(4)
    
    !FIND THE.true.?) RH(T,P)
    !NOTE: PR= PRESSURE CORRESPONDING TO GUESSED RH
    !      DPR= (dP_ / dRH)
    !      the values (1-exp(-RH))**i are stored in the array ERMI(i)
    do20: do LOO=1,100
    
      RH=RHN
      if (RH<=1D-8) RH=1D-8
      if (RH> 1.9D0) RH=1.9D0
      RH2= RH*RH
      Y=   RH*B/4.D0
      ER=  exp(-RH)
      Y3=  (1D0-Y)**3
      ALY= 11.D0*Y
      BETY=44.33333333333333D0*Y*Y
      F1= (1.D0+ALY+BETY)/Y3
      F2= 4.D0*Y*(BB/B-3.5D0)
      ERMI(0)=1D0
      ERMI(1)=1D0-ER
      do I=2,9; ERMI(I)=ERMI(I-1)*ERMI(1); enddo
      PR=0.0D0
      DPR=0.0D0
      do I=1,36
        S=  GI(I)/TAUI(LI(I))*ERMI(KI(I)-1)
        PR= PR+S
        DPR=DPR+(2D0+RH*(KI(I)*ER-1D0)/ERMI(1))*S
      end do
      do I=37,40
        DEL=  RH/RHOI(I)-One
        RHOI2=RHOI(I)*RHOI(I)
        TAU=  T/TTTI(I)-One
        QHEX=(-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
        if (QHEX>-150.0D0) then
          Q10= GI(I) *DEL**LI(I) &
          &          *exp( -ALPI(I)*DEL**LI(I) -BETI(I)*TAU*TAU )
        else
          Q10=0.0D0
        end if
        QM= LI(I)/DEL-KI(I)*ALPI(I)*DEL**(KI(I)-1)
        S= Q10*QM*RH2/RHOI(I)
        PR= PR+S
        DPR= DPR &
        &  + S*(2.0D0/RH+QM/RHOI(I)) &
        &  - RH2 /RHOI2 *Q10 &
        &    *( LI(I)/DEL/DEL +KI(I)*(KI(I)-1)*ALPI(I)*DEL**(KI(I)-2) )
      end do
      PR=RH*(RH*ER*PR+RT*(F1+F2))
      DPR= RH*ER*DPR &
      &  + RT * ((1.D0 +2.D0*ALY +3D0*BETY)/Y3 &
      &         + 3.D0*Y*F1/(1D0-Y) &
      &         + 2.D0*F2)
      !
      if (DPR<=0.0D0) then
        if (P<=PS) then; RHN=RHN*0.95D0
        else           ; RHN=RHN*1.05D0
        end if
      else
        if(DPR<0.01D0) DPR=0.01D0
        S=(P-PR)/DPR
        if(ABS(S)>0.1D0) S=0.1D0*S/ABS(S)
        RHN=RH+S
      end if
      DP_=ABS(One-PR/P)
      DR_=ABS(One-RHN/RH)
      if (DP_<1D-5 .and. DR_<1D-12) exit do20 !goto 30
    
    enddo do20 !20 continue
    
    RHIT(ITER)=RHN !30
    DPRIT(ITER)=DPR
  
  enddo do100
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
  
  return
end subroutine Eos_H2O_rho

end module M_Eos_H2O_Rho_Psat
