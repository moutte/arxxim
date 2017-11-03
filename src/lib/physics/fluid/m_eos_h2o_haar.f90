module M_Eos_H2O_Haar
!! was M_Haar_H2O

!-------------------------------------------------------------------------------
! HAAR      : Eos_H2O_Haar
!-------------------------------------------------------------------------------
! Purpose : Haar and al EOS to compute Gibbs Free energy and Volume of H2O
! Reference : Haar-Gallager-Kell, '84
!
! Source : CdeCapitani Theriak [ fsol.f90 ]
!
!-------------------------------------------------------------------------------
! Arxim Integration : J.Moutte
!-------------------------------------------------------------------------------

use M_Kinds
implicit none
private

!real(dp), parameter :: R_jK= 8.314510D0

!// public functionS
public:: Eos_H2O_Haar !! CalcGH2O_Haar

contains

subroutine Eos_H2O_Haar(TdgK,Pbar,G_H2Ort,V_H2O_m3)
!--
!-- CdeCapitani, SupCrt
!-- compute Gibbs free energy of Water G_H2O, using EoS (Haar-Gallager-Kell, '84)
!-- same as used in Eos_H2O_Rho
!--
  use M_Dtb_Const,only: R_jK, Tref,DtbConv_Benson,S0_Hydrogen,S0_Oxygen
  use M_Eos_H2O_Rho_Psat
  use M_Eos_Utils
  !
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(out):: G_H2Ort,V_H2O_m3
  !
  real(dp):: &
    G_H2O,   &
    TAUI(0:6),ERMI(0:9), &
    P, T, R, RR, T0, RT, S, &
    B, BB, ARK, BRK, OFT, &
    cubB, cubC, cubD, X1,X2,X2I,X3, &
    PS, PR, DPR, GREF, &
    Y3,F1,F2,&
    VOL, RH, RHN, RH2, ER, DEL, TAU, &
    ALY, BETY, QHEX, Q10, QM, RHOI2, &
    DP_, DR_,  X, Y, AA, TR, W, AID
  integer:: LOO, I
  !
  real(dp)::RHOI(37:40)=(/0.319D0, 0.310D0, 0.310D0,  1.550D0  /)
  real(dp)::TTTI(37:40)=(/640.0D0, 640.0D0, 641.6D0,  270.0D0  /)
  real(dp)::ALPI(37:40)=(/34.0D0,  40.0D0,  30.0D0,   1050.0D0 /)
  real(dp)::BETI(37:40)=(/2.0D4,   2.0D4,   4.0D4,    25.0D0   /)
  !
  integer::KI(1:40)=(/&
  &  1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5, &
  &  6,6,6,6, 7,7,7,7, 9,9,9,9, 3,3,1,5, 2,2,2,4  /)
  integer::LI(1:40)=(/&
  &  1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, &
  &  1,2,4,6, 1,2,4,6, 1,2,4,6, 0,3,3,3, 0,2,0,0  /)
  real(dp)::CI(1:18)=(/&
  &  .19730271018D2,     .209662681977D2,  -.483429455355D0, &
  &  .605743189245D1,  22.56023885D0,      -9.87532442D0,    -.43135538513D1,&
  &  .458155781D0,      -.47754901883D-1,   .41238460633D-2, -.27929052852D-3,&
  &  .14481695261D-4,   -.56473658748D-6,   .16200446D-7,    -.3303822796D-9,&
  &  .451916067368D-11, -.370734122708D-13, .137546068238D-15 /)
  !
  !GI ARE IN (bar cc / g)=    10 * (J / g)
  real(dp)::GI(1:40)=(/&
  &  -.53062968529023D4,  .22744901424408D5, .78779333020687D4,&
  &  -.69830527374994D3,  .17863832875422D6,-.39514731563338D6,&
  &   .33803884280753D6, -.13855050202703D6,-.25637436613260D7,&
  &   .48212575981415D7, -.34183016969660D7, .12223156417448D7,&
  &   .11797433655832D8, -.21734810110373D8, .10829952168620D8,&
  &  -.25441998064049D7, -.31377774947767D8, .52911910757704D8,&
  &  -.13802577177877D8, -.25109914369001D7, .46561826115608D8,&
  &  -.72752773275387D8,  .41774246148294D7, .14016358244614D8,&
  &  -.31555231392127D8,  .47929666384584D8, .40912664781209D7,&
  &  -.13626369388386D8,  .69625220862664D7,-.10834900096447D8,&
  &  -.22722827401688D7,  .38365486000660D7, .68833257944332D5,&
  &   .21757245522644D6, -.26627944829770D5,-.70730418082074D6,&
  &  -.225D1,            -1.68D1,            .055D1,           &
  &  -93.0D1 /)
  !
  T=TdgK
  P=Pbar
  !
  RR= 8.31441D0
  ! supposed to be R_JK ??
  ! Nota: in RobieHemingway (from CODATA), R_jK= 8.314510D0 J/K/Mole
  ! R=  4.6152D0  =10 * 8.31441 / 18.0152 = R in bar.cc/g/K

  R= 4.61522D0 != value in code supcrt92/h2O92
  RT= R*T
  !GREF=-729.0855D0*18.0152D0*4.184D0
  !-> conversion cal/gram_H2O -> J/Mole

  !--- GREF CALCULATED WITH THIS ROUTINE AT 1 BAR AND 25 DEG. C
  GREF= -54955.2356146119409D0
  ! (value used in Haar_Ghiorso : gref= -54955.23970014D0)

  T0= 647.073D0
  ! The values (T/T0)**i are stored in the array TAUI(i)
  TAUI(0)=One
  TAUI(1)=T/T0
  do I=2,6
    TAUI(I)=TAUI(I-1)*TAUI(1)
  end do
  !
  B= - 0.3540782D0 *log(TAUI(1)) &
     + 0.7478629D0 &
     + 0.007159876D0 /TAUI(3) &
     - 0.003528426D0 /TAUI(5)
  BB=  1.1278334D0 &
     - 0.5944001D0  /TAUI(1) &
     - 5.010996D0   /TAUI(2) &
     + 0.63684256D0 /TAUI(4)
  !
  !NAME='STEAM'
  PS= 220.55D0
  if (T<=647.25D0) call Eos_H2O_psat(T,PS)    !if (P >= PS) NAME='STEAM'
  !
  !-- SET INITIAL GUESS FOR RHO USING THB-FIT TO REDLICH-KWONG
  ARK=  1.279186D8-2.241415D4 *T   !RedlichKwong Param. A
  BRK=  1.428062D1+6.092237D-4*T   !RedlichKwong Param. B
  OFT=  ARK/(P*sqrt(T))
  cubB= -10.0D0*RR*T/P
  cubC= OFT -BRK*BRK +BRK*cubB
  cubD= -BRK*OFT
  !
  call CUBIC(cubB,cubC,cubD,X1,X2,X2I,X3)
  !
  if(X2I/=Zero) then
    VOL=X1
  else
    if (P<PS)  then ; VOL=DMAX1(X1,X2,X3)
    else            ; VOL=DMIN1(X1,X2,X3)
    end if
  end if
  if (VOL<=Zero) then  ; RHN=1.9D0
  else                 ; RHN=One/VOL*18.0152D0
  end if
  !-- FIND THE TRUE(?) RH(T,P)
  !-- NOTE: PR-  PRESSURE CORRESPONDING TO GUESSED RH
  !-- DPR-  (dP_ / dRH)
  !-- the values (1-exp(-RH))**i are stored in the array ERMI(i)
  do LOO=1,100
    RH=RHN
    if(RH<=1D-8)  RH=1D-8
    if(RH> 1.9D0) RH=1.9D0
    RH2= RH*RH
    Y=   RH*B/4.D0
    ER=  exp(-RH)
    Y3=  (One-Y)**3
    ALY= 11.D0*Y
    BETY=133.D0/3.D0*Y*Y
    F1=  (One+ALY+BETY)/Y3
    F2=  4.D0*Y*(BB/B-3.5D0)
    !
    ERMI(0)=One
    ERMI(1)=One-ER
    do I=2,9; ERMI(I)=ERMI(I-1)*ERMI(1); end do
    !
    PR=Zero
    DPR=Zero
    do I=1,36
      S=    GI(I)/TAUI(LI(I))*ERMI(KI(I)-1)
      PR=   PR+S
      DPR=  DPR+(2D0+RH*(KI(I)*ER-One)/ERMI(1))*S
    end do
    !
    do I=37,40
      DEL=   RH/RHOI(I)-One
      RHOI2= RHOI(I)*RHOI(I)
      TAU=   T/TTTI(I)-One
      QHEX= (-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
      if (QHEX>-150.0D0) then
        Q10= GI(I) *DEL**LI(I) *exp(-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
      else
        Q10= Zero
      end if
      QM=  LI(I)/DEL-KI(I)*ALPI(I)*DEL**(KI(I)-1)
      S=   Q10*QM*RH2/RHOI(I)
      PR=  PR+S
      DPR= DPR &
         + S*(Two/RH +QM/RHOI(I)) &
         - RH2/RHOI2*Q10 &
           *(LI(I)/DEL/DEL+KI(I)*(KI(I)-1)*ALPI(I)*DEL**(KI(I)-2))
    end do
    PR= RH*(RH*ER*PR+RT*(F1+F2))
    DPR=RH*ER*DPR+RT*((One+Two*ALY+3.D0*BETY)/Y3 +3.D0*Y*F1/(One-Y)+2.D0*F2)
    if (DPR <= Zero) then
      if (P<=PS) then; RHN=RHN*0.95D0
      else;            RHN=RHN*1.05D0
      end if
    else
      if (DPR < 0.01D0) DPR=0.01D0
      S=(P-PR)/DPR
      if (abs(S) > 0.1D0) S=0.1D0*S/abs(S)
      RHN=RH+S
    end if
    DP_=abs(One-PR/P)
    DR_=abs(One-RHN/RH)
    if (DP_ < 1.D-5.and.DR_ < 1.D-5) exit
  end do
  !
  RH=     RHN
  Y=      RH*B/4.D0
  X=      One-Y
  ER=     exp(-RH)
  ERMI(0)=One
  ERMI(1)=One-ER
  do I=2,9; ERMI(I)=ERMI(I-1)*ERMI(1); end do
  !
  !--- base function --
  AA= RT*(-log(X) &
    - 130.0D0/3.0D0 /X     & !
    + 169.0D0/6.0D0 /X/X   & !
    + 4.D0*Y*(BB/B-3.5D0)  &
    + 91.D0/6.D0           & !
    + log(RH*RT/1.01325D0))
  !
  !--- residual function --
  do I=1,36
    AA= AA &
    & + GI(I) /KI(I) /TAUI(LI(I)) *ERMI(KI(I))
  end do
  do I=37,40
    DEL=  RH/RHOI(I)-One
    TAU=  T/TTTI(I)-One
    QHEX= (-ALPI(I)*DEL**KI(I)-BETI(I)*TAU*TAU)
    if (QHEX > -150.0D0) AA=AA+GI(I)*DEL**LI(I)*exp(QHEX)
  end do
  !
  !--- ideal gas function --
  TR=  T/1.0D2
  W=   TR**(-3)
  AID= One +( CI(1)/TR +CI(2) )*log(TR)
  do I=3,18
    AID= AID +CI(I)*W
    W=   W*TR
  end do
  !
  AA= AA -RT*AID
  !-- CALCULATE G-  AA + P/RH  AND  V-  1/RH
  !-- G CORRECTED TO DH-T*S FOR WATER AT 1 BAR AND 25 DEG. C
  !-- ACCORDING TO ROBIE ET Al. 1978
  !
  V_H2O_m3= 18.0152D0 /RH *1.0D-6  !molar volume in m3
  !
  G_H2O= (AA+P/RH)*1.80152D0 - GREF -306685.5925D0
  !*1.80152D0=  conversion from Joule/Gram to Joule/mole
  !
  !Robie/codata data for entropies of elements at 298.15:
  !! H2: 130.68 J/K/mole
  !! O2: 205.15 J/K/mole
  !! C:    5.74 J/K/mole
  !! -> for H2O, to transform G(BermBrown) to G(BensHelg),
  !!____add 298.15*(130.68 + 205.15/2)=  69544.98 J
  !
  if(DtbConv_Benson) G_H2O=   G_H2O + Tref *(S0_Hydrogen *Two + S0_Oxygen)
  G_H2Ort= G_H2O /R_jK /TdgK
  !
  return
end subroutine Eos_H2O_Haar

subroutine WHAAR2(P,T,Grt,VH2O)
!-- for tests, not used !!
!--
  use M_Dtb_Const,only: R_JK,Tref,DtbConv_Benson,S0_Hydrogen,S0_Oxygen
  !
  real(dp),intent(in)    :: P,T
  real(dp),intent(out)   :: Grt,VH2O
  !
  character(len=16):: NAME
  !
  real(dp):: TAUI(0:6),ERMI(0:9),GI(40),CI(18)
  integer :: KI(40),LI(40)
  real(dp):: RHOI(37:40),TTTI(37:40),ALPI(37:40),BETI(37:40)
  real(dp):: GH2O,R,RR,RT,GREF,T0,PS
  real(dp):: B,BB,ARK,BRK,X,Y,Y3,VOL,RHN,RH,RH2,F1,F2,PR,DPR,ER,S,AA
  real(dp):: RHOI2,TAU,DEL,QHEX,Q10,QM,DP_,DR,TR,W,AID
  real(dp):: OFT,BUK,CUK,DUK,X1,X2,X2I,X3
  real(dp):: ALY,BETY
  integer:: I,LOO
  !
  !-- GI ARE IN (bar cc / g)  -  10 * (J / g)
  data GI/-.53062968529023D4,.22744901424408D5,.78779333020687D4,&
          -.69830527374994D3,.17863832875422D6,-.39514731563338D6,&
          .33803884280753D6,-.13855050202703D6,-.25637436613260D7,&
          .48212575981415D7,-.34183016969660D7, .12223156417448D7,&
          .11797433655832D8,-.21734810110373D8, .10829952168620D8,&
          -.25441998064049D7,-.31377774947767D8,.52911910757704D8,&
          -.13802577177877D8,-.25109914369001D7, .46561826115608D8,&
          -.72752773275387D8,.41774246148294D7,.14016358244614D8,&
          -.31555231392127D8,.47929666384584D8,.40912664781209D7,&
          -.13626369388386D8, .69625220862664D7,-.10834900096447D8,&
          -.22722827401688D7,.38365486000660D7,.68833257944332D5,&
          .21757245522644D6,-.26627944829770D5,-.70730418082074D6,&
          -.225D1,-1.68D1,.055D1,-93.0D1/
  data KI/4*1,4*2,4*3,4*4,4*5,4*6,4*7,4*9,2*3,1,5,3*2,4/
  data LI/1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,0,3*3,0,2,0,0/
  data CI/.19730271018D2,    .209662681977D2,    -.483429455355D0,&
          .605743189245D1,   22.56023885D0,      -9.87532442D0,    -.43135538513D1,&
          .458155781D0,-.47754901883D-1,.41238460633D-2,-.27929052852D-3,&
          .14481695261D-4,-.56473658748D-6,.16200446D-7,-.3303822796D-9,&
          .451916067368D-11,-.370734122708D-13,.137546068238D-15/
  data RHOI/0.319D0,0.310D0,0.310D0,1.550D0/
  data TTTI/640.0D0,640.0D0,641.6D0,270.0D0/
  data ALPI/34.0D0,40.0D0,30.0D0,1050.0D0/
  data BETI/2.0D4,2.0D4,4.0D4,25.0D0/
  !
  !!R=  4.6152D0
  R= 4.61522D0 != value in code supcrt92/h2O92
  RT= R*T
  ! GREF=-729.0855D0*18.0152D0*4.184D0
  !     = GREF CALCULATED WITH THIS ROUTINE AT 1 BAR AND 25 DEG. C
  GREF= -54955.2356146119409D0
  T0=    647.073D0
  !
  !The values (T/T0)**i are stored in the array TAUI(i)
  TAUI(0)=1.D0
  TAUI(1)=T/T0
  do I=2,6
    TAUI(I)=TAUI(I-1)*TAUI(1)
  end do
  !
  B  = -0.3540782D0*Dlog(TAUI(1)) &
       +0.7478629D0 &
       +0.007159876D0 /TAUI(3) &
       -0.003528426D0 /TAUI(5)
  BB = 1.1278334D0 &
       -0.5944001D0  /TAUI(1) &
       -5.010996D0   /TAUI(2) &
       +0.63684256D0 /TAUI(4)
  !
  !!NAME='STEAM'
  PS=220.55D0
  if (T <= 647.25D0) then
     call PSAT2(T,PS)
     if (P >= PS) NAME='STEAM'
  end if

  !-- SET INITIAL GUESS FOR RHO USING THB-FIT TO REDLICH-KWONG
  ARK= 1.279186D8 -2.241415D4*T
  BRK= 1.428062D1 +6.092237D-4*T
  RR=  8.31441D0
  OFT= ARK/(P*sqrt(T))
  BUK=-10D0*RR*T/P
  CUK=OFT-BRK*BRK+BRK*BUK
  DUK=-BRK*OFT

  call CUBIC(BUK,CUK,DUK,X1,X2,X2I,X3)

  if (X2I /= 0.0D0) then
    VOL=X1
  else
    if (P < PS)  then
      VOL=DMAX1(X1,X2,X3)
    else
      VOL=DMIN1(X1,X2,X3)
    end if
  end if

  if (VOL <= 0.0D0) then
    RHN=1.9D0
  else
    RHN=1D0/VOL*18.0152D0
  end if

  !-- FIND THE TRUE(?) RH(T,P)
  !-- NOTE: PR - PRESSURE CORRESPONDING TO GUESSED RH
  !-- DPR - (DP_ / dRH)
  !-- the values (1-exp(-RH))**i are stored in the array ERMI(i)
  do LOO=1,100

    RH=RHN
    if (RH <= 1D-8) RH=1D-8
    if (RH > 1.9D0) RH=1.9D0

    RH2=RH*RH
    Y=RH*B/4.D0
    ER=Dexp(-RH)
    Y3=(1D0-Y)**3
    ALY=11.D0*Y
    BETY=44.33333333333333D0*Y*Y

    F1=(1.D0+ALY+BETY)/Y3
    F2=4.D0*Y*(BB/B-3.5D0)
    !
    ERMI(0)=1D0
    ERMI(1)=1D0-ER
    do I=2,9
      ERMI(I)=ERMI(I-1)*ERMI(1)
    end do
    !
    PR=  0.0D0
    DPR= 0.0D0
    !
    do I=1,36
      S=GI(I)/TAUI(LI(I))*ERMI(KI(I)-1)
      PR=PR+S
      DPR=DPR+(2D0+RH*(KI(I)*ER-1D0)/ERMI(1))*S
    end do
    !
    do I=37,40
      DEL=RH/RHOI(I)-1.0D0
      RHOI2=RHOI(I)*RHOI(I)
      TAU=T/TTTI(I)-1.0D0
      QHEX=(-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
      if (QHEX > -150.0D0) then
        Q10=GI(I)*DEL**LI(I)*Dexp(-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
      else
        Q10=0.0D0
      end if
      QM=  LI(I)/DEL-KI(I)*ALPI(I)*DEL**(KI(I)-1)
      S=   Q10*QM*RH2/RHOI(I)
      PR=  PR+S
      DPR= DPR &
      &  + S*( 2.0D0/RH + QM/RHOI(I) ) &
      &  - RH2/RHOI2*Q10 &
      &    *( LI(I)/DEL/DEL +KI(I)*(KI(I)-1)*ALPI(I)*DEL**(KI(I)-2) )
    end do
    !
    PR=  RH*(RH*ER*PR+RT*(F1+F2))
    DPR= RH*ER*DPR &
    &  + RT*( (1.D0  + 2.D0*ALY + 3.D0*BETY)/Y3 &
    &       + 3.D0*Y*F1/(1.D0-Y) &
    &       + 2.D0*F2)
    !
    if (DPR <= 0.0D0) then
      if (P <= PS) then
        RHN=RHN*0.95D0
      else
        RHN=RHN*1.05D0
      end if
    else
      if (DPR < 0.01D0) DPR=0.01D0
      S=(P-PR)/DPR
      if (abs(S) > 0.1D0) S=0.1D0*S/abs(S)
      RHN=RH+S
    end if
    DP_= abs(1.0D0-PR/P)
    DR=  abs(1.0D0-RHN/RH)

    if (DP_ < 1D-5.and.DR < 1D-5) exit !goto 30

  end do !20

  RH=RHN
  Y=RH*B/4.D0
  X=1D0-Y
  ER=Dexp(-RH)
  ERMI(0)=1D0
  ERMI(1)=1D0-ER
  !
  do I=2,9
    ERMI(I)=ERMI(I-1)*ERMI(1)
  end do
  !
  !--- CALCULATE BASE function
  AA= RT*(-Dlog(X) &
     -43.33333333333333D0/X &
     +28.16666666666667D0/X/X &
     +4D0*Y*(BB/B-3.5D0) &
     +15.16666666666667D0 &
     +Dlog(RH*RT/1.01325D0))
  !
  !--- CALCULATE RESIDUAL function
  do I=1,36
    AA=AA+GI(I)/KI(I)/TAUI(LI(I))*ERMI(KI(I))
  end do
  do I=37,40
    DEL=RH/RHOI(I)-1.0D0
    TAU=T/TTTI(I)-1.0D0
    QHEX=(-ALPI(I)*DEL**KI(I)-BETI(I)*TAU*TAU)
    if (QHEX > -150.0D0) then
      AA=AA+GI(I)*DEL**LI(I)*Dexp(QHEX)
    end if
  end do
  !---/

  !--- CALCULATE IDEAL GAS function (-AID, HELMHOLTZ function, IDEAL)
  TR=  T/1.0D2
  W=   TR**(-3)
  AID= 1.D0 + (CI(1)/TR+CI(2))*Dlog(TR)
  do I=3,18
    AID= AID +CI(I)*W
    W=   W*TR
  end do
  !---/
  !
  AA=AA-RT*AID
  !
  !--- CALCULATE G - AA + P/RH  AND  V - 1/RH
  !--  G CORRECTED TO DH-T*S FOR WATER AT 1 BAR AND 25 DEG. C
  !--  ACCORDING TO ROBIE ET AL. 1978
  GH2O= (AA+P/RH)*1.80152D0 -GREF -306685.5925D0
  VH2O= 1/RH
  !
  !Robie/COdata data for entropies of elements at 298.15:
  !! H2: 130.68 J/K/mole
  !! O2: 205.15 J/K/mole
  !! C:    5.74 J/K/mole
  !! -> for H2O, to transform G(BermanBrown) to G(BensonHelgeson),
  !!____add 298.15*(130.68 + 205.15/2)=  69544.98 J
  !
  if(DtbConv_Benson) GH2O= GH2O + Tref *(S0_Hydrogen *Two + S0_Oxygen)
  Grt=  GH2O /T/R_JK
  !
  return
end subroutine WHAAR2

subroutine PSAT2(T,PS)
  real(dp):: T,PS,A(8),W,WSQ,V,FF
  integer :: I
  !
  data A &
  /-7.8889166D0,  2.5514255D0,  -6.716169D0, 33.239495D0, &
  & -105.38479D0, 174.35319D0, -148.39348D0, 48.631602D0/
  !
  if (T <= 314.00D0) then
    PS=Dexp(6.3573118D0-8858.843D0/T+607.56335D0/(T**0.6D0))
  else
    V=T/647.25D0
    W=abs(1.0D0-V)
    WSQ=sqrt(W)
    FF=0.0D0
    do 11,I=1,8
      FF=FF+A(I)*W
    11 W=W*WSQ
    PS=220.93D0*Dexp(FF/V)
  end if
  !
  return
end subroutine PSAT2

subroutine CUBIC(B,C,D,X1,X2,X2I,X3)
  real(dp):: B,C,D
  real(dp):: X1,X2,X2I,X3
  real(dp):: PI,Q,P,R,PHI3,FF
  !
  PI=3.14159263538979D0
  X2=0.0D0
  X2I=0.0D0
  X3=0.0D0
  !
  if (C == 0.0D0.and.D == 0.0D0) then
    X1=-B
    return
  end if
  !
  Q=((2.D0*B*B*B)/(27.D0)-(B*C)/(3.D0)+D)/2.D0
  P=(3.D0*C-B*B)/(9.D0)
  FF=abs(P)
  R=sqrt(FF)
  FF=R*Q
  if (FF < 0.0D0) R=-R
  FF=Q/(R*R*R)
  !
  if (P > 0.0D0) then
    PHI3=Dlog(FF+sqrt(FF*FF+1.D0))/3.D0
    X1=-R*(Dexp(PHI3)-Dexp(-PHI3))-B/(3.D0)
    X2I=1.D0
  else
    if (Q*Q+P*P*P > 0.0D0) then
      PHI3=Dlog(FF+sqrt(FF*FF-1.D0))/3.D0
      X1=-R*(Dexp(PHI3)+Dexp(-PHI3))-B/(3.D0)
      X2I=1.D0
    else
      PHI3=dataN(sqrt(1.D0-FF*FF)/FF)/3.D0
      X1=-2.D0*R*DCOS(PHI3)-B/(3.D0)
      X2=2.D0*R*DCOS(PI/3.D0-PHI3)-B/(3.D0)
      X2I=0.D0
      X3=2.D0*R*DCOS(PI/3.D0+PHI3)-B/(3.D0)
    end if
  end if
  !
  return
end subroutine CUBIC

end module M_Eos_H2O_Haar
