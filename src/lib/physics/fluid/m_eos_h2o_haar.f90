MODULE M_Eos_H2O_Haar
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

USE M_Kinds
IMPLICIT NONE
PRIVATE

!REAL(dp), PARAMETER :: R_jK= 8.314510D0

!// PUBLIC FUNCTIONS
PUBLIC:: Eos_H2O_Haar !! CalcGH2O_Haar

CONTAINS

SUBROUTINE Eos_H2O_Haar(TdgK,Pbar,G_H2Ort,V_H2O_m3)
!--
!-- CdeCapitani, SupCrt
!-- compute Gibbs free energy of Water G_H2O, using EoS (Haar-Gallager-Kell, '84)
!-- same as used in Eos_H2O_Rho
!--
  USE M_Dtb_Const,ONLY: R_jK, Tref,DtbConv_Benson,S0_Hydrogen,S0_Oxygen
  USE M_Eos_H2O_Rho_Psat
  USE M_Eos_Utils
  !
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(OUT):: G_H2Ort,V_H2O_m3
  !
  REAL(dp):: &
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
  INTEGER:: LOO, I
  !
  REAL(dp)::RHOI(37:40)=(/0.319D0, 0.310D0, 0.310D0,  1.550D0  /)
  REAL(dp)::TTTI(37:40)=(/640.0D0, 640.0D0, 641.6D0,  270.0D0  /)
  REAL(dp)::ALPI(37:40)=(/34.0D0,  40.0D0,  30.0D0,   1050.0D0 /)
  REAL(dp)::BETI(37:40)=(/2.0D4,   2.0D4,   4.0D4,    25.0D0   /)
  !
  INTEGER::KI(1:40)=(/&
  &  1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5, &
  &  6,6,6,6, 7,7,7,7, 9,9,9,9, 3,3,1,5, 2,2,2,4  /)
  INTEGER::LI(1:40)=(/&
  &  1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, &
  &  1,2,4,6, 1,2,4,6, 1,2,4,6, 0,3,3,3, 0,2,0,0  /)
  REAL(dp)::CI(1:18)=(/&
  &  .19730271018D2,     .209662681977D2,  -.483429455355D0, &
  &  .605743189245D1,  22.56023885D0,      -9.87532442D0,    -.43135538513D1,&
  &  .458155781D0,      -.47754901883D-1,   .41238460633D-2, -.27929052852D-3,&
  &  .14481695261D-4,   -.56473658748D-6,   .16200446D-7,    -.3303822796D-9,&
  &  .451916067368D-11, -.370734122708D-13, .137546068238D-15 /)
  !
  !GI ARE IN (bar cc / g)=    10 * (J / g)
  REAL(dp)::GI(1:40)=(/&
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
  DO I=2,6
    TAUI(I)=TAUI(I-1)*TAUI(1)
  ENDDO
  !
  B= - 0.3540782D0 *LOG(TAUI(1)) &
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
  IF (T<=647.25D0) CALL Eos_H2O_psat(T,PS)    !IF (P >= PS) NAME='STEAM'
  !
  !-- SET INITIAL GUESS FOR RHO USING THB-FIT TO REDLICH-KWONG
  ARK=  1.279186D8-2.241415D4 *T   !RedlichKwong Param. A
  BRK=  1.428062D1+6.092237D-4*T   !RedlichKwong Param. B
  OFT=  ARK/(P*SQRT(T))
  cubB= -10.0D0*RR*T/P
  cubC= OFT -BRK*BRK +BRK*cubB
  cubD= -BRK*OFT
  !
  CALL CUBIC(cubB,cubC,cubD,X1,X2,X2I,X3)
  !
  IF(X2I/=Zero) THEN
    VOL=X1
  ELSE
    IF (P<PS)  THEN ; VOL=DMAX1(X1,X2,X3)
    ELSE            ; VOL=DMIN1(X1,X2,X3)
    ENDIF
  ENDIF
  IF (VOL<=Zero) THEN  ; RHN=1.9D0
  ELSE                 ; RHN=One/VOL*18.0152D0
  ENDIF
  !-- FIND THE TRUE(?) RH(T,P)
  !-- NOTE: PR-  PRESSURE CORRESPONDING TO GUESSED RH
  !-- DPR-  (dP_ / dRH)
  !-- the values (1-EXP(-RH))**i are stored in the array ERMI(i)
  DO LOO=1,100
    RH=RHN
    IF(RH<=1D-8)  RH=1D-8
    IF(RH> 1.9D0) RH=1.9D0
    RH2= RH*RH
    Y=   RH*B/4.D0
    ER=  EXP(-RH)
    Y3=  (One-Y)**3
    ALY= 11.D0*Y
    BETY=133.D0/3.D0*Y*Y
    F1=  (One+ALY+BETY)/Y3
    F2=  4.D0*Y*(BB/B-3.5D0)
    !
    ERMI(0)=One
    ERMI(1)=One-ER
    DO I=2,9; ERMI(I)=ERMI(I-1)*ERMI(1); ENDDO
    !
    PR=Zero
    DPR=Zero
    DO I=1,36
      S=    GI(I)/TAUI(LI(I))*ERMI(KI(I)-1)
      PR=   PR+S
      DPR=  DPR+(2D0+RH*(KI(I)*ER-One)/ERMI(1))*S
    ENDDO
    !
    DO I=37,40
      DEL=   RH/RHOI(I)-One
      RHOI2= RHOI(I)*RHOI(I)
      TAU=   T/TTTI(I)-One
      QHEX= (-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
      IF (QHEX>-150.0D0) THEN
        Q10= GI(I) *DEL**LI(I) *EXP(-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
      ELSE
        Q10= Zero
      ENDIF
      QM=  LI(I)/DEL-KI(I)*ALPI(I)*DEL**(KI(I)-1)
      S=   Q10*QM*RH2/RHOI(I)
      PR=  PR+S
      DPR= DPR &
         + S*(Two/RH +QM/RHOI(I)) &
         - RH2/RHOI2*Q10 &
           *(LI(I)/DEL/DEL+KI(I)*(KI(I)-1)*ALPI(I)*DEL**(KI(I)-2))
    ENDDO
    PR= RH*(RH*ER*PR+RT*(F1+F2))
    DPR=RH*ER*DPR+RT*((One+Two*ALY+3.D0*BETY)/Y3 +3.D0*Y*F1/(One-Y)+2.D0*F2)
    IF (DPR <= Zero) THEN
      IF (P<=PS) THEN; RHN=RHN*0.95D0
      ELSE;            RHN=RHN*1.05D0
      ENDIF
    ELSE
      IF (DPR < 0.01D0) DPR=0.01D0
      S=(P-PR)/DPR
      IF (ABS(S) > 0.1D0) S=0.1D0*S/ABS(S)
      RHN=RH+S
    ENDIF
    DP_=ABS(One-PR/P)
    DR_=ABS(One-RHN/RH)
    IF (DP_ < 1.D-5.AND.DR_ < 1.D-5) EXIT
  ENDDO
  !
  RH=     RHN
  Y=      RH*B/4.D0
  X=      One-Y
  ER=     EXP(-RH)
  ERMI(0)=One
  ERMI(1)=One-ER
  DO I=2,9; ERMI(I)=ERMI(I-1)*ERMI(1); ENDDO
  !
  !--- base function --
  AA= RT*(-LOG(X) &
    - 130.0D0/3.0D0 /X     & !
    + 169.0D0/6.0D0 /X/X   & !
    + 4.D0*Y*(BB/B-3.5D0)  &
    + 91.D0/6.D0           & !
    + LOG(RH*RT/1.01325D0))
  !
  !--- residual function --
  DO I=1,36
    AA= AA &
    & + GI(I) /KI(I) /TAUI(LI(I)) *ERMI(KI(I))
  ENDDO
  DO I=37,40
    DEL=  RH/RHOI(I)-One
    TAU=  T/TTTI(I)-One
    QHEX= (-ALPI(I)*DEL**KI(I)-BETI(I)*TAU*TAU)
    IF (QHEX > -150.0D0) AA=AA+GI(I)*DEL**LI(I)*EXP(QHEX)
  ENDDO
  !
  !--- ideal gas function --
  TR=  T/1.0D2
  W=   TR**(-3)
  AID= One +( CI(1)/TR +CI(2) )*LOG(TR)
  DO I=3,18
    AID= AID +CI(I)*W
    W=   W*TR
  ENDDO
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
  !Robie/CODATA DATA for entropies of elements at 298.15:
  !! H2: 130.68 J/K/mole
  !! O2: 205.15 J/K/mole
  !! C:    5.74 J/K/mole
  !! -> for H2O, to transform G(BermBrown) to G(BensHelg),
  !!____add 298.15*(130.68 + 205.15/2)=  69544.98 J
  !
  IF(DtbConv_Benson) G_H2O=   G_H2O + Tref *(S0_Hydrogen *Two + S0_Oxygen)
  G_H2Ort= G_H2O /R_jK /TdgK
  !
  RETURN
ENDSUBROUTINE Eos_H2O_Haar

SUBROUTINE WHAAR2(P,T,Grt,VH2O)
!-- for tests, not used !!
!--
  USE M_Dtb_Const,ONLY: R_JK,Tref,DtbConv_Benson,S0_Hydrogen,S0_Oxygen
  !
  REAL(dp),INTENT(IN)    :: P,T
  REAL(dp),INTENT(OUT)   :: Grt,VH2O
  !
  CHARACTER(LEN=16):: NAME
  !
  REAL(dp):: TAUI(0:6),ERMI(0:9),GI(40),CI(18)
  INTEGER :: KI(40),LI(40)
  REAL(dp):: RHOI(37:40),TTTI(37:40),ALPI(37:40),BETI(37:40)
  REAL(dp):: GH2O,R,RR,RT,GREF,T0,PS
  REAL(dp):: B,BB,ARK,BRK,X,Y,Y3,VOL,RHN,RH,RH2,F1,F2,PR,DPR,ER,S,AA
  REAL(dp):: RHOI2,TAU,DEL,QHEX,Q10,QM,DP_,DR,TR,W,AID
  REAL(dp):: OFT,BUK,CUK,DUK,X1,X2,X2I,X3
  REAL(dp):: ALY,BETY
  INTEGER:: I,LOO
  !
  !-- GI ARE IN (bar cc / g)  -  10 * (J / g)
  DATA GI/-.53062968529023D4,.22744901424408D5,.78779333020687D4,&
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
  DATA KI/4*1,4*2,4*3,4*4,4*5,4*6,4*7,4*9,2*3,1,5,3*2,4/
  DATA LI/1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,1,2,4,6,0,3*3,0,2,0,0/
  DATA CI/.19730271018D2,    .209662681977D2,    -.483429455355D0,&
          .605743189245D1,   22.56023885D0,      -9.87532442D0,    -.43135538513D1,&
          .458155781D0,-.47754901883D-1,.41238460633D-2,-.27929052852D-3,&
          .14481695261D-4,-.56473658748D-6,.16200446D-7,-.3303822796D-9,&
          .451916067368D-11,-.370734122708D-13,.137546068238D-15/
  DATA RHOI/0.319D0,0.310D0,0.310D0,1.550D0/
  DATA TTTI/640.0D0,640.0D0,641.6D0,270.0D0/
  DATA ALPI/34.0D0,40.0D0,30.0D0,1050.0D0/
  DATA BETI/2.0D4,2.0D4,4.0D4,25.0D0/
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
  DO I=2,6
    TAUI(I)=TAUI(I-1)*TAUI(1)
  END DO
  !
  B  = -0.3540782D0*DLOG(TAUI(1)) &
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
  IF (T <= 647.25D0) THEN
     CALL PSAT2(T,PS)
     IF (P >= PS) NAME='STEAM'
  END IF

  !-- SET INITIAL GUESS FOR RHO USING THB-FIT TO REDLICH-KWONG
  ARK= 1.279186D8 -2.241415D4*T
  BRK= 1.428062D1 +6.092237D-4*T
  RR=  8.31441D0
  OFT= ARK/(P*DSQRT(T))
  BUK=-10D0*RR*T/P
  CUK=OFT-BRK*BRK+BRK*BUK
  DUK=-BRK*OFT

  CALL CUBIC(BUK,CUK,DUK,X1,X2,X2I,X3)

  IF (X2I /= 0.0D0) THEN
    VOL=X1
  ELSE
    IF (P < PS)  THEN
      VOL=DMAX1(X1,X2,X3)
    ELSE
      VOL=DMIN1(X1,X2,X3)
    END IF
  END IF

  IF (VOL <= 0.0D0) THEN
    RHN=1.9D0
  ELSE
    RHN=1D0/VOL*18.0152D0
  END IF

  !-- FIND THE TRUE(?) RH(T,P)
  !-- NOTE: PR - PRESSURE CORRESPONDING TO GUESSED RH
  !-- DPR - (DP_ / dRH)
  !-- the values (1-EXP(-RH))**i are stored in the array ERMI(i)
  DO LOO=1,100

    RH=RHN
    IF (RH <= 1D-8) RH=1D-8
    IF (RH > 1.9D0) RH=1.9D0

    RH2=RH*RH
    Y=RH*B/4.D0
    ER=DEXP(-RH)
    Y3=(1D0-Y)**3
    ALY=11.D0*Y
    BETY=44.33333333333333D0*Y*Y

    F1=(1.D0+ALY+BETY)/Y3
    F2=4.D0*Y*(BB/B-3.5D0)
    !
    ERMI(0)=1D0
    ERMI(1)=1D0-ER
    DO I=2,9
      ERMI(I)=ERMI(I-1)*ERMI(1)
    END DO
    !
    PR=  0.0D0
    DPR= 0.0D0
    !
    DO I=1,36
      S=GI(I)/TAUI(LI(I))*ERMI(KI(I)-1)
      PR=PR+S
      DPR=DPR+(2D0+RH*(KI(I)*ER-1D0)/ERMI(1))*S
    END DO
    !
    DO I=37,40
      DEL=RH/RHOI(I)-1.0D0
      RHOI2=RHOI(I)*RHOI(I)
      TAU=T/TTTI(I)-1.0D0
      QHEX=(-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
      IF (QHEX > -150.0D0) THEN
        Q10=GI(I)*DEL**LI(I)*DEXP(-ALPI(I)*DEL**LI(I)-BETI(I)*TAU*TAU)
      ELSE
        Q10=0.0D0
      END IF
      QM=  LI(I)/DEL-KI(I)*ALPI(I)*DEL**(KI(I)-1)
      S=   Q10*QM*RH2/RHOI(I)
      PR=  PR+S
      DPR= DPR &
      &  + S*( 2.0D0/RH + QM/RHOI(I) ) &
      &  - RH2/RHOI2*Q10 &
      &    *( LI(I)/DEL/DEL +KI(I)*(KI(I)-1)*ALPI(I)*DEL**(KI(I)-2) )
    END DO
    !
    PR=  RH*(RH*ER*PR+RT*(F1+F2))
    DPR= RH*ER*DPR &
    &  + RT*( (1.D0  + 2.D0*ALY + 3.D0*BETY)/Y3 &
    &       + 3.D0*Y*F1/(1.D0-Y) &
    &       + 2.D0*F2)
    !
    IF (DPR <= 0.0D0) THEN
      IF (P <= PS) THEN
        RHN=RHN*0.95D0
      ELSE
        RHN=RHN*1.05D0
      END IF
    ELSE
      IF (DPR < 0.01D0) DPR=0.01D0
      S=(P-PR)/DPR
      IF (DABS(S) > 0.1D0) S=0.1D0*S/DABS(S)
      RHN=RH+S
    END IF
    DP_= DABS(1.0D0-PR/P)
    DR=  DABS(1.0D0-RHN/RH)

    IF (DP_ < 1D-5.AND.DR < 1D-5) EXIT !GOTO 30

  END DO !20

  RH=RHN
  Y=RH*B/4.D0
  X=1D0-Y
  ER=DEXP(-RH)
  ERMI(0)=1D0
  ERMI(1)=1D0-ER
  !
  DO I=2,9
    ERMI(I)=ERMI(I-1)*ERMI(1)
  END DO
  !
  !--- CALCULATE BASE FUNCTION
  AA= RT*(-DLOG(X) &
     -43.33333333333333D0/X &
     +28.16666666666667D0/X/X &
     +4D0*Y*(BB/B-3.5D0) &
     +15.16666666666667D0 &
     +DLOG(RH*RT/1.01325D0))
  !
  !--- CALCULATE RESIDUAL FUNCTION
  DO I=1,36
    AA=AA+GI(I)/KI(I)/TAUI(LI(I))*ERMI(KI(I))
  END DO
  DO I=37,40
    DEL=RH/RHOI(I)-1.0D0
    TAU=T/TTTI(I)-1.0D0
    QHEX=(-ALPI(I)*DEL**KI(I)-BETI(I)*TAU*TAU)
    IF (QHEX > -150.0D0) THEN
      AA=AA+GI(I)*DEL**LI(I)*DEXP(QHEX)
    END IF
  END DO
  !---/

  !--- CALCULATE IDEAL GAS FUNCTION (-AID, HELMHOLTZ FUNCTION, IDEAL)
  TR=  T/1.0D2
  W=   TR**(-3)
  AID= 1.D0 + (CI(1)/TR+CI(2))*DLOG(TR)
  DO I=3,18
    AID= AID +CI(I)*W
    W=   W*TR
  END DO
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
  !Robie/CODATA DATA for entropies of elements at 298.15:
  !! H2: 130.68 J/K/mole
  !! O2: 205.15 J/K/mole
  !! C:    5.74 J/K/mole
  !! -> for H2O, to transform G(BermBrown) to G(BensHelg),
  !!____add 298.15*(130.68 + 205.15/2)=  69544.98 J
  !
  IF(DtbConv_Benson) GH2O= GH2O + Tref *(S0_Hydrogen *Two + S0_Oxygen)
  Grt=  GH2O /T/R_JK
  !
  RETURN
END SUBROUTINE WHAAR2

SUBROUTINE PSAT2(T,PS)
  REAL(dp):: T,PS,A(8),W,WSQ,V,FF
  INTEGER :: I
  !
  DATA A &
  /-7.8889166D0,  2.5514255D0,  -6.716169D0, 33.239495D0, &
  & -105.38479D0, 174.35319D0, -148.39348D0, 48.631602D0/
  !
  IF (T <= 314.00D0) THEN
    PS=DEXP(6.3573118D0-8858.843D0/T+607.56335D0/(T**0.6D0))
  ELSE
    V=T/647.25D0
    W=DABS(1.0D0-V)
    WSQ=DSQRT(W)
    FF=0.0D0
    DO 11,I=1,8
      FF=FF+A(I)*W
    11 W=W*WSQ
    PS=220.93D0*DEXP(FF/V)
  END IF
  !
  RETURN
END SUBROUTINE PSAT2

SUBROUTINE CUBIC(B,C,D,X1,X2,X2I,X3)
  REAL(dp):: B,C,D
  REAL(dp):: X1,X2,X2I,X3
  REAL(dp):: PI,Q,P,R,PHI3,FF
  !
  PI=3.14159263538979D0
  X2=0.0D0
  X2I=0.0D0
  X3=0.0D0
  !
  IF (C == 0.0D0.AND.D == 0.0D0) THEN
    X1=-B
    RETURN
  END IF
  !
  Q=((2.D0*B*B*B)/(27.D0)-(B*C)/(3.D0)+D)/2.D0
  P=(3.D0*C-B*B)/(9.D0)
  FF=DABS(P)
  R=DSQRT(FF)
  FF=R*Q
  IF (FF < 0.0D0) R=-R
  FF=Q/(R*R*R)
  !
  IF (P > 0.0D0) THEN
    PHI3=DLOG(FF+DSQRT(FF*FF+1.D0))/3.D0
    X1=-R*(DEXP(PHI3)-DEXP(-PHI3))-B/(3.D0)
    X2I=1.D0
  ELSE
    IF (Q*Q+P*P*P > 0.0D0) THEN
      PHI3=DLOG(FF+DSQRT(FF*FF-1.D0))/3.D0
      X1=-R*(DEXP(PHI3)+DEXP(-PHI3))-B/(3.D0)
      X2I=1.D0
    ELSE
      PHI3=DATAN(DSQRT(1.D0-FF*FF)/FF)/3.D0
      X1=-2.D0*R*DCOS(PHI3)-B/(3.D0)
      X2=2.D0*R*DCOS(PI/3.D0-PHI3)-B/(3.D0)
      X2I=0.D0
      X3=2.D0*R*DCOS(PI/3.D0+PHI3)-B/(3.D0)
    END IF
  END IF
  !
  RETURN
END SUBROUTINE CUBIC

ENDMODULE M_Eos_H2O_Haar
