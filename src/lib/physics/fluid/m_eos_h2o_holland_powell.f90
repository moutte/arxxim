MODULE M_Eos_H2O_Holland_Powell
!! was M_Holland_Powell_H2O

!----------------------------------------------------------------------
! H2O HP91      : Eos_H2O_HolPow91
! H2O HOLPOW    : Eos_H2O_HolPow
!----------------------------------------------------------------------
! Purpose : Holland and Powel EOS for Pure H2O - Water and Steam
!
! References :
!        [HolPow] Holland and Powell (1990), from older CdeCapitani
!        [HP91]   Holland and Powell (1991,1998), from CdeCapitani  
!
! Source : CdeCapitani Theriak [ fsol.f90 ] 
!
!------------------------------------------------------------------------
! Arxim Integration : J.Moutte
!------------------------------------------------------------------------

  USE M_Kinds
  IMPLICIT NONE
  PRIVATE

  REAL(dp), PARAMETER :: R_jK= 8.314510D0

  !// Public Functions
  PUBLIC::  Eos_H2O_HolPow     ! was CalcGH2O_HolPow
  PUBLIC::  Eos_H2O_HolPow91   ! was CalcGH2O_HP91

  PRIVATE:: CalcGT_H2O_HP91

CONTAINS

!---

SUBROUTINE Eos_H2O_HolPow(TdgK,Pbar,Grt) 
!--------------------------------------------------------
! Holland and Powell (1990), from older CdeCapitani !
!--------------------------------------------------------
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(OUT):: Grt
  REAL(dp) &
  & A1,A2,A3,A4,A5,    &
  & B1,B2,B3,B4,B5,B6, &
  & C1,C2,C3,C4,C5,C6, &
  & T,T0,TT,SQT,Pkb,Pk2,Pk3,SQP, &
  & RTLNF,             &
  & K1,K2,K3,K4,       &
  & CPRDT,CPRTDT,      &
  & H0,S0,TT0,SQT0
  !
  H0=-241.81D0; S0= 188.80D-3
  K1= 0.0401D0; K2= 0.8656D-5; K3=487.5D0; K4=-0.2512D0
  T0= 298.15D0; TT0=88893.4225D0
  !
  DATA A1,A2,A3,A4,A5 &
  /-40.338D0,  1.6474D0,   -0.0062115D0, 2.0068D0,      0.0562929D0/
  DATA B1,B2,B3,B4,B5,B6 &
  /0.117372D0, Zero,      Zero,        -0.00046710D0, Zero,      Zero/
  DATA C1,C2,C3,C4,C5,C6 &
  /-7.3681D-6, 1.10295D-7, -9.8774D-7,   -2.4819D-5,    8.2948D-6,  8.33667D-8/
  !
  T=    TdgK
  SQT0= SQRT(T0)
  TT=   T*T
  SQT=  SQRT(T)
  !
  CPRDT= K1 *(T    -T0    )      &
  &    + K2 *(TT   -TT0   ) /Two &
  &    - K3 *(One/T-One/T0)      &
  &    + K4 *(SQT  -SQT0  ) *Two
  CPRTDT= K1*LOG(T/T0)                &
  &     + K2*(T       -T0       )     &
  &     - K3*(One/(TT)-One/(TT0))/Two &
  &     - K4*(One/SQT -One/SQT0 )*Two
  !
  Pkb=  Pbar/1.0D3
  Pk2=  Pkb*Pkb
  Pk3=  Pk2*Pkb
  SQP=  SQRT(Pkb)
  !
  RTLnF=   A1 +A2*Pkb +A3*Pk2 +A4/Pkb +A5/Pk2 &
  &    + ( B1 +B2*Pkb         +B3/Pkb +B4/Pk2 +B5/SQP +B6/Pk3 ) *T &
  &    + ( C1 +C2*Pkb         +C5/Pkb +C3/Pk2 +C4/SQP +C6/Pk3 ) *T*T
  !
  Grt=   (H0 +CPRDT -T*(S0+CPRTDT) +RTLnF) *1.0D3 /R_jk/TdgK
  
  RETURN
ENDSUBROUTINE Eos_H2O_HolPow

!---

SUBROUTINE Eos_H2O_HolPow91(TdgK,Pbar,Grt,V_H2O_m3)
  !-------------------------------------------------------
  !.Holland and Powell (1991,1998), from CdeCapitani ! 
  !-------------------------------------------------------
  USE M_Eos_Utils
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(OUT):: Grt,V_H2O_m3
  !
  REAL(dp):: &
    T,SQT,RT,DT,&
    A,AGas,G,GT,&
    delP,Pkb,PSA,VMRK,&
    cubB,cubC,cubD,X1,X2,X2I,X3,&
    V1,V2,VREF,PREF,AF, & !RTLNP,&
    V_H2Ocm3, &
    A1,A2,A3,A4,A5,A6,A7,A8,A9,&
    dGGas1,dGGas2,dGGas3,dGGas
    !GA1,GA2,GA3,GF
  REAL(dp):: &
  & AVir=1.9853D-3, BVir=-8.9090D-2, CVir=8.0331D-2, &
  & R=   8.3142D-3, &
  & P0=  Two,       TK=  695.0D0,    TA=  673.0D0, &
  & A0=  1113.4D0,  B=   1.465D0 
  DATA A1,A2,A3,A4,A5,A6,A7,A8,A9 &
   /-0.88517D0,  4.5300D-3, -1.3183D-5, &
    -0.22291D0, -3.8022D-4,  1.7791D-7, &
     5.84870D0, -2.1370D-2,  6.8133D-5/
  !
  T=     TdgK
  RT=    R*T
  Pkb=   Pbar/1000.0D0 !-> Pkb in kbar
  SQT=   SQRT(T)
  
  G=     Zero
  delP=  Pkb-P0
  dGGas= Zero; dGGas1=Zero; dGGas2=Zero; dGGas3=Zero
  V1=    Zero; V2=    Zero
  
  Pref=  0.001D0 !=0.001 kbar
  
  PSA= - 13.6270D-3        &
       + 7.29395D-7  *T**2 &
       - 2.34622D-9  *T**3 &
       + 4.83607D-15 *T**5
  
  IF (T<TA) THEN
    DT=  TA-T
    !A=   A0 +A1*(TA-T) +A2*(TA-T)**2 +A3*(TA-T)**3
    !AGas=A0 +A7*(TA-T) +A8*(TA-T)**2 +A9*(TA-T)**3
    A=   A0 + DT*( A1 + DT*(A2 + DT*A3) )
    AGas=A0 + DT*( A7 + DT*(A8 + DT*A9) )
  ELSE
    DT=  T-TA
    A=   A0 + DT*( A4 + DT*(A5 + DT*A6) )
    AGas=Zero
  ENDIF
  
  VREF=RT/PREF
  
  !volume at Pkb and T
  IF (T<TA.AND.Pkb<PSA) THEN; AF=AGas
  ELSE;                       AF=A
  ENDIF
  !
  cubB=-RT/Pkb
  cubC=-(B*RT+B*B*Pkb-AF/SQT)/Pkb
  cubD=-AF*B/SQT/Pkb
  CALL CUBIC(cubB,cubC,cubD,X1,X2,X2I,X3)
  
  IF (X2I/=Zero) THEN
    VMRK=X1
  ELSE
    IF (Pkb<PSA.AND.T<TK)  THEN
      VMRK=DMAX1(X1,X2,X3)
    ELSE
      VMRK=DMIN1(X1,X2,X3)
      IF (VMRK<B) VMRK=DMAX1(X1,X2,X3)
    ENDIF
  ENDIF
  
  V_H2Ocm3=VMRK
  
  IF (Pkb>P0) &
  & V_H2Ocm3=  VMRK       &
  &    + AVir *delP       &
  &    + BVir *SQRT(delP) &
  &    + CVir *delP**0.25D0
  
  !CALL GAGA(AF,B,Pkb,VMRK,T,GA3)
  
  dGGas3=  &
  &   VMRK *Pkb &
  & - VREF *PREF &
  & - RT           *LOG((VMRK-B)      /(VREF-B)     ) &
  & + (AF/(SQT*B)) *LOG( VMRK*(VREF+B)/VREF/(VMRK+B))
  
  IF (Pkb>P0) &
  & dGGas3= dGGas3 &
  &  + AVir/Two       *delP*delP   &
  &  + BVir*Two/3.0D0 *delP**1.5D0 &
  &  + CVir*0.8D0     *delP**1.25D0
  
  !volume of gas at T and PSAT
  IF (T<TK.AND.Pkb>PSA) THEN
  
    IF (T>TA) THEN; AF=A
    ELSE;           AF=AGas
    ENDIF
    cubB= -   RT                 /PSA
    cubC= -(B*RT+B*B*PSA-AF/SQT) /PSA
    cubD= -              AF*B/SQT/PSA
    !
    CALL CUBIC(cubB,cubC,cubD,X1,X2,X2I,X3)
    !
    IF (X2I/=Zero) THEN; V1=X1
    ELSE;                V1=MAX(X1,X2,X3)
    ENDIF
    
    !CALL GAGA(AF,B,PSA,V1,T,GA1)
    dGGas1=  V1*PSA &
    &     -  VREF*PREF &
    &     -  RT *LOG((V1-B)/(VREF-B)) &
    &     + (AF/(SQT*B)) *LOG(V1 *(VREF+B)/(VREF*(V1+B)))
    
    !volume of liquid at a T and PSAT
    cubB= -   RT /PSA
    cubC= -(B*RT+B*B*PSA-A/SQT) /PSA
    cubD= -            A*B/SQT  /PSA
    !
    CALL CUBIC(cubB,cubC,cubD,X1,X2,X2I,X3)
    !
    IF (X2I/=Zero) THEN; V2=X1
    ELSE;                V2=DMIN1(X1,X2,X3)
    ENDIF
    
    !CALL GAGA(A,B,PSA,V2,T,GA2)
    dGGas2=  VMRK*Pkb &
    &      - V2*PSA &
    &      - RT*LOG((VMRK-B) /(V2-B)) &
    &      + (A/(SQT*B)) *LOG( VMRK *(V2+B) /(V2*(VMRK+B)) )
    IF (Pkb>P0) &
    & dGGas2=  dGGas2 &
    &       + AVir /Two       *delP**2     &
    &       + BVir *Two/3.0D0 *delP**1.5D0 &
    &       + CVir *0.8D0     *delP**1.25D0
    !GF=RT*(GA1-GA2+GA3)
    dGGas=1.0D3*(dGGas2+dGGas1)
  
  ELSE
    
    !GF=RT*GA3
    dGGas=1.0D3*dGGas3
    
  ENDIF
  !
  CALL CalcGT_H2O_HP91(T,GT)
  !
  G=  GT +dGGas
  !
  V_H2O_m3= V_H2Ocm3 /1.0D6
  Grt=      G /RT
  
  RETURN
END SUBROUTINE Eos_H2O_HolPow91

!---

SUBROUTINE CalcGT_H2O_HP91(Tk,GT)
!--------------------------------------------------------
!.compute G(H2O) at (T,Pref)
!.Holland and Powell (1991,1998), from CdeCapitani codes
!--------------------------------------------------------
  IMPLICIT NONE
  REAL(dp),INTENT(IN) :: Tk
  REAL(dp),INTENT(OUT):: GT
  !
  REAL(dp):: K1,K2,K3,K4
  REAL(dp):: H0,S0,T0,TT0
  REAL(dp):: TT,SQT,CPRDT,CPRTDT,SQT0
  DATA K1,K2,K3,K4 /0.0401D0, 0.8656D-5,  487.5D0, -0.2512D0/
  DATA H0,S0       /-241.81D0,188.80D-3/
  DATA T0,TT0      /298.15D0, 88893.4225D0/
  !
  SQT0=  SQRT(T0)
  TT=    Tk*Tk
  SQT=   SQRT(Tk)
  !
  CPRDT= K1*(Tk     -T0)         &
  &    + K2*(TT     -TT0)   /Two &
  &    - K3*(One/Tk -One/T0)     &
  &    + K4*(SQT    -SQT0)  *Two
  CPRTDT=K1 *LOG(Tk   /T0)      &
  &    + K2 *(Tk      -T0)      &
  &    - K3 *(One/TT  -One/TT0) /Two &
  &    - K4 *(One/SQT -One/SQT0)*Two
  !
  GT= ( H0 +CPRDT -Tk*(S0 +CPRTDT) ) *1.0D3
  RETURN
ENDSUBROUTINE CalcGT_H2O_HP91

!---

END MODULE M_Eos_H2O_Holland_Powell
