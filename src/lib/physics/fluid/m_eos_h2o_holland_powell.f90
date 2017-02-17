module M_Eos_H2O_Holland_Powell
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

  use M_Kinds
  implicit none
  private

  real(dp), parameter :: R_jK= 8.314510D0

  !// Public Functions
  public::  Eos_H2O_HolPow     ! was CalcGH2O_HolPow
  public::  Eos_H2O_HolPow91   ! was CalcGH2O_HP91

  private:: CalcGT_H2O_HP91

contains

!---

subroutine Eos_H2O_HolPow(TdgK,Pbar,Grt,LnFug) 
!--------------------------------------------------------
! Holland and Powell (1990), from older CdeCapitani !
!--------------------------------------------------------
  implicit none
  !
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(out):: Grt,LnFug
  !
  real(dp) &
  & A1,A2,A3,A4,A5,    &
  & B1,B2,B3,B4,B5,B6, &
  & C1,C2,C3,C4,C5,C6, &
  & T,T0,TT,SQT,Pkb,Pk2,Pk3,SQP, &
  & RTLNF,             &
  & K1,K2,K3,K4,       &
  & CPRDT,CPRTDT,      &
  & H0,S0,TT0,SQT0,    &
  & S0Ele
  !
  H0=-241.81D0  ; S0= 188.80D-3
  K1= 0.0401D0  ; K2= 0.8656D-5    ; K3=487.5D0  ; K4=-0.2512D0
  T0= 298.15D0  ; TT0=88893.4225D0
  !
  data A1,A2,A3,A4,A5 &
  /-40.338D0,  1.6474D0,   -0.0062115D0, 2.0068D0,      0.0562929D0/
  data B1,B2,B3,B4,B5,B6 &
  /0.117372D0, Zero,      Zero,        -0.00046710D0, Zero,      Zero/
  data C1,C2,C3,C4,C5,C6 &
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
  CPRTDT= K1*log(T/T0)                &
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
  S0Ele= 102.575D0 + 65.34D0 *2.0D0
  Grt= (H0 +CPRDT -T*(S0+CPRTDT) +RTLnF)*1.0D3 
  Grt= Grt + S0Ele*298.15D0 !--------------------------Benson Convention
  Grt= Grt/R_jk/TdgK
  LnFug= RTLnF*1.0D3 /R_jk/TdgK
  !
  return
end subroutine Eos_H2O_HolPow

!---

subroutine Eos_H2O_HolPow91(TdgK,Pbar,Grt,V_m3,LnFug)
  !-------------------------------------------------------
  !.Holland and Powell (1991,1998), from CdeCapitani ! 
  !-------------------------------------------------------
  use M_Eos_Utils
  implicit none
  !
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(out):: Grt,V_m3,LnFug
  !
  real(dp):: &
  & T,SQT,RT,DT,&
  & A,AGas,G,GT,&
  & delP,Pkb,PSA,VMRK,&
  & cubB,cubC,cubD,X1,X2,X2I,X3,&
  & V1,V2,VREF,PREF,AF, & !RTLNP,&
  & Volum,S0Ele, &
  & A1,A2,A3,A4,A5,A6,A7,A8,A9,&
  & dGGas1,dGGas2,dGGas3,dGGas
  ! GA1,GA2,GA3,GF
  real(dp):: &
  & AVir=1.9853D-3, BVir=-8.9090D-2, CVir=8.0331D-2, &
  & R=   8.3142D-3, &
  & P0=  Two,       TK=  695.0D0,    TA=  673.0D0, &
  & A0=  1113.4D0,  B=   1.465D0 
  data A1,A2,A3,A4,A5,A6,A7,A8,A9 &
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
  dGGas= Zero
  dGGas1=Zero; dGGas2=Zero; dGGas3=Zero
  V1=    Zero; V2=    Zero
  
  Pref=  0.001D0 !=0.001 kbar
  
  PSA= - 13.6270D-3        &
  &    + 7.29395D-7  *T**2 &
  &    - 2.34622D-9  *T**3 &
  &    + 4.83607D-15 *T**5
  
  if (T<TA) then
    DT=  TA-T
    !A=   A0 +A1*(TA-T) +A2*(TA-T)**2 +A3*(TA-T)**3
    !AGas=A0 +A7*(TA-T) +A8*(TA-T)**2 +A9*(TA-T)**3
    A=   A0 + DT*( A1 + DT*(A2 + DT*A3) )
    AGas=A0 + DT*( A7 + DT*(A8 + DT*A9) )
  else
    DT=  T-TA
    A=   A0 + DT*( A4 + DT*(A5 + DT*A6) )
    AGas=Zero
  end if
  
  VREF=RT/PREF
  
  !volume at Pkb and T
  if (T<TA.and.Pkb<PSA) then; AF=AGas
  else;                       AF=A
  end if
  !
  cubB=-RT/Pkb
  cubC=-(B*RT+B*B*Pkb-AF/SQT)/Pkb
  cubD=-AF*B/SQT/Pkb
  call CUBIC(cubB,cubC,cubD,X1,X2,X2I,X3)
  
  if (X2I/=Zero) then
    VMRK=X1
  else
    if (Pkb<PSA.and.T<TK)  then
      VMRK=max(X1,X2,X3)
    else
      VMRK=DMIN1(X1,X2,X3)
      if (VMRK<B) VMRK=max(X1,X2,X3)
    end if
  end if
  
  Volum=VMRK
  
  if (Pkb>P0) &
  & Volum=  VMRK          &
  &    + AVir *delP       &
  &    + BVir *sqrt(delP) &
  &    + CVir *delP**0.25D0
  
  !call GAGA(AF,B,Pkb,VMRK,T,GA3)
  
  dGGas3= VMRK*Pkb - VREF*PREF &
  &     - RT           *log((VMRK-B)      /(VREF-B)     ) &
  &     + (AF/(SQT*B)) *log( VMRK*(VREF+B)/VREF/(VMRK+B))
  
  if (Pkb>P0) &
  & dGGas3= dGGas3 &
  &  + AVir/Two       *delP**2     &
  &  + BVir*Two/3.0D0 *delP**1.5D0 &
  &  + CVir*0.8D0     *delP**1.25D0
  
  !volume of gas at T and PSAT
  if (T<TK.and.Pkb>PSA) then
  
    if (T>TA) then  ; AF=A
    else            ; AF=AGas
    end if
    cubB= -   RT                 /PSA
    cubC= -(B*RT+B*B*PSA-AF/SQT) /PSA
    cubD= -              AF*B/SQT/PSA
    !
    call CUBIC(cubB,cubC,cubD,X1,X2,X2I,X3)
    !
    if (X2I/=Zero) then  ; V1=X1
    else                 ; V1=MAX(X1,X2,X3)
    end if
    
    !call GAGA(AF,B,PSA,V1,T,GA1)
    dGGas1=  V1*PSA &
    &     -  VREF*PREF &
    &     -  RT *log((V1-B)/(VREF-B)) &
    &     + (AF/(SQT*B)) *log(V1 *(VREF+B)/(VREF*(V1+B)))
    
    !volume of liquid at a T and PSAT
    cubB= -   RT /PSA
    cubC= -(B*RT+B*B*PSA-A/SQT) /PSA
    cubD= -            A*B/SQT  /PSA
    !
    call CUBIC(cubB,cubC,cubD,X1,X2,X2I,X3)
    !
    if (X2I/=Zero) then; V2=X1
    else;                V2=DMIN1(X1,X2,X3)
    end if
    
    !call GAGA(A,B,PSA,V2,T,GA2)
    dGGas2=  VMRK*Pkb &
    &      - V2*PSA &
    &      - RT*log((VMRK-B) /(V2-B)) &
    &      + (A/(SQT*B)) *log( VMRK *(V2+B) /(V2*(VMRK+B)) )
    if (Pkb>P0) &
    & dGGas2=  dGGas2 &
    &       + AVir /Two       *delP**2     &
    &       + BVir *Two/3.0D0 *delP**1.5D0 &
    &       + CVir *0.8D0     *delP**1.25D0
    !GF=RT*(GA1-GA2+GA3)
    dGGas= dGGas2+dGGas1 !kiloJoule
  
  else
    
    !GF=RT*GA3
    dGGas= dGGas3 !kiloJoule
    
  end if
  !
  call CalcGT_H2O_HP91(T,GT) ! GT in Joules
  !
  G=  GT +dGGas*1.0D3
  !
  V_m3= Volum *1.0D-5
  !
  S0Ele= 102.575D0 + 65.34D0 *2.0D0
  G=     G + S0Ele*298.15D0 !--------------------------Benson Convention
  Grt=   G /T /R_jK
  LnFug= dGGas*1.0D3 /T /R_jK
  !
  return
end subroutine Eos_H2O_HolPow91

!---

subroutine CalcGT_H2O_HP91(Tk,GT)
!-----------------------------------------------------------------------
!.compute G(H2O) at (T,Pref)
!.Holland and Powell (1991,1998), from CdeCapitani codes
!-----------------------------------------------------------------------
  implicit none
  real(dp),intent(in) :: Tk
  real(dp),intent(out):: GT
  !
  real(dp):: K1,K2,K3,K4
  real(dp):: H0,S0,T0,TT0
  real(dp):: TT,SQT,CPRDT,CPRTDT,SQT0
  data K1,K2,K3,K4 /0.0401D0, 0.8656D-5,  487.5D0, -0.2512D0/
  data H0,S0       /-241.81D0,188.80D-3/
  data T0,TT0      /298.15D0, 88893.4225D0/
  !
  SQT0=  SQRT(T0)
  TT=    Tk*Tk
  SQT=   SQRT(Tk)
  !
  CPRDT= K1*(Tk     -T0)         &
  &    + K2*(TT     -TT0)   /Two &
  &    - K3*(One/Tk -One/T0)     &
  &    + K4*(SQT    -SQT0)  *Two
  CPRTDT=K1 *log(Tk   /T0)      &
  &    + K2 *(Tk      -T0)      &
  &    - K3 *(One/TT  -One/TT0) /Two &
  &    - K4 *(One/SQT -One/SQT0)*Two
  !
  GT= H0 +CPRDT -Tk*(S0 +CPRTDT) 
  GT= GT*1.0D3 !---------------------------------------------------Joule
  !
  return
end subroutine CalcGT_H2O_HP91

!---

end module M_Eos_H2O_Holland_Powell
