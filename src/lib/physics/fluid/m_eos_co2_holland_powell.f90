MODULE M_Eos_CO2_Holland_Powell
!! was M_Holland_Powell_CO2

!-----------------------------------------------------------------------
! CO2 HP91      : Eos_CO2_HolPow91
!-----------------------------------------------------------------------
! Purpose : Holland and Powel EOS for Pure CO2 - Liquid and Vapor!
! Reference : Holland and Powell (1991,1998)
!
! Source : CdeCapitani Theriak [ fsol.f90 ] 
!
!-----------------------------------------------------------------------
! Arxim Integration : J.Moutte
!-----------------------------------------------------------------------

  USE M_Kinds

  IMPLICIT NONE

  PRIVATE

  !// PUBLIC FUNCTIONS
  PUBLIC::  Eos_CO2_HolPow91 ! was CalcGCO2_HP91

  PRIVATE:: CalcGT_CO2_HP91
  
CONTAINS

SUBROUTINE Eos_CO2_HolPow91(TdgK,Pbar,G,V_m3)
!-----------------------------------------------------------------------
!.Holland and Powell (1991,1998), from CdeCapitani codes
!-----------------------------------------------------------------------
  USE M_Eos_Utils,ONLY: CUBIC
  !
  REAL(dp),INTENT(IN) :: TdgK,Pbar
  REAL(dp),INTENT(OUT):: G,V_m3
  !
  REAL(dp):: T,P0,R,RT
  REAL(dp):: AVir,AVir0,AVirT,BVir,BVir0,BVirT,CVir
  REAL(dp):: B,A,A0,A1,A2,A3,GT,SQT,P,VMRK
  REAL(dp):: CubB,CubC,CubD
  REAL(dp):: X1,X2,X2I,X3
  REAL(dp):: dGGas3,dGGas,VREF,Pref,Vcm3
  !
  !C,C0,C1,D,D0,D1,RTLNP
  AVir0= 5.40776D-3  ;  AVirT=-1.59046D-6
  BVir0=-1.78198D-1  ;  BVirT= 2.45317D-5
  CVir=Zero
  !
  R=     8.3142D-3   ;  P0=    5.0D0
  A0=    741.2D0     ;  B=     3.057D0
  A1=   -0.10891D0   ;  A2=   -3.4203D-4    ;  A3=Zero
  !
  !Chp91 DATA C0,C1/-2.26924D-1,7.73793D-5/
  !Chp91 DATA D0,D1/1.33790D-2,-1.01740D-5/
  
  T=     TdgK
  P=     Pbar/1000.0D0 !P in kbar
  Pref=  0.001D0       !Pref= 0.001 kbar
  
  G=     Zero
  dGGas= Zero
  dGGas3=Zero
  RT=    R*T
  SQT=   SQRT(T)
  
  AVir=  AVir0 + AVirT*T !                            Chp91 AVir=D0+D1*T
  BVir=  BVir0 + BVirT*T !                            Chp91 BVir=C0+C1*T
  A=     A0+ A1*T +A2*T*T +A3*T*T*T
  
  VRef=  RT/Pref !                       reference volume at 1 Bar and T
  
  !                                               calc volume at P and T
  CubB= -RT/P
  CubC= -(B*RT+B*B*P-A/SQT)/P
  CubD= -A*B/SQT/P
  
  CALL CUBIC(CubB,CubC,CubD,X1,X2,X2I,X3)
  
  IF (X2I/=Zero) THEN
    VMRK=X1
  ELSE
    VMRK=MIN(X1,X2,X3)
    IF (VMRK<B) VMRK=MAX(X1,X2,X3)
  ENDIF
  
  IF (P>P0) THEN
    Vcm3= VMRK + AVir*(P-P0) + BVir*SQRT(P-P0) + CVir*(P-P0)**0.25D0
  ELSE
    Vcm3= VMRK
  ENDIF
  
  V_m3= Vcm3 /1.0D6
  
  dGGas3= VMRK*P-VREF*Pref &
  &     - RT*LOG((VMRK-B)/(VREF-B)) &
  &     + A/SQT/B *LOG(VMRK*(VREF+B)/VREF/(VMRK+B))
  
  IF (P>P0) &
  & dGGas3=  dGGas3 &
  &       + AVir/Two       *(P-P0)**2 &
  &       + BVir*Two/3.0D0 *(P-P0)**1.5D0 &
  &       + CVir*0.8D0     *(P-P0)**1.25D0
  
  dGGas= dGGas3 *1.0D3
  
  CALL CalcGT_CO2_HP91(T,GT)
  
  !RTLNP=1D3*RT*LOG(P*1D3)
  G=    GT +dGGas
  
  RETURN
ENDSUBROUTINE Eos_CO2_HolPow91

!---

SUBROUTINE CalcGT_CO2_HP91(T,GT) 
!-----------------------------------------------------------------------
!Holland and Powell (1991), from CdeCapitani codes
!-----------------------------------------------------------------------
  REAL(dp),INTENT(IN) :: T
  REAL(dp),INTENT(OUT):: GT
  !
  REAL(dp):: CPRDT,CPRTDT
  REAL(dp):: TT,SQT,SQT0
  REAL(dp):: H0=-393.51D0, S0=  213.70D-3
  REAL(dp):: K1= 0.0878D0, K2= -0.2644D-5,   K3= 706.4D0, K4=-0.9989D0
  REAL(dp):: T0= 298.15D0, TT0= 88893.4225D0
  
  SQT0=SQRT(T0)
  TT=T*T
  SQT=SQRT(T)
  !
  CPRDT=  K1 *(T     -T0    )     &
  &     + K2 *(TT    -TT0   )/Two &
  &     - K3 *(One/T -One/T0)     &
  &     + K4 *(SQT   -SQT0  )*Two
  CPRTDT= K1 *LOG(T    /T0)       &
  &     + K2 *(T       -T0      ) &
  &     - K3 *(One/TT  -One/TT0 )/Two &
  &     - K4 *(One/SQT -One/SQT0)*Two
  !
  GT= (H0+CPRDT-T*(S0+CPRTDT)) *1.0D3
  
  RETURN
ENDSUBROUTINE CalcGT_CO2_HP91

END MODULE M_Eos_CO2_Holland_Powell
