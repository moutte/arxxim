module M_Eos_CO2_Holland_Powell
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

  use M_Kinds

  implicit none

  private

  !// public functionS
  public::  Eos_CO2_HolPow91 ! was CalcGCO2_HP91

  private:: CalcGT_CO2_HP91
  
contains

subroutine Eos_CO2_HolPow91(TdgK,Pbar,Grt,V_m3,LnFug)
!-----------------------------------------------------------------------
!.Holland and Powell (1991,1998), from CdeCapitani codes
!-----------------------------------------------------------------------
  use M_Eos_Utils,only: CUBIC
  !
  real(dp),intent(in) :: TdgK,Pbar
  real(dp),intent(out):: Grt,V_m3,LnFug
  !
  real(dp):: T,P0,R,RT
  real(dp):: AVir,AVir0,AVirT,BVir,BVir0,BVirT,CVir
  real(dp):: B,A,A0,A1,A2,A3,GT,SQT,P,VMRK
  real(dp):: CubB,CubC,CubD
  real(dp):: X1,X2,X2I,X3
  real(dp):: dGGas3,dGGas,VREF,Pref,Vol
  real(dp):: G,S0Ele
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
  !Chp91 data C0,C1/-2.26924D-1,7.73793D-5/
  !Chp91 data D0,D1/1.33790D-2,-1.01740D-5/
  
  T=     TdgK
  P=     Pbar/1000.0D0 !P in kbar
  Pref=  0.001D0       !Pref= 0.001 kbar
  
  G=     Zero
  dGGas= Zero
  dGGas3=Zero
  RT=    R*T
  SQT=   sqrt(T)
  
  AVir=  AVir0 + AVirT*T !----------------------------Chp91 AVir=D0+D1*T
  BVir=  BVir0 + BVirT*T !----------------------------Chp91 BVir=C0+C1*T
  A=     A0+ A1*T +A2*T*T +A3*T*T*T
  
  VRef=  RT/Pref !-----------------------reference volume at 1 Bar and T
  
  !-----------------------------------------------calc volume at P and T
  CubB= -RT/P
  CubC= -(B*RT+B*B*P-A/SQT)/P
  CubD= -A*B/SQT/P
  
  call CUBIC(CubB,CubC,CubD,X1,X2,X2I,X3)
  
  if (X2I/=Zero) then
    VMRK=X1
  else
    VMRK=MIN(X1,X2,X3)
    if (VMRK<B) VMRK=MAX(X1,X2,X3)
  end if
  
  if (P>P0) then
    Vol= VMRK + AVir*(P-P0) + BVir*sqrt(P-P0) + CVir*(P-P0)**0.25D0
  else
    Vol= VMRK
  end if
  
  V_m3= Vol *1.0D-5
  
  dGGas3= VMRK*P-VREF*Pref &
  &     - RT*log((VMRK-B)/(VREF-B)) &
  &     + A/SQT/B *log(VMRK*(VREF+B)/VREF/(VMRK+B))
  
  if (P>P0) &
  & dGGas3=  dGGas3 &
  &       + AVir/Two       *(P-P0)**2 &
  &       + BVir*Two/3.0D0 *(P-P0)**1.5D0 &
  &       + CVir*0.8D0     *(P-P0)**1.25D0
  
  dGGas= dGGas3*1.0D3
  
  call CalcGT_CO2_HP91(T,GT)
  
  !RTLNP=1D3*RT*log(P*1D3)
  G=    GT +dGGas
  
  ! Benson convention
  S0Ele= 5.74D0  + 102.575D0 *2.0D0
  G= G + S0Ele*298.15D0
  Grt= G/8.314502D0/TdgK
  LnFug= dGGas/8.314502D0/TdgK
  
  return
end subroutine Eos_CO2_HolPow91

!---

subroutine CalcGT_CO2_HP91(T,GT) 
!-----------------------------------------------------------------------
!Holland and Powell (1991), from CdeCapitani codes
!-----------------------------------------------------------------------
  real(dp),intent(in) :: T
  real(dp),intent(out):: GT
  !
  real(dp):: CPRDT,CPRTDT
  real(dp):: TT,SQT,SQT0
  real(dp):: H0=-393.51D0, S0=  213.70D-3
  real(dp):: K1= 0.0878D0, K2= -0.2644D-5,   K3= 706.4D0, K4=-0.9989D0
  real(dp):: T0= 298.15D0, TT0= 88893.4225D0
  
  SQT0=sqrt(T0)
  TT=T*T
  SQT=sqrt(T)
  !
  CPRDT=  K1 *(T     -T0    )     &
  &     + K2 *(TT    -TT0   )/Two &
  &     - K3 *(One/T -One/T0)     &
  &     + K4 *(SQT   -SQT0  )*Two
  CPRTDT= K1 *log(T    /T0)       &
  &     + K2 *(T       -T0      ) &
  &     - K3 *(One/TT  -One/TT0 )/Two &
  &     - K4 *(One/SQT -One/SQT0)*Two
  !
  GT= (H0+CPRDT-T*(S0+CPRTDT)) *1.0D3
  
  return
end subroutine CalcGT_CO2_HP91

end module M_Eos_CO2_Holland_Powell
