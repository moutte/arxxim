module M_Mixmodel_Duan
!--
!-- interface between Mixture model and optimization routine
!--
  use M_Kinds

  implicit none
  
  private
  
  public:: EoS_H2O_CO2_Duan92
  public:: EoS_CO2_Duan92
  public:: EoS_H2O_Duan92
  
  public:: EoS_H2O_CO2_Duan06
  public:: EoS_H2O_Duan06_LP
  public:: EoS_CO2_Duan06_LP
  public:: EoS_H2O_Duan06_HP
  public:: EoS_CO2_Duan06_HP
  
  !! if(SOLNAM=='HC-DU92') then
  !!   F1PLO=0.0D0
  !!   F2PLO=0.0D0
  !!   call EoS_H2O_CO2_Duan92(T,P,X,VAT,AC1,AC2,F1PLO,F2PLO)
  !!   A(1)=AC1
  !!   A(2)=AC2
  !!   return
  !! end if
  !!
  !! if (FALL=='DU92H2O') call EoS_H2O_Duan92(T,PGAS,G,V,F1)
  !! if (FALL=='DU92CO2') call EoS_CO2_Duan92(T,PGAS,G,V,F1)
  
  real(dp):: BVC,CVC2,DVC4,EVC5,FVC2,BETA1,GAMMAVC2
  real(dp):: R,RT,MWT
  real(dp):: BST1,CST1,DST1,EST1,FST1,BESTST1,GAMST1
  real(dp):: BST2,CST2,DST2,EST2,FST2,BESTST2,GAMST2
  
  !! real*8 BVC,CVC2,DVC4,EVC5,FVC2,BETA1,GAMMAVC2,R,RT,MWT, &
  !! BST1,CST1,DST1,EST1,FST1,BESTST1,GAMST1, &
  !! BST2,CST2,DST2,EST2,FST2,BESTST2,GAMST2
  !! common /DTHINGS/ BVC,CVC2,DVC4,EVC5,FVC2,BETA1,GAMMAVC2,R,RT, &
  !! MWT, &
  !! BST1,CST1,DST1,EST1,FST1,BESTST1,GAMST1, &
  !! BST2,CST2,DST2,EST2,FST2,BESTST2,GAMST2
  
contains

subroutine EoS_H2O_Duan92( &  !
& T,P, &                      !IN
& Grt,V_m3,LnPhi)             !OUT
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--
  real(dp),intent(in) :: T,P
  real(dp),intent(out):: Grt,V_m3,LnPhi
  !---------------------------------------------------------------------
  real(dp):: &
  & PR,TR,TR2,TR3, &
  & VC,VC2,VC4,VC5,EU3,Z,V,INTE, &
  & V2,V4,V5,EXPON, &
  & G,H0,S0,S0Ele, &
  & K1,K2,K3,K4,K6, &
  & CPDT,CPTDT,T0,TT,TT0,SQT,SQT0
  integer:: NP,NSTAB
  !---------------------------------------------------------------------
  T0=298.15D0
  TT=T*T
  TT0=T0*T0
  SQT=sqrt(T)
  SQT0=sqrt(T0)

  R=0.08314467D0 ! R_JK=8.314502D0 
  RT=R*T
  
  EU3=1.0D0/3.0D0
  
  ! PR=P/221.19D0
  TR=T/647.25D0
  
  TR2=TR*TR
  TR3=TR2*TR
  
  VC=R*647.25D0/221.19D0
  VC2=VC*VC
  VC4=VC2*VC2
  VC5=VC4*VC
  
  BVC= ( 8.64449220D-02 -3.96918955D-01/TR2  -5.73334886D-02/TR3)*VC
  CVC2=(-2.93893000D-04 -4.15775512D-03/TR2  +1.99496791D-02/TR3)*VC2
  DVC4=( 1.18901426D-04 +1.55212063D-04/TR2  -1.06855859D-04/TR3)*VC4
  EVC5=(-4.93197687D-06 -2.73739155D-06/TR2  +2.65571238D-06/TR3)*VC5
  FVC2=(8.96079018D-03/TR3)*VC2
  
  BETA1=4.02D+00
  GAMMAVC2=2.57D-02*VC2
  MWT=0.0180154D0  ! molar weight kg/mole
  
  BST1=2.0D0*BVC
  CST1=3.0D0*CVC2
  DST1=5.0D0*DVC4
  EST1=6.0D0*EVC5
  FST1=2.0D0*FVC2
  BESTST1=BETA1
  GAMST1=3.0D0*GAMMAVC2
  
  BST2=0.0D0
  CST2=0.0D0
  DST2=0.0D0
  EST2=0.0D0
  FST2=0.0D0
  BESTST2=0.0D0
  GAMST2=0.0D0
 
  !----------------------------------------------------- find the volume
  call DuanZ( &
  & T,P,      &
  & NP,NSTAB, &
  & Z,V,INTE)
  
  !--------------------------------------------------------- EQ9 for H2O
  V2=V**2
  V4=V2**2
  V5=V4*V
  EXPON=exp(-GAMMAVC2/V2)
  lnPhi= -log(Z) &
  &   + BST1/V &
  &   + CST1/(2.0D0*V2) &
  &   + DST1/(4.0D0*V4) &
  &   + EST1/(5.0D0*V5) &
  &   +((FST1*BETA1+BESTST1*FVC2)/(2.0D0*GAMMAVC2))*(1.0D0-EXPON) &
  &   +((FST1*GAMMAVC2 +GAMST1*FVC2 -FVC2*BETA1*(GAMST1-GAMMAVC2)) &
  &      /(2.0D0*GAMMAVC2**2)) &
  &    *(1.0D0-((GAMMAVC2/V2)+1.0D0)*EXPON) &
  &   -(((GAMST1-GAMMAVC2)*FVC2)/(2.0D0*GAMMAVC2**2)) &
  &    *(2.-(GAMMAVC2**2/V4+2.0D0*GAMMAVC2/V2+2.0D0)*EXPON)
      
  !----------------------------------------------------- G for ideal gas
  H0=-241816.00D0 
  S0=188.7200D0
  K1=115.45000D0
  K4=-3799.900D0
  K3=-2871300.00D0
  K6=51055.16D0
  K2=-0.002062D0
  !
  CPDT= K1*(T-T0) &
  &   + K2*(TT-TT0)/2.0D0 &
  &   - K3*(1.0D0/T-1.0D0/T0) &
  &   + K4*2D0*(SQT-SQT0) &
  &   + K6*log(T/T0)
  CPTDT=K1*log(T/T0) &
  &    +K2*(T-T0) &
  &    -K3*(1.0D0/(TT)-1.0D0/(TT0))/2.0D0 &
  &    -K4*2D0*(1.0D0/SQT-1.0D0/SQT0) &
  &    -K6*(1.0D0/T-1.0D0/T0)
  G=H0+CPDT-T*(S0+CPTDT)
  
  !----------------------------------------------------Benson convention
  S0Ele= 102.575D0 + 65.34D0 *2.0D0 !H2O, joule
  G= G + S0Ele*T0

  !------------------------------------------------ correct for fugacity
  G=       G + R*100.0D0 *T *(lnPhi + log(P))
  Grt=     G/8.314502D0/T
  V_m3=    V*1.0D-3

  ! Michelsen Mollerup eq 86
  ! for pure fluid,
  ! RT.ln(fug(T,P) / P0)= mu(T,P) - mu*(T,P0)
  ! mu*(T,P0) is for pure perfect gas at T,P0
  ! mu*(T,P) - mu*(T,P0)= RT.ln(P / P0)
  !
  ! P0= 1
  ! mu(T,P)= mu*(T,P0) + RT.ln(fug(T,P))
  !        = mu*(T,P0) + RT.ln(phi(T,P)) + RT.ln(P)
  
  return
end subroutine EoS_H2O_Duan92

subroutine EoS_CO2_Duan92( &
& T,P, &
& Grt,V_m3,LnPhi)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--
  implicit none
  !---------------------------------------------------------------------
  real(dp),intent(in) :: T,P
  real(dp),intent(out):: Grt,V_m3,lnPhi
  !---------------------------------------------------------------------
  real(dp):: &
  & TR,TR2,TR3,  &
  & VC,VC2,VC4,VC5, &
  & EU3,Z,V,INTE,   &
  & V2,V4,V5,EXPON,G,H0,S0,S0Ele, &
  & K1,K2,K3,K4,K6, &
  & CPDT,CPTDT,T0,TT,TT0,SQT,SQT0
  integer:: NP,NSTAB
  !---------------------------------------------------------------------
  T0=298.15D0
  TT=T*T
  TT0=T0*T0
  SQT=sqrt(T)
  SQT0=sqrt(T0)
  
  R=0.08314467D0
  RT=R*T
  EU3=1.0D0/3.0D0
  
  TR=T/304.20D0
  TR2=TR*TR
  TR3=TR2*TR
  
  VC=R*304.20D0/73.825D0
  VC2=VC*VC
  VC4=VC2*VC2
  VC5=VC4*VC
  
  !! common !!
  BVC= ( 8.99288497D-02 - 4.94783127D-01/TR2 + 4.77922245D-02/TR3)*VC
  CVC2=( 1.03808883D-02 - 2.82516861D-02/TR2 + 9.49887563D-02/TR3)*VC2
  DVC4=( 5.20600880D-04 - 2.93540971D-04/TR2 - 1.77265112D-03/TR3)*VC4
  EVC5=(-2.51101973D-05 + 8.93353441D-05/TR2 + 7.88998563D-05/TR3)*VC5
  FVC2=(                                     - 1.66727022D-02/TR3)*VC2
  BETA1=1.398D0
  GAMMAVC2=2.96D-02*VC2
  
  MWT=0.0440098D0
  
  BST2=2.0D0*BVC
  CST2=3.0D0*CVC2
  DST2=5.0D0*DVC4
  EST2=6.0D0*EVC5
  FST2=2.0D0*FVC2
  BESTST2=BETA1
  GAMST2=3.0D0*GAMMAVC2
  
  BST1=0.0D0
  CST1=0.0D0
  DST1=0.0D0
  EST1=0.0D0
  FST1=0.0D0
  BESTST1=0.0D0
  GAMST1=0.0D0
  !----------------------------------------------------- find the volume
  call DuanZ( &
  & T,P,      &
  & NP,NSTAB, &
  & Z,V,INTE)
  !
  !--------------------------------------------------------- EQ9 for CO2
  V2=V*V
  V4=V2*V2
  V5=V4*V
  EXPON= exp(-GAMMAVC2/V2)
  lnPhi= -log(Z) &
  &   +  BST2/V &
  &   +  CST2/(2.0D0*V2) &
  &   +  DST2/(4.0D0*V4) &
  &   +  EST2/(5.0D0*V5) &
  &   + ((FST2*BETA1+BESTST2*FVC2)/2.0D0/GAMMAVC2)*(1.0D0-EXPON) &
  &   + ((FST2*GAMMAVC2+GAMST2*FVC2-FVC2*BETA1*(GAMST2-GAMMAVC2)) &
  &     /(2.0D0*GAMMAVC2**2)) &
  &     *(1.0D0-(GAMMAVC2/V2+1.0D0)*EXPON) &
  &   - (((GAMST2-GAMMAVC2)*FVC2)/(2.0D0*GAMMAVC2**2)) &
  &     *(2.-(GAMMAVC2**2/V4+2.0D0*GAMMAVC2/V2+2.0D0)*EXPON)  
  !
  !----------------------------------------------------- G for ideal gas
  H0= -393510.010
  S0=  213.6770D0
  K1=  93.0000D0
  K4= -1340.900D0
  K3=  123800.000D0
  K6=  6336.20D0
  K2= -0.002876D0
  !
  CPDT= K1*(T-T0) &
  &   + K2*(TT-TT0)/2.0D0 &
  &   - K3*(1.0D0/T-1.0D0/T0) &
  &   + K4*2D0*(SQT-SQT0) &
  &   + K6*log(T/T0)
  !
  CPTDT= K1*log(T/T0) &
  &    + K2*(T-T0) &
  &    - K3*(1.0D0/(TT)-1.0D0/(TT0))/2.0D0 &
  &    - K4*2D0*(1.0D0/SQT-1.0D0/SQT0) &
  &    - K6*(1.0D0/T-1.0D0/T0)
  !
  G=H0+CPDT-T*(S0+CPTDT)
  !-----------------------------------------------------/G for ideal gas
  
  !----------------------------------------------------Benson convention
  S0Ele= 5.74D0  + 102.575D0 *2.0D0
  G= G + S0Ele*T0

  !------------------------------------------------ correct for fugacity
  ! mu(T,P)= mu*(T,P0) + RT.ln(fug(T,P))
  !        = mu*(T,P0) + RT.ln(phi(T,P)) + RT.ln(P)
  G=       G + R*100.0D0 *T *(lnPhi + log(P))
  !FugCoef=exp(lnPhi)
  !G= G + R*100.0D0 *T *(lnPhi +log(P))
  Grt= G/8.314502D0/T
  V_m3=V*1.0D-3
  
  return
end subroutine EoS_CO2_Duan92

subroutine EoS_H2O_CO2_Duan92( & !
& T,P,XF,                    & ! IN
& UpdatePhi,                 & ! IN
& LnPhiH2O_Pur,LnPhiCO2_Pur, & ! INOUT
& V_H2O_Pur,   V_CO2_Pur,    & ! INOUT
& V_m3,LnActH2O,LnActCO2)      ! OUT
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--
  real(dp),intent(in) :: T,P
  real(dp),intent(in) :: XF(2)
  logical, intent(in) :: UpdatePhi
  real(dp),intent(inout):: LnPhiH2O_Pur,LnPhiCO2_Pur
  real(dp),intent(inout):: V_H2O_Pur,V_CO2_Pur
  real(dp),intent(out):: V_m3
  real(dp),intent(out):: LnActH2O,LnActCO2
  !--------------------------------------------------------------------- 
  real(dp):: VATP,INTATP,ZATP
  real(dp):: LnPhiH2O,LnPhiCO2,G
  integer :: NP,NSTAB
  !---------------------------------------------------------------------
  if(UpdatePhi) then
    call EoS_H2O_Duan92(T,P,   G,V_H2O_Pur,LnPhiH2O_Pur)
    call EoS_CO2_Duan92(T,P,   G,V_CO2_Pur,LnPhiCO2_Pur)
  end if
  ! print *,FugH2O,FugCO2

  call MixRules92(T,P,           XF)
  call DuanZ     (T,P,           NP,NSTAB,ZATP,VATP,INTATP)
  call EQUA9     (VATP,ZATP,     LnPhiH2O,LnPhiCO2)

  write(91,'(A,2I7)') "NP,NSTAB",NP,NSTAB
  
  V_m3=VATP*1.0D-3
  LnActH2O= log(XF(1)) +LnPhiH2O -LnPhiH2O_Pur
  LnActCO2= log(XF(2)) +LnPhiCO2 -LnPhiCO2_Pur
  
  ! Michelsen Mollerup eq 119 - 122
  ! fug_i_mix(T,P,x)= phi_i_mix(T,P,x).P.x_i
  !                 = fug_i_pur(T,P).gamma_i(T,P,x).x_i
  !
  ! gamma_i(T,P,x)= activ_i(T,P,x) / x_i
  !               = fug_i_mix(T,P,x) / fug_i_pur(T,P)
  ! 
  ! activ_i(T,P,x)= x_i * fug_i_mix(T,P,x) / fug_i_pur(T,P) 

  return
end subroutine EoS_H2O_CO2_Duan92

subroutine EQUA9( &
& V,Z, &
& lnPhi1,lnPhi2)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!-- 
!--
  real(dp),intent(in) :: V,Z
  real(dp),intent(out):: lnPhi1,lnPhi2
  !---------------------------------------------------------------------
  real(dp):: V2,V4,V5,EXPON
  !---------------------------------------------------------------------
  V2=V*V
  V4=V2*V2
  V5=V4*V
  
  EXPON=exp(-GAMMAVC2/V2)
  
  lnPhi1=-log(Z) &
  &   + BST1/V &
  &   + CST1/(2.0D0*V2) &
  &   + DST1/(4.0D0*V4) &
  &   + EST1/(5.0D0*V5) &
  &   +((FST1*BETA1+BESTST1*FVC2)/(2.0D0*GAMMAVC2))*(1.0D0-EXPON) &
  &   +((FST1*GAMMAVC2+GAMST1*FVC2-FVC2*BETA1*(GAMST1-GAMMAVC2))  &
  &     /(2.0D0*GAMMAVC2**2)) &
  &    *(1.0D0-((GAMMAVC2/V2)+1.0D0)*EXPON) &
  &   -(((GAMST1-GAMMAVC2)*FVC2)/(2.0D0*GAMMAVC2**2)) &
  &    *(2.-(GAMMAVC2**2/V4+2.0D0*GAMMAVC2/V2+2.0D0)*EXPON)

  lnPhi2=-log(Z) &
  &   + BST2/V  &
  &   + CST2/(2.0D0*V2) &
  &   + DST2/(4.0D0*V4) &
  &   + EST2/(5.0D0*V5) &
  &   +((FST2*BETA1+BESTST2*FVC2)/(2.0D0*GAMMAVC2))*(1.0D0-EXPON) &
  &   +((FST2*GAMMAVC2+GAMST2*FVC2-FVC2*BETA1*(GAMST2-GAMMAVC2))  &
  &     /(2.0D0*GAMMAVC2**2)) &
  &    *(1.0D0-((GAMMAVC2/V2)+1.0D0)*EXPON) &
  &   -(((GAMST2-GAMMAVC2)*FVC2)/(2.0D0*GAMMAVC2**2)) &
  &    *(2.-(GAMMAVC2**2/V4+2.0D0*GAMMAVC2/V2+2.0D0)*EXPON)
  
  return
end subroutine EQUA9

subroutine MixRules92(T,P,XF)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--

  real(dp),intent(in):: T,P
  real(dp),intent(in):: XF(2)
  
  !! real(dp):: BVC,CVC2,DVC4,EVC5,FVC2,BETA1,GAMMAVC2,R,RT,MWT, &
  !! BST1,CST1,DST1,EST1,FST1,BESTST1,GAMST1, &
  !! BST2,CST2,DST2,EST2,FST2,BESTST2,GAMST2
  !! common /DTHINGS/ &
  !! & BVC,CVC2,DVC4,EVC5,FVC2,BETA1,GAMMAVC2,R,RT,MWT, &
  !! & BST1,CST1,DST1,EST1,FST1,BESTST1,GAMST1, &
  !! & BST2,CST2,DST2,EST2,FST2,BESTST2,GAMST2
  
  !------------------------------------------------- 1=H2O, 2=CO2, 3=CH4
  real(dp):: &
  & A1(2),A2(2),A3(2),A4(2),A5(2),A6(2), &
  & A7(2),A8(2),A9(2),A10(2),A11(2),A12(2), &
  & ALPHA(2),BETA(2),GAMMA(2), &
  & TC(2),PC(2), &
  & BB(2),CC(2),DD(2),EE(2),FF(2), &
  & BS(2),CS(2),DS(2),ES(2),FS(2), &
  & PM(2), &
  & VC(2),VC2(2),VC4(2),VC5(2), &
  & TR(2),TR2(2),TR3(2),PR(2),EU3
  
  integer:: I,J,K,II,L,M,N
  
  real(dp):: &
  & K1(2,2),K2(2,2,2),K3(2,2,2), &
  & BIJ(2,2),VCIJ(2,2),CIJK(2,2,2),VCIJK(2,2,2), &
  & DIJKLM(2,2,2,2,2),VCIJKLM(2,2,2,2,2), &
  & EIJKLMN(2,2,2,2,2,2),VCIJKLMN(2,2,2,2,2,2), &
  & FIJ(2,2),GAMIJK(2,2,2)
  
  data A1  / 8.64449220D-02, 8.99288497D-02/
  data A2  /-3.96918955D-01,-4.94783127D-01/
  data A3  /-5.73334886D-02, 4.77922245D-02/
  data A4  /-2.93893000D-04, 1.03808883D-02/
  data A5  /-4.15775512D-03,-2.82516861D-02/
  data A6  / 1.99496791D-02, 9.49887563D-02/
  data A7  / 1.18901426D-04, 5.20600880D-04/
  data A8  / 1.55212063D-04,-2.93540971D-04/
  data A9  /-1.06855859D-04,-1.77265112D-03/
  data A10 /-4.93197687D-06,-2.51101973D-05/
  data A11 /-2.73739155D-06, 8.93353441D-05/
  data A12 / 2.65571238D-06, 7.88998563D-05/
  
  data ALPHA / 8.96079018D-03,-1.66727022D-02/
  data BETA  / 4.02D+00,       1.398D+00/
  data GAMMA / 2.57D-02,       2.96D-02/
  
  data TC /647.25D0,304.20D0/
  data PC /221.19D0,73.825D0/
  data PM /0.0180154D0,0.0440098D0/
  
  R=0.08314467D0
  RT=R*T
  EU3=1.0D0/3.0D0
  II=1
  
  do I=1,2
  
    PR(I)=P/PC(I)
    
    TR(I)=T/TC(I)
    TR2(I)=TR(I)**2
    TR3(I)=TR(I)**3

    VC(I)=R*TC(I)/PC(I)
    VC2(I)=VC(I)**2
    VC4(I)=VC(I)**4
    VC5(I)=VC(I)**5
    
    BB(I)= A1(I)  +A2(I)/TR2(I)  +A3(I)/TR3(I)
    CC(I)= A4(I)  +A5(I)/TR2(I)  +A6(I)/TR3(I)
    DD(I)= A7(I)  +A8(I)/TR2(I)  +A9(I)/TR3(I)
    EE(I)= A10(I) +A11(I)/TR2(I) +A12(I)/TR3(I)
    FF(I)= ALPHA(I)/TR3(I)
    
    if (BB(I) < 0.0D0) then
      BB(I)=-BB(I)
      BS(I)=-1.0D0
    else
      BS(I)=1.0D0
    end if
    
    if (CC(I) < 0.0D0) then
      CC(I)=-CC(I)
      CS(I)=-1.0D0
    else
      CS(I)=1.0D0
    end if
    
    if (DD(I) < 0.0D0) then
      DD(I)=-DD(I)
      DS(I)=-1.0D0
    else
      DS(I)=1.0D0
    end if
    
    if (EE(I) < 0.0D0) then
      EE(I)=-EE(I)
      ES(I)=-1.0D0
    else
      ES(I)=1.0D0
    end if
    
    if (FF(I) < 0.0D0) then
      FF(I)=-FF(I)
      FS(I)=-1.0D0
    else
      FS(I)=1.0D0
    end if
    
  end do
  
  call MAKEK1_92(T,K1,K2,K3)
  
  do I=1,2
    do J=1,2
      BIJ(I,J)=(((BS(I)*BB(I)**EU3 +BS(J)*BB(J)**EU3)/2.0D0)**3)*K1(I,J)
      VCIJ(I,J)=((VC(I)**EU3+VC(J)**EU3)/2.0D0)**3
    end do
  end do
  
  BVC=0.0D0
  do I=1,2
    do J=1,2
      BVC=BVC+XF(I)*XF(J)*BIJ(I,J)*VCIJ(I,J)
    end do
  end do
  
  BST1=0.0D0
  BST2=0.0D0
  do J=1,2
    BST1=BST1+2.0D0*XF(J)*BIJ(1,J)*VCIJ(1,J)
    BST2=BST2+2.0D0*XF(J)*BIJ(J,2)*VCIJ(J,2)
  end do
  
  do I=1,2
    do J=1,2
      do K=1,2
        CIJK(I,J,K)=( (( CS(I)*CC(I)**EU3 &
        &              + CS(J)*CC(J)**EU3 &
        &              +       CC(K)**EU3 )/3.0D0)**3) &
        &           *K2(I,J,K)
        VCIJK(I,J,K)=( (VC(I)**EU3 +VC(J)**EU3 +VC(K)**EU3)/3.0D0 )**3
      end do
    end do
  end do
  
  CVC2=0.0D0
  do I=1,2
    do J=1,2
      do K=1,2
        CVC2=CVC2 +XF(I)*XF(J)*XF(K)*CIJK(I,J,K)*VCIJK(I,J,K)**2
      end do
    end do
  end do
  
  CST1=0.0D0
  CST2=0.0D0
  do J=1,2
    do K=1,2
      CST1=CST1 +3.0D0*XF(J)*XF(K)*CIJK(1,J,K)*VCIJK(1,J,K)**2
      CST2=CST2 +3.0D0*XF(J)*XF(K)*CIJK(J,2,K)*VCIJK(J,2,K)**2
    end do
  end do
  
  do I=1,2
    do J=1,2
      do K=1,2
        do L=1,2
          do M=1,2
            DIJKLM(I,J,K,L,M)= &
            & (( DS(I)*DD(I)**EU3 &
            &  + DS(J)*DD(J)**EU3 &
            &  + DS(K)*DD(K)**EU3 &
            &  + DS(L)*DD(L)**EU3 &
            &  + DS(M)*DD(M)**EU3 )/5.0D0)**3
            VCIJKLM(I,J,K,L,M)= &
            & (( VC(I)**EU3 &
            &  + VC(J)**EU3 &
            &  + VC(K)**EU3 &
            &  + VC(L)**EU3 &
            &  + VC(M)**EU3 )/5.0D0)**3
          end do
        end do
      end do
    end do
  end do
  DVC4=0.0D0
  do I=1,2
    do J=1,2
      do K=1,2
        do L=1,2
          do M=1,2
            DVC4= DVC4 &
            &   + XF(I)*XF(J)*XF(K)*XF(L)*XF(M) &
            &          *DIJKLM(I,J,K,L,M) &
            &          *VCIJKLM(I,J,K,L,M)**4
          end do
        end do
      end do
    end do
  end do
  DST1=0.0D0
  DST2=0.0D0
  do J=1,2
    do K=1,2
      do L=1,2
        do M=1,2
          DST1=DST1 &
          &   +5.0D0*XF(J)*XF(K)*XF(L)*XF(M) &
          &         *DIJKLM(1,J,K,L,M)*VCIJKLM(1,J,K,L,M)**4
          DST2=DST2 &
          &   +5.0D0*XF(J)*XF(K)*XF(L)*XF(M) &
          &         *DIJKLM(J,2,K,L,M)*VCIJKLM(J,2,K,L,M)**4
        end do
      end do
    end do
  end do
  !---
  do I=1,2
    do J=1,2
      do K=1,2
        do L=1,2
          do M=1,2
            do N=1,2
              EIJKLMN(I,J,K,L,M,N)= &
              & ((ES(I)*EE(I)**EU3 &
              &  +ES(J)*EE(J)**EU3 &
              &  +ES(K)*EE(K)**EU3 &
              &  +ES(L)*EE(L)**EU3 &
              &  +ES(M)*EE(M)**EU3 &
              &  +ES(N)*EE(N)**EU3)/6.0D0)**3
              VCIJKLMN(I,J,K,L,M,N)= &
              & ((VC(I)**EU3 &
              &  +VC(J)**EU3 &
              &  +VC(K)**EU3 &
              &  +VC(L)**EU3 &
              &  +VC(M)**EU3 &
              &  +VC(N)**EU3)/6.0D0)**3
            end do
          end do
        end do
      end do
    end do
  end do
  
  EVC5=0.0D0
  do I=1,2
    do J=1,2
      do K=1,2
        do L=1,2
          do M=1,2
            do N=1,2
              EVC5= EVC5 &
              &   + XF(I)*XF(J)*XF(K)*XF(L)*XF(M)*XF(N) &
              &    *EIJKLMN(I,J,K,L,M,N)*VCIJKLMN(I,J,K,L,M,N)**5
            end do
          end do
        end do
      end do
    end do
  end do
  
  EST1=0.0D0
  EST2=0.0D0
  do J=1,2
    do K=1,2
      do L=1,2
        do M=1,2
          do N=1,2
            EST1=EST1 &
            &   +6.0D0*XF(J)*XF(K)*XF(L)*XF(M)*XF(N) &
            &         *EIJKLMN(1,J,K,L,M,N)*VCIJKLMN(1,J,K,L,M,N)**5
            EST2=EST2 &
            &   +6.0D0*XF(J)*XF(K)*XF(L)*XF(M)*XF(N) &
            &         *EIJKLMN(J,2,K,L,M,N)*VCIJKLMN(J,2,K,L,M,N)**5
          end do
        end do
      end do
    end do
  end do
  
  do I=1,2
    do J=1,2
      FIJ(I,J)=((FS(I)*FF(I)**EU3+FS(J)*FF(J)**EU3)/2.0D0)**3
    end do
  end do
  
  FVC2=0.0D0
  do I=1,2
    do J=1,2
      FVC2=FVC2+XF(I)*XF(J)*FIJ(I,J)*VCIJ(I,J)**2
    end do
  end do
  
  FST1=0.0D0
  FST2=0.0D0
  do J=1,2
    FST1=FST1+2.0D0*XF(J)*FIJ(1,J)*VCIJ(1,J)**2
    FST2=FST2+2.0D0*XF(J)*FIJ(J,2)*VCIJ(J,2)**2
  end do
  
  BETA1=0.0D0
  do I=1,2
    BETA1=BETA1+XF(I)*BETA(I)
  end do
  
  BESTST1=BETA(1)
  BESTST2=BETA(2)
  
  do I=1,2
    do J=1,2
      do K=1,2
        GAMIJK(I,J,K)=(((GAMMA(I)**EU3 +GAMMA(J)**EU3 +GAMMA(K)**EU3)/3.0D0)**3) &
        &             *K3(I,J,K)
      end do
    end do
  end do
  
  GAMMAVC2=0.0D0
  do I=1,2
    do J=1,2
      do K=1,2
        GAMMAVC2= GAMMAVC2 &
        &       + XF(I)*XF(J)*XF(K)*GAMIJK(I,J,K)*VCIJK(I,J,K)**2
      end do
    end do
  end do
  
  GAMST1=0.0D0
  GAMST2=0.0D0
  do J=1,2
    do K=1,2
      GAMST1=GAMST1 &
      &     +3.0D0*XF(J)*XF(K)*GAMIJK(1,J,K)*VCIJK(1,J,K)**2
      GAMST2=GAMST2 &
      &     +3.0D0*XF(J)*XF(K)*GAMIJK(J,2,K)*VCIJK(J,2,K)**2
    end do
  end do
  
  MWT=XF(1)*PM(1)+XF(2)*PM(2)
  
  return
end subroutine MixRules92

subroutine MAKEK1_92(T,K1,K2,K3)
!--
!-- called by MixRules92
!--
  implicit none
  
  real(dp),intent(in) :: T
  real(dp),intent(out):: K1(2,2),K2(2,2,2),K3(2,2,2)
  
  real(dp):: T2
  
  !---  bis jetzt nur index 1 und 2
  !---- 1=H2O, 2=CO2, 3=CH4
    
  T2=T*T
  
  K1(:,:)=0.0D0
  K2(:,:,:)=0.0D0
  K3(:,:,:)=0.0D0
  
  K1(1,1)=1.0D0
  K1(2,2)=1.0D0
  
  K2(1,1,1)=1.0D0
  K2(2,2,2)=1.0D0
  
  K3(1,1,1)=1.0D0
  K3(2,2,2)=1.0D0
  
  !----- binary H2O-CO2
  
  if (T<=373.15) then
    
    K1(1,2)=0.20611D0+0.0006D0*T
    K1(2,1)=K1(1,2)
    
    K2(1,1,2)=0.8023278D0-0.0022206D0*T+184.76824D0/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    
    K3(1,1,2)=1.80544D0-0.0032605D0*T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
    
  end if
  
  if (T > 373.15D0.and.T<=495.15D0) then
    
    K1(1,2)=-10084.5042 -4.27134485*T +256477.783D0/T &
    &       +0.00166997474D0*T2 +1816.78D0*log(T)
    K1(2,1)=K1(1,2)
    
    K2(1,1,2)=9.000263D0-0.00623494D0*T-2307.7125D0/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    
    K3(1,1,2)=-74.1163D0+0.1800496D0*T-1.40904946D-4*T2+10130.5246D0/T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
    
  end if
  
  if (T > 495.15.and.T<=623.15) then
    
    K1(1,2)=-0.3568D0+7.8888D-4*T+333.399D0/T
    K1(2,1)=K1(1,2)
    
    K2(1,1,2)=-19.97444D0+0.0192515D0*T+5707.4229D0/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    
    K3(1,1,2)=12.1308D0-0.0099489D0*T-3042.09583D0/T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  
  end if
  
  if (T > 623.15) then
    
    K1(1,2)=-4.53122D0+0.0042113D0*T+1619.7D0/T
    K1(2,1)=K1(1,2)
    
    K2(1,1,2)=-163.4855D0+0.190552D0*T-7.228514D-5*T2+46082.885D0/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    
    K3(1,1,2)=1.7137D0-6.7136D-4*T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  
  end if
  
  !!!!      if (T > 623.15) then
  !!      K1(1,2)=0.5D0
  !!      K1(2,1)=K1(1,2)
  !!      K2(1,1,2)=1.0D0
  !!      K2(1,2,1)=K2(1,1,2)
  !!      K2(2,1,1)=K2(1,1,2)
  !!      K2(1,2,2)=K2(1,1,2)
  !!      K2(2,1,2)=K2(1,1,2)
  !!      K2(2,2,1)=K2(1,1,2)
  !!      K3(1,1,2)=1.0D0
  !!      K3(1,2,1)=K3(1,1,2)
  !!      K3(2,1,1)=K3(1,1,2)
  !!      K3(1,2,2)=K3(1,1,2)
  !!      K3(2,1,2)=K3(1,1,2)
  !!      K3(2,2,1)=K3(1,1,2)
  !!!!      end if
  !-----
  !!      call FAKEK(T,1,FF)
  !!      K1(1,2)=FF
  !!      K1(2,1)=K1(1,2)
  !!      call FAKEK(T,2,FF)
  !!      K2(1,1,2)=FF
  !!      K2(1,2,1)=K2(1,1,2)
  !!      K2(2,1,1)=K2(1,1,2)
  !!      K2(1,2,2)=K2(1,1,2)
  !!      K2(2,1,2)=K2(1,1,2)
  !!      K2(2,2,1)=K2(1,1,2)
  !!      call FAKEK(T,3,K3)
  !!      K3(1,1,2)=FF
  !!      K3(1,2,1)=K3(1,1,2)
  !!      K3(2,1,1)=K3(1,1,2)
  !!      K3(1,2,2)=K3(1,1,2)
  !!      K3(2,1,2)=K3(1,1,2)
  !!      K3(2,2,1)=K3(1,1,2)
  
  return
end subroutine MAKEK1_92

subroutine DuanZ( &
& T,P, &
& NP,NSTAB, &
& ZST,VOLST,INTST)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--
  implicit none
  !
  real(dp),intent(in) :: T,P
  integer, intent(out) :: NP,NSTAB 
  real(dp),intent(out):: ZST,VOLST,INTST
  ! logical, intent(out):: Error
  !
  real(dp):: Z,DELTA,START,XX,FF
  real(dp):: ZP(5),INTP(5),VOLP(5)
  integer :: I
  logical :: FOUND
  !
  ! Error= .false.
  !
  NP=0
  START=0.001D0
  DELTA=0.001D0
  
  FOUND=.false.
  
  do I=1,5
  
    call FINDPT( &
    & T,P,START,DELTA, &
    & FOUND,Z,XX)
    
    if (FOUND) then
      NP=NP+1
      VOLP(NP)=XX
      ZP(NP)=Z
      call INTEGRALF(XX,FF)
      INTP(NP)=FF
      START=XX+DELTA
      FOUND=.false.
    else
      exit !goto 888
    end if
    
  end do
  
  ! 888 continue
  
  if (NP==1) NSTAB=1
  
  if (NP==2) then
    if (INTP(2) > INTP(1)) then  ;  NSTAB=2
    else                         ;  NSTAB=1
    end if
  end if
  
  if (NP > 2) then
    write (UNIT=6,FMT='(2(F20.10),'' more phases than expected'')') T,P
    if (INTP(NP) > INTP(1)) then
      NSTAB=NP
    else
      NSTAB=1
    end if
    !     NSTAB=NP
  end if
  
  if (NP==0) then
    write (UNIT=6,FMT= &
    '(2(F20.10),'' less phases than expected'')') T,P
    ! Error= .true.
    ! ErrorMessage= "less phases than expected"
    NSTAB=1
    VOLP(1)=1000.0D0
    INTP(1)=0.0d0
  end if
  
  VOLST=VOLP(NSTAB)
  INTST=INTP(NSTAB)
  ZST=ZP(NSTAB)
  
  return
end subroutine DuanZ

subroutine FINDPT( &
& T,P,START,DELTA, &
& FOUND,Z,XX)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!-- called by DuanZ
!--
  implicit none
  !---------------------------------------------------------------------
  real(dp),intent(in) :: T,P
  real(dp),intent(in) :: START,DELTA
  logical, intent(out):: FOUND
  real(dp),intent(out):: Z,XX
  !---------------------------------------------------------------------
  real(dp):: &
  & V,P1,P2,P3,X1,X2,X3, &
  & VMAXI,PMINI,PDEL,DX,XDEL,VXMINI,AB1,AB2,VVOR
  logical*4 COD1,COD2,COD3
  logical*4 JUM
  !---------------------------------------------------------------------
  JUM=.false.
  FOUND=.false.
  XX=START
  VMAXI=500.0D0
  PMINI=0.00001D0
  VXMINI=1.0D-5
  DX=DELTA
  V=START
  X1=V

  call JUSTF(V,Z,P1)
  call ABLEITF(V,AB1)
  
  COD1=(P > P1)
  
  P2=P1
  AB2=AB1
  COD2=COD1
  
  do1: do while(V<=VMAXI)
    
    DX=V/5.0D0
    VVOR=V
    
    do3: do
    
      V=VVOR+DX      
      X2=V
      
      call JUSTF  (V,     Z,P2)
      call ABLEITF(V,     AB2)
      
      COD2=(P > P2)
      
      if (    COD1 &
      & .and. COD2 &
      & .and.(AB1 > 0.0D0) &
      & .and.(AB2 < 0.0D0) &
      & .and. DX  > 1D-5) then
        DX=DX/2.0D0
        cycle do3 ! goto 3
      end if
      
      exit do3
    
    end do do3
    
    !--- ...or.COD1 ist damit positive Ableitungen uebersprungen werden
    if ((COD1.eqv.COD2).or.COD1) then
      X1=X2
      P1=P2
      AB1=AB2
      COD1=COD2
    else
      JUM=.true.
      exit do1 !goto 10
    end if
    !---

  end do do1
  
  ! 10 continue
  
  if (.not.JUM) return
  
  do11: do
  
    PDEL=abs(P2-P1)
    XDEL=abs(X2-X1)
    
    if ((PDEL < PMINI).and.(XDEL < VXMINI)) exit do11 !goto 15
    
    V=(X1+X2)/2.0D0
    call JUSTF(V,Z,P3)
    X3=V
    COD3=(P > P3)
    !-----
    if (COD3.EQV.COD1) then
      X1=X3
      P1=P3
      COD1=COD3
    else
      X2=X3
      P2=P3
      COD2=COD3
    end if
    
    !goto 11
  end do do11
  
  !-----
  ! 15 continue
  V=(X1+X2)/2.0D0
  call JUSTF(V,Z,P1)
  X1=V
  
  FOUND=.true.
  XX=X1

  return
end subroutine FINDPT

subroutine INTEGRALF(V,PST)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!-- called by DuanZ
!--
  implicit none
  !---------------------------------------------------------------------
  real(dp),intent(in) :: V
  real(dp),intent(out):: PST
  !---------------------------------------------------------------------
  real(dp):: V2,V4,V5
  !---------------------------------------------------------------------
  V2=V**2
  V4=V**4
  V5=V**5
  
  PST=(exp(-GAMMAVC2/V2)) &
  &   *(FVC2/V2 +BETA1*FVC2/GAMMAVC2 +FVC2/GAMMAVC2)/2.0D0 &
  &  +log(V) &
  &  -BVC/V -CVC2/(2.0D0*V2) -DVC4/(4.0D0*V4) -EVC5/(5.0D0*V5)
  PST=PST*RT
  
  return
end subroutine INTEGRALF

subroutine JUSTF(V,Z,P1)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!-- called by FINDPT
!--
  implicit none
  !---------------------------------------------------------------------
  real(dp),intent(in) :: V
  real(dp),intent(out):: Z,P1
  !---------------------------------------------------------------------
  real(dp):: V2,V4,V5
  !--
  V2=V**2
  V4=V2**2
  V5=V4*V
  !
  Z=  1.0D0   &
  & + BVC /V  &
  & + CVC2/V2 &
  & + DVC4/V4 &
  & + EVC5/V5 &
  & +(FVC2/V2)*(BETA1+GAMMAVC2/V2)*exp(-GAMMAVC2/V2)
  !
  P1=Z*RT/V
  !
  return
end subroutine JUSTF

subroutine ABLEITF(V,PST)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!-- called by FINDPT
!--
  implicit none
  !---------------------------------------------------------------------
  real(dp),intent(in) :: V
  real(dp),intent(out):: PST
  !---------------------------------------------------------------------
  real(dp):: V2,V3,V4
  !--
  V2=V*V
  V3=V2*V
  V4=V2*V2
  PST=(exp(-GAMMAVC2/V2)) &
  &  *( -3.0D0 *BETA1 *FVC2 &
  &     +2.0D0 *BETA1 *FVC2*GAMMAVC2/V2 &
  &     -5.0D0 *FVC2 *GAMMAVC2/V2 &
  &     +(2.0D0*FVC2*GAMMAVC2**2)/V4 ) /V4 &
  &  +(-V -2.0D0*BVC -3.0D0*CVC2/V -5.0D0*DVC4/V3 -6.0D0*EVC5/V4)/V3
  PST=PST*RT
  !
  return
end subroutine ABLEITF

subroutine EoS_H2O_CO2_Duan06(   &
& T,P,XF,              &
& V_m3,                &
& ACTIV1,   ACTIV2,    &
& F1PLO,    F2PLO,     &
& F1PHI,    F2PHI,     &
& F1PHI2000,F2PHI2000, &
& F1PLO2000,F2PLO2000)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--
!-- Duan et al, 2006
!-- T= 273 - 533 K
!-- P= 1 - 1000 bar
!--
  implicit none
  
  real(dp),intent(in) :: T,P,XF(2)
  real(dp),intent(out):: V_m3,ACTIV1,ACTIV2
  real(dp),intent(out):: F1PLO,F2PLO,F1PLO2000,F2PLO2000
  real(dp),intent(out):: F1PHI,F2PHI,F1PHI2000,F2PHI2000
  !
  real(dp):: &
  & VAT,AC1,AC2,P2000, &
  & FH2OP,FCO2P, &
  & FH2OX,FCO2X, &
  & F1XHI,    F2XHI, &
  & F1XHI2000,F2XHI2000, &
  & F1XLO2000,F2XLO2000, &
  & F1XLO,F2XLO
      
  P2000=2000.0D0

  if (P<=P2000) then
  
    call DUAN06EQLO( &
    & T,P,XF,        &
    & VAT,AC1,AC2,   &
    & F1PLO,F2PLO,   &
    & F1XLO,F2XLO)
    
    !VOLUM=VAT*100.0D0
    V_m3=VAT*1.0D-3
    
    FH2OP=F1PLO
    FCO2P=F2PLO
    ACTIV1=F1XLO*XF(1)/F1PLO
    ACTIV2=F2XLO*XF(2)/F2PLO
    
  else
  
    call DUAN06EQHI( &
    & T,P,XF,        &
    & VAT,AC1,AC2,   &
    & F1PHI,F2PHI,   &
    & F1XHI,F2XHI)
    
    call DUAN06EQHI( &
    & T,P2000,XF,    &
    & VAT,AC1,AC2,   &
    & F1PHI2000,F2PHI2000, &
    & F1XHI2000,F2XHI2000)
    
    call DUAN06EQLO( &
    & T,P2000,XF,    &
    & VAT,AC1,AC2,   &
    & F1PLO2000,F2PLO2000, &
    & F1XLO2000,F2XLO2000)
    
    FH2OP= F1PHI *F1PLO2000 /F1PHI2000
    FCO2P= F2PHI *F2PLO2000 /F2PHI2000
    !   
    FH2OX= F1XHI *F1XLO2000 /F1XHI2000
    FCO2X= F2XHI *F2XLO2000 /F2XHI2000
    !
    ACTIV1=XF(1)*FH2OX/FH2OP
    ACTIV2=XF(2)*FCO2X/FCO2P
    !
    V_m3=VAT*1.0D-3
    !
    !! ACTIV1=(F1XHI*F1XLO2000/F1XHI2000) /(F1PHI*F1PLO2000/F1PHI2000)
    !! ACTIV1=ACTIV1*XF(1)
    !! ACTIV2=(F2XHI*F2XLO2000/F2XHI2000)/(F2PHI*F2PLO2000/F2PHI2000)
    !! ACTIV2=ACTIV2*XF(2)
    
  end if

end subroutine EoS_H2O_CO2_Duan06

subroutine DUAN06EQHI( &
& T,P,XF, &
& VAT,AH2O,ACO2, &
& FugH2O,FugCO2, &
& ACH,ACC)
! & Error)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--
  implicit none
  
  real(dp),intent(in)   :: T,P,XF(2)
  real(dp),intent(inout):: FugH2O,FugCO2
  real(dp),intent(out)  :: VAT,AH2O,ACO2,ACH,ACC
  ! logical, intent(out):: Error
  !
  real(dp):: &
  & VATP,INTATP,ZATP, &
  & lnPhi1,lnPhi2,V01,V02, &
  & G
  integer:: NP,NSTAB
  ! logical:: Error_DuanZ
  
  ! Error= .false.
  
  if (FugH2O==0.0D0 .and. FugCO2==0.0D0) then
    call EoS_H2O_Duan06_HP(T,P,G, V01,FugH2O)
    call EoS_CO2_Duan06_HP(T,P,G, V02,FugCO2)
  end if
  
  call MixRulesHP(T,P,      XF)
  call DuanZ     (T,P,      NP,NSTAB,ZATP,VATP,INTATP)!,Error_DuanZ)
  !if(Error_DuanZ) then
  !  Error= .true.  
  !  return
  !end if
  call EQUA9(VATP,ZATP,     lnPhi1,lnPhi2)
  
  VAT=VATP
  ACH=exp(lnPhi1)
  ACC=exp(lnPhi2)
  
  AH2O=ACH*XF(1)/FugH2O
  ACO2=ACC*XF(2)/FugCO2
  
  return
end subroutine DUAN06EQHI

subroutine DUAN06EQLO( &
& T,P,XF, &
& VAT,AH2O,ACO2, &
& FugH2O,FugCO2, &
& ACH,ACC)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!-- called by EoS_H2O_CO2_Duan06
!--
  implicit none
  !---------------------------------------------------------------------
  real(dp),intent(in)    :: T,P,XF(2)
  real(dp),intent(inout) :: FugH2O,FugCO2
  real(dp),intent(out)   :: VAT,AH2O,ACO2,ACH,ACC
  !---------------------------------------------------------------------
  real(dp):: VATP,INTATP,ZATP,lnPHI1,lnPHI2,G
  real(dp):: V01,V02
  integer :: NP,NSTAB
  !---------------------------------------------------------------------
  if (FugH2O==0.0D0 .and. FugCO2==0.0D0) then
    call EoS_H2O_Duan06_LP(T,P,G,V01,FugH2O)    
    call EoS_CO2_Duan06_LP(T,P,G,V02,FugCO2)    
  end if
  
  call MixRulesLP(T,P,           XF)
  call DuanZ     (T,P,           NP,NSTAB,ZATP,VATP,INTATP)
  call EQUA9     (VATP,ZATP,     lnPHI1,lnPHI2)
  
  VAT=VATP
  ACH=exp(lnPHI1)
  ACC=exp(lnPHI2)
  
  AH2O=ACH*XF(1)/FugH2O
  ACO2=ACC*XF(2)/FugCO2
  
  return
end subroutine DUAN06EQLO

subroutine EoS_H2O_Duan06_HP( &
& T,P, &
& G,VOLUM,FUGCF)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--
  implicit none
  
  real(dp):: T,P
  real(dp):: G,VOLUM,FUGCF
  !
  real(dp):: TR,TR2,TR3
  real(dp):: VC,VC2,VC4,VC5,EU3,Z,V,INTE
  real(dp):: V2,V4,V5,EXPON,PHI1,H0,S0
  real(dp):: K1,K2,K3,K4,K6
  real(dp):: CPDT,CPTDT,T0,TT,TT0,SQT,SQT0
  integer :: NP,NSTAB
  
  T0=298.15D0
  TT=T*T
  TT0=T0*T0
  SQT=sqrt(T)
  SQT0=sqrt(T0)
  
  R=0.08314467D0
  RT=R*T
  EU3=1.0D0/3.0D0
  ! PR=P/221.19D0
  
  TR=T/647.25D0
  TR2=TR*TR
  TR3=TR2*TR
  
  VC=R*647.25D0/221.19D0
  VC2=VC*VC
  VC4=VC2*VC2
  VC5=VC4*VC
  
  BVC=(4.68071541D-02-2.81275941D-01/TR2-2.43926365D-01/TR3)*VC
  CVC2=(1.10016958D-02-3.86603525D-02/TR2+9.30095461D-02/TR3)*VC2
  DVC4=(-1.15747171D-05+4.19873848D-04/TR2-5.82739501D-04/TR3)*VC4
  EVC5=(1.00936000D-06-1.01713593D-05/TR2+1.63934213D-05/TR3)*VC5
  FVC2=(-4.49505919D-02/TR3)*VC2
  BETA1=-3.15028174D-01
  GAMMAVC2=1.25000000D-02*VC2
  MWT=0.0180154D0
  
  BST1=2.0D0*BVC
  CST1=3.0D0*CVC2
  DST1=5.0D0*DVC4
  EST1=6.0D0*EVC5
  FST1=2.0D0*FVC2
  BESTST1=BETA1
  GAMST1=3.0D0*GAMMAVC2
  
  BST2=0.0D0
  CST2=0.0D0
  DST2=0.0D0
  EST2=0.0D0
  FST2=0.0D0
  BESTST2=0.0D0
  GAMST2=0.0D0
      
  !------------------------------------------------------find the volume
  call DuanZ( &
  & T,P,      &
  & NP,NSTAB, &
  & Z,V,INTE)
    
  !----------------------------------------------------------EQ9 for H2O
  V2=V*V
  V4=V2*V2
  V5=V4*V
  EXPON=exp(-GAMMAVC2/V2)
  PHI1=-log(Z) +BST1/V +CST1/(2.0D0*V2) +DST1/(4.0D0*V4)+EST1/(5.0D0*V5)
  PHI1=PHI1+((FST1*BETA1+BESTST1*FVC2)/(2.0D0*GAMMAVC2))*(1.0D0-EXPON)
  PHI1=PHI1 &
  &   +((FST1*GAMMAVC2+GAMST1*FVC2-FVC2*BETA1*(GAMST1-GAMMAVC2)) &
  &     /(2.0D0*GAMMAVC2**2)) &
  &    *(1.0D0-((GAMMAVC2/V2)+1.0D0)*EXPON) &
  &   -(((GAMST1-GAMMAVC2)*FVC2)/(2.0D0*GAMMAVC2**2)) &
  &    *(2.-(GAMMAVC2**2/V4+2.0D0*GAMMAVC2/V2+2.0D0)*EXPON)
  
  !------------------------------------------------------G for ideal gas
  H0=-241816.00D0
  S0=188.7200D0
  K1=115.45000D0
  K4=-3799.900D0
  K3=-2871300.00D0
  K6=51055.16D0
  K2=-0.002062D0
  CPDT=K1*(T-T0)+K2*(TT-TT0)/2.0D0 &
  &   -K3*(1.0D0/T-1.0D0/T0)+2D0*K4*(SQT-SQT0) &
  &   +K6*log(T/T0)
  CPTDT=K1*log(T/T0)+K2*(T-T0) &
  &    -K3*(1.0D0/(TT)-1.0D0/(TT0))/2.0D0 &
  &    -2D0*K4*(1.0D0/SQT-1.0D0/SQT0) &
  &    -K6*(1.0D0/T-1.0D0/T0)
  G=H0+CPDT-T*(S0+CPTDT)
      
  !-------------------------------------------------correct for fugacity
  FUGCF=exp(PHI1)
  G=G + RT*100.0D0 *log(FUGCF*P)
  VOLUM=V*100.0D0
  
  return
end subroutine EoS_H2O_Duan06_HP

subroutine EoS_H2O_Duan06_LP(T,P,G,VOLUM,FUGCF)
  implicit none

  real(dp),intent(in) :: T,P
  real(dp),intent(out):: G,VOLUM,FUGCF
  !---------------------------------------------------------------------
  real(dp):: &
  & TR,TR2,TR3,VC,VC2,VC4,VC5,EU3,Z,V,INTE, &
  & V2,V4,V5,EXPON,PHI1,&
  & H0,S0,K1,K2,K3,K4,K6, &
  & CPDT,CPTDT, &
  & T0,TT,TT0,SQT,SQT0
  integer :: NP,NSTAB
  !---------------------------------------------------------------------
  EU3=1.0D0/3.0D0
  T0=298.15D0
  TT=T*T
  TT0=T0*T0
  SQT=sqrt(T)
  SQT0=sqrt(T0)
  
  R=0.08314467D0
  RT=R*T
  
  ! PR=P/221.19D0
  TR=T/647.25D0
  TR2=TR*TR
  TR3=TR2*TR
  
  VC=R*647.25D0/221.19D0
  VC2=VC*VC
  VC4=VC2*VC2
  VC5=VC4*VC
  
  BVC= ( 4.38269941D-02 -1.68244362D-01/TR2 -2.36923373D-01/TR3)*VC
  CVC2=( 1.13027462D-02 -7.67764181D-02/TR2 +9.71820593D-02/TR3)*VC2
  DVC4=( 6.62674916D-05 +1.06637349D-03/TR2 -1.23265258D-03/TR3)*VC4
  EVC5=(-8.93953948D-06 -3.88124606D-05/TR2 +5.61510206D-05/TR3)*VC5
  FVC2=(7.51274488D-03/TR3)*VC2
  BETA1=2.51598931D+00
  GAMMAVC2=3.94000000D-02*VC2
  MWT=0.0180154D0
  
  BST1=2.0D0*BVC
  CST1=3.0D0*CVC2
  DST1=5.0D0*DVC4
  EST1=6.0D0*EVC5
  FST1=2.0D0*FVC2
  BESTST1=BETA1
  GAMST1=3.0D0*GAMMAVC2
  BST2=0.0D0
  CST2=0.0D0
  DST2=0.0D0
  EST2=0.0D0
  FST2=0.0D0
  BESTST2=0.0D0
  GAMST2=0.0D0

  !----------------------------------------------------find the volume--
  call DuanZ( &
  & T,P,      &
  & NP,NSTAB, &
  & Z,V,INTE)
  
  !--------------------------------------------------------EQ9 for H2O--
  V2=V*V
  V4=V2*V2
  V5=V4*V
  
  EXPON=exp(-GAMMAVC2/V2)
  
  PHI1=-log(Z)+BST1/V+CST1/(2.0D0*V2)+DST1/(4.0D0*V4)+EST1/(5.0D0*V5)
  PHI1=PHI1 &
  &   +((FST1*BETA1+BESTST1*FVC2)/(2.0D0*GAMMAVC2))*(1.0D0-EXPON)
  PHI1=PHI1 &
  &   +( (FST1*GAMMAVC2+GAMST1*FVC2-FVC2*BETA1*(GAMST1-GAMMAVC2)) &
  &       /(2.0D0*GAMMAVC2**2) )*(1.0D0-((GAMMAVC2/V2)+1.0D0)*EXPON) &
  &   -( ((GAMST1-GAMMAVC2)*FVC2) / (2.0D0*GAMMAVC2**2) ) &
  &    *(2.-(GAMMAVC2**2/V4+2.0D0*GAMMAVC2/V2+2.0D0)*EXPON)
  
  !----------------------------------------------------G for ideal gas--
  H0=-241816.00D0
  S0=188.7200D0
  K1=115.45000D0
  K4=-3799.900D0
  K3=-2871300.00D0
  K6=51055.16D0
  K2=-0.002062D0
  
  CPDT=K1*(T-T0) &
  &   +K2*(TT-TT0)/2.0D0 &
  &   -K3*(1.0D0/T-1.0D0/T0) &
  &   +2D0*K4*(SQT-SQT0) &
  &   +K6*log(T/T0)
  
  CPTDT=K1*log(T/T0) &
  &    +K2*(T-T0) &
  &    -K3*(1.0D0/(TT)-1.0D0/(TT0))/2.0D0 &
  &    -2D0*K4*(1.0D0/SQT-1.0D0/SQT0) &
  &    -K6*(1.0D0/T-1.0D0/T0)
  
  G=H0+CPDT-T*(S0+CPTDT)

  !-----------------------------------------------correct for fugacity--
  FUGCF=exp(PHI1)
  G=G+RT*100.0D0*log(FUGCF*P)
  VOLUM=V*100.0D0

  return
end subroutine EoS_H2O_Duan06_LP

subroutine EoS_CO2_Duan06_HP(T,P,G,VOLUM,FUGCF)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--
  implicit none
  !---------------------------------------------------------------------
  real(dp),intent(in) :: T,P 
  real(dp),intent(out):: G,VOLUM,FUGCF
  real(dp):: &
  TR,TR2,TR3,VC,VC2,VC4,VC5,EU3, &
  & Z,V,INTE, &
  V2,V4,V5,EXPON,PHI2,H0,S0,K1,K2,K3,K4,K6, &
  CPDT,CPTDT,T0,TT,TT0,SQT,SQT0
  integer :: NP,NSTAB
  !---------------------------------------------------------------------
  T0=298.15D0
  TT=T*T
  TT0=T0*T0
  SQT=sqrt(T)
  SQT0=sqrt(T0)
  !---
  R=0.08314467D0
  RT=R*T
  EU3=1.0D0/3.0D0
  !      PR=P/221.19D0
  TR=T/304.1282D0
  TR2=TR*TR
  TR3=TR2*TR
  
  VC=R*304.1282D0/73.773D0
  VC2=VC*VC
  VC4=VC2*VC2
  VC5=VC4*VC
  
  BVC=(5.72573440D-03+7.94836769D+00/TR2-3.84236281D+01/TR3)*VC
  CVC2=(3.71600369D-02-1.92888994D+00/TR2+6.64254770D+00/TR3)*VC2
  DVC4=(-7.02203950D-06+1.77093234D-02/TR2-4.81892026D-02/TR3)*VC4
  EVC5=(3.88344869D-06-5.54833167D-04/TR2+1.70489748D-03/TR3)*VC5
  FVC2=(-4.13039220D-01/TR3)*VC2
  
  BETA1=-8.47988634D+00
  GAMMAVC2=2.80000000D-02*VC2
  MWT=0.0440098D0
  BST2=2.0D0*BVC
  CST2=3.0D0*CVC2
  DST2=5.0D0*DVC4
  EST2=6.0D0*EVC5
  FST2=2.0D0*FVC2
  BESTST2=BETA1
  GAMST2=3.0D0*GAMMAVC2
  BST1=0.0D0
  CST1=0.0D0
  DST1=0.0D0
  EST1=0.0D0
  FST1=0.0D0
  BESTST1=0.0D0
  GAMST1=0.0D0
  
  !---- find the volume
  call DuanZ( &
  & T,P,      &
  & NP,NSTAB, &
  & Z,V,INTE)
  
  !---- EQ9 for CO2
  V2=V*V
  V4=V2*V2
  V5=V4*V
  EXPON=exp(-GAMMAVC2/V2)
  PHI2=-log(Z)+BST2/V+CST2/(2.0D0*V2)+DST2/(4.0D0*V4)+EST2/(5.0D0*V5)
  PHI2=PHI2 &
  &   +((FST2*BETA1+BESTST2*FVC2)/(2.0D0*GAMMAVC2))*(1.0D0-EXPON)
  PHI2=PHI2 &
  &   +((FST2*GAMMAVC2+GAMST2*FVC2-FVC2*BETA1*(GAMST2-GAMMAVC2))/(2.0D0*GAMMAVC2**2)) &
  &   *(1.0D0-((GAMMAVC2/V2)+1.0D0)*EXPON) &
  &   -(((GAMST2-GAMMAVC2)*FVC2)/(2.0D0*GAMMAVC2**2)) &
  &    *(2.-(GAMMAVC2**2/V4+2.0D0*GAMMAVC2/V2+2.0D0)*EXPON)
  
  !---- G for ideal gas
  H0=-393510.010
  S0=213.6770D0
  K1=93.0000D0
  K4=-1340.900D0
  K3=123800.000D0
  K6=6336.20D0
  K2=-0.002876D0
  
  CPDT=K1*(T-T0) &
  &   +K2*(TT-TT0)/2.0D0 &
  &   -K3*(1.0D0/T-1.0D0/T0) &
  &   +2D0*K4*(SQT-SQT0) &
  &   +K6*log(T/T0)
  
  CPTDT=K1*log(T/T0) &
  &    +K2*(T-T0) &
  &    -K3*(1.0D0/(TT)-1.0D0/(TT0))/2.0D0 &
  &    -2D0*K4*(1.0D0/SQT-1.0D0/SQT0) &
  &    -K6*(1.0D0/T-1.0D0/T0)
  
  G=H0+CPDT-T*(S0+CPTDT)
  
  !---- correct for fugacity
  FUGCF=exp(PHI2)
  G=G+RT*100.0D0*log(FUGCF*P)
  VOLUM=V*100.0D0
  
  return
end subroutine EoS_CO2_Duan06_HP

subroutine EoS_CO2_Duan06_LP(T,P,G,VOLUM,FUGCF)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--
  implicit none
  !---------------------------------------------------------------------
  real(dp):: &
  TR,TR2,TR3,VC,VC2,VC4,VC5,EU3,T,P,Z,V,INTE, &
  V2,V4,V5,EXPON,PHI2,G,H0,S0,K1,K2,K3,K4,K6, &
  CPDT,CPTDT,T0,TT,TT0,SQT,SQT0,FUGCF,VOLUM
  integer :: NP,NSTAB
  !---------------------------------------------------------------------
  T0=298.15D0
  TT=T*T
  TT0=T0*T0
  SQT=sqrt(T)
  SQT0=sqrt(T0)
  
  R=0.08314467D0
  RT=R*T
  EU3=1.0D0/3.0D0
  !      PR=P/221.19D0
  TR=T/304.1282D0
  TR2=TR*TR
  TR3=TR2*TR
  VC=R*304.1282D0/73.773D0
  VC2=VC*VC
  VC4=VC2*VC2
  VC5=VC4*VC
  BVC=(1.14400435D-01-9.38526684D-01/TR2+7.21857006D-01/TR3)*VC
  CVC2=(8.81072902D-03+6.36473911D-02/TR2-7.70822213D-02/TR3)*VC2
  DVC4=(9.01506064D-04-6.81834166D-03/TR2+7.32364258D-03/TR3)*VC4
  EVC5=(-1.10288237D-04+1.26524193D-03/TR2-1.49730823D-03/TR3)*VC5
  FVC2=(7.81940730D-03/TR3)*VC2
  BETA1=-4.22918013D+00
  GAMMAVC2=1.58500000D-01*VC2
  MWT=0.0440098D0
  BST2=2.0D0*BVC
  CST2=3.0D0*CVC2
  DST2=5.0D0*DVC4
  EST2=6.0D0*EVC5
  FST2=2.0D0*FVC2
  BESTST2=BETA1
  GAMST2=3.0D0*GAMMAVC2
  BST1=0.0D0
  CST1=0.0D0
  DST1=0.0D0
  EST1=0.0D0
  FST1=0.0D0
  BESTST1=0.0D0
  GAMST1=0.0D0
  
  !---- find the volume
  call DuanZ( &
  & T,P,      &
  & NP,NSTAB, &
  & Z,V,INTE)
  
  !---- EQ9 for CO2
  V2=V*V
  V4=V2*V2
  V5=V4*V
  EXPON=exp(-GAMMAVC2/V2)
  PHI2=-log(Z)+BST2/V+CST2/(2.0D0*V2)+ &
  DST2/(4.0D0*V4)+EST2/(5.0D0*V5)
  PHI2=PHI2+((FST2*BETA1+BESTST2*FVC2)/ &
  (2.0D0*GAMMAVC2))*(1.0D0-EXPON)
  PHI2=PHI2+((FST2*GAMMAVC2+GAMST2*FVC2- &
  FVC2*BETA1*(GAMST2-GAMMAVC2))/(2.0D0*GAMMAVC2**2))* &
  (1.0D0-((GAMMAVC2/V2)+1.0D0)*EXPON)-(((GAMST2-GAMMAVC2)* &
  FVC2)/(2.0D0*GAMMAVC2**2))*(2.-(GAMMAVC2**2/V4+ &
  2.0D0*GAMMAVC2/V2+2.0D0)*EXPON)
  
  !---- G for ideal gas
  H0=-393510.010
  S0=213.6770D0
  K1=93.0000D0
  K4=-1340.900D0
  K3=123800.000D0
  K6=6336.20D0
  K2=-0.002876D0
  CPDT=K1*(T-T0)+K2*(TT-TT0)/2.0D0 &
  &   -K3*(1.0D0/T-1.0D0/T0) &
  &   +2D0*K4*(SQT-SQT0) &
  &   +K6*log(T/T0)
  CPTDT=K1*log(T/T0) &
  &    +K2*(T-T0) &
  &    -K3*(1.0D0/(TT)-1.0D0/(TT0))/2.0D0 &
  &    -2D0*K4*(1.0D0/SQT-1.0D0/SQT0) &
  &    -K6*(1.0D0/T-1.0D0/T0)
  G=H0+CPDT-T*(S0+CPTDT)
  
  !---- correct for fugacity
  FUGCF=exp(PHI2)
  G=G+RT*100.0D0*log(FUGCF*P)
  VOLUM=V*100.0D0
  
  return
end subroutine EoS_CO2_Duan06_LP

subroutine MixRulesHP(T,P,XF)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!-- called by DUAN06EQHI
!--
  implicit none
  
  !---- 1=H2O, 2=CO2, 3=CH4
  real(dp):: &
  A1(2),A2(2),A3(2),A4(2),A5(2),A6(2), &
  A7(2),A8(2),A9(2),A10(2),A11(2),A12(2), &
  ALPHA(2),BETA(2),GAMMA(2), &
  TC(2),PC(2), &
  BB(2),CC(2),DD(2),EE(2),FF(2), &
  BS(2),CS(2),DS(2),ES(2),FS(2), &
  PM(2), &
  VC(2),VC2(2),VC4(2),VC5(2), &
  TR(2),TR2(2),TR3(2),PR(2),P,T,XF(2),EU3
  integer:: I,J,K,II,L,M,N
  real(dp):: &
  K1(2,2),K2(2,2,2),K3(2,2,2), &
  BIJ(2,2),VCIJ(2,2),CIJK(2,2,2),VCIJK(2,2,2), &
  DIJKLM(2,2,2,2,2),VCIJKLM(2,2,2,2,2), &
  EIJKLMN(2,2,2,2,2,2),VCIJKLMN(2,2,2,2,2,2), &
  FIJ(2,2),GAMIJK(2,2,2)
  
  data A1  / 4.68071541D-02, 5.72573440D-03/
  data A2  /-2.81275941D-01, 7.94836769D+00/
  data A3  /-2.43926365D-01,-3.84236281D+01/
  data A4  / 1.10016958D-02, 3.71600369D-02/
  data A5  /-3.86603525D-02,-1.92888994D+00/
  data A6  / 9.30095461D-02, 6.64254770D+00/
  data A7  /-1.15747171D-05,-7.02203950D-06/
  data A8  / 4.19873848D-04, 1.77093234D-02/
  data A9  /-5.82739501D-04,-4.81892026D-02/
  data A10 / 1.00936000D-06, 3.88344869D-06/
  data A11 /-1.01713593D-05,-5.54833167D-04/
  data A12 / 1.63934213D-05, 1.70489748D-03/
  
  data ALPHA /-4.49505919D-02,-4.13039220D-01/
  data BETA  /-3.15028174D-01,-8.47988634D+00/
  data GAMMA / 1.25000000D-02, 2.80000000D-02/
  
  data TC /647.25D0,304.1282D0/
  data PC /221.19D0,73.773D0/
  data PM /0.0180154D0,0.0440098D0/
  
  R=0.08314467D0
  RT=R*T
  EU3=1.0D0/3.0D0
  II=1
      
  do I=1,2

    PR(I)=P/PC(I)
    TR(I)=T/TC(I)
    TR2(I)=TR(I)*TR(I)
    TR3(I)=TR2(I)*TR(I)
    VC(I)=R*TC(I)/PC(I)
    VC2(I)=VC(I)*VC(I)
    VC4(I)=VC2(I)*VC2(I)
    VC5(I)=VC4(I)*VC(I)

    BB(I)=A1(I)+A2(I)/TR2(I)+A3(I)/TR3(I)
    CC(I)=A4(I)+A5(I)/TR2(I)+A6(I)/TR3(I)
    DD(I)=A7(I)+A8(I)/TR2(I)+A9(I)/TR3(I)
    EE(I)=A10(I)+A11(I)/TR2(I)+A12(I)/TR3(I)
    FF(I)=ALPHA(I)/TR3(I)
       
    if (BB(I)<0.0D0) then
      BB(I)=-BB(I)
      BS(I)=-1.0D0
    else
      BS(I)=1.0D0
    end if
    
    if (CC(I)<0.0D0) then
      CC(I)=-CC(I)
      CS(I)=-1.0D0
    else
      CS(I)=1.0D0
    end if
    
    if (DD(I)<0.0D0) then
      DD(I)=-DD(I)
      DS(I)=-1.0D0
    else
      DS(I)=1.0D0
    end if
    
    if (EE(I)<0.0D0) then
      EE(I)=-EE(I)
      ES(I)=-1.0D0
    else
      ES(I)=1.0D0
    end if
    
    if (FF(I)<0.0D0) then
      FF(I)=-FF(I)
      FS(I)=-1.0D0
    else
      FS(I)=1.0D0
    end if
    
  end do

  call MAKEK1_HI(T,K1,K2,K3)

  do I=1,2
    do J=1,2
      BIJ(I,J)=(((BS(I)*BB(I)**EU3+BS(J)*BB(J)**EU3) &
           /2.0D0)**3)*K1(I,J)
      VCIJ(I,J)=((VC(I)**EU3+VC(J)**EU3)/2.0D0)**3
    end do
  end do
  
  BVC=0.0D0
  do I=1,2
    do J=1,2
      BVC=BVC+XF(I)*XF(J)*BIJ(I,J)*VCIJ(I,J)
    end do
  end do
  
  BST1=0.0D0
  BST2=0.0D0
  do J=1,2
    BST1=BST1+2.0D0*XF(J)*BIJ(1,J)*VCIJ(1,J)
    BST2=BST2+2.0D0*XF(J)*BIJ(J,2)*VCIJ(J,2)
  end do

  do I=1,2
    do J=1,2
      do K=1,2
        CIJK(I,J,K)=(((CS(I)*CC(I)**EU3+CS(J)*CC(J)**EU3+CC(K)**EU3) &
              /3.0D0)**3)*K2(I,J,K)
        VCIJK(I,J,K)=((VC(I)**EU3+VC(J)**EU3+VC(K)**EU3)/3.0D0)**3
      end do
    end do
  end do
  
  CVC2=0.0D0
  do I=1,2
    do J=1,2
      do K=1,2
        CVC2=CVC2+XF(I)*XF(J)*XF(K)*CIJK(I,J,K)*VCIJK(I,J,K)**2
      end do
    end do
  end do
  
  CST1=0.0D0
  CST2=0.0D0
  do J=1,2
    do K=1,2
      CST1=CST1+3.0D0*XF(J)*XF(K)*CIJK(1,J,K)*VCIJK(1,J,K)**2
      CST2=CST2+3.0D0*XF(J)*XF(K)*CIJK(J,2,K)*VCIJK(J,2,K)**2
    end do
  end do
      
  do I=1,2
   do J=1,2
    do K=1,2
     do L=1,2
      do M=1,2
        DIJKLM(I,J,K,L,M)= &
          ((DS(I)*DD(I)**EU3 &
          + DS(J)*DD(J)**EU3 &
          + DS(K)*DD(K)**EU3 &
          + DS(L)*DD(L)**EU3 &
          + DS(M)*DD(M)**EU3) /5.0D0)**3
        VCIJKLM(I,J,K,L,M)= &
          ((VC(I)**EU3      &
          + VC(J)**EU3      &
          + VC(K)**EU3      &
          + VC(L)**EU3      &
          + VC(M)**EU3)/5.0D0)**3
      end do
     end do
    end do
   end do
  end do
  
  DVC4=0.0D0
  do I=1,2
   do J=1,2
    do K=1,2
     do L=1,2
      do M=1,2
        DVC4= DVC4 &
        &   + XF(I)*XF(J)*XF(K)*XF(L)*XF(M) &
        &     *DIJKLM(I,J,K,L,M)*VCIJKLM(I,J,K,L,M)**4
      end do
     end do
    end do
   end do
  end do
  DST1=0.0D0
  DST2=0.0D0
  do J=1,2
   do K=1,2
    do L=1,2
     do M=1,2
  DST1=DST1+5.0D0*XF(J)*XF(K)*XF(L)*XF(M)*DIJKLM(1,J,K,L,M)* &
      VCIJKLM(1,J,K,L,M)**4
  DST2=DST2+5.0D0*XF(J)*XF(K)*XF(L)*XF(M)*DIJKLM(J,2,K,L,M)* &
      VCIJKLM(J,2,K,L,M)**4
     end do
    end do
   end do
  end do
      
  do I=1,2
   do J=1,2
    do K=1,2
     do L=1,2
      do M=1,2
       do N=1,2
  EIJKLMN(I,J,K,L,M,N)=((ES(I)*EE(I)**EU3+ES(J)*EE(J)**EU3+ &
                         ES(K)*EE(K)**EU3+ES(L)*EE(L)**EU3+ &
                         ES(M)*EE(M)**EU3+ &
                         ES(N)*EE(N)**EU3)/6.0D0)**3
  VCIJKLMN(I,J,K,L,M,N)=((VC(I)**EU3+VC(J)**EU3+ &
                          VC(K)**EU3+VC(L)**EU3+ &
                          VC(M)**EU3+VC(N)**EU3)/6.0D0)**3
       end do
      end do
     end do
    end do
   end do
  end do
  EVC5=0.0D0
  do I=1,2
   do J=1,2
    do K=1,2
     do L=1,2
      do M=1,2
       do N=1,2
        EVC5=EVC5+XF(I)*XF(J)*XF(K)*XF(L)*XF(M)*XF(N)* &
        EIJKLMN(I,J,K,L,M,N)*VCIJKLMN(I,J,K,L,M,N)**5
       end do
      end do
     end do
    end do
   end do
  end do
  
  EST1=0.0D0
  EST2=0.0D0
  do J=1,2
   do K=1,2
    do L=1,2
     do M=1,2
      do N=1,2
       EST1=EST1+6.0D0*XF(J)*XF(K)*XF(L)*XF(M)*XF(N)* &
       EIJKLMN(1,J,K,L,M,N)*VCIJKLMN(1,J,K,L,M,N)**5
       EST2=EST2+6.0D0*XF(J)*XF(K)*XF(L)*XF(M)*XF(N)* &
       EIJKLMN(J,2,K,L,M,N)*VCIJKLMN(J,2,K,L,M,N)**5
      end do
     end do
    end do
   end do
  end do
      
  do I=1,2
    do J=1,2
      FIJ(I,J)=((FS(I)*FF(I)**EU3+FS(J)*FF(J)**EU3)/2.0D0)**3
    end do
  end do
  
  FVC2=0.0D0
  do I=1,2
    do J=1,2
      FVC2=FVC2+XF(I)*XF(J)*FIJ(I,J)*VCIJ(I,J)**2
    end do
  end do
  
  FST1=0.0D0
  FST2=0.0D0
  do J=1,2
    FST1=FST1+2.0D0*XF(J)*FIJ(1,J)*VCIJ(1,J)**2
    FST2=FST2+2.0D0*XF(J)*FIJ(J,2)*VCIJ(J,2)**2
  end do
      
  BETA1=0.0D0
  do I=1,2
    BETA1=BETA1+XF(I)*BETA(I)
  end do
  BESTST1=BETA(1)
  BESTST2=BETA(2)
      
  do I=1,2
   do J=1,2
    do K=1,2
      GAMIJK(I,J,K)=(((GAMMA(I)**EU3+GAMMA(J)**EU3+ &
      GAMMA(K)**EU3)/3.0D0)**3)*K3(I,J,K)
    end do
   end do
  end do
  
  GAMMAVC2=0.0D0
  do I=1,2
   do J=1,2
    do K=1,2
      GAMMAVC2=GAMMAVC2+XF(I)*XF(J)*XF(K)*GAMIJK(I,J,K)* &
      VCIJK(I,J,K)**2
    end do
   end do
  end do
  
  GAMST1=0.0D0
  GAMST2=0.0D0
  do J=1,2
    do K=1,2
      GAMST1=GAMST1+3.0D0*XF(J)*XF(K)*GAMIJK(1,J,K)* &
      VCIJK(1,J,K)**2
      GAMST2=GAMST2+3.0D0*XF(J)*XF(K)*GAMIJK(J,2,K)* &
      VCIJK(J,2,K)**2
    end do
  end do
  
  MWT=XF(1)*PM(1)+XF(2)*PM(2)
  
  return
end subroutine MixRulesHP

subroutine MAKEK1_HI(T,K1,K2,K3)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!-- called by MixRulesHP
!--
  implicit none
  
  real(dp):: T,T2
  real(dp):: K1(2,2),K2(2,2,2),K3(2,2,2)
  
  !--- bis jetzt nur index 1 und 2
  !---- 1=H2O, 2=CO2, 3=CH4
  K1(:,:)=0.0D0
  K2(:,:,:)=0.0D0
  K3(:,:,:)=0.0D0
     
  T2=T*T
  
  K1(1,1)=1.0D0
  K1(2,2)=1.0D0
  
  K2(1,1,1)=1.0D0
  K2(2,2,2)=1.0D0
  K3(1,1,1)=1.0D0
  K3(2,2,2)=1.0D0

  !----- binary H2O-CO2
  if (T<=373.15) then
    K1(1,2)=0.20611+0.0006*T
    K1(2,1)=K1(1,2)
    
    K2(1,1,2)=0.8023278-0.0022206*T+184.76824/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    
    K3(1,1,2)=1.80544-0.0032605*T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  end if

  if (T>373.15.and.T<=495.15) then
    K1(1,2)=-10084.5042-4.27134485*T &
    &       +256477.783/T+0.00166997474*T2+1816.78*log(T)
    K1(2,1)=K1(1,2)
    
    K2(1,1,2)=9.000263-0.00623494*T-2307.7125/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    
    K3(1,1,2)=-74.1163+0.1800496*T-1.40904946E-04*T2+10130.5246/T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  end if
  !--
  if (T>495.15.and.T<=623.15) then
    K1(1,2)=-0.3568+7.8888E-04*T+333.399/T
    K1(2,1)=K1(1,2)
    
    K2(1,1,2)=-19.97444+0.0192515*T+5707.4229/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    
    K3(1,1,2)=12.1308-0.0099489*T-3042.09583/T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  end if
  !--
  if (T>623.15 .and. T<=672.15) then
    K1(1,2)=-4.53122+0.0042113*T+1619.7/T
    K1(2,1)=K1(1,2)
    
    K2(1,1,2)=-163.4855+0.190552*T-7.228514E-05*T2+46082.885/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    
    K3(1,1,2)=1.7137-6.7136E-04*T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  end if
  !--
  !!!!      if (T>672.15) then
    K1(1,2)=9.034D0-7.9212D-3*T+2.3285D-6*T2-2.4221D+03/T
    K1(2,1)=K1(1,2)
    
    K2(1,1,2)=-1.068D0+1.8756D-3*T-4.9371D-7*T2+6.6180D2/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    
    K3(1,1,2)=1.0D0
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  !!!!      end if
  
  return
end subroutine MAKEK1_HI

subroutine MixRulesLP(T,P,XF)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!-- called by DUAN06EQLO
!--
  implicit none
  
  !! real*8 BVC,CVC2,DVC4,EVC5,FVC2,BETA1,GAMMAVC2,R,RT,MWT, &
  !! BST1,CST1,DST1,EST1,FST1,BESTST1,GAMST1, &
  !! BST2,CST2,DST2,EST2,FST2,BESTST2,GAMST2
  !! common /DTHINGS/ BVC,CVC2,DVC4,EVC5,FVC2,BETA1,GAMMAVC2,R,RT, &
  !! MWT, &
  !! BST1,CST1,DST1,EST1,FST1,BESTST1,GAMST1, &
  !! BST2,CST2,DST2,EST2,FST2,BESTST2,GAMST2
      
  !---- 1=H2O, 2=CO2, 3=CH4
  real(dp):: &
  A1(2),A2(2),A3(2),A4(2),A5(2),A6(2), &
  A7(2),A8(2),A9(2),A10(2),A11(2),A12(2), &
  ALPHA(2),BETA(2),GAMMA(2), &
  TC(2),PC(2), &
  BB(2),CC(2),DD(2),EE(2),FF(2), &
  BS(2),CS(2),DS(2),ES(2),FS(2), &
  PM(2), &
  VC(2),VC2(2),VC4(2),VC5(2), &
  TR(2),TR2(2),TR3(2),PR(2),P,T,XF(2),EU3
  integer:: I,J,K,II,L,M,N
  real(dp):: &
  K1(2,2),K2(2,2,2),K3(2,2,2), &
  BIJ(2,2),VCIJ(2,2), &
  CIJK(2,2,2),VCIJK(2,2,2), &
  DIJKLM(2,2,2,2,2),VCIJKLM(2,2,2,2,2), &
  EIJKLMN(2,2,2,2,2,2),VCIJKLMN(2,2,2,2,2,2), &
  FIJ(2,2),GAMIJK(2,2,2)
  
  data A1  / 4.38269941D-02, 1.14400435D-01/
  data A2  /-1.68244362D-01,-9.38526684D-01/
  data A3  /-2.36923373D-01, 7.21857006D-01/
  data A4  / 1.13027462D-02, 8.81072902D-03/
  data A5  /-7.67764181D-02, 6.36473911D-02/
  data A6  / 9.71820593D-02,-7.70822213D-02/
  data A7  / 6.62674916D-05, 9.01506064D-04/
  data A8  / 1.06637349D-03,-6.81834166D-03/
  data A9  /-1.23265258D-03, 7.32364258D-03/
  data A10 /-8.93953948D-06,-1.10288237D-04/
  data A11 /-3.88124606D-05, 1.26524193D-03/
  data A12 / 5.61510206D-05,-1.49730823D-03/
  
  data ALPHA / 7.51274488D-03, 7.81940730D-03/
  data BETA  / 2.51598931D+00,-4.22918013D+00/
  data GAMMA / 3.94000000D-02, 1.58500000D-01/
  
  data TC /647.25D0,304.1282D0/
  data PC /221.19D0,73.773D0/
  data PM /0.0180154D0,0.0440098D0/
  
  R=0.08314467D0
  RT=R*T
  EU3=1.0D0/3.0D0
  II=1
      
  do I=1,2
    PR(I)=P/PC(I)
    TR(I)=T/TC(I)
    TR2(I)=TR(I)*TR(I)
    TR3(I)=TR2(I)*TR(I)
    VC(I)=R*TC(I)/PC(I)
    VC2(I)=VC(I)*VC(I)
    VC4(I)=VC2(I)*VC2(I)
    VC5(I)=VC4(I)*VC(I)
    !
    BB(I)=A1(I)+A2(I)/TR2(I)+A3(I)/TR3(I)
    CC(I)=A4(I)+A5(I)/TR2(I)+A6(I)/TR3(I)
    DD(I)=A7(I)+A8(I)/TR2(I)+A9(I)/TR3(I)
    EE(I)=A10(I)+A11(I)/TR2(I)+A12(I)/TR3(I)
    FF(I)=ALPHA(I)/TR3(I)
    
    if (BB(I).LT.0.0D0) then
      BB(I)=-BB(I)
      BS(I)=-1.0D0
    else
      BS(I)=1.0D0
    end if
    if (CC(I).LT.0.0D0) then
      CC(I)=-CC(I)
      CS(I)=-1.0D0
    else
      CS(I)=1.0D0
    end if
    if (DD(I).LT.0.0D0) then
      DD(I)=-DD(I)
      DS(I)=-1.0D0
    else
      DS(I)=1.0D0
    end if
    if (EE(I).LT.0.0D0) then
      EE(I)=-EE(I)
      ES(I)=-1.0D0
    else
      ES(I)=1.0D0
    end if
    if (FF(I).LT.0.0D0) then
     FF(I)=-FF(I)
     FS(I)=-1.0D0
    else
     FS(I)=1.0D0
    end if
    
  end do
      
  call MAKEK1_LO(T,K1,K2,K3)
      
  do I=1,2
    do J=1,2
      BIJ(I,J)=(((BS(I)*BB(I)**EU3+BS(J)*BB(J)**EU3)/2.0D0)**3)*K1(I,J)
      VCIJ(I,J)=((VC(I)**EU3+VC(J)**EU3)/2.0D0)**3
    end do
  end do
  BVC=0.0D0
  do I=1,2
   do J=1,2
  BVC=BVC+XF(I)*XF(J)*BIJ(I,J)*VCIJ(I,J)
   end do
  end do
  BST1=0.0D0
  BST2=0.0D0
  do J=1,2
   BST1=BST1+2.0D0*XF(J)*BIJ(1,J)*VCIJ(1,J)
   BST2=BST2+2.0D0*XF(J)*BIJ(J,2)*VCIJ(J,2)
  end do
      
  do I=1,2
   do J=1,2
    do K=1,2
  CIJK(I,J,K)=(((CS(I)*CC(I)**EU3+CS(J)*CC(J)**EU3+CC(K)**EU3) &
              /3.0D0)**3)*K2(I,J,K)
  VCIJK(I,J,K)=((VC(I)**EU3+VC(J)**EU3+VC(K)**EU3)/3.0D0)**3
    end do
   end do
  end do
  CVC2=0.0D0
  do I=1,2
   do J=1,2
    do K=1,2
  CVC2=CVC2+XF(I)*XF(J)*XF(K)*CIJK(I,J,K)*VCIJK(I,J,K)**2
    end do
   end do
  end do
  CST1=0.0D0
  CST2=0.0D0
  do J=1,2
   do K=1,2
    CST1=CST1+3.0D0*XF(J)*XF(K)*CIJK(1,J,K)*VCIJK(1,J,K)**2
    CST2=CST2+3.0D0*XF(J)*XF(K)*CIJK(J,2,K)*VCIJK(J,2,K)**2
   end do
  end do
      
  do I=1,2
   do J=1,2
    do K=1,2
     do L=1,2
      do M=1,2
        DIJKLM(I,J,K,L,M)= &
        & ((DS(I)*DD(I)**EU3+DS(J)*DD(J)**EU3+ &
            DS(K)*DD(K)**EU3+DS(L)*DD(L)**EU3+ &
            DS(M)*DD(M)**EU3)/5.0D0)**3
        VCIJKLM(I,J,K,L,M)= &
        & ((VC(I)**EU3+VC(J)**EU3+VC(K)**EU3+ &
            VC(L)**EU3+VC(M)**EU3)/5.0D0)**3
      end do
     end do
    end do
   end do
  end do
  
  DVC4=0.0D0
  do I=1,2
   do J=1,2
    do K=1,2
     do L=1,2
      do M=1,2
       DVC4=DVC4 &
       &   + XF(I)*XF(J)*XF(K)*XF(L)*XF(M) &
       &    *DIJKLM(I,J,K,L,M)*VCIJKLM(I,J,K,L,M)**4
      end do
     end do
    end do
   end do
  end do
  
  DST1=0.0D0
  DST2=0.0D0
  do J=1,2
   do K=1,2
    do L=1,2
     do M=1,2
      DST1=DST1+5.0D0*XF(J)*XF(K)*XF(L)*XF(M)*DIJKLM(1,J,K,L,M)* &
      VCIJKLM(1,J,K,L,M)**4
      DST2=DST2+5.0D0*XF(J)*XF(K)*XF(L)*XF(M)*DIJKLM(J,2,K,L,M)* &
      VCIJKLM(J,2,K,L,M)**4
     end do
    end do
   end do
  end do
      
  do I=1,2
   do J=1,2
    do K=1,2
     do L=1,2
      do M=1,2
       do N=1,2
        EIJKLMN(I,J,K,L,M,N)= &
        & ((ES(I)*EE(I)**EU3 +ES(J)*EE(J)**EU3 &
        & + ES(K)*EE(K)**EU3 +ES(L)*EE(L)**EU3 &
        & + ES(M)*EE(M)**EU3 +ES(N)*EE(N)**EU3)/6.0D0)**3
        VCIJKLMN(I,J,K,L,M,N)= &
        & ((VC(I)**EU3+VC(J)**EU3 &
        & + VC(K)**EU3+VC(L)**EU3 &
        & + VC(M)**EU3+VC(N)**EU3)/6.0D0)**3
       end do
      end do
     end do
    end do
   end do
  end do
  EVC5=0.0D0
  do I=1,2
   do J=1,2
    do K=1,2
     do L=1,2
      do M=1,2
       do N=1,2
  EVC5=EVC5+XF(I)*XF(J)*XF(K)*XF(L)*XF(M)*XF(N)* &
  EIJKLMN(I,J,K,L,M,N)*VCIJKLMN(I,J,K,L,M,N)**5
       end do
      end do
     end do
    end do
   end do
  end do
  EST1=0.0D0
  EST2=0.0D0
  do J=1,2
   do K=1,2
    do L=1,2
     do M=1,2
      do N=1,2
 EST1=EST1+6.0D0*XF(J)*XF(K)*XF(L)*XF(M)*XF(N)* &
 EIJKLMN(1,J,K,L,M,N)*VCIJKLMN(1,J,K,L,M,N)**5
 EST2=EST2+6.0D0*XF(J)*XF(K)*XF(L)*XF(M)*XF(N)* &
 EIJKLMN(J,2,K,L,M,N)*VCIJKLMN(J,2,K,L,M,N)**5
      end do
     end do
    end do
   end do
  end do
      
  do I=1,2
    do J=1,2
      FIJ(I,J)=((FS(I)*FF(I)**EU3+FS(J)*FF(J)**EU3)/2.0D0)**3
    end do
  end do
  FVC2=0.0D0
  do I=1,2
    do J=1,2
      FVC2=FVC2+XF(I)*XF(J)*FIJ(I,J)*VCIJ(I,J)**2
    end do
  end do
  
  FST1=0.0D0
  FST2=0.0D0
  do J=1,2
    FST1=FST1+2.0D0*XF(J)*FIJ(1,J)*VCIJ(1,J)**2
    FST2=FST2+2.0D0*XF(J)*FIJ(J,2)*VCIJ(J,2)**2
  end do
  
  BETA1=0.0D0
  do I=1,2
    BETA1=BETA1+XF(I)*BETA(I)
  end do
  BESTST1=BETA(1)
  BESTST2=BETA(2)
      
  do I=1,2
    do J=1,2
      do K=1,2
        GAMIJK(I,J,K)=(((GAMMA(I)**EU3+GAMMA(J)**EU3+ &
        GAMMA(K)**EU3)/3.0D0)**3)*K3(I,J,K)
      end do
    end do
  end do
  
  GAMMAVC2=0.0D0
  do I=1,2
    do J=1,2
      do K=1,2
        GAMMAVC2=GAMMAVC2+XF(I)*XF(J)*XF(K)*GAMIJK(I,J,K)* &
        VCIJK(I,J,K)**2
      end do
    end do
  end do
  
  GAMST1=0.0D0
  GAMST2=0.0D0
  do J=1,2
    do K=1,2
      GAMST1=GAMST1+3.0D0*XF(J)*XF(K)*GAMIJK(1,J,K)* &
      VCIJK(1,J,K)**2
      GAMST2=GAMST2+3.0D0*XF(J)*XF(K)*GAMIJK(J,2,K)* &
      VCIJK(J,2,K)**2
    end do
  end do
  
  MWT=XF(1)*PM(1)+XF(2)*PM(2)
  
  return
end subroutine MixRulesLP

subroutine MAKEK1_LO(T,K1,K2,K3)
!--
!-- code (MODIFIED) from THERIAK / C DeCapitani
!--

  implicit none
  
  real(dp),intent(in) :: T
  real(dp),intent(out):: K1(2,2),K2(2,2,2),K3(2,2,2)
  
  real(dp):: T2
  
  !--- bis jetzt nur index 1 und 2
  !---- 1=H2O, 2=CO2, 3=CH4
  K1(:,:)=0.0D0
  K2(:,:,:)=0.0D0
  K3(:,:,:)=0.0D0
  
  T2=T*T
  
  K1(1,1)=1.0D0
  K1(2,2)=1.0D0
  K2(1,1,1)=1.0D0
  K2(2,2,2)=1.0D0
  K3(1,1,1)=1.0D0
  K3(2,2,2)=1.0D0
  
  !----- binary H2O-CO2
  if (T<=373.15) then
    K1(1,2)=0.20611+0.0006*T
    K1(2,1)=K1(1,2)
    K2(1,1,2)=0.8023278-0.0022206*T+184.76824/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    K3(1,1,2)=1.80544-0.0032605*T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  end if
  !--
  if (T > 373.15.and.T<=495.15) then
    K1(1,2)=-10084.5042-4.27134485*T+256477.783/T &
    &       +0.00166997474*T2+1816.78*log(T)
    K1(2,1)=K1(1,2)
    K2(1,1,2)=9.000263-0.00623494*T-2307.7125/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    K3(1,1,2)=-74.1163+0.1800496*T-1.40904946E-04*T2+10130.5246/T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  end if
  !--
  if (T > 495.15.and.T<=623.15) then
    K1(1,2)=-0.3568+7.8888E-04*T+333.399/T
    K1(2,1)=K1(1,2)
    K2(1,1,2)=-19.97444+0.0192515*T+5707.4229/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    K3(1,1,2)=12.1308-0.0099489*T-3042.09583/T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  end if
  !--
  if (T > 623.15.and.T<=672.15) then
    K1(1,2)=-4.53122+0.0042113*T+1619.7/T
    K1(2,1)=K1(1,2)
    K2(1,1,2)=-163.4855+0.190552*T-7.228514E-05*T2+46082.885/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    K3(1,1,2)=1.7137-6.7136E-04*T
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  end if
  
  if (T > 672.15) then
    K1(1,2)=3.131D0-5.0624D-3*T+1.8641D-6*T2-31.409D0/T
    K1(2,1)=K1(1,2)
    K2(1,1,2)=-46.646D0+4.2877D-2*T-1.0892D-5*T2+1.5782D4/T
    K2(1,2,1)=K2(1,1,2)
    K2(2,1,1)=K2(1,1,2)
    K2(1,2,2)=K2(1,1,2)
    K2(2,1,2)=K2(1,1,2)
    K2(2,2,1)=K2(1,1,2)
    K3(1,1,2)=0.9D0
    K3(1,2,1)=K3(1,1,2)
    K3(2,1,1)=K3(1,1,2)
    K3(1,2,2)=K3(1,1,2)
    K3(2,1,2)=K3(1,1,2)
    K3(2,2,1)=K3(1,1,2)
  end if
  
  return
end subroutine MAKEK1_LO

end module M_Mixmodel_Duan
