module M_Eos_KerJac
!--
!-- module for calculations of H2O-CO2 fluids according to Kerrick & Jacobs EoS
!-- modified from CdeCapitani F77 sources (most) and from SupCrt92 sources (some)
!--
  use M_Kinds
  use M_Trace
  use M_Dtb_Const,only: T_CK
  
  implicit none

  private

  public:: EoS_H2O_KerJac      !(J_,Pbar,Tk,G,FugCoef,VOLUM)
  public:: EoS_CO2_KerJac      !(J_,Pbar,Tk,G,FugCoef,VOLUM)
  public:: EoS_H2O_CO2_KerJac  !(Pbar,TdgK,XCO2,FugC_CO2,FugC_H2O,ActCO2,ActH2O)

  real(dp),parameter:: &
  & R_JK=  8.314502D0, & !GasConstant, J/MOLE*K
  & R_CbK= 83.14502D0    !GasConstant, CC*BARS/MOLE*K !1 cm3=  0.1 J/bar

  logical,parameter:: develop= .true. ! .false.  ! !

  !! type t_Sol_KerJac
  !!   character(len=7)::Nam
  !!   real(dp)        ::&
  !!   & H0,S0,&
  !!   & K1,K2,K3,K4,K6,&
  !!   & B,C,D,E
  !! end type t_Sol_KerJac

contains

subroutine EoS_KerJac_CalcPur( &
& iCod,Tk,Pbar, &
& Grt,FugCoef,V_m3)
!--
!-- from THERIAK sources (C. deCapitani), modified
!--
  integer, intent(in) :: iCod
  real(dp),intent(in) :: Pbar,Tk
  real(dp),intent(out):: Grt,FugCoef,V_m3
  !
  real(dp):: H0,S0,K1,K2,K3,K4,K6
  real(dp):: G,Vol,CPDT,CPTDT,Tk0,TT,TT0,SQT,SQT0

  if(iDebug>2) write(fTrc,'(A)') &
  & "===========================================< EoS_KerJac_CalcPur =="

  Tk0=  298.15D0
  TT=   Tk*Tk
  TT0=  Tk0*Tk0
  SQT=  SQRT(Tk)
  SQT0= SQRT(Tk0)

  select case(iCod)

  case(1)
  ! H2O
  ! data in _DtbThr_JUN92_cc.txt
  ! WATER H(2)O(1) W WATR
  ! ST       -228538.00     -241816.00       188.7200          0.000        1.00000
  ! C1        115.45000      -3799.900   -2871300.000             0.        0.00000
  ! C2         51055.16      -0.002062     0.00000000          0.000        0.00000
  ! C1               K1             K4             K3
  ! C2               K6             K2
  ! Cp= K1*T +K2*T*T -K3/T +K4*SQRT(T) +K6*log(T)
    H0= -241816.00D0
    S0=     188.7200D0
    K1=     115.45000D0
    K2=      -0.002062D0
    K3=-2871300.00D0
    K4=   -3799.900D0
    K6=   51055.16D0

  case(2)
  ! CO2
  ! data in _DtbThr_JUN92_cc.txt
  ! CARBON DIOXIDE C(1)O(2) CO2 CARB
  ! ST       -394342.00     -393510.01       213.6770          0.000        1.00000
  ! C1         93.00000      -1340.900     123800.000             0.        0.00000
  ! C2          6336.20      -0.002876     0.00000000          0.000        0.00000
  ! C1               K1             K4             K3
  ! C2               K6             K2
  ! Cp= K1*T +K2*T*T -K3/T +K4*SQRT(T) +K6*log(T)
    H0= -393510.010D0
    S0=     213.677D0
    K1=      93.0D0
    K2=      -0.002876D0
    K3=  123800.0D0
    K4=   -1340.9D0
    K6=    6336.2D0

  end select

  !Cp= K1 +K2.T +K3/Sqr(T) + K4/SqRt(T) + K6/T
  CpdT= K1*(Tk     -Tk0    )       &
  &   + K2*(TT     -TT0    )/2.0D0 &
  &   - K3*(One/Tk -One/Tk0)       &
  &   + K4*(SQT    -SQT0   )*2.0D0 &
  &   + K6*log(Tk/Tk0)

  CpTdT= K1*log(Tk/Tk0)                  &
  &    + K2*(Tk       -Tk0      )        &
  &    - K3*(One/(TT) -One/(TT0)) /2.0D0 &
  &    - K4*(One/SQT  -One/SQT0 ) *2.0D0 &
  &    - K6*(One/Tk   -One/Tk0  )

  ! compute G at T=Tk and P=Pref (=1bar)
  G= H0 +CpdT -Tk*(S0+CpTdT)

  call EoS_KerJac_CalcPhiPur(iCod,Tk,Pbar,FugCoef,Vol)

  V_m3= Vol !!!§§§ apply what conversion factor ???

  ! bring fluid from (Tk,Pref) to (Tk,Pbar)
  ! G(T,P)= G(T,P0) + Integr(P0toP,(@G/@P)dp)
  Grt= G/R_JK/Tk + log(FugCoef*Pbar)

  if(iDebug>2) write(fTrc,'(A)') &
  & "===========================================</EoS_KerJac_CalcPur =="
  
  return
end subroutine EoS_KerJac_CalcPur

subroutine EoS_H2O_KerJac( &
& Tk,Pbar, &
& Grt,FugCoef,V_m3)
!--
!-- from THERIAK sources (C. deCapitani), modified
!--
  real(dp),intent(in) :: Pbar,Tk
  real(dp),intent(out):: Grt,FugCoef,V_m3
  !
  real(dp):: H0,S0,K1,K2,K3,K4,K6
  real(dp):: G,Vol,CPDT,CPTDT,Tk0,TT,TT0,SQT,SQT0

  if(iDebug>2) write(fTrc,'(/,A)') &
  & "===============================================< EoS_H2O_KerJac =="

  Tk0=  298.15D0
  TT=   Tk*Tk
  TT0=  Tk0*Tk0
  SQT=  SQRT(Tk)
  SQT0= SQRT(Tk0)

  ! data in _DtbThr_JUN92_cc.txt
  ! WATER H(2)O(1) W WATR
  ! ST       -228538.00     -241816.00       188.7200          0.000        1.00000
  ! C1        115.45000      -3799.900   -2871300.000             0.        0.00000
  ! C2         51055.16      -0.002062     0.00000000          0.000        0.00000
  ! C1               K1             K4             K3
  ! C2               K6             K2
  ! Cp= K1*T +K2*T*T -K3/T +K4*SQRT(T) +K6*log(T)
  
  H0= -241816.00D0
  S0=     188.7200D0
  K1=     115.45000D0
  K2=      -0.002062D0
  K3=-2871300.00D0
  K4=   -3799.900D0
  K6=   51055.16D0

  !Cp= K1 +K2.T +K3/Sqr(T) + K4/SqRt(T) + K6/T
  CpdT= K1*(Tk     -Tk0    )       &
  &   + K2*(TT     -TT0    )/2.0D0 &
  &   - K3*(One/Tk -One/Tk0)       &
  &   + K4*(SQT    -SQT0   )*2.0D0 &
  &   + K6*log(Tk/Tk0)
  !
  CpTdT= K1*log(Tk/Tk0)                  &
  &    + K2*(Tk       -Tk0      )        &
  &    - K3*(One/(TT) -One/(TT0)) /2.0D0 &
  &    - K4*(One/SQT  -One/SQT0 ) *2.0D0 &
  &    - K6*(One/Tk   -One/Tk0  )

  ! compute G at T=Tk and P=Pref (=1bar)
  G= H0 +CpdT -Tk*(S0+CpTdT)
  !
  if(develop) then ;  call EoS_H2O_KerJac_PhiPur(Tk,Pbar,FugCoef,Vol)
  else             ;  call EoS_KerJac_CalcPhiPur(1,Tk,Pbar,FugCoef,Vol)
  end if

  V_m3= Vol !!!§§§ apply what conversion factor ???

  ! bring fluid from (Tk,Pref) to (Tk,Pbar)
  ! G(T,P)= G(T,P0) + Integr(P0toP,(@G/@P)dp)
  Grt= G/R_JK/Tk + log(FugCoef*Pbar)

  if(iDebug>2) write(fTrc,'(A,/)') &
  & "===============================================</EoS_H2O_KerJac =="
  
  return
end subroutine EoS_H2O_KerJac

subroutine EoS_CO2_KerJac( & !CdCapitani
& Tk,Pbar, &
& Grt,FugCoef,V_m3)

  real(dp),intent(in) :: Pbar,Tk
  real(dp),intent(out):: Grt,FugCoef,V_m3
  !
  real(dp):: H0,S0,K1,K2,K3,K4,K6
  real(dp):: G,Vol,CPDT,CPTDT,Tk0,TT,TT0,SQT,SQT0

  if(iDebug>2) write(fTrc,'(/,A)') &
  & "===============================================< EoS_CO2_KerJac =="

  Tk0=  298.15D0
  TT=   Tk*Tk
  TT0=  Tk0*Tk0
  SQT=  SQRT(Tk)
  SQT0= SQRT(Tk0)

  ! data in _DtbThr_JUN92_cc.txt
  ! CARBON DIOXIDE C(1)O(2) CO2 CARB
  ! ST       -394342.00     -393510.01       213.6770          0.000        1.00000
  ! C1         93.00000      -1340.900     123800.000             0.        0.00000
  ! C2          6336.20      -0.002876     0.00000000          0.000        0.00000
  ! C1               K1             K4             K3
  ! C2               K6             K2
  ! Cp= K1*T +K2*T*T -K3/T +K4*SQRT(T) +K6*log(T)
  
  H0= -393510.010D0
  S0=     213.677D0
  K1=      93.0D0
  K2=      -0.002876D0
  K3=  123800.0D0
  K4=   -1340.9D0
  K6=    6336.2D0

  !Cp= K1 +K2.T +K3/Sqr(T) + K4/SqRt(T) + K6/T
  CpdT= K1*(Tk     -Tk0    )       &
  &   + K2*(TT     -TT0    )/2.0D0 &
  &   - K3*(One/Tk -One/Tk0)       &
  &   + K4*(SQT    -SQT0   )*2.0D0 &
  &   + K6*log(Tk/Tk0)

  CpTdT= K1*log(Tk/Tk0)                  &
  &    + K2*(Tk       -Tk0      )        &
  &    - K3*(One/(TT) -One/(TT0)) /2.0D0 &
  &    - K4*(One/SQT  -One/SQT0 ) *2.0D0 &
  &    - K6*(One/Tk   -One/Tk0  )

  ! compute G at T=Tk and P=Pref (=1bar)
  G= H0 +CpdT -Tk*(S0+CpTdT)

  if(develop) then ;  call EoS_CO2_KerJac_PhiPur(Tk,Pbar,FugCoef,Vol)
  else             ;  call EoS_KerJac_CalcPhiPur(2,Tk,Pbar,FugCoef,Vol)
  end if

  V_m3= Vol !!!§§§ apply what conversion factor ???

  ! bring fluid from (Tk,Pref) to (Tk,Pbar)
  ! G(T,P)= G(T,P0) + Integr(P0toP,(@G/@P)dp)
  Grt= G/R_JK/Tk + log(FugCoef*Pbar)

  if(iDebug>2) write(fTrc,'(A,/)') &
  & "===============================================</EoS_CO2_KerJac =="

  return
end subroutine EoS_CO2_KerJac

subroutine EoS_H2O_CO2_KerJac( & !
& TdgK,Pbar,XCO2,              & !in
& FugC_H2O,FugC_CO2,           & !out
& PhiH2O_Mix,PhiCO2_Mix,       & !out
& ActH2O,ActCO2)                 !out
!--
!-- Reference: Kerrick, D. M. And Jacobs, G. K.,Am.Jour.Sci.
!--   a modified Redlich-Kwong equation for H2O, CO2,and H2O-CO2 mixtures
!--   at elevated pressures and tempuratures
!--
!-- Variables And Rules Are As Follows:
!--
!--   actCO2,actH2O...activity of CO2 and H2O, respectively
!--   bCO2,bH2O,bMix..covolume; cc/mole
!--   cCO2,cH2O,cMix..attractive term in MRK equation; BAR*(CC**2)**2SQRT(T)/MOLE**2
!--   dCO2,dH2O,dMix..attractive term in MRK equation; BAR*(CC**3)*SQRT(T)/MOLE**3
!--   eCO2,eH2O,eMix..attractive term in MRK equation; BAR*(CC**4)*SQRT(T)/MOLE**4
!--   vCO2,vH2O,vMix..molar volume; cc/mole
!--   zCO2,zH2O,zMix..compressibility
!--   CIJ,DIJ,EIJ...  cross coefficients of c,d,e
!--   FKCMix,FKWMix.  fugacity coefficient of CO2 and H2O in the fluid mixture
!--   FKCPur,FKWPur.  fugacity coefficients of pure CO2 and pure H2O, respectively
!--   XC,XCO2.......  mole fraction of co2 in the fluid mixture
!--   XW,XH2O.......  mole fraction of h2o in the fluid mixture
!--   Y.............  b/v4; variable in hard sphere-equation
!--
  real(dp),intent(in) :: Pbar,TdgK,XCO2
  real(dp),intent(out):: FugC_CO2,FugC_H2O
  real(dp),intent(out):: PhiH2O_Mix,PhiCO2_Mix
  real(dp),intent(out):: ActCO2,ActH2O
  
  !
  real(dp):: P,T,RT,TC,PK
  real(dp):: XC,XW
  real(dp):: bMix,cMix,dMix,eMix
  real(dp):: CIJ,DIJ,EIJ
  real(dp):: zCO2,vCO2,zH2O,vH2O,vMix,zMix
  ! real(dp):: PhiCO2_Mix, PhiH2O_Mix
  
  !coeff's of KerJac equ. for CO2 (cCO2) and H2O (cH2O)
  real(dp):: & 
  & bH2O,cCO2,dCO2,eCO2, &
  & bCO2,cH2O,dH2O,eH2O
  
  P=Pbar
  T=TdgK

  !CALCULATION OF parameterS useD IN THE program.
  PK=    P /1.0D3 !-> PK is in kbar
  TC=    T -T_CK
  RT=    R_CbK*T *SQRT(T)

  XC=    XCO2      !=X(CO2)
  XW=    One -XC !=X(H2O)

  ActCO2= XC
  ActH2O= XW

  bH2O=  29.0D0        !Covolume Of H2O
  cH2O=  (   290.78D0  -0.30276D0*T +0.00014774D0*T*T)*1.0D6
  dH2O=  (  -8374.0D0   +19.437D0*T   -0.008148D0*T*T)*1.0D6
  eH2O=  (  76600.0D0    -133.9D0*T     +0.1071D0*T*T)*1.0D6

  bCO2=  58.0D0        !Covolume Of CO2
  cCO2=  (    28.31D0  +0.10721D0*T -0.00000881D0*T*T)*1.0D6
  dCO2=  (   9380.0D0     -8.53D0*T   +0.001189D0*T*T)*1.0D6
  eCO2=  (-368654.0D0    +715.9D0*T     +0.1534D0*T*T)*1.0D6
  
  if (TC >= 325.0D0 .and. TC <= 1050.0D0) then
    bMix=   bCO2*XC + bH2O*XW
    CIJ=    SQRT(cCO2*cH2O)
    DIJ=    SQRT(dCO2*dH2O)
    EIJ=    SQRT(eCO2*eH2O)
    cMix=   cCO2*XC*XC + cH2O*XW*XW + 2.0D0*CIJ*XC*XW
    dMix=   dCO2*XC*XC + dH2O*XW*XW + 2.0D0*DIJ*XC*XW
    eMix=   eCO2*XC*XC + eH2O*XW*XW + 2.0D0*EIJ*XC*XW
  end if

  ! compute Z of pure H2O
  call ZPur(1,bH2O,cH2O,dH2O,eH2O,TdgK,Pbar, & !in
  &         zH2O,vH2O)

  ! compute Z of pure CO2
  call ZPur(2,bCO2,cCO2,dCO2,eCO2,TdgK,Pbar, & !in
  &         zCO2,vCO2)

  ! compute Z of H2O-CO2 mixture
  if (TC>=325.0D0 .and. TC<=1050.0D0) &
  & call ZMixCalc( &
  & TdgK,Pbar,XC,XW,vCO2,vH2O, & !in
  & bMix,cMix,dMix,eMix, vMix,zMix)
  
  if(iDebug>2) write(fTrc,'(7G15.6)') &
  & XCO2,zH2O,zCO2,zMix,vH2O,vCO2,vMix

  ! compute fugacity coefficients of CO2 and H2O
  FugC_H2O= FugCoef_Pure(bH2O,cH2O,dH2O,eH2O,TdgK,Pbar, zH2O,vH2O)
  FugC_CO2= FugCoef_Pure(bCO2,cCO2,dCO2,eCO2,TdgK,Pbar, zCO2,vCO2)

  ! compute Activities, Fugacities in the mixture
  if (TC>=325.0D0 .and. TC<=1050.0D0) then

    ! fug.coeff H2O in mixture
    PhiH2O_Mix= FugCoef_Mix( &
    & 1,T,P, &
    & XC,XW, &
    & bH2O,cH2O,dH2O,eH2O, &
    & bMix,cMix,dMix,eMix, &
    & CIJ,DIJ,EIJ,vMix,zMix)

    ! fug.coeff CO2 in mixture
    PhiCO2_Mix= FugCoef_Mix( &
    & 2,T,P, &
    & XC,XW, &
    & bCO2,cCO2,dCO2,eCO2, &
    & bMix,cMix,dMix,eMix, &
    & CIJ,DIJ,EIJ,vMix,zMix)

    ActH2O=  XW *PhiH2O_Mix /FugC_H2O
    !PhiH2O_Mix*XW /FGH2O with PhiH2O_Mix=fug.coeff H2O in mixture

    ActCO2=  XC *PhiCO2_Mix /FugC_CO2
    !PhiCO2_Mix*XC /FGCO2 with PhiCO2_Mix=fug.coeff CO2 in mixture

  else

    ActCO2=  XCO2
    ActH2O=  XW

  end if

  return
end subroutine EoS_H2O_CO2_KerJac

subroutine EoS_KerJac_CalcPhiPur( &
& iCod,TdgK,Pbar, &
& FugCoef,Volum)
!--
!-- program by Kerrick & Jacobs (1981)
!-- for calculation of free energy of CO2, H2O and for their mixtures
!--
!-- mixtures are restricted to 325 - 1050 C,
!-- because values of c, d, or e become negative
!-- and thus can't take the square root
!--
!-- VARIABLES AND RULES ARE AS FOLLOWS:
!-- TC,T             temperature; celsius,kelvin, respectively
!-- ACO2,aH2O        activity of CO2 and H2O, respectively
!-- cCO2-for CO2; cH2O-for H2O; cMix: for mixture
!-- bCO2,bH2O,bMix   covolume; cc/mole
!-- cCO2,cH2O,cMix   attractive term in MRK equation; BAR*(CC**2)**2SQRT(T)/MOLE**2
!-- dCO2,dH2O,dMix   attractive term in MRK equation; BAR*(CC**3)*SQRT(T)/MOLE**3
!-- eCO2,eH2O,eMix   attractive term in MRK equation; BAR*(CC**4)*SQRT(T)/MOLE**4
!-- CIJ,DIJ,EIJ      cross coefficients of c,d,e
!-- FKCMix,FKWMix    fugacity coefficient of CO2 and H2O in the fluid mixture
!-- FKCPur,FKWPur    fugacity coefficients of pure CO2 and pure H2O, respectively
!-- vCO2,vH2O,vMix   molar volumes; CC/MOLE
!-- XC,XCO2          Mole Fraction Of CO2 In The Fluid Mixture
!-- XW,XH2O          Mole Fraction Of H2O In The Fluid Mixture
!-- Y B/V4           Variable In Hard Sphere-Equation
!-- zCO2,zH2O,zMix   Compressibilities
!--

  integer, intent(in) :: iCod
  real(dp),intent(in) :: Pbar,TdgK
  real(dp),intent(out):: FugCoef,Volum
  !
  real(dp):: P,T,T2
  real(dp):: Z,V,PhiPur
  !
  !coeff's of KerJac equ. for CO2 (cCO2) and H2O (cH2O)
  real(dp):: bH2O,cCO2,dCO2,eCO2
  real(dp):: bCO2,cH2O,dH2O,eH2O
  
  P=  Pbar !bars
  T=  TdgK !Kelvin
  T2= T*T
  
  select case(iCod)
  
  case(1)
    bH2O=  29.0D0   !Covolume Of H2O
    cH2O= (290.78D0    -0.30276D0*T +0.00014774D0*T2)*1.0D6   ! cH2O
    dH2O= (-8374.0D0   +19.437D0 *T -0.008148D0  *T2)*1.0D6   ! dH2O
    eH2O= (76600.0D0   -133.9D0  *T +0.1071D0    *T2)*1.0D6   ! eH2O
    
    call ZPur( &
    & 1,bH2O,cH2O,dH2O,eH2O,T,P,& !In
    & Z,V)       !Out
    
    PhiPur= FugCoef_Pure(bH2O,cH2O,dH2O,eH2O,T,P,Z,V)

  case(2)
    bCO2=  58.0D0   !Covolume Of CO2
    cCO2= (28.31D0     +0.10721D0*T -0.00000881D0*T2)*1.0D6   ! cCO2
    dCO2= (9380.0D0    -8.53D0   *T +0.001189D0  *T2)*1.0D6   ! dCO2
    eCO2= (-368654.0D0 +715.9D0  *T +0.1534D0    *T2)*1.0D6   ! eCO2
    
    call ZPur( &
    & 2,bCO2,cCO2,dCO2,eCO2,T,P,& !In
    & Z,V)       !Out
    
    PhiPur= FugCoef_Pure(bCO2,cCO2,dCO2,eCO2,T,P,Z,V)

  end select
  
  FugCoef=PhiPur
  Volum=V

  return
end subroutine EoS_KerJac_CalcPhiPur

subroutine EoS_H2O_KerJac_PhiPur( &
& TdgK,Pbar, &
& FugCoef,Volum)
!--
  real(dp),intent(in) :: Pbar,TdgK
  real(dp),intent(out):: FugCoef,Volum
  !
  real(dp):: T,T2
  real(dp):: TC,PK,VI
  real(dp):: Z,V
  real(dp):: bH2O,cH2O,dH2O,eH2O
  
  if(iDebug>2) write(fTrc,'(/,A)') &
  & "=========================================< EoS_H2O_KerJac_PhiPur=="
  
  T=  TdgK !Kelvin
  T2= T*T
  !
  bH2O=  29.0D0
  cH2O= (290.78D0    -0.30276D0*T +0.00014774D0*T2)*1.0D6
  dH2O= (-8374.0D0   +19.437D0 *T -0.008148D0  *T2)*1.0D6
  eH2O= (76600.0D0   -133.9D0  *T +0.1071D0    *T2)*1.0D6
  !
  !--- compute initial guess volume for Solve
  TC=   TdgK -T_CK     !-> Celsius
  PK=   Pbar /1000.0D0 !-> kilobar
  !
  if (PK>=One)                                      VI=  22.0D0
  if (PK>=0.90D0 .and. PK<One)                      VI=  24.2D0
  if (PK>=0.60D0 .and. PK<0.9D0)                    VI=  31.2D0
  if (PK>=0.21D0 .and. PK<0.60D0 .and. TC>=550.0D0) VI=  75.0D0
  if (PK>=0.21D0 .and. PK<0.60D0 .and. TC< 550.0D0) VI=  35.0D0
  if (PK>=0.10D0 .and. PK<0.21D0 .and. TC< 400.0D0) VI=  15.0D0
  if (PK>=0.10D0 .and. PK<0.21D0 .and. TC>=400.0D0) VI= 100.0D0
  if (PK>=0.005D0.and. PK<0.10D0)                   VI= 500.0D0
  if (PK< 0.005D0)                                  VI=1000.0D0
  !---/
  
  !!! if(iDebug>2) write(fTrc,'(5G15.6)') bH2O,cH2O,dH2O,eH2O,VI

  call Solve( &
  & TdgK,Pbar, &
  & bH2O,cH2O,dH2O,eH2O,VI, &
  & Z,V)
  !
  if(iDebug>3) write(fTrc,'(2(A,G15.6,/))') "Z=",Z,"V=",V
  !
  FugCoef=FugCoef_Pure(bH2O,cH2O,dH2O,eH2O,T,Pbar,Z,V)
  Volum=V

  if(iDebug>2) write(fTrc,'(A,/)') &
  & "=========================================</EoS_H2O_KerJac_PhiPur=="
  
  return
end subroutine EoS_H2O_KerJac_PhiPur

subroutine EoS_CO2_KerJac_PhiPur( &
& TdgK,Pbar, &
& FugCoef,Volum)
!--
  real(dp),intent(in) :: Pbar,TdgK
  real(dp),intent(out):: FugCoef,Volum
  !
  real(dp):: P,T,T2
  real(dp):: Z,V
  real(dp):: bCO2,cCO2,dCO2,eCO2
  real(dp):: TC,PK,VI
  
  if(iDebug>2) write(fTrc,'(/,A)') &
  & "=========================================< EoS_CO2_KerJac_PhiPur=="
  
  P=  Pbar !bars
  T=  TdgK !Kelvin
  T2= T*T
  !
  bCO2=  58.0D0  !Covolume Of CO2
  cCO2= (28.31D0     +0.10721D0*T -0.00000881D0*T2)*1.0D6   ! cCO2
  dCO2= (9380.0D0    -8.53D0   *T +0.001189D0  *T2)*1.0D6   ! dCO2
  eCO2= (-368654.0D0 +715.9D0  *T +0.1534D0    *T2)*1.0D6   ! eCO2
  !
  !--- compute initial guess volume for Solve
  TC=   TdgK -T_CK     !-> Celsius
  PK=   Pbar /1000.0D0 !-> kilobar
  !
  if (PK >= One)                        VI=  35.0D0
  if (PK >= 0.10D0  .and. PK < One)     VI= 100.0D0
  if (PK >= 0.005D0 .and. PK < 0.10D0)  VI= 500.0D0
  if (PK <  0.005D0)                    VI=5000.0D0
  !---/
  call Solve( &
  & TdgK,Pbar, &
  & bCO2,cCO2,dCO2,eCO2,VI, &
  & Z,V)
  !
  if(iDebug>3) write(fTrc,'(2(A,G15.6,/))') "Z=",Z,"V=",V
  !
  FugCoef=FugCoef_Pure(bCO2,cCO2,dCO2,eCO2,T,P,Z,V)
  Volum=V

  if(iDebug>2) write(fTrc,'(A,/)') &
  & "=========================================</EoS_CO2_KerJac_PhiPur=="
  
  return
end subroutine EoS_CO2_KerJac_PhiPur

subroutine ZPur( & !
& iCod,          & !
& B,C,D,E,       & !
& T,P,           & !
& Zz,Vv)           !
!--
!-- compute the volume of CO2 and H2O at T,P
!-- an initial guess of the volume (vi) is chosen for a given pressure range.
!-- this value is then used in the routine Solve,
!-- which solves for the exact volume by means of an iterative newton - raphson technique.
!--
  integer, intent(in) :: iCod
  real(dp),intent(in) :: T,P !P=bar, T=K
  real(dp),intent(in) :: B,C,D,E
  real(dp),intent(out):: Zz,Vv !zCO2,zH2O,vCO2,vH2O
  real(dp):: PK,TC
  real(dp):: VI

  !--- compute initial guess volume for Solve
  TC=   T -T_CK     !-> Celsius
  PK=   P /1000.0D0 !-> kilobar
  !
  select case(iCod)
  
  case(1) ! H2O
    !B=  bH2O ; C=  cH2O ; D=  dH2O ; E=  eH2O
    if (PK>=One)                                      VI=  22.0D0
    if (PK>=0.90D0 .and. PK<One)                      VI=  24.2D0
    if (PK>=0.60D0 .and. PK<0.9D0)                    VI=  31.2D0
    if (PK>=0.21D0 .and. PK<0.60D0 .and. TC>=550.0D0) VI=  75.0D0
    if (PK>=0.21D0 .and. PK<0.60D0 .and. TC< 550.0D0) VI=  35.0D0
    if (PK>=0.10D0 .and. PK<0.21D0 .and. TC< 400.0D0) VI=  15.0D0
    if (PK>=0.10D0 .and. PK<0.21D0 .and. TC>=400.0D0) VI= 100.0D0
    if (PK>=0.005D0.and. PK<0.10D0)                   VI= 500.0D0
    if (PK< 0.005D0)                                  VI=1000.0D0
  !_
  case(2) ! CO2
    !B=  bCO2 ; C=  cCO2 ; D=  dCO2 ; E=  eCO2
    if (PK >= One)                        VI=  35.0D0
    if (PK >= 0.10D0  .and. PK < One)     VI= 100.0D0
    if (PK >= 0.005D0 .and. PK < 0.10D0)  VI= 500.0D0
    if (PK <  0.005D0)                    VI=5000.0D0
  !_
  end select
  !---/

  call Solve( &
  & T,P, &
  & B,C,D,E,VI, &
  & Zz,Vv)

  if(iDebug>3) write(fTrc,'(2(A,G15.6,/))') "Z=",Zz,"V=",Vv
  
  return
end subroutine ZPur

subroutine ZMixCalc( &
& T,P, &
& XC,XW,vCO2,vH2O, &
& bMix,cMix,dMix,eMix, &
& VMix,ZMix)
!--
!-- calculate the volume of a mixture of CO2 and H2O at P,T,Xi
!--
!-- the molar volumes of CO2 and H2O as calculated in ZPur
!-- are used to define the initial estimate
!--
!-- routine Solve is then used to calculate the volume of the mixture.
!--
  real(dp),intent(in) :: P,T,XC,XW,vCO2,vH2O,bMix,cMix,dMix,eMix
  real(dp),intent(out):: VMix,ZMix
  real(dp):: VInit !,Z,V

  VInit= vCO2*XC + vH2O*XW
  ! VInit is VIdeal = Sum(X_i*V_i)

  call Solve( &
  & T,P,&
  & bMix,cMix,dMix,eMix,VInit,&
  & ZMix,VMix)

  return
end subroutine ZMixCalc

subroutine Solve( & 
& T,P,          & ! IN
& B,C,D,E,VInit,& ! IN
& Z,V)            ! OUT
!--
!-- NEWTON-RAPHSON method called By ZPur And Zmix
!-- to Calculate The Volume Of CO2, H2O, And A Mixture Of CO2-H2O
!-- at Pressure, Temperature, and XCO2
!--
!-- Newton-Raphson Method Summarized As Follows:
!--   f(x)-  0;  x(k+1)-  x(k) - f(x)/df(x)
!--   where df(x) is the partial differential of f(x) with respect to x.
!--   x(k+1) is calculated for "k" iterations,
!--   until no change in x(k+1) occurs.

  real(dp),intent(in) :: P,T
  real(dp),intent(in) :: B,C,D,E
  real(dp),intent(in) :: VInit
  real(dp),intent(out):: Z,V
  !
  real(dp):: VI,X,Y
  real(dp):: BI,BI2,PN,PR,PA1,PA2,F
  real(dp):: D1,D2,D3,D4,D5,D6,D7,D8,D9
  real(dp):: DPR,DPA,DF,DifF
  integer ::K

  VI=VInit

  do K=1,50

    Y=    B/VI/4.0D0
    X=    One - Y
    BI=   VI + B
    BI2=  BI*BI

    ! DEFINITION OF THE F(X) FOR Solve:
    ! F(X)=  0=  P(REPULSIVE) - P(ATTRACTIVE) - P
    ! where F(X) IS A REARRANGEMENT OF KERRICK AND JACOBS'(1980) EQUATION (14)
    PN=   One +Y +Y**2 -Y**3
    PR=   R_CbK*T *PN /VI /X**3
    PA1=  C + D/VI + E/VI/VI
    PA2=  PA1/SQRT(T)/VI/BI
    F=    PR - PA2 - P

    ! DEFINITION OF THE DifFERENTIAL OF F(X) FOR NEWARP:
    ! DF(X)=  dP_(REPULSIVE) - dP_(ATTRACTIVE)
    ! D1=   (-3.0D0*B) / (4.0D0*(VI**3)*X**4)
    ! D2=   -One / ((VI**2)*(X**3))
    ! D3=    One / (VI*(X**3))
    ! D4=   -B / (4.0D0*VI**2)
    ! D5=   -2.0D0*(B**2) / (16.0D0*(VI**3))
    ! D6=    3.0D0*(B**3) / (64.0D0*(VI**4))
    ! DPR=  ((PN*(D1+D2)) + (D3*(D4+D5+D6)))*R*T
    ! D7=   (-One/(VI*BI2)) + (-One/(VI**2*BI))
    ! D8=   One / (VI*BI)
    ! D9=   (-D/VI**2) + ((-2.0D0*E)/(VI**3))
    ! DPA=  (PA1*D7+D8*D9) / T12
    ! DF=   DPR - DPA

    D1=   -B    /VI**3 /X**4 *3.0D0/4.0D0
    D2=   -One  /VI**2 /X**3
    D3=    One  /VI    /X**3
    D4=   -B    /VI**2      /4.0D0
    D5=   -B**2 /VI**3      *2.0D0/16.0D0
    D6=    B**3 /VI**4      *3.0D0/64.0D0
    DPR=  (PN*(D1+D2)   +  D3*(D4+D5+D6)) *R_CbK*T
    D7=   -One/VI/BI2 -  One/VI**2/BI
    D8=   One /VI /BI
    D9=   -    D/VI**2  -  2.0D0*E/VI**3
    DPA=  (PA1*D7+D8*D9) / SQRT(T)
    DF=   DPR - DPA

    ! CALCULATION OF V(K+1)
    V=     VI - F/DF
    DifF=  ABS(V-VI)
    
    ! CONTINUATION OR end OF ITERATIONS FOR NEWARP
    if (DifF<0.01D0) exit

    if (V>1.0D6) V=1.0D6
    if (V<9.9D0) V=1.0D1

    VI= V

  enddo

  Z=  V*P /R_CbK/T

  return
end subroutine Solve

real(dp) function FugCoef_Pure(B,C,D,E,T,P,Z,V) !modified
!--
!-- calculate fugacity coefficients of pure co2 or pure h2o at T,P
!-- kerrick and jacobs
!-- for a derivation of the pure fugacity coefficient expression (fcp)
!--
  real(dp),intent(in):: P,T
  real(dp),intent(in):: B,C,D,E
  real(dp),intent(in):: Z,V
  !
  real(dp):: RT,Y,FCP
  !
  RT= R_CbK*T *SQRT(T)

  !select case(iCod)
  !  case(1); B=bH2O; C=cH2O; D=dH2O; E=eH2O !H2O
  !  case(2); B=bCO2; C=cCO2; D=dCO2; E=eCO2 !CO2
  !end select

  Y=    B / (4.0D0*V)

  FCP= (8.0D0*Y -9.0D0*Y*Y +3.0D0*Y**3)       &
  &    / ((One-Y)**3)                         &
  &   - log(Z)                                &
  &   - (C/(RT*(V+B)))      -(D/(RT*V*(V+B))) &
  &   - (E/(RT*V*V*(V+B)))  +((C/(RT*B))*(log(V/(V+B))))    &
  &   - (D/(RT*B*V))        +((D/(RT*B*B))*(log((V+B)/V)))  &
  &   - (E/(RT*2.0D0*B*V*V))+(E/(RT*B*B*V))   &
  &   -((E/(RT*B**3))*(log((V+B)/V)))
  FCP=  exp(FCP)

  ! if(iDebug>2) write(fTrc,'(D12.3)') FCP

  !if(J_==1) FKWPur=  FCP
  !if(J_==2) FKCPur=  FCP

  FugCoef_Pure= FCP

  return
end function FugCoef_Pure

real(dp) function FugCoef_Mix(&
& iCod, &
& T,P,  &
& XC,XW,&
& bPur,cPur,dPur,ePur, &
& bMix,cMix,dMix,eMix, &
& CIJ,DIJ,EIJ,vMix,zMix)
!--
!-- fugacity coefficients of CO2 and H2O in a H2O-CO2 mixture.
!-- these together with the fugacity coefficients of pure CO2 and Pure H2O
!-- can be used to calculate activities of CO2 And H2O in the mixture
!-- for each pressure, temperature, and XCO2.
!--
!-- see Kerrick and Jacobs for a derivation
!-- of fugacity coefficient expression (Fcm-FugCMix)
!--
  integer, intent(in) :: iCod
  real(dp),intent(in) :: T,P
  real(dp),intent(in) :: XC,XW
  real(dp),intent(in) :: bPur,cPur,dPur,ePur
  real(dp),intent(in) :: bMix,cMix,dMix,eMix
  real(dp),intent(in) :: CIJ,DIJ,EIJ,vMix,zMix

  real(dp):: RT
  real(dp):: B,C,D,E
  real(dp):: B1,C1,D1,E1
  real(dp):: X1,X2,FCM,VBV,VB,V,Y,Z

  RT= R_CbK*T*SQRT(T)
  B=  bMix; C=  cMix; D=  dMix; E=  eMix
  V=  vMix; Z=  zMix
  Y=  B /V /4.0D0
  
  B1=bPur; C1=cPur; D1=dPur; E1=ePur
  !select case(iCod)
  !  case(1) ; B1=bH2O; C1=cH2O; D1=dH2O; E1=eH2O; X1=XW; X2=XC !H2O
  !  case(2) ; B1=bCO2; C1=cCO2; D1=dCO2; E1=eCO2; X1=XC; X2=XW !CO2
  !end select
  select case(iCod)
    case(1) ; X1=XW; X2=XC !H2O
    case(2) ; X1=XC; X2=XW !CO2
  end select

  VB= V+B
  VBV=VB/V

  FCM=  (4.0D0*Y -3.0D0*Y*Y) /(One-Y)**2                            &
  &  + ((B1/B)*((4.0D0*Y-2.0D0*Y*Y)/((One-Y)**3)))                  &
  &  - (((2.0D0*C1*X1+2.0D0*CIJ*X2)/(RT*B))       *log(VBV))        &
  &  - (                     (C*B1)/(RT*B*VB))                      &
  &  + ((                    (C*B1)/(RT*B*B))     *log(VBV))        &
  &  - ((2.0D0*D1*X1+2.0D0*DIJ*X2+D)/(RT*B*V))                      &
  &  + (((2.0D0*X1*D1+2.0D0*DIJ*X2+D)/(RT*B*B))   *log(VBV))        &
  &  + ((D*B1)/(RT*V*B*VB))                                         &
  &  + ((2.0D0*B1*D) /(RT*B*B*VB))                                  &
  &  - (((2.0D0*B1*D)/(RT*B**3     ))*(log(VBV)))                   &
  &  - ((2.0D0*E1*X1 +2.0D0*EIJ*X2 +2.0D0*E)/(RT*2.0D0*B*V*V))      &
  &  + ((2.0D0*E1*X1 +2.0D0*EIJ*X2 +2.0D0*E)/(RT*B*B*V))            &
  &  - (((2.0D0*E1*X1+2.0D0*EIJ*X2 +2.0D0*E)/(RT*(B**3)))*log(VBV)) &
  &  +         (E*B1) /(RT*2.0D0*B*V*V*VB)                          &
  &  - ( (3.0D0*E*B1) /(RT*2.0D0*B*B*V*VB))                         &
  &  + (((3.0D0*E*B1) /(RT*(B**4))      )*log(VBV))                 &
  &  - ( (3.0D0*E*B1) /(RT*(B**3)*VB))                              &
  &  - log(Z)

  FugCoef_Mix= exp(FCM)

  !if (J==1) PhiH2O_Mix=  FCM
  !if (J==2) PhiCO2_Mix=  FCM

  return
end function FugCoef_Mix

end module M_Eos_KerJac

!! subroutine ZPur_H2O( &
!! & T,P, &
!! & B,C,D,E, &
!! & Zz,Vv)
!! !--
!! !-- compute the volume of H2O at T,P
!! !-- an initial guess of the volume (vi) is chosen for a given pressure range.
!! !-- this value is then used in the routine Solve,
!! !-- which solves for the exact volume by means of an iterative newton - raphson technique.
!! !--
!!   real(dp),intent(in) :: T,P !P=bar, T=K
!!   real(dp),intent(in) :: B,C,D,E
!!   real(dp),intent(out):: Zz,Vv
!!   !
!!   real(dp):: PK,TC
!!   real(dp):: VI
!!   
!!   !--- compute initial guess volume for Solve
!!   TC=   T -T_CK     !-> Celsius
!!   PK=   P /1000.0D0 !-> kilobar
!!   !
!!   if (PK>=One)                                      VI=  22.0D0
!!   if (PK>=0.90D0 .and. PK<One)                      VI=  24.2D0
!!   if (PK>=0.60D0 .and. PK<0.9D0)                    VI=  31.2D0
!!   if (PK>=0.21D0 .and. PK<0.60D0 .and. TC>=550.0D0) VI=  75.0D0
!!   if (PK>=0.21D0 .and. PK<0.60D0 .and. TC< 550.0D0) VI=  35.0D0
!!   if (PK>=0.10D0 .and. PK<0.21D0 .and. TC< 400.0D0) VI=  15.0D0
!!   if (PK>=0.10D0 .and. PK<0.21D0 .and. TC>=400.0D0) VI= 100.0D0
!!   if (PK>=0.005D0.and. PK<0.10D0)                   VI= 500.0D0
!!   if (PK< 0.005D0)                                  VI=1000.0D0
!!   !---/
!! 
!!   call Solve( &
!!   & T,P, &
!!   & B,C,D,E,VI, &
!!   & Zz,Vv)
!! 
!!   if(iDebug>2) write(fTrc,'(2(A,D12.3))') "ZPur_H2O, Z=",Zz," / V=",Vv
!! 
!!   return
!! end subroutine ZPur_H2O

!! subroutine ZPur_CO2( &
!! & T,P, &
!! & B,C,D,E, &
!! & Zz,Vv)
!! !--
!! !-- compute the volume of CO2 at T,P
!! !-- an initial guess of the volume (vi) is chosen for a given pressure range.
!! !-- this value is then used in the routine Solve,
!! !-- which solves for the exact volume by means of an iterative newton - raphson technique.
!! !--
!!   real(dp),intent(in) :: T,P !P=bar, T=K
!!   real(dp),intent(in) :: B,C,D,E
!!   real(dp),intent(out):: Zz,Vv
!!   !
!!   real(dp):: PK,TC
!!   real(dp):: VI
!! 
!!   !--- compute initial guess volume for Solve
!!   TC=   T -T_CK     !-> Celsius
!!   PK=   P /1000.0D0 !-> kilobar
!!   !
!!   if (PK >= One)                        VI=  35.0D0
!!   if (PK >= 0.10D0  .and. PK < One)     VI= 100.0D0
!!   if (PK >= 0.005D0 .and. PK < 0.10D0)  VI= 500.0D0
!!   if (PK <  0.005D0)                    VI=5000.0D0
!!   !---/
!! 
!!   call Solve( &
!!   & T,P, &
!!   & B,C,D,E,VI, &
!!   & Zz,Vv)
!! 
!!   if(iDebug>2) write(fTrc,'(2(A,D12.3))') "ZPur_CO2, Z=",Zz," / V=",Vv
!! 
!!   return
!! end subroutine ZPur_CO2


