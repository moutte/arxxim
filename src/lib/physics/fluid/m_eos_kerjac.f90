module M_Eos_KerJac
!--
!-- module for calculations of H2O-CO2 fluids according to Kerrick & Jacobs EoS
!-- MODIFIED from CdeCapitani F77 sources (most) and from SupCrt92 sources (some)
!--
  use M_Kinds
  use M_Trace,only: fTrc,iDebug
  use M_Dtb_Const,only: T_CK
  
  implicit none

  private
  
  public:: EoS_H2O_KerJac_PhiPur
  public:: EoS_CO2_KerJac_PhiPur
  ! public:: EoS_KerJac_CalcPur  !(i,TdgK,Pbar)    (Grt,FugCoef,V_m3)
  public:: EoS_H2O_KerJac      !(TdgK,Pbar)      (Grt,FugCoef,V_m3)
  public:: EoS_CO2_KerJac      !(TdgK,Pbar)      (Grt,FugCoef,V_m3)
  public:: EoS_H2O_CO2_KerJac  !(TdgK,Pbar,XCO2) (FugC_CO2,FugC_H2O,ActCO2,ActH2O)

  real(dp),parameter:: & !
  & R_JK=  8.314502D0, & !GasConstant, J/MOLE/K
  & R_CbK= 83.14502D0    !GasConstant, CC*BARS/MOLE/K !1 cm3=  0.1 J/bar

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
& Grt,V_m3,LnPhi)
!--
!-- from THERIAK sources (C. deCapitani), MODIFIED
!-- H20= iCod=1; CO2= iCod=2
!--
  integer, intent(in) :: iCod
  real(dp),intent(in) :: Pbar,Tk
  real(dp),intent(out):: Grt,V_m3,LnPhi
  !
  real(dp):: H0,S0,K1,K2,K3,K4,K6
  real(dp):: S0Ele
  real(dp):: G,Vcm3,CPDT,CPTDT
  real(dp):: Tk0,TT,TT0,SQT,SQT0

  if(iDebug>2) write(fTrc,'(A)') &
  & "------------------------------------------------EoS_KerJac_CalcPur"

  Tk0=  298.15D0
  TT=   Tk*Tk
  TT0=  Tk0*Tk0
  SQT=  sqrt(Tk)
  SQT0= sqrt(Tk0)

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
  ! Cp= K1*T +K2*T*T -K3/T +K4*sqrt(T) +K6*log(T)
    H0= -241816.00D0
    S0=     188.7200D0
    K1=     115.45000D0
    K2=      -0.002062D0
    K3=-2871300.00D0
    K4=   -3799.900D0
    K6=   51055.16D0
    S0Ele= 65.34D0 *2.0D0 + 102.575D0

  case(2)
  ! CO2
  ! data in _DtbThr_JUN92_cc.txt
  ! CARBON DIOXIDE C(1)O(2) CO2 CARB
  ! ST       -394342.00     -393510.01       213.6770          0.000        1.00000
  ! C1         93.00000      -1340.900     123800.000             0.        0.00000
  ! C2          6336.20      -0.002876     0.00000000          0.000        0.00000
  ! C1               K1             K4             K3
  ! C2               K6             K2
  ! Cp= K1*T +K2*T*T -K3/T +K4*sqrt(T) +K6*log(T)
    H0= -393510.010D0
    S0=     213.677D0
    K1=      93.0D0
    K2=      -0.002876D0
    K3=  123800.0D0
    K4=   -1340.9D0
    K6=    6336.2D0
    S0Ele= 5.74D0  + 102.575D0 *2.0D0

  end select

  !Cp= K1 +K2.T +K3/Sqr(T) + K4/sqrt(T) + K6/T
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

  ! Benson convention
  G= G + S0Ele*298.15D0

  call EoS_KerJac_PhiPur(iCod,Tk,Pbar,        Vcm3,LnPhi)

  ! bring fluid from (Tk,Pref) to (Tk,Pbar)
  ! G(T,P)= G(T,P0) + Integr(P0toP,(@G/@P)dp)
  Grt= G/R_JK/Tk + LnPhi + log(Pbar)
  V_m3= Vcm3*1.D-6 !!!§§§ apply what conversion factor ???


  if(iDebug>2) write(fTrc,'(A)') &
  & "-----------------------------------------------/EoS_KerJac_CalcPur"
  
  return
end subroutine EoS_KerJac_CalcPur

subroutine EoS_H2O_KerJac( &
& Tk,Pbar, &
& Grt,V_m3,LnPhi)
!--
!-- from THERIAK sources (C. deCapitani), MODIFIED
!--
  real(dp),intent(in) :: Pbar,Tk
  real(dp),intent(out):: Grt,V_m3,LnPhi
  !
  real(dp):: H0,S0,S0Ele
  real(dp):: K1,K2,K3,K4,K6
  real(dp):: G,Vcm3,CPDT,CPTDT,Tk0,TT,TT0,SQT,SQT0

  if(iDebug>2) write(fTrc,'(/,A)') &
  & "----------------------------------------------------EoS_H2O_KerJac"

  Tk0=  298.15D0
  TT=   Tk*Tk
  TT0=  Tk0*Tk0
  SQT=  sqrt(Tk)
  SQT0= sqrt(Tk0)

  ! data in _DtbThr_JUN92_cc.txt
  ! WATER H(2)O(1) W WATR
  ! ST       -228538.00     -241816.00       188.7200          0.000        1.00000
  ! C1        115.45000      -3799.900   -2871300.000             0.        0.00000
  ! C2         51055.16      -0.002062     0.00000000          0.000        0.00000
  ! C1               K1             K4             K3
  ! C2               K6             K2
  ! Cp= K1*T +K2*T*T -K3/T +K4*sqrt(T) +K6*log(T)
  
  H0= -241816.00D0
  S0=     188.7200D0
  K1=     115.45000D0
  K2=      -0.002062D0
  K3=-2871300.00D0
  K4=   -3799.900D0
  K6=   51055.16D0

  !Cp= K1 +K2.T +K3/Sqr(T) + K4/sqrt(T) + K6/T
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
  if(develop) then ;  call EoS_H2O_KerJac_PhiPur(Tk,Pbar,    Vcm3,LnPhi)
  else             ;  call EoS_KerJac_PhiPur(1,Tk,Pbar,      Vcm3,LnPhi)
  end if

  ! Benson convention
  S0Ele= 102.575D0 + 65.34D0 *2.0D0
  G= G + S0Ele*Tk0

  ! bring fluid from (Tk,Pref) to (Tk,Pbar)
  ! G(T,P)= G(T,P0) + Integr(P0toP,(@G/@P)dp)
  Grt= G/R_JK/Tk + LnPhi + log(Pbar)
  ! for pure fluid, Fugacity(T,P)= FugCoeff(T,P) *Pressure
  V_m3= Vcm3*1.D-6 !!!§§§ apply what conversion factor ???
  
  if(iDebug>2) write(fTrc,'(A,/)') &
  & "---------------------------------------------------/EoS_H2O_KerJac"
  
  return
end subroutine EoS_H2O_KerJac

subroutine EoS_CO2_KerJac( & !CdeCapitani
& Tk,Pbar, &
& Grt,V_m3,LnPhi)
!-----------------------------------------------------------------------
!-- from THERIAK sources (C. deCapitani), MODIFIED
!-----------------------------------------------------------------------
  real(dp),intent(in) :: Pbar,Tk
  real(dp),intent(out):: Grt,V_m3,LnPhi
  !---------------------------------------------------------------------
  real(dp):: H0,S0,S0Ele
  real(dp):: K1,K2,K3,K4,K6
  real(dp):: G,Vcm3,CPDT,CPTDT
  real(dp):: Tk0,TT,TT0,SQT,SQT0
  !---------------------------------------------------------------------
  if(iDebug>2) write(fTrc,'(/,A)') &
  & "----------------------------------------------------EoS_CO2_KerJac"

  Tk0=  298.15D0
  TT=   Tk**2
  TT0=  Tk0**2
  SQT=  sqrt(Tk)
  SQT0= sqrt(Tk0)

  ! data in _DtbThr_JUN92_cc.txt
  ! CARBON DIOXIDE C(1)O(2) CO2 CARB
  ! ST       -394342.00     -393510.01       213.6770          0.000        1.00000
  ! C1         93.00000      -1340.900     123800.000             0.        0.00000
  ! C2          6336.20      -0.002876     0.00000000          0.000        0.00000
  ! C1               K1             K4             K3
  ! C2               K6             K2
  ! Cp= K1*T +K2*T*T -K3/T +K4*sqrt(T) +K6*log(T)
  
  H0= -393510.010D0
  S0=     213.677D0
  K1=      93.0D0
  K2=      -0.002876D0
  K3=  123800.0D0
  K4=   -1340.9D0
  K6=    6336.2D0

  !Cp= K1 +K2.T +K3/Sqr(T) + K4/sqrt(T) + K6/T
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

  !---------------------------------compute G at T=Tk and P=Pref (=1bar)
  G= H0 +CpdT -Tk*(S0+CpTdT)
  !----------------------------------------------------Benson convention
  S0Ele= 5.74D0  + 102.575D0 *2.0D0
  G= G + S0Ele*Tk0
  !
  if(develop) then ;  call EoS_CO2_KerJac_PhiPur(Tk,Pbar,    Vcm3,LnPhi)
  else             ;  call EoS_KerJac_PhiPur(2,Tk,Pbar,      Vcm3,LnPhi)
  end if
  !------------------------------bring fluid from (Tk,Pref) to (Tk,Pbar)
  !----------------------------G(T,P)= G(T,P0) + Integr(P0toP,(@G/@P)dp)
  Grt= G/R_JK/Tk + LnPhi + log(Pbar)

  V_m3= Vcm3*1.D-6 !!!§§§ apply what conversion factor ???

  if(iDebug>2) write(fTrc,'(A,/)') &
  & "---------------------------------------------------/EoS_CO2_KerJac"

  return
end subroutine EoS_CO2_KerJac

subroutine EoS_H2O_CO2_KerJac( & !
& TdgK,Pbar,XCO2,              & !in
& Update_Pure,                 & !in
& LnPhiH2O_Pur,LnPhiCO2_Pur,   & !inout
& V_H2O_Pur,V_CO2_Pur,         & !inout !volumes in cm3 !!!
& LnPhiH2O_Mix,LnPhiCO2_Mix,   & !out
& LnActH2O,    LnActCO2,       & !out
& V_m3)                          !out
!--
!-- Reference: Kerrick, D. M. And Jacobs, G. K.,Am.Jour.Sci.
!--   a MODIFIED Redlich-Kwong equation for H2O, CO2,and H2O-CO2 mixtures
!--   at elevated pressures and tempuratures
!--
!--   LnactCO2,LnactH2O=
!--     log_activity of CO2 and H2O, respectively
!--   LnPhiH2O_Mix,LnPhiCO2_Mix=
!--     log_fugacity coefficient of CO2 and H2O in the fluid mixture
!--   LnPhiH2O_Pur,LnPhiCO2_Pur=
!--     log_fugacity coefficients of pure CO2 and pure H2O, respectively
!--   XC,XCO2.........mole fraction of co2 in the fluid mixture
!--   XW,XH2O.........mole fraction of h2o in the fluid mixture
!--
!--   bCO2,bH2O,bMix..covolume; cc/mole
!--   cCO2,cH2O,cMix..attractive term in MRK equation; BAR*(CC**2)**2SQRT(T)/MOLE**2
!--   dCO2,dH2O,dMix..attractive term in MRK equation; BAR*(CC**3)*sqrt(T)/MOLE**3
!--   eCO2,eH2O,eMix..attractive term in MRK equation; BAR*(CC**4)*sqrt(T)/MOLE**4
!--   vCO2,vH2O,vMix..molar volume; cc/mole
!--   zCO2,zH2O,zMix..compressibility
!--   CIJ,DIJ,EIJ.....cross coefficients of c,d,e
!--   Y...............b/v4; variable in hard sphere-equation
!--
  real(dp),intent(in) :: Pbar,TdgK,XCO2
  logical, intent(in) :: Update_Pure
  real(dp),intent(inout):: LnPhiH2O_Pur,LnPhiCO2_Pur
  real(dp),intent(inout):: V_H2O_Pur,V_CO2_Pur
  real(dp),intent(out):: LnPhiH2O_Mix,LnPhiCO2_Mix
  real(dp),intent(out):: LnActCO2,LnActH2O
  real(dp),intent(out):: V_m3
  !---------------------------------------------------------------------
  real(dp):: P,T,RT,TC,T2,PK,Vcm3
  real(dp):: XC,XW
  real(dp):: bMix,cMix,dMix,eMix
  real(dp):: CIJ,DIJ,EIJ
  real(dp):: zCO2,zH2O
  real(dp):: vMix,zMix
  !---------------------------------------------------------------------
  !coeff's of KerJac equ. for CO2 (cCO2) and H2O (cH2O)
  real(dp):: bCO2,cCO2,dCO2,eCO2
  real(dp):: bH2O,cH2O,dH2O,eH2O
  !
  zMix= 1.0D0
  !
  P=  Pbar
  T=  TdgK
  T2= T*T
  !
  !calculation of parameters used in the program.
  PK=    P /1.0D3 !-> PK is in kbar
  TC=    T -T_CK
  RT=    R_CbK*T *sqrt(T)
  !
  XC=    XCO2      !=X(CO2)
  XW=    One -XC   !=X(H2O)
  !
  LnActCO2= log(XC)
  LnActH2O= log(XW)
  !
  bH2O=  29.0D0        !Covolume Of H2O
  cH2O=  (   290.78D0  -0.30276D0*T +0.00014774D0*T2)*1.0D6
  dH2O=  (  -8374.0D0   +19.437D0*T   -0.008148D0*T2)*1.0D6
  eH2O=  (  76600.0D0    -133.9D0*T     +0.1071D0*T2)*1.0D6
  !
  bCO2=  58.0D0        !Covolume Of CO2
  cCO2=  (    28.31D0  +0.10721D0*T -0.00000881D0*T2)*1.0D6
  dCO2=  (   9380.0D0     -8.53D0*T   +0.001189D0*T2)*1.0D6
  eCO2=  (-368654.0D0    +715.9D0*T     +0.1534D0*T2)*1.0D6
  !
  if(Update_Pure) then
    !
    call ZPur( & !------------------------------------Z,Vcm3 of pure H2O
    & 1,bH2O,cH2O,dH2O,eH2O,TdgK,Pbar, & !in
    & zH2O,V_H2O_Pur)
    !
    LnPhiH2O_Pur= LnFugCoef_Pure(bH2O,cH2O,dH2O,eH2O,TdgK,zH2O,V_H2O_Pur)
    !
    call ZPur( & !------------------------------------Z,Vcm3 of pure CO2
    & 2,bCO2,cCO2,dCO2,eCO2,TdgK,Pbar, & !in
    & zCO2,V_CO2_Pur)
    !
    LnPhiCO2_Pur= LnFugCoef_Pure(bCO2,cCO2,dCO2,eCO2,TdgK,zCO2,V_CO2_Pur)
    !
  end if
  !--------------------------------------------cubic param's for mixture
  if (TC >= 325.0D0 .and. TC <= 1050.0D0) then
    bMix=   bCO2*XC + bH2O*XW
    CIJ=    sqrt(cCO2*cH2O)
    DIJ=    sqrt(dCO2*dH2O)
    EIJ=    sqrt(eCO2*eH2O)
    cMix=   cCO2*XC*XC + cH2O*XW*XW + 2.0D0*CIJ*XC*XW
    dMix=   dCO2*XC*XC + dH2O*XW*XW + 2.0D0*DIJ*XC*XW
    eMix=   eCO2*XC*XC + eH2O*XW*XW + 2.0D0*EIJ*XC*XW
  end if
  !-------------------------------------------------------------------//
  
  !--------------------------------------------Z,Vcm3 of H2O-CO2 mixture
  if (TC>=325.0D0 .and. TC<=1050.0D0) &
  & call ZMixCalc( & !
  & TdgK,Pbar,XC,XW,V_CO2_Pur,V_H2O_Pur, & !in
  & bMix,cMix,dMix,eMix,                 & !in
  & vMix,zMix)
  
  if(iDebug>2) write(fTrc,'(7G15.6)') &
  & XCO2,zH2O,zCO2,zMix,V_H2O_Pur,V_CO2_Pur,vMix

  !--------------------------------Activities, Fugacities in the mixture
  if (TC>=325.0D0 .and. TC<=1050.0D0) then
    !-------------------------------------------fug.coeff H2O in mixture
    LnPhiH2O_Mix= LnFugCoef_Mix( &
    & 1,T,XC,XW,                 &
    & bH2O,cH2O,dH2O,eH2O,       &
    & bMix,cMix,dMix,eMix,       &
    & CIJ,DIJ,EIJ,vMix,zMix)
    !-------------------------------------------fug.coeff CO2 in mixture
    LnPhiCO2_Mix= LnFugCoef_Mix( &
    & 2,T,XC,XW,                 &
    & bCO2,cCO2,dCO2,eCO2,       &
    & bMix,cMix,dMix,eMix,       &
    & CIJ,DIJ,EIJ,vMix,zMix)
    !
    Vcm3= vMix
    ! e.g. Michelsen-Mollerup eq.119
    LnActH2O=  log(XW) +LnPhiH2O_Mix -LnPhiH2O_Pur
    LnActCO2=  log(XC) +LnPhiCO2_Mix -LnPhiCO2_Pur
  else
    LnActH2O=  log(XW)
    LnActCO2=  log(XC)
    Vcm3=      XC*V_CO2_Pur + XW*V_H2O_Pur
    zMix =     1.0D0
  end if
  !------------------------------//Activities, Fugacities in the mixture
  
  V_m3= Vcm3*1.D-6
  
  ! Z= P*V/R/T
  
  return
end subroutine EoS_H2O_CO2_KerJac

subroutine EoS_KerJac_PhiPur( &
& iCod,TdgK,Pbar, &
& Vcm3,LnPhi)
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
!-- dCO2,dH2O,dMix   attractive term in MRK equation; BAR*(CC**3)*sqrt(T)/MOLE**3
!-- eCO2,eH2O,eMix   attractive term in MRK equation; BAR*(CC**4)*sqrt(T)/MOLE**4
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
  real(dp),intent(out):: Vcm3,LnPhi
  !
  real(dp):: P,T,T2
  real(dp):: Z,V
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
    
    LnPhi= LnFugCoef_Pure(bH2O,cH2O,dH2O,eH2O,T,Z,V) +log(Pbar)

  case(2)
    bCO2=  58.0D0   !Covolume Of CO2
    cCO2= (28.31D0     +0.10721D0*T -0.00000881D0*T2)*1.0D6   ! cCO2
    dCO2= (9380.0D0    -8.53D0   *T +0.001189D0  *T2)*1.0D6   ! dCO2
    eCO2= (-368654.0D0 +715.9D0  *T +0.1534D0    *T2)*1.0D6   ! eCO2
    
    call ZPur( &
    & 2,bCO2,cCO2,dCO2,eCO2,T,P,& !In
    & Z,V)       !Out
    
    LnPhi= LnFugCoef_Pure(bCO2,cCO2,dCO2,eCO2,T,Z,V)

  end select
  
  Vcm3=    V

  return
end subroutine EoS_KerJac_PhiPur

subroutine EoS_H2O_KerJac_PhiPur( &
& TdgK,Pbar, &
& Vcm3,LnPhi)
!--
  real(dp),intent(in) :: Pbar,TdgK
  real(dp),intent(out):: Vcm3,LnPhi
  !
  real(dp):: T,T2
  real(dp):: TC,PK,VI
  real(dp):: Z,V
  real(dp):: bH2O,cH2O,dH2O,eH2O
  
  if(iDebug>2) write(fTrc,'(/,A)') &
  & "---------------------------------------------EoS_H2O_KerJac_PhiPur"
  
  T=  TdgK !Kelvin
  T2= T*T
  !
  bH2O=  29.0D0
  cH2O= (290.78D0    -0.30276D0*T +0.00014774D0*T2)*1.0D6
  dH2O= (-8374.0D0   +19.437D0 *T -0.008148D0  *T2)*1.0D6
  eH2O= (76600.0D0   -133.9D0  *T +0.1071D0    *T2)*1.0D6
  !
  !-------------------------------compute initial guess volume for Solve
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
  !--------------------------------------------------------------------/
  
  !!! if(iDebug>2) write(fTrc,'(5G15.6)') bH2O,cH2O,dH2O,eH2O,VI

  call Solve( &
  & TdgK,Pbar, &
  & bH2O,cH2O,dH2O,eH2O,VI, &
  & Z,V)
  !
  if(iDebug>3) write(fTrc,'(2(A,G15.6,/))') "Z=",Z,"Vcm3=",V
  !
  LnPhi= LnFugCoef_Pure(bH2O,cH2O,dH2O,eH2O,T,Z,V)
  Vcm3=  V

  if(iDebug>2) write(fTrc,'(A,/)') &
  & "--------------------------------------------/EoS_H2O_KerJac_PhiPur"
  
  return
end subroutine EoS_H2O_KerJac_PhiPur

subroutine EoS_CO2_KerJac_PhiPur( &
& TdgK,Pbar, &
& Vcm3,LnPhi)
!--
  real(dp),intent(in) :: Pbar,TdgK
  real(dp),intent(out):: Vcm3,LnPhi
  !
  real(dp):: T,T2
  real(dp):: Z,V
  real(dp):: bCO2,cCO2,dCO2,eCO2
  real(dp):: TC,PK,VI
  
  if(iDebug>2) write(fTrc,'(/,A)') &
  & "---------------------------------------------EoS_CO2_KerJac_PhiPur"
  
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
  LnPhi= LnFugCoef_Pure(bCO2,cCO2,dCO2,eCO2,T,Z,V)
  Vcm3=V

  if(iDebug>2) write(fTrc,'(A,/)') &
  & "--------------------------------------------/EoS_CO2_KerJac_PhiPur"
  
  return
end subroutine EoS_CO2_KerJac_PhiPur

subroutine ZPur( & !
& iCod,          & !in
& B,C,D,E,       & !in
& T,P,           & !in
& Z,Vcm3)          !out
!--
!-- compute the volume of CO2 and H2O at T,P
!-- an initial guess of the volume (vi) is chosen for a given pressure range.
!-- this value is then used in the routine Solve,
!-- which solves for the exact volume by means of an iterative newton - raphson technique.
!--
  integer, intent(in) :: iCod ! 1=H2O 2=CO2
  real(dp),intent(in) :: T,P  ! P=bar, T=K
  real(dp),intent(in) :: B,C,D,E
  real(dp),intent(out):: Z,Vcm3 !zCO2,zH2O,vCO2,vH2O
  !---------------------------------------------------------------------
  real(dp):: PK,TC
  real(dp):: VI    ! V initial guess
  !-------------------------------compute initial guess volume for Solve
  TC=   T -T_CK     !-> Celsius
  PK=   P /1000.0D0 !-> kilobar
  select case(iCod)
  case(1) ! H2O
    if (PK>=One)                                      VI=  22.0D0
    if (PK>=0.90D0 .and. PK<One)                      VI=  24.2D0
    if (PK>=0.60D0 .and. PK<0.9D0)                    VI=  31.2D0
    if (PK>=0.21D0 .and. PK<0.60D0 .and. TC>=550.0D0) VI=  75.0D0
    if (PK>=0.21D0 .and. PK<0.60D0 .and. TC< 550.0D0) VI=  35.0D0
    if (PK>=0.10D0 .and. PK<0.21D0 .and. TC< 400.0D0) VI=  15.0D0
    if (PK>=0.10D0 .and. PK<0.21D0 .and. TC>=400.0D0) VI= 100.0D0
    if (PK>=0.005D0.and. PK<0.10D0)                   VI= 500.0D0
    if (PK< 0.005D0)                                  VI=1000.0D0
  case(2) ! CO2
    if (PK >= One)                        VI=  35.0D0
    if (PK >= 0.10D0  .and. PK < One)     VI= 100.0D0
    if (PK >= 0.005D0 .and. PK < 0.10D0)  VI= 500.0D0
    if (PK <  0.005D0)                    VI=5000.0D0
  end select
  !--------------------------------------------------------------------/
  call Solve( &
  & T,P, &
  & B,C,D,E,VI, &
  & Z,Vcm3)
  !
  if(iDebug>3) write(fTrc,'(2(A,G15.6,/))') "Z=",Z,"Vcm3=",Vcm3
  !
  return
end subroutine ZPur

subroutine ZMixCalc(   &
& T,P,                 &
& XC,XW,vCO2,vH2O,     &
& bMix,cMix,dMix,eMix, &
& Vcm3,Z)
!--
!-- calculate the volume of a mixture of CO2 and H2O at P,T,Xi
!--
!-- the molar volumes of CO2 and H2O as calculated in ZPur
!-- are used to define the initial estimate
!--
!-- routine Solve is then used to calculate the volume of the mixture.
!--
  real(dp),intent(in) :: P,T,XC,XW,vCO2,vH2O,bMix,cMix,dMix,eMix
  real(dp),intent(out):: Vcm3,Z
  real(dp):: VInit !,Z,V

  VInit= vCO2*XC + vH2O*XW
  ! VInit is VIdeal = Sum(X_i*V_i)

  call Solve( &
  & T,P,&
  & bMix,cMix,dMix,eMix,VInit,&
  & Z,Vcm3)

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
!-- volume unit = cm3
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
  integer :: K

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
    PA2=  PA1/sqrt(T)/VI/BI
    F=    PR - PA2 - P

    ! DEFINITION OF THE DIFFERENTIAL OF F(X) FOR NEWARP:
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
    DPA=  (PA1*D7+D8*D9) / sqrt(T)
    DF=   DPR - DPA

    ! CALCULATION OF V(K+1)
    V=     VI - F/DF
    DifF=  ABS(V-VI)
    
    ! CONTINUATION OR end OF ITERATIONS FOR NEWARP
    if (DifF<0.01D0) exit

    if (V>1.0D6) V=1.0D6
    if (V<9.9D0) V=1.0D1

    VI= V

  end do

  Z=  V*P /R_CbK/T

  return
end subroutine Solve

real(dp) function LnFugCoef_Pure(B,C,D,E,T,Z,V) !MODIFIED
!--
!-- calculate fugacity coefficient of pure co2 or pure h2o at T,P
!-- cf kerrick and jacobs
!-- for a derivation of the pure fugacity coefficient expression (fcp)
!--
  real(dp),intent(in):: T
  real(dp),intent(in):: B,C,D,E
  real(dp),intent(in):: Z,V     !V is in cm3
  !
  real(dp):: RT,Y,FCP
  !
  RT= R_CbK*T *sqrt(T)

  !select case(iCod)
  !  case(1); B=bH2O; C=cH2O; D=dH2O; E=eH2O !H2O
  !  case(2); B=bCO2; C=cCO2; D=dCO2; E=eCO2 !CO2
  !end select

  Y= B/(4.0D0*V)

  FCP= (8.0D0*Y -9.0D0*Y*Y +3.0D0*Y**3)       &
  &    / ((One-Y)**3)                         &
  &   - log(Z)                                &
  &   - C/(RT    *(V+B))    - D/(RT*V*(V+B))  &
  &   - E/(RT*V*V*(V+B))    + (C/(RT*B))*(log(V/(V+B)))    &
  &   - D/(RT*B*V)          + (D/(RT*B*B))*(log((V+B)/V))  &
  &   - E/(RT*2.0D0*B*V*V)  + E/(RT*B*B*V)    &
  &   - (E/(RT*B**3))*(log((V+B)/V))

  LnFugCoef_Pure= FCP

  return
end function LnFugCoef_Pure

real(dp) function LnFugCoef_Mix(&
& iCod,                &
& T,XC,XW,             &
& bPur,cPur,dPur,ePur, &
& bMix,cMix,dMix,eMix, &
& CIJ,DIJ,EIJ,         &
& vMix,zMix)
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
  real(dp),intent(in) :: T
  real(dp),intent(in) :: XC,XW
  real(dp),intent(in) :: bPur,cPur,dPur,ePur
  real(dp),intent(in) :: bMix,cMix,dMix,eMix
  real(dp),intent(in) :: CIJ,DIJ,EIJ
  real(dp),intent(in) :: vMix,zMix

  real(dp):: RT
  real(dp):: B,C,D,E
  real(dp):: B1,C1,D1,E1
  real(dp):: X1,X2,FCM,VBV,VB,V,Y,Z

  RT= R_CbK*T*sqrt(T)
  B=  bMix; C=  cMix; D=  dMix; E=  eMix
  V=  vMix; Z=  zMix
  Y=  B/V/4.0D0
  
  B1=bPur; C1=cPur; D1=dPur; E1=ePur
  
  !select case(iCod)
  !  case(1) ; B1=bH2O; C1=cH2O; D1=dH2O; E1=eH2O; X1=XW; X2=XC !H2O
  !  case(2) ; B1=bCO2; C1=cCO2; D1=dCO2; E1=eCO2; X1=XC; X2=XW !CO2
  !end select
  
  select case(iCod)
  case(1) ; X1=XW  ; X2=XC !H2O
  case(2) ; X1=XC  ; X2=XW !CO2
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

  LnFugCoef_Mix= FCM

  !if (J==1) PhiH2O_Mix=  FCM
  !if (J==2) PhiCO2_Mix=  FCM

  return
end function LnFugCoef_Mix

end module M_Eos_KerJac

!!----------------------------------------------------------------------
!! this routine is included in Zpur
!!----------------------------------------------------------------------
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

!!----------------------------------------------------------------------
!! this routine is included in Zpur
!!----------------------------------------------------------------------
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


