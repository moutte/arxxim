MODULE M_T_DtbH2OHkf
!--
!-- implements TYPE T_DtbH2OHkf,
!-- container of Solvent properties of water for the HKF EoS,
!-- and related routines
!--
!-- tools for storing/computing Solvent parameters (- density, dielectric, etc.) fot HKF EoS
!-- routines retrieved from SupCrt92, sometime with help of Theriak code
!
!-- main issue for application of HKF EoS is derivation of
!-- solvent (-H2O) properties:
!--   density:    computed by RhoEtc
!--   g-FUNCTION: computed by GShock9
!--   dielectric: computed by JohnsonNorton91
!
!-- as by-product: compute A and B debye hueckel parameter
!-- for given (P,T) conditions
!
!-- NB: SupCrt tools for checking validity of P,T domain not re-implemented !!!
!--
  
  USE M_Kinds
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC:: T_DtbH2OHkf
  PUBLIC:: DtbH2OHkf_Calc
  
  !DATA satmin / 0.01d0, 0.006117316772d0 /
  !DATA aa, ba, ca / 0.549824d3,  0.65995d0, -0.4973d-4 /
  !DATA VPtTta, VPrTtb, Stran / 0.23348d2, 0.2372d2, 0.342d0 /
  
  REAL(dp),PARAMETER,PUBLIC:: &
  ! limits of application of HKF model
  & TMax= 1000._dp, & ! TdgC
  & PMax= 5000._dp    ! Pbar
  
  REAL(dp),PARAMETER,PUBLIC:: &
  ! ZPrTr,YPrTr from SupCrt:
  ! "assignment of [Z,Y]PrTr to Johnson-Norton (1991) values
  !  for distribution version"
  & ZPrTr= -0.1278034682D-1,& !value of Z (Born coeff') at Pref,Tref
  & YPrTr= -0.5798650444D-4   !value of Y (Born coeff') at Pref,Tref
   
  REAL(dp),PARAMETER,PUBLIC:: &
  ! parameters for HKF model for aqu.species:
  & Eta=  166027._dp,&  !=Avogadro*electron^2/2
  !                     !=6.022E23*(4.803E-10)^2 /2
  !                     !=1.66027E5 Angst.Cal/mole !JOH92,p905
  & Theta=228._dp,   &  != 228 K
  ! Theta= experimental constant for non solv. term in EoS of aqu. species
  & Psi=  2600._dp,  &  !=2600 bars
  ! Psi=  experimental constant for non solv. term in EoS of aqu. species 
  & gShockRef= Zero     !value gShock at Pr,Tr
  ! anion= Zero
  ! cation= 0.94d0
  ! can be used for a0 term in D-H  (rx= crystallographic radius)
  !   for cations, r(j)= rx(j) + Z(j).(0.94 + g)
  !   for anions,  r(j)= rx(j) + Z(j). g) 
    
  TYPE:: T_DtbH2OHkf ! H2O Properties
  ! contains, for a given (Pbar,TdgK),
  ! all DATA on H2O necessary for application of the HKF EoS for aqueous species,
  ! also produces values of a and b Debye-Hueckel parameters
    REAL(dp)::&
    & Pbar,TdgK,&
    & Rho, &
    & Alfa,dAlfdT, &                    !thermal expansion.
    & Beta,dBetdT,dBetdP,&              !compressibility
    & GShok,dGshdP,dGshdT,d2GshdT2,&    !G Shock FUNCTION
    & Eps,dEpsdP,dEpsdT,d2EpsdT2,&      !dielectric
    & Q,X,Y,Z, &                        !Born coeff's
    & dhA,dhB                           !DebyeHueckelParam's
  ENDTYPE T_DtbH2OHkf 

CONTAINS

SUBROUTINE DtbH2OHkf_Calc(TdgK,Pbar,pW) 
!--
!-- computes H2O properties (saved in variable pW, of TYPE T_DtbH2OHkf,
!-- for a given (Pbar,TdgK) point
!-- calls SUBROUTINEs CalcRhoW, CalcGSHOK, CalcEpsJN91 adapted from SupCrt92
!-- (routines implemented outside the module, could be nested inside,
!--  as they are of no use outside)
!--

  USE M_Fluid_Calc,ONLY: Eos_H2O_Rho
  USE M_Trace,     ONLY: iDebug,fTrc
  
  REAL(dp),         INTENT(IN) :: TdgK,Pbar
  TYPE(T_DtbH2OHkf),INTENT(OUT):: pW
  
  pW%TdgK=TdgK
  pW%Pbar=Pbar
  
  CALL Eos_H2O_Rho( & !-> density in g/cm3 !!
  & TdgK,Pbar,& !IN
  & pW%Rho, pW%Alfa, pW%Beta, pW%dAlfdT, pW%dBetdT, pW%dBetdP) !OUT
  
  CALL CalcGSHOK( &
  & TdgK,Pbar,pW%Rho,pW%Beta,  pW%Alfa,  pW%dAlfdT, & !IN
  & pW%Gshok, pW%dGShdP, pW%dGShdT, pW%d2GShdT2)      !OUT
  
  CALL CalcEpsJN91( &
  & TdgK, pW%Rho, pW%Beta, pW%Alfa, pW%dAlfdT, & !IN
  & pW%Eps, pW%dEpsdP, pW%dEpsdT, pW%d2EpsdT2)   !OUT
  
  !CALL CalcBorn92( &
  !& pW%Eps,pW%dEpsdP,pW%dEpsdT,pW%d2EpsdT2, & !IN
  !& pW%Z_, pW%Q_, pW%Y_, pW%X_)               !OUT
  
  pW%Z=   One/pW%Eps -One !not the original definition of Z;=Zer_= - Zeps -1
  pW%Y=   pW%dEpsdT /pW%Eps /pW%Eps
  pW%Q=   pW%dEpsdP /pW%Eps /pW%Eps
  pW%X=   (pW%d2EpsdT2 -2.0D0 /pW%Eps *pW%dEpsdT*pW%dEpsdT) /pW%Eps /pW%Eps
  
  ! formulas for pW%Rho in g/cm3 !!
  ! with Rho in kg/m3, should use Rho/1.D3 !!
  pW%dhA= 1.824829238D6 *SQRT(pW%Rho) /SQRT((pW%Eps*TdgK)**3)
  !with IonSizeParam in NM, dhA=0.5095 at 25°C
  pW%dhB= 50.29158649D0 *SQRT(pW%Rho) /SQRT(pW%Eps*TdgK)
  !with IonSizeParam in NM, dhB=0.3284 at 25°C
  
ENDSUBROUTINE DtbH2OHkf_Calc

SUBROUTINE CalcGShok( & !SupCrt92
& TdgK,Pbar,Dgcm3,beta,alpha,daldT, &
& g,dgdP,dgdT,d2gdT2)
!--
!-- Computes g, dgdP, dgdT, and d2gdT2
!-- using equations given by Shock et al. (1992) 
!-- g is a P-T dependent Solvent function (Tanger & Helgeson, 1988)
!-- g and its derivatives are necessary for calculation of
!--   Re(j), effective electrostatic radius of species j, 
!--   W(j), conventional Born coefficient
!--

  USE M_Dtb_Const,ONLY:T_CK
  
  REAL(dp),INTENT(IN) :: TdgK,Pbar,Dgcm3,beta,alpha,daldT
  REAL(dp),INTENT(OUT):: g,dgdP,dgdT,d2gdT2
  
  !beta,dgdP: /bar; alpha,dgdT: /K; daldT,d2gdT2: /K^2
  REAL(dp)::&
  & T,P,D,&
  & a, b,&
  & dgdD,dgdD2,&
  & dadT,dadTT,&
  & dbdT,dbdTT,&
  & dDdT,dDdP,dDdTT,&
  & Db,dDbdT,dDbdTT,&
  & f,ft,dftdT,dftdTT,&
  & fp,dfpdP,&
  & dfdP,dfdT,d2fdT2
  
  REAL(dp),DIMENSION(6):: c= &
  (/ -0.2037662D+01,  0.5747000D-02, -0.6557892D-05, &
  &   0.6107361D+01, -0.1074377D-01,  0.1268348D-04  /)
  REAL(dp),DIMENSION(3):: cc= &
  (/  0.3666666D+02, -0.1504956D-9,   0.5017997D-13 /)
  
  T=TdgK-T_CK; P=Pbar; D=Dgcm3
  
  IF (D >= One) RETURN
  
  a= c(1) + c(2)*T + c(3)*T**2
  b= c(4) + c(5)*T + c(6)*T**2
  g= a*(One - D)**b
  
  dgdD= - a *b *(One - D)**(b - One)
  dgdD2=  a *b *(b - One)  *(One - D)**(b - Two)
  
  dadT=   c(2) + c(3)*Two*T
  dadTT=         c(3)*Two
  
  dbdT=   c(5) + c(6)*Two*T
  dbdTT=         c(6)*Two
  
  dDdT= - D *alpha
  dDdP=   D *beta
  dDdTT=- D *(daldT - alpha**2)
  
  Db=    (One - D) **b
  dDbdT= -b *(One - D)**(b-One) *dDdT + LOG(One - D) *Db  *dbdT
  dDbdTT= -(b * (One - D)**(b-One) * dDdTT &
        + (One - D)**(b-One) * dDdT * dbdT + b * dDdT * (-(b-One) * (One - D)**(b-Two) * dDdT &
        + LOG(One - D) * (One - D)**(b-One) * dbdT)) &
        + LOG(One - D) * (One - D)**b * dbdTT - (One - D)**b * dbdT * dDdT / (One - D) &
        + LOG(One - D) * dbdT * dDbdT
  dgdP= dgdD * dDdP
  dgdT= a*dDbdT + Db*dadT
  d2gdT2= a*dDbdTT + dDbdT*dadT*Two + Db*dadTT
  IF ((T < 155.0D0) .OR. (P > 1000.0D0) .OR. (T > 355.0D0)) RETURN
  !
  ft=       ((T - 155.0D0)/300.0D0)**4.8D0 &
  & + cc(1)*((T - 155.0D0)/300.0D0)**16
  !
  dftdT= 4.8D0 /300.0D0       *((T - 155.0D0)/300.0D0)**3.8D0 &
  &    + 16.0D0/300.0D0 *cc(1)*((T - 155.0D0)/300.0D0)**15
  
  dftdTT= 3.8D0 *4.8D0 /300.0D0**2       *((T - 155.0D0)/300.0D0)**2.8D0 &
  &     + 15.0D0*16.0D0/300.0D0**2 *cc(1)*((T - 155.0D0)/300.0D0)**14
  !
  fp= cc(2)*(1000.0D0 - P)**3 &
  & + cc(3)*(1000.0D0 - P)**4
  !
  dfpdP= -3.0D0*cc(2)*(1000.0D0 - P)**2 &
  &    - 4.0D0 *cc(3)*(1000.0D0 - P)**3 
  !
  f=      ft *fp
  dfdP=   ft *dfpdP
  dfdT=   fp *dftdT
  d2fdT2= fp *dftdTT
  g=      g - f
  dgdP=   dgdP   - dfdP
  dgdT=   dgdT   - dfdT
  d2gdT2= d2gdT2 - d2fdT2
  
  RETURN
ENDSUBROUTINE CalcGShok

!! SUBROUTINE CalcBorn92(Eps,dEps_dP,dEps_dT,d2Eps_dT2,Zeps,Qeps,Yeps,Xeps)
!! !                     >InPut________________________>OutPut
!! !Compute the Z, Q, Y, and X Born FUNCTIONs from their eps, dedP, dedT, and d2edT2 counterparts.
!! !produced by Johnson-Norton (1991) equation (= JN91)
!! !from SupCrt92
!! !(not used, rewritten within calling procedures)
!!   IMPLICIT NONE !REAL(dp) (a-h,o-z)
!!   REAL(dp)::&
!!     Eps,dEps_dP,dEps_dT,d2Eps_dT2,& !IN
!!     Zeps,Qeps,Yeps,Xeps             !OUT
!!   !
!!   Qeps= dEps_dP/Eps/Eps
!!   Xeps= (d2Eps_dT2-2.0D0/Eps*dEps_dT**2)/Eps/Eps
!!   Yeps= dEps_dT/Eps/Eps
!!   Zeps=-One/Eps
!! ENDSUBROUTINE CalcBorn92

SUBROUTINE CalcEpsJN91( &
& TdgK,Dgcm3,Beta,Alpha,dAldT, &
& Eps,dEdP,dEdT,d2EdT2)
!--
!-- from SupCrt (Johnson & Norton 91), also found in deCapitani
!-- Compute (eps, dedP, dedT, d2edT2)(T,D) 
!-- using equations given by Johnson and Norton (1991)
!-- -- fit parameters regressed from least squares fit to dielectric DATA
!-- consistent with the HK74 equation at low temperatures
!-- and with the Pitz83 equation at high temperatures
!--
!-- Units:
!--   T: K; D: g/cm3;
!--   beta,dedP: /bar;
!--   alpha,dedT: /K;
!--   daldT,d2edT2: /K^2
!--

  USE M_Dtb_Const, ONLY:TRef
  
  REAL(dp),INTENT(IN) :: TdgK,Dgcm3,Beta,Alpha,dAldT
  REAL(dp),INTENT(OUT):: Eps,dEdP,dEdT,d2EdT2
  
  REAL(dp):: a(10)
  REAL(dp):: c(5), dcdT(5), dc2dTT(5)
  REAL(dp):: T, D, Tn
  INTEGER :: J,K
  !SAVE
  !REAL(dp),PARAMETER:: Tref=298.15_R8_
  
  T=TdgK
  D=Dgcm3
  
  a=(/ &
  &  0.1470333593E+02, 0.2128462733E+03,-0.1154445173E+03, &
  &  0.1955210915E+02,-0.8330347980E+02, 0.3213240048E+02, &
  & -0.6694098645E+01,-0.3786202045E+02, 0.6887359646E+02, &
  & -0.2729401652E+02 /)
  
  Tn= T / Tref
  c(1)= One
  dCdT(1)= Zero
  dC2dTT(1)= Zero
  
  C(2)= a(1)/Tn
  dCdT(2)= -a(1)*Tref/T**2
  dC2dTT(2)= Two*a(1)*Tref/T**3
  
  C(3)= a(2)/Tn &
      + a(3) &
      + a(4)*Tn           
  dCdT(3)= -a(2)*Tref/T**2 &
         +  a(4)/Tref           
  dC2dTT(3)= Two*a(2)*Tref/T**3
  
  C(4)= a(5)/Tn &
      + a(6)*Tn &
      + a(7)*Tn**2
  dCdT(4)= -a(5)    *Tref/T**2 &
         +  a(6)         /Tref &
         +  a(7)*Two*T   /Tref**2
  dC2dTT(4)= Two*a(5)*Tref/T**3 &
           + Two*a(7)/Tref**2
  
  C(5)= a(8)/Tn**2 &
      + a(9)/Tn &
      + a(10)
  dCdT(5)= -Two*a(8)*Tref**2/T**3 &
         - a(9)*Tref/T**2
  dC2dTT(5)= 6.0d0*a(8)*Tref**2/T**4 &
           + Two*a(9)*Tref/T**3
  
  eps= Zero
  DO k=1,5
    eps= eps  + c(k)*D**(k-1)
  ENDDO !Eps=C1+C2*D^2+..+C5*D^4
  
  dedP= Zero
  DO j= 0,4
    dedP= dedP + j*c(j+1)*D**j
  ENDDO
  dedP= beta * dedP
  
  dedT= Zero
  DO j= 0,4
    dedT= dedT &
    &   + D**j *( dcdT(j+1) - j*alpha*c(j+1) )
  ENDDO
  
  d2edT2= Zero
  DO j= 0,4
    d2edT2= d2edT2 &
    & + D**j &
    & * (dc2dTT(j+1) &
    &   -j *(alpha*dcdT(j+1) + c(j+1)*daldT) &
    &   -j *alpha *(dcdT(j+1) - j*alpha*c(j+1)))
  ENDDO
  
  RETURN
ENDSUBROUTINE CalcEpsJN91

ENDMODULE M_T_DtbH2OHkf

!SUBROUTINE Omeg92(& 
!  WRef_,Chg_,aname_,&
!  g,dgdP,dgdT,d2gdT2,&
!  W,dWdP,dWdT,d2WdT2)
!! from, SupCrt92 !not used, rewritten within calling procedure
!! 
!! conventional Born coefficient (W) of the current aqueous species, 
!! and its derivatives, dwdP, dwdP, and dw2dT2
!! as a FUNCTION of g, dgdP, dgdT, d2gdT2, wref, and Charge 
!! using equations given by Johnson et al. (1991)
!!
!!IN=
!!  G-FUNCTION and its derivatives= g,dgdP,dgdT,d2gdT2
!!  Species properties= WRef_,Chg_,aname_
!!OUT=
!!  Omega and its derivatives= W,dWdP,dWdT,d2WdT2
!!
!  USE M_T_DtbH2OHkf, ONLY: eta,gref
!  IMPLICIT NONE 
!  CHARACTER*20,INTENT(IN) :: aname_
!  REAL(dp),   INTENT(IN) :: WRef_,Chg_, g,dgdP,dgdT,d2gdT2
!  REAL(dp),   INTENT(OUT)::             w,dwdP,dwdT,d2wdT2 
!  REAL(dp)::reref, re, sqZ, Z3, Z4  
!  !
!  IF ((Chg_== Zero) .OR. (aname_== 'H+')) THEN !NEUTRAL AQUEOUS SPECIES OR H+ 
!    w= wref_
!    dwdP= Zero
!    dwdT= Zero
!    d2wdT2= Zero
!    RETURN
!  ELSE !__________________________________________!CHARGED AQUEOUS SPECIES OTHER THAN H+
!    sqZ= Chg_*Chg_
!    reref= sqZ / (wref_/eta + Chg_/(3.082d0 + gref))
!    re= reref + ABS(Chg_) * g
!    w= eta*(sqZ/re   - Chg_/(3.082d0 + g))
!    Z3= ABS(Chg_)*sqZ/re**2 - Chg_/(3.082d0 + g)**2
!    Z4= sqZ*sqZ/re**3 - Chg_/(3.082d0 + g)**3
!    dwdP= -eta * Z3 * dgdP
!    dwdT= -eta * Z3 * dgdT
!    d2wdT2= 2.0d0 * eta * Z4 * dgdT**2 - eta * Z3 * d2gdT2
!  END IF
!ENDSUBROUTINE Omeg92

