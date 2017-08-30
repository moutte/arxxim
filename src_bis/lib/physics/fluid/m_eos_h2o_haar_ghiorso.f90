module M_Eos_H2O_Haar_Ghiorso
!! was M_Haar_Ghiorso_H2O

!-----------------------------------------------------------------------
! HAAR.GHIORSO : Eos_H2O_Haar_Ghiorso
!-----------------------------------------------------------------------
! Purpose : Haar and al EOS to compute Properties of H2O
!           Gibbs energy, Enthalpy, Entropy and Molar Volume
!
! Original : L.Mastin [ haarH2O.f ]
! http://volcanomodels.sr.unh.edu/workshops/workshop-2002/module_library/haarH2O.f
!
!-----------------------------------------------------------------------
! Arxim Integration : J.Moutte
!-----------------------------------------------------------------------

use M_Kinds
implicit none

private

!// public functionS
public::  Eos_H2O_Haar_Ghiorso !!CalcGH2O_Haar_Ghiorso

contains

subroutine Eos_H2O_Haar_Ghiorso ( &
& TdgK,Pbar, &
& gH2Ort,hH2O,sH2O,vH2O_m3 )
  !
  use M_Dtb_Const,only: R_jK, DtbConv_Benson,S0_Hydrogen,S0_Oxygen,Tref
  !
  real(dp),intent(in)  :: TdgK,Pbar
  real(dp),intent(out) :: gH2Ort,hH2O,sH2O,vH2O_m3
  !
  real(dp) :: gH2O,cp,cv,rhof,dvdth2o,blkgas,sndspeed
  !
  call CalcGH2O_Haar_Ghiorso_Detail &
  & (TdgK,Pbar, &
  &  gH2O,hH2O,sH2O,vH2O_m3,&
  &  rhof,cp,cv,blkgas,dvdth2o,sndspeed)
  !
  !Robie/COdata data for entropies of elements at 298.15:
  !! H2: 130.68 J/K/mole = 2*S0_Electron
  !! O2: 205.15 J/K/mole
  !! C:    5.74 J/K/mole
  !! -> for H2O, to transform G(BermBrown) to G(BensHelg),
  !!____add 298.15*(130.68 + 205.15/2)=  69544.98 J
  !
  if(DtbConv_Benson) gH2O= gH2O + Tref *(S0_Hydrogen *Two + S0_Oxygen)
  !
  gH2Ort= GH2O /R_jK /TdgK
  !
end subroutine Eos_H2O_Haar_Ghiorso

!---

subroutine CalcGH2O_Haar_Ghiorso_Detail &
& (TdgK,Pbar, &
&  gH2O,hH2O,sH2O,vH2O_m3,&
&  rhof,cp,cv,blkgas,dvdth2o,sndspeed)
!--
!-- CALCULATES THERMODYNAMIC PROPERTIES OF WATER AND STEAM
!-- USING THE EQUATIONS OF HAAR ET AL. (1984).  
!--
!-- program AUTHORS:
!--  (1) MARK GHIORSO, DEPT. OF GEOLOGY & GEOPHYSICS UNIVERSITY OF WASHINGTON
!--  Seattle, WA 98195-1310
!--  ghiorso@u.washington.edu 
!--  (2) LARRY G. MASTIN, U.S. GEOlogical SURVEY CASCADES VOLCANO OBSERVATORY
!--  lgmastin@usgs.gov
!--  tel. 360-993-8925 (USA)
!-- http://vulcan.wr.usgs.gov/Projects/Mastin
!-- DATE:        WRITTEN, 1998, REVISED DECEMBER 2002
!-- LANGUAGE:    FORTRAN 90 (98% fortran 77)
!-- THIS program WAS ORIGINALLY  WRITTEN BY MARK GHIORSO IN C.
!-- IT WAS TRANSLATED INTO FORTRAN BY LARRY MASTIN, JUNE 1998,
!-- AND doNATED TO THE UNIVERSITY OF NEW HAMPSHIRE CONDUIT MODELING WORKSHOP WEB PAGE DECEMBER, 2002.
!-- 
!-- INPUT parameterS:
!--   TdgK temperature, Kelvin
!--   Pbar pressure, bar
!-- 
!-- OUTPUT parameterS:
!--   cp       specific heat at constant pressure (J/kg K)
!--   cv       specific heat at constant volume (J/kg K)
!--   gH2O     Gibbs free energy (J/kg)
!--   hH2O     enthalpy relative to the elements (H and O) at STP (J/kg)
!--   sH2O     entropy relative to the elements at STP (J/kg K)
!--   rhof     density (kg/m3)
!--   vH2Ocm3  molar volume (cm3/mole)
!--   blkgas   bulk modulus (Pascals)
!--   dvdth2o  change in molar volume with temperature (m3/kg K)
!--   sndspeed sound speed (m/s)
!--
  use M_Eos_H2O_Rho_Psat, only: Eos_H2O_psat
  use M_Eos_Utils, only: kubik, f_Square, f_Cube, f_Quart, f_Quint
  !
  real(dp),intent(in)  :: TdgK,Pbar
  real(dp),intent(out) :: gH2O,hH2O,sH2O,vH2O_m3
  real(dp),intent(out) :: cp,cv,rhof,dvdth2o,blkgas,sndspeed 
  !
  real(dp):: vH2Ocm3
  real(dp):: &
  r=  4.6152D0, & 
  !4.6152= 10 * 8.31441 / 18.0152 = R in bar.cc/ gramH2O /K
  rr= 8.31441D0
  !Nota: in RobieHemingway (from COdata), R_jK= 8.314510D0 J/K/Mole
  !
  real(dp):: &
  gref=  -54955.23970014D0, &
  !-54955.2356146121147 25 c and 1 bar 
  href=  -34099.89230644D0
  !! sref=  69.94917790942D0  !not used ??
  !! 69.9146D0 is used instead ??
  !
  real(dp):: &
  t0= 647.25D0, &
  ps= 220.55D0, &
  P0= 1.01325D0 != 1 atm
  !
  real(dp):: T,P
  !
  real(dp):: alpha,beta,gamma
  real(dp):: &
  & taui(0:6), ermi(0:9), &
  & vol, rhn, &
  & dp_, dr_, rh, pr, &
  & dpr, q10, qm, x_, tr, dtrdt, dpdrh, dpdt, drhdt, d2pdt2, &
  & d2pdtdrh, d2pdrh2, temp
  !
  real(dp):: ark, brk, oft, buk, cuk, duk 
  real(dp):: x1,x2,x2i,x3
  integer :: i, icount
  !
  !TERMS useD IN BASE function
  real(dp):: b,  dbdt,  d2bdt2,  d3bdt3
  real(dp):: bb, dbbdt, d2bbdt2, d3bbdt3
  real(dp):: &
  & y, dydt, dydrh, d2ydt2, d2ydtdrh, d2ydrh2, &
  & d3ydt3, d3ydt2drh, d3ydtdrh2, d3ydrh3, &
  & yy,yy2,yy3,yy4
  !
  !BASE function AND ITS DERIVATIVES
  real(dp):: &
  & Z, dZdt, dZdrh, d2Zdt2, d2Zdtdrh, d2Zdrh2,& 
  & d3Zdt3, d3Zdt2drh,d3Zdtdrh2, d3Zdrh3
  real(dp):: &
  & Ab, dAbdt, dAbdrh, d2Abdt2, d2Abdtdrh, &
  & d2Abdrh2, d3Abdt3, d3Abdt2drh, d3Abdtdrh2, d3Abdrh3
  !
  !TERMS useD IN RESIDUAL function
  real(dp):: del, tau, abc
  !
  !RESIDUAL function AND ITS DERIVATIVES
  real(dp):: &
  & Ar, dArdt, dArdrh, d2Ardt2, d2Ardtdrh, &
  & d2Ardrh2, d3Ardt3, d3Ardt2drh, d3Ardtdrh2, d3Ardrh3
  !
  !IDEAL GAS function AND ITS DERIVATIVES
  real(dp):: Ai, dAidt, d2Aidt2, d3Aidt3
  !
  !HELMHOLTZ function AND ITS DERIVATIVES
  real(dp):: &
  & A, dAdt, dAdrh, d2Adt2, d2Adtdrh, d2Adrh2, &
  & d3Adt3, d3Adt2drh, d3Adtdrh2, d3Adrh3
  !
  !final THERMODYNAMIC valueS (besides those as intent(out))
  real(dp):: &
  & cpH2O,cvH2O,dcpdtH2O, &
  & dvdpH2O,d2vdt2H2O,&
  & d2vdtdpH2O,d2vdp2H2O
  !
  !OTHER TERMS
  real(dp):: &
  & eps1,expQ,Q,dQdt,dQdrh,d2Qdtdrh,&
  & d2Qdrh2,d2Qdt2,d3Qdt3,&
  & d3Qdt2drh,d3Qdtdrh2,d3Qdrh3
  !
  real(dp):: gi(1:40)=(/ &
  & -.53062968529023D4,  .22744901424408D5,  .78779333020687D4, &
  & -.69830527374994D3,  .17863832875422D6, -.39514731563338D6, &
  &  .33803884280753D6, -.13855050202703D6, -.25637436613260D7, &
  &  .48212575981415D7, -.34183016969660D7,  .12223156417448D7, &
  &  .11797433655832D8, -.21734810110373D8,  .10829952168620D8, &
  & -.25441998064049D7, -.31377774947767D8,  .52911910757704D8, &
  & -.13802577177877D8, -.25109914369001D7,  .46561826115608D8, &
  & -.72752773275387D8,  .41774246148294D7,  .14016358244614D8, &
  & -.31555231392127D8,  .47929666384584D8,  .40912664781209D7, &
  & -.13626369388386D8,  .69625220862664D7, -.10834900096447D8, &
  & -.22722827401688D7,  .38365486000660D7,  .68833257944332D5, &
  &  .21757245522644D6, -.26627944829770D5, -.70730418082074D6, &
  & -.225D1,            -1.68D1,             .055D1,            &
  & -93.0D1 /)                                                  
  !
  real(dp)::ci(1:18)=(/ &
  & .19730271018D2,      .209662681977D2,   -.483429455355D0,   &
  & .605743189245D1,   22.56023885D0,        -9.87532442D0,     &
  &-.43135538513D1,      .458155781D0,        -.47754901883D-1, &
  & .41238460633D-2,    -.27929052852D-3,    .14481695261D-4,   &
  &-.56473658748D-6,     .16200446D-7,      -.3303822796D-9,    &
  & .451916067368D-11,  -.370734122708D-13,  .137546068238D-15 /)
  !
  integer::ki(1:40)=(/ &
  & 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, &
  & 6, 6, 6, 6, 7, 7, 7, 7, 9, 9, 9, 9, 3, 3, 1, 5, 2, 2, 2, 4  /)
  integer:: li(1:40)=(/ &
  & 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, &
  & 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 0, 3, 3, 3, 0, 2, 0, 0  /)
  !
  real(dp):: rhoi(37:40)= (/ 0.319D0, 0.310D0, 0.310D0,  1.55D0 /)
  real(dp):: ttti(37:40)= (/ 640.0D0, 640.0D0, 641.6D0, 270.0D0 /)
  real(dp):: alpi(37:40)= (/  34.0D0,  40.0D0,  30.0D0,1050.0D0 /)
  real(dp):: beti(37:40)= (/   2.0D4,   2.0D4,   4.0D4,  25.0D0 /)
  !
  real(dp)::bi(0:5) =&
  (/ 0.7478629D0,  -0.3540782D0,  Zero,         &
  &  0.007159876D0, Zero,        -0.003528426D0 /)
  real(dp)::bbi(0:5)=&
  (/ 1.1278334D0,  -0.5944001D0, -5.010996D0,   &
  &  Zero,          0.63684256D0, Zero          /)

  !R     SPECifIC GAS CONSTANT
  !GREF  GIBBS FREE ENERGY
  !HREF  ENTHALPY
  !SREF  ENTROPY
  !TO    TEMPERATURE AT CRITICAL POINT
  !P     PRESSURE AT CRITICAL POINT
  !RR    UNIVERSAL GAS CONSTANT (J/MOLE K)
  !ALPHA CONSTANT useD TO CALCULATE BASE function 
  !BETA  CONSTANT useD TO CALCULATE BASE function
  !GAMMA CONSTANT useD TO CALCULATE BASE function
  !PO    PRESSURE AT 1 ATM (IN BARS)
  
  p=  Pbar
  t=  TdgK
  !
  alpha= 11.0D0
  beta=  133.0D0 /3.0D0
  gamma= 7.0D0   /2.0D0
  !
  !The values of reduced temperature (T/T0)**i are stored in the  array TAUI(i)
  taui(0)=  One
  taui(1)=  t/t0
  do i=2,6
    taui(i)=  taui(i-1)*taui(1)
  enddo
  b=    bi(1)*log(taui(1)) + bi(0)  + bi(3)/taui(3)  + bi(5)/taui(5)
  bb=   bbi(0) + bbi(1)/taui(1) + bbi(2)/taui(2) + bbi(4)/taui(4)
  !
  if (t<=647.25D0) call Eos_H2O_psat(t,ps) !ps=  psat2(t)
  !set initial guess for rho using thb-fit to redlich-kwong
  !
  !ark=  redlich-kwong constant a 
  !brk=  redlich-kwong constant b
  !
  !the redlich-kwong equation of state is:
  !
  !  p=  R*T/(v-brk) - ark/(v*(v+brk)*sqrt(t))
  !
  !where T is absolute temp, and v is molar volume (cm3/mole)
  !
  ark=  1.279186e8 - 2.241415e4 * t
  brk=  1.428062e1 + 6.092237e-4 * t
  !
  oft=  ark/(p*sqrt(t))
  buk=  -10.0*rr*t/p
  cuk=  oft - brk*brk + brk*buk
  duk=  - brk*oft
  
  call kubik(buk,cuk,duk,x1,x2,x2i,x3)
  !
  !CALCULATE MOLAR VOLUME
  !
  if (x2i/=Zero) then
    vol=  x1
  else
    if (p<ps) then ;  vol=  MAX(x1, MAX(x2, x3))
    else           ;  vol=  MIN(x1, MIN(x2, x3))
    end if
  end if
  !
  !CALCULATE SPECifIC VOLUME
  !
  if (vol<=Zero) then  ;  rhn=  1.9
  else                 ;  rhn=  (One/vol)*18.0152
  end if
  !
  !CALCULATE parameterS THAT GO INTO BASE function:
  !
  dp_=  9.99e+10
  dr_=  9.99e+10
  !
  icount=1
  eps1=  10.0*epsilon(real(8))
  !
  !do while ((icount<=100).and.
  ! *  ((dp_>eps1).or.
  ! *   (dr_>eps1)))
  !
  do while (icount<=100)
    
    if (dp_<1.0D-06) exit
    if (dr_<1.0D-06) exit
    
    rh=  rhn
    if (rh<=Zero) rh=  1.D-8
    if (rh> 1.9) rh=  1.9D0
    
    !CALCULATE Y, WHICH IS useD IN BASE function  (HAAR ET AL., P. 272)
    y=  rh*b/4.0
    
    !CALCULATE TERMS useD IN RESIDUAL function (EQ. A.5 OF HAAR ET AL.)
    ermi(0)=  One
    ermi(1)=  One -exp(-rh)
    do i=2,9; ermi(i)=  ermi(i-1)*ermi(1); enddo
    
    pr=   Zero
    dpr=  Zero
    
    !CALCULATE RESIDUAL function
    do i=1,36
      pr= pr &
      & + gi(i)/taui(li(i))*ermi(ki(i)-1)
      dpr= dpr &
      &  + (2.0+rh*(ki(i)*exp(-rh)-One)/ermi(1)) &
      &     *gi(i)/taui(li(i))*ermi(ki(i)-1)
    enddo
    !
    do i=37,40
      del=  rh/rhoi(i) - One
      tau=  t/ttti(i)  - One
      abc=  -alpi(i)*del**ki(i)-beti(i)*tau*tau
      if (abc>-100.0D0) then  ;  q10=  gi(i)*del**li(i) * exp(abc)
      else                    ;  q10=  Zero
      end if
      !
      qm=   li(i)/del - ki(i)*alpi(i)*del**(ki(i)-1)
      pr=   pr + q10*qm*rh*rh/rhoi(i)
      dpr=  dpr &
      &  + (q10*qm*rh*rh/rhoi(i)) * (2.0/rh+qm/rhoi(i)) &
      &  - rh*rh/(rhoi(i)*rhoi(i))*q10 &
      &         * ( li(i)/del/del + ki(i)*(ki(i)-1)*alpi(i)*del**(ki(i)-2) )
    enddo
    !
    pr= rh*(rh*exp(-rh)*pr  &
    & + r*t*((One + alpha*y + beta*y*y)/f_Cube(One-y) &
    & + 4.0*y*(bb/b - gamma))) 
    dpr= rh*exp(-rh)*dpr &
    &  + r*t*( (One + 2.0*alpha*y + 3.0*beta*y*y)/f_Cube(One-y) &
    &  + 3.0*y*(One + alpha*y + beta*y*y)/f_Quart(One-y) &
    &  + 2.0*4.0*y*(bb/b - gamma))
    
    if (dpr<=Zero) then
      if (p<=ps) then  ;  rhn=  rhn*0.95
      else             ;  rhn=  rhn*1.05
      end if
      else
      if (dpr<0.01) dpr=  0.01
      x_=  (p - pr)/dpr
      if (ABS(x_)>0.1) x_=  0.1*x_/ABS(x_)
      rhn=  rh + x_
    end if
    
    dp_=  ABS(One - pr/p)
    dr_=  ABS(One - rhn/rh)
    icount=icount+1
    
  enddo
  !
  !write (6,1) t, p, ps, rhn
  !!c1  format (5x,'t=',f8.2,' p=',f8.2,' ps=',f8.2,' rhn=',f8.4)
  !
  !CALCULATE DERIVATIVES W/R TO TEMPERATURE
  rh=  rhn
  dbdt= bi(1)/t &
  &   - 3.D0*bi(3)/(taui(3)*t)  &
  &   - 5.D0*bi(5)/(taui(5)*t)
  !
  d2bdt2=          -bi(1)/(t*t) &
  &     + 4.D0*3.D0*bi(3)/(taui(3)*t*t) &
  &     + 6.D0*5.D0*bi(5)/(taui(5)*t*t)
  !
  d3bdt3= 2.D0          *bi(1)/(t*t*t) &
  &     - 5.D0*4.D0*3.D0*bi(3)/(taui(3)*t*t*t) &
  &     - 7.D0*6.D0*5.D0*bi(5)/(taui(5)*t*t*t)
  !
  dbbdt= -      bbi(1)/(taui(1)*t) &
  &      - 2.D0*bbi(2)/(taui(2)*t) &
  &      - 4.D0*bbi(4)/(taui(4)*t)
  !
  d2bbdt2= 2.D0     *bbi(1)/(taui(1)*t*t) &
  &      + 3.D0*2.D0*bbi(2)/(taui(2)*t*t) &
  &      + 5.D0*4.D0*bbi(4)/(taui(4)*t*t)
  !
  d3bbdt3= - 3.D0   *2.D0*bbi(1)/(taui(1)*t*t*t)   &
  &      - 4.D0*3.D0*2.D0*bbi(2)/(taui(2)*t*t*t)  &
  &      - 6.D0*5.D0*4.D0*bbi(4)/(taui(4)*t*t*t) 

  y=          rh *b     /4.D0
  dydt=       rh *dbdt  /4.D0
  dydrh=      b         /4.0
  d2ydt2=     rh *d2bdt2/4.D0
  d2ydtdrh=   dbdt      /4.0
  d2ydrh2=    Zero
  d3ydt3=     rh *d3bdt3/4.D0
  d3ydt2drh=  d2bdt2    /4.D0
  d3ydtdrh2=  Zero
  d3ydrh3=    Zero
  !
  ermi(0)=  One
  ermi(1)=  One-exp(-rh)
  do i=2,9
    ermi(i)=  ermi(i-1)*ermi(1)
  enddo
  !
  yy=  One-y
  yy2= yy*yy
  yy3= yy2*yy
  yy4= yy2*yy2
  !
  !calculate base function Ab and its derivatives
  Z=      - log(yy) &
  &       - (beta-One)        /yy &
  &       + (alpha+beta+One)  /(yy2*2.0D0) &
  &       + 4.0D0*(bb/b-gamma)*y &
  &       - (alpha-beta+3.0D0)/2.0D0 &
  &       + log(rh*r*t/P0)
  dZdt=     One/t &
  &       + 4.0D0*(b*dbbdt-bb*dbdt)*y/(b*b) &
  &       + 4.0D0*(bb/b-gamma)*dydt &
  &       + dydt/(yy) &
  &       + (alpha+beta+One)*dydt/yy3 &
  &       - (beta-One)*dydt/yy2
  dZdrh=    One/rh &
  &       + 4.0D0*(bb/b-gamma)*dydrh &
  &       + dydrh/(yy) &
  &       + (alpha+beta+One)*dydrh/yy3 &
  &       - (beta-One)*dydrh/yy2
  d2Zdt2= -One/(t*t) &
  &       + 4.0D0*(((b*d2bbdt2-bb*d2bdt2)*(b*b) &  
  &       - 2.0D0*(b*dbbdt-bb*dbdt)*b*dbdt)*y)/yy4 &
  &       + f_Square(dydt/(yy)) &
  &       + 3.0D0*(alpha+beta+One)*f_Square(dydt)/yy4 &
  &       - 2.0D0*(beta-One)*f_Square(dydt)/yy3 &
  &       + 8.0D0*(b*dbbdt-bb*dbdt) *dydt/(b*b) &
  &       + 4.0D0*(bb/b-gamma)      *d2ydt2 &
  &       +                          d2ydt2 /yy &
  &       + (alpha+beta+One)        *d2ydt2 /yy3 &
  &       - (beta-One)              *d2ydt2 /yy2
  d2Zdtdrh=  4.0D0*(bb/b-gamma)*d2ydtdrh + d2ydtdrh/(yy) &
  &       + (alpha+beta+One)*d2ydtdrh/yy3 &
  &       - (beta-One)*d2ydtdrh/yy2 &
  &       + 4.0D0*(b*dbbdt-bb*dbdt)*dydrh/(b*b) &
  &       + dydt*dydrh/yy2 &
  &       + 3.0D0*(alpha+beta+One)*dydt*dydrh/yy4 &
  &       - 2.0D0*(beta-One)*dydt*dydrh/yy3
  d2Zdrh2= -One/(rh*rh) &
  &      + f_Square(dydrh/(yy)) &
  &      + 3.0D0*(alpha+beta+One)*f_Square(dydrh)/yy4  &
  &      - 2.0D0*(beta-One)*f_Square(dydrh)/yy3 &
  &      + 4.0D0*(bb/b-gamma)*d2ydrh2 + d2ydrh2/(yy) &
  &      + (alpha+beta+One)*d2ydrh2/yy3 &
  &      - (beta-One)*d2ydrh2/yy2
  d3Zdt3=  2.0D0/f_Cube(t) + 4.0D0*(((-2.0D0*(b*dbbdt-bb*dbdt)*f_Square(dbdt) &
  &     + (b*d3bbdt3+d2bbdt2*dbdt-dbbdt*d2bdt2-bb*d3bdt3)*f_Square(b)   &
  &     - 2.0D0*(b*dbbdt-bb*dbdt)*b*d2bdt2)*y                           &
  &     + ((b*d2bbdt2-bb*d2bdt2)*f_Square(b)                            &
  &     - 2.0D0*(b*dbbdt-bb*dbdt)*b*dbdt)*dydt)*yy4                     &
  &     - 4.0D0*((b*d2bbdt2-bb*d2bdt2)*f_Square(b)                      &
  &     - 2.0D0*(b*dbbdt-bb*dbdt)*b*dbdt)*f_Cube(b)*y*dbdt)/(yy4*yy4)   &
  &     + 2.0D0*f_Cube(dydt/(yy)) &
  &     + 12.0D0*(alpha+beta+One)*f_Cube(dydt)/f_Quint(yy) &
  &     - 6.0D0*(beta-One)*f_Cube(dydt)/yy4 &
  &     + 8.0D0*((b*d2bbdt2-bb*d2bdt2)*f_Square(b)              &
  &     - 2.0D0*(b*dbbdt-bb*dbdt)*b*dbdt)*dydt/yy4     &
  &     + 12.0D0*(b*dbbdt-bb*dbdt)*d2ydt2/f_Square(b)                   &
  &     + 3.0D0*dydt*d2ydt2/yy2                           &
  &     + 9.0D0*(alpha+beta+One)*dydt*d2ydt2/yy4       &
  &     - 6.0D0*(beta-One)*dydt*d2ydt2/yy3                &
  &     + 4.0D0*(bb/b-gamma)*d3ydt3 + d3ydt3/(yy)                &
  &     + (alpha+beta+One)*d3ydt3/yy3                     &
  &     - (beta-One)*d3ydt3/yy2;
  d3Zdt2drh= 4.0D0*(((b*d2bbdt2-bb*d2bdt2)*f_Square(b) &
  &        - 2.0D0*(b*dbbdt-bb*dbdt)*b*dbdt)*dydrh)/yy4 &
  &        + 4.0D0*(bb/b-gamma)*d3ydt2drh + d3ydt2drh/(yy) &
  &        + (alpha+beta+One)*d3ydt2drh/yy3 &
  &        - (beta-One)*d3ydt2drh/yy2 &
  &        + 8.0D0*(b*dbbdt-bb*dbdt)*d2ydtdrh/f_Square(b)               &
  &        + 2.0D0*dydt*d2ydtdrh/yy2                      &
  &        + 6.0D0*(alpha+beta+One)*dydt*d2ydtdrh/yy4  &
  &        - 4.0D0*(beta-One)*dydt*d2ydtdrh/yy3           &
  &        + 2.0D0*f_Square(dydt)*dydrh/yy3                   &
  &        + 12.0D0*(alpha+beta+One)*f_Square(dydt)*dydrh/f_Quint(yy)  &
  &        - 6.0D0*(beta-One)*f_Square(dydt)*dydrh/yy4         &
  &        + d2ydt2*dydrh/yy2                                   &
  &        + 3.0D0*(alpha+beta+One)*d2ydt2*dydrh/yy4         &
  &        - 2.0D0*(beta-One)*d2ydt2*dydrh/yy3;
  d3Zdtdrh2=  2.0D0*dydt*f_Square(dydrh)/yy3                        &
  &        + 12.0D0*(alpha+beta+One)*dydt*f_Square(dydrh)/f_Quint(yy)  &
  &        - 6.0D0*(beta-One)*dydt*f_Square(dydrh)/yy4         &
  &        + 4.0D0*(bb/b-gamma)*d3ydtdrh2 + d3ydtdrh2/(yy)             &
  &        + (alpha+beta+One)*d3ydtdrh2/yy3 &
  &        - (beta-One)*d3ydtdrh2/yy2 &
  &        + 2.0D0*d2ydtdrh*dydrh/yy2 &
  &        + 6.0D0*(alpha+beta+One)*d2ydtdrh*dydrh/yy4 &
  &        - 4.0D0*(beta-One)*d2ydtdrh*dydrh/yy3 &
  &        + 4.0D0*(b*dbbdt-bb*dbdt)*d2ydrh2/f_Square(b) &
  &        + dydt*d2ydrh2/yy2                                   &
  &        + 3.0D0*(alpha+beta+One)*dydt*d2ydrh2/yy4         &
  &        - 2.0D0*(beta-One)*dydt*d2ydrh2/yy3
  d3Zdrh3=    2.0D0/f_Cube(rh) + 2.0D0*f_Cube(dydrh/(yy))                  &
  &        + 12.0D0*(alpha+beta+One)*f_Cube(dydrh)/f_Quint(yy)         &
  &        - 6.0D0*(beta-One)*f_Cube(dydrh)/yy4                &
  &        + 3.0D0*dydrh*d2ydrh2/yy2                            &
  &        + 9.0D0*(alpha+beta+One)*dydrh*d2ydrh2/yy4        &
  &        - 6.0D0*(beta-One)*dydrh*d2ydrh2/yy3                 &
  &        + 4.0D0*(bb/b-gamma)*d3ydrh3 + d3ydrh3/(yy)                 &
  &        + (alpha+beta+One)*d3ydrh3/yy3                       &
  &        - (beta-One)*d3ydrh3/yy2
  !
  Ab=          r*t*Z
  dAbdt=       r*Z + r*t*dZdt
  dAbdrh=      r*t*dZdrh
  d2Abdt2=     2.0D0*r*dZdt + r*t*d2Zdt2
  d2Abdtdrh=   r*dZdrh + r*t*d2Zdtdrh
  d2Abdrh2=    r*t*d2Zdrh2
  d3Abdt3=     3.0D0*r*d2Zdt2 + r*t*d3Zdt3
  d3Abdt2drh=  2.0D0*r*d2Zdtdrh + r*t*d3Zdt2drh
  d3Abdtdrh2=  r*d2Zdrh2 + r*t*d3Zdtdrh2
  d3Abdrh3=    r*t*d3Zdrh3
  !
  !calculate residual function Ar and its derivatives
  Ar=          Zero
  dArdt=       Zero
  dArdrh=      Zero
  d2Ardt2=     Zero
  d2Ardtdrh=   Zero
  d2Ardrh2=    Zero
  d3Ardt3=     Zero
  d3Ardt2drh=  Zero
  d3Ardtdrh2=  Zero
  d3Ardrh3=    Zero
  !
  do i=1,36
    Ar=          Ar + gi(i)/ki(i)/taui(li(i))*ermi(ki(i))
    dArdt=       dArdt + (-li(i)*gi(i)/ki(i)/(taui(li(i))*t) * ermi(ki(i)))
    dArdrh=      dArdrh + gi(i)/taui(li(i))*ermi(ki(i)-1)*exp(-rh)
    d2Ardt2=     d2Ardt2 + (li(i)+One)*li(i)*gi(i)/ki(i) /(taui(li(i))*t*t)*ermi(ki(i))
    d2Ardtdrh=   d2Ardtdrh + (-li(i)*gi(i)/ (taui(li(i))*t)*ermi(ki(i)-1)*exp(-rh))
    if (ki(i)>1) then
      d2Ardrh2= d2Ardrh2 &
      &       + gi(i)/taui(li(i)) &
      &         *((ki(i)-One)*ermi(ki(i)-2) *exp(-rh)-ermi(ki(i)-1)) *exp(-rh)
    else
      d2Ardrh2=  d2Ardrh2 + (-gi(i)/taui(li(i))*exp(-rh))
    end if
    !
    d3Ardt3=     d3Ardt3 - (li(i)+2.0D0)*(li(i)+One)*li(i)*gi(i)/ki(i) &
    &                     /(taui(li(i))*t*t*t) *ermi(ki(i))
    d3Ardt2drh=  d3Ardt2drh +(li(i)+One)*li(i)*gi(i) &
    &                       /(taui(li(i))*t*t) *ermi(ki(i)-1) *exp(-rh)
    !
    if (ki(i)>1) then
      d3Ardtdrh2=  d3Ardtdrh2 - li(i)*gi(i)/(taui(li(i))*t) &
      &         *((ki(i)-One)*ermi(ki(i)-2) &
      &         *exp(-rh)-ermi(ki(i)-1))*exp(-rh)
    else
      d3Ardtdrh2=  d3Ardtdrh2 + li(i)*gi(i)/(taui(li(i))*t)*exp(-rh)
    end if
    if (ki(i)>2) then
      d3Ardrh3=  d3Ardrh3 + gi(i)/taui(li(i)) &
      &       *(((ki(i)-2.0D0)*ermi(ki(i)-3)*exp(-rh)  &
      &       - 3.0D0*ermi(ki(i)-2)) &
      &       *(ki(i)-One)*exp(-rh)+ermi(ki(i)-1))*exp(-rh)
    else if (ki(i)>1) then 
      d3Ardrh3=  d3Ardrh3 - gi(i)/taui(li(i)) *(4.0D0*exp(-rh)-One)*exp(-rh)
    else
      d3Ardrh3=  d3Ardrh3 + gi(i)/taui(li(i))*exp(-rh)
    end if
  enddo
  !
  do3740: do i=37,40
    del=    rh/rhoi(i) - One
    tau=    t/ttti(i) - One
    Q=      -alpi(i)*del**ki(i) - beti(i)*tau*tau
    dQdt=   -beti(i)*2.0D0*tau/ttti(i)
    if (ki(i)==0) then  ;  dQdrh=  Zero
    else                ;  dQdrh=  -alpi(i)*ki(i)*del**(ki(i)-1)/rhoi(i)
    end if
    !
    d2Qdt2=     -beti(i)*2.0D0/f_Square(ttti(i))
    d2Qdtdrh=   Zero
    !
    if (ki(i)==0 .or. ki(i)==1) then; d2Qdrh2=  Zero
    else; d2Qdrh2=  -alpi(i)*ki(i)*(ki(i)-One)*del**(ki(i)-2) /f_Square(rhoi(i))
    end if
    !
    d3Qdt3=  Zero
    d3Qdt2drh=  Zero
    d3Qdtdrh2=  Zero
    if (ki(i)==0 .or. ki(i)==1 .or. ki(i)==2) then; d3Qdrh3=  Zero
    else; d3Qdrh3=  -alpi(i)*ki(i)*(ki(i)-One)*(ki(i)-2.0D0) *del**(ki(i)-3)/f_Cube(rhoi(i))
    end if
    !
    if (Q>-100.0D0) then; expQ= dexp(Q)
    else;                 expQ= Zero
    end if
    Ar=  Ar + gi(i)*del**li(i)*expQ
    dArdt=  dArdt + gi(i)*del**li(i)*expQ*dQdt
    !
    if (li(i)==0) then
      dArdrh= dArdrh &
      &     + gi(i)*expQ*dQdrh
    else
      dArdrh= dArdrh &
      &     + gi(i)*li(i)*del**(li(i)-1)*expQ/rhoi(i) &
      &     + gi(i)*del**li(i)*expQ*dQdrh
    end if
    d2Ardt2= d2Ardt2 &
    &      + gi(i)*del**li(i)*expQ*(f_Square(dQdt)+d2Qdt2)
    !
    if (li(i)==0) then
      d2Ardtdrh=   d2Ardtdrh + gi(i)*expQ*dQdt*dQdrh  + gi(i)*expQ*d2Qdtdrh
    else
      d2Ardtdrh= d2Ardtdrh &
      &        + gi(i)*li(i) *del**(li(i)-1)*expQ*dQdt/rhoi(i) &
      &        + gi(i)*del**li(i)*expQ*(dQdt*dQdrh+d2Qdtdrh)
    end if
    !
    if (li(i)==0) then
      d2Ardrh2= d2Ardrh2 &
      &       + gi(i)*expQ*f_Square(dQdrh) &
      &       + gi(i)*expQ*d2Qdrh2
    else if (li(i)==1) then
      d2Ardrh2= d2Ardrh2 + 2.0D0*gi(i)*expQ*dQdrh/rhoi(i) &
      &       + gi(i)*del*expQ*f_Square(dQdrh) &
      &       + gi(i)*del*expQ*d2Qdrh2
    else 
      d2Ardrh2= d2Ardrh2 &
      &       + gi(i)*li(i)*(li(i)-One)*del**(li(i)-2) *expQ/f_Square(rhoi(i)) &
      &       + 2.0D0*gi(i)*li(i)*(del**(li(i)-1))*expQ*dQdrh/rhoi(i) &
      &       + gi(i)*del**li(i)*(expQ*f_Square(dQdrh)+expQ*d2Qdrh2)
    end if
    !
    d3Ardt3= d3Ardt3 &
    &      + gi(i)*del**li(i)*expQ *(f_Cube(dQdt)+3.0D0*dQdt*d2Qdt2+d3Qdt3)
    !
    if (li(i)==0) then
      d3Ardt2drh=  d3Ardt2drh &
      &         + gi(i)*(expQ*f_Square(dQdt)*dQdrh  &
      &                + expQ*(d2Qdt2*dQdrh + dQdt*d2Qdtdrh)) &
      &         + gi(i)*(expQ*dQdt*d2Qdtdrh + expQ*d3Qdt2drh)
    else
      d3Ardt2drh=  d3Ardt2drh &
      &         + gi(i)*li(i)*del**(li(i)-1) *(expQ*f_Square(dQdt)+expQ*d2Qdt2)/rhoi(i) &
      &         + gi(i)*del**li(i)*(expQ*dQdt*(dQdt*dQdrh+d2Qdtdrh) &
      &                           + expQ*(d2Qdt2*dQdrh+dQdt*d2Qdtdrh+d3Qdt2drh))
    end if
    !
    if (li(i)==1) then
      d3Ardtdrh2= d3Ardtdrh2 &
      &         + 2.0D0*gi(i)*li(i)*(expQ*dQdt*dQdrh + expQ*d2Qdtdrh)/rhoi(i) &
      &         + gi(i)*del*(expQ*dQdt*f_Square(dQdrh)  &
      &                    + expQ*2.0D0*dQdrh*d2Qdtdrh  &
      &                    + expQ*dQdt*d2Qdrh2        &
      &                    + expQ*d3Qdtdrh2)
    else
      d3Ardtdrh2= d3Ardtdrh2 &
      &         + gi(i)*li(i)*(li(i)-One) *del**(li(i)-2) &
      &                *expQ *dQdt/f_Square(rhoi(i)) &
      &         + 2.0D0*gi(i)*li(i)*del**(li(i)-1) &
      &                *(expQ*dQdt*dQdrh + expQ*d2Qdtdrh)/rhoi(i) &
      &         + gi(i)*del**li(i)*(expQ*dQdt*f_Square(dQdrh) &
      &                           + expQ*2.0D0*dQdrh*d2Qdtdrh &
      &                           + expQ*dQdt*d2Qdrh2       &
      &                           + expQ*d3Qdtdrh2)
    end if
    !
    if (li(i)==0) then
      d3Ardrh3= d3Ardrh3 &
      &       + gi(i)*(expQ*f_Cube(dQdrh) &
      &              + expQ*2.0D0*dQdrh*d2Qdrh2 &            
      &              + expQ*dQdrh*d2Qdrh2 &
      &              + expQ*d3Qdrh3)
    else if (li(i)==1) then
      d3Ardrh3=  d3Ardrh3 &
      &       + 3.0D0*gi(i)*(expQ*f_Square(dQdrh)+ expQ*d2Qdrh2)/rhoi(i) &
      &       + gi(i)*del*(expQ*f_Cube(dQdrh)      &
      &                  + expQ*2.0D0*dQdrh*d2Qdrh2 & 
      &                  + expQ*dQdrh*d2Qdrh2     &
      &                  + expQ*d3Qdrh3)
    else if (li(i)==2) then
      d3Ardrh3= d3Ardrh3                                 &
      &       + 3.0D0*gi(i)*2.0D0*expQ*dQdrh/f_Square(rhoi(i)) &
      &       + 3.0D0*gi(i)*2.0D0*del*(expQ*f_Square(dQdrh) + expQ*d2Qdrh2)/rhoi(i)  &
      &       + gi(i)*f_Square(del)*(expQ*f_Cube(dQdrh)          &
      &                          + expQ*2.0D0*dQdrh*d2Qdrh2  &
      &                          + expQ*dQdrh*d2Qdrh2        &
      &                          + expQ*d3Qdrh3)
    else 
      d3Ardrh3=  d3Ardrh3 &
      &       + gi(i)*li(i)*(li(i)-One)*(li(i)-2.0D0) *del**(li(i)-3)*expQ/f_Cube(rhoi(i))    &
      &       + 3.0D0*gi(i)*li(i)*(li(i)-One)*del**(li(i)-2) *expQ*dQdrh/f_Square(rhoi(i))    &
      &       + 3.0D0*gi(i)*li(i)*del**(li(i)-1) *(expQ*f_Square(dQdrh) + expQ*d2Qdrh2)/rhoi(i) &
      &       + gi(i)*del**li(i)*(expQ*f_Cube(dQdrh) &
      &                         + expQ*2.0D0*dQdrh*d2Qdrh2 &
      &                         + expQ*dQdrh*d2Qdrh2 &
      &                         + expQ*d3Qdrh3)
    end if
    !
  enddo do3740
  !
  !calculate ideal gas function Ai and derivatives
  tr=      t/1.0D2
  dtrdt=   One/100.0
  !
  Z=      One + (ci(1)/tr + ci(2))*log(tr)
  dZdt=   (-ci(1)*dtrdt/f_Square(tr))*log(tr) + (ci(1)/tr + ci(2))*dtrdt/tr
  d2Zdt2= (2.0*ci(1)*f_Square(dtrdt)/f_Cube(tr))*log(tr)  & 
  &     + (-ci(1)*dtrdt/f_Square(tr))*dtrdt/tr           &
  &     + (-ci(1)*dtrdt/f_Square(tr))*dtrdt/tr           &
  &     - (ci(1)/tr + ci(2))*f_Square(dtrdt/tr)
  d3Zdt3= (-3.0*2.0*ci(1)*f_Cube(dtrdt)/f_Quart(tr))*log(tr)  &
  &     + (2.0*ci(1)*f_Square(dtrdt)/f_Cube(tr))*dtrdt/tr   &
  &     + (2.0*ci(1)*f_Square(dtrdt)/f_Cube(tr))*dtrdt/tr   &
  &     - (-ci(1)*dtrdt/f_Square(tr))*f_Square(dtrdt/tr)     &
  &     + (2.0*ci(1)*f_Square(dtrdt)/f_Cube(tr))*dtrdt/tr   &
  &     - (-ci(1)*dtrdt/f_Square(tr))*f_Square(dtrdt/tr)     &
  &     - (-ci(1)*dtrdt/f_Square(tr))*f_Square(dtrdt/tr)     &
  &     + 2.0*(ci(1)/tr + ci(2))*f_Cube(dtrdt/tr)
  !
  do i=3,18
    Z=       Z      +                         ci(i) *tr**(i-6)
    dZdt=    dZdt   + (i-6.0)                *ci(i) *tr**(i-7) *dtrdt
    d2Zdt2=  d2Zdt2 + (i-6.0)*(i-7.0)        *ci(i) *tr**(i-8) *f_Square(dtrdt)
    d3Zdt3=  d3Zdt3 + (i-6.0)*(i-7.0)*(i-8.0)*ci(i) *tr**(i-9) *f_Cube(dtrdt)
  enddo
  !
  Ai=         -r*t*Z
  dAidt=      - r*Z - r*t*dZdt
  d2Aidt2=    - 2.0*r*dZdt - r*t*d2Zdt2
  d3Aidt3=    - 3.0*r*d2Zdt2 - r*t*d3Zdt3
  !
  !calculate the total helmholtz A function and its deriviatives
  A=          Ab         + Ar         + Ai
  dAdt=       dAbdt      + dArdt      + dAidt
  dAdrh=      dAbdrh     + dArdrh
  d2Adt2=     d2Abdt2    + d2Ardt2    + d2Aidt2
  d2Adtdrh=   d2Abdtdrh  + d2Ardtdrh
  d2Adrh2=    d2Abdrh2   + d2Ardrh2
  d3Adt3=     d3Abdt3    + d3Ardt3    + d3Aidt3
  d3Adt2drh=  d3Abdt2drh + d3Ardt2drh
  d3Adtdrh2=  d3Abdtdrh2 + d3Ardtdrh2
  d3Adrh3=    d3Abdrh3   + d3Ardrh3
  !
  !calculate g=  a + p/rh and v=  1/rh
  !p=            rh*rh*dAdrh
  dpdrh=        2.0*rh*dAdrh + rh*rh*d2Adrh2
  dpdt=         rh*rh*d2Adtdrh
  drhdt=        -dpdt/dpdrh
  d2pdt2=       rh*rh*d3Adt2drh
  d2pdtdrh=     2.0*rh*d2Adtdrh + rh*rh*d3Adtdrh2
  d2pdrh2=      2.0*dAdrh + 4.0*rh*d2Adrh2 + rh*rh*d3Adrh3
  !
  !calculate other thermodynamic properties
  !("temp" is a temporary variable, not temperature ...)
  gH2O=        A + p/rh
  hH2O=        A + p/rh - t*dAdt
  sH2O=        - dAdt
  cpH2O=       -t*d2Adt2 + (t/(rh*rh))*f_Square(dpdt)/(dpdrh)
  cvH2O=       -t *d2Adt2
  vH2Ocm3=     One/rh
  dvdtH2O=     -(One/f_Square(rh))*drhdt
  dvdpH2O=     -(One/f_Square(rh))/dpdrh
  temp=        (-d2pdt2 - 2.D0*d2pdtdrh*drhdt  - d2pdrh2*f_Square(drhdt))/dpdrh
  d2vdt2H2O=   2.D0*f_Square(drhdt)/f_Cube(rh) - temp/f_Square(rh)
  temp=        (-d2pdtdrh/dpdrh - d2pdrh2*drhdt/dpdrh)/dpdrh
  d2vdtdpH2O=  -2.D0*dpdt/(f_Cube(rh)*f_Square(dpdrh)) - temp/f_Square(rh)
  d2vdp2H2O=   (One/f_Square(rh))*d2pdrh2/f_Cube(dpdrh) + (2.0/f_Cube(rh))/f_Square(dpdrh)
  temp=        - d2Adt2 &
  &            - t*d3Adt3 &
  &            + f_Square(dpdt/rh)/dpdrh &
  &            + t*(2.D0*dpdrh*dpdt*d2pdt2 - f_Square(dpdt)*d2pdtdrh) /f_Square(rh*dpdrh);
  dcpdtH2O=    temp + t*(d2vdt2H2O)*dpdt
  !
  !CONVERT FROM SPECifIC (PER G) valueS TO MOLAR valueS
  gH2O=        gH2O       *18.0152D0 /10.0D0
  hH2O=        hH2O       *18.0152D0 /10.0D0
  sH2O=        sH2O       *18.0152D0 /10.0D0
  cpH2O=       cpH2O      *18.0152D0 /10.0D0
  cvH2O=       cvH2O      *18.0152D0 /10.0D0
  dcpdtH2O=    dcpdtH2O   *18.0152D0 /10.0D0
  vH2Ocm3=     vH2Ocm3    *18.0152D0
  dvdtH2O=     dvdtH2O    *18.0152D0
  dvdpH2O=     dvdpH2O    *18.0152D0
  d2vdt2H2O=   d2vdt2H2O  *18.0152D0
  d2vdtdpH2O=  d2vdtdpH2O *18.0152D0
  d2vdp2H2O=   d2vdp2H2O  *18.0152D0
  !
  !SUBSTRACT OUT BERMAN'S value (CHECK THIS!)
  !_______________________Berman 1988____________________Haar 1977
  gH2O=        gH2O     - gref - 285829.96D0 - (298.15D0*69.9146D0)
  hH2O=        hH2O     - href - 285829.96D0                        
  vH2O_m3=     vH2Ocm3 /1.0D6
  !
  ! deltaH_0_f= -285829.96D0
  ! S_0_ref=    +69.9146D0
  !
  !CONVERT TO SI UNITS
  cp=       cpH2O /0.0180152D0
  cv=       cvH2O /0.0180152D0
  rhof=     1.0D6 *0.0180152D0 /vH2Ocm3 !-> kg/m3
  gamma=    cpH2O /cvH2O
  blkgas=   rhof *dpdrh *100.D0
  !! error with compaq !! sndspeed= sqrt((cp/cv)*blkgas/rhof)
  sndspeed= zero
  !
  !CONVERT DVDTH2O FROM CM3/MOLE TO M3/(KG K)
  dvdtH2O=  dvdth2o /0.0180152D0 /1.0d+06
  !
  return
end subroutine CalcGH2O_Haar_Ghiorso_Detail

end module M_Eos_H2O_Haar_Ghiorso
