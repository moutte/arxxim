module M_Dtb_Const
!--
!-- constant parameters used by M_T_DtbAquHkf, M_T_DtbMinThr, ...
!--
  use M_Kinds
  implicit none
  !
  private
  !
  logical,public:: DtbConv_Benson= .true.  !.false.  !
  !
  real(dp),parameter,public:: Avogadro=  6.0221415D+23  !Avogadro's number in / mol, , COdata_2002
  real(dp),parameter,public:: Boltzmann= 1.3806505D-23  !Boltzman constant in  J / K , COdata_2002
  real(dp),parameter,public:: Electron=  1.60217653D-19 !elementary charge, C, COdata_2002
  real(dp),parameter,public:: EpsVacuum= 8.8538D-12     !permittivity in vacuum, C.C/J/m
  !
  real(dp),parameter,public:: &
  & R_cK= 1.987216D0,&  !R constant, in cal/K/mole
  & R_jK= 8.314510D0    !R constant, in J/K/Mole, RobieHemingway (from COdata)
  !                     !should be =Avogadro*Boltzmann
  !
  real(dp),parameter,public:: &
  & F_JV=       96485.309D0, & !Faraday,   in J/V/Mole, RobieHemingway (from COdata) 
  & CalToJoule= 4.183999D0,  & !1 Cal = 4.183999 Joule
  & JouleToCal= 0.2390058D0    !
  !
  ! 1 bar=   10^5 Pascal
  ! 1 Joule= (1 Pascal).(1 m^3) = (10^-5 bar).(10^6 cm^3) = 10 bar.cm^3
  ! 1 bar=   0.1 Joule/cm^3 = 0.02390058 Cal/cm^3
  ! 1 cm3=   0.1 J/bar= 0.02390058 Cal/bar
  !
  real(dp),parameter,public:: &
  & OneAtm= 1.01325D0, & !1 atm in bar
  & OneMPa= 10.0D0       !1 megaPascal = 10 Bar
  !
  real(dp),parameter,public:: &
  & Pref=   1.01325D0, & !1 atm in Bars, or is it 1 bar ????
  & Tref=   298.150D0, & !in degK,(=25°C)
  & T_CK=   273.150D0    !from T_C to T_K
  !
  real(dp),parameter,public:: &
  & PminHSV=  1.0D-3,   &
  & PmaxHSV=  12.0D3,   &
  & TCminHSV= 0.00D0,   &
  & TCmaxHSV= 1.60D3
  !
  real(dp),parameter,public:: &
  !-- Robie/COdata data for entropies of elements at 298.15 --
  !-- H2: 130.68 J/K/mole
  !-- O2: 205.15 J/K/mole
  !-- C:    5.74 J/K/mole
  !-- -> for H2O, to transform G(BermBrown) to G(BensHelg),
  !--    add 298.15*(130.68 + 205.15/2)-  69544.98 J
  & S0_Hydrogen= 65.34D0, &
  & S0_Oxygen=   102.575D0
 
  !
end module M_Dtb_Const


 
