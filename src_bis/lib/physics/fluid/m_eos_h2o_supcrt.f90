module M_Eos_H2O_Supcrt
!! was M_Supcrt_H2O
!--
!-- interface to routine h2o92
!-- from program supcrt92 by J.W.Johnson, LLNL, 1991
!--

public:: CalcGH2O_Supcrt

contains

subroutine CalcGH2O_Supcrt(TdgK,Pbar,Grt,H,S,V)
  use M_Dtb_Const,only: R_jK, DtbConv_Benson,S0_Hydrogen,S0_Oxygen
  implicit none
  !
  real(8),intent(in) :: TdgK,Pbar
  real(8),intent(out):: Grt,H,S,V
  !
  integer, parameter:: NPROP2=  46
  integer:: specs(10)
  real(8):: states(4),X
  real(8):: props(NPROP2)
  logical:: error
  !
  data specs  / 2,2,2,5,1,0,0,0,0,0 /
  data states / 4*0.0d0 /
  !
  specs(6)=   0
  specs(7)=   2
  specs(8)=   0
  specs(9)=   4
  !
  !calculate H2O properties at standard state of TdgK, P bar
  states(1)=  TdgK - 273.15d0 !298.15D0 -273.15d0
  states(2)=  Pbar !1.01325D0
  !
  call Data_Consts
  call Data_HGKcon
  call Data_LVScon
  !
  call H2O92(specs,states,props,error)
  !
  X=   props(3)*4.184D0
  if(.not. DtbConv_Benson) X= X -298.15D0 *(S0_Hydrogen *2.0D0 + S0_Oxygen)
  Grt= X /TdgK /R_jK !!8.314510D0
  H=   props(9)*4.184D0
  S=   props(5)*4.184D0
  V=   18.0152d0 /states(3) /1.0D6 !mwH2O/states(3)
  !
  return
end subroutine CalcGH2O_Supcrt

subroutine H2O92(specs,states,props,error)
!H2O92
!Computes state, thermodynamic, transport, and electroststic properties of fluid H2O at T,[P,D] 
!using equations and data given by Haar et al. (1984), Levelt Sengers et al. (1983), Johnson and Norton (1991), 
!Watson et al. (1980), Sengers and Kamgar-Parsi (1984), Sengers et al. (1984), Helgeson and Kirkham (1974), 
!Uematsu and Franck (1980), and Pitzer (1983). 
!Author: James W. Johnson Earth Sciences Dept., L-219 Lawrence Livermore National Laboratory Livermore, CA 94550
!johnson@s05.es.llnl.gov
!Abandoned:  8 November 1991
!specs  - 
!       !Input unit, triple point, saturation,
!and option specs:
!       !it, id, ip, ih, itripl, isat, iopt, useLVS, epseqn, icrit;
!       !note that the returned value of isat may differ from its input value 
!       !and that icrit need not be specified prior to invocation.
!!states - State variables:
!       !temp, pres, dens(1), dens(2);
!       !note that the first three of these must be specified prior to invocation
!       !and that, in the case of saturation, vapor density is returned in dens(1), liquid in dens(2).
!props - Thermodynamic, transport, electrostatic, and combined property values:
!       ! A, G, S, U, H, Cv, Cp, Speed, 
!       ! alpha, beta, diel, visc, tcond, surten, tdiff, Prndtl, visck, albe,
!       ! ZBorn, YBorn, QBorn, daldT, XBorn
!error =
!       !logical argument that indicates success (.false.) or failure (.true.) of the call,
!       !the latter value in response to out-of-bounds specs or states variables.
!
  implicit real(8) (a-h,o-z)
  integer, parameter:: NPROP=   23
  integer, parameter:: NPROP2=  46
  !
  integer:: specs(10)
  real(8):: states(4)
  real(8):: props(NPROP2)
  logical:: error
  !
  real(8):: Dens(2)
  real(8):: wpliq(NPROP), wprops(NPROP)
  logical:: useLVS !crtreg, valid, 
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  common /wpvals/ wprops, wpliq
  save
  !
  call unit(specs(1),specs(2),specs(3),specs(4),specs(5))
  !         (it,id,ip,ih,itripl)
  !
  if (.not.(valid(specs(1),specs(2),specs(3),specs(4), &
  &               specs(5),specs(6),specs(7),specs(8),specs(9), &
  &               states(1),states(2),states(3)))) then
    error=  .true.
    !!print *,"not valid"  ; pause
    return
  else
    error=  .false.
  end if
  !
  if (crtreg(specs(6),specs(7),specs(1), states(1),states(2),states(3))) then
    specs(10)=  1
    useLVS= (specs(8) == 1)
  else
    specs(10)=  0
    useLVS=  .false.
  end if
  !
  if (useLVS) then
    Dens(1)=  states(3)
    call LVSeqn(specs(6),specs(7),specs(5), states(1),states(2),Dens,specs(9))
    Dens(1)=  Dens(1) / 1.0d3
    if (specs(6)==1) Dens(2)=Dens(2)/1.0d3
  else
    Dens(1)=  states(3) / 1.0d3
    call HGKeqn(specs(6),specs(7),specs(5), states(1),states(2),Dens,specs(9))
  end if
  !
  call load(1,wprops,props)
  !
  if (specs(6) == 1) then
    tempy=  Dens(1)
    Dens(1)=  Dens(2)
    Dens(2)=  tempy
    call load(2,wpliq,props)
  end if
  !
  states(1)=  TdegUS(specs(1),states(1))
  states(2)=  states(2) * fp
  states(3)=  Dens(1) / fd
  !
  if (specs(6)==1) states(4)= Dens(2) /fd
  !
  return
end subroutine H2O92
     
! module M_H2O_tolers
!   TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
! end module M_H2O_tolers
! 
! module M_H2O_units
!   ft, fd, fvd, fvk, fs, fp, fh, fst, fc
! end module M_H2O_units

logical function valid(it,id,ip,ih,itripl,isat,iopt,useLVS,epseqn,Temp,Pres,Dens)
!-- .true.
!-- if unit and equation specifications are valid
!-- and input state conditions 
!-- fall within the HGK equation's region of validity;
!-- .false. otherwise.

  implicit real(8) (a-h,o-z)
  
  integer  useLVS, epseqn
  
  !!logical  valspc, valTD, valTP
  common /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  common /crits/  Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon
  common /tpoint/ Utr, Str, Htr, Atr, Gtr,Ttr, Ptripl, Dltrip, Dvtrip
  common /HGKbnd/ Ttop, Tbtm, Ptop, Pbtm, Dtop, Dbtm
  common /liqice/ sDli1, sPli1, sDli37, sPli37, sDIB30, Tli13, Pli13, Dli13, TnIB30, DnIB30
  save
  
  ! ensure validity of input specifications
  if (.not. valspc(it,id,ip,ih,itripl,isat,iopt, useLVS,epseqn)) then
    valid=  .false.
    !!print *,"valspc false" ; pause
    return
  end if
  
  ! convert to  degC, bars, g/cm3 ***
  T=  TdegK(it,Temp) - 273.15d0
  D=  Dens *fd
  P=  Pres /fp *1.0d1
  Ttripl=  Ttr - 273.15d0
  Tcrit=   Tc  - 273.15d0
  Pcrit=   Pc  * 1.0d1
  !
  if (isat == 0) then
    if(iopt==1) then ; valid=  valTD(T,D,isat,epseqn)
    else             ; valid=  valTP(T,P)
    end if
  else
    if(iopt==1) then  ; valid=  ((T+FPTOL >= Ttripl) .and. (T-FPTOL <= Tcrit))
    else              ; valid=  ((P+FPTOL >= Ptripl) .and. (P-FPTOL <= Pcrit))
    end if
  end if
  
  return
end function valid

logical function valspc(it,id,ip,ih,itripl,isat,iopt, useLVS,epseqn)
!-- returns .true. if  it, id, ip, ih, itripl, isat, iopt, useLVS, 
!-- and epseqn values all define valid input;
!-- returns .false. otherwise.
  
  integer  useLVS, epseqn
  save
  
  valspc=  (1 <= it)     .and. (it     <= 4) .and. &
  &        (1 <= id)     .and. (id     <= 4) .and. &
  &        (1 <= ip)     .and. (ip     <= 5) .and. &
  &        (1 <= ih)     .and. (ih     <= 6) .and. &
  &        (0 <= itripl) .and. (itripl <= 1) .and. &
  &        (0 <= isat)   .and. (isat   <= 1) .and. &
  &        (1 <= iopt)   .and. (iopt   <= 2) .and. &
  &        (0 <= useLVS) .and. (useLVS <= 1) .and. & 
  &        (1 <= epseqn) .and. (epseqn <= 5)
  
  return
end function valspc

logical function valTD(T,D,isat,epseqn)
!Returns .true. if  T-D  defines liquid or vapor H2O within validity limits of the HGK equation of state;
!returns .false. otherwise.
  implicit real(8) (a-h,o-z)
  
  integer epseqn
  
  common /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
  common /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
  common /RTcurr/ rt
  common /crits/  Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon
  common /tpoint/ Utr, Str, Htr, Atr, Gtr, Ttr, Ptripl, Dltrip, Dvtrip
  common /HGKbnd/ Ttop, Tbtm, Ptop, Pbtm, Dtop, Dbtm
  common /liqice/ sDli1, sPli1, sDli37, sPli37, sDIB30, Tli13, Pli13, Dli13, TnIB30, DnIB30
  common /coefs/  a(20), q(20), x(11)
  common /satur/  Dliq, Dvap, DH2O, iphase
  save
  
  equivalence (TmnLVS, x(1))
  if(   (T-FPTOL > Ttop) &
    .or.(T+FPTOL < Tbtm) &
    .or.(D-FPTOL > Dtop) &
    .or.(D+FPTOL < Dbtm)) &
  then
    valTD=  .false.
    return
  end if
  
  Tcrit=   Tc  - 273.15d0
  Ttripl=  Ttr - 273.15d0
  
  if ( (T+FPTOL >= Tcrit) &
  .or. ((T >= TnIB30) .and. (D >= Dltrip))) then
    
    Dlimit=  sDIB30 * (T-TnIB30) + Dtop
    valTD=   (D-FPTOL <= Dlimit)
    
  else
    
    if (D-FPTOL <= Dltrip) then
      if (T >= Ttripl) then
        valTD=  .true.
        Tk=  T + 273.15d0
        if (Tk < TmnLVS) then
          rt=  gascon * Tk
          call pcorr(0,Tk,Ps,Dl,Dv,epseqn)
        else
          istemp=  1
          DH2O=  0.0d0
          P=  Pfind(istemp,Tk,DH2O)
          call denLVS(istemp,Tk,P)
          Dv=  Dvap / 1.0d3
          Dl=  Dliq / 1.0d3
        end if
        if ((D >= Dv) .and. (D<= Dl)) isat=  1
      else
        P=  Psublm(T)
        PMPa=  P / 1.0d1
        Tk=  T + 273.15d0
        Dguess=  PMPa / Tk / 0.4d0
        rt=  gascon * Tk
        call bb(Tk)
        call denHGK(Dsublm,PMPa,Dguess,Tk,dPdD)
        valTD=  (D-FPTOL <= Dsublm)
      end if
    else
      if (D <= Dli13) then
        Dlimit=  sDli1 * (T-Tli13) + Dli13
        valTD=  (D+FPTOL >= Dlimit)
      else
        Dlimit=  sDli37 * (T-Tli13) + Dli13
        valTD=  (D-FPTOL <= Dlimit)
      end if
    end if
    
  end if
  
  return
end function valTD

logical function valTP(T,P)
!.true. if  T-P  defines liquid or vapor H2O within validity limits of the HGK equation of state;
!.false. otherwise.

  implicit real(8) (a-h,o-z)
  
  common /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
  common /crits/  Tcrit, rhoC, Pc, Pcon, Ucon, Scon, dPcon
  common /tpoint/ Utr, Str, Htr, Atr, Gtr, Ttr, Ptripl, Dltrip, Dvtrip
  common /HGKbnd/ Ttop, Tbtm, Ptop, Pbtm, Dtop, Dbtm
  common /liqice/ sDli1, sPli1, sDli37, sPli37, sDIB30, Tli13, Pli13, Dli13, TnIB30, DnIB30
  save
  
  if ((T-FPTOL > Ttop) .or. (T+FPTOL < Tbtm) .or. (P-FPTOL > Ptop) .or. (P+FPTOL < Pbtm)) then
    valTP=  .false.
    return
  else
    valTP=  .true.
  end if

  if (P >= Pli13) then
    Plimit=  sPli37 * (T-Tli13) + Pli13
    valTP=  (P-FPTOL <= Plimit)
  else
    if (P >= Ptripl) then
      Plimit=  sPli1 * (T-Tli13) + Pli13
      valTP=  (P+FPTOL >= Plimit)
    else
      Psubl=  Psublm(T)
      valTP=  (P-FPTOL <= Psubl)
    end if
  end if
  
  return
end function valTP

real(8) function Psublm(Temp)
!-- Psublimation(T)  computed from the equation given by 
!-- Washburn (1924): Monthly Weather Rev., v.52, pp.488-490.

  implicit real(8) (a-h,o-z)
  save
  
  T=  Temp + 2.731d2
  PmmHg= power(1.0d1, &
  &  (-2.4455646d3/T + 8.2312d0*DLOG10(T) - 1.677006d-2*T + 1.20514d-5*T*T - 6.757169d0))
  
  !convert mmHg to bars ***
  Psublm=  PmmHg * 1.33322d-3
  
  return
end function Psublm

!consts - Constants
subroutine Data_Consts
!-- block data consts
  
  implicit real(8) (a-h,o-z)
  
  parameter (MAXISO=  21, MAXINC=  75, NPLOTF=  8)
  
  integer:: rterm, wterm, reacf, pronf, tabf
  integer:: plotf(NPLOTF),tempf, mapiso(2,3), mapinc(2,3), mapv3(2,3)
  integer:: rec1m1, rec1m2, rec1m3, rec1m4, rec1aa, rec1gg 
  real(8):: mwH2O, satmin(2)
  real(8)::&
  & dsvar(MAXINC,MAXISO),  Vw (MAXINC,MAXISO), &
  & bew  (MAXINC,MAXISO),  alw(MAXINC,MAXISO), dalw(MAXINC,MAXISO), &
  & Sw   (MAXINC,MAXISO),  Cpw(MAXINC,MAXISO), Hw  (MAXINC,MAXISO), &
  & Gw   (MAXINC,MAXISO),  &
  & Zw   (MAXINC,MAXISO),  Qw (MAXINC,MAXISO), Yw  (MAXINC,MAXISO), &
  & Xw   (MAXINC,MAXISO)
  
  logical lvdome(MAXINC,MAXISO), H2Oerr(MAXINC,MAXISO),EQ3run, lv1bar
  !
  character*4  incvar(2,3)
  character*10 isov(2,3), incv(2,3), var3(2,3), isosat(2)
  character*12 isovar(2,3)
  character*20 namecf, namerf, nametf, namepf(NPLOTF)
  !
  !integer rterm, wterm, reacf, pronf, tabf, plotf(6)
  common /io/     rterm, wterm, iconf, reacf, pronf, tabf, plotf
  common /io2/    tempf
  common /stvars/ isosat, isovar, incvar
  common /TPDmap/ mapiso, mapinc, mapv3
  common /headmp/ isov, incv, var3
  common /refval/ mwH2O, R, Pref, Tref, ZPrTr, YPrTr
  common /rlimit/ nmin1,  nmin2,  nmin3,  nmin4,  ngas, &
  &               naqs, rec1m1, rec1m2, rec1m3, rec1m4, rec1gg, rec1aa
  common /tranm2/ ntrm2
  common /aqscon/ eta, theta, psi, anion, cation, gref
  common /qtzcon/ aa, ba, ca, VPtTta, VPrTtb, Stran
  common /satend/ satmin
  common /defval/ DPMIN,  DPMAX,  DPINC, DTMIN, DTMAX, DTINC, DTSMIN, DTSMAX, DTSINC 
  common /null/   XNULLM, XNULLA
  common /badtd/  lvdome, H2Oerr
  common /fnames/ namecf, namerf, nametf, namepf
  common /EQ36/   EQ3run
  common /lv1b/   lv1bar
  common /H2Ogrd/ dsvar, Vw, bew, alw, dalw, Sw, Cpw, Hw, Gw, Zw, Qw, Yw, Xw
  common /H2Oss/  Dwss, Vwss, bewss, alwss, dalwss, Swss, Cpwss, Hwss, Gwss, &
  &               Zwss, Qwss, Ywss, Xwss
  save
  
  data EQ3run, lv1bar /2*.false. /
  !8=  NPLOTF
  data namepf / 8*'                    ' /
  !13*MAXISO*MAXINC=  20475
  data dsvar, Vw, bew, alw, dalw, Sw, Cpw, Hw, Gw, Zw, Qw, Yw, Xw &
  &    / 20475*0.0d0 /
  !2*MAXISO*MAXINC=  3150 
  data lvdome, H2Oerr / 3150*.false. /
  data Dwss, Vwss, bewss, alwss, dalwss, Swss, Cpwss, Hwss, Gwss, Zwss, Qwss, Ywss, Xwss &
  &    / 13*0.0d0 /
  data XNULLM, XNULLA         / 999999.0d0,999.0d0   /
  data DPMIN,  DPMAX,  DPINC  / 500.0d0,   5000.0d0,  500.0d0 /
  data DTMIN,  DTMAX,  DTINC  / 100.0d0,   1000.0d0,  50.0d0  /
  data DTSMIN, DTSMAX, DTSINC /   0.0d0,   350.0d0,   25.0d0  /
  data satmin                 / 0.01d0,    0.006117316772d0   /
  data aa, ba, ca             / 0.549824d3,0.65995d0,-0.4973d-4 /
  data VPtTta, VPrTtb, Stran  / 0.23348d2, 0.2372d2,  0.342d0   /
  data eta, theta, psi        / 0.166027d6,0.228d3,   0.26d4    /
  data anion, cation, gref    / 0.0d0,     0.94d0,    0.0d0     /
  data mwH2O, R               / 18.0152d0, 1.9872d0  /
  data Pref, Tref             /  0.1d1,    0.29815d3 /
  !ZPrTr, YPrTr calculated in SUBR getH2O as f(epseqn)
  data rterm, wterm, iconf, reacf, pronf, tabf, tempf / 5,     6,     40,    41,    42,   43,    44 /
  !8=  NPLOTF
  data plotf / 61, 62, 63, 64, 65, 66, 67, 68 /
  data isovar / 'CHORES(g/cc)', 'BARS(bars)  ', 3*'THERMS(degC)', 'BARS(bars)  ' /
  data incvar / 2*'TEMP', 'DENS'  , 2*'PRES',   'TEMP' /
  data isosat / 'TEMP(degC)',   'PRES(bars)' /
  data isov   / 'DH2O(g/cc)', 'PRES(bars)', 3*'TEMP(degC)', 'PRES(bars)' /
  data incv   / 2*'TEMP(degC)', 'DH2O(g/cc)', 'PRES(bars)', 'PRES(bars)', 'TEMP(degC)' /

  data var3   / 'PRES(bars)', 'DH2O(g/cc)', 'PRES(bars)', 3*'DH2O(g/cc)' /
  
  data mapiso / 3, 2, 1, 1, 1, 2 /
  data mapinc / 1, 1, 3, 2, 2, 1 /
  data mapv3  / 2, 3, 2, 3, 4, 4 /
  
!end block data consts
end subroutine Data_Consts

subroutine Data_HGKcon
!Constant parameters for the H2O equation of state given by  Haar, Gallagher, & Kell (1984): 
!bp, bq=      b(j), B(j) from Table A.1, p.272
!g1, g2, gf=  alpha, beta, gamma from eq (A.2), p.272
!g, ii, jj=   g(i), k(i), l(i) from eq (A.5), p.272.
!Note that  tz < tcHGK.
!Tolerence limits required in various real & inexact comparisons are set and stored in common /tolers/.
  implicit real(8) (a-h,o-z)
  common /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
  common /nconst/ g(40), ii(40), jj(40), nc
  common /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
  common /bconst/ bp(10), bq(10)
  common /addcon/ atz(4), adz(4), aat(4), aad(4)
  common /HGKcrt/ tcHGK, dcHGK, pcHGK
  common /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
  common /HGKbnd/ Ttop, Tbtm, Ptop, Pbtm, Dtop, Dbtm
  common /liqice/ sDli1, sPli1, sDli37, sPli37, sDIB30, Tli13, Pli13, Dli13, TnIB30, DnIB30
  common /tpoint/ Utripl, Stripl, Htripl, Atripl, Gtripl, Ttripl, Ptripl, Dltrip, Dvtrip
  save
  data    Ttripl, Ptripl, Dltrip, Dvtrip &
       /  2.731600d2, &
          0.611731d-2,&
          0.999778d0, &
          0.485467d-5 /
  data   Ttop,    Tbtm,    Ptop,    Pbtm,    Dtop,    Dbtm &
      /  2.25d3, -2.0d1,   3.0d4,   1.0d-3, &
         0.1380746d1,&
         0.8587455d-7 /
  data  sDli1, sPli1, sDli37, sPli37, sDIB30,&
        Tli13, Pli13, Dli13, TnIB30, DnIB30 &
      / -0.5847974d-2,&
        -0.1381808d3,&
         0.1832440d-2,&
         0.1745368d3,&
        -0.1683754d-3,&
        -0.1500000d2,&
         0.2074100d4,&
         0.1087556d1,&
         0.1450000d3,&
         0.1026316d1 /
  data    TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL &
       !!JM_2009!! / 1.0d-6, 1.0d-6, 1.0d-9, 1.0d-5, -673.5d0, 1.0d-7 /
       / 1.0d-6, 1.0d-6, 1.0d-9, 1.0d-5, -100.0d0, 1.0d-7 /
  data tcHGK, dcHGK, pcHGK / .647126d3, .322d3, .22055d2 /
  data atz /.64d3,  .64d3,  .6416d3, .27d3/
  data adz /.319d0, .319d0, .319d0,  .155d1/
  data aat /.2d5,   .2d5,   .4d5,    .25d2/
  data aad /.34d2,  .4d2,   .3d2,    .105d4/
  data wm, gascon, tz, aa, uref, sref &
      /  .1801520000d2,  .46152200d0, .647073d3, &
         .1d1,          -.43284550d4, .761808d1  /
  data g1, g2, gf /.11d2, .44333333d2, .35d1/
  data bp / .7478629d0,  -.3540782d0,  2*.0d0,  .7159876d-2, 0d0, -.3528426d-2, 3*.0d0/
  data bq / .11278334d1,  .0d0, -.5944001d0, -.5010996d1, .0d0,.63684256d0,  4*.0d0/
  data nc / 36 /
  data g /-.53062968529023d3,  .22744901424408d4, .78779333020687d3, &
  &      -.69830527374994d2,  .17863832875422d5, -.39514731563338d5, &
  &       .33803884280753d5, -.13855050202703d5, -.25637436613260d6, &
  &       .48212575981415d6, -.34183016969660d6,  .12223156417448d6, &
  &       .11797433655832d7, -.21734810110373d7,  .10829952168620d7, &
  &      -.25441998064049d6, -.31377774947767d7,  .52911910757704d7, &
  &      -.13802577177877d7, -.25109914369001d6,  .46561826115608d7, &
  &      -.72752773275387d7,  .41774246148294d6,  .14016358244614d7, &
  &      -.31555231392127d7,  .47929666384584d7,  .40912664781209d6, &
  &      -.13626369388386d7,  .69625220862664d6, -.10834900096447d7, &
  &      -.22722827401688d6,  .38365486000660d6,  .68833257944332d4, &
  &       .21757245522644d5, -.26627944829770d4, -.70730418082074d5, &
  &      -.22500000000000d0, -.16800000000000d1,                     &
  &       .5500000000000d-1, -.93000000000000d2/
  data ii / 4*0, 4*1, 4*2, 4*3, 4*4, 4*5, 4*6, 4*8, 2*2, 0, 4, 3*2, 4/
  data jj / 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, &
            2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, &
            1, 4, 4, 4, 0, 2, 0, 0/
end subroutine Data_HGKcon

subroutine Data_LVScon
!Constant parameters for the H2O critical region equation of state given by 
!Levelt Sengers, Kamgar-Parsi, Balfour,& Sengers (1983).
  implicit real(8)  (a-h,o-z)
  common /crits/ Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon
  common /coefs/ a(20), q(20), x(11)
  save
  data Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon / &
     647.067d0, 322.778d0, 22.046d0, &
     0.034070660379837018423130834983d0, 22046.0d0,&
     0.034070660379837018423130834983d3,&
     0.000000327018783663660700780197d0 /
  data a / &
    -0.017762d0,  5.238000d0,  0.000000d0, -2.549150d1, &
     6.844500d0,  0.325000d0,  1.440300d0,  0.000000d0, &
     1.375700d0,  2.366660d1,  4.820000d0,  0.294200d0, &
    -1.123260d1, -2.265470d1, -1.788760d1, -4.933200d0, &
     1.109430391161019373812391218008d0,&
    -1.981395981400671095301629432211d0,&
     0.246912528778663959151808173743d0,&
    -0.843411332867484343974055795059d0 / 
  data q / &
    -0.006000d0, -0.003000d0,  0.000000d0,  6.470670d2,&
     3.227780d2,  2.204600d1,  0.267000d0, -1.600000d0,&
     0.491775937675717720291497417773d0,    0.108500d0,&
     0.586534703230779473334597524774d0,&
    -1.026243389120214352553706598564d0,&
     0.612903225806451612903225804745d0,    0.500000d0,&
    -0.391500d0,  0.825000d0,  0.741500d0,&
     0.103245882826119154987166286332d0,&
     0.160322434159191991394857495360d0,&
    -0.169859514687100893997445721324d0 /
  data x / &
    6.430000d2,  6.453000d2,  6.950000d2,&
    1.997750d2,  4.200400d2,&
    2.09945691135940719075293945960d1,&
    2.15814057875264119875397458907d1,&
    3.0135d1, 4.0484d1,&
    175777517046267847932127026995d0,&
    380293646126229135059562456934d0 /
    
!     equivalence (cc,     a(1) ),  (pointA, q(1) ),  (Tmin1,  x(1)),
!    1            (p3,     a(2) ),  (pointB, q(2) ),  (Tmin2,  x(2)),
!    2            (delroc, a(3) ),  (delpc,  q(3) ),  (Tmax,   x(3)),
!    3            (p2,     a(4) ),  (Tc,     q(4) ),  (Dmin,   x(4)),
!    4            (p1,     a(5) ),  (rhoc,   q(5) ),  (Dmax,   x(5)),
!    5            (beta,   a(6) ),  (Pc,     q(6) ),  (Pmin1,  x(6)),
!    6            (xko,    a(7) ),  (dPcdTc, q(7) ),  (Pmin2,  x(7)),
!    7            (delTc,  a(8) ),  (slopdi, q(8) ),  (Pmax1,  x(8)),
!    8            (besq,   a(9) ),  (p11,    q(9) ),  (Pmax2,  x(9)),
!    9            (aa,     a(10)),  (alpha,  q(10)),  (sl1,    x(10)),
!    0            (delta,  a(11)),  (p00,    q(11)),  (sl2,    x(11)),
!    1            (k1,     a(12)),  (p20,    q(12)),
!    2            (muc,    a(13)),  (p40,    q(13)),
!    3            (mu1,    a(14)),  (deli,   q(14)),
!    4            (mu2,    a(15)),  (alh1,   q(15)),
!    5            (mu3,    a(16)),  (beti,   q(16)),
!    6            (s00,    a(17)),  (gami,   q(17)),
!    7            (s20,    a(18)),  (p01,    q(18)),
!    8            (s01,    a(19)),  (p21,    q(19)),
!    9            (s21,    a(20)),  (p41,    q(20))

end subroutine Data_LVScon

subroutine unit(it,id,ip,ih,itripl)
! Sets internal parameters according to user-specified choice of units.
! Internal program units are degK(T), and gm/cm**3(D); 
! all other properties are computed in dimensionless form and dimensioned at output time.
! NOTE:  conversion factors for j/g ---> cal/(g,mole) (ffh (4 & 5)) 
! are consistent with those given in Table 1, Helgeson & Kirkham (1974a) for thermal calories, 
! and differ slightly with those given by Haar et al (1984) for international calories.  
  implicit real(8) (a-h,o-z)
  real(8)  fft(4), ffd(4), ffvd(4), ffvk(4), ffs(4), ffp(5), ffh(6), ffst(4), ffcd(4), ffch(6)
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  save
  !
  data fft  /1.0d0,  1.0d0, 0.555555556d0, 0.555555556d0 /
  data ffd  /1.0d-3, 1.0d0, 1.80152d-2,    1.6018d-2     /
  data ffvd /1.0d0,  1.0d1, 0.555086816d0, 0.671968969d0 /
  data ffvk /1.0d0,  1.0d4, 1.0d4,         1.076391042d1 /
  data ffs  /1.0d0,  1.0d2, 1.0d2,         3.280833d0    /
  data ffp  /1.0d0,  1.0d1, 9.869232667d0, 1.45038d2,     1.01971d1/
  data ffh  /1.0d0,  1.0d0, 1.80152d1,     2.3901d-1,     4.305816d0, 4.299226d-1/
  data ffst /1.0d0,  1.0d3, 0.555086816d2, 0.2205061d1   /
  data ffcd /1.0d0,  1.0d-2,1.0d-2,        0.3048d0      /
  data ffch /1.0d-3, 1.0d0, 1.0d0,         0.23901d0,     0.23901d0,  0.947244d-3 /
  !
  ft=   fft(it)
  fd=   ffd(id)
  fvd=  ffvd(id)
  fvk=  ffvk(id)
  fs=   ffs(id)
  fp=   ffp(ip)
  fh=   ffh(ih)
  fst=  ffst(id)
  fc=   ffcd(id) * ffch(ih)
  if (itripl == 1)  call tpset
  return
end subroutine unit

logical function crtreg(isat,iopt,it,T,P,D)
!-- Returns .true. if input state conditions fall within the critical region of H2O; 
!-- otherwise returns .false..
!-- T, P, D, input in user-specified units, are returned in degK, MPa, kg/m3.
  
  implicit real(8) (a-h,o-z)
  
  logical  llim, ulim
  common /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
  common /coefs/ a(20), q(20), x(11)
  common /units/ ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  save
  
  equivalence &
  & (Tmin1,  x(1)),(Tmin2, x(2)),(Tmax,  x(3)), &
  & (Dmin,   x(4)),(Dmax,  x(5)),&
  & (Pbase1, x(6)),(Pbase2,x(7)),(PTmins, x(10)),(PTmaxs,x(11))
  
  T=  TdegK(it,T)
  
  if (isat == 0) then
    if (iopt == 1) then
      D= D * fd * 1.0d3
      crtreg= ((T >= Tmin1) .and. (T <= Tmax) .and. (D >= Dmin) .and. (D <= Dmax))
    else
      P= P / fp
      if ((T<Tmin1) .or. (T>Tmax)) then
        crtreg=  .false.
      else
        Pmin= Pbase1 + PTmins * (T - Tmin1)
        Pmax= Pbase2 + PTmaxs * (T - Tmin2)
        llim= (P >= Pmin)
        ulim= (P <= Pmax)
        if (llim.and.ulim) then
          crtreg= .true.
        else
          if (llim .and. (T <= Tmin2)) then
            isat1=  1
            ddummy= 0.0d0
            Pstest= Pfind(isat1,T,ddummy)
            crtreg= (P <= Pstest)
          else
            crtreg=  .false.
          end if
        end if
      end if
    end if
  else
    if(iopt==1) then;         crtreg=  (T >= Tmin1)
    else            ; P=P/fp; crtreg=  (P >= Pbase1)
    end if
  end if
  
  return
end function crtreg

subroutine HGKeqn(isat,iopt,itripl,Temp,Pres,Dens,epseqn)
! Computes thermodynamic and transport properties of H2O 
! from the equation of state given by Haar, Gallagher, & Kell (1984).
  implicit real(8) (a-h,o-z)
  parameter (NPROP=  23)
  integer epseqn
  real(8)  Dens(2), wprops(NPROP), wpliq(NPROP)
  common /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
  common /wpvals/ wprops, wpliq
  common /RTcurr/ rt
  save
  !
  rt=  gascon * Temp
  !
  call HGKsat(isat,iopt,itripl,Temp,Pres,Dens,epseqn)
  if (isat == 0) then
    call bb(Temp)
    call calcv3(iopt,itripl,Temp,Pres,Dens(1),epseqn)
    call thmHGK(Dens(1),Temp)
    call dimHGK(isat,itripl,Temp,Pres,Dens(1),epseqn)
  else
    wpliq(1:NPROP)=  wprops(1:NPROP)
    call dimHGK(2,itripl,Temp,Pres,Dens(2),epseqn)
  end if
  !
  return
end subroutine HGKeqn

subroutine HGKsat(isat,iopt,itripl,Temp,Pres,Dens,epseqn)
!-- If  isat-1, computes  Psat(T) or Tsat(P) (iopt-1,2), liquid and vapor densities, 
!-- and associated thermodynamic and transport properties.
!-- If  isat-0, checks whether  T-D or T-P (iopt-1,2) falls on or within  TOL  of the liquid-vapor surface; 
!-- if so, sets isat -- 1 and computes properties.
  implicit real(8) (a-h,o-z)
  
  real(8)  Dens(2)
  integer  epseqn
  
  common /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
  common /HGKcrt/ tcHGK, dcHGK, pcHGK
  common /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  common /tpoint/ Utr, Str, Htr, Atr, Gtr, Ttripl, Ptripl, Dltrip, Dvtrip
  common /crits/  Tc, rhoC, Pc, Pcon, Ucon, Scon, dPcon
  save
  
  if (isat == 1) then
    if (iopt==1) then
      call pcorr(itripl,Temp,Pres,Dens(1),Dens(2),epseqn)
    else
      call tcorr(itripl,Temp,Pres,Dens(1),Dens(2),epseqn)
    end if
  else
    if ((Temp>Tc).or.(Temp<Ttripl).or.((iopt==2).and.(Pres>Pc))) then
      return
    else
      call pcorr(itripl,Temp,Ptemp,dltemp,dvtemp,epseqn)
      if ( ((iopt==2) .and.(DABS(Pres-Ptemp)<=PTOL))&
      .or. ((iopt==1) &
      &    .and.( (DABS(Dens(1)-dltemp)<=DTOL) &
      &       .or.(DABS(Dens(1)-dvtemp)<=DTOL)  )  ) )&
      then
        isat=  1
        Pres=  Ptemp
        Dens(1)=  dltemp
        Dens(2)=  dvtemp
      end if
    end if
  end if
  
  return
end subroutine HGKsat

subroutine calcv3(iopt,itripl,Temp,Pres,Dens,epseqn)
!Compute the dependent state variable.
  implicit real(8) (a-h,o-z)
  
  integer  epseqn
  
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  common /qqqq/   q0, q5
  common /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
  common /fcts/   ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, dpdd, cjtt, cjth
  common /RTcurr/ rt
  save
  
  if (iopt == 1) then
    call resid(Temp,Dens)
    call base(Dens,Temp)
    call ideal(Temp)
    Pres=   rt * Dens * z + q0
  else
    if (Temp < tz) then
      call pcorr(itripl,Temp,ps,dll,dvv,epseqn)
    else
      ps=    2.0d4
      dll=   0.0d0
    end if
    !
    if (Pres > ps) then
      dguess=dll
    else
      dguess=Pres/Temp/0.4d0
    end if
    !
    call denHGK(Dens,Pres,dguess,Temp,dpdd)
    call ideal(Temp)
    !
  end if
  
  return
end subroutine calcv3
 
subroutine thmHGK(d,t)
!-- Computes thermodynamic functions in dimensionless units from the HGK equation of state:
!-- Helmholtz, Gibbs, internal energy, and enthalpy functions (ad, gd, ud, hd) are per RT;
!-- entropy and heat capacities (sd, cvd, cpd) are per R.
  implicit real(8) (a-h,o-z)
  
  common /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
  common /qqqq/   qp, qdp
  common /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
  common /resf/   ar, gr, sr, ur, hr, cvr, dpdtr
  common /idf/    ai, gi, si, ui, hi, cvi, cpi
  common /fcts/   ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, dpdd, cjtt, cjth
  common /RTcurr/ rt
  save
  
  !z=     zb + qp/rt/d
  dpdd=  rt * (zb + yb * dzb) + qdp
  ad=    ab + ar + ai - uref/t + sref !Helmholtz/RT
  gd=    ad + zb + qp/rt/d            !Gibbs/RT
  ud=    ub + ur + ui - uref/t        !U/RT
  dpdt=  rt * d * dpdtb + dpdtr
  cvd=   cvb + cvr + cvi              !Cv/R
  cpd=   cvd + t*dpdt*dpdt/(d*d*dpdd*gascon) !Cp/R
  hd=    ud + zb + qp/rt/d            !Enthalpy/RT
  sd=    sb + sr + si - sref          !Entropy/R
  dvdt=  dpdt / dpdd / d / d
  cjtt=  1.0d0 / d - t * dvdt
  cjth=  -cjtt / cpd / gascon
  
  return
end subroutine thmHGK

subroutine bb(t)
!-- Computes molecular parameters b, the "excluded volume" (eq A.3),
!-- and B, the second virial coefficient (eq A.4), in cm3/g (b1,b2) 
!-- and their first and second derivatives with respect to temperature (b1t,b1tt,b2t,b2tt).
  implicit real(8) (a-h,o-z)
  
  real(8) v(10)
  
  common /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
  common /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
  common /bconst/ bp(10), bq(10)
  save
  
  v(1)=  1.0d0
  do i=2,10; v(i)=  v(i-1) * tz / t; end do
  b1=    bp(1) + bp(2) * Dlog(1.0 / v(2))
  b2=    bq(1)
  b1t=   bp(2) * v(2) / tz
  b2t=   0.0d0
  b1tt=  0.0d0
  b2tt=  0.0d0
  do i=3,10
    b1=    b1   + bp(i) * v(i-1)
    b2=    b2   + bq(i) * v(i-1)
    b1t=   b1t  - (i-2) * bp(i) * v(i-1) / t
    b2t=   b2t  - (i-2) * bq(i) * v(i-1) / t
    b1tt=  b1tt + bp(i) * (i-2)*(i-2) * v(i-1) / t / t
    b2tt=  b2tt + bq(i) * (i-2)*(i-2) * v(i-1) / t / t
  end do
  b1tt=  b1tt - b1t / t
  b2tt=  b2tt - b2t / t
  
  return
end subroutine bb

subroutine base(d,t)
!-- Computes Abase, Gbase, Sbase, Ubase, Hbase, Cvbase -- all per RT (dimensionless) --
!-- as well as Pbase & dP/dT -- both per (DRT) -- 
!-- for the base function (ab, gb, sb, ub, hb, cvb, pb, dpdtb).
!-- See Haar, Gallagher & Kell (1979), eq(1).
  implicit real(8) (a-h,o-z)
  
  common /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
  common /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
  common /aconst/ wm, gascon, tz, a, z, dz, y, uref, sref
  save
  
  y=      .25d0 * b1 * d
  x=      1.0d0 - y
  z0=     (1.0d0 + g1*y + g2*y*y) / (x*x*x)
  z=      z0 + 4.0d0*y*(b2/b1 - gf)
  dz0=    (g1 + 2.0d0*g2*y)/(x*x*x) + 3.0d0*(1.0d0 + g1*y + g2*y*y)/(x*x*x*x)
  dz=     dz0 + 4.0d0*(b2/b1 - gf)
  pb=     z
  ab=     -log(x) &
  &       - (g2 - 1.0d0)/x &
  &       + 28.16666667d0/x/x &
  &       + 4.0d0*y*(b2/b1 - gf) &
  &       + 15.166666667d0 &
  &       + log(d*t*gascon/.101325d0)
  gb=     ab + z
  ub=     -t*b1t*(z - 1.0d0 - d*b2)/b1 - d*t*b2t
  sb=     ub - ab
  hb=     z + ub
  bb2tt=  t * t * b2tt
  cvb=    2.0d0*ub &
  &     + (z0 - 1.0d0)*(((t*b1t/b1)*(t*b1t/b1)) &
  &     - t*t*b1tt/b1) &
  &     - d*(bb2tt - gf*b1tt*t*t) &
  &     - (t*b1t/b1)*(t*b1t/b1)*y*dz0
  dpdtb=  pb/t + d*(dz*b1t/4.0d0 + b2t - b2/b1*b1t)
  return
end subroutine base

subroutine resid(t,d)
!-- Computes residual contributions to pressure (q),
!-- the Helmloltz function (ar) ,!
!-- dP/dD (q5), 
!-- Gibbs function (gr), 
!-- entropy (sr),
!-- internal energy (ur),
!-- enthalpy (hr), 
!-- isochoric heat capacity (cvr),
!-- dP/dT.  
!-- The first 36 terms of the residual function represent a global least-squares fit to experimental data 
!-- outside the critical region, 
!-- terms 37-39 affect only the immediate vicinity of the critical point, 
!-- and the last term (40) contributes only in the high pressure, low temperature region.
!-- ideal - Computes thermodynamic properties for H2O in the ideal gas state using equations given by Woolley (1979).
  implicit real(8) (a-h,o-z)
  
  real(8) qr(11), qt(10), qzr(9), qzt(9)
  
  common /resf/   ar, gr, sr, ur, hr, cvr, dpdtr
  common /qqqq/   q, q5
  common /nconst/ g(40), ii(40), jj(40), n
  common /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
  common /addcon/ atz(4), adz(4), aat(4), aad(4)
  common /RTcurr/ rt
  common /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
  save
  
  equivalence (qr(3), qzr(1)), (qt(2), qzt(1))
  
  qr(1)=  0.0d0
  q5=     0.0d0
  q=      0.0d0
  ar=     0.0d0
  dadt=   0.0d0
  cvr=    0.0d0
  dpdtr=  0.0d0
  e=      Dexp(-aa * d)
  q10=    d * d * e
  q20=    1.0d0 - e
  qr(2)=  q10
  v=      tz / t
  qt(1)=  t / tz
  do i=2,10
    qr(i+1)=  qr(i) * q20
    qt(i)=    qt(i-1) * v
  end do 
  do i=1,n
    k=  ii(i) + 1
    l=  jj(i)
    zz= k
    if(k==1) then; qp= g(i) * aa * qr(2)    * qzt(l)
    else         ; qp= g(i) * aa * qzr(k-1) * qzt(l)
    end if
    q=      q + qp
    q5=     q5 + aa*(2.0/d - aa*(1.0 - e*(k-1)/q20))*qp
    ar=     ar + g(i)*qzr(k)*qzt(l)/q10/zz/rt
    dfdt=   power(q20,DBLE(k))*(1-l)*qzt(l+1)/tz/k
    d2f=    l * dfdt
    dpt=    dfdt*q10*aa*k/q20
    dadt=   dadt  + g(i)*dfdt
    dpdtr=  dpdtr + g(i)*dpt
    cvr=    cvr   + g(i)*d2f/gascon
  end do
  qp=   0.0d0
  q2a=  0.0d0
  do j=37,40
    if(g(j)==0.0d0) exit
    k=      ii(j)
    km=     jj(j)
    ddz=    adz(j-36)
    del=    d/ddz - 1.0d0
    if (DABS(del) < 1.0d-10)  del=  1.0d-10
    !
    ex1=    -aad(j-36) * power(del,DBLE(k))
    if (ex1 < -100.0D0) then ; dex=  0.0d0
    else                     ; dex=  Dexp(ex1) *power(del,DBLE(km))
    end if
    !
    att=    aat(j-36)
    tx=     atz(j-36)
    tau=    t/tx - 1.0d0
    ex2=    -att * tau * tau
    !
    if (ex2 <= EXPTOL) then ; tex=  0.0d0
    else                    ; tex=  Dexp(ex2)
    end if
    !
    q10=    dex * tex
    qm=     km/del - k*aad(j-36)*power(del,DBLE(k-1))
    fct=    qm * d*d * q10 / ddz
    q5t=    fct*(2.0d0/d + qm/ddz) &
    &     - (d/ddz)*(d/ddz)*q10 * &
    &       (km/del/del + k*(k-1)*aad(j-36) * power(del,DBLE(k-2)))
    q5=     q5 + q5t*g(j)
    qp=     qp + g(j)*fct
    dadt=   dadt  - 2.0d0*g(j)*att*tau* q10 /tx
    dpdtr=  dpdtr - 2.0d0*g(j)*att*tau* fct /tx
    q2a=    q2a + t*g(j)*(4.0d0*att*ex2 + 2.0d0*att)*q10/tx/tx
    ar=     ar  + q10*g(j)/rt
  enddo
  sr=   -dadt / gascon
  ur=   ar + sr
  cvr=  cvr + q2a/gascon
  q=    q + qp
  
  return
end subroutine resid

subroutine ideal(t)
!-- Computes thermodynamic properties for H2O in the ideal gas state &
!-- using equations given by Woolley (1979) 
  implicit real(8) (a-h,o-z)
  
  real(8) c(18)
  common /idf/ ai, gi, si, ui, hi, cvi, cpi
  save
  
  data c / .197302710180d02,   .209662681977d2, -.483429455355d0,  &
  &        .605743189245d01,   .225602388500d2, -.987532442000d1,  &
  &       -.431355385130d01,   .458155781000d0, -.477549018830d-1, &
  &        .412384606330d-2,  -.279290528520d-3, .144816952610d-4, &
  &       -.564736587480d-6,   .162004460000d-7,-.330382279600d-9, &
  &        .451916067368d-11, -.370734122710d-13, &
  &        .137546068238d-15/
  
  tt=   t / 1.0d2
  tl=   Dlog(tt)
  gi=   -(c(1)/tt + c(2)) * tl
  hi=   (c(2) + c(1)*(1.0d0 - tl)/tt)
  cpi=  c(2) - c(1)/tt
  do i=3,18
    emult=  power(tt,DBLE(i-6))
    gi=   gi - c(i) * emult
    hi=   hi + c(i) * (i-6) * emult
    cpi=  cpi + c(i) * (i-6) * (i-5) * emult
  end do 
  ai=   gi - 1.0d0
  ui=   hi - 1.0d0
  cvi=  cpi - 1.0d0
  si=   ui - ai
  
  return
end subroutine ideal

real(8) function dalHGK(d,t,alpha)
!-- Computes/Returns (d(alpha)/dt)p(d,t,alpha) 
!-- for the Haar et al. (1983) equation of state
  implicit real(8) (a-h,o-z)
  
  real(8) tempi(4), densi(4), betai(4), alphai(4), g(40), k, l, km, lm, kp, lp
  
  integer          ll(40), kk(40)
  
  common /aconst/ wm, gascon, tz, a, z, dz, y, uref, sref
  common /ellcon/ g1, g2, gf, b1, b2, b1t, b2t, b1tt, b2tt
  common /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
  common /resf/   ar, gr, sr, ur, hr, cvr, dpdtr
  common /qqqq/   q, q5
  common /nconst/ g, kk, ll, n
  common /addcon/ tempi, densi, betai, alphai
  common /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
  save
  
  !evaluate derivatives for the base function
  y=      .25d0 * b1 * d
  x=      1.0d0 - y
  dydtp=  (d/4.0d0)*(b1t - b1*alpha)
  dbdd=   gascon*t  *((b1/4.0d0/x) * (1.0d0 - (g2-1.0d0)/x + (g1+g2+1.0d0)/x/x) + b2 - b1*gf + 1.0d0/d)
  db2dd=  gascon*t *((b1*b1/16.0d0/x/x) &
          * (1.0d0 - 2.0d0*(g2-1.0d0)/x + 3.0d0*(g1+g2+1.0d0)/x/x) - 1.0d0/d/d)
  db2ddt=  gascon*t &
           * ((b1t/4.0d0/x/x) &
           * (1.0d0 - (g2-1.0d0)*(1.0d0+y)/x + (g1+g2+1.0d0)*(1.0d0+2.0d0*y)/x/x) + b2t - gf*b1t) &
           + dbdd/t 
  db2dtp=  dbdd/t + gascon * t &
           * ((b1*dydtp/4.0d0/x/x/x) * (1.0d0 - g2 + 2.0d0*(g1+g2+1.0d0)/x) &
           + ((x*b1t + b1*dydtp)/4.0d0/x/x) &
           * (1.0d0 - (g2-1.0d0)/x + (g1+g2+1.0d0)/x/x) &
           + b2t - gf*b1t + alpha/d )
  db3ddt=  db2dd/t + gascon*t &
           * ( (b1*b1*dydtp/8.0d0/x/x/x/x) * (1.0d0 - g2 + 3.0d0*(g1+g2+1.0d0)/x) &
           + (b1*(x*b1t + b1*dydtp)/8.0d0/x/x/x) &
           * (1.0d0 - 2.0d0*(g2-1.0d0)/x + 3.0d0*(g1+g2+1.0d0)/x/x) - 2.0d0*alpha/d/d )
  db3dtt=  (db2ddt - dbdd/t)/t + gascon*t* ( &
           (b1t*dydtp/2.0d0/x/x/x/x) * (1.0d0 - g2 + &
           (g1+g2+1.0d0)*(2.0d0+y)/x) + &
           ((x*b1tt + 2.0d0*b1t*dydtp)/4.0d0/x/x/x) * (1.0d0 - &
           (g2-1.0d0)*(1+y)/x + (g1+g2+1.0d0)*(1.0d0+2.0d0*y)/x/x) &
           + b2tt - gf*b1tt ) + (t*db2dtp - dbdd)/t/t
  !evaluate derivatives for the residual function
  ! drdd=    q/d/d
  ! dr2dd=   (q5 - 2.0d0/d*q)/d/d
  ! dr2ddt=  dpdtr/d/d
  e1=   Dexp(-a * d)
  e2=   1.0d0 - e1
  tzt=  tz / t
  drdd=    0.0d0
  dr2dd=   0.0d0
  dr2ddt=  0.0d0
  dr2dtp=  0.0d0
  dr3ddt=  0.0d0
  dr3dtt=  0.0d0
  !evaluate terms 1-36
  do i=1,n
    k=  DBLE(kk(i)) + 1.0d0
    l=  DBLE(ll(i)) - 1.0d0
    km=  k - 1.0d0
    lm=  l - 1.0d0
    kp=  k + 1.0d0
    lp=  l + 1.0d0
    xtzt=  power(tzt,l)
    drdd=    drdd + g(i) * xtzt*power(e2,km)*e1
    dr2dd=   dr2dd + g(i) * e1*xtzt*power(e2,km) * (km*e1/e2 - 1.0d0)
    dr2ddt=  dr2ddt &
           - g(i)*e1*l*power(e2,km)*power(tzt,lp)/tz
    dr2dtp=  dr2dtp &
           + g(i)*e1*power(e2,km)*xtzt * ( d*alpha - l/t - km*e1*d*alpha/e2 )
    dr3ddt=  dr3ddt &
           + g(i)*(km*d*alpha*e1*e1*xtzt &
             * power(e2,k-3.0d0) &
             + e1*xtzt*power(e2,km)* (km*e1/e2 - 1.0d0) * (d*alpha -  l/t - km*d*alpha*e1/e2) )
    dr3dtt=  dr3dtt &
           + g(i)*l*e1*power(e2,km)*power(tzt,lp)/tz &
             *(lp/t + d*alpha*km*e1/e2 - d*alpha)
  end do
  !evaluate terms 37-40
  do i=37,40
    k=   DBLE(kk(i))
    l=   DBLE(ll(i))
    km=  k - 1.0d0
    lm=  l - 1.0d0
    kp=  k + 1.0d0
    lp=  l + 1.0d0
    ai=  alphai(i-36)
    bi=  betai(i-36)
    di=  densi(i-36)
    ti=  tempi(i-36)
    tau=  t/ti - 1.0d0
    del=  d/di - 1.0d0
    if (DABS(del) < 1.0d-10)  del=  1.0d-10
    !
    ex1=  -ai * power(del,k)
    if (ex1 < EXPTOL) then ; dex=  0.0d0
    else                   ; dex=  Dexp(ex1)
    end if
    ex2=   -bi * tau * tau
    !
    if (ex2 <= EXPTOL) then ; tex=  0.0d0
    else                    ; tex=  Dexp(ex2)
    end if
    ex12=   dex * tex
    !
    qm=     l/del - k*ai*power(del,km)
    xdell=  power(del,l)
    xdelk=  power(del,k)
    drdd=   drdd + g(i)*xdell*ex12/di*qm
    dr2dd=  dr2dd &
          + g(i)*xdell*ex12/di/di * (qm*qm - l/di/di - ai*k*km*power(del,k-2.0d0))
    dr2ddt= dr2ddt &
          - g(i)*2.0d0*bi*tau*ex12*xdell/ti/di*qm
    dr2dtp= dr2dtp &
          + g(i)/di &
            *(d*alpha*xdell*ex12/di/del/del*(l+ai*k*km*xdelk) &
            + qm * ( ex12 &
              * ( xdell* (k*ai*d*alpha*power(del,km)/di - 2.0d0*bi*tau/ti) &
              - l*d*alpha*power(del,lm)/di) ) )
    dr3ddt= dr3ddt &
          + g(i)/di/di *( xdell*ex12* (2.0d0*qm* &
      (l*d*alpha/di/del/del + ai*k*km*d*alpha* &
      power(del,k-2.0d0)/di) - 2.0d0*l*d*alpha/di/del &
      /del/del + ai*k*km*(k-2.0d0)*power(del,k-3.0d0)* &
      d*alpha/di) + (qm*qm - l/del/del - ai*k*km* & 
      power(del,k-2.0d0)) *(ex12*xdell*( ai*k* & 
      power(del,k-1.0d0)*d*alpha/di - 2.0d0*bi*tau/ti ) -  &
      ex12*l*power(del,l-1.0d0)*d*alpha/di) )
    dr3dtt=  &
      dr3dtt - 2.0d0*g(i)*bi/ti/di * ( tau*xdell*ex12*d* &
      alpha/del/del/di * (l + ai*k*km*power(del,k)) + &
      qm*( xdell*ex12*( ai*k*d*alpha*tau*power(del,km)/di & 
      + (1.0d0 - 2.0d0*bi*tau*tau)/ti -  &
      tau*l*d*alpha/di/del ) ) )
  end do
  !compute (d(alpha)/dT)P
  dalHGK=  &
  &  ((db3dtt + dr3dtt)*(2.0d0*(dbdd + drdd) + &
  &  d*(db2dd + dr2dd)) - &
  &  (db2ddt + dr2ddt)*(2.0d0*(db2dtp + dr2dtp) + &
  &  d*(db3ddt + dr3ddt) - d*alpha*(db2dd + dr2dd))) /&
  &  (2.0d0*(dbdd + drdd) + d*(db2dd + dr2dd)) / &
  &  (2.0d0*(dbdd + drdd) + d*(db2dd + dr2dd))
     
  return
end function dalHGK

subroutine denHGK(d,p,dguess,t,dpdd)
!-- Computes density (d in g/cm3) and dP/dD (dPdd) as f(p(MPa),t(degK))
!-- from an initial density guess (dguess).
  implicit real(8) (a-h,o-z)
  common /qqqq/   q0, q5
  common /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
  common /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
  common /RTcurr/ rt
  save
  i=   0
  d=   dguess
10 i=  i + 1
  if (d <= 0.0d0)  d= 1.0d-8
  if (d >  1.9d0)  d= 1.9d0
  !
  call resid(t,d)
  call base(d,t)
  !
  pp=    rt * d * pb + q0
  dpdd=  rt * (z + y * dz) + q5
  !*** if  dpdd < 0  assume d in 2-phase region and adjust accordingly ***
  if (dpdd > 0.0d0) GO TO 20
  !
  if (dguess >= 0.2967d0)  d= d *1.02d0
  if (dguess < 0.2967d0)   d= d *0.98d0
  if (i <= 10) GO TO 10
  !
  20 dpdx=  dpdd * 1.1d0
  if (dpdx < 0.1d0)  dpdx=  0.1d0
  dp=    DABS(1.0d0 - pp/p)
  if ((dp     < 1.0d-8) .or. &
     ((dguess > 0.3d0) .and. (dp < 1.0d-7)) .or.  &
     ((dguess > 0.7d0) .and. (dp < 1.0d-6))) return
  x=     (p - pp) / dpdx
  if (DABS(x) > 0.1d0)  x=  x * 0.1d0 / DABS(x)
  d=  d + x
  if (d <= 0.0d0)  d=  1.0d-8
  if (i <= 30) GO TO 10
  !
  return
end subroutine denHGK

real(8) function PsHGK(t)
!-- Returns an approximation to Psaturation(T) that agrees to within 0.02% of that predicted 
!-- by the HGK surface for temperatures up to within roughly a degree of the critical point.
  implicit real(8) (a-h,o-z)
  real(8)  a(8)
  save
  data a /-.78889166d1,  .25514255d1, -.6716169d1,  .33239495d2,&
          -.10538479d3,  .17435319d3, -.14839348d3, .48631602d2/
  if (T <= 314.0d0) then
    pl=     6.3573118d0 - 8858.843d0/t + 607.56335d0 * power(t,-0.6d0)
    PsHGK=  0.1d0 * Dexp(pl)
  else
    v=  t / 647.25d0
    w=  DABS(1.0d0 - v)
    b=  0.0d0
    do i=1,8
      z=  i
      b=  b + a(i)*power(w,(z + 1.0d0)/2.0d0)
    enddo  
    q=  b / v
    PsHGK=  22.093d0 * Dexp(q)
  end if
  return
end function PsHGK

real(8) function TsHGK(p)
!-- Returns Tsaturation(P)
  implicit real(8) (a-h,o-z)
  save

  TsHGK=  0.0d0
  if (p > 22.05d0)  return

  k=   0
  pl=  2.302585d0 + Dlog(p)
  tg=  372.83d0 &
  + pl*(27.7589d0 + pl*(2.3819d0 + pl*(0.24834d0 + pl*0.0193855d0)))
  
  1 if (tg < 273.15d0)  tg=  273.15d0
  
  if (tg > 647.00d0)  tg=  647.00d0
  
  if (k >= 8) then
    TsHGK=  tg
  else
    k=   k + 1
    pp=  PsHGK(tg)
    dp=  TdPsdT(tg)
    if (ABS(1.0d0 - pp/p) < 1.0d-5) then
      TsHGK=  tg
    else
      tg=  tg * (1.0d0 + (p - pp)/dp)
      GO TO 1
    end if
  end if

  return
end function TsHGK

real(8) function TdPsdT(t)
   !Returns  T*(dPsat/dT)
   implicit real(8) (a-h,o-z)
   real(8) a(8)
   save
   data a /-.78889166d1,  .25514255d1, -.6716169d1,  .33239495d2, &
           -.10538479d3,  .17435319d3, -.14839348d3, .48631602d2/
   v=  t / 647.25d0
   w=  1.0 - v
   b=  0.0d0
   c=  0.0d0
   do 4 i=1,8
        z=  i
        y=  a(i) * power(w,(z + 1.0d0)/2.0d0)
        c=  c + y/w*(0.5d0 - 0.5d0*z - 1.0d0/v)
   4 b=  b + y
   q=       b / v
   TdPsdT=  22.093d0 * Dexp(q) * c
   return
end function TdPsdT

subroutine corr(itripl,t,p,dl,dv,delg,epseqn)
!Computes liquid and vapor densities (dliq & dvap) and  (Gl-Gv)/RT  (delg) 
!for T-P conditions on or near the saturation surface.
  implicit real(8) (a-h,o-z)
  integer epseqn
  common /qqqq/   q00, q11
  common /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
  common /fcts/   ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, dpdd, cjtt, cjth
  common /basef/  ab, gb, sb, ub, hb, cvb, pb, dpdtb
  common /RTcurr/ rt
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  common /HGKcrt/ tcHGK, dcHGK, pcHGK
  save
  call bb(t)
  dguess=  dl
  if (dl <= 0.0d0)  dguess=  1.11d0 - 0.0004d0*t
  call denHGK(dl,p,dguess,t,dpdd)
  call ideal(t)
  call thmHGK(dl,t)
  !save liquid properties
  call dimHGK(1,itripl,t,p,dl,epseqn)
  gl=     gd
  dguess= dv
  if (dv<=0.0d0) dguess=  p / rt
  call denHGK(dv,p,dguess,t,dpdd)
  if (dv<5.0d-7) dv=  5.0d-7
  call ideal(t)
  call thmHGK(dv,t)
  !vapor properties will be available in common /fcts/ (dimensionless) after pcorr's final call of corr (delg < 10d-4)
  gv=    gd
  delg=  gl - gv
  return
end subroutine corr

!** pcorr - Computes Psaturation(T) (p) and liquid and vapor
!           densities (dl & dv) from refinement of an initial
!           approximation (PsHGK(t)) in accord with  Gl=  Gv.

subroutine pcorr(itripl,t,p,dl,dv,epseqn)
!Computes Psaturation(T) (p) and liquid and vapor densities (dl & dv) 
!from refinement of an initial approximation (PsHGK(t)) in accord with  Gl=  Gv.
  implicit real(8) (a-h,o-z)
  integer epseqn
  common /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
  save
  p=   PsHGK(t)
  dl=  0.0d0; dv=  0.0d0
  do
    call corr(itripl,t,p,dl,dv,delg,epseqn)
    dp=  delg * gascon * T / (1.0d0/dv - 1.0d0/dl)
    p=   p + dp
    if (DABS(delg)<1.0d-4) exit
  enddo
  return
end subroutine pcorr

subroutine tcorr(itripl,t,p,dl,dv,epseqn)
!-- Computes Tsaturation(P) (t) and liquid and vapor densities (dl & dv)
!-- from refinement of an initial approximation (TsHGK(p)) in accord with  Gl-  Gv.
  implicit real(8) (a-h,o-z)
  
  integer  epseqn
  
  common /aconst/ wm, gascon, tz, aa, zb, dzb, yb, uref, sref
  common /RTcurr/ rt
  save
  
  t= TsHGK(p)
  if (t == 0.0d0) return
  dl=  0.0d0
  dv=  0.0d0
  do
    rt=  t * gascon
    call corr(itripl,t,p,dl,dv,delg,epseqn)
    dp=  delg * gascon * t / (1.0d0/dv - 1.0d0/dl)
    t=  t * (1.0d0 - dp/TdPsdT(t))
    if (DABS(delg) < 1.0d-4) exit
  end do
  
  return
end subroutine tcorr

subroutine LVSeqn(isat,iopt,itripl,T,P,Dens,epseqn)
!Computes thermodynamic and transport properties of critical region H2O (369.85-419.85 degC, 0.20-0.42 gm/cm3)
!from the fundamental equation given by Levelt Sengers, et al (1983): J.Phys.Chem.Ref.Data, V.12, No.1, pp.1-28.
  implicit real(8) (a-h,o-z)
  parameter (NPROP=  23)
  real(8)  wprops(NPROP), wpliq(NPROP), Dens(2)
  logical           cpoint
  integer           epseqn
  common /coefs/  a(20), q(20), x(11)
  common /crits/  Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
  common /therm/  AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw, heat, Speed
  common /satur/  Dliq, Dvap, DH2O, iphase
  common /param/  r1, th1
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  common /wpvals/ wprops, wpliq
  save
  cpoint=  .false.
  DH2O=  Dens(1)
  10 call LVSsat(iopt,isat,T,P,DH2O)
  if ((isat /= 0) .or. (iopt /= 1))  call denLVS(isat,T,P)
  if (isat == 0) then
    Dens(1)=  DH2O
  else
    Dens(1)=  Dliq
    Dens(2)=  Dvap
  end if
  if (isat == 0) then
    call thmLVS(isat,T,r1,th1)
    call dimLVS(isat,itripl,th1,T,P*1.0d1,dl,dv,wprops,epseqn)
    if (cpoint) then
      call cpswap
      Dens(1)=  cdens
      Dens(2)=  cdens
      isat=  1
      iopt=  ioptsv
    end if
  else
    th1=  -1.0d0
    call thmLVS(isat,T,r1,th1)
    call dimLVS(isat,itripl,th1,T,P*1.0d1,dl,dv,wprops,epseqn)
    th1=   1.0d0
    call thmLVS(isat,T,r1,th1)
    call dimLVS(isat,itripl,th1,T,P*1.0d1,dl,dv,wpliq,epseqn)
    if (dl == dv) then
      cpoint=  .true.
      cdens=  dl
      T=  647.0670000003d0
      P=   22.0460000008d0
      ioptsv=  iopt
      iopt=  2
      isat=  0
      GO TO 10 
    end if
  end if
end subroutine LVSeqn

subroutine cpswap
!Load critical point A, G, U, H, S, Vs, Di, ZB, 
!albe values from wpliq into wprops and 
!approximations to critical Cv, Cp, alpha, beta, 
!visc, tcond, Prndtl, tdiff, visck, YB, QB, XB, daldT, st values from wprops into wpliq.
  implicit real(8) (a-h,o-z)
  parameter (NPROP=  23)
  integer:: &
    aw,  gw,  sw,  uw,  hw,  cvw, cpw,  vsw,  alw, bew, &
    diw, viw, tcw, stw, tdw, Prw, vikw, albew, &
    ZBw, YBw, QBw, dalwdT, XBw
  real(8)  wprops(NPROP), wpliq(NPROP)
  common /wpvals/ wprops, wpliq
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  save
  data aw,gw,sw,uw,hw,cvw,cpw,vsw,alw,bew,diw,viw,tcw,stw,tdw,Prw,vikw,albew,ZBw,YBw,QBw,dalwdT,XBw &
     /  1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,  17,   18, 19, 20, 21,    22, 23 /
  wprops(aw)=     wpliq(aw)
  wprops(gw)=     wpliq(gw)
  wprops(sw)=     wpliq(sw)
  wprops(uw)=     wpliq(uw)
  wprops(hw)=     wpliq(hw)
  wprops(diw)=    wpliq(diw)
  wprops(ZBw)=    wpliq(ZBw)
  wprops(stw)=    wpliq(stw)
  wpliq(cvw)=     wprops(cvw)
  wpliq(cpw)=     wprops(cpw)
  wpliq(alw)=     wprops(alw)
  wpliq(bew)=     wprops(bew)
  wpliq(YBw)=     wprops(YBw)
  wpliq(QBw)=     wprops(QBw)
  wpliq(XBw)=     wprops(XBw)
  wpliq(tcw)=     wprops(tcw)
  wpliq(tdw)=     wprops(tdw)
  wpliq(Prw)=     wprops(Prw)
  wpliq(dalwdT)=  wprops(dalwdT)
  wpliq(albew)=   wprops(albew)
  wprops(vsw)=    0.429352766443498d2 * fs
  wprops(viw)=    1.0d6
  wprops(vikw)=   1.0d6
  wpliq(vsw)=     wprops(vsw)
  wpliq(viw)=     wprops(viw)
  wpliq(vikw)=    wprops(vikw)
end subroutine cpswap     

subroutine LVSsat(iopt,isat,T,P,D)
!If  isat=1,  computes  Psat(T) or Tsat(P) (iopt=1,2).
!If  isat=0,  checks whether  T-D or T-P (iopt=1,2) falls on or within TOL 
!of the liq-vap surface; if so, isat <- 1  and  T <- Tsat.
  implicit real(8) (a-h,o-z)
  common /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
  common /crits/  Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
  save
  data ERRTOL, TCTOL / 1.0d-12, 1.0d-2 /
  if (isat == 1) then
    if (iopt == 1) then
       P=  Pfind(isat,T,D)
    end if
    T=  TsLVS(isat,P)
  else
    if (iopt == 1) then
       P=  Pfind(isat,T,D)
    end if
    if (P-ERRTOL > Pc) then
      return
    else
      call backup
      Tsat=  TsLVS(isat,P)
      if (DABS(Tsat-T) < TCTOL) then
        T=  Tsat
        isat=  1
      else
        call restor
      end if
    end if
  end if
  return
end subroutine LVSsat

subroutine denLVS(isat,T,P)
!Calculates  DH2O(T,P)  or  Dvap,Dliq(T,P) from the Levelt Sengers, et al (1983) critical region equation of state.
  implicit real(8) (a-h,o-z)
  real(8) s(2), sd(2)
  common /coefs/ a(20), q(20), x(11)
  common /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
  common /satur/ Dliq, Dvap, DH2O, iphase
  common /therm/ AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw, heat, Speed
  common /param/ r1, th1
  common /deri2/ dPdD, dPdT
  save
  equivalence (Dmin, x(4)), (Dmax, x(5)), (pw11, q(9)), (xk0,  a(7)), (xk1,  a(12))
  if (isat == 0) then
    DH2O=  rhoc
    do i=1,20
      Pnext=  Pfind(isat,T,DH2O)
      Pdif=   Pnext - P
      if (iphase == 2) then
        if (DABS(Pdif) <= 0.0d0) then
          return
        else
        end if
        if (Pdif < 0.0d0) then
          DH2O=  Dmax
        else
          DH2O=  Dmin
        end if
      else
        delD=   -Pdif/dPdD
        DH2O=  DH2O + delD
        if (DH2O < Dmin)  DH2O=  Dmin
        if (DH2O > Dmax)  DH2O=  Dmax
        if (DABS(delD/DH2O) < 1.0d-6)  return
      end if
    end do
  else
    Tw=    -Tc/T
    dTw=   1.0d0 + Tw
    call ss(r1,th1,s,sd)
    rho1=  1.0d0+pw11*dTw+a(1)*(s(1)+s(2))
    rho2=  xk0*power(r1,a(6)) + xk1*power(r1,q(16))
    Dvap=  rhoc * (rho1 - rho2)
    Dliq=  rhoc * (rho1 + rho2)
    return
  end if
  return
end subroutine denLVS

!TsLVS - Returns saturation T(P) 
real(8) function TsLVS(isat,P)
  implicit real(8) (a-h,o-z)
  common /therm/ AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw, heat, Speed
  common /satur/ Dliq, Dvap, DH2O, iphase
  common /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
  common /deri2/ dPdD, dPdT
  save
  !
  TsLVS2=  Tc - 1.0d0
  D=  rhoc
  do 10 i=1,20
    Pnext=  Pfind(isat,TsLVS2,D)
    dT=  (Pnext - P)/dPdT
    TsLVS2=  TsLVS2 - dT
    if (TsLVS2 > Tc) then
      TsLVS2=  Tc
    else
      if (DABS(dT/TsLVS2) < 1.0d-8) then
        GO TO 20
      else
      end if
    end if
  10 continue
  20 TsLVS=  TsLVS2
  return
end function TsLVS

!Pfind
!Returns P(T,D)
!Computes (dP/dD)T when invoked by SUB Dfind (isat=0)
!and (dP/dT)D when invoked by SUB TsLVS (isat=1)
!Also computes 1st & 2nd partial derivatives the singular part of the potential (Delta P tilde)
!that are used in SUB thmLVS.
real(8) function Pfind(isat,T,D)
  implicit real(8) (a-h,o-z)
  real(8) s(2), xk(2), sd(2)
  common /coefs/ a(20), q(20), x(11)
  common /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
  common /satur/ Dliq, Dvap, DH2O, iphase
  common /therm/ AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw, heat, Speed
  common /param/ r1, th1
  common /tolers/ TTOL, PTOL, DTOL, XTOL, EXPTOL, FPTOL
  common /deriv/ amu, s, Pw, Tw, dTw, dM0dT, dP0dT, d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
  common /deri2/ dPdD, dPdT
  common /abc2/  r, th
  save
  equivalence (Pw1, a(5)),   (Pw2, a(4)),   (Pw3,  a(2)),&
              (amc, a(13)),  (am1, a(14)),  (am2,  a(15)),&
              (am3, a(16)),  (p00, q(11)),  (p20,  q(12)),&
              (p40, q(13)),  (p01, q(18)),  (p21,  q(19)),&
              (p41, q(20)),  (aa,  a(10)),  (xk0,  a(7)),&
              (xk1, a(12)),  (pw11,q(9)),   (alpha,q(10)),&
              (alhi,q(15)),  (besq,a(9))
  xk(1)=  xk0
  xk(2)=  xk1
  if (DABS(T-Tc) < FPTOL)  T=  Tc
  Tee=    (T-Tc)/Tc
  Tw=     -Tc/T
  dTw=    1.0d0 + Tw
  if (isat == 0) then
    rho=  D / rhoc
    call conver(rho,Tee,amu,th1,r1,rho1,s,rhodi,err)
  else
    th1=  -1.0d0
    th=   th1
    r1=   dTw/(1.0d0-besq)
    r=    r1
    call ss(r1,th1,s,sd)
    rho=  th1 * (xk0*power(r1,a(6)) &
        + xk1*power(r1,q(16))) &
        + a(1)*(s(1)+s(2))
    rho=  1.0d0+pw11*dTw+rho
    amu=  0.0d0
    D=  rho * rhoc
  end if
  tt1=  th1*th1
  tt2=  tt1*tt1
  Pw0=   1.0d0+dTw*(Pw1+dTw*(Pw2+dTw*Pw3))
  if (isat == 0) then
    Pwmu=  amu*rhodi
  else
    Pwmu=  0.0d0
  end if
  p0th=  p00+p20*tt1+p40*tt2
  p1th=  p01+p21*tt1+p41*tt2
  dPw0=  xk0*p0th*power(r1,2.0d0-alpha)
  dPw1=  xk1*p1th*power(r1,2.0d0-alhi)
  dPw=   aa*(dPw0+dPw1)
  Pw=    Pw0 + Pwmu + dPw
  Pfind=  Pw * Pcon * T
  if (DABS(th1) < 1.0d0) then
    iphase=  1
  else
    iphase=  2
    dP0dT=  Pw1+dTw*(2.0d0*Pw2+3.0d0*Pw3*dTw)
    dM0dT=  am1+dTw*(2.0d0*am2+3.0d0*am3*dTw)
    Uw=     dP0dT-rho*dM0dT+pw11*amu+s(1)+s(2)
    dPdTcd=  Uw + rho*dM0dT
    dPwdTw=  Pw - Tw*dPdTcd
    dPdT=    Pcon * dPwdTw
  end if
  call aux(r1,th1,d2PdT2,d2PdMT,d2PdM2,aa,xk,sd,Cvcoex)
  if (iphase == 1) dPdD=  dPcon * D * T / d2PdM2
  return
end function Pfind

!aux - Calculates some second derivatives of the anomalous part of the equation of state.
subroutine aux(r1,th1,d2PdT2,d2PdMT,d2PdM2,aa,xk,sd,Cvcoex)
  implicit real(8) (a-h,o-z)
  real(8) xk(2), s(2), sd(2), w(2), y(2), z(2), coex(2)
  common /coefs/ a(20), q(20), x(11)
  save
  equivalence (cc,   a(1)),  (beta, a(6)),  (besq,a(9)),&
              (delta,a(11)), (alpha,q(10)), (s00, a(17)),&
              (s20,  a(18)), (s01,  a(19)), (s21, a(20))
  !
  deli=   0.0d0
  s(1)=    s00+s20*th1*th1
  s(2)=    s01+s21*th1*th1
  sd(1)=   2.0*th1*s20
  sd(2)=   2.0*th1*s21
  ww=      0.0d0
  yy=      0.0d0
  zz=      0.0d0
  gamma=   beta*(delta-1.0d0)
  tt1=     th1*th1
  ter=     2.0d0*beta*delta-1.0d0
  g=       (1.0+(besq*ter-3.0)*tt1 - besq*(ter-2.0)*tt1*tt1)
  Cvcoex=  0.0d0
  do i=1,2
    alhi=     alpha - deli
    beti=     beta + deli
    gami=     gamma - deli
    if (r1 /= 0.0d0) then
      w(i)=   (1.0-alhi)*(1.0-3.0*tt1)*s(i) - beta*delta*(1.0-tt1)*th1*sd(i)
      w(i)=   (w(i)*power(r1,-alhi))/g
      w(i)=   w(i) * xk(i)
      ww=     ww + w(i)
      y(i)=   beti*(1.0d0-3.0d0*tt1)*th1 - beta*delta*(1.0d0-tt1)*th1
      y(i)=   (y(i)*power(r1,beti-1.0d0)) * xk(i) / g
      yy=     yy + y(i)
      z(i)=   1.0d0-besq*(1.0d0-(2.0d0*beti))*tt1
      z(i)=   (z(i)*power(r1,-gami)) * xk(i) / g
      zz=     zz + z(i)
      a1=  (beta*(delta-3.0d0)-3.0d0*deli-besq*alhi*gami) &
         / (2.0d0*besq*besq*(2.0d0-alhi)*(1.0d0-alhi)*alhi)
      a2=  1 &
         +((beta*(delta-3.0d0)-3.0d0*deli-besq*alhi*ter) &
          /(2.0d0*besq*(1.0d0-alhi)*alhi))
      a2=  -a2
      !
      a4=  1.0d0+((ter-2.0d0)/(2.0d0*alhi))
      f1=  a1 + a2 + a4
      !
      coex(i)=  ((2.0d0-alhi)*(1.0d0-alhi)*power(r1,-alhi)*f1*xk(i))
      Cvcoex=   Cvcoex + coex(i)
    end if
    deli=  0.5d0
  end do
  d2PdT2=  aa * ww
  d2PdMT=  yy + aa*cc*ww
  d2PdM2=  zz/aa + 2.0d0*cc*yy + cc*cc*aa*ww
  return
end subroutine aux

!Transforms  T,D  to  parametric variables  r,theta according to the revised and scaled equations.
subroutine conver(rho,Tee,amu,th1,r1,rho1s,s1,rhodi,error1)
  implicit real(8) (a-h,o-z)
  real(8) s1(2), sd(2)
  common /coefs/ a(20), q(20), x(11)
  common /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
  common /abc2/  r, th
  save
  equivalence &
    (beta,a(6)),  (delta,a(11)),  (xk1,  a(12)),&
    (cc,  a(1)),  (alhi, q(15)),  (alpha,q(10)),&
    (besq,a(9)),  (p11,  q(9)),   (deli, q(14)),&
    (p1w, q(18)), (p2w,  q(19)),  (p4w,  q(20)),&
    (aa,  a(10)), (xk0,  a(7)),   (s00,  a(17)),&
    (s20, a(18)), (betai,q(16))
  Tstar=   Tee + 1.0d0
  dtstin=  1.0d0 - (1.0d0 / Tstar)
  r1=      dtstin
  if (dtstin <= 0.0d0)  then
    r1=   dtstin/(1.0d0-besq)
    th1=  1.0d0
  else
    th1=  0.0d0
  end if
  call ss(r1,th1,s1,sd)
  rhodi=   1.0d0 + p11*dtstin
  rhodit=  rhodi + cc*s1(1) + cc*s1(2)
  drho=    rho - rhodit
  amu=     0.0d0
  if (dtstin <= 0.0d0) then
    rho1co=  xk0*power(r1,beta) + xk1*power(r1,betai)
    twofaz=  rho1co
    if (DABS(drho) <= twofaz) then
      rho1s=   DSIGN(rho1co,drho) + cc*s1(1)
      th1=     DSIGN(1.00d0,drho)
      error1=  1.0d0
      r=  r1
      th=  th1
      return
    end if
  end if
  if (drho == 0.0d0) then
    th1=    0.0d0
    r1=     dtstin
    rho1s=  cc*s1(1)
  end if 
  !rule for first pass ***
  y1=    dtstin      
  den1=  rho - rhodit
  call rtheta(r1,th1,den1,y1)
  !
  tt=    th1*th1
  amu=   aa*power(r1,beta*delta)*th1*(1.0d0-tt)
  y1=    dtstin + cc*amu
  !
  call ss(r1,th1,s1,sd)
  !
  rhoweg=  xk1*power(r1,betai)*th1 + cc*s1(2)      
  rho1s=   den1 + cc*s1(1) + rhoweg
  error1=  rho - rhodi - rho1s
  r=   r1
  th=  th1
  !
  if (DABS(error1) < 1.0d-5) then
       return
  end if
  !
  !rule for second pass ***
  !
  den12=  rho - rhodi - cc*s1(1) + rhoweg
  !
  if (den12 == den1) den12=  den1 - 1.0d-6
  !
  call rtheta(r1,th1,den12,y1)
  !
  tt=   th1*th1
  amu=  aa*power(r1,beta*delta)*th1*(1.0d0-tt)
  y1=   dtstin + cc*amu
  !
  call ss(r1,th1,s1,sd)
  !
  rhoweg=  xk1*power(r1,betai)*th1 + cc*s1(2)
  rho1s2=  den12 + cc*s1(1) + rhoweg
  error2=  rho - rhodi - rho1s2
  !
  if (DABS(error2) <= 1.0d-5) then
    r=   r1
    th=  th1
    error1=  error2
    rho1s=   rho1s2
    return      
  end if
  !_____________________________!rule for nth pass ***
  den2=    den12
  do isig=1,10
    slope=   (error2-error1)/(den2-den1)
    hold=    den2
    den2=    den1 - (error1/slope)
    !
    call rtheta(r1,th1,den2,y1)
    tt=   th1*th1
    amu=  aa*power(r1,beta*delta)*th1*(1.0d0-tt)
    y1=   dtstin + cc*amu
    call ss(r1,th1,s1,sd)
    rhoweg=  xk1*power(r1,betai)*th1 + cc*s1(2)
    rho1s=   den2 + cc*s1(1) + rhoweg
    error1=  error2
    error2=  rho - rhodi - rho1s
    r=   r1
    th=  th1
    if (DABS(error2) < 1.0d-6) return
    den1=  hold
  end do
  !
  return
end subroutine conver

!rtheta - Fits data for  1.0 < theta < 1.000001.
!Solves:
!  rho=  em*theta*(r**beta)
!  Tee=  r*(1.0d0-besq*theta*theta)
!Routine given by Moldover (1978): Jour. Res. NBS, v. 84, n. 4, p. 329 - 334.
subroutine rtheta(r,theta,rho,Tee)
  implicit real(8) (a-h,o-z)
  common /coefs/ a(20), q(20), x(11)
  save
  equivalence (beta,a(6)), (em,a(7)), (besq,a(9))
  if (em <= 0.0d0  .or.  besq <= 1.0d0) GO TO 600
  absrho=  DABS(rho)
  if (absrho < 1.0d-12) GO TO 600
  bee=  sqrt(besq)
  if (DABS(Tee) < 1.0d-12) GO TO 495
  if (Tee < 0.0d0) then
    z=  1.0d0-(1.0d0-bee)*Tee/(1.0d0-besq) * power(em/absrho,1.0d0/beta)
  else
    z=  power(1.0d0+Tee*power(em/bee/absrho,1.0d0/beta), -beta)
  end if
  if (z > 1.00234d0*bee) GO TO 496
  c=  -rho*bee/em/power(DABS(Tee),beta)
  z=  DSIGN(z,rho)
  do n=1,16
    z2=  z*z
    z3=  1.0d0 - z2
    dz=  z3*(z+c*power(DABS(z3),beta))/(z3+2.0d0*beta*z2)
    z=   z - dz
    if (DABS(dz/z) < 1.0d-12) GO TO 498
  end do
  601  if (DABS(theta) > 1.0001d0) theta=  theta/DABS(theta)
  return
  498  theta=  z/bee
  r=      Tee/(1.0d0-z*z)
  r=      DABS(r)
  return
  495  theta=  DSIGN(1.0d0,rho)/bee
  r=      power(rho/(em*theta),1.0d0/beta)
  return
  496  theta=  DSIGN(1.0d0,rho)
  r=      Tee/(1.0d0-besq)
  r=      DABS(r)
  return
  600  if (DABS(Tee) < 1.0d-12) GO TO 601
  if (Tee < 0.0d0) GO TO 496
  theta=  1.0d-12
  r=      Tee
  return
end subroutine rtheta

subroutine ss(r,th,s,sd)!terms of the summation that defines dPotl/dT
!                       !and the 1st derivative of the theta (s) square polynomial.
  implicit real(8) (a-h,o-z)
  real(8) s(2), sd(2), sx(2)
  common /coefs/ a(20), q(20), x(11)
  common /abc1/  dPdM
  save
  equivalence (alpha,q(10)),  (beta,a(6)),  (besq,a(9)),&
              (delta,a(11)),  (deli,q(14)), (alhi,q(15)),&
              (beti, q(16)),  (gami,q(17)), (p00, q(11)),&
              (p01,  q(18)),  (s00, a(17)), (s20, a(18)),&
              (s01,  a(19)),  (s21,  a(20))
  tt=     th*th
  sx(1)=   s00 + s20*tt
  sd(1)=  2.0d0*s20*th
  sx(2)=   s01 + s21*tt
  sd(2)=  2.0d0*s21*th
  s(1)=   sx(1)*a(10)*a(7)*power(r,1.0d0-alpha)
  s(2)=   sx(2)*a(10)*a(12)*power(r,1.0d0-alhi)
  dPdM=   power(r,beta)*a(7)*th  + a(1)*power(r,1.0d0-alpha)* a(10)*a(7)*sx(1) + &
          power(r,beti)*a(12)*th + a(1)*power(r,1.0d0-alhi)* a(10)*a(12)*sx(2)
  return
end subroutine ss

subroutine thmLVS(isat,T,r1,th1)!thermodynamic and transport properties of critical region H2O
!                               !using Levelt Sengers et al (1983) EoS
  implicit real(8) (a-h,o-z)
  real(8) s(2), xk(2), sd(2)
  common /coefs/ a(20), q(20), x(11)
  common /crits/ Tc, rhoc, Pc, Pcon, Ucon, Scon, dPcon
  common /therm/ AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw,  heat, Speed
  common /satur/ Dliq, Dvap, DH2O, iphase
  common /deriv/ amu, s, Pw, Tw, dTw, dM0dT, dP0dT, d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
  common /deri2/ dPdD, dPdT
  common /abc1/  dPdM
  common /abc3/  dPdTcd
  save
  equivalence (pw2, a(4)),   (pw3, a(2)),  (besq,  a(9)),&
              (amc, a(13)),  (am1, a(14)), (am2,   a(15)),&
              (aa,  a(10)),  (xk0, a(7)),  (am3,   a(16)),&
              (xk1, a(12)),  (pw11,q(9)),  (alpha, q(10)),&
              (alhi,q(15)),  (pw1, a(5))
  d2P0dT=  2.0d0*pw2 + 6.0d0*pw3*dTw
  d2M0dT=  2.0d0*am2 + 6.0d0*am3*dTw
  dP0dT=   pw1+dTw*(2.0d0*pw2+3.0d0*pw3*dTw)
  dM0dT=   am1+dTw*(2.0d0*am2+3.0d0*am3*dTw)
  if (isat == 0) then
    rho=     DH2O / rhoc
    Uw=      dP0dT-rho*dM0dT+pw11*amu+s(1)+s(2)
  else
    rho=     th1 * (xk0*power(r1,a(6)) + xk1*power(r1,q(16))) + a(1)*(s(1)+s(2))
    rho=     1.0d0+pw11*dTw+rho
    Uw=      dP0dT-rho*dM0dT+pw11*amu+s(1)+s(2)
    DH2O=    rho * rhoc
    dPdT2=   Pw - Tw*(Uw+rho*dM0dT)
    heat=    1.0d3*T*(Pcon*dPdT2)*(1.0d0/Dvap-1.0d0/Dliq)
    call ss(r1,th1,s,sd)
    call aux(r1,th1,d2PdT2,d2PdMT,d2PdM2,aa,xk,sd,Cvcoex)
    if (r1 /= 0.0d0) then
      dPdD=  dPcon * DH2O * T / d2PdM2
    end if
  end if
  if (r1 /= 0.0d0) then
    dPdTcd=  dP0dT+pw11*(amu-rho/d2PdM2)+s(1)+s(2) -  d2PdMT*rho/d2PdM2
    dPwdTw=  Pw - Tw*dPdTcd
    dPdTal=  Pcon * dPwdTw
    CviTw2=  d2P0dT - rho*d2M0dT + d2PdT2 - (pw11+d2PdMT)*(pw11+d2PdMT)/d2PdM2
    Cvw=     CviTw2 * Tw*Tw
    Cpw=     Cvw + d2PdM2*dPwdTw*dPwdTw / (rho*rho)
    betaw=   1.0d0 / (DH2O*dPdD)
    alphw=   betaw * dPdTal
    Speed=    1.0d3 * sqrt(Cpw/Cvw*dPdD)
  else
    Cvw=    1.0d0
    Cpw=    1.0d0
    betaw=  1.0d0
    alphw=  1.0d0
    Speed=  0.0d0
  end if
  Hw=  Pw - Tw*Uw
  Sw=  Hw - rho*(amu+amc+dTw*(am1+dTw*(am2+dTw*am3)))
  Scond=   Scon/DH2O
  U=       Uw * Ucon/DH2O
  H=       Hw * Scond * T
  entrop=  Sw * Scond
  AE=      U - T * entrop
  GE=      H - T * entrop
  Cv=      Cvw * Scond
  Cp=      Cpw * Scond
  return
end subroutine thmLVS

!dalLVS - Computes/returns (d(alpha)/dt)p(D,T,alpha) for the Levelt Sengers et al. (1983) equation of state
!Note that D (kg/m**3),T (degK), P (MPa), alpha (degK**-1).
real(8) function dalLVS(D,T,P,alpha)
  implicit real(8) (a-h,o-z)
  real(8) sss(2), xk(2), s(2), dsdT(2), sp(2), dspdT(2),&
                   k(2), calpha(2), cbeta(2), cgamma(2),&
                   u(2), v(2), w(2), dudT(2), dvdT(2), dwdT(2)
  common /coefs/ aa(20), qq(20), xx(11)
  common /crits/ Tc, Dc, Pc, Pcon, Ucon, Scon, dPcon
  common /deriv/ amu, sss, Pw, Tw, dTw, dM0dT, dP0dT,&
                 d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
  common /deri2/ dPdD, dPdT
  common /abc1/  dPdM
  common /abc2/  r,th
  common /abc3/  dPdTcd
  save
  equivalence (a,   aa(10)), (c,   aa(1)),  (delta,  aa(11)),&
              (bsq, aa(9)),  (P11, qq(9)),  (Delta1, qq(14)),&
              (P1,  aa(5)),  (P2,  aa(4)),  (P3,     aa(2)),&
              (s00, aa(17)), (s01, aa(19)), (s20,    aa(18)),&
              (s21, aa(20))
  if (r == 0.0d0) then
       dalLVS=  1.0d6
       return
  end if
  k(1)=       aa(7)
  k(2)=       aa(12)
  calpha(1)=  qq(10)
  calpha(2)=  qq(15)
  cbeta(1)=   aa(6)
  cbeta(2)=   qq(16)
  cgamma(1)=  cbeta(1)*(delta - 1.0d0)
  cgamma(2)=  cgamma(1) - Delta1
  delT=       (T - Tc) / T
  s(1)=       s00 + s20*th**2
  s(2)=       s01 + s21*th**2
  sp(1)=      2.0d0*s20*th
  sp(2)=      2.0d0*s21*th
  ! Compute drdT and d0dT from solution of the linear system
  !                      ax=  b
  ! d(dPdM)/dT=  -D/Dc*alpha - P11*Tc/T**2=  ar1*drdT + a01*d0dT=  b1
  ! d(delT)/dT=            Tc/T**2=          ar2*drdT + a02*d0dT=  b2
  b1=  -D/Dc*alpha - P11*Tc/T/T
  b2=   Tc/T**2   
  ar1=  0.0d0                 
  a01=  0.0d0
  do i=  1,2
    ar1=  ar1 + k(i) * (cbeta(i)*th*power(r,cbeta(i)-1.0d0) &
        + a*c*(1.0d0 - calpha(i))*power(r,-calpha(i))*s(i))
    a01=  a01 + k(i) * (power(r,cbeta(i)) + a*c*sp(i)       &
        * power(r,1.0d0-calpha(i)))
  end do
  ar2=  1.0d0 - bsq*th**2 - a*c*cbeta(1)*delta &
      * (1.0d0 - th**2)*th*power(r,(cbeta(1)*delta - 1.0d0))
  a02=  3.0d0*a*c*th**2*power(r,cbeta(1)*delta) &
      - 2.0d0*bsq*r*th - a*c*power(r,cbeta(1)*delta)
  !solve the linear system with simplistic GE w/ partial pivoting
  if (DABS(ar1) > DABS(ar2)) then 
    amult=  -ar2 / ar1
    d0dT=   (b2 + amult*b1) / (a02 + amult*a01)
    drdT=   (b1 - a01*d0dT) / ar1
  else
    amult=  -ar1 / ar2
    d0dT=   (b1 + amult*b2) / (a01 + amult*a02)
    drdT=   (b2 - a02*d0dT) / ar2
  end if
  !Compute theta polynomials and their tempertaure derivatives
  dsdT(1)=    2.0d0*s20*th*d0dT
  dsdT(2)=    2.0d0*s21*th*d0dT
  dspdT(1)=   2.0d0*s20*d0dT
  dspdT(2)=   2.0d0*s21*d0dT
  !
  q=      1.0d0 + (bsq*(2.0d0*cbeta(1)*delta - 1.0d0) - 3.0d0) &
        * th**2 - bsq*(2.0d0*cbeta(1)*delta - 3.0d0)*th**4
  !
  dqdT=   2.0d0*(bsq*(2.0d0*cbeta(1)*delta - 1.0d0) - 3.0d0) * th * d0dT &
        - 4.0d0*bsq*(2.0d0*cbeta(1)*delta - 3.0d0) * th**3 * d0dT
  !
  do i=  1,2
    u(i)=     (1.0d0 - bsq*(1.0d0 - 2.0d0*cbeta(i))*th**2) / q
    dudT(i)=  (-2.0d0*bsq*(1.0d0 - 2.0d0*cbeta(i))*th*d0dT - &
              u(i)*dqdT) / q
    v(i)=     ((cbeta(i) - cbeta(1)*delta)*th + &
              (cbeta(1)*delta - 3.0d0*cbeta(i))*th**3) / q
    dvdT(i)=  ((cbeta(i) - cbeta(1)*delta)*d0dT + &
              3.0d0*(cbeta(1)*delta - 3.0d0*cbeta(i))* &
              th**2*d0dT - v(i)*dqdT) / q
    w(i)=     ((1.0d0 - calpha(i))*(1.0d0 - 3.0d0*th**2)* &
              s(i) - cbeta(1)*delta*(th - th**3)*sp(i)) / q
    dwdT(i)=  ((1.0d0 - calpha(i))*((1.0d0 - 3.0d0*th**2)* &
              dsdT(i) - 6.0d0*th*s(i)*d0dT) - cbeta(1)* &
              delta*((th - th**3)*dspdT(i) + sp(i)* &
              (d0dT - 3.0d0*th**2*d0dT)) - w(i)*dqdT) / q
  end do
  !Compute dP0dTT, ddelMT, dPdTT, dPdMMT, dPdMTT, dPPTT
  dP0dTT=  Tc/T**2 * (2.0d0*P2 + 6.0d0*P3*delT) 
  ddelMT=  a*power(r,cbeta(1)*delta)* (cbeta(1)*delta*th/r* &
           (1.0d0 - th**2)*drdT + (1.0d0 - 3.0d0*th**2)*d0dT)
  dPdTT=   0.0d0
  dPdMMT=  0.0d0
  dPdMTT=  0.0d0
  do i=  1,2
    dPdTT=   dPdTT + a*k(i) * (power(r,1.0d0-calpha(i))* &
             dsdT(i) + s(i)*(1.0d0 - calpha(i))* &
             power(r,-calpha(i))*drdT)
    !
    dPdMMT=  dPdMMT + k(i) * ((power(r,-cgamma(i))*dudT(i) - & 
             u(i)*cgamma(i)*power(r,-1.0d0-cgamma(i))*drdT) / &
             a + 2.0d0*c*(power(r,cbeta(i)-1.0d0)*dvdT(i) +  &
             v(i)*(cbeta(i) - 1.0d0)*power(r,cbeta(i)-2.0d0)* &
             drdT) + a*c**2*(power(r,-calpha(i))*dwdT(i) -  &
             calpha(i)*w(i)*power(r,-1.0d0-calpha(i))*drdT))
    !
    dPdMTT=  dPdMTT + k(i) * (power(r,cbeta(i)-1.0d0)*dvdT(i) + &
             v(i)*(cbeta(i) - 1.0d0)*power(r,cbeta(i)-2.0d0)* &
             drdT + a*c*(power(r,-calpha(i))*dwdT(i) -  &
             calpha(i)*power(r,-1.0d0-calpha(i))*drdT*w(i)))
    !
  end do
  !
  dPPTT=  dP0dTT + dPdTT + P11*ddelMT - D/Dc*dPdMTT/d2PdM2 &
        + (P11 + d2PdMT)*(D/Dc*alpha/d2PdM2 & 
        + D/Dc*dPdMMT/d2PdM2**2)
  !
  pterm=  P/Pc + dPdTcd
  !compute (d(alpha)/dT)P
  !
  dalLVS=   Tc*Dc**2/D**2/T**2 * (-2.0d0/T*d2PdM2*pterm  &
          + 2.0d0*alpha*d2PdM2*pterm + pterm*dPdMMT &
          + d2PdM2*dPPTT)
  return
end function dalLVS

!backup - Save Pfind common values during saturation check.
subroutine backup
  implicit real(8) (a-h,o-z)
  real(8) s(2), xk(2)
  common /satur/ Dliq, Dvap, DH2O, iphase
  common /param/ r1, th1
  common /deriv/ amu, s, Pw, Tw, dTw, dM0dT, dP0dT, d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
  common /deri2/ dPdD, dPdT
  common /store/ sav2, sav3, sav4, sav5, sav6, sav7, sav8, & 
                 sav9, sav10, sav11, sav12, sav13, sav14, sav15, & 
                 sav16, sav17, sav18, sav19, isav1
  save
  isav1=  iphase
  sav2=   r1
  sav3=   th1
  sav4=   amu
  sav5=   s(1)
  sav6=   s(2)
  sav7=   Pw
  sav8=   Tw
  sav9=   dTw
  sav10=  dM0dT
  sav11=  dP0dT
  sav12=  d2PdM2
  sav13=  d2PdMT
  sav14=  d2PdT2
  sav15=  p0th
  sav16=  p1th
  sav17=  xk(1)
  sav18=  xk(2)
  sav19=  dPdD
  return
end subroutine backup

!restor - Restore Pfind common values after saturation check.
subroutine restor
  implicit real(8) (a-h,o-z)
  real(8) s(2), xk(2)
  common /satur/ Dliq, Dvap, DH2O, iphase
  common /param/ r1, th1
  common /deriv/ amu, s, Pw, Tw, dTw, dM0dT, dP0dT, d2PdM2, d2PdMT, d2PdT2, p0th, p1th, xk
  common /deri2/ dPdD, dPdT
  common /store/ sav2, sav3, sav4, sav5, sav6, sav7, sav8, &
                 sav9, sav10, sav11, sav12, sav13, sav14, sav15, & 
                 sav16, sav17, sav18, sav19, isav1
  save
  !
  iphase=  isav1
  r1=      sav2
  th1=     sav3
  !
  amu=     sav4
  s(1)=    sav5
  s(2)=    sav6
  Pw=      sav7
  Tw=      sav8
  dTw=     sav9
  dM0dT=   sav10
  dP0dT=   sav11
  d2PdM2=  sav12
  d2PdMT=  sav13
  d2PdT2=  sav14
  p0th=    sav15
  p1th=    sav16
  xk(1)=   sav17
  xk(2)=   sav18
  dPdD=    sav19
  !
  return
end subroutine restor

subroutine load(phase,ptemp,props)
!Load thermodynamic and transport property values from ptemp into props.
  implicit real(8) (a-h,o-z)
  parameter (NPROP=   23, NPROP2=  46)
  real(8)  ptemp(NPROP), props(NPROP2)
  integer           phase, key(NPROP,2)
  save
  data  key &
     /  1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23, &
       25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, & 
        2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, & 
       26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46  /
  do i=  1,NPROP 
    props(key(i,phase))=  ptemp(i)
  end do
  return
end subroutine load

!Dimension triple point  U, S, H, A, G  values (in J/g from
!Table 2, Helgeson & Kirkham, 1974a) into user-specified units.
subroutine tpset
  implicit real(8) (a-h,o-z)
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  common /tpoint/ Utripl, Stripl, Htripl, Atripl, Gtripl, Ttripl, Ptripl, Dltrip, Dvtrip
  save
  data       Utr,        Str,       Htr,        Atr,        Gtr &
       / -15766.0d0,  3.5144d0, -15971.0d0, -12870.0d0, -13073.0d0 /
  Utripl=  Utr * fh
  Stripl=  Str * fh
  Htripl=  Htr * fh
  Atripl=  Atr * fh
  Gtripl=  Gtr * fh
end subroutine tpset

!Convert  U, S, H, A, G  values computed with reference to zero triple point properties
!(Haar et al., 1984; Levelt Sengers et al., 1983)
!into values referenced to triple point properties given by  Helgeson and Kirkham, 1974a.

subroutine triple(T,wpzero)
  implicit real(8) (a-h,o-z)
  parameter (NPROP=  23)
  real(8)  wpzero(NPROP)
  integer  A, G, S, U, H
  common /tpoint/ Utr, Str, Htr, Atr, Gtr, Ttripl, Ptripl, Dltrip, Dvtrip
  save
  data    A, G, S, U, H  / 1, 2, 3, 4, 5 /
  wpzero(S)=  wpzero(S) + Str
  TS=  T*wpzero(S) - Ttripl*Str
  wpzero(G)=  wpzero(H) - TS + Gtr            
  wpzero(A)=  wpzero(U) - TS + Atr 
  wpzero(H)=  wpzero(H) + Htr
  wpzero(U)=  wpzero(U) + Utr
end subroutine triple

!Returns  base**exp  utilizing the intrinsic FORTRAN exponentiation function
!in such a manner so as to insure computation of machine-independent values for all defined exponentiations
!Attempted undefined exponentiations produce an error message and cause program termination.
real(8) function power(base,exp)
  implicit real(8) (a-h,o-z)
  integer rterm, wterm, reacf, pronf, tabf, plotf(6)
  common /io/ rterm, wterm, iconf, reacf, pronf, tabf, plotf
  save
  data TOL / 1.0d-7 /
  if (base > 0.0d0) then
    power=  base**exp
  else
    if (DABS(base) > TOL) then
      if (DBLE(INT(exp)) /= exp) then
        write(wterm,10) base, exp
        10 format(/,' neg base ** real exp is complex',/,' base,exp: ',2e20.13,/)
        stop
      else
        if (MOD(exp,2.0d0) == 0.0d0) then
          power=   (-base)**exp
        else
          power=  -((-base)**exp)
        end if
      end if
    else
      if (exp > 0.0d0) then
        power=  0.0d0
      else
        write(wterm,20) base, exp
        20 format(/,' zero base ** (exp <= 0) is undefined', /,' base,exp: ',2e20.13)
        stop
      end if
    end if
  end if
  return
end function power

!TdegK - Returns input temperature  t  converted from  user-specified units to degK.
real(8) function TdegK(it,t)
  implicit real(8) (a-h,o-z)
  save
  GO TO (1,2,3,4), it
  1 TdegK=  t
  return
  2 TdegK=  t + 273.15d0
  return
  3 TdegK=  t / 1.8d0
  return
  4 TdegK=  (t + 459.67d0) / 1.8d0
  return
end function TdegK

real(8) function TdegUS(it,t)
!Returns input temperature  t  converted from degK to user-specified units.
  implicit real(8) (a-h,o-z)
  save
  GO TO (1,2,3,4), it
  1 TdegUS=  t
  return
  2 TdegUS=  t - 273.15d0
  return
  3 TdegUS=  t * 1.8d0
  return
  4 TdegUS=  t * 1.8d0 - 459.67d0
  return
end function TdegUS

!dim[HGK,LVS] - Dimensioning routines for H2O88.
!dimHGK - Dimensions thermodynamic and transport property values computed from the HGK equation of state 
!per user-specified choice of units.     
subroutine dimHGK(isat,itripl,t,p,d,epseqn)
  implicit real(8) (a-h,o-z)
  parameter (NPROP=  23)
  real(8)  wprops(NPROP), wpliq(NPROP)
  integer::&
    aw, gw, sw, uw, hw, cvw, cpw, vsw, alw, bew, &
    diw, viw, tcw, stw, tdw, Prw, vikw, albew, &
    ZBw, YBw, QBw, dalwdT, XBw
  integer  epseqn    
  common /units/  ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  common /fcts/   ad, gd, sd, ud, hd, cvd, cpd, dpdt, dvdt, dpdd, cjtt, cjth
  common /aconst/ wm, gascon, tz, aa, z, dz, y, uref, sref
  common /RTcurr/ rt    
  common /wpvals/ wprops, wpliq
  save
  data aw, gw, sw, uw, hw, cvw, cpw, vsw, alw, bew, diw, viw, &
       tcw, stw, tdw, Prw, vikw, albew, ZBw, YBw, QBw, & 
       dalwdT, XBw &
     /  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12, &
       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 /
  wprops(aw)=    ad * rt * fh
  wprops(gw)=    gd * rt * fh
  wprops(sw)=    sd * gascon * fh * ft
  wprops(uw)=    ud * rt * fh
  wprops(hw)=    hd * rt * fh
  wprops(cvw)=   cvd * gascon * fh * ft
  wprops(cpw)=   cpd * gascon * fh * ft
  wprops(vsw)=   sqrt(DABS(cpd*dpdd*1.0d3/cvd)) * fs
  wprops(bew)=   1.0d0 / (d * dpdd * fp)
  wprops(alw)=   d * dvdt
  wprops(dalwdT)=  dalHGK(d,t,wprops(alw))

  pbars=  p*1.0d1
  dkgm3=  d * 1.0d3   
  betaPa=  wprops(bew)*fp / 1.0d6
  betab=   wprops(bew)*fp / 1.0d1
  CpJKkg=  wprops(cpw)/fh/ft * 1.0d3

  wprops(viw)=   viscos(t,pbars,dkgm3,betaPa) * fvd
  wprops(tcw)=   thcond(t,pbars,dkgm3,wprops(alw),betaPa) * fc * ft
  if ((isat==0).or.(isat==2)) then; wprops(stw)=  0.0d0
                              else; wprops(stw)=  surten(t) * fst; end if
  call Born92(&
    t,pbars,dkgm3/1.0d3,betab,&
    wprops(alw),wprops(dalwdT), wprops(diw),&
    wprops(ZBw),wprops(QBw),wprops(YBw),wprops(XBw),epseqn)
  wprops(tdw)=    wprops(tcw)/fc/ft  / (dkgm3 * CpJKkg) * fvk
  if (wprops(tcw) /= 0.0d0) then; wprops(Prw)=  wprops(viw)/fvd * CpJKkg / (wprops(tcw)/fc/ft)
                            else; wprops(Prw)=  0.0d0; end if
  wprops(vikw)=   wprops(viw)/fvd / dkgm3 * fvk
  wprops(albew)=  wprops(alw) / wprops(bew)
  if (itripl == 1) call triple(t,wprops)
end subroutine dimHGK

subroutine dimLVS(isat,itripl,theta,T,Pbars,dl,dv,tprops,epseqn)
!Dimension critical region properties per user-specs and load into tprops.
  implicit real(8) (a-h,o-z)
  parameter (NPROP=  23)
  real(8)  tprops(NPROP)
  integer  &
    aw, gw, sw, uw, hw, cvw, cpw, vsw, alw, bew, diw, viw,&
    tcw, stw, tdw, Prw, vikw, albew, &
    ZBw, YBw, QBw, dalwdT, XBw
  integer  epseqn
  common /therm/   AE, GE, U, H, Entrop, Cp, Cv, betaw, alphw, heat, Speed
  common /satur/   Dliq, Dvap, DH2O, iphase
  common /units/   ft, fd, fvd, fvk, fs, fp, fh, fst, fc
  common /abc2/    r, th
  save
  data &
    aw, gw, sw, uw, hw, cvw, cpw,  vsw, alw, bew, diw, viw, &
    tcw, stw, tdw, Prw, vikw,albew,ZBw, YBw, QBw, &
    dalwdT, XBw &
  / 1,  2,  3,  4,  5,  6,   7,   8,  9,  10, 11, 12, &  
    13, 14, 15, 16, 17, 18,  19,  20, 21, 22, 23 /
  if (isat == 1) then
      dv=    Dvap
      dl=    Dliq
  end if
  tprops(aw)=   AE * fh
  tprops(gw)=   GE * fh
  tprops(sw)=   Entrop * fh * ft
  tprops(uw)=   U * fh
  tprops(hw)=   H * fh
  tprops(cvw)=  Cv * fh * ft
  tprops(cpw)=  Cp * fh * ft
  tprops(vsw)=  Speed * fs
  tprops(bew)=  betaw / fp
  tprops(alw)=  alphw
  th=  theta
  tprops(dalwdT)=  dalLVS(DH2O,T,Pbars/1.0d1,tprops(alw))
  CpJKkg=   Cp * 1.0d3
  betaPa=   betaw / 1.0d6
  betab=    betaw / 1.0d1
  if (DABS(theta) /= 1.0d0) then
    dkgm3=  DH2O
    tprops(stw)=  0.0d0
  else
    if (theta < 0.0d0) then
      dkgm3=  Dvap
      tprops(stw)=  0.0d0
    else
      dkgm3=  Dliq
      dkgm3=  Dliq
      tprops(stw)=  surten(T) * fst
    end if
  end if
  call Born92(T,Pbars,dkgm3/1.0d3,betab,tprops(alw),tprops(dalwdT), &
              tprops(diw),tprops(ZBw),tprops(QBw),tprops(YBw), &
              tprops(XBw),epseqn)
  tprops(viw)=    viscos(T,Pbars,dkgm3,betaPa) * fvd
  tprops(tcw)=    thcond(T,Pbars,dkgm3,tprops(alw),betaPa) * fc * ft
  tprops(tdw)=    tprops(tcw)/fc/ft  / (dkgm3 * CpJKkg) * fvk
  tprops(Prw)=    tprops(viw)/fvd * CpJKkg / (tprops(tcw)/fc/ft)
  tprops(vikw)=   tprops(viw)/fvd / dkgm3 * fvk
  tprops(albew)=  tprops(alw) / tprops(bew)
  if (itripl == 1)  call triple(T,tprops)
end subroutine dimLVS      

!tran88 - 
!Set of FORTRAN77 functions that compute transport properties of fluid H2O.
!Input state parameters should be computed 
!from the Haar et al. (1984) and Levelt Sengers et al. (1983) equations of state 
!in order to facilitate comparision with published tabular values referenced below for each function.

!programmer:  James W. Johnson
!abandoned:   20 January 1988
!viscos -
!Returns dynamic viscosity of H2O in kg/m*s (= Pa*s)
!if Tk, Pbars falls within the validity region (specified by the initial if statement) of the Watson et al. (1980) equation
!otherwise returns zero
!See equations 3.1-2 and 4.1-5 and Tables 1, 6, and 8 from Sengers and Kamgar-Parsi (1984).
real(8) function viscos(Tk,Pbars,Dkgm3,betaPa)
  implicit real(8) (a-h,o-z)
  parameter (Tstar=  647.270d0)
  parameter (Dstar=  317.763d0)
  parameter (Pstar=  22.1150d6)
  parameter (ustar=  1.0d-6)
  real(8)  a(4), b(6,7)
  save
  data a / 0.0181583d0,  0.0177624d0,  0.0105287d0, -0.0036744d0 /
  data b / 0.5132047d0,  0.3205656d0,  0.0d0,        0.0d0, &
          -0.7782567d0,  0.1885447d0,  0.2151778d0,  0.7317883d0, &
           1.2410440d0,  1.4767830d0,  0.0d0,        0.0d0, &
          -0.2818107d0, -1.0707860d0, -1.2631840d0,  0.0d0, &
           0.0d0,        0.0d0,        0.1778064d0,  0.4605040d0, &
           0.2340379d0, -0.4924179d0,  0.0d0,        0.0d0, &
          -0.0417661d0,  0.0d0,        0.0d0,        0.1600435d0, &
           0.0d0,        0.0d0,        0.0d0,       -0.01578386d0, &
           0.0d0,        0.0d0,        0.0d0,        0.0d0, &
           0.0d0,        0.0d0,        0.0d0,      -0.003629481d0, &
           0.0d0,        0.0d0 /
  data TOL /1.0d-2/
  viscos=  0.0d0
  TdegC=   Tk - 273.15d0
  if ((Pbars > 5000.0d0+TOL) .or. &
     ((Pbars > 3500.0d0+TOL).and.(TdegC > 150.0d0+TOL)).or. & 
     ((Pbars > 3000.0d0+TOL).and.(TdegC > 600.0d0+TOL)) .or. &
     (TdegC  > 900.0d0+TOL))  return
  T=  Tk / Tstar
  D=  Dkgm3 / Dstar
  sum=  0.0d0
  do i=0,3; sum=  sum + a(i+1)/T**i; end do
  u0=  ustar * sqrt(T) / sum
  sum=  0.0d0
  do i=0,5
    do j=0,6; sum=  sum + b(i+1,j+1) * (1.0d0/T-1)**i * (D-1)**j; end do
  end do
  u1=  Dexp(D*sum)
  if ((0.997d0 <= T) .and. (T <= 1.0082d0) .and. &
      (0.755d0 <= D) .and. (D <= 1.2900d0)) then
    xt=  Pstar/Dstar**2 * betaPa * Dkgm3**2
    if (xt < 22.0d0) then; u2=  1.0d0
                     else; u2=  0.922 * power(xt,0.0263d0); end if
  else
    u2=  1.0d0
  end if
  viscos=  u0 * u1 * u2
  return
end function viscos

real(8) function thcond(Tk,Pbars,Dkgm3,alph,betaPa)
!Returns thermal conductivity of H2O in J/m*deg*s (=W/m*deg)
!if  Tk, Pbars  falls within the validity region (specified by the initial if statement)
!  of the Sengers et al. (1984) equation; 
!returns zero otherwise
!See equations 3.2-14 and tables 2-5 and I.5-6 from the above reference.
  implicit real(8) (a-h,o-z)
  parameter (Tstar=  647.270d0)
  parameter (Dstar=  317.763d0)
  parameter (Pstar=  22.1150d6)
  parameter (ustar=  1.0d-6)
  parameter (C=      3.7711d-8)
  real(8)  aL(4), au(4), bL(6,5), bu(5,6), L0, L1, L2
  save
  data aL / 0.2022230d1,  0.1411166d2,  0.5255970d1, -0.2018700d1 /
  data au / 0.0181583d0,  0.0177624d0,  0.0105287d0, -0.0036744d0 /
  data bL / 1.329304600d0, -0.404524370d0,  0.244094900d0, &
            0.018660751d0, -0.129610680d0,  0.044809953d0, &
            1.701836300d0, -2.215684500d0,  1.651105700d0, &
           -0.767360020d0,  0.372833440d0, -0.112031600d0, &
            5.224615800d0, -1.012411100d1,  4.987468700d0, &
           -0.272976940d0, -0.430833930d0,  0.133338490d0, &
            8.712767500d0, -9.500061100d0,  4.378660600d0, &
           -0.917837820d0,  0.0d0,          0.0d0, &
           -1.852599900d0,  0.934046900d0,  0.0d0, &
            0.0d0,          0.0d0,          0.0d0  /
  data bu / 0.5019380d0,  0.2356220d0, -0.2746370d0,  0.1458310d0, &
           -0.0270448d0,  0.1628880d0,  0.7893930d0, -0.7435390d0, &
            0.2631290d0, -0.0253093d0, -0.1303560d0,  0.6736650d0, &
           -0.9594560d0,  0.3472470d0, -0.0267758d0,  0.9079190d0, &
            1.2075520d0, -0.6873430d0,  0.2134860d0, -0.0822904d0, &
           -0.5511190d0,  0.0670665d0, -0.4970890d0,  0.1007540d0, &
            0.0602253d0,  0.1465430d0, -0.0843370d0,  0.1952860d0, &
           -0.0329320d0, -0.0202595d0  /
  data TOL /1.0d-2/
  thcond=  0.0d0
  TdegC=   Tk - 273.15d0
  if ((Pbars > 4000.0d0+TOL) .or. &
     ((Pbars > 2000.0d0+TOL).and.(TdegC > 125.0d0+TOL)).or. &
     ((Pbars > 1500.0d0+TOL).and.(TdegC > 400.0d0+TOL)).or. &
     (TdegC  > 800.0d0+TOL))  return
  T=  Tk / Tstar
  D=  Dkgm3 / Dstar
  sum=  0.0d0
  do i=0,3; sum=  sum + aL(i+1)/T**i; end do
  L0=  sqrt(T) / sum
  sum=  0.0d0
  do i=0,4
    do j=0,5; sum=  sum + bL(j+1,i+1) * (1.0d0/T-1)**i * (D-1)**j; end do
  end do
  L1=  Dexp(D*sum)
  sum=  0.0d0
  do i=0,3; sum=  sum + au(i+1)/T**i; end do
  u0=  ustar * sqrt(T) / sum
  sum=  0.0d0
  do i=0,5
    do j=0,4; sum=  sum + bu(j+1,i+1) * (1.0d0/T-1)**i * (D-1)**j; end do
  end do
  u1=  Dexp(D*sum)
  xt=    Pstar/Dstar**2 * betaPa * Dkgm3**2
  dPdT=  Tstar/Pstar * alph/betaPa
  L2=  C / (u0*u1) * (T/D)**2 * dPdT**2 * power(xt,0.4678d0) &
       * sqrt(D) * Dexp(-18.66d0*(T-1)**2 - (D-1)**4)
  thcond=  L0 * L1 + L2
  return
end function thcond

real(8) function surten(Tsatur)
!Returns the surface tension of vapor-saturated liquid H2O in MPa*cm (converted from N/m)
!as computed from the Vargaftik et al. (1983) equation
!See equations 10.1-2, Kestin et al. (1984); compare also equation C.5 and table 11, Haar et al. (1984).
  implicit real(8) (a-h,o-z)
  parameter (Ttripl=  273.16d0)
  parameter (Tcrit=   647.067d0)
  parameter (Tstar=   647.27d0)
  parameter (Tcstar=  0.999686d0)
  parameter (v=   1.256d0)
  parameter (B=  -0.625d0)
  parameter (stref=  0.2358d0)
  parameter (FPTOL=  1.0d-10)
  save
  if ((Tsatur < Ttripl) .or. (Tsatur > Tcrit)) then
    surten=  0.0d0
    return
  end if
  if (Tsatur>=Tcrit-FPTOL) then; Tnorm=  0.0d0
                           else; Tnorm=  (Tcstar - Tsatur/Tstar) / Tcstar; end if
  surten=  stref * power(Tnorm,v) * (1.0d0 + B*Tnorm)
  return
end function surten

!Born92 - Computes the Z, Q, Y, and X Born functions at TK, Pbars.
!  epseqn=  1 use Helgeson-Kirkham (1974) equation 
!  epseqn=  2 use Pitzer (1983) equation
!  epseqn=  3 use Uematsu-Franck (1980) equation
!  epseqn=  4 use Johnson-Norton (1991) equation
!  epseqn=  5 use Archer-Wang (1990) equation
!                   
subroutine Born92(TK,Pbars,Dgcm3,betab,alphaK,daldT,eps,Z,Q,Y,X,epseqn)
  implicit real(8) (a-h,o-z)
  character,parameter::T_=char(9)
  parameter (TMAX=  1000.0d0, PMAX=  5000.0d0, TOL=  1.0d-3)
  integer  epseqn
  save
  eps=  0.0d0
  Z=    0.0d0
  Y=    0.0d0
  Q=    0.0d0
  X=    0.0d0

  TdegC=  TK - 273.15d0
!The following line can be commented out to facilitate probably unreliable, yet potentially useful, predictive calculations 
!at state conditions beyond the validity limits of the aqueous species equation of state.
  if ((TdegC > TMAX+TOL) .or. (Pbars > PMAX+TOL)) return
! if (epseqn == 1) then
!   call HK74(TK,Dgcm3,betab,alphaK,daldT,eps,dedP,dedT,d2edT2)
!   call epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
!   return
! end if

! if (epseqn == 2) then
!   call Pitz83(TK,Dgcm3,betab,alphaK,daldT,eps,dedP,dedT,d2edT2)
!   call epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
!   return
! end if

! if (epseqn == 3) then
!   call UF80(TK,Dgcm3,betab,alphaK,daldT,eps,dedP,dedT,d2edT2)
!   call epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
!   return
! end if
  if (epseqn==4) then
    call JN91(TK,Dgcm3,betab,alphaK,daldT,eps,dedP,dedT,d2edT2)
    call epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
    return
  end if
! if (epseqn == 5) then
!   Dkgm3=  Dgcm3 * 1.0d3
!   PMPa=   Pbars / 1.0d1
!   betam=  betab * 1.0d1
!   call AW90(TK,PMPa,Dkgm3,betam,alphaK,daldT,eps,dedP,dedT,d2edT2)
!***convert  dedP  FROM  MPa**-1  TO  bars**-1
!   dedP=  dedP / 1.0d1
!   call epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
!   return
! end if
  write(11,'(A)') "BORN92"
  write(11,'(F7.1,A1,E15.8,A1)',advance="no") & 
             TK,  t_,Dgcm3,t_
  write(11,'(E15.8,A1,E15.8,A1,E15.8,A1,E15.8, A1,E15.8,A1,E15.8,A1,E15.8,A1,E15.8,A1)') &
             eps,  t_,dEdP, t_,dEdT, t_,d2EdT2,t_,Q,    t_,X,    t_,Y,    t_,Z,    t_     

end subroutine Born92

subroutine JN91(T,D,beta,alpha,daldT,eps,dedP,dedT,d2edT2)
!Compute (eps, dedP, dedT, d2edT2)(T,D) using equations given by Johnson and Norton (1991); 
!fit parameters regressed from least squares fit to dielectric data consistent
!with the HK74 equation and low temperatures,
!and with the Pitz83 equation at high temperatures.
! Units:
!   T ............... K
!   D ............... g/cm**3
!   beta, dedP ...... bar**(-1)
!   alpha, dedT ..... K**(-1)
!   daldT, d2edT2 ... K**(-2)
  implicit real(8) (a-h,o-z)
  real(8)  a(10), c(5), dcdT(5), dc2dTT(5)
  save
  data Tref / 298.15d0 /
  data a / &
            0.1470333593E+02, & 
            0.2128462733E+03, & 
           -0.1154445173E+03, & 
            0.1955210915E+02, & 
           -0.8330347980E+02, & 
            0.3213240048E+02, & 
           -0.6694098645E+01, & 
           -0.3786202045E+02, & 
            0.6887359646E+02, & 
           -0.2729401652E+02 /
  Tn=  T / Tref
  !
  c(1)=       1.0d0
  dcdT(1)=    0.0d0
  dc2dTT(1)=  0.0d0
  !
  c(2)=       a(1)/Tn
  dcdT(2)=    -a(1)*Tref/T**2
  dc2dTT(2)=  2.0d0*a(1)*Tref/T**3
  !
  c(3)=       a(2)/Tn + a(3) + a(4)*Tn           
  dcdT(3)=    -a(2)*Tref/T**2 + a(4)/Tref           
  dc2dTT(3)=  2.0d0*a(2)*Tref/T**3
  !
  c(4)=       a(5)/Tn + a(6)*Tn + a(7)*Tn**2
  dcdT(4)=    -a(5)*Tref/T**2 + a(6)/Tref + 2.0d0*a(7)*T/Tref**2
  dc2dTT(4)=  2.0d0*a(5)*Tref/T**3 + 2.0d0*a(7)/Tref**2
  !
  c(5)=       a(8)/Tn**2 + a(9)/Tn + a(10)
  dcdT(5)=    -2.0d0*a(8)*Tref**2/T**3 - a(9)*Tref/T**2
  dc2dTT(5)=  6.0d0*a(8)*Tref**2/T**4 + 2.0d0*a(9)*Tref/T**3
  !
  eps=  0.0d0
  do k=1,5
    eps=  eps + c(k)*D**(k-1)
  end do
  !
  dedP=  0.0d0
  do j=  0,4; dedP=  dedP + j*c(j+1)*D**j; end do
  dedP=  beta * dedP
  !
  dedT=  0.0d0
  do j=  0,4; dedT=  dedT + D**j*(dcdT(j+1) - j*alpha*c(j+1)); end do
  d2edT2=  0.0d0
  do j=  0,4
    d2edT2=  d2edT2 &
           + D**j*(dc2dTT(j+1) &
           - j*(alpha*dcdT(j+1) &
           + c(j+1)*daldT) &
           - j*alpha*(dcdT(j+1) &
           - j*alpha*c(j+1)))
  end do
end subroutine JN91

subroutine epsBrn(eps,dedP,dedT,d2edT2,Z,Q,Y,X)
!Compute the Z, Q, Y, and X Born functions from their eps, dedP, dedT, and d2edT2 counterparts.
  implicit real(8) (a-h,o-z)
  save
  Z=  -1.0d0/eps
  Q=   1.0d0/eps**2 * dedP
  Y=   1.0d0/eps**2 * dedT
  X=   1.0d0/eps**2 * d2edT2 - 2.0d0*eps*Y**2
end subroutine epsBrn

end module M_Eos_H2O_Supcrt

