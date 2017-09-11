module M_T_DtbMinHkf
!--
!-- module for data structures dedicated to thermodynamic data on min./gas species 
!-- Hkf approach (apparent G a la Benson-Helgeson, minerals have constant volumes, etc.)
!--
  use M_Kinds
  implicit none
  !
  private
  !
  public:: T_DtbMinHkf
  public:: DtbMinHkf_Calc
  public:: DtbMinHkf_Zero
  !
  type:: T_DtbMinHkf
  !--
  !-- data structure for thermo data for Min.Species according to SupCrt database
  !--
    character(len=15):: Num
    character(len=23):: Name
    character(len=5) :: Abbr
    character(len=71):: Formula
    character(len=3) :: Typ !MIN or GAS <-not used yet !!!??
    integer :: Div=1 !!,nO_=0
    integer :: nTran, phaser
    real(dp):: &
    & G0R,     &
    & H0R,     &
    & S0_,     &
    & V0R,     &
    & Tmax,    &
    & WeitKg,  &
    & S0Ele !ref. state entropy
    real(dp),dimension(3):: Ttran, Htran, Vtran
    real(dp),dimension(3):: MK1, MK2, MK3, MK4, dPdTtr !MK=MaierKelly Coeff's
  end type T_DtbMinHkf
  !
contains

subroutine DtbMinHkf_Calc(M,TdgK,Pbar,S)
!--
!-- computes the standard molal thermodynamic properties of mineral M at Pbar,TdgK
!-- using equations from Helgeson et al. (1978).
!--
  use M_Dtb_Const,only: CalToJoule,R_cK,R_jK,Tref,Pref
  use M_Numeric_Const,    only: Ln10
  use M_T_Species,only: T_Species
  !
  type(T_DtbMinHkf),intent(in)   :: M !thermodyn' model
  real(dp),         intent(in)   :: TdgK,Pbar
  type(T_Species),  intent(inout):: S !species
  !
  integer :: PhTran,PhasRegion
  real(dp):: CprdT,CprdlT,Spttrm,Hpttrm,Gpttrm,VdP
  real(dp):: GR, HR, SR, VR, CpR !-> stores values computed for given P,T
  !
  Spttrm=  Zero
  Hpttrm=  Zero
  Gpttrm=  Zero
  !--------------------------------------------------- Phase Transition terms --
  PhasRegion=f_MinHkf_PhaseRegion(M,Pbar,Pref,TdgK)
  do PhTran=1,PhasRegion-1
    Spttrm=  Spttrm + M%Htran(PhTran) /M%Ttran(PhTran)
    Hpttrm=  Hpttrm + M%Htran(PhTran)
    Gpttrm=  Gpttrm + M%Htran(PhTran) *(One - TdgK /M%Ttran(PhTran))
    !!-> when there is phase transition, transition enthalpy is used  ...
  end do
  !--------------------------------------------------/ Phase Transition terms --
  !
  select case(trim(M%Typ))
  
  case("MIN")
    call MinHkf_Vol_Calc( &
    & M,Pbar,Pref,TdgK,PhasRegion, &
    & VdP,Vr) ! VdP in calorie, Vr in cm^3
    !
    !--- from bar/cm^3 to Joule
    !--- Pascal*m^3-Joule -> (1D-5 bar)*(1D+6 cm^3) - Joule
    !--- -> Joule - bar*cm^3 /10
    VdP= VdP /10.0D0
    !--- from cm^3 to m^3
    Vr= Vr*1.0D-6
  
  case("GAS") ! ideal gas
    Vdp= R_jk *TdgK *log(Pbar /Pref) !=Joule
    Vr=  R_jK *TdgK /Pbar /1.D5 ! Joule/Pascal= m^3
  
  end select
  !
  call MinHkf_Cp_Calc( &
  & M,TdgK,Tref, &    ! in
  & Cpr,CprdT,CprdlT) ! out
  !
  Sr=   M%S0_ + Spttrm + CprdlT 
  Hr=   M%H0R + Hpttrm + CprdT  + VdP
  Gr=   M%G0R + Gpttrm + CprdT  + VdP - M%S0_*(TdgK-Tref) - TdgK*CprdlT
  !-> here, enthalpy Hr is not used for calculation of Gr !!!
  !   (a difference with calcul'n on Thr-like data)
  !
  S%WeitKg= M%WeitKg
  !
  S%G0rt=  GR /R_jk /TdgK
  S%H0=    HR
  S%S0=    SR
  S%V0=    VR
  S%Cp0=   CpR
  S%LnFug= Zero
  !
  !!S%VMol= S%V0
  !! Rho=  S%WeitKg / S%V0 !-> density of mineral, in kg/m3
end subroutine DtbMinHkf_Calc

subroutine DtbMinHkf_Zero(M)
  type(T_DtbMinHkf),intent(out)::M
  !
  M%Div=   1
  M%G0R=   Zero
  M%H0R=   Zero
  M%S0_=   Zero
  M%V0R=   Zero
  M%S0Ele= Zero
  !
  M%MK1(1:3)= Zero
  M%MK2(1:3)= Zero
  M%MK3(1:3)= Zero
end subroutine DtbMinHkf_Zero

!! integer function DtbMinHkf_Index(Str,V) 
!! !.-> index of mineral named Str in vDtbMinHkf
  !! character(*),                  intent(in):: Str
  !! type(T_DtbMinHkf),dimension(:),intent(in):: V
  !! integer     ::I
  !! DtbMinHkf_Index=0
  !! I=0
  !! do
    !! I=I+1 !; if(idebug>1) write(fTrc,'(A)') vEle(I)%SpName
    !! if(trim(Str)==trim(V(I)%Name)) then;
      !! DtbMinHkf_Index=I
      !! exit
    !! end if
    !! if(I==size(V)) exit
  !! end do !if Str not found -> DtbMinThr_Index=0
!! end function DtbMinHkf_Index

function Cp(T,a,b,c)
!.standard molal heat capacity at T.
  real(dp)::Cp,T,a,b,c
  Cp=  a + b*T + c/T/T
end function Cp

function CpdT(T1,T2,a,b,c)
!.integral Cp.dT evaluated from T1 to T2.
  real(dp)::CpdT,T1,T2,a,b,c
  CpdT= a*(T2     - T1    ) &
  &   + b*(T2*T2  - T1*T1 )/2.0d0 &
  &   - c*(One/T2 - One/T1)
end function CpdT

function CpdlnT(T1,T2,a,b,c)
!.integral Cp.dlnT evaluated from T1 to T2.
  real(dp)::CpdlnT,T1,T2,a,b,c
  CpdlnT= a*(log(T2)   - log(T1)  ) &
  &     + b*(T2        - T1       ) &
  &     - c*(One/T2/T2 - One/T1/T1)/2.0D0
end function CpdlnT

subroutine MinHkf_Vol_Calc( &
& M,P,P_ref,T,PhaseRegion, &
& VdP,Vr)
!--
!-- computes VMin_(P,T), VMin_*dP_, and (if necesary) PtranT_.
!--
  use M_Dtb_Const, only: CalToJoule
  !
  type(T_DtbMinHkf),intent(in) :: M
  real(dp),         intent(in) :: P,P_ref ! bar
  real(dp),         intent(in) :: T       ! dgK
  integer,          intent(in) :: PhaseRegion
  real(dp),         intent(out):: VdP     ! bar.cm^3
  real(dp),         intent(out):: Vr      ! cm^3
  !
  real(dp),dimension(2)::PtranT
  integer :: I
  !
  Vr= M%V0R ! cm^3
  do i=1,PhaseRegion-1 ; Vr=Vr + M%Vtran(i) ; end do
  !
  VdP= Vr*(P - P_ref) !-> VdP in cm^3*bar
  !
  !return if Pressure integration does not cross phase transition boundaries
  if ((M%nTran    ==0) &
  .or.(M%dPdTtr(1)==Zero) &
  .or.(M%Ttran(1) >=T)) return 
  ! 
  if ((M%nTran==1).and.(PhaseRegion==2))  return
  if ((M%nTran==2).and.(PhaseRegion==3))  return
  if ((M%nTran==2).and.(PhaseRegion==2).and.(M%Ttran(2)>T)) return
  !
  !-- take account of cross-boundary pressure integration 
  if ((M%nTran    == 1) &
 .or.((PhaseRegion==1).and.(M%Ttran(2)>T))) then
    PtranT(1)=  P_ref + (T - M%Ttran(1))*M%dPdTtr(1)
    VdP= M%V0R     *(P - P_ref)  &
    &  + M%Vtran(1)*(PtranT(1) - P_ref) !-> VdP in cm^3*bar
    return
  end if
  !
  !-- nTran-  2 and T >- Ttran(2)
  PtranT(2)=  P_ref + (T - M%Ttran(2))*M%dPdTtr(2)
  if (PhaseRegion==2) then
    VdP=(M%V0R + M%Vtran(1))*(P -         P_ref) &
    &  + M%Vtran(2)         *(PtranT(2) - P_ref)
  else
    PtranT(1)= P_ref + (T - M%Ttran(1))*M%dPdTtr(1)
    VdP= M%V0R     *(P - P_ref) &
    &  + M%Vtran(1)*(PtranT(1) - P_ref) &
    &  + M%Vtran(2)*(PtranT(2) - P_ref)
  end if
  return
end subroutine MinHkf_Vol_Calc

integer function f_MinHkf_PhaseRegion(M,P,P_ref,T)
  type(T_DtbMinHkf),intent(in)::M
  real(dp),      intent(in)::P,P_ref,T
  !
  !-- Returns phase region for mineral imin at P, T; 
  !-- and, as a side effect, TtranP(1..MXTRAN) as f(P).
  !-- f_MinHkf_PhaseRegion-  1 ... TtranP(1) > T  [or imin lacks transn]
  !-- f_MinHkf_PhaseRegion-  2 ... TtranP(1) - T  [- TtranP(2)]
  !-- f_MinHkf_PhaseRegion-  3 ... TtranP(2) - T  [- TtranP(3)]
  !-- f_MinHkf_PhaseRegion-  4 ... TtranP(3) - T
  !
  real(dp),dimension(3)::TtranP
  !
  f_MinHkf_PhaseRegion=  1 !-- phase region 1
  if (M%nTran==0) return
  if (M%dPdTtr(1)==Zero) then  ;  TtranP(1)=M%Ttran(1)
  else                         ;  TtranP(1)=M%Ttran(1)+(P-P_ref)/M%dPdTtr(1)
  end if
  if (T<=TtranP(1)) return
  f_MinHkf_PhaseRegion=  2 !-- phase region 2
  if (M%nTran==1) return
  if (M%dPdTtr(2)==Zero) then  ;  TtranP(2)=M%Ttran(2)
  else                         ;  TtranP(2)=M%Ttran(2)+(P-P_ref)/M%dPdTtr(2)
  end if
  if (T<=TtranP(2)) return
  f_MinHkf_PhaseRegion=  3 !-- phase region 3
  if (M%nTran== 2)   return
  if (M%dPdTtr(3)==Zero) then  ;  TtranP(3)=M%Ttran(3)
  else                         ;  TtranP(3)=M%Ttran(3)+(P-P_ref)/M%dPdTtr(3)
  end if
  if (T<=TtranP(3)) return
  f_MinHkf_PhaseRegion=  4 !-- phase region 4
  return
end function f_MinHkf_PhaseRegion

integer function f_MinHkf_CpRegion(M,T) !in M_T_DtbMinHkf
!.-> effective phase region for temperature integration of Cpr(T)
!.for mineral imin (i.e., the phase region specified by T at 1 bar).
  type(T_DtbMinHkf)::M
  real(dp)::T
  integer ::I
  !
  f_MinHkf_CpRegion=1
  do i=1,M%nTran
    if (T > M%TTran(i)) f_MinHkf_CpRegion=f_MinHkf_CpRegion+1
  end do
end function f_MinHkf_CpRegion

subroutine MinHkf_Cp_Calc(M,TdgK,T_ref,Cpr,CprdT,CprdlT)
!--
!-- computes the standard molal heat capacity and heat capacity temperature integrals, 
!-- evaluated from Tref to TdgK at 1 bar.
!--
  !
  type(T_DtbMinHkf),intent(in) :: M
  real(dp),         intent(in) :: TdgK,T_ref
  real(dp),         intent(out):: Cpr,CprdT,CprdlT 
  !
  !iCpRegion=f_MinHkf_CpRegion(M,TdgK)
  !
  select case(f_MinHkf_CpRegion(M,TdgK)) !iCpRegion)
  case(1)
    Cpr=     Cp(TdgK,          M%MK1(1),M%MK1(2),M%MK1(3))
    CprdT=   CpdT(T_ref,  TdgK,M%MK1(1),M%MK1(2),M%MK1(3))
    CprdlT=  CpdlnT(T_ref,TdgK,M%MK1(1),M%MK1(2),M%MK1(3))
  case(2)
    Cpr=   Cp(TdgK,                 M%MK2(1),M%MK2(2),M%MK2(3))
    !
    CprdT= CpdT(T_ref,M%Ttran(1),   M%MK1(1),M%MK1(2),M%MK1(3)) &
    &    + CpdT(M%Ttran(1),TdgK,    M%MK2(1),M%MK2(2),M%MK2(3))
    !
    CprdlT= CpdlnT(T_ref,M%Ttran(1),M%MK1(1),M%MK1(2),M%MK1(3)) &
    &     + CpdlnT(M%Ttran(1),TdgK,M%MK2(1),M%MK2(2),M%MK2(3))
  case(3)
    Cpr=     Cp(TdgK, M%MK3(1),M%MK3(2),M%MK3(3))
    !
    CprdT= CpdT(T_ref,M%Ttran(1),       M%MK1(1),M%MK1(2),M%MK1(3)) &
    &    + CpdT(M%Ttran(1),M%Ttran(2),  M%MK2(1),M%MK2(2),M%MK2(3)) &
    &    + CpdT(M%Ttran(2),TdgK,        M%MK3(1),M%MK3(2),M%MK3(3))
    !
    CprdlT= CpdlnT(T_ref,M%Ttran(1),     M%MK1(1),M%MK1(2),M%MK1(3)) &
    &     + CpdlnT(M%Ttran(1),M%Ttran(2),M%MK2(1),M%MK2(2),M%MK2(3)) &
    &     + CpdlnT(M%Ttran(2),TdgK,      M%MK3(1),M%MK3(2),M%MK3(3))
  case(4)
    Cpr= Cp(  TdgK,                      M%MK4(1),M%MK4(2),M%MK4(3))
    !
    CprdT= CpdT(T_ref,     M%Ttran(1),M%MK1(1),M%MK1(2),M%MK1(3)) &
    &    + CpdT(M%Ttran(1),M%Ttran(2),M%MK2(1),M%MK2(2),M%MK2(3)) &
    &    + CpdT(M%Ttran(2),M%Ttran(3),M%MK3(1),M%MK3(2),M%MK3(3)) &
    &    + CpdT(M%Ttran(3),TdgK,      M%MK4(1),M%MK4(2),M%MK4(3)) 
    !
    CprdlT= CpdlnT(T_ref,     M%Ttran(1),M%MK1(1),M%MK1(2),M%MK1(3)) &
    &     + CpdlnT(M%Ttran(1),M%Ttran(2),M%MK2(1),M%MK2(2),M%MK2(3)) &
    &     + CpdlnT(M%Ttran(2),M%Ttran(3),M%MK3(1),M%MK3(2),M%MK3(3)) &
    &     + CpdlnT(M%Ttran(3),TdgK,      M%MK4(1),M%MK4(2),M%MK4(3))
  end select
  !
  return
end subroutine MinHkf_Cp_Calc

!!!   subroutine CalcGasHKF(s_,TdgK) !gas-> calculate G,H,S at (TdgK, Pref)
!!!     use M_Dtb_Const,only:CalToJoule,R_jk,Tref
!!!     use M_Numeric_Const,only:Ln10
!!!     type(T_GasHkf),intent(inout)::S_
!!!     real(dp),      intent(in)   ::TdgK
!!!     !
!!!     integer::nEl
!!!     real(dp)::Cpr,CprdT,CprdlT
!!!     !Cptrms(phase,i,CpRegion,T,Cpr,CprdT,CprdlT)
!!!     !Computes the standard molal heat capacity and heat capacity temperature integrals, 
!!!     !evaluated from Tref to T at 1 bar.
!!!     Cpr=       Cp    (TdgK,     S_%MK(1),S_%MK(2),S_%MK(3))
!!!     CprdT=     CpdT  (Tref,TdgK,S_%MK(1),S_%MK(2),S_%MK(3))
!!!     CprdlT=    CpdlnT(Tref,TdgK,S_%MK(1),S_%MK(2),S_%MK(3))
!!!     !H(Pr,T)=  Hf(Pr,Tr) + Intg<Tr,T>(Cp.dT)
!!!     !S(Pr,T)=  S(Pr,Tr)  + Intg<Tr,T>(Cp.dT/T)
!!!     !G(Pr,T)=  Gf(Pr,Tr) - S(Pr,Tr).(T-Tr) + Intg<Tr,T>(Cp.dT) - T.Intg<Tr,T>(Cp.dT/T)
!!!     S_%Vr=   S_%V0R
!!!     S_%Sr=   S_%S0_ + CprdlT
!!!     S_%Hr=   S_%H0R  + CprdT
!!!     S_%Gr=   S_%G0R  - S_%S0_*(TdgK - Tref) + CprdT - TdgK*CprdlT
!!!     !conversions Cal -> Joule
!!!     S_%Gr=   S_%Gr*CalToJoule
!!!     S_%Hr=   S_%Hr*CalToJoule
!!!     S_%Sr=   S_%Sr*CalToJoule
!!!     !
!!!     nEl=size(vEle)
!!!     S_%S0Ele= dot_product(S_%Stoik(1:nEl),vEle(1:nEl)%S0) /S_%Div
!!!     S_%WeitKg=dot_product(S_%Stoik(1:nEl)/S_%Div,vEle(1:nEl)%WeitKg) /S_%Div
!!!     !
!!!     S_%logK=-S_%Gr/R_jk/TdgK/Ln10
!!!   end subroutine CalcGasHKF

end module M_T_DtbMinHkf

