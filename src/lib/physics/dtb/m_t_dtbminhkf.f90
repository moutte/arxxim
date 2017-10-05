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
    real(dp),dimension(3):: Ttran, Htran, Vtran, dPdTtr
    real(dp),dimension(4):: MK1, MK2, MK3, MK4
  end type T_DtbMinHkf
  !
contains

subroutine DtbMinHkf_Calc(M,TdgK,Pbar,S)
!--
!-- computes the standard molal thermodynamic properties
!-- of mineral M at Pbar,TdgK,
!-- using equations from Helgeson et al. (1978).
!--
  use M_Dtb_Const,    only: R_cK,R_jK,Tref,Pref
  use M_Numeric_Const,only: Ln10
  use M_T_Species,    only: T_Species
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
  !---------------------------------------------- Phase Transition terms
  PhasRegion=f_MinHkf_PhaseRegion(M,Pbar,Pref,TdgK)
  do PhTran=1,PhasRegion-1
    Spttrm=  Spttrm + M%Htran(PhTran) /M%Ttran(PhTran)
    Hpttrm=  Hpttrm + M%Htran(PhTran)
    Gpttrm=  Gpttrm + M%Htran(PhTran) *(One - TdgK /M%Ttran(PhTran))
    !!-> when there is phase transition, transition enthalpy is used  ...
  end do
  !---------------------------------------------/ Phase Transition terms
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
  M%MK1(1:4)= Zero
  M%MK2(1:4)= Zero
  M%MK3(1:4)= Zero
end subroutine DtbMinHkf_Zero

function Cp(T,cf)
!.standard molal heat capacity at T.
  real(dp):: Cp,T,cf(:)
  Cp=  cf(1) + cf(2)*T + cf(3)/T/T + cf(4)/sqrt(T)
end function Cp

function CpdT(T1,T2,cf)
!.integral Cp.dT evaluated from T1 to T2.
  real(dp):: CpdT,T1,T2,cf(:)
  CpdT= cf(1)*(T2       - T1      ) &
  &   + cf(2)*(T2*T2    - T1*T1   )/2.0D0 &
  &   - cf(3)*(One/T2   - One/T1  ) &
  &   + cf(4)*(sqrt(T2) - sqrt(T1))*2.0D0
end function CpdT

function CpdlnT(T1,T2,cf)
!.integral Cp.dlnT evaluated from T1 to T2.
  real(dp):: CpdlnT,T1,T2,cf(:)
  CpdlnT= cf(1)*(log(T2)      - log(T1)     ) &
  &     + cf(2)*(T2           - T1          ) &
  &     - cf(3)*(One/T2/T2    - One/T1/T1   )/2.0D0 &
  &     - cf(4)*(One/sqrt(T2) - One/sqrt(T1))*2.0D0
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
!-- Returns phase region for mineral imin at P, T; 
!-- and, as a side effect, TtranP(1..MXTRAN) as f(P).
!-- f_MinHkf_PhaseRegion-  1 ... TtranP(1) > T  [or imin lacks transn]
!-- f_MinHkf_PhaseRegion-  2 ... TtranP(1) - T  [- TtranP(2)]
!-- f_MinHkf_PhaseRegion-  3 ... TtranP(2) - T  [- TtranP(3)]
!-- f_MinHkf_PhaseRegion-  4 ... TtranP(3) - T
  type(T_DtbMinHkf),intent(in)::M
  real(dp),      intent(in)::P,P_ref,T
  !
  real(dp),dimension(3)::TtranP
  !
  f_MinHkf_PhaseRegion=  1 !-- phase region 1
  if (M%nTran==0) return
  if (M%dPdTtr(1)==Zero) then  ;  TtranP(1)=M%Ttran(1)
  else                         ;  TtranP(1)=M%Ttran(1)+(P-P_ref)/M%dPdTtr(1)
  end if
  if (T<=TtranP(1)) return
  !
  f_MinHkf_PhaseRegion=  2 !-- phase region 2
  if (M%nTran==1) return
  if (M%dPdTtr(2)==Zero) then  ;  TtranP(2)=M%Ttran(2)
  else                         ;  TtranP(2)=M%Ttran(2)+(P-P_ref)/M%dPdTtr(2)
  end if
  if (T<=TtranP(2)) return
  !
  f_MinHkf_PhaseRegion=  3 !-- phase region 3
  if (M%nTran== 2)   return
  if (M%dPdTtr(3)==Zero) then  ;  TtranP(3)=M%Ttran(3)
  else                         ;  TtranP(3)=M%Ttran(3)+(P-P_ref)/M%dPdTtr(3)
  end if
  if (T<=TtranP(3)) return
  !
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
  type(T_DtbMinHkf),intent(in) :: M
  real(dp),         intent(in) :: TdgK,T_ref
  real(dp),         intent(out):: Cpr,CprdT,CprdlT 
  !
  !iCpRegion=f_MinHkf_CpRegion(M,TdgK)
  !
  select case(f_MinHkf_CpRegion(M,TdgK)) !iCpRegion)
  case(1)
    Cpr=     Cp(TdgK,          M%MK1(1:4))
    CprdT=   CpdT(T_ref,  TdgK,M%MK1(1:4))
    CprdlT=  CpdlnT(T_ref,TdgK,M%MK1(1:4))
  case(2)
    Cpr=   Cp(TdgK,                 M%MK2(1:4))
    !
    CprdT= CpdT(T_ref,M%Ttran(1),   M%MK1(1:4)) &
    &    + CpdT(M%Ttran(1),TdgK,    M%MK2(1:4))
    !
    CprdlT= CpdlnT(T_ref,M%Ttran(1),M%MK1(1:4)) &
    &     + CpdlnT(M%Ttran(1),TdgK, M%MK2(1:4))
  case(3)
    Cpr=     Cp(TdgK, M%MK3(1:4))
    !
    CprdT= CpdT(T_ref,M%Ttran(1),       M%MK1(1:4)) &
    &    + CpdT(M%Ttran(1),M%Ttran(2),  M%MK2(1:4)) &
    &    + CpdT(M%Ttran(2),TdgK,        M%MK3(1:4))
    !
    CprdlT= CpdlnT(T_ref,M%Ttran(1),     M%MK1(1:4)) &
    &     + CpdlnT(M%Ttran(1),M%Ttran(2),M%MK2(1:4)) &
    &     + CpdlnT(M%Ttran(2),TdgK,      M%MK3(1:4))
  case(4)
    Cpr= Cp(TdgK, M%MK4(1:4))
    !
    CprdT= CpdT(T_ref,     M%Ttran(1),M%MK1(1:4)) &
    &    + CpdT(M%Ttran(1),M%Ttran(2),M%MK2(1:4)) &
    &    + CpdT(M%Ttran(2),M%Ttran(3),M%MK3(1:4)) &
    &    + CpdT(M%Ttran(3),TdgK,      M%MK4(1:4)) 
    !
    CprdlT= CpdlnT(T_ref,     M%Ttran(1),M%MK1(1:4)) &
    &     + CpdlnT(M%Ttran(1),M%Ttran(2),M%MK2(1:4)) &
    &     + CpdlnT(M%Ttran(2),M%Ttran(3),M%MK3(1:4)) &
    &     + CpdlnT(M%Ttran(3),TdgK,      M%MK4(1:4))
  end select
  !
  return
end subroutine MinHkf_Cp_Calc

end module M_T_DtbMinHkf

