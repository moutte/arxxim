MODULE M_T_DtbMinHkf
!--
!-- module for data structures dedicated to thermodynamic data on min./gas species 
!-- Hkf approach (apparent G a la Benson-Helgeson, minerals have constant volumes, etc.)
!--
  USE M_Kinds
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: T_DtbMinHkf
  PUBLIC:: DtbMinHkf_Calc
  PUBLIC:: DtbMinHkf_Zero
  !
  TYPE:: T_DtbMinHkf
  !--
  !-- data structure for thermo data for Min.Species according to SupCrt database
  !--
    CHARACTER(LEN=15):: Num
    CHARACTER(LEN=23):: Name
    CHARACTER(LEN=71):: Formula
    CHARACTER(LEN=3) :: Typ !MIN or GAS <-not used yet !!!??
    INTEGER :: Div=1 !!,nO_=0
    INTEGER :: nTran, phaser
    REAL(dp):: &
    & G0R,     &
    & H0R,     &
    & S0_,     &
    & V0R,     &
    & Tmax,    &
    & WeitKg,  &
    & S0Ele !ref. state entropy
    REAL(dp),DIMENSION(3):: Ttran, Htran, Vtran
    REAL(dp),DIMENSION(3):: MK1, MK2, MK3, MK4, dPdTtr !MK=MaierKelly Coeff's
  ENDTYPE T_DtbMinHkf
  !
CONTAINS

SUBROUTINE DtbMinHkf_Calc(M,TdgK,Pbar,S)
!--
!-- computes the standard molal thermodynamic properties of mineral M at Pbar,TdgK
!-- using equations from Helgeson et al. (1978).
!--
  USE M_Dtb_Const,ONLY: CalToJoule,R_cK,R_jK,Tref,Pref
  USE M_Numeric_Const,    ONLY: Ln10
  USE M_T_Species,ONLY: T_Species
  !
  TYPE(T_DtbMinHkf),INTENT(IN)   :: M !thermodyn' model
  REAL(dp),         INTENT(IN)   :: TdgK,Pbar
  TYPE(T_Species),  INTENT(INOUT):: S !species
  !
  INTEGER :: PhTran,PhasRegion
  REAL(dp):: CprdT,CprdlT,Spttrm,Hpttrm,Gpttrm,VdP
  REAL(dp):: GR, HR, SR, VR, CpR !-> stores values computed for given P,T
  !
  Spttrm=  Zero
  Hpttrm=  Zero
  Gpttrm=  Zero
  !--------------------------------------------------- Phase Transition terms --
  PhasRegion=f_MinHkf_PhaseRegion(M,Pbar,Pref,TdgK)
  DO PhTran=1,PhasRegion-1
    Spttrm=  Spttrm + M%Htran(PhTran) /M%Ttran(PhTran)
    Hpttrm=  Hpttrm + M%Htran(PhTran)
    Gpttrm=  Gpttrm + M%Htran(PhTran) *(One - TdgK /M%Ttran(PhTran))
    !!-> when there is phase transition, transition enthalpy is USEd  ...
  ENDDO
  !--------------------------------------------------/ Phase Transition terms --
  !
  SELECT CASE(TRIM(M%Typ))
  
  CASE("MIN")
    CALL MinHkf_Vol_Calc( &
    & M,Pbar,Pref,TdgK,PhasRegion, &
    & VdP,Vr) ! VdP in calorie, Vr in cm^3
    !
    !--- from bar/cm^3 to Joule
    !--- Pascal*m^3-Joule -> (1D-5 bar)*(1D+6 cm^3) - Joule
    !--- -> Joule - bar*cm^3 /10
    VdP= VdP /10.0D0
    !--- from cm^3 to m^3
    Vr= Vr*1.0D-6
  
  CASE("GAS") ! ideal gas
    Vdp= R_jk *TdgK *LOG(Pbar /Pref) !=Joule
    Vr=  R_jK *TdgK /Pbar /1.D5 ! Joule/Pascal= m^3
  
  END SELECT
  !
  CALL MinHkf_Cp_Calc( &
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
ENDSUBROUTINE DtbMinHkf_Calc

SUBROUTINE DtbMinHkf_Zero(M)
  TYPE(T_DtbMinHkf),INTENT(OUT)::M
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
ENDSUBROUTINE DtbMinHkf_Zero

!~ INTEGER FUNCTION DtbMinHkf_Index(Str,V) 
!~ !.-> index of mineral named Str in vDtbMinHkf
  !~ CHARACTER(*),                  INTENT(IN):: Str
  !~ TYPE(T_DtbMinHkf),DIMENSION(:),INTENT(IN):: V
  !~ INTEGER     ::I
  !~ DtbMinHkf_Index=0
  !~ I=0
  !~ DO
    !~ I=I+1 !; IF(iDebug>0) WRITE(fTrc,'(A)') vEle(I)%SpName
    !~ IF(TRIM(Str)==TRIM(V(I)%Name)) THEN;
      !~ DtbMinHkf_Index=I
      !~ EXIT
    !~ ENDIF
    !~ IF(I==SIZE(V)) EXIT
  !~ ENDDO !IF Str not found -> DtbMinThr_Index=0
!~ ENDFUNCTION DtbMinHkf_Index

FUNCTION Cp(T,a,b,c)
!.standard molal heat capacity at T.
  REAL(dp)::Cp,T,a,b,c
  Cp=  a + b*T + c/T/T
ENDFUNCTION Cp

FUNCTION CpdT(T1,T2,a,b,c)
!.integral Cp.dT evaluated from T1 to T2.
  REAL(dp)::CpdT,T1,T2,a,b,c
  CpdT= a*(T2     - T1    ) &
  &   + b*(T2*T2  - T1*T1 )/2.0d0 &
  &   - c*(One/T2 - One/T1)
END FUNCTION CpdT

FUNCTION CpdlnT(T1,T2,a,b,c)
!.integral Cp.dlnT evaluated from T1 to T2.
  REAL(dp)::CpdlnT,T1,T2,a,b,c
  CpdlnT= a*(LOG(T2)   - LOG(T1)  ) &
  &     + b*(T2        - T1       ) &
  &     - c*(One/T2/T2 - One/T1/T1)/2.0D0
ENDFUNCTION CpdlnT

SUBROUTINE MinHkf_Vol_Calc( &
& M,P,P_ref,T,PhaseRegion, &
& VdP,Vr)
!--
!-- computes VMin_(P,T), VMin_*dP_, and (IF necesary) PtranT_.
!--
  USE M_Dtb_Const, ONLY: CalToJoule
  !
  TYPE(T_DtbMinHkf),INTENT(IN) :: M
  REAL(dp),         INTENT(IN) :: P,P_ref ! bar
  REAL(dp),         INTENT(IN) :: T       ! dgK
  INTEGER,          INTENT(IN) :: PhaseRegion
  REAL(dp),         INTENT(OUT):: VdP     ! bar.cm^3
  REAL(dp),         INTENT(OUT):: Vr      ! cm^3
  !
  REAL(dp),DIMENSION(2)::PtranT
  INTEGER :: I
  !
  Vr= M%V0R ! cm^3
  DO i=1,PhaseRegion-1 ; Vr=Vr + M%Vtran(i) ; ENDDO
  !
  VdP= Vr*(P - P_ref) !-> VdP in cm^3*bar
  !
  !RETURN IF Pressure integration does not cross phase transition boundaries
  IF ((M%nTran    ==0) &
  .OR.(M%dPdTtr(1)==Zero) &
  .OR.(M%Ttran(1) >=T)) RETURN 
  ! 
  IF ((M%nTran==1).AND.(PhaseRegion==2))  RETURN
  IF ((M%nTran==2).AND.(PhaseRegion==3))  RETURN
  IF ((M%nTran==2).AND.(PhaseRegion==2).AND.(M%Ttran(2)>T)) RETURN
  !
  !-- take account of cross-boundary pressure integration 
  IF ((M%nTran    == 1) &
 .OR.((PhaseRegion==1).AND.(M%Ttran(2)>T))) THEN
    PtranT(1)=  P_ref + (T - M%Ttran(1))*M%dPdTtr(1)
    VdP= M%V0R     *(P - P_ref)  &
    &  + M%Vtran(1)*(PtranT(1) - P_ref) !-> VdP in cm^3*bar
    RETURN
  ENDIF
  !
  !-- nTran-  2 and T >- Ttran(2)
  PtranT(2)=  P_ref + (T - M%Ttran(2))*M%dPdTtr(2)
  IF (PhaseRegion==2) THEN
    VdP=(M%V0R + M%Vtran(1))*(P -         P_ref) &
    &  + M%Vtran(2)         *(PtranT(2) - P_ref)
  ELSE
    PtranT(1)= P_ref + (T - M%Ttran(1))*M%dPdTtr(1)
    VdP= M%V0R     *(P - P_ref) &
    &  + M%Vtran(1)*(PtranT(1) - P_ref) &
    &  + M%Vtran(2)*(PtranT(2) - P_ref)
  ENDIF
  RETURN
ENDSUBROUTINE MinHkf_Vol_Calc

INTEGER FUNCTION f_MinHkf_PhaseRegion(M,P,P_ref,T)
  TYPE(T_DtbMinHkf),INTENT(IN)::M
  REAL(dp),      INTENT(IN)::P,P_ref,T
  !
  !-- Returns phase region for mineral imin at P, T; 
  !-- and, as a side effect, TtranP(1..MXTRAN) as f(P).
  !-- f_MinHkf_PhaseRegion-  1 ... TtranP(1) > T  [or imin lacks transn]
  !-- f_MinHkf_PhaseRegion-  2 ... TtranP(1) - T  [- TtranP(2)]
  !-- f_MinHkf_PhaseRegion-  3 ... TtranP(2) - T  [- TtranP(3)]
  !-- f_MinHkf_PhaseRegion-  4 ... TtranP(3) - T
  !
  REAL(dp),DIMENSION(3)::TtranP
  !
  f_MinHkf_PhaseRegion=  1 !== phase region 1
  IF (M%nTran==0) RETURN
  IF (M%dPdTtr(1)==Zero) THEN  ;  TtranP(1)=M%Ttran(1)
  ELSE                         ;  TtranP(1)=M%Ttran(1)+(P-P_ref)/M%dPdTtr(1)
  END IF
  IF (T<=TtranP(1)) RETURN
  f_MinHkf_PhaseRegion=  2 !== phase region 2
  IF (M%nTran==1) RETURN
  IF (M%dPdTtr(2)==Zero) THEN  ;  TtranP(2)=M%Ttran(2)
  ELSE                         ;  TtranP(2)=M%Ttran(2)+(P-P_ref)/M%dPdTtr(2)
  END IF
  IF (T<=TtranP(2)) RETURN
  f_MinHkf_PhaseRegion=  3 !== phase region 3
  IF (M%nTran== 2)   RETURN
  IF (M%dPdTtr(3)==Zero) THEN  ;  TtranP(3)=M%Ttran(3)
  ELSE                         ;  TtranP(3)=M%Ttran(3)+(P-P_ref)/M%dPdTtr(3)
  END IF
  IF (T<=TtranP(3)) RETURN
  f_MinHkf_PhaseRegion=  4 !== phase region 4
  RETURN
ENDFUNCTION f_MinHkf_PhaseRegion

INTEGER FUNCTION f_MinHkf_CpRegion(M,T) !in M_T_DtbMinHkf
!.-> effective phase region for temperature integration of Cpr(T)
!.for mineral imin (i.e., the phase region specIFied by T at 1 bar).
  TYPE(T_DtbMinHkf)::M
  REAL(dp)::T
  INTEGER ::I
  !
  f_MinHkf_CpRegion=1
  DO i=1,M%nTran
    IF (T > M%TTran(i)) f_MinHkf_CpRegion=f_MinHkf_CpRegion+1
  ENDDO
ENDFUNCTION f_MinHkf_CpRegion

SUBROUTINE MinHkf_Cp_Calc(M,TdgK,T_ref,Cpr,CprdT,CprdlT)
!--
!-- computes the standard molal heat capacity and heat capacity temperature integrals, 
!-- evaluated from Tref to TdgK at 1 bar.
!--
  !
  TYPE(T_DtbMinHkf),INTENT(IN) :: M
  REAL(dp),         INTENT(IN) :: TdgK,T_ref
  REAL(dp),         INTENT(OUT):: Cpr,CprdT,CprdlT 
  !
  !iCpRegion=f_MinHkf_CpRegion(M,TdgK)
  !
  SELECT CASE(f_MinHkf_CpRegion(M,TdgK)) !iCpRegion)
  CASE(1)
    Cpr=     Cp(TdgK,          M%MK1(1),M%MK1(2),M%MK1(3))
    CprdT=   CpdT(T_ref,  TdgK,M%MK1(1),M%MK1(2),M%MK1(3))
    CprdlT=  CpdlnT(T_ref,TdgK,M%MK1(1),M%MK1(2),M%MK1(3))
  CASE(2)
    Cpr=   Cp(TdgK,                 M%MK2(1),M%MK2(2),M%MK2(3))
    !
    CprdT= CpdT(T_ref,M%Ttran(1),   M%MK1(1),M%MK1(2),M%MK1(3)) &
    &    + CpdT(M%Ttran(1),TdgK,    M%MK2(1),M%MK2(2),M%MK2(3))
    !
    CprdlT= CpdlnT(T_ref,M%Ttran(1),M%MK1(1),M%MK1(2),M%MK1(3)) &
    &     + CpdlnT(M%Ttran(1),TdgK,M%MK2(1),M%MK2(2),M%MK2(3))
  CASE(3)
    Cpr=     Cp(TdgK, M%MK3(1),M%MK3(2),M%MK3(3))
    !
    CprdT= CpdT(T_ref,M%Ttran(1),       M%MK1(1),M%MK1(2),M%MK1(3)) &
    &    + CpdT(M%Ttran(1),M%Ttran(2),  M%MK2(1),M%MK2(2),M%MK2(3)) &
    &    + CpdT(M%Ttran(2),TdgK,        M%MK3(1),M%MK3(2),M%MK3(3))
    !
    CprdlT= CpdlnT(T_ref,M%Ttran(1),     M%MK1(1),M%MK1(2),M%MK1(3)) &
    &     + CpdlnT(M%Ttran(1),M%Ttran(2),M%MK2(1),M%MK2(2),M%MK2(3)) &
    &     + CpdlnT(M%Ttran(2),TdgK,      M%MK3(1),M%MK3(2),M%MK3(3))
  CASE(4)
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
  END SELECT
  !
  RETURN
ENDSUBROUTINE MinHkf_Cp_Calc

!!!   SUBROUTINE CalcGasHKF(s_,TdgK) !gas-> calculate G,H,S at (TdgK, Pref)
!!!     USE M_Dtb_Const,ONLY:CalToJoule,R_jk,Tref
!!!     USE M_Numeric_Const,ONLY:Ln10
!!!     TYPE(T_GasHkf),INTENT(INOUT)::S_
!!!     REAL(dp),      INTENT(IN)   ::TdgK
!!!     !
!!!     INTEGER::nEl
!!!     REAL(dp)::Cpr,CprdT,CprdlT
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
!!!     nEl=SIZE(vEle)
!!!     S_%S0Ele= DOT_PRODUCT(S_%Stoik(1:nEl),vEle(1:nEl)%S0) /S_%Div
!!!     S_%WeitKg=DOT_PRODUCT(S_%Stoik(1:nEl)/S_%Div,vEle(1:nEl)%WeitKg) /S_%Div
!!!     !
!!!     S_%logK=-S_%Gr/R_jk/TdgK/Ln10
!!!   ENDSUBROUTINE CalcGasHKF

ENDMODULE M_T_DtbMinHkf

