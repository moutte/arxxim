module M_T_DtbMinThr
!--
!-- definition of structured data type, and related functions,
!-- for minerals tabulated according to CdeCapitani / Theriak (Berman) database
!--
  use M_Kinds
  use M_Trace,only: fTrc,iDebug,T_
  !
  implicit none
  !
  private
  !
  public:: T_DtbMinThr
  public:: DtbMinThr_Index
  public:: DtbMinThr_Zero 
  public:: DtbMinThr_Calc
  !
  type:: T_DtbMinThr
    character(len=15):: Num
    character(len=23):: Name
    character(len=71):: Formula
    character(len=3) :: Typ ! MIN,GAS
    character(len=15):: Special
    !
    integer :: Div= 1
    !
    character(len=7) :: CodGas !used for gases: "REDKWON","VDW",...
    !
    integer :: nLanda    ! number of landau transitions
    integer :: codVol    ! code for volume function
    real(dp)::&          ! for landau transitions
    & TQ1B(1:3),TRE(1:3), ASPK(1:3),BSPK(1:3),DHTR(1:3),&
    & TEQ(1:3), DVTR(1:3),DVDT(1:3),DVDP(1:3)
    !
    real(dp):: G0R,H0R,S0_,V0R
    real(dp):: S0Ele,WeitKg
    
    real(dp):: AA0,AAT,BB0,BBT      !for VDW/RDK/... cubic EoS of Gases (VanDerWaals,RedlichKwong,...)
    real(dp):: Tcrit,Pcrit,Acentric !for SRK/PRS/... cubic EoS of gases
    !
    real(dp):: K1,K2,K3,K4,K5,K6,K7,K8,K9  !coeff's for Cp formulas
    real(dp):: D1,D2,D3,D4,D5,D6,D7,D8,D9  !disorder /Berman (HP:D1,D2,D3)
    real(dp):: TD0,TDMax,VAdj              !disorder /Berman 
    !
    real(dp):: VTA,VTB,VPA,VPB             !volume function
    real(dp):: TKri,SMA                    !
    ! a la Berman: use VTA,VTB,VPA,VPB
    ! a la HP:     use VTA,VTB,VPA,TKRI,SMA
    logical :: Dis
  end type T_DtbMinThr
  !
contains

subroutine DtbMinThr_Zero(M)
  type(T_DtbMinThr),intent(out)::M
  !
  M%Div=    1
  M%codVol= 1
  M%nLanda= 0
  M%Dis=.false. !; M%VdW=.false.; M%RdK=.false.; M%Min=.false. !M%TL1=.false.
  M%Name=    "ZZZ"
  M%Special= "ZZZ"
  M%Typ=     "ZZZ"
  M%CodGas=  "IDEAL"
  M%AA0=Zero;   M%AAT=Zero;   M%BB0=Zero; M%BBT=Zero !for VdW/RdK EoS of Gases
  M%Tcrit=Zero; M%Pcrit=One;  M%Acentric=Zero
  M%G0R=Zero;   M%H0R=Zero;   M%S0_=Zero; M%V0R=Zero
  M%S0Ele=Zero; M%WeitKg=Zero
  M%K1= Zero;   M%K2= Zero;   M%K3= Zero; M%K4= Zero
  M%K5= Zero;   M%K6= Zero;   M%K6= Zero; M%K7= Zero; M%K8=Zero; M%K9=Zero
  M%D1= Zero;   M%D2= Zero;   M%D3= Zero; M%D4= Zero
  M%D5= Zero;   M%D6= Zero;   M%D6= Zero; M%D7= Zero; M%D8=Zero
  M%TD0=Zero;   M%TDMax=Zero; M%VAdj=Zero
  M%VTA=Zero;   M%VTB=Zero;   M%VPA=Zero; M%VPB=Zero
  M%TKri=Zero;  M%SMA=Zero;   M%VPA=Zero
  !
end subroutine DtbMinThr_Zero

integer function DtbMinThr_Index(Str,V) 
!--
!-- index of mineral named Str in vDtbMinThr,
!-- used in ReadSolution
!--
  character(*),     intent(in):: Str
  type(T_DtbMinThr),intent(in):: V(:)
  !
  integer     ::I
  !
  DtbMinThr_Index=0
  !
  I=0
  do
    I=I+1 !; if(iDebug>0) write(fTrc,'(A)') vEle(I)%SpName
    if(trim(Str)==trim(V(I)%Name)) then;
      DtbMinThr_Index=I
      exit
    end if
    if(I==size(V)) exit
  enddo
  !if Str not found -> DtbMinThr_Index=0
  !
end function DtbMinThr_Index

subroutine DtbMinThr_Calc(M,TdgK,Pbar,S)

  use M_Numeric_Const, only: Ln10
  use M_T_Species,only: T_Species
  use M_Dtb_Const,only: R_jk,Pref,Tref,DtbConv_Benson
  use M_Fluid_Calc,only: Eos_H2O_Haar_Ghiorso
  !
  type(T_DtbMinThr),intent(in) :: M
  real(dp),         intent(in) :: Pbar,TdgK
  type(T_Species),  intent(inout):: S
  !
  real(dp):: DH,DS,DG,DCp,DV
  real(dp):: HR,SR,GR,CpR,Volum,LnFug
  integer ::I
  !
  S%WeitKg= M%WeitKg
  !
  select case(trim(M%Special))
  
  case("H2O_HGK")
    call Eos_H2O_Haar_Ghiorso(TdgK,Pbar,S%G0rt,S%H0,S%S0,S%V0)
    !call CalcGH2O_Supcrt(TdgK,Pbar,S%G0rt,S%H0,S%S0,S%V0)
    return
  
  end select
  !
  HR=    M%H0R
  SR=    M%S0_
  Volum= M%V0R
  !
  GR=    M%H0R - TdgK*M%S0_
  if(DtbConv_Benson) GR= GR + Tref *M%S0Ele !!!Benson Convention!!!
  !
  !-- first heat the mineral or gas from Tref to T current
  !-- at constant P(-Pref)
  call Heating(M,Tref,TdgK,DH,DS,DCp)
  !
  HR=    HR + DH     !H(T,P0)=H(T0,P0) + Sum<T0->T><(@H@T)dT>
  SR=    SR + DS     !S(T,P0)=S(T0,P0) + Sum<T0->T><(@S@T)dT>
  CpR=   DCp
  DG=    DH - TdgK*DS
  GR=    GR + DG
  !
  !-- Now add volume term to G(T,P0) to obtain G(T,P)
  call Compressing(M,Tref,TdgK,Pref,Pbar,DH,DS,DG,Volum,LnFug)
  !
  HR= HR + DH
  SR= SR + DS
  GR= GR + DG
  !
  do I=1,M%NLANDA
    call Capitani_Landau(M,I,TdgK,Pbar,DH,DS,DG,DV,DCp)
    HR=    HR  + DH
    SR=    SR  + DS
    GR=    GR  + DG
    Volum= Volum + DV
    CpR=   CpR + DCp
  enddo
  !
  if (M%Dis.and.M%TD0/=Zero .and. M%TDMax/=Zero .and. TdgK>M%TD0) then
    call Berman_Disorder(M,TdgK,Pbar,DH,DS,DG,DV,DCp)
    GR=    GR    + DG
    HR=    HR    + DH
    SR=    SR    + DS
    CpR=   CpR   + DCp
  end if
  !
  !-- caveat programmator: here, Molar Volumes in J/bar
  !-- -> Volum(cm3)-Volum(J/bar)*10.0D0
  !--    Volum(m3)- Volum(J/bar)*1.0D-5
  !
  S%V0=  Volum *1.0D-5 !-> conversion from J/bar/Mole to M^3/Mole
  S%H0=  HR
  S%S0=  SR
  S%Cp0= CpR
  S%G0rt=GR /R_jk /TdgK 
  !
  if(S%Typ=="GAS") S%LnFug= LnFug
  !
  return
end subroutine DtbMinThr_Calc

subroutine Compressing( & !
& M, &            ! IN:  species
& T0,T, P0,P, &   ! IN:  T,P integration limits
& DH,DS,DG,VR,&   ! OUT: variations in H,S,G,V
& LnFug)          ! OUT
!--
!-- H(T,P)- H(T,P0) + Sum-P0->P>(@H@P)dP - H(T,P0) + Sum-P0_P>-(V - T(@V@T))dP>
!-- S(T,P)- S(T,P0) + Sum-P0->P>(@S@P)dP - S(T,P0) - Sum-P0_P>-(@V@T) dP>
!-- (@G@P)- V
!-- G(T,P)- G(T,P0) + Sum-P0->P>-(@G@P)dP> - G(T,P0) + Sum-P0_P>-VdP>
!--
  use M_Dtb_Const,only: R_jk
  type(T_DtbMinThr),intent(in):: M !mineral
  real(dp),intent(in) :: T0,T,P0,P
  real(dp),intent(out):: DH,DS,DG,VR
  real(dp),intent(out):: LnFug
  !
  real(dp):: Integr_VdP,Integr_dVdTdP
  !
  real(dp):: DGGas
  real(dp):: PK,DPK,A00,A22,KAPS,DKDT,KTT,VMA,FV,TKR,TEF
  real(dp):: Q2,Q4,Q6,Q20,Q40,Q60
  real(dp):: FVLA,FVDP,FHLA,FSLA,DGLA
  !
  DH= Zero ; DS= Zero ; DG= Zero ; VR= Zero
  LnFug= Zero
  !
  !------------------------------------------------------------------GAS
  if (M%Typ=="GAS") then
    ! adapted from Theriak, deCapitani
    VR= R_jk*T /P *1.D-5 !P(bar)-> RT/P(J/bar) (1.M^3 =1.D+5.J/bar)
    DGGas= Zero
    !
    if(trim(M%CodGas)/="IDEAL") then
      call DtbGasThr_Calc( &
      & M,P0,P,T, &     ! IN
      & DGGas,VR)       ! OUT
      ! volume conversion cm^3 to m^3
      VR= VR*1.D-6
    end if
    !
    DG= R_jk*T*log(P/P0) + DGGas
    LnFug= log(P/P0) + DGGas /T /R_jk
    !!
    !! should also compute DH and DS !!!
    !!
    return !------------------------------------------------------return
  end if
  !-----------------------------------------------------------------/GAS
  
  !------------------------------------------------------------------MIN
  select case(M%codVol)
  
  case(1)
    ! mineral volume according to Berman volume function
    ! adapted from Theriak, deCapitani
    VR= M%V0R &
    & + M%VTA*(T-T0) &
    & + M%VTB*(T-T0)*(T-T0) &
    & + M%VPA*(P-P0) &
    & + M%VPB*(P-P0)*(P-P0)
    !
    Integr_VdP= (P-P0)*(M%V0R + M%VTA*(T-T0) + M%VTB*(T-T0)**2) &
    &         +  M%VPA*(P*P /2.0D0 -P  *P0           +P0*P0/2.0D0) &
    &         +  M%VPB*(P**3/3.0D0 -P*P*P0 +P*P0*P0  -P0**3/3.0D0)
    !
    Integr_dVdTdP= (M%VTA + 2*M%VTB*(T-T0)) *(P-P0)
    !
    DH= Integr_VdP + T*Integr_dVdTdP 
    DS=            -   Integr_dVdTdP
    DG= Integr_VdP
    !
  !case(2) then
  !  fVVol= M%V0R*exp(M%VAA*(T-T0)+M%VAB*(TT-TT0)/2.0D0)
  !  if (M%VB>Zero) then
  !    Integr_VdP= (fVVol/M%VB)*(One-exp(-M%VB*(P-P0)))
  !    Volum= fVVol*exp(-M%VB*(P-P0))
  !    fVVol= Volum-M%V0R
  !  else
  !    Integr_VdP= fVVol*(P-P0)
  !    Volum= fVVol
  !    fVVol= Volum-M%V0R
  !  end if
  !  GR= GR + Integr_VdP
  !end if
  !if (M%codVol==3) then
  !  fVVol=M%VL0+M%VLA*T
  !  if (M%VL2==Zero.and.M%VLN/=Zero) then
  !    FF= 1.0D-6/(0.7551D0+2.76D0*FV/M%VLN)
  !    Integr_VdP= fVVol &
  !              * ((P-P0) - FF*((P**2-P0**2)/2.0D0 - FF*(P**3-P0**3)/3.0D0))
  !    Volum= fVVol*(One-FF*(P-FF*P**2/2.0D0))
  !  else
  !    Integr_VdP= fVVol*(P-P0) &
  !              + M%VL2*(P**2-P0**2)
  !    Volum= FV+2.0D0*M%VL2*P
  !  end if
  !  fVVol=Volum-M%V0R
  !  GR= GR+ Integr_VdP
  !end if
  !
  case(4)
  if(M%VTB/=Zero) then
    ! adapted from Theriak, deCapitani
    ! integration of VdP from P=0.001 kBar to P kBar
    !
    ! currently no computation of DH nor DS
    !
    PK=   P*1.0D-3 !P in kBar
    DPK=  PK-1D-3  !deltaP en kBar= PKbar - 0.001
    A00=  M%VTA
    A22=  M%D1 !in HoPo1998 A22  is always = 10
    KAPS= M%D2 !in HoPo1998 KAPS is always = 4
    DKDT= M%D3 !in HoPo1998 dkdt is always =-1.5D-4
    !
    KTT= M%VTB*(One+DKDT*(T-T0))
    VMA= M%VPA
    FV=  M%V0R &
    &    *(One + A00*(T-T0) - 20.0D0*A00*(SQRT(T)-SQRT(T0)))
    !VR=FV*((One-4.0D0*PK/(KTT+4.0D0*PK))**0.25D0)
    VR= FV *(KTT/(4.0D0*DPK+KTT))**0.25D0
    !
    DG= 1.0D3 *FV *KTT /3.0D0 *( (One + 4.0D0*DPK/KTT)**0.75D0 -One )
    !
    if (M%TKri/=Zero.and.M%SMA/=Zero) then
      TKR= M%TKri + VMA /M%SMA *P
      if (T<TKR) then  ; TEF=T
      else             ; TEF=TKR
      end if
      Q40=  One -T0/M%TKri
      Q4=   One -TEF/TKR
      Q20=  SQRT(Q40)
      Q2=   SQRT(Q4)
      Q60=  Q20*Q40
      Q6=   Q2*Q4
      !
      DGLA= M%SMA*((TEF-TKR)*Q2 + TKR*Q6/3.0D0)
      !
      FVLA= VMA *Q20 &
      &   *( One + A00*(TEF-T0) - 20.0D0*A00*(SQRT(TEF)-SQRT(T0)) )
      FVDP= FVLA*KTT/3.0D0 *((One+4.0D0*DPK/KTT)**0.75D0 -One)
      FVDP= FVDP*1D3
      FHLA= M%SMA *M%TKri *(Q20 - Q60/3.0D0)
      FSLA= M%SMA *Q20
      !
      DG= DG + FHLA -T*FSLA +FVDP +DGLA
      !
      VR= VR + FVLA
      !
    end if !(TKri/=Zero.and.SMA/=Zero)
  end if !if (VOLVO==4.and.VTB/=Zero)
  
  end select
  !-----------------------------------------------------------------/MIN
  !
end subroutine Compressing

subroutine Heating(M,T0,T,DH,DS,DCp)
  type(T_DtbMinThr),intent(in):: M !mineral
  real(dp),intent(in) :: T0,T
  real(dp),intent(out):: DH,DS,DCp
  !
  real(dp):: Cp0, TT, TT0, SQT, SQT0
  real(dp):: Cp,CpdT,CpTdT
  !
  !dH= (@H@T)dT + (@H@P)dP =   (Cp)dT + (V - T(@V@T))dP
  !dS= (@S@T)dT + (@S@P)dP = (Cp/T)dT - (@V@T)dP 
  !dG= 
  !
  !HKF (Maier-Kelley): Cp= a + b*T + c/T**2; a=K1; b=K2; c=K3
  !K1= Cp0   Cp_ref
  !K2= CpT   + _ *T
  !K3= CpT_2 + _ /T/T
  !K4= CpT_r + _ /SQRT(T)
  !K5= CpT2  + _ *T*T
  !K6= CpT_1 + _ /T
  !K7= CpTr  + _ *SQRT(T)
  !K8= CpT_3 + _ /T/T/T
  !K9= CpT3  + _ *T*T*T
  !
  TT=  T*T       ;   TT0=  T0*T0
  SQT= SQRT(T)   ;   SQT0= SQRT(T0)
  !
  !--- Cp at temp. T-T0
  Cp0=   M%K1       &
       + M%K2 *T0   &
       + M%K3 /TT0  &
       + M%K4 /SQT0 &
       + M%K5 *TT0  &
       + M%K6 /T0   &
       + M%K7 *SQT0 &
       + M%K8 /TT0/T0 &
       + M%K9 *TT0*T0
  !--- Cp at temp. T
  Cp=    M%K1       & 
       + M%K2 *T    &
       + M%K3 /TT   &
       + M%K4 /SQT  &
       + M%K5 *TT   &
       + M%K6 /T    &
       + M%K7 *SQT  &
       + M%K8 /TT/T &
       + M%K9 *TT*T
  !--- sum_{from T0 to T} (@H@T)dT - sum_{from T0 to T}(CpdT)- CpRdT
  CpdT=  M%K1 *(T       -T0     ) & !_______
       + M%K2 *(TT      -TT0    )/2.0D0 &
       - M%K3 *(One/T   -One/T0 ) &
       + M%K4 *(SQT     -SQT0   )*2.D0 &
       + M%K5 *(TT*T    -TT0*T0 )/3.0D0 &
       + M%K6 *log(T/T0) &
       + M%K7 *(T*SQT   -T0*SQT0)*2.0D0/3.0D0 &
       - M%K8 *(One/TT  -One/TT0)/2.0D0 &
       + M%K9 *(TT*TT   -TT0*TT0)/4.0D0
  !--- sum_{from T0 to T} (@S@T)dT - sum_{from T0 to T} -(Cp/T)dT>- CpRTdT
  CpTdT= M%K1 *log(T/T0)             & !____
       + M%K2 *(T       -T0        ) &
       - M%K3 *(One/TT  -One/TT0   )/2.0D0 &
       - M%K4 *(One/SQT -One/SQT0  )*2D0 &
       + M%K5 *(TT      -TT0       )/2.0D0 &
       - M%K6 *(One/T   -One/T0    ) &
       + M%K7 *(SQT     -SQT0      )*2D0 &
       - M%K8 *(One/TT/T-One/TT0/T0)/3.0D0 &
       + M%K9 *(TT*T    -TT0*T0    )/3.0D0
  DH=  CpdT
  DS=  CpTdT
  DCp= Cp
  !
end subroutine Heating

subroutine Capitani_Landau(M,K,TdgK,Pbar,FGTR,FHTR,FSTR,FVTR,FCPTR)
  !adapted from Theriak, deCapitani
  type(T_DtbMinThr),intent(in):: M !mineral
  integer,       intent(in)   :: K !which transition
  real(dp),      intent(in)   :: Pbar,TdgK
  real(dp),      intent(out)  :: FGTR,FHTR,FSTR,FVTR,FCPTR
  !parameters in tMin
  !TQ1B(1:3),TRE(1:3), ASPK(1:3),BSPK(1:3),DHTR(1:3),&
  !TEQ(1:3), DVTR(1:3),DVDT(1:3),DVDP(1:3)
  real(dp):: Pr,Tr
  real(dp):: A11,B11,C11,L22,CTR,T1,T2
  real(dp):: DHSPK,DSSPK,TQP,DSTR
  !~ real(dp) :: A11,B11,C11,D11,CTR,T9,TR9,DHSPK,DSSPK,GSPK, AASP,ABSP,CCTR,TEQK
  !
  Pr=    One !bar
  Tr=    298.15D0 !dgK
  !
  !~ TEQK=  M%TEQ(K)*(Pbar-Pr) +M%TQ1B(K)
  !~ CTR=   M%TQ1B(K) -TEQK
  !~ TR9=   M%TRE(K)  -CTR
  !~ if(TdgK>TEQK) then; T9=TEQK
  !~ else              ; T9=TdgK
  !~ end if
  !~ D11=   M%BSPK(K)*M%BSPK(K)
  !~ AASP=  M%ASPK(K)*M%ASPK(K)
  !~ ABSP=  M%ASPK(K)*M%BSPK(K)
  !~ CCTR=  CTR*CTR
  !~ A11=   AASP*CTR  +ABSP*CCTR*2.0D0 +D11*CCTR*CTR
  !~ B11=   AASP      +ABSP*CTR *4.0D0 +D11*CCTR*3.0D0
  !~ C11=              ABSP*2.0D0      +D11*CTR *3.0D0
  !~ DHSPK= A11*(T9     -TR9    )       &
       !~ + B11*(T9*T9  -TR9*TR9)/2.0D0 &
       !~ + C11*(T9**3  -TR9**3 )/3.0D0 &
       !~ + D11*(T9**4  -TR9**4 )/4.0D0
  !~ DSSPK= A11*(log(T9)-log(TR9))       &
       !~ + B11*(T9     -TR9     )       &
       !~ + C11*(T9*T9  -TR9*TR9 )/2.0D0 &
       !~ + D11*(T9**3  -TR9**3  )/3.0D0
  !~ if ((T9+CTR) < M%TRE(K)) then
    !~ DHSPK= Zero
    !~ DSSPK= Zero
  !~ end if
  !~ !
  !~ GSPK= DHSPK -T9*DSSPK
  !~ if (TdgK>TEQK) GSPK=GSPK-(M%DHTR(K)/M%TQ1B(K)+DSSPK)*(TdgK-TEQK)
  !~ GSPK= GSPK &
  !~ &   + M%DVDT(K)*(Pbar      -Pr   )*(T9-Tr) &
  !~ &   + M%DVDP(K)*(Pbar*Pbar -Pr*Pr)/2.0D0  &
  !~ &   - M%DVDP(K)*(Pbar-Pr)
  !~ !
  !~ GR= GR+GSPK
  L22= M%BSPK(K)**2
  CTR= M%TEQ(K)*(Pbar-1.0D0)
  TQP= M%TQ1B(K) + CTR
  T2=  MIN(TdgK,TQP)
  !
  if (M%BSPK(K)/=Zero) then  ;  T1= -M%ASPK(K)/M%BSPK(K) + CTR
  else                       ;  T1= T2
  end if
  !
  if (TdgK<T1) then
    !
    DHSPK= Zero ; DSSPK= Zero
    FCPTR= Zero ; FHTR=  Zero ; FSTR=  Zero ; FGTR=  Zero
    !
  else
    !
    A11= -CTR*T1**2
    B11=  T1*(2.0D0*CTR+T1)
    C11= -(CTR+2.0D0*T1)
    DHSPK= A11*(T2-T1) &
    &    + B11*(T2**2-T1**2)/2.0D0 &
    &    + C11*(T2**3-T1**3)/3.0D0 &
    &    +     (T2**4-T1**4)/4.0D0
    DSSPK= A11*(log(T2/T1)) &
    &    + B11*(T2-T1) &
    &    + C11*(T2**2-T1**2)/2.0D0 &
    &    +     (T2**3-T1**3)/3.0D0
    FHTR=  DHSPK *L22
    FSTR=  DSSPK *L22
    FGTR= (DHSPK-TdgK*DSSPK)*L22
    !
    if (TdgK < TQP) then
      FCPTR= (TdgK-CTR) *(TdgK-T1)**2 *L22
    else
      FCPTR= Zero
      DSTR=  M%DHTR(K) /M%TQ1B(K)
      FHTR=  FHTR +DSTR*TQP
      FSTR=  FSTR +DSTR
      FGTR=  FGTR +DSTR*TQP -TdgK*DSTR
    end if
    !
  end if
  !
  FGTR= FGTR &
  &   + M%DVDT(K) *(Pbar-1.0D0)  *(T2-298.15D0) &
  &   + M%DVDP(K) *(Pbar*Pbar-1.0D0)/2.0D0 &
  &   - M%DVDP(K) *(Pbar-1.0D0)
  !
  FVTR= M%DVDT(K) *(T2-298.15D0) &
  &   + M%DVDP(K) *(Pbar-1.0D0)
  !
  !~ HR=   HR    +FHTR
  !~ SR=   SR    +FSTR
  !~ CPR=  CPR   +FCPTR
  !~ VOLUM=VOLUM +FVTR
  !~ GR=   GR    +GSPK
  return
end subroutine Capitani_Landau

subroutine  Berman_Disorder(M,TdgK,Pbar,Hdis,Sdis,Gdis,Vdis,Cpdis)
!--
!-- for data with Berman format,
!-- Add disorder contributions to Cp, H, S, G, V
!-- adapted from Theriak, deCapitani
!--
  type(T_DtbMinThr),intent(in) :: M !mineral
  real(dp),         intent(in) :: TdgK,Pbar
  real(dp),         intent(out):: Hdis,Sdis,Gdis,Vdis,Cpdis
  real(dp):: TD, CpRdT, CpRTdT
  TD= MIN(TdgK,M%TDMax) !DMIN1(T,M%TDMax) DMIN1: arg=Real, res=Real
  
  CpDis=  M%D1 &
  &     + M%D2 *TD &
  &     + M%D3 /(TD*TD) &
  &     + M%D4 /SQRT(TD) &
  &     + M%D5 *TD*TD &
  &     + M%D6 /TD &
  &     + M%D7 *SQRT(TD) &
  &     + M%D8 /(TD**3) &
  &     + M%D9 *(TD**3)
  
  CpRdT=  M%D1 *(TD          -M%TD0      ) &
  &     + M%D2 *(TD*TD       -M%TD0*M%TD0)/2.0D0 &
  &     - M%D3 *(One/TD      -One/M%TD0  ) &
  &     + M%D4 *(SQRT(TD)    -SQRT(M%TD0))*2D0 &
  &     + M%D5 *(TD**3       -M%TD0**3   )/3.0D0 &
  &     + M%D6 *log(TD/M%TD0) &
  &     + M%D7 *(TD*SQRT(TD) -M%TD0*SQRT(M%TD0))*2.0D0/3.0D0 &
  &     - M%D8 *(One/(TD*TD) -One/(M%TD0*M%TD0))/2.0D0 &
  &     + M%D9 *(TD**4       -M%TD0**4         )/4.0D0
  
  CpRTdT= M%D1 *log(TD/M%TD0) &
  &     + M%D2 *(TD          -M%TD0            ) &
  &     - M%D3 *(One/(TD*TD) -One/(M%TD0*M%TD0))/2.0D0 &
  &     - M%D4 *(One/SQRT(TD)-One/SQRT(M%TD0)  )*2D0 &
  &     + M%D5 *(TD*TD       -M%TD0*M%TD0      )/2.0D0 &
  &     - M%D6 *(One/TD      -One/M%TD0        ) &
  &     + M%D7 *(SQRT(TD)    -SQRT(M%TD0)      )*2D0 &
  &     - M%D8 *(One/(TD**3) -One/(M%TD0**3)   )/3.0D0 &
  &     + M%D9 *(TD**3       -M%TD0**3         )/3.0D0
  !
  if (ABS(M%VAdj)>10.0D0) then  ;  VDis= CpRdT/(10.0D0*M%VAdj)
  else                          ;  VDis= Zero
  end if
  !
  HDis= CpRdT
  SDis= CpRTdT
  GDis= CpRdT -TdgK*CpRTdT +VDis*(Pbar-One)
  !
end subroutine Berman_Disorder

subroutine CUBIC(A2,A1,A0,X1,X2,X2I,X3)
!--
!-- find roots of X^3 +A2.X^2 +A1.X +A0 - 0
!--
  use M_Numeric_Const,only:Pi
  real(dp),intent(in) :: A2,A1,A0
  real(dp),intent(out):: X1,X2,X2I,X3
  !
  real(dp)::Q,P,R,PHI3,FF
  !
  X2=  Zero
  X2I= Zero
  X3=  Zero
  !
  if (A1==Zero .and. A0==Zero) then !X^2.(X+A2) = 0 
    X1=-A2
    return
  end if
  !
  Q=(2.0D0*A2*A2*A2/27.D0 - A2*A1/3.D0 +A0) /2.D0
  P=(3.0D0*A1             - A2*A2         ) /9.D0
  FF=ABS(P)
  R=SQRT(FF)
  FF=R*Q
  if (FF < Zero) R=-R
  FF=Q/R/R/R
  
  if (P > Zero) then
    PHI3= log(FF+SQRT(FF*FF+ One))    /3.D0
    X1=  -R*(exp(PHI3)-exp(-PHI3)) -A2/3.D0
    X2I=  One
  else
    if (Q*Q+P*P*P > Zero) then
      PHI3= log(FF+SQRT(FF*FF-1.D0))/3.D0
      X1=  -R*(exp(PHI3)+exp(-PHI3)) -A2/3.D0
      X2I=  One
    else
      PHI3= ATAN(SQRT(1.D0-FF*FF)/FF)/3.D0
      X1=  -2.0D0*R*COS(PHI3)         -A2/3.D0
      X2=   2.0D0*R*COS(PI/3.D0-PHI3) -A2/3.D0
      X2I=  Zero
      X3=   2.0D0*R*COS(PI/3.D0+PHI3) -A2/3.D0
    end if
  end if
  !
  return
end subroutine CUBIC

subroutine DtbGasThr_Calc( &
& M, &
& Pref,Pbar,TdgK, &
& DGGas,Vol)
!--
!-- adapted from Theriak / CdeCapitani
!-- compute deltaG of isothermal (at TdgK) compression from Pref to Pbar
!-- deltaG - integr VdP
!--
  use M_Dtb_Const, only:R_jk
  !
  type(T_DtbMinThr),intent(in) :: M
  real(dp),         intent(in) :: Pref,Pbar,TdgK
  real(dp),         intent(out):: DGGas,Vol
  !
  real(dp):: A2,A1,A0
  real(dp):: VRef,VGas,VLiq
  real(dp):: AA,BB,X1,X2,X2I,X3
  !
  call CubicEoS_CalcParamsAB(M,TdgK,AA,BB)
  !
  !--------- FIRST compute VRef at PRef using the Cubic EoS in volume --
  call CubicEqu_CalcCoeffs(M,TdgK,Pref,AA,BB,A2,A1,A0)
  call CUBIC(A2,A1,A0,X1,X2,X2I,X3)
  if (X2I/=Zero) then
    VRef=X1            ! if one root, one-phase
  else
    VRef=MAX(X1,X2,X3) ! if 3 real roots, VGas is highest root
  end if
  !
  !-------- then compute Volume at Pbar using the Cubic EoS in volume --
  call CubicEqu_CalcCoeffs(M,TdgK,Pbar,AA,BB,A2,A1,A0)
  call CUBIC(A2,A1,A0,X1,X2,X2I,X3)
  if (X2I/=Zero) then !if one root, one-phase
    Vol=X1
  else ! if 3 real roots, VGas is highest root, VLiq the lowest
    VGas= MAX(X1,X2,X3)
    VLiq= MIN(X1,X2,X3)
    DGGas= -One
    !
    if (VLiq>BB) &
    & DGGas= CubicEoS_RTLnPhi(M,AA,BB,TdgK,Pbar,VGas) &
    &      - CubicEoS_RTLnPhi(M,AA,BB,TdgK,Pbar,VLiq)
    !
    ! stable phase is phase with lowest Phi
    if (DGGas >Zero) then; Vol=VLiq
    else                 ; Vol=VGas
    end if
  end if
  !
  DGGas= CubicEoS_RTLnPhi(M,AA,BB,TdgK,Pbar,Vol)
  !
  return
end subroutine DtbGasThr_Calc

real(dp) function CubicEoS_RTLnPhi(M,AA,BB,T,P,V)
  use M_Dtb_Const,only: R_jk
  
  type(T_DtbMinThr),intent(in):: M
  real(dp),         intent(in):: AA,BB
  real(dp),         intent(in):: T,P,V
  !
  real(dp):: XX,RTbar
  real(dp):: SQT2
  !
  ! compute R*T with gas constant in bar*cm3 /mol/K
  RTbar= R_jk *T *10.0D0
  !
  XX= Zero
  !
  select case(trim(M%CodGas))
  
  case("REDKWON") ! Redlich-Kwong
    XX= V*P - RTbar &
    & - RTbar *log(P/RTbar) &
    & - RTbar *log(V-BB) &
    & - AA/BB *log((V+BB)/V) /SQRT(T)
  
  case("VDWAALS") ! Van der Waals
    XX= V*P &
    & - RTbar *log(V-BB) &
    & - AA/V
  
  case("SOAVE") ! Soave
    XX= V*P - RTbar - RTbar *log(P/RTbar) &
    & - RTbar *log(V-BB) &
    & - AA/BB *log((V+BB)/V)
  
  case("PENGROB","PRSV") ! Peng-Robinson
    SQT2=  SQRT(Two)
    XX= V*P - RTbar - RTbar *log(P/RTbar) &
    & - RTbar *log(V-BB) &
    & - AA/BB *log((V+BB*(One+Sqt2))/(V+BB*(One-Sqt2))) /Sqt2/Two
  
  end select
  !
  ! conversion from from Cm3.Bar to Joule
  CubicEoS_RTLnPhi= XX *1.D-1 
  !
  return
end function CubicEoS_RTLnPhi

subroutine CubicEqu_CalcCoeffs( &
& M,     &  !IN
& T,P,   &  !IN
& AA,BB, &  !IN
& A2,A1,A0) !OUT
!--
!-- coeff's for cubic equation X^3 +A2.X^2 +A1.X +A0 - 0
!--
  use M_Dtb_Const,only: R_jk
  type(T_DtbMinThr),intent(in) :: M
  real(dp),         intent(in) :: T,P
  real(dp),         intent(in) :: AA,BB
  real(dp),         intent(out):: A2,A1,A0
  !
  real(dp):: RTbar,A3
  !
  RTbar= R_jk *T *10.0D0
  A3= P
  A2=-RTbar/A3
  A1= Zero
  A0= Zero
  select case(trim(M%CodGas))
  
  case("REDKWON")
  !Coeffs in the cubic RdK equation in volume
  !P0.SQT.V3 - R.TdgK.SQT.V2 + (A3 - R.TdgK.SQT.A2 -A2.A2.P.SQT).V - A3.A2=0
  !V3 - (R.TdgK/P).V2 + (A3/(P0.SQT) - A2.R.TdgK/P - A2.A2).V - A3.A2/(P0.SQT)= 0
    A3=  P *SQRT(T) 
    A2= -RTbar /P
    A1= -BB*BB +AA/A3 +BB*A2
    A0= -AA*BB /A3
  
  case("VDWAALS")
    A3=  P
    A2= -BB -RTbar/A3
    A1=  AA /A3
    A0= -AA*BB /A3
  
  case("SOAVE")
    A3=  P
    A2= -RTbar /A3
    A1= -BB*BB -BB*RTbar/A3 +AA/A3 
    A0= -AA*BB /A3
  
  case("PENGROB") !Peng-Robinson,1976
    A3=  P
    A2= -RTbar/A3 +BB
    A1= -BB*BB*3.D0 -2.D0*BB*RTbar/A3 +AA/A3
    A0= (BB*BB           +BB*RTbar/A3 -AA/A3) *BB
  
  end select
  !
  A3= One
  !
  return
end subroutine CubicEqu_CalcCoeffs

subroutine CubicEoS_CalcParamsAB(M,TdgK,AA,BB)
  use M_Dtb_Const,only: R_jk
  type(T_DtbMinThr),intent(in) :: M
  real(dp),intent(in) :: TdgK
  real(dp),intent(out):: AA,BB
  !
  real(dp):: Rbar,RTbar,MM,Alpha
  real(dp):: K1 !should be a species-dependent parameter ..!!
  !
  Rbar=  R_jk *10.0D0
  RTbar= R_jk *TdgK *10.0D0 !R*T with gas constant in bar*cm3 /mol/K
  !
  select case(trim(M%CodGas))
  
  case("REDKWON","VDWAALS")
    AA= M%AA0 +TdgK*M%AAT
    BB= M%BB0 +TdgK*M%BBT
  
  case("PENGROB") !Peng-Robinson,1976
    !AA= OmegaA * (R*Tcrit)**2 /Pcrit ->
    AA= 0.45724D0 *Rbar*Rbar *M%Tcrit*M%Tcrit /M%Pcrit
    MM= 0.37464D0 &
    & + 1.54226D0 *M%Acentric &
    & - 0.26992D0 *M%Acentric*M%Acentric
    alpha= One + MM *(One -SQRT(TdgK /M%Tcrit))
    alpha= alpha*alpha
    AA= alpha*AA
    !BB= OmegaB * R*Tcrit/Pcrit ->
    BB= 0.07780D0 *Rbar *M%Tcrit /M%Pcrit
  
  case("PRSV") !Peng-Robinson-Stryjek-Vera
    K1= -0.06635D0 !provisonal !! (value for H2O, should be read from file) !!
    AA= 0.457235D0 *Rbar*Rbar *M%Tcrit*M%Tcrit /M%Pcrit
    !AA= OmegaA * (R*Tcrit)**2 /Pcrit
    MM= 0.378893D0 &
    & + 1.4897153D0  *M%Acentric    &
    & - 0.17131848D0 *M%Acentric**2 &
    & + 0.0196554D0  *M%Acentric**3
    MM= MM + K1*(One +SQRT(TdgK /M%Tcrit))*(0.7D0 -TdgK /M%Tcrit)
    alpha= One + MM *(One -SQRT(TdgK /M%Tcrit))
    alpha= alpha*alpha
    AA= alpha*AA
    BB= 0.077796D0 *Rbar *M%Tcrit /M%Pcrit
    !BB= OmegaB * R*Tcrit/Pcrit
  
  case("SOAVE") !Soave,1972
    AA= 0.42748D0 *Rbar*Rbar *M%Tcrit*M%Tcrit /M%Pcrit
    !1 /(9*2**(1/3) - 9) -> 0.42748
    MM= 0.480D0 &
    & + 1.574D0 *M%Acentric &
    & - 0.176D0 *M%Acentric*M%Acentric
    ! MM= 0.48508D0 &
    ! & + 1.55171D0*M%Acentric &
    ! & - 0.15613D0*M%Acentric*M%Acentric
    alpha= One + MM *(One - SQRT(TdgK /M%Tcrit))
    alpha= alpha*alpha
    AA= alpha*AA
    !
    !(2**(1/3) - 1)/3 -> 0.086640349965D0
    BB= & !Bpure_Soave
    & 0.08664D0 *Rbar *M%Tcrit /M%Pcrit
  
  end select
  !
end subroutine CubicEoS_CalcParamsAB

end module M_T_DtbMinThr

!TWQ/Theriak/databases (Berman/Capitani)
!code= contains
!---------- ----------------
!ST= G0R,H0R,S0_,V0R
!R-K= AA0,AAT,BB0,BBT
!VDW= AA0,AAT,BB0,BBT
!CP1 or C1= K1,K4,K3,K8
!CP3 or C3= K1,K2,K3,K4
!CP2 or C2= K6,K2,K5,K7,K9
!D1= D1,D4,D3,D8,D6
!D2= D2,D5,TD0,TDMax,VAdj
!AQ1= K1,K2,K8,K9,K7
!AQ2= D1,D2,D3,D4
!AQP= D1,D2
!SPC
!COM
!TR1 or T1= TQ1B(I),TRE(I),ASPK(I),BSPK(I),DHTR(I), I=1,NLANDA
!TR2 or T2= TEQ(I),DVTR(I),DVDT(I),DVDP(I)
!TL1= TKri,SMA
!V1= VTA,VTB,VPA,VPB
!  VTA=VTA/100000.0D0*V0R
!  VTB=VTB/100000.0D0*V0R
!  VPA=VPA/100000.0D0*V0R
!  VPB=VPB/100000000.0D0*V0R
!  VO1=.true.
!  VO2=.false.
!  VO3=.false.
!V2= VAA,VAB,VB
!  VO2=.true.
!  VO1=.false.
!  VO3=.false.
!V3= VL0,VLA,VLN,VL2
!  VO3=.true.
!  VO1=.false.
!  VO2=.false.
!VHP= VTA,VTB,TKri,SMA,VPA
!  VOLVO=4
!  VO3=.false.
!  VO1=.false.
!  VO2=.false.
!VH2= D1,D2,D3
!  VOLVO=4
!  VO3=.false.
!  VO1=.false.
!  VO2=.false.

