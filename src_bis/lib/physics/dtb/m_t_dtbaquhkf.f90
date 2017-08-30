module M_T_DtbAquHkf
!--
!-- data type T_DtbAquHkf, to store parameters and results for aqu.species,
!-- and related routines for calculation of G,S,H,V at a given (P,T)
!-- program should read ref.state values from file to G0R,H0R,S0_, ...
!-- and compute (CalcAquSpXxx) G,S,H,V, ... for a given value of the H2O property (pW) 
!-- two 'versions' of CalcAquSpXxx are implemented:
!-- -- DtbAquHkf_Calc: apparent G according to Benson-Helgeson convention 
!-- -- DtbAquHkf_CalcThr: apparent G according to Berman-Brown convention
!-- the results differ by one term corresponding to Tref.Sum(Nu_i.S(Pref,Tref,element_i))
!--
  use M_Kinds
  use M_Trace
  !
  implicit none
  !
  private
  !
  public:: T_DtbAquHkf
  public:: DtbAquHkf_Calc
  public:: DtbAquHkf_CalcThr
  !
  !--- data structure for Aqu.Species following SupCrt, HKF EoS
  type:: T_DtbAquHkf
    character(len=15):: num
    character(len=23):: name
    character(len=71):: Formula
    integer :: Div=1
    integer :: Chg
    real(dp):: &
      G0R,H0R,S0_,V0R,&  !ref.state values
      WRef, &            !omega (Born)
      A1,A2,A3,A4,&      !EoS volume parameters (non-solvation part)
      C1,C2, &           !EoS Cp parameters (non-solvation part)
      S0Ele, WeitKg
  end type T_DtbAquHkf
  !---/

contains 
 
!~ subroutine DtbAquHkf_Init(M)
  !~ type(T_DtbAquHkf),intent(out)::M
  !~ M%Div=1
  !~ M%G0R=Zero; M%H0R= Zero; M%S0_=  Zero; M%V0R=Zero
  !~ M%C1= Zero; M%C2=  Zero
  !~ M%Chg=0;    M%Wref=Zero
  !~ M%A1= Zero; M%A2=  Zero; M%A3=   Zero; M%A4= Zero
  !~ M%S0Ele= Zero
  !~ !M%CpS=Zero; M%Gs=  Zero; M%Hs=   Zero; M%Ss= Zero; M%Vs= Zero
  !~ !M%CpR=Zero; M%GR=  Zero; M%HR=   Zero; M%SR= Zero; M%VR= Zero
  !~ !M%VR1=Zero; M%VR2=Zero
!~ end subroutine DtbAquHkf_Init

!Hkf data vs Thr Data
!hkf	Gf	Hf	Sref	a1	a2	a3	a4	c1	c2	omega	charge
!thr	G	H	M	A1	A2	A3	A4	K1	K2	K9	K8
!hkf	-115609	-126834	-77.7	-3.3802	-17.007	14.518	-2.0758	10.7	-8.06	2.753	3
!thr	-116543	-128681	-80.8	-0.334	-1711.0	14.992	-20716	10.7	-80600	3.33	3
!

!A1H=A1C/1E-1; A2H=A2C/1E2; A3H=A3C; A4H=A4C/1E4; C1H=C1C; C2H=C2C/1E5
!consistent with instructions in DtbIo:
!  S_%A1= S_%A1*1.0D-1
!  S_%A2= S_%A2*1.0D02
!  S_%A4= S_%A4*1.0D04
!  S_%C2= S_%C2*1.0D04
!  S_%wref= S_%wref*1.0D05

subroutine DtbAquHkf_Calc( &
& M,  & !IN:  database model
& pW, & !IN:  T,P and dielectric properties of Solvent
& S)    !OUT: species
!--
!-- computes the standard partial molal thermodynamic properties (Benson Convention) 
!-- of aqueous species s_ at T,P
!-- using Tanger And Helgeson (1988), Shock Et Al. (1991), and Johnson Et Al. (1991).
!--
  use M_Numeric_Const,only: Ln10
  use M_T_Element,    only: T_Element
  use M_T_Species,    only: T_Species
  use M_Dtb_Const,    only: Tref, Pref, CalToJoule, R_jk, T_CK, DtbConv_Benson
  
  !-- parameters for HKF EoS --
  use M_T_DtbH2OHkf, only: T_DtbH2OHkf
  use M_T_DtbH2OHkf, only: Theta, Psi
  use M_T_DtbH2OHkf, only: gShockRef, Eta, ZPrTr, YPrTr
  
  type(T_DtbAquHkf),intent(in) :: M
  type(T_DtbH2OHkf),intent(in) :: pW
  !
  type(T_Species),  intent(out):: S
  !
  integer :: Zj,Zj2
  real(dp):: P,T
  real(dp):: Q_,X_,Y_,Z_
  real(dp):: ReH,ReHRef ! values for H+
  real(dp):: Rej,RejRef 
  real(dp):: w,dwdP,dwdT,d2wdT2
  real(dp):: X1,X2, PpPsi, P0pPsi, TmTht, T0mTht
  !
  ! results for a given T,P, 
  real(dp)::  Vs,Ss,Cps,Hs,Gs  !solvation term
  real(dp)::  GR,HR,SR,VR,CpR  !total
  !
  real(dp):: Rej_Cation, Rej_Anion
  integer :: Nu_Cation,  Nu_Anion
  
  RejRef= Zero
  Rej_Cation=  1.91D0 != value for Na+ at Tref,Pref
  Rej_Anion=   1.81D0 != value for Cl- at Tref,Pref
  Nu_Cation=   1 ! NaCl stoichiometry
  Nu_Anion=    1 ! NaCl stoichiometry
  !
  P=     pW%Pbar
  T=     pW%TdgK
  PpPsi= P+Psi    ;  P0pPsi= Psi+Pref
  TmTht= T-Theta  ;  T0mTht= Tref-Theta
  !
  !--- Q_,X_,Y_,Z_: Born coefficents
  Q_= pW%Q !Q_= @(1/Eps)/@P -> = pW%dEpsdP /pW%Eps /pW%Eps
  Y_= pW%Y !Y_= @(1/Eps)/@T -> = pW%dEpsdT /pW%Eps /pW%Eps
  X_= pW%X !X_= @Y/@T       -> = (pW%d2EpsdT2 -2.0D0/pW%Eps*pW%dEpsdT**2)/pW%Eps/pW%Eps
  !        !N_= @Q/@P  (not used)
  Z_= pW%Z !Z_= One /pW%Eps - 1 !not the original definition of Z
  !---/
  !
  Zj=  M%Chg
  Zj2= Zj*Zj
  !
  !------------------- compute Omega (here W) and related derivatives --
  !------------------------ (equival't to routine Omeg92 of SUPCRT92) --
  !
  !-------------------------------- for NEUTRAL AQUEOUS SPECIES OR H+ --
  if ((Zj==0) .or. (M%name== 'H+')) then
  
    ReHRef=  3.082D0
    RejRef=  Zero
    W=       M%WRef
    dWdP=    Zero
    dWdT=    Zero
    d2WdT2=  Zero

  !------------------------ for CHARGED AQUEOUS SPECIES OTHER THAN H+ --
  else
   
    RejRef= Zj2 / (M%WRef/eta + Zj/(3.082d0 + gShockRef))
    ! actually, gShockRef is set to Zero ...
    ReHRef= 3.082D0
    Rej=    RejRef + ABS(Zj) *pW%GShok
    ReH=    ReHRef + pW%GShok
    
    W=             eta*Zj2/Rej         - eta*Zj/ReH
    X1=        ABS(Zj)*Zj2/Rej/Rej     - Zj/ReH/ReH
    X2=            Zj2*Zj2/Rej/Rej/Rej - Zj/ReH/ReH/ReH
    
    dWdP=  -eta *X1 *pW%dGshdP
    dWdT=  -eta *X1 *pW%dGshdT
    d2WdT2= eta *X2 *pW%dGshdT**2 *2.0D0 &
    &     - eta *X1 *pW%d2GshdT2
    
  end if
  !---------------------------------------------------/ compute Omega --
  !
  !if(M%chg/=0) then; M%RejRef= M%chg*M%chg /(M%WRef/Eta + M%chg/3.082D0)
  !else             ; M%RejRef= Zero
  !end if
  !
  Vs=-w*Q_ + Z_*dwdP ! solvation part
  VR= Vs         &
    + M%A1       & !non-solvation part
    + M%A2/PpPsi &
    + M%A3/TmTht &
    + M%A4/PpPsi/TmTht
  !
  Ss= w*Y_ - Z_*dwdT - M%wref*YPrTr !solvation part
  SR= Ss                 &
    + M%S0_              &
    + M%C1  *log(T/Tref) &
    - M%C2               &
      *(One/TmTht - One/T0mTht + log(Tref*TmTht/T/T0mTht)/Theta) &
      /Theta &
    + (M%A3*(P-Pref) + M%A4*log(PpPsi/P0pPsi)) &
       /TmTht/TmTht
  !
  CpS= w*T*X_ &
    + T*Two*Y_*dwdT &
    - T*Z_*d2wdT2
  CpR= M%C1 &
    + M%C2/TmTht/TmTht &
    - (M%A3*(P-Pref) + M%A4*log(PpPsi/P0pPsi)) *Two *T /TmTht**3
  CpR= CpS + CpR
  !
  Hs= w *Z_                &
    + w *Y_ *T             &
    - T *Z_ *dwdT          &
    + M%wref*(ZPrTr + One) &
    - M%wref*Tref*YPrTr
  HR= M%H0R &
    + M%C1*(T-Tref) &
    - M%C2*(One/TmTht - One/T0mTht) &
    + M%A1*(P-Pref) &
    + M%A2*log(PpPsi/P0pPsi) &
    +(M%A3*(P-Pref) + M%A4*log(PpPsi/P0pPsi))*(T + TmTht)/TmTht/TmTht
  HR= Hs + HR
  !GZterm
  Gs= w*Z_ &
    - M%wref *(-ZPrTr - One) &
    + M%wref *YPrTr *(T-Tref)
  GR= M%G0R &
    - M%S0_ *(T-Tref) &
    - M%C1  *(T*log(T/Tref)-T+Tref) &
    + M%A1  *(P-Pref) &
    + M%A2  *log(PpPsi/P0pPsi) &
    - M%C2  *( (T-Tref)/T0mTht - T*log(Tref*TmTht/T/T0mTht)/Theta ) /Theta &
    +(M%A3*(P-Pref) + M%A4*log(PpPsi/P0pPsi))/TmTht
  GR= Gs + GR
  !
  S%V0=  VR *CalToJoule*1.0D-5 !Cal/Bar to m3  
  !! S%Vs=  Vs *CalToJoule*1.0D-5 !Cal/Bar ro m3
  S%S0=  SR *CalToJoule 
  !! S%Ss=  Ss *CalToJoule  
  S%Cp0= CpR*CalToJoule
  !! S%CpS= CpS*CalToJoule
  S%H0=  HR *CalToJoule
  !! S%Hs=  Hs *CalToJoule
  GR=    GR *CalToJoule
  !! S%Gs=  Gs *CalToJoule
  !
  if(.not. DtbConv_Benson) GR= GR - Tref*M%S0ele
  S%G0rt= GR /R_jk /T  
  !!S%logK= - S%G0rt/Ln10
  !M%dlogKdT=   M%Haqs / (2.302585d0 * R_ck * T**2) 
  !M%dlogKdP= -0.23901488d-1 * M%VR / (2.302585d0 * R * T)
  !
  S%WeitKg= M%WeitKg
  !
  !!S%VMol=S%V0
  !
  S%AquSize= 3.72D0
  ! 3.72 A is default value corresponding to Na+=1.91 + Cl-=1.81
  !
  if(M%name=='H+') RejRef= ReHRef
  !
  ! following calculations of size parameters for aqu'species
  ! are taken from TOUGHREACT user's guide, appendix F
  ! (derived from Reed, 1982), applying Helgeson eqn 125,
  ! they assume that NaCl is the dominant electrolyte,
  ! which explains values 1.91 (= RejRef for Na+) and 1.81 (= RejRef for Cl-)
  !
  ! cation -> use RejRef of Cl-, i.e. 1.81 Angstrom
  if(Zj>0) S%AquSize= &
  & 2.D0*(Nu_Cation*RejRef + abs(Zj) *Rej_Anion) &
  &     /(Nu_Cation        + abs(Zj))
  !
  ! anion -> 
  if(Zj<0) S%AquSize= &
  & 2.D0*(Nu_Anion *RejRef + abs(Zj) *Rej_Cation) &
  &     /(Nu_Anion         + abs(Zj))
  !
  if(M%name=='H+') S%AquSize= 9.0D0
  !
  !~ !! in the case of CaCl2-dominant medium:
  !~ if(Zj>0) & ! cation -> use Rej(Cl-)
  !~ & S%AquSize= 2.D0*(RejRef + Rej_Anion *abs(Zj))/(One + abs(Zj))
  !~ !
  !~ if(Zj<0) & ! anion -> use Rej(Ca+2)
  !~ & S%AquSize= 2.D0*(2*RejRef + Rej_Ca *abs(Zj))/(2 + abs(Zj))
  !
  return
end subroutine DtbAquHkf_Calc

subroutine DtbAquHkf_CalcThr( &
& M, & !IN:  database model
& pW,& !IN:  T,P and dielectric properties of Solvent
& S)   !OUT: species
!--
!-- derived from AQUA, by CdeCapitani, 2004, with modifications
!-- calculates Hkf aqu.species using (H,M) data
!-- (DtbAquHkf_Calc uses (G,M))
!-- before call to DtbAquHkf_CalcThr, must have computed prop' of H2O -> pW
!--
  use M_Dtb_Const,only: CalToJoule, Tref, Pref, T_CK, R_jk, DtbConv_Benson
  use M_Numeric_Const,   only:Ln10
  use M_T_DtbH2OHkf
  use M_T_Species, only: T_Species
  !
  type(T_DtbAquHkf),  intent(in) :: M
  type(T_DtbH2OHkf),intent(in) :: pW
  type(T_Species),    intent(out):: S
  !
  real(dp):: P,T
  real(dp):: &
  GShok,dGShdP,dGShdT,d2GShdT2,& ! from T_DtbH2OHkf
  Z_,Y_,Q_,X_,&                  ! from T_DtbH2OHkf
  Zj,RejRef,Rej,ReH,WjRef,Wj,&   ! species properties
  !ReRef= "effective electrostat. radius at Tref,Pref"
  PpPsi,P0pPsi,TmTht,T0mTht,&
  X1,X2, &
  fHCp, fSCp, fGVol, fGCp, fG0
  !
  ! results for a given T,P, 
  real(dp):: Vs,Ss,Cps,Hs,Gs !solvation term
  real(dp):: GR,HR,SR,VR,CpR !total
  !
  ! INPUT= 
  ! AQ1:  K1,K2,K8, K9, K7
  ! Zj, RejRef=rx+Z*0.94 for Z>0
  ! AQ2:  D1,D2,D3,D4
  !
  P=pW%Pbar          ;  T=pW%TdgK
  PpPsi= P  + Psi    ;  P0pPsi= Pref + Psi
  TmTht= T  - Theta  ;  T0mTht= Tref - Theta
  GShok=  pW%GShok   ;  dGShdP=   pW%dGShdP
  dGShdT= pW%dGShdT  ;  d2GShdT2= pW%d2GShdT2
  !
  Q_= pW%Q ! Q_= pW%dEpsdP /pW%Eps /pW%Eps
  X_= pW%X ! X_= (pW%d2EpsdT2 -2.0D0/pW%Eps*pW%dEpsdT**2)/pW%Eps/pW%Eps
  Y_= pW%Y ! Y_= pW%dEpsdT /pW%Eps /pW%Eps
  Z_= pW%Z ! Z_= One /pW%Eps - 1 !not the original definition of Z
  
  ! RejRef= M%Wref
  ! in Original database, W_ contains Omega for neutral species,
  ! and Re,j,PrTr for charged species
  ! but, in order to use slop98.dat as input for theriak-like calculation,
  ! the relation JOH92/56 is re-introduced
  Zj=     M%Chg
  RejRef= Zero
  WjRef=  M%Wref !default values, valid for neutral species
  Wj=     M%Wref !default values, valid for neutral species
  X1=     Zero   !default values, valid for neutral specie
  X2=     Zero   !default values, valid for neutral species
  if (Zj/=Zero) then
    RejRef= Zj*Zj /(WjRef/Eta + Zj/3.082D0)   ! JOH92/56
    Rej=    RejRef  +ABS(Zj)*GShok            ! JOH92/48
    ReH=    3.082D0 +GShok
    Wj=     Eta *Zj *(Zj/Rej - One/ReH)       ! Omega_j, JOH92/55
    X1=    -Eta *( ABS(Zj**3)/Rej/Rej - Zj/ReH/ReH ) !dWdP=X1*dgdP; dWdT=X1*dgdT
    X2=     Eta *2.0D0*Zj *( (Zj/Rej)**3 - One/(ReH**3) ) !d2WdT2= dGShdT**2*X2+d2GShdT2*X1
  end if
  !
  CpR= M%C1 + M%C2/TmTht/TmTht
  ! CpS = solvation contrib'
  CpS= Wj*T*X_ &
  &  + Two*T*Y_*dGShdT*X1 & 
  &  - T*Z_*(dGShdT**2*X2+d2GShdT2*X1)
  CpR= CpR + CpS
  !
  VR= M%A1 &
    + M%A2/PpPsi &
    + M%A3/TmTht &
    + M%A4/PpPsi/TmTht
  ! Vs = solvation contrib'
  Vs= - Wj*Q_ + Z_*dGShdP*X1
  VR=   VR + Vs
  !
  fHCp= M%C1*(T-Tref) &
      - M%C2*(One/TmTht - One/T0mTht)
  Hs=   Wj*Z_ +Wj*T*Y_ -WjRef*Tref*YPrTr +WjRef*(ZPrTr +One) -T*Z_*dGShdT*X1
  HR=   M%H0R + fHCp + Hs
  !
  fSCp= M%C1*log(T/Tref) &
      - M%C2*( (Tref-T)/TmTht/T0mTht + log(Tref*TmTht/T0mTht/T)/Theta )/Theta
  Ss=   Wj*Y_ -Z_*dGShdT*X1 -WjRef*YPrTr
  SR=   M%S0_ + fSCp + Ss
  !
  fG0=   M%H0R - T* M%S0_
  fGCp=  fHCp-T*fSCp
  fGVol= M%A1*(P-Pref) &
       + M%A2*log(PpPsi/P0pPsi) &
       + One/TmTht*( M%A3*(P-Pref) + M%A4*log(PpPsi/P0pPsi) )
  Gs=    Wj*Z_ +WjRef*(ZPrTr+1) +WjRef*YPrTr*(T-Tref) !-WjRef*E01 replaced by +WjRef*(ZPrTr+1) 
  GR=    fG0+fGCp+fGVol+Gs
  !
  !--- unit conversions
  !! S%Vs=  Vs *CalToJoule*1.0D-5 ! Cal/Bar to m3
  !! S%Ss=  Ss *CalToJoule  
  !! S%CpS= CpS*CalToJoule
  !! S%Hs=  Hs *CalToJoule
  !! S%Gs=  Gs *CalToJoule
  !
  S%V0=  VR *CalToJoule*1.0D-5 ! *1E-5= Joule/Bar to m3  
  S%S0=  SR *CalToJoule 
  S%Cp0= CpR*CalToJoule
  S%H0=  HR *CalToJoule
  GR=    GR *CalToJoule
  !---/ unit conversions
  !
  S%WeitKg=M%WeitKg
  !
  !!!conversion from Berman-Brown to Benson-Helgeson Convention!!!
  !! S%S0Ele= M%S0Ele
  if(DtbConv_Benson) GR= GR + Tref*M%S0Ele
  S%G0rt= GR /R_jk /T 
  !!S%logK= -S%G0rt/Ln10
  !
  if(M%name== 'H+') then
    !! S%AquSize= 9.0D0  ! why ???
    S%AquSize= 3.08D0 !! Helgeson'81, p1304
  else
    S%AquSize= Rej !Ref !*1.D-10
  end if
  !
  if (iDebug>0) then
    write(fTrc,'(A,A1,2(F7.1,A1))',advance="no") M%Name,T_,T-T_CK,T_,P,T_
    write(fTrc,'(5(F15.4,A1))') S%V0,T_, S%S0,T_, S%Cp0,T_, S%H0,T_, -S%G0rt/Ln10,T_ 
  end if
  !
  return
end subroutine DtbAquHkf_CalcThr

end module M_T_DtbAquHkf

