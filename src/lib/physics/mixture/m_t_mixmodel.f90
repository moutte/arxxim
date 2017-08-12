module M_T_MixModel
!--
!-- implementation of solution models, and related potentials, activities, ...
!-- developed for minerals,
!-- to be used also for other mixtures (fluids, melts) that follow a mole fraction model
!--
  use M_Kinds
  implicit none
  private
  !
  public:: T_MixModel,T_Margul
  public:: MaxSite,MaxPole,MaxAtom,MaxMulti,MaxMarg
  !
  public:: MixModel_Param_Update
  public:: MixModel_Index
  public:: MixModel_Zero
  public:: MixModel_Name
  public:: MixModel_XPoleToXSite

  public:: MixModel_Activities
  public:: MixModel_GibbsMixRT
  public:: MixModel_GibbsIdeal

  public:: MixModel_Site_ActivIdeal

  !public:: MixModel_Site_LnActivsIdeal
  !public:: MixModel_Pole_LnActivsIdeal

  !public:: MixModel_Site_SConf
  !public:: MixModel_Pole_SConf

  !public:: MixModel_Site_XsMargules
  !public:: MixModel_Pole_XsMargules

  !public:: MixModel_Site_LnGammaMargules
  !public:: MixModel_Pole_LnGammaMargules

  !public:: MixModel_Margules_Monome
  !
  !public:: MixModel_SConf
  !public:: MixModel_XsMargules
  !
  !real(dp),parameter::& !from ModConst
  !& R_jK= 8.314510D0,&  !R constant, in J/K/Mole, RobieHemingway (from COdata)
  !& Pref= 1.0000D0,&  !bar
  !& Tref= 298.150D0    !in degK,=25°C
  
  integer,parameter,public:: & !
  & Mix_Molecular= 1, & !
  & Mix_Felspar=   2, & !
  & Mix_Site=      3, & !
  & Mix_Special=   4, & !
  !
  & Mix_Exchange= 11, & !
  & Mix_Surface=  12

  !-----------------------------fixed dimensions for all T_MixModel's --
  !
  integer,parameter:: & !fixed dimensions for T_MixModel
  & MaxPole=  6,   & !max number of endmembers (or components) to define a mixture
  & MaxSite=  6,   & !max number of sites in the mixture
  & MaxAtom= 16,   & !max dimension of composition vector of mixture, used in T_MixPhase
  & MaxMulti= 6,   & !max multiplicity
  & MaxMarg= 24      !max nr of Margules terms
  !
  !----------------------------/fixed dimensions for all T_MixModel's --
  !
  !---------- type for Margules parameter for up to Quaternary mixing --
  !
  !! type:: T_Margul
  !!   integer:: NPoleMarg
  !!   !number (2:4) of "poles" (end-members) involved in the Margules term
  !!   !most commonly, NPole=2 (binary interaction), sometimes NPole=3
  !!   integer,dimension(1:4)::&
  !!     vIPole, &
  !!     ! for molecular models,
  !!     !   vIPole lists the indices, in the e-m list of the mixture,
  !!     !   of the e-m ("pole") involved in Margules term
  !!     ! for site models,
  !!     !   vIPole lists the indices, in the element list of the site,
  !!     !   of the elements involved in Margules term
  !!     !
  !!     vPower !, & !
  !!     !power(i)= power of "pole" IPole(i) in the the Margules term
  !!     !ISite    !for site models
  !! end type T_Margul
  
  type:: T_Margul
    integer,dimension(1:MaxAtom):: vDegree
    ! new format,
    !   vDegree(:) contains all degrees,
    !     listed in the same order as
    !     the e-m list in case of molecular model,
    !     or as the site fractions in case of site model
    !   if a e-m or a site fraction is not involved, its degree is zero
    ! -> vDegree(:) contains information that was contained in vIPole and in vPower
    real(dp)::WG,WH,WS,WV,WCp
    integer ::Kohler
    !caveat!! deCapitani allows the Kohler parameter to take fractional values
  end type T_Margul
  !
  !---------------------------------------------------------------------
  !
  !WG=   WH &
  !&   + WCP*(T-T0) &
  !&   -(WS +WCP*log(T/T0))*T &
  !&   + WV *P
  !WCp commonly Zero -> WG= WH -T*WS +P.WV
  !
  !---------------------------------- about using Margules parameters --
  != one Margules term M of the excess G calculated as
  !  Y=M%WG
  !<old version>
  !  do I=1,M%NPoleMarg !NPole= number of poles involved in term M
  !    Y=Y*X(M%iPol(I))**(M%Power(I))
  !  !! IPole=list of indices (in the solution) of the poles involved in M
  !  end do
  !  example:
  !    term 112 -> Power(1)=2, Power(2)=1 -> W112 *X(IPol(1))**2 *X(Ipol(2))
  !
  !<new version>
  !  do I=1,S%NPole != number of poles in the mixture
  !    if(M%vDegree(I)>0) Y=Y*X(I)**M%vDegree(I))
  !    !! no more use of vIPole(:)
  !  end do
  !
  !--------------------------------------------------------------------
  !
  !-----------------------------------------------------------T_MixModel
  type:: T_MixModel

    character(len=23):: Name    !same length as name in tMinThr
    character(len=3) :: Typ     !"MIN","GAS",...
    ! character(len=15):: Model   !type of EoS: IDEAL, MOLECULAR, SITE, FELSPAR, GUGGENHEIM, ...
    integer          :: Model

    integer:: NPole             !Nr of EndMembers, must be <=MaxPole
    integer:: NSite             !Nr of Sites, must be <=MaxSite
    integer:: NMarg             !Nr of Margules terms

    ! logical:: vHasPole(1:MaxPole)

    character(len=23):: vNamPole(1:MaxPole)
    ! names of end-members
    ! (same name must be found also in pure species database)

    integer:: vIPole(1:MaxPole)
    ! vIPole(i)= index of end member i in current species list
    !   species list is generally vSpc in M_Global_Vars,
    !   thus,
    !     vIPole(i)   is Species_Index(vSpc, vNamPole(:)
    !     vNamPole(i) is vSpc(vIPole(i))%Name

    integer:: vMulti(1:MaxSite)
    ! site multiplicity
    ! for molecular models, multiplicity is vMulti(1)
    ! vMulti(i)<=MaxMulti forall i
    !
    !! real(dp):: GuggA0,GuggA1
    !! ! for Guggenheim binary assymetric model (e.g. Mg-calcite)
    !
    ! fields GuggA0,GuggA1 have been deleted,
    ! because implementation of a "Guggenheim mode" for a binary mixture
    ! is possible using fields already available in T_Margul
    !
    ! G_Xs= x *(1-x) *( A0 + A1*(x - (1-x))
    ! equivalent with Margules:
    ! G_XsMix= A0*X1*X2 + A1*X1*X1*X2 -A1*X1*X2*X2
    ! >>
    ! W12= A0  ;  W112= A1  ;  W221= -A1
    !
    type(T_Margul):: vMarg(1:MaxMarg)
    ! array of Margules parameters
    !
    !! vIMargSite still useful, although information is also contained in vPower
    integer:: vIMargSite(1:MaxMarg)
    ! for molecular (or one-site models), vIMargSite(1:NMarg)= 1
    ! for multisite models, vIMargSite(i) tells which site vMarg(i) is for (??)

    !----------------------------------------------------for site models 
    integer:: NAtom
    ! dimension of composition vector, <=MaxAtom
    !
    real(dp):: vPoleCoeff(1:MaxPole)
    ! normalization factors for site models
    !
    character(len=3)::vNamSite(1:MaxSite)
    ! names of sites, e.g. M1_, M2_, M23, A__, OH_, ...
    !
    character(len=6):: vNamAtom(1:MaxAtom)
    ! 6 chars= 3 for element name, 3 for site name
    ! names of site fractions, e.g. Mg_M1_, Mg_M2_, Si_T__, OH_OH_, ...
    !
    integer:: vAtomMulti(1:MaxAtom)
    ! gives the multiplicity of the site containing the atom
    integer:: vIAtomSite(1:MaxAtom)
    ! gives the index of the site containing the atom
    integer:: tPoleAtom(1:MaxPole,1:MaxAtom)
    ! table of end-members w/r site fractions
    ! example
    ! SITE
    !   A    1  K__Na_     !iSite=1; K:iEle=1, Na:iEle=2
    !   M2A  1  Al_Mg_Fe_
    !   T1   2  Al_Si_
    ! end
    ! >> composition vector will contain
    !   (xK_A, xNa_A, xAl_M2A, xMg_M2A, xFe_M2A, xAl_T1, xSi_T1)
    ! list of poles and site occupancy table:
    ! POLE
    !   !    A(1) M2A(1) T1(2)  ! M2B(1) T2(2)
    !   Mu   K    Al     AlSi   ! Al     SiSi
    !   Pa   Na   Al     AlSi   ! Al     SiSi
    !   MCel K    Mg     SiSi   ! Al     SiSi
    !   FCel K    Fe     SiSi   ! Al     SiSi
    ! end
    ! ->>
    ! tPoleAtom=
    !      !K_A   Na_A   Al_M2A   Mg_M2A   Fe_M2A   Al_T1   Si_T1
    ! Mu    1     0      1        0        0        1       1
    ! Pa    0     1      1        0        0        1       1
    ! MCel  1     0      0        1        0        0       2
    ! FCel  1     0      0        0        1        0       2
    ! vIAtomSite=
    !       1     1      2        2        2        3       3
    ! vAtomMulti=
    !       1     1      1        1        1        2       2
    !
    !---------------------------------------------------/for site models
    !
  end type T_MixModel
  !----------------------------------------------------------/T_MixModel
  !
  !!!-----------------------------------------------------------examples
  !!!
  !!!  SITE
  !!!    M2 2 MG_FE_
  !!!    M1 1 MG_FE_AL_ !may have up to MaxElSite elements
  !!!    T2 2 AL_SI_
  !!!  endSITE
  !!!  POLE
  !!!    !__________       M2____    M1_        T_____
  !!!    !possible=        MG_FE_    MG_FE_AL_  AL_SI_
  !!!    !vMulti=          2         1          2
  !!!    !phlogopite       Mg2       Mg         AlSi3
  !!!    PHLOGOPITE        MG_MG_    MG_        AL_SI_
  !!!    ANNITE            FE_FE_    FE_        AL_SI_
  !!!    SIDEROPHYLLITE    FE_FE_    AL_        AL_AL_
  !!!    EASTONITE         MG_MG_    AL_        AL_AL_
  !!!  endPOLE
  !!!  .../...
  !!!endSOLUTION
  !!!
  !!!
  !!!
  !!!-> example of sol.composition:
  !!! COMPOSITION
  !!!   M2 1.00 1.00 !sum=2.00, =vMulti(1)
  !!!   M1 0.50 0.50 !sum=1.00, =vMulti(2)
  !!!   T2 1.00 1.00 !sum=2.00, =vMulti(3)
  !!! end
  !!!
  !!!----------------------------------------------------------/examples
  !

contains

subroutine MixModel_Zero(S)

  type(T_MixModel),intent(out)::S
  
  S%Name=     "XXX"
  S%Model=    1 ! "MOLECULAR"
  S%NPole=    0
  S%NSite=    1
  S%vMulti(1)=1
  S%nMarg=    0
  !S%vHasPole(:)= .false.
  S%vNamPole(:)= "Z"
  S%vIPole(:)=   0
  S%vMulti(:)=   1
  S%NAtom= 1
  S%vNamSite(:)=     "Z"
  S%vNamAtom(:)=     "EleSit"
  S%vPoleCoeff(:)=   One
  S%vAtomMulti(:)=   1
  S%vIAtomSite(:)=    1
  S%tPoleAtom(:,:)=  0
  
end subroutine MixModel_Zero

integer function MixModel_Index(V,Str)
!-- position of solution model named Str in V(1:nMixModel)
!-----------------------------------------------------------------------
  type(T_MixModel),dimension(:),intent(in):: V
  character(*),                 intent(in):: Str
  integer     ::I
  MixModel_Index=0
  !---------------------------------------------------------------------
  if(size(V)==0) return
  !
  I=0
  do
    I=I+1
    if(trim(Str)==trim(V(I)%Name)) then
      MixModel_Index=I
      exit
    end if
    if(I==size(V)) exit
  end do
  !
  return
end function MixModel_Index

subroutine MixModel_Name(S,Str,Err)
  type(T_MixModel),intent(in) :: S
  character(len=*),intent(out):: Str
  logical,         intent(out):: Err
  !
  Err= .false.
  !
  select case(S%Model)
  case(Mix_Molecular)  ; Str= "MOLECULAR"
  case(Mix_Felspar)    ; Str= "FELSPAR"
  case(Mix_Site)       ; Str= "SITE"
  case(Mix_Special)    ; Str= "SPECIAL"
  !
  case(Mix_Exchange)   ; Str= "EXCHANGE"
  case(Mix_Surface)    ; Str= "SURFACE"
  !
  case default
    Str= "NONE"
    Err= .true.
  end select
  !
  return
end subroutine MixModel_Name

subroutine MixModel_XPoleToXSite(MM,vX,vY)

  type(T_MixModel),intent(in) :: MM
  real(dp),        intent(in) :: vX(:) ! end member mole fractions
  real(dp),        intent(out):: vY(:) ! site fractions
  !---------------------------------------------------------------------
  integer:: iP,iEl,J
  !
  vY(:)= Zero
  do iP=1,MM%NPole
    do iEl=1,MM%NAtom
      J= MM%tPoleAtom(iP,iEl)
      !print *,J
      if(J /= 0) vY(iEl)= vY(iEl) + vX(iP) *real(J) /real(MM%vAtomMulti(iEl))
      !print *,vY
    end do
    !pause
  end do
  !
  return
end subroutine MixModel_XPoleToXSite

real(dp) function MixModel_Site_ActivIdeal( &
& MM,   & !IN= Mix.model, generally =vMixModel(Fas%iModel)
& IPol, & !IN= index of end-member in models
& vXatom) !IN= phase composition
!--
!-- for site mixing models
!-- returns ideal activity of end-member IPol,
!-- following solution model MM
!--
  type(T_MixModel),intent(in):: MM
  real(dp),        intent(in):: vXAtom(:)
  integer,         intent(in):: IPol
  !---------------------------------------------------------------------
  real(dp):: Y
  integer :: iEl
  !---------------------------------------------------------------------
  Y= One
  !
  do iEl= 1,MM%NAtom
    !if(MM%tPoleAtom(iPol,iEl)/=0) Y= Y *vXAtom(iEl)**MM%tPoleAtom(iPol,iEl)
    if(MM%tPoleAtom(iPol,iEl)/=0) Y= Y *vXAtom(iEl)**MM%vAtomMulti(iEl)
  end do
  !
  MixModel_Site_ActivIdeal= Y *MM%vPoleCoeff(iPol)
  !
  !if(iDebug>1) write(fTrc,'(I3,1X,G15.6,1X,A)') &
  !& iPol,MixModel_Site_ActivIdeal,trim(MM%vNamPole(iPol))
  !
end function MixModel_Site_ActivIdeal

subroutine MixModel_Pole_LnActivsIdeal( & !
& MM,     &  !IN= Mix.model
& vX,     &  !IN= phase composition
& vLpole, &  !IN
& vLnAct)    !OUT
!-----------------------------------------------------------------------
!-- for end-member mixing models
!-- returns log(ideal activity) of all 'active' end-members
!-----------------------------------------------------------------------
  type(T_MixModel),intent(in) :: MM
  real(dp),        intent(in) :: vX(:)
  logical,         intent(in) :: vLpole(:)
  real(dp),        intent(out):: vLnAct(:)
  !---------------------------------------------------------------------
  integer:: iP
  !---------------------------------------------------------------------
  vLnAct(:)= Zero
  !
  select case(MM%Model) ! (trim(MM%Model))
  !
  case(Mix_Molecular) ! ("IDEAL","POLE","MOLECULAR")
    do iP=1,MM%NPole
      if(vLPole(iP)) &
      & vLnAct(iP)= vLnAct(iP) + log(vX(iP))*MM%vMulti(1)
    end do
  !
  case(Mix_Felspar) ! ("FELSPAR")
    if(vLPole(1)) vLnAct(1)= log(vX(1)) !+log(1-vX(3)*vX(3))
    if(vLPole(2)) vLnAct(2)= log(vX(2)) !+log(1-vX(3)*vX(3))
    if(vLPole(3)) then
      vLnAct(3)= log(vX(3))+2*log((One+vX(3))/2.0D0)
      if(vLPole(1)) vLnAct(1)= log(vX(1)) +log(One-vX(3)**2)
      if(vLPole(2)) vLnAct(2)= log(vX(2)) +log(One-vX(3)**2)
    end if
  !
  end select
  !
end subroutine MixModel_Pole_LnActivsIdeal

subroutine MixModel_Site_LnActivsIdeal( & !
& MM,     &  !IN= Mix.model
& vXatom, &  !IN= phase composition
& vLpole, &  !IN
& Ok,     &  !OUT
& vLnAct)    !OUT
!-----------------------------------------------------------------------
!-- for site mixing models
!-- returns log(ideal activity) of all 'active' end-members
!-----------------------------------------------------------------------
  type(T_MixModel),intent(in) :: MM
  real(dp),        intent(in) :: vXAtom(:)
  logical,         intent(in) :: vLpole(:)
  logical,         intent(out):: Ok
  real(dp),        intent(out):: vLnAct(:)
  !---------------------------------------------------------------------
  real(dp):: Y
  integer :: iP,iEl,J
  !---------------------------------------------------------------------
  Ok= .true.
  !
  do iP= 1,MM%NPole

    if(vLPole(iP)) then

      Y= Zero
      do iEl= 1,MM%NAtom
        J= MM%tPoleAtom(iP,iEl)
        if(J/=0) then
          if(vXAtom(iEl)>Zero) then
            Y= Y + log(vXAtom(iEl)) *MM%vAtomMulti(iEl) !*J
          else
          !! normally, should never happen, if vLPole is .true. ...
            Ok= .false.
          end if
        end if
      end do
      
      ! add correction term, to have activ=1 for pure end member
      vLnAct(iP)= Y + log(MM%vPoleCoeff(iP))

    end if

  end do
  !
end subroutine MixModel_Site_LnActivsIdeal

real(dp) function MixModel_Site_SConf( & !
& MM,   & !IN, mixing model
& vXatom)    !IN, composition
!-----------------------------------------------------------------------
!-- compute configurational entropy for site mixing --
!-----------------------------------------------------------------------
  use M_Dtb_Const,only: R_jk
  !
  type(T_MixModel),intent(in):: MM
  real(dp),        intent(in):: vXatom(:)
  !---------------------------------------------------------------------
  real(dp):: Sconf,X
  integer :: iEl
  !---------------------------------------------------------------------
  SConf= Zero
  !
  ! sum of x*ln(x) on the different sites,
  ! taking account also of site multiplicity
  do iEl=1,MM%NAtom
    X= vXAtom(iEl)
    if(X > Zero) &
    & Sconf= Sconf + X *log(X) *MM%vAtomMulti(iEl)
  end do
  !
  MixModel_Site_Sconf= -R_jk*Sconf
  !
end function MixModel_Site_SConf

real(dp) function MixModel_Pole_SConf( & !
& MM,   & !IN, mixing model
& vLpole,  & !IN, pole is present / absent
& vX)        !IN, composition
!--
!-- compute configurational entropy for end member mixing --
!--
  use M_Dtb_Const,only: R_jk
  !
  type(T_MixModel),intent(in):: MM
  logical,         intent(in):: vLpole(:)
  real(dp),        intent(in):: vX(:)
  !
  real(dp):: Sconf
  integer :: iP
  !
  SConf= Zero
  !
  select case(MM%Model) ! (trim(MM%Model))

  case(Mix_Molecular) ! ("IDEAL","POLE","MOLECULAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR"
    do iP=1,MM%NPole
      if(vLPole(iP)) Sconf= Sconf + vX(iP)*log(vX(iP))
    end do
    Sconf= SConf *MM%vMulti(1)
    !
  case(Mix_Felspar) ! ("FELSPAR")
  !! "FELSPAR" is actually a special case of "MOLECULAR" mixing:
  !! configurational entropy with "Al-avoidance"
    ! ALBITE
    if(vLPole(1)) Sconf= Sconf + vX(1)*log(vX(1))
    ! K-FELSPAR
    if(vLPole(2)) Sconf= Sconf + vX(2)*log(vX(2))
    ! ANORTHITE
    if(vLPole(3)) then
      Sconf= Sconf + vX(3) *( log(vX(3)) +log((One +vX(3))/Two) *Two )
      if(vLPole(1)) Sconf= Sconf + vX(1)*log(One -vX(3)*vX(3))
      if(vLPole(2)) Sconf= Sconf + vX(2)*log(One -vX(3)*vX(3))
    end if
    !
  end select
  !
  MixModel_Pole_Sconf= -R_jk*Sconf
  !
end function MixModel_Pole_SConf

real(dp) function MixModel_SConf( & !
& MM,      & !IN, mixing model
& vLPole,  & !IN
& vX)        !IN
!--
!-- compute configurational entropy --
!--
  type(T_MixModel),intent(in):: MM
  logical,         intent(in):: vLPole(:)
  real(dp),        intent(in):: vX(:)
  !
  real(dp),allocatable:: vXAtom(:)
  !
  select case(MM%Model) ! (trim(MM%Model))

  case(Mix_Molecular,Mix_Felspar) ! ("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR"
    MixModel_SConf= MixModel_Pole_SConf(MM,vLPole,vX) ! vX is Fas%vXPole
  !
  case(Mix_Site) ! ("SITE")
    allocate(vXAtom(MM%NAtom))
    vXAtom(:)= vX(:)
    MixModel_SConf= MixModel_Site_SConf(MM,vXAtom) ! vX is Fas%vXAtom
    deallocate(vXAtom)
  !
  end select
  !
end function MixModel_SConf

real(dp) function MixModel_Margules_Monome(M,vX)
  type(T_Margul),intent(in):: M
  real(dp),      intent(in):: vX(:)
  !
  real(dp):: Y
  integer :: i
  !
  Y= M%WG
  do i=1,size(vX)
    if(M%vDegree(i)>0) &
    & Y= Y *vX(i)**M%vDegree(i)
  end do
  !
  MixModel_Margules_Monome= Y
  !
end function MixModel_Margules_Monome

real(dp) function MixModel_Pole_LnGammaMargules( & !
& MM,      & !IN, mixing model
& iP,      & !IN
& vMonome, & !
& vXpole)    !IN, composition
!--
!-- for MOLECULAR model, compute RTln(Gamma) related to Margules Terms
!-- (according to expression by deCapitani-Kirschen)
!--
  type(T_MixModel),intent(in):: MM
  integer,         intent(in):: iP
  real(dp),        intent(in):: vMonome(:),vXpole(:)
  !
  type(T_Margul)::M
  real(dp):: LnGam,S_i,Y,Z
  integer :: iM,I
  integer :: dSi_dXj
  !
  ! RT.LnGam(j)=      &
  !   sum(i=1..NMarg) &
  !   W_i1..ip . X_i1..X_ip (S_i)^(-k_i) &
  !   ( (1 - p_i + k_i) + q_j / X_j - (k_i/S_i).dSi_dXj ) )
  ! P_i = degree of monomial X_i1... = sum(M%vDegree(:))
  ! S_i = sum of mole fractions involved in term = sum(vX(M%iPol(1:M%NPole)))
  ! Q_j = number of indices, among i1..ip, equal to j,
  !     = M%vDegree(j)
  !
  LnGam=Zero
  !
  do iM=1,MM%NMarg

    M= MM%vMarg(iM)
    Y= vMonome(iM)
    Z= M%vDegree(iP) /vXpole(iP) +One -sum(M%vDegree(:))
    !
    if(M%Kohler/=0) then
      S_i= Zero
      do I=1,MM%NPole
        if(M%vDegree(I)>0) S_i= S_i +vXpole(I)
      end do
      dSi_dXj=0
      if(M%vDegree(iP)>0) dSi_dXj= 1
      !
      Y= Y / S_i**M%Kohler
      Z= Z + M%Kohler
      if(dSi_dXj/=0) Z= Z - M%Kohler *dSi_dXj/S_i
    end if
    !
    LnGam= LnGam + Y*Z
  end do
  !
  MixModel_Pole_LnGammaMargules= LnGam ! gives actually R*T*Ln(Gamma) !!
  !
end function MixModel_Pole_LnGammaMargules

real(dp) function MixModel_Site_LnGammaMargules( & !
& MM,      & !IN, mixing model
& iP,      & !IN
& vMonome, & !
& vXatom)    !IN, composition
!--
!-- for SITE model, compute RTln(Gamma) related to Margules Terms
!-- ( should be called only when vLpole(iP) is true )
!--
  type(T_MixModel),intent(in):: MM
  integer,         intent(in):: iP
  real(dp),        intent(in):: vXatom(:),vMonome(:)
  !
  type(T_Margul)::M
  real(dp):: X != RT.ln(gamma)
  integer :: iM,J,iEl,iS
  !
  ! n * RT.LnGam_m=   &
  !   sum(i=1..NMarg) &
  !   [ W_i1..ip . X_i1..X_ip  * ( Q_m /X_m - 2 ) ] !for degree 3 Margules !!
  ! n   = number of sites on which mixing occurs
  ! q_m = number of indices, among i1..ip, equal to m,
  !     = M%vDegree(m)
  !
  !! LnAct= Zero
  X= Zero
  !
  do iEl= 1,MM%NAtom
    J= MM%tPoleAtom(iP,iEl)
    if(J/=0) then
      !! LnAct= LnAct + log(vXAtom(iEl)) *J
      do iM=1,MM%NMarg
        iS= MM%vIMargSite(iM)
        if( iS == MM%vIAtomSite(iEl) ) then
          M= MM%vMarg(iM)
          X=   X &
          &  + vMonome(iM) &
          &    *(M%vDegree(iEl) /vXatom(iEl) +One -sum(M%vDegree(:))) &
          &    /real(MM%vAtomMulti(iEl))
        end if
      end do
    end if
  end do
  !
  MixModel_Site_LnGammaMargules= X ! is R*T*ln(Gamma) !!
  !
end function MixModel_Site_LnGammaMargules

real(dp) function MixModel_Pole_XsMargules( & !
& MM,     & !IN, mixing model
& vLPole, & !IN
& vXpole)   !IN, composition
!--
!-- for MOLECULAR model
!--
  type(T_MixModel),intent(in):: MM
  logical,         intent(in):: vLPole(:)
  real(dp),        intent(in):: vXpole(:)
  !
  type(T_Margul):: M
  real(dp):: Y,G_XSMix
  integer :: iM,iP
  !
  G_XSMix= Zero
  !
  do iM=1,MM%NMarg ! number of Margules terms
    !
    M= MM%vMarg(iM)  ! the Margules structure
    Y= M%WG          ! the interaction energy for this Margules term
    !
    do iP=1,MM%NPole
      if(M%vDegree(iP)>0) then != pole iP involved in monome Y
        if(vLPole(iP)) then
          Y= Y *vXPole(iP)**M%vDegree(iP)
        else
          Y= Zero
        end if
      end if
    end do
    !
    G_XSMix= G_XSMix +Y
  end do
  !
  MixModel_Pole_XsMargules= G_XSMix
  !
  return
end function MixModel_Pole_XsMargules

real(dp) function MixModel_Site_XsMargules( & !
& MM,    & !IN, mixing model
& vXatom)  !IN
!--
!-- for SITE model
!--
  type(T_MixModel),intent(in):: MM
  real(dp),        intent(in):: vXAtom(:)
  !
  type(T_Margul):: M
  real(dp):: Y,G_XSMix
  integer :: iM,iEl
  !
  G_XSMix= Zero
  !
  do iM=1,MM%NMarg ! number of Margules terms
    !
    M= MM%vMarg(iM)  ! the Margules structure
    Y= M%WG          ! the interaction energy for this Margules term
    !
    do iEl=1,MM%NAtom
      if(M%vDegree(iEl)>0) & != site fraction iEl is involved in monome Y
      & Y= Y *vXAtom(iEl)**M%vDegree(iEl)
    end do
    !
    G_XSMix= G_XSMix +Y
  end do
  !
  MixModel_Site_XsMargules= G_XSMix
  !
  return
end function MixModel_Site_XsMargules

real(dp) function MixModel_XsMargules( & !
& MM,      & !IN, mixing model
& vLPole,  & !IN
& vX)        !IN
!--
  type(T_MixModel),intent(in):: MM
  logical,         intent(in):: vLPole(:)
  real(dp),        intent(in):: vX(:)
  !
  MixModel_XsMargules= Zero
  !
  if(MM%NMarg==0) return
  !
  select case(MM%Model) ! (trim(MM%Model))
  !
  case(Mix_Molecular,Mix_Felspar) ! ("IDEAL","POLE","MOLECULAR","FELSPAR")
  ! should consider "IDEAL","POLE" as obsolete,
  ! and recommand using "MOLECULAR" or "FELSPAR"
    MixModel_XsMargules= &
    & MixModel_Pole_XsMargules(MM,vLPole,vX) ! vX is Fas%vXPole
  !
  case(Mix_Site) ! ("SITE")
    MixModel_XsMargules= &
    & MixModel_Site_XsMargules(MM,vX) ! vX is Fas%vXAtom
  !
  end select
  !
  return
end function MixModel_XsMargules

subroutine MixModel_Param_Update(TdgK,Pbar,S)
  real(dp),        intent(in)   :: TdgK,Pbar
  type(T_MixModel),intent(inout):: S
  !
  if(S%NMarg>0) call MixModel_Margul_Wg_Calc(TdgK,Pbar,S)
  !
  return
end subroutine MixModel_Param_Update

subroutine MixModel_Margul_Wg_Calc(TdgK,Pbar,S)
!--
!-- from S%vMarg(:)%WH,WS,..., calc. values of S%vMarg(:)%WG at given P,T
!--
  use M_Dtb_Const,only: Tref
  !
  real(dp),        intent(in)   :: TdgK,Pbar
  type(T_MixModel),intent(inout):: S
  !
  type(T_Margul) ::M
  integer        ::iM
  !
  if(S%NMarg>0) then

    do iM=1,S%NMarg
      M=S%vMarg(iM)
      S%vMarg(iM)%WG= M%WH + M%WCP*(TdgK-Tref) &
      &             -(M%WS + M%WCP*log(TdgK/Tref))*TdgK &
      &             + M%WV*Pbar
    end do

  end if
  !
  return
end subroutine MixModel_Margul_Wg_Calc

subroutine MixModel_Activities( & !
& TdgK,Pbar, & ! in
& MM,        & ! in: mixing model
& vXPole,    & ! in
& vLPole,    & ! in
& Ok, Msg,   & ! out
& vLGam,     & ! out
& vLIdeal,   & ! out
& vLnAct)      ! out
!--
!-- calculate activities of end-members in phase F at given T,P
!--
  use M_Dtb_Const,only: R_jk
  use M_Trace
  use M_MixModel_Special
  !
  real(dp),        intent(in) :: Pbar, TdgK
  type(T_MixModel),intent(in) :: MM
  logical,         intent(in) :: vLPole(:)
  real(dp),        intent(in) :: vXPole(:)
  logical,         intent(out):: Ok
  character(*),    intent(out):: Msg
  real(dp),        intent(out):: vLGam(:)
  real(dp),        intent(out):: vLIdeal(:)
  real(dp),        intent(out):: vLnAct(:)
  !
  integer :: iP,iM
  real(dp):: P
  real(dp):: vMonome(MM%NMarg)
  real(dp):: vXAtom(MM%NAtom)
  !real(dp),allocatable:: vMonome(:)
  !real(dp),allocatable:: vXAtom(:)
  !
  P=Pbar !for future use ??
  !
  Ok= .true.
  Msg= "Ok"
  !
  !F%vLPole(1:MM%NPole)= vXPole(1:MM%NPole)>Zero
  !
  vLGam(:)=   Zero !default
  vLIdeal(:)= Zero !default
  !
  select case(MM%Model)
  !
  case(Mix_Special) ! "SPECIAL"
    !
    call MixModel_Special_Activities( &
    & MM%Name,TdgK,Pbar,vXPole,vLpole, &
    & vLIdeal,vLGam)
    !
  case(Mix_Molecular,Mix_Felspar) ! "IDEAL","POLE","MOLECULAR","FELSPAR"
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR" or "FELSPAR"
    !
    call MixModel_Pole_LnActivsIdeal(MM,vXPole,vLpole,vLIdeal)
    !
    !---------------------------------------------------- Margules Terms 
    if(MM%NMarg>0) then
      !allocate(vMonome(MM%NMarg))
      !
      do iM=1,MM%NMarg
        vMonome(iM)= MixModel_Margules_Monome(MM%vMarg(iM),vXPole)
      end do
      do iP=1,MM%NPole
        if(vLPole(iP)) &
        & vLGam(iP)= MixModel_Pole_LnGammaMargules(MM,iP,vMonome,vXPole) &
        &            /R_jk/TdgK
      end do
      !
      !deallocate(vMonome)
      !
    end if
    !----------------------------------------------------/Margules Terms 
    !
  case(Mix_Site)  ! "SITE"
    !
    !allocate(vXAtom(MM%NAtom))
    call MixModel_XPoleToXSite(MM,vXPole,vXAtom)
    !
    call MixModel_Site_LnActivsIdeal(MM,vXAtom,vLpole,Ok,vLIdeal)
    !
    !---------------------------------------------------- Margules Terms 
    if(MM%NMarg>0) then
      !allocate(vMonome(MM%NMarg))
      !
      do iM=1,MM%NMarg
        vMonome(iM)= MixModel_Margules_Monome(MM%vMarg(iM),vXAtom)
      end do
      do iP=1,MM%NPole
        if(vLPole(iP)) &
        & vLGam(iP)= MixModel_Site_LnGammaMargules(MM,iP,vMonome,vXAtom) &
        &            /R_jk/TdgK
      end do
      !
      !deallocate(vMonome)
      !
    end if
    !----------------------------------------------------/Margules Terms 
    !deallocate(vXAtom)
    !
  ! case default
  !   Ok= .false.
  !   Msg= trim(MM%Model)//"= invalid MM%Model in MixModel_CalcActivities"
  !   !
  end select
  !
  vLnAct(1:MM%NPole)= vLIdeal(1:MM%NPole) + vLGam(1:MM%NPole)
  !
  if(iDebug>3) then !----------------------------------------------trace
    write(fTrc,'(A)') "MixModel_CalcActivs -> X,ActIdeal,Gamma,Activ"
    write(fTrc,'(2A)') "MixModel=",MM%Name
    do iP=1,MM%NPole
      if(vLPole(iP)) then
        write(fTrc,'()')
        write(fTrc,'(A,I2,A1,A15,A1,4(A4,G15.6,A1))') &
        & "POLE",iP,              T_,&
        & trim(MM%vNamPole(iP)),  T_,&
        & "Frc=",vXPole(iP),          T_,&
        & "XId=",exp(vLIdeal(iP)),T_,&
        & "Gam=",exp(vLGam(iP)),  T_,&
        & "Act=",exp(vLnAct(iP)), T_
      end if
    end do
  end if !--------------------------------------------------------/trace
  !
end subroutine MixModel_Activities

real(dp) function MixModel_GibbsIdeal( & !
& TdgK,Pbar, & !IN
& MM,        & !IN, mixing model
& vLPole,    & !IN
& vXPole)      !IN
!--
!-- !!! todo - check implementation of site models !!!
!--
  use M_Dtb_Const,only: R_jk
  use M_MixModel_Special
  !
  real(dp),        intent(in):: TdgK,Pbar
  type(T_MixModel),intent(in):: MM
  logical,         intent(in):: vLPole(:)
  real(dp),        intent(in):: vXPole(:)
  !
  real(dp),allocatable:: vXAtom(:)
  real(dp),allocatable:: vLIdeal(:),vLGam(:)
  real(dp):: G_IdMix
  !
  ! G_IdMix: ideal part of free energy of mixing of solution MM for a given P,T
  !
  select case(MM%Model)
  !
  case(Mix_Molecular,Mix_Felspar) !("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommend using "MOLECULAR" or "FELSPAR" only
    G_IdMix= -TdgK *MixModel_Pole_SConf(MM,vLPole,vXPole)
    !
  case(Mix_Site) ! ("SITE")
    allocate(vXAtom(MM%NAtom))
    call MixModel_XPoleToXSite(MM,vXPole,vXAtom)
    !
    G_IdMix= -TdgK *MixModel_Site_SConf(MM,vXAtom)
    !
    deallocate(vXAtom)
  !
  case(Mix_Special) ! ("SPECIAL")
    allocate(vLIdeal(size(vXPole)))
    allocate(vLGam(size(vXPole)))
    !
    call MixModel_Special_Activities( &
    & MM%Name,TdgK,Pbar,vXPole,vLpole, &
    & vLIdeal,vLGam)
    !
    G_IdMix= sum(vXPole(:)*vLIdeal(:), MASK=vLPole(:)) *R_jk  *TdgK
    !
    deallocate(vLIdeal)
    deallocate(vLGam)
  !
  end select
  !
  MixModel_GibbsIdeal= G_IdMix
  !
end function MixModel_GibbsIdeal

real(dp) function MixModel_GibbsMixRT( & !
& TdgK,Pbar, & !IN
& MM,        & !IN, mixing model
& vLPole,    & !IN
& vXPole)      !IN
!--
!-- !!! check implementation of site models !!!
!--
  use M_Dtb_Const,only: R_jk
  use M_MixModel_Special
  !
  real(dp),        intent(in):: TdgK,Pbar
  type(T_MixModel),intent(in):: MM
  logical,         intent(in):: vLPole(:)
  real(dp),        intent(in):: vXPole(:)
  !
  real(dp),allocatable:: vXAtom(:)
  real(dp),allocatable:: vLIdeal(:),vLGam(:)
  real(dp):: G_IdMixRT,G_XSMixRT
  !
  ! GMixing= G_IdealMixing + G_XSMixing
  !   G_IdMix: ideal part of free energy of mixing of solution MM for a given P,T
  !   G_XSMix: excess free energy of mixing of MM MM for a given P,T
  !
  ! Gibbs free energy of the mixture at P,T is GMeca+GMix
  !
  G_XSMixRT=Zero
  !
  select case(MM%Model)
  !
  case(Mix_Molecular,Mix_Felspar) !("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR" or "FELSPAR"
    G_IdMixRT= - MixModel_Pole_SConf(MM,vLPole,vXPole) /R_jk
    !
    if(MM%NMarg>0) G_XSMixRT= &
    & MixModel_Pole_XsMargules(MM,vLPole,vXPole) /R_jk /TdgK
  !
  case(Mix_Site) ! ("SITE")
    allocate(vXAtom(MM%NAtom))
    call MixModel_XPoleToXSite(MM,vXPole,vXAtom)
    !
    G_IdMixRT= - MixModel_Site_SConf(MM,vXAtom) /R_jk
    !
    if(MM%NMarg>0) G_XSMixRT= &
    & MixModel_Site_XsMargules(MM,vXAtom) /R_jk /TdgK
    !
    deallocate(vXAtom)
  !
  case(Mix_Special) ! ("SPECIAL")
    allocate(vLIdeal(size(vXPole)))
    allocate(vLGam(size(vXPole)))
    !
    call MixModel_Special_Activities( &
    & MM%Name,TdgK,Pbar,vXPole,vLpole, &
    & vLIdeal,vLGam)
    !
    G_IdMixRT= sum(vXPole(:)*vLIdeal(:), MASK=vLPole(:))
    G_XSMixRT= sum(vXPole(:)*vLGam(:),   MASK=vLPole(:))
    !
    deallocate(vLIdeal)
    deallocate(vLGam)
  !
  end select
  !
  !!! GMeca is not computed here,
  !!! because we are concerned only with mixing properties
  !
  MixModel_GibbsMixRT= G_IdMixRT +G_XSMixRT
  !
end function MixModel_GibbsMixRT

end module M_T_MixModel


