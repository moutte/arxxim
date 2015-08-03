MODULE M_T_MixModel
!--
!-- implementation of solution models, and related potentials, activities, ...
!-- developed for minerals,
!-- to be used also for other mixtures (fluids, melts) that follow a mole fraction model
!--
  USE M_Kinds
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC:: T_MixModel,T_Margul
  PUBLIC:: MaxSite,MaxPole,MaxAtom,MaxMulti,MaxMarg
  !
  PUBLIC:: MixModel_Margul_Wg_Calc
  PUBLIC:: MixModel_Index
  PUBLIC:: MixModel_Zero
  PUBLIC:: MixModel_XPoleToXSite

  PUBLIC:: MixModel_Activities
  PUBLIC:: MixModel_GibbsMixRT
  PUBLIC:: MixModel_GibbsIdeal

  PUBLIC:: MixModel_Site_ActivIdeal

  !PUBLIC:: MixModel_Site_LnActivsIdeal
  !PUBLIC:: MixModel_Pole_LnActivsIdeal

  !PUBLIC:: MixModel_Site_SConf
  !PUBLIC:: MixModel_Pole_SConf

  !PUBLIC:: MixModel_Site_XsMargules
  !PUBLIC:: MixModel_Pole_XsMargules

  !PUBLIC:: MixModel_Site_LnGammaMargules
  !PUBLIC:: MixModel_Pole_LnGammaMargules

  !PUBLIC:: MixModel_Margules_Monome
  !
  !PUBLIC:: MixModel_SConf
  !PUBLIC:: MixModel_XsMargules
  !
  !REAL(dp),PARAMETER::& !from ModConst
  !& R_jK= 8.314510D0,&  !R constant, in J/K/Mole, RobieHemingway (from CODATA)
  !& Pref= 1.0000D0,&  !bar
  !& Tref= 298.150D0    !in degK,=25°C

  !----------------------------fixed dimensions for all T_MixModel's --
  !
  INTEGER,PARAMETER:: & !fixed dimensions for T_MixModel
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
  !! TYPE:: T_Margul
  !!   INTEGER:: NPoleMarg
  !!   !number (2:4) of "poles" (end-members) involved in the Margules term
  !!   !most commonly, NPole=2 (binary interaction), sometimes NPole=3
  !!   INTEGER,DIMENSION(1:4)::&
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
  !! ENDTYPE T_Margul
  
  TYPE:: T_Margul
    INTEGER,DIMENSION(1:MaxAtom):: vDegree
    ! new format,
    !   vDegree(:) contains all degrees,
    !     listed in the same order as
    !     the e-m list in case of molecular model,
    !     or as the site fractions in case of site model
    !   if a e-m or a site fraction is not involved, its degree is zero
    ! -> vDegree(:) contains information that was contained in vIPole and in vPower
    REAL(dp)::WG,WH,WS,WV,WCp
    INTEGER ::Kohler
    !caveat!! deCapitani allows the Kohler parameter to take fractional values
  ENDTYPE T_Margul
  !
  !---------------------------------------------------------------------
  !
  !WG=   WH &
  !&   + WCP*(T-T0) &
  !&   -(WS +WCP*LOG(T/T0))*T &
  !&   + WV *P
  !WCp commonly Zero -> WG= WH -T*WS +P.WV
  !
  !---------------------------------- about using Margules parameters --
  != one Margules term M of the excess G calculated as
  !  Y=M%WG
  !<old version>
  !  DO I=1,M%NPoleMarg !NPole= number of poles involved in term M
  !    Y=Y*X(M%iPol(I))**(M%Power(I))
  !  !! IPole=list of indices (in the solution) of the poles involved in M
  !  ENDDO
  !  example:
  !    term 112 -> Power(1)=2, Power(2)=1 -> W112 *X(IPol(1))**2 *X(Ipol(2))
  !
  !<new version>
  !  DO I=1,S%NPole != number of poles in the mixture
  !    IF(M%vDegree(I)>0) Y=Y*X(I)**M%vDegree(I))
  !    !! no more use of vIPole(:)
  !  ENDDO
  !
  !--------------------------------------------------------------------
  !
  !-----------------------------------------------------------T_MixModel 
  TYPE:: T_MixModel

    CHARACTER(LEN=23):: Name    !same length as name in tMinThr
    CHARACTER(LEN=3) :: Typ     !"MIN","GAS",...
    CHARACTER(LEN=15):: Model   !TYPE of EoS: IDEAL, MOLECULAR, SITE, FELSPAR, GUGGENHEIM, ...

    INTEGER:: NPole             !Nr of EndMembers, must be <=MaxPole
    INTEGER:: NSite             !Nr of Sites, must be <=MaxSite
    INTEGER:: NMarg             !Nr of Margules terms

    ! LOGICAL:: vHasPole(1:MaxPole)

    CHARACTER(LEN=23):: vNamPole(1:MaxPole)
    ! names of end-members
    ! (same name must be found also in pure species database)

    INTEGER:: vIPole(1:MaxPole)
    ! vIPole(i)= index of end member i in current species list
    !   species list is generally vSpc in M_Global_Vars,
    !   thus,
    !     vIPole(i)   is Species_Index(vSpc, vNamPole(:)
    !     vNamPole(i) is vSpc(vIPole(i))%Name

    INTEGER:: vMulti(1:MaxSite)
    ! site multiplicity
    ! for molecular models, multiplicity is vMulti(1)
    ! vMulti(i)<=MaxMulti forall i
    !
    !! REAL(dp):: GuggA0,GuggA1
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
    TYPE(T_Margul):: vMarg(1:MaxMarg)
    ! array of Margules parameters
    !
    !! vIMargSite still useful, although information is also contained in vPower
    INTEGER:: vIMargSite(1:MaxMarg)
    ! for molecular (or one-site models), vIMargSite(1:NMarg)= 1
    ! for multisite models, vIMargSite(i) tells which site vMarg(i) is for (??)

    !----------------------------------------------------for site models 
    INTEGER:: NAtom
    ! dimension of composition vector, <=MaxAtom
    !
    REAL(dp):: vPoleCoeff(1:MaxPole)
    ! normalization factors for site models
    !
    CHARACTER(LEN=3)::vNamSite(1:MaxSite)
    ! names of sites, e.g. M1_, M2_, M23, A__, OH_, ...
    !
    CHARACTER(LEN=6):: vNamAtom(1:MaxAtom)
    ! 6 chars= 3 for element name, 3 for site name
    ! names of site fractions, e.g. Mg_M1_, Mg_M2_, Si_T__, OH_OH_, ...
    !
    INTEGER:: vAtomMulti(1:MaxAtom)
    ! gives the multiplicity of the site containing the atom
    INTEGER:: vIAtomSite(1:MaxAtom)
    ! gives the index of the site containing the atom
    INTEGER:: tPoleAtom(1:MaxPole,1:MaxAtom)
    ! table of end-members w/r site fractions
    ! example
    ! SITE
    !   A    1  K__Na_     !iSite=1; K:iEle=1, Na:iEle=2
    !   M2A  1  Al_Mg_Fe_
    !   T1   2  Al_Si_
    ! END
    ! >> composition vector will contain
    !   (xK_A, xNa_A, xAl_M2A, xMg_M2A, xFe_M2A, xAl_T1, xSi_T1)
    ! list of poles and site occupancy table:
    ! POLE
    !   !    A(1) M2A(1) T1(2)  ! M2B(1) T2(2)
    !   Mu   K    Al     AlSi   ! Al     SiSi
    !   Pa   Na   Al     AlSi   ! Al     SiSi
    !   MCel K    Mg     SiSi   ! Al     SiSi
    !   FCel K    Fe     SiSi   ! Al     SiSi
    ! END
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
  END TYPE T_MixModel
  !----------------------------------------------------------/T_MixModel
  !
  !!!-----------------------------------------------------------examples
  !!!
  !!!  SITE
  !!!    M2 2 MG_FE_
  !!!    M1 1 MG_FE_AL_ !may have up to MaxElSite elements
  !!!    T2 2 AL_SI_
  !!!  ENDSITE
  !!!  POLE
  !!!    !__________       M2____    M1_        T_____
  !!!    !possible=        MG_FE_    MG_FE_AL_  AL_SI_
  !!!    !vMulti=          2         1          2
  !!!    !phlogopite       Mg2       Mg         AlSi3
  !!!    PHLOGOPITE        MG_MG_    MG_        AL_SI_
  !!!    ANNITE            FE_FE_    FE_        AL_SI_
  !!!    SIDEROPHYLLITE    FE_FE_    AL_        AL_AL_
  !!!    EASTONITE         MG_MG_    AL_        AL_AL_
  !!!  ENDPOLE
  !!!  .../...
  !!!ENDSOLUTION
  !!!
  !!!
  !!!
  !!!-> example of sol.composition:
  !!! COMPOSITION
  !!!   M2 1.00 1.00 !sum=2.00, =vMulti(1)
  !!!   M1 0.50 0.50 !sum=1.00, =vMulti(2)
  !!!   T2 1.00 1.00 !sum=2.00, =vMulti(3)
  !!! END
  !!!
  !!!----------------------------------------------------------/examples
  !

CONTAINS

SUBROUTINE MixModel_Zero(S)

  TYPE(T_MixModel),INTENT(OUT)::S
  
  S%Name=     "XXX"
  S%Model=    "IDEAL"
  S%NPole=    0
  S%NSite=    1
  S%vMulti(1)=1
  S%nMarg=    0
  ! S%vHasPole(:)= .FALSE.
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
  
END SUBROUTINE MixModel_Zero

INTEGER FUNCTION MixModel_Index(V,Str)
!-- position of solution model named Str in V(1:nMixModel)
!-----------------------------------------------------------------------
  TYPE(T_MixModel),DIMENSION(:),INTENT(IN):: V
  CHARACTER(*),                 INTENT(IN):: Str
  INTEGER     ::I
  MixModel_Index=0
  !---------------------------------------------------------------------
  IF(SIZE(V)==0) RETURN
  !
  I=0
  DO
    I=I+1
    IF(TRIM(Str)==TRIM(V(I)%Name)) THEN
      MixModel_Index=I   ;    EXIT
    ENDIF
    IF(I==SIZE(V)) EXIT
  ENDDO
  !
  RETURN
END FUNCTION MixModel_Index

SUBROUTINE MixModel_XPoleToXSite(MM,vX,vY)

  TYPE(T_MixModel),INTENT(IN) :: MM
  REAL(dp),        INTENT(IN) :: vX(:) ! end member mole fractions
  REAL(dp),        INTENT(OUT):: vY(:) ! site fractions
  !---------------------------------------------------------------------
  INTEGER:: iP,iEl,J
  !
  vY(:)= Zero
  DO iP=1,MM%NPole
    DO iEl=1,MM%NAtom
      J= MM%tPoleAtom(iP,iEl)
      !PRINT *,J
      IF(J /= 0) vY(iEl)= vY(iEl) + vX(iP) *REAL(J) /REAL(MM%vAtomMulti(iEl))
      !PRINT *,vY
    END DO
    !PAUSE
  END DO
  !
  RETURN
END SUBROUTINE MixModel_XPoleToXSite

REAL(dp) FUNCTION MixModel_Site_ActivIdeal( &
& MM,   & !IN= Mix.model, generally =vMixModel(Fas%iModel)
& IPol, & !IN= index of end-member in models
& vXatom) !IN= phase composition
!--
!-- for site mixing models
!-- returns ideal activity of end-member IPol,
!-- following solution model MM
!--
  TYPE(T_MixModel),INTENT(IN):: MM
  REAL(dp),        INTENT(IN):: vXAtom(:)
  INTEGER,         INTENT(IN):: IPol
  !---------------------------------------------------------------------
  REAL(dp):: Y
  INTEGER :: iEl
  !---------------------------------------------------------------------
  Y= One
  !
  DO iEl= 1,MM%NAtom
    !IF(MM%tPoleAtom(iPol,iEl)/=0) Y= Y *vXAtom(iEl)**MM%tPoleAtom(iPol,iEl)
    IF(MM%tPoleAtom(iPol,iEl)/=0) Y= Y *vXAtom(iEl)**MM%vAtomMulti(iEl)
  ENDDO
  !
  MixModel_Site_ActivIdeal= Y *MM%vPoleCoeff(iPol)
  !
  !IF(iDebug>1) WRITE(fTrc,'(I3,1X,G15.6,1X,A)') &
  !& iPol,MixModel_Site_ActivIdeal,TRIM(MM%vNamPole(iPol))
  !
ENDFUNCTION MixModel_Site_ActivIdeal

SUBROUTINE MixModel_Pole_LnActivsIdeal( & !
& MM,     &  !IN= Mix.model
& vX,     &  !IN= phase composition
& vLpole, &  !IN
& vLnAct)    !OUT
!--
!-- for end-member mixing models
!-- returns log(ideal activity) of all 'active' end-members
!--
  TYPE(T_MixModel),INTENT(IN) :: MM
  REAL(dp),        INTENT(IN) :: vX(:)
  LOGICAL,         INTENT(IN) :: vLpole(:)
  REAL(dp),        INTENT(OUT):: vLnAct(:)
  !---------------------------------------------------------------------
  INTEGER:: iP
  !---------------------------------------------------------------------
  vLnAct(:)= Zero
  !
  SELECT CASE(TRIM(MM%Model))
  !
  CASE("IDEAL","POLE","MOLECULAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR"
    DO iP=1,MM%NPole
      IF(vLPole(iP)) &
      & vLnAct(iP)= vLnAct(iP) + LOG(vX(iP))*MM%vMulti(1)
    ENDDO
  !
  CASE("FELSPAR")
    IF(vLPole(1)) vLnAct(1)= LOG(vX(1)) !+LOG(1-vX(3)*vX(3))
    IF(vLPole(2)) vLnAct(2)= LOG(vX(2)) !+LOG(1-vX(3)*vX(3))
    IF(vLPole(3)) THEN
      vLnAct(3)= LOG(vX(3))+2*LOG((One+vX(3))/2.0D0)
      IF(vLPole(1)) vLnAct(1)= LOG(vX(1)) +LOG(One-vX(3)**2)
      IF(vLPole(2)) vLnAct(2)= LOG(vX(2)) +LOG(One-vX(3)**2)
    ENDIF
  !
  END SELECT
  !
END SUBROUTINE MixModel_Pole_LnActivsIdeal

SUBROUTINE MixModel_Site_LnActivsIdeal( & !
& MM,     &  !IN= Mix.model
& vXatom, &  !IN= phase composition
& vLpole, &  !IN
& Ok,     &  !OUT
& vLnAct)    !OUT
!--
!-- for site mixing models
!-- returns log(ideal activity) of all 'active' end-members
!--
  !---------------------------------------------------------------------
  TYPE(T_MixModel),INTENT(IN) :: MM
  REAL(dp),        INTENT(IN) :: vXAtom(:)
  LOGICAL,         INTENT(IN) :: vLpole(:)
  LOGICAL,         INTENT(OUT):: Ok
  REAL(dp),        INTENT(OUT):: vLnAct(:)
  !---------------------------------------------------------------------
  REAL(dp):: Y
  INTEGER :: iP,iEl,J
  !---------------------------------------------------------------------
  Ok= .TRUE.
  !
  DO iP= 1,MM%NPole

    IF(vLPole(iP)) THEN

      Y= Zero
      DO iEl= 1,MM%NAtom
        J= MM%tPoleAtom(iP,iEl)
        IF(J/=0) THEN
          IF(vXAtom(iEl)>Zero) THEN
            Y= Y + LOG(vXAtom(iEl)) *MM%vAtomMulti(iEl) !*J
          ELSE
          !! normally, should never happen, if vLPole is .true. ...
            Ok= .FALSE.
          ENDIF
        ENDIF
      ENDDO
      
      ! add correction term, to have activ=1 for pure end member
      vLnAct(iP)= Y + LOG(MM%vPoleCoeff(iP))

    ENDIF

  ENDDO
  !
END SUBROUTINE MixModel_Site_LnActivsIdeal

REAL(dp) FUNCTION MixModel_Site_SConf( & !
& MM,   & !IN, mixing model
& vXatom)    !IN, composition
!--
!------------------------------------ compute configurational entropy --
!---------------------------------------------------- for site mixing --
!--
  USE M_Dtb_Const,ONLY: R_jk
  !
  TYPE(T_MixModel),INTENT(IN):: MM
  REAL(dp),        INTENT(IN):: vXatom(:)
  !
  REAL(dp):: Sconf,X
  INTEGER :: iEl
  !
  SConf= Zero
  !
  ! sum of x*ln(x) on the different sites,
  ! taking account also of site multiplicity
  DO iEl=1,MM%NAtom
    X= vXAtom(iEl)
    IF(X > Zero) &
    & Sconf= Sconf + X *LOG(X) *MM%vAtomMulti(iEl)
  ENDDO
  !
  MixModel_Site_Sconf= -R_jk*Sconf
  !
ENDFUNCTION MixModel_Site_SConf

REAL(dp) FUNCTION MixModel_Pole_SConf( & !
& MM,   & !IN, mixing model
& vLpole,  & !IN, pole is present / absent
& vX)        !IN, composition
!--
!-- compute configurational entropy --
!-- for end member mixing --
!--
  USE M_Dtb_Const,ONLY: R_jk
  !
  TYPE(T_MixModel),INTENT(IN):: MM
  LOGICAL,         INTENT(IN):: vLpole(:)
  REAL(dp),        INTENT(IN):: vX(:)
  !
  REAL(dp):: Sconf
  INTEGER :: iP
  !
  SConf= Zero
  !
  SELECT CASE(TRIM(MM%Model))

  CASE("IDEAL","POLE","MOLECULAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR"
    DO iP=1,MM%NPole
      IF(vLPole(iP)) Sconf= Sconf + vX(iP)*LOG(vX(iP))
    ENDDO
    Sconf= SConf *MM%vMulti(1)
    !
  CASE("FELSPAR")
  !! "FELSPAR" is actually a special case of "MOLECULAR" mixing
    ! ALBITE
    IF(vLPole(1)) Sconf= Sconf + vX(1)*LOG(vX(1))
    ! K-FELSPAR
    IF(vLPole(2)) Sconf= Sconf + vX(2)*LOG(vX(2))
    ! ANORTHITE
    IF(vLPole(3)) THEN
      Sconf= Sconf + vX(3) *( LOG(vX(3)) +LOG((One +vX(3))/Two) *Two )
      IF(vLPole(1)) Sconf= Sconf + vX(1)*LOG(One -vX(3)*vX(3))
      IF(vLPole(2)) Sconf= Sconf + vX(2)*LOG(One -vX(3)*vX(3))
    ENDIF
    !
  END SELECT
  !
  MixModel_Pole_Sconf= -R_jk*Sconf
  !
ENDFUNCTION MixModel_Pole_SConf

REAL(dp) FUNCTION MixModel_SConf( & !
& MM,   & !IN, mixing model
& vLPole,  & !IN
& vX)        !IN
!--
!-- compute configurational entropy --
!--
  TYPE(T_MixModel),INTENT(IN):: MM
  LOGICAL,         INTENT(IN):: vLPole(:)
  REAL(dp),        INTENT(IN):: vX(:)
  !
  REAL(dp),ALLOCATABLE:: vXAtom(:)
  !
  SELECT CASE(TRIM(MM%Model))

  CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR"
    MixModel_SConf= MixModel_Pole_SConf(MM,vLPole,vX) ! vX is Fas%vXPole

  CASE("SITE")
    ALLOCATE(vXAtom(MM%NAtom))
    vXAtom(:)= vX(:)
    MixModel_SConf= MixModel_Site_SConf(MM,vXAtom) ! vX is Fas%vXAtom
    DEALLOCATE(vXAtom)

  END SELECT
  !
ENDFUNCTION MixModel_SConf

REAL(dp) FUNCTION MixModel_Margules_Monome(M,vX)
  TYPE(T_Margul),INTENT(IN):: M
  REAL(dp),      INTENT(IN):: vX(:)
  !
  REAL(dp):: Y
  INTEGER :: i
  !
  Y= M%WG
  DO i=1,SIZE(vX)
    IF(M%vDegree(i)>0) &
    & Y= Y *vX(i)**M%vDegree(i)
  ENDDO
  !
  MixModel_Margules_Monome= Y
  !
END FUNCTION MixModel_Margules_Monome

REAL(dp) FUNCTION MixModel_Pole_LnGammaMargules( & !
& MM,      & !IN, mixing model
& iP,      & !IN
& vMonome, & !
& vXpole)    !IN, composition
!--
!-- for MOLECULAR model, compute RTln(Gamma) related to Margules Terms
!-- (according to expression by deCapitani-Kirschen)
!--
  TYPE(T_MixModel),INTENT(IN):: MM
  INTEGER,         INTENT(IN):: iP
  REAL(dp),        INTENT(IN):: vMonome(:),vXpole(:)
  !
  TYPE(T_Margul)::M
  REAL(dp):: LnGam,S_i,Y,Z
  INTEGER :: iM,I
  INTEGER :: dSi_dXj
  !
  ! RT.LnGam(j)=      &
  !   SUM(i=1..NMarg) &
  !   W_i1..ip . X_i1..X_ip (S_i)^(-k_i) &
  !   ( (1 - p_i + k_i) + q_j / X_j - (k_i/S_i).dSi_dXj ) )
  ! P_i = degree of monomial X_i1... = SUM(M%vDegree(:))
  ! S_i = sum of mole fractions involved in term = SUM(vX(M%iPol(1:M%NPole)))
  ! Q_j = number of indices, among i1..ip, equal to j,
  !     = M%vDegree(j)
  !
  LnGam=Zero
  !
  DO iM=1,MM%NMarg

    M= MM%vMarg(iM)
    Y= vMonome(iM)
    Z= M%vDegree(iP) /vXpole(iP) +One -SUM(M%vDegree(:))
    !
    IF(M%Kohler/=0) THEN
      S_i= Zero
      DO I=1,MM%NPole
        IF(M%vDegree(I)>0) S_i= S_i +vXpole(I)
      END DO
      dSi_dXj=0
      IF(M%vDegree(iP)>0) dSi_dXj= 1
      !
      Y= Y / S_i**M%Kohler
      Z= Z + M%Kohler
      IF(dSi_dXj/=0) Z= Z - M%Kohler *dSi_dXj/S_i
    ENDIF
    !
    LnGam= LnGam + Y*Z
  ENDDO
  !
  MixModel_Pole_LnGammaMargules= LnGam ! gives actually R*T*Ln(Gamma) !!
  !
END FUNCTION MixModel_Pole_LnGammaMargules

REAL(dp) FUNCTION MixModel_Site_LnGammaMargules( & !
& MM,      & !IN, mixing model
& iP,      & !IN
& vMonome, & !
& vXatom)    !IN, composition
!--
!-- for SITE model, compute RTln(Gamma) related to Margules Terms
!-- ( should be called only when vLpole(iP) is true )
!--
  TYPE(T_MixModel),INTENT(IN):: MM
  INTEGER,         INTENT(IN):: iP
  REAL(dp),        INTENT(IN):: vXatom(:),vMonome(:)
  !
  TYPE(T_Margul)::M
  REAL(dp):: X != RT.ln(gamma)
  INTEGER :: iM,J,iEl,iS
  !
  ! n * RT.LnGam_m=   &
  !   SUM(i=1..NMarg) &
  !   [ W_i1..ip . X_i1..X_ip  * ( Q_m /X_m - 2 ) ] !for degree 3 Margules !!
  ! n   = number of sites on which mixing occurs
  ! q_m = number of indices, among i1..ip, equal to m,
  !     = M%vDegree(m)
  !
  !! LnAct= Zero
  X= Zero
  !
  DO iEl= 1,MM%NAtom
    J= MM%tPoleAtom(iP,iEl)
    IF(J/=0) THEN
      !! LnAct= LnAct + LOG(vXAtom(iEl)) *J
      DO iM=1,MM%NMarg
        iS= MM%vIMargSite(iM)
        IF( iS == MM%vIAtomSite(iEl) ) THEN
          M= MM%vMarg(iM)
          X=   X &
          &  + vMonome(iM) &
          &    *(M%vDegree(iEl) /vXatom(iEl) +One -SUM(M%vDegree(:))) &
          &    /REAL(MM%vAtomMulti(iEl))
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !
  MixModel_Site_LnGammaMargules= X ! is R*T*ln(Gamma) !!
  !
END FUNCTION MixModel_Site_LnGammaMargules

REAL(dp) FUNCTION MixModel_Pole_XsMargules( & !
& MM,     & !IN, mixing model
& vLPole, & !IN
& vXpole)   !IN, composition
!--
!-- for MOLECULAR model
!--
  TYPE(T_MixModel),INTENT(IN):: MM
  LOGICAL,         INTENT(IN):: vLPole(:)
  REAL(dp),        INTENT(IN):: vXpole(:)
  !
  TYPE(T_Margul):: M
  REAL(dp):: Y,G_XSMix
  INTEGER :: iM,iP
  !
  G_XSMix= Zero
  !
  DO iM=1,MM%NMarg ! number of Margules terms
    !
    M= MM%vMarg(iM)  ! the Margules structure
    Y= M%WG             ! the interaction energy for this Margules term
    !
    DO iP=1,MM%NPole
      IF(M%vDegree(iP)>0) THEN != pole iP involved in monome Y
        IF(vLPole(iP)) THEN
          Y= Y *vXPole(iP)**M%vDegree(iP)
        ELSE
          Y= Zero
        ENDIF
      ENDIF
    ENDDO
    !
    G_XSMix= G_XSMix +Y
  ENDDO
  !
  MixModel_Pole_XsMargules= G_XSMix
  !
  RETURN
END FUNCTION MixModel_Pole_XsMargules

REAL(dp) FUNCTION MixModel_Site_XsMargules( & !
& MM,    & !IN, mixing model
& vXatom)  !IN
!--
!-- for SITE model
!--
  TYPE(T_MixModel),INTENT(IN):: MM
  REAL(dp),        INTENT(IN):: vXAtom(:)
  !
  TYPE(T_Margul):: M
  REAL(dp):: Y,G_XSMix
  INTEGER :: iM,iEl
  !
  G_XSMix= Zero
  !
  DO iM=1,MM%NMarg ! number of Margules terms
    !
    M= MM%vMarg(iM)  ! the Margules structure
    Y= M%WG          ! the interaction energy for this Margules term
    !
    DO iEl=1,MM%NAtom
      IF(M%vDegree(iEl)>0) & != site fraction iEl is involved in monome Y
      & Y= Y *vXAtom(iEl)**M%vDegree(iEl)
    ENDDO
    !
    G_XSMix= G_XSMix +Y
  ENDDO
  !
  MixModel_Site_XsMargules= G_XSMix
  !
  RETURN
END FUNCTION MixModel_Site_XsMargules

REAL(dp) FUNCTION MixModel_XsMargules( & !
& MM,      & !IN, mixing model
& vLPole,  & !IN
& vX)        !IN
!--
  TYPE(T_MixModel),INTENT(IN):: MM
  LOGICAL,         INTENT(IN):: vLPole(:)
  REAL(dp),        INTENT(IN):: vX(:)
  !
  MixModel_XsMargules= Zero
  !
  IF(MM%NMarg==0) RETURN
  !
  SELECT CASE(TRIM(MM%Model))
  !
  CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR" or "FELSPAR"
    MixModel_XsMargules= &
    & MixModel_Pole_XsMargules(MM,vLPole,vX) ! vX is Fas%vXPole
  !
  CASE("SITE")
    MixModel_XsMargules= &
    & MixModel_Site_XsMargules(MM,vX) ! vX is Fas%vXAtom
  !
  END SELECT
  !
  RETURN
END FUNCTION MixModel_XsMargules

SUBROUTINE MixModel_Margul_Wg_Calc(TdgK,Pbar,S)
!--
!-- from S%vMarg(:)%WH,WS,..., calc. values of S%vMarg(:)%WG at given P,T
!--
  USE M_Dtb_Const,ONLY: Tref
  !
  REAL(dp),        INTENT(IN)   :: TdgK,Pbar
  TYPE(T_MixModel),INTENT(INOUT):: S
  !
  TYPE(T_Margul) ::M
  INTEGER        ::iM
  !
  IF(S%NMarg>0) THEN

    DO iM=1,S%NMarg
      M=S%vMarg(iM)
      S%vMarg(iM)%WG= M%WH + M%WCP*(TdgK-Tref) &
      &             -(M%WS +M%WCP*LOG(TdgK/Tref))*TdgK &
      &             + M%WV*Pbar
    ENDDO

  ENDIF
  !
  RETURN
END SUBROUTINE MixModel_Margul_Wg_Calc

SUBROUTINE MixModel_Activities( & !
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
  USE M_Dtb_Const,ONLY: R_jk
  USE M_Trace
  USE M_MixModel_Special
  !
  REAL(dp),        INTENT(IN) :: Pbar, TdgK
  TYPE(T_MixModel),INTENT(IN) :: MM
  LOGICAL,         INTENT(IN) :: vLPole(:)
  REAL(dp),        INTENT(IN) :: vXPole(:)
  LOGICAL,         INTENT(OUT):: Ok
  CHARACTER(*),    INTENT(OUT):: Msg
  REAL(dp),        INTENT(OUT):: vLGam(:)
  REAL(dp),        INTENT(OUT):: vLIdeal(:)
  REAL(dp),        INTENT(OUT):: vLnAct(:)
  !
  INTEGER :: iP,iM
  REAL(dp):: P
  REAL(dp),ALLOCATABLE:: vMonome(:)
  REAL(dp),ALLOCATABLE:: vXAtom(:)
  !
  P=Pbar !for future use ??
  !
  Ok= .TRUE.
  Msg= "Ok"
  !
  !F%vLPole(1:MM%NPole)= vXPole(1:MM%NPole)>Zero
  !
  vLGam(:)=   Zero !default
  vLIdeal(:)= Zero !default
  !
  SELECT CASE(TRIM(MM%Model))
  !
  CASE("SPECIAL")
    !
    CALL MixModel_Special_Activities( &
    & MM%Name,TdgK,Pbar,vXPole,vLpole, &
    & vLIdeal,vLGam)
    !
  CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR" or "FELSPAR"
    !
    CALL MixModel_Pole_LnActivsIdeal(MM,vXPole,vLpole,vLIdeal)
    !
    !------------------------------------------------- Margules Terms --
    IF(MM%NMarg>0) THEN
      ALLOCATE(vMonome(MM%NMarg))
      !
      DO iM=1,MM%NMarg
        vMonome(iM)= MixModel_Margules_Monome(MM%vMarg(iM),vXPole)
      ENDDO
      DO iP=1,MM%NPole
        IF(vLPole(iP)) &
        & vLGam(iP)= MixModel_Pole_LnGammaMargules(MM,iP,vMonome,vXPole) &
        &            /R_jk/TdgK
      ENDDO
      !
      DEALLOCATE(vMonome)
      !
    ENDIF
    !-------------------------------------------------/Margules Terms --
    !
  CASE("SITE")
    !
    ALLOCATE(vXAtom(MM%NAtom))
    CALL MixModel_XPoleToXSite(MM,vXPole,vXAtom)
    !
    CALL MixModel_Site_LnActivsIdeal(MM,vXAtom,vLpole,Ok,vLIdeal)
    !
    !------------------------------------------------- Margules Terms --
    IF(MM%NMarg>0) THEN
      ALLOCATE(vMonome(MM%NMarg))
      !
      DO iM=1,MM%NMarg
        vMonome(iM)= MixModel_Margules_Monome(MM%vMarg(iM),vXAtom)
      ENDDO
      DO iP=1,MM%NPole
        IF(vLPole(iP)) &
        & vLGam(iP)= MixModel_Site_LnGammaMargules(MM,iP,vMonome,vXAtom) &
        &            /R_jk/TdgK
      ENDDO
      !
      DEALLOCATE(vMonome)
      !
    ENDIF
    !-------------------------------------------------/Margules Terms --
    DEALLOCATE(vXAtom)
    !
  CASE DEFAULT
    Ok= .FALSE.
    Msg= TRIM(MM%Model)//"= invalid MM%Model in MixModel_CalcActivities"
    !
  END SELECT
  !
  vLnAct(1:MM%NPole)= vLIdeal(1:MM%NPole) + vLGam(1:MM%NPole)
  !
  IF(iDebug>3) THEN !=========================================< trace ==
    WRITE(fTrc,'(A)') "MixModel_CalcActivs -> X,ActIdeal,Gamma,Activ"
    WRITE(fTrc,'(2A)') "MixModel=",MM%Name
    DO iP=1,MM%NPole
      IF(vLPole(iP)) THEN
        WRITE(fTrc,'()')
        WRITE(fTrc,'(A,I2,A1,A15,A1,4(A4,G15.6,A1))') &
        & "POLE",iP,              T_,&
        & TRIM(MM%vNamPole(iP)),  T_,&
        & "Frc=",vXPole(iP),          T_,&
        & "XId=",EXP(vLIdeal(iP)),T_,&
        & "Gam=",EXP(vLGam(iP)),  T_,&
        & "Act=",EXP(vLnAct(iP)), T_
      ENDIF
    ENDDO
  ENDIF !====================================================</ trace ==
  !
END SUBROUTINE MixModel_Activities

REAL(dp) FUNCTION MixModel_GibbsIdeal( & !
& TdgK,Pbar, & !IN
& MM,        & !IN, mixing model
& vLPole,    & !IN
& vXPole)      !IN
!--
!-- !!! todo - check implementation of site models !!!
!--
  USE M_Dtb_Const,ONLY: R_jk
  USE M_MixModel_Special
  !
  REAL(dp),        INTENT(IN):: TdgK,Pbar
  TYPE(T_MixModel),INTENT(IN):: MM
  LOGICAL,         INTENT(IN):: vLPole(:)
  REAL(dp),        INTENT(IN):: vXPole(:)
  !
  REAL(dp),ALLOCATABLE:: vXAtom(:)
  REAL(dp),ALLOCATABLE:: vLIdeal(:),vLGam(:)
  REAL(dp):: G_IdMix
  !
  ! G_IdMix: ideal part of free energy of mixing of solution MM for a given P,T
  !
  SELECT CASE(TRIM(MM%Model))
  !
  CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommend using "MOLECULAR" or "FELSPAR" only
    G_IdMix= -TdgK *MixModel_Pole_SConf(MM,vLPole,vXPole)
    !
  CASE("SITE")
    ALLOCATE(vXAtom(MM%NAtom))
    CALL MixModel_XPoleToXSite(MM,vXPole,vXAtom)
    !
    G_IdMix= -TdgK *MixModel_Site_SConf(MM,vXAtom)
    !
    DEALLOCATE(vXAtom)
  !
  CASE("SPECIAL")
    ALLOCATE(vLIdeal(SIZE(vXPole)))
    ALLOCATE(vLGam(SIZE(vXPole)))
    !
    CALL MixModel_Special_Activities( &
    & MM%Name,TdgK,Pbar,vXPole,vLpole, &
    & vLIdeal,vLGam)
    !
    G_IdMix= SUM(vXPole(:)*vLIdeal(:), MASK=vLPole(:)) *R_jk  *TdgK
    !
    DEALLOCATE(vLIdeal)
    DEALLOCATE(vLGam)
  !
  END SELECT
  !
  MixModel_GibbsIdeal= G_IdMix
  !
ENDFUNCTION MixModel_GibbsIdeal

REAL(dp) FUNCTION MixModel_GibbsMixRT( & !
& TdgK,Pbar, & !IN
& MM,        & !IN, mixing model
& vLPole,    & !IN
& vXPole)      !IN
!--
!-- !!! check implementation of site models !!!
!--
  USE M_Dtb_Const,ONLY: R_jk
  USE M_MixModel_Special
  !
  REAL(dp),        INTENT(IN):: TdgK,Pbar
  TYPE(T_MixModel),INTENT(IN):: MM
  LOGICAL,         INTENT(IN):: vLPole(:)
  REAL(dp),        INTENT(IN):: vXPole(:)
  !
  REAL(dp),ALLOCATABLE:: vXAtom(:)
  REAL(dp),ALLOCATABLE:: vLIdeal(:),vLGam(:)
  REAL(dp):: G_IdMixRT,G_XSMixRT
  !
  ! GMixing= G_IdealMixing + G_XSMixing
  !   G_IdMix: ideal part of free energy of mixing of solution MM for a given P,T
  !   G_XSMix: excess free energy of mixing of MM MM for a given P,T
  !
  ! Gibbs free energy of the mixture at P,T is GMeca+GMix
  !
  G_XSMixRT=Zero
  !
  SELECT CASE(TRIM(MM%Model))
  !
  CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR" or "FELSPAR"
    G_IdMixRT= - MixModel_Pole_SConf(MM,vLPole,vXPole) /R_jk
    !
    IF(MM%NMarg>0) G_XSMixRT= &
    & MixModel_Pole_XsMargules(MM,vLPole,vXPole) /R_jk /TdgK
  !
  CASE("SITE")
    ALLOCATE(vXAtom(MM%NAtom))
    CALL MixModel_XPoleToXSite(MM,vXPole,vXAtom)
    !
    G_IdMixRT= - MixModel_Site_SConf(MM,vXAtom) /R_jk
    !
    IF(MM%NMarg>0) G_XSMixRT= &
    & MixModel_Site_XsMargules(MM,vXAtom) /R_jk /TdgK
    !
    DEALLOCATE(vXAtom)
  !
  CASE("SPECIAL")
    ALLOCATE(vLIdeal(SIZE(vXPole)))
    ALLOCATE(vLGam(SIZE(vXPole)))
    !
    CALL MixModel_Special_Activities( &
    & MM%Name,TdgK,Pbar,vXPole,vLpole, &
    & vLIdeal,vLGam)
    !
    G_IdMixRT= SUM(vXPole(:)*vLIdeal(:), MASK=vLPole(:))
    G_XSMixRT= SUM(vXPole(:)*vLGam(:),   MASK=vLPole(:))
    !
    DEALLOCATE(vLIdeal)
    DEALLOCATE(vLGam)
  !
  END SELECT
  !
  !!! GMeca is not computed here,
  !!! because we are concerned only with mixing properties
  !
  MixModel_GibbsMixRT= G_IdMixRT +G_XSMixRT
  !
ENDFUNCTION MixModel_GibbsMixRT

ENDMODULE M_T_MixModel


