MODULE M_T_MixPhase
  USE M_Kinds
  USE M_Trace,ONLY: fTrc, iDebug, T_, Stop_
  USE M_T_MixModel
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: T_MixPhase
  !
  PUBLIC:: MixPhase_Index
  PUBLIC:: MixPhase_Zero
  PUBLIC:: MixPhase_CalcActivs
  PUBLIC:: MixPhase_NormCompo
  PUBLIC:: MixPhase_Weit
  PUBLIC:: MixPhase_GibbsRT
  PUBLIC:: MixPhase_Constraint
  PUBLIC:: MixPhase_Volume
  !
  PUBLIC:: MixPhase_CalcMixing   !used in Dtb_Test_Solution

  TYPE:: T_MixPhase
  ! implements a mixture phase (= symetric model)(= )
  ! includes a pointer to the mixing model,
  ! the composition data, and stores activity data
    !----------------------------------------------------- fixed data --
    CHARACTER(LEN=23):: Name
    INTEGER          :: iModel !index of mixing model in vMixModel
    !-----------------------------------------------------/fixed data --
    !
    !-------------------------------------------------- variable data --
    !----------------------------------------------- composition data --
    LOGICAL :: vLPole(1:MaxPole)
    ! vLPole(I)= T/F == end member I is present /absent
    REAL(dp):: vXPole(1:MaxPole)
    ! comp'n in mole fractions of the end members
    !REAL(dp):: vXAtom(1:MaxAtom) !comp'n in atom fractions, size=Model%NAtom
    !
    !----------------- values at current (T,P,X), for the end-members --
    REAL(dp),DIMENSION(1:MaxPole):: & !
    & vLnAct,  & !Ln(Activities of end-members 1:Model%nPole)
    & vLIdeal, & !Ln(ideal Activities of end-members 1:Model%nPole)
    & vLGam      !Ln(activ.coeff's)
    !! & vLMarg     !Ln(activ.coeff's related to Margules parameters)
    !
    !------------------------- values at current (T,P), for the phase --
    REAL(dp):: Grt,H,S,V,Cp
    !--------------------------------------------------/variable data --
    !
  ENDTYPE T_MixPhase

  TYPE:: T_MixPhaseData
  !-- variable data characteristic of a mixture phase (- symetric model)
  !-- includes composition data, activity data, etc.
    !----------------------------------------------- composition data --
    LOGICAL :: vLPole(1:MaxPole) ! end member is present/absent
    REAL(dp):: vXPole(1:MaxPole) ! comp'n in mole fractions of the end members
    !REAL(dp):: vXAtom(1:MaxAtom) ! comp'n in atom fractions, size=Model%NAtom
    !
    !----------------- values at current (T,P,X), for the end-members --
    REAL(dp),DIMENSION(1:MaxPole):: & !
    & vLnAct,  & ! Ln(Activities of end-members 1:Model%nPole)
    & vLIdeal, & ! Ln(ideal Activities of end-members 1:Model%nPole)
    & vLGam,   & ! Ln(activ.coeff's)
    & vLMarg     ! Ln(activ.coeff's related to Margules parameters)
    !
    !------------------------- values at current (T,P), for the phase --
    REAL(dp):: Grt,H,S,V,Cp
    !
  ENDTYPE T_MixPhaseData

CONTAINS

INTEGER FUNCTION MixPhase_Index(V,Str)
!--
!-- position of mixture phase named Str in V(1:SIZE(V))
!--
  TYPE(T_MixPhase),INTENT(IN)::V(:)
  CHARACTER(*),    INTENT(IN)::Str
  !
  INTEGER::I
  !
  MixPhase_Index=0
  IF(SIZE(V)==0) RETURN
  !
  I=0
  DO
    I=I+1 !; IF(iDebug>0) WRITE(fTrc,'(A)') vCpn(I)%SpName
    IF(TRIM(Str)==TRIM(V(I)%Name)) THEN
      MixPhase_Index= I
      EXIT
    ENDIF
    IF(I==SIZE(V)) EXIT
  ENDDO
  !if Str not found, I=0
  !
  RETURN
ENDFUNCTION MixPhase_Index

SUBROUTINE MixPhase_Zero(S)
  TYPE(T_MixPhase),INTENT(OUT):: S
  S%Name=      "Z"
  S%iModel=  0
  !
  S%vLPole=  .FALSE.
  S%vXPole=  Zero
  !S%vXAtom=  Zero
  !~ S%vLnAct=  Zero
  S%vLIdeal= Zero
  S%vLGam=   Zero
  !! S%vLMarg=  Zero
  !
ENDSUBROUTINE MixPhase_Zero

SUBROUTINE MixPhase_NormCompo(Mix)
!--
!-- normalize end-member fractions to tot-100%
!-- should not always normalize ???
!--
  TYPE(T_MixPhase),INTENT(INOUT):: Mix
  !
  REAL(dp),DIMENSION(1:MaxPole)::vX
  ! = composition: mole fractions of the end members
  REAL(dp):: Y,Eps
  !
  Eps= EPSILON(Y)
  !
  WHERE(.NOT. Mix%vLPole) Mix%vXPole=Zero
  !
  vX(1:MaxPole)=Mix%vXPole(1:MaxPole)
  Y= SUM(vX(1:MaxPole),MASK=Mix%vLPole(1:MaxPole))
  !
  IF(Y>Eps) Mix%vXPole(1:MaxPole)=vX(1:MaxPole)/Y
  !
ENDSUBROUTINE MixPhase_NormCompo

REAL(dp) FUNCTION MixPhase_Weit( &
& vSpc,      & !IN
& MM,        & !IN, mixing model
& Mix)         !IN, mixture phase, composition normalized
  USE M_T_Species, ONLY: T_Species
  !
  TYPE(T_Species), INTENT(IN),DIMENSION(:)::vSpc
  TYPE(T_MixModel),INTENT(IN)   :: MM
  TYPE(T_MixPhase),INTENT(IN)   :: Mix !phase
  !
  REAL(dp):: Weit
  INTEGER :: I
  !
  Weit=Zero
  DO I=1,MM%NPole
    IF(Mix%vLPole(I)) &
    & Weit= Weit &
    &     + Mix%vXPole(I) *vSpc(MM%vIPole(I))%WeitKg
  ENDDO
  !
  MixPhase_Weit= Weit
  !
ENDFUNCTION MixPhase_Weit

REAL(dp) FUNCTION MixPhase_Volume( & !not complete,
!& TdgK,Pbar, & !IN
& vSpc,      & !IN, database
& SM,        & !IN, solution model
& Fas)         !IN, solution phase, composition normalized
  USE M_T_Species, ONLY: T_Species
  !
  !REAL(dp),        INTENT(IN)   :: TdgK,Pbar
  TYPE(T_Species), INTENT(IN):: vSpc(:)
  TYPE(T_MixModel),INTENT(IN):: SM
  TYPE(T_MixPhase),INTENT(IN):: Fas !phase
  !
  REAL(dp),DIMENSION(1:MaxPole):: vX !composition: mole fractions of the end members
  REAL(dp):: V, V_XS
  INTEGER :: I
  !
  vX(1:SM%NPole)=Fas%vXPole(1:SM%NPole)
  V= Zero
  !
  DO I=1,SM%NPole
    IF(Fas%vLPole(I)) &
    & V= V + vX(I) *vSpc(SM%vIPole(I))%V0
  ENDDO
  !
  !V= DOT_PRODUCT(vX(1:SM%NPole),vSpc(SM%vIPole(1:SM%NPole))%V0)
  !------------------------------ XS volume calc. not implemented !!! --
  V_XS=Zero
  !
  MixPhase_Volume= V + V_XS
  !
ENDFUNCTION MixPhase_Volume

SUBROUTINE MixPhase_CalcMixing( & !
& TdgK,Pbar, & !IN
& MM,        & !IN
& Fas,       & !IN
& GMix,G_IdMix,G_XsMix)  !OUT
!--
!-- Output -
!--   G_IdMix: ideal part of free energy of mixing of mixture Phase at T,P
!--   GMix: free energy of mixing of phase Phase at T,P
!--   Gibbs free energy of the solution at T,P - GMeca + GMix
!--
  USE M_Dtb_Const,ONLY: R_jk
  !
  REAL(dp),        INTENT(IN) :: Pbar, TdgK
  TYPE(T_MixModel),INTENT(IN) :: MM
  TYPE(T_MixPhase),INTENT(IN) :: Fas
  !REAL(dp),        INTENT(IN) :: vLPole(:),vXPole(:)
  REAL(dp),        INTENT(OUT):: GMix,G_IdMix,G_XsMix
  !
  ! REAL(dp),ALLOCATABLE:: vXAtom(:)
  !
  G_IdMix= MixModel_GibbsIdeal( & !
  & TdgK,Pbar,  & !IN
  & MM,         & !IN, mixing model
  & Fas%vLPole, & !IN
  & Fas%vXPole)   !IN
  !
  GMix= MixModel_GibbsMixRT( & !
  & TdgK,Pbar,  & !IN
  & MM,         & !IN, mixing model
  & Fas%vLPole, & !IN
  & Fas%vXPole)   !IN
  GMix= GMix *R_jk*TdgK
  !
  G_XsMix= GMix - G_IdMix
  !
  RETURN
ENDSUBROUTINE MixPhase_CalcMixing

REAL(dp) FUNCTION MixPhase_GibbsRT( & !
& TdgK,Pbar, & !IN
& vSpc,      & !IN, database
& MM,        & !IN, mixing model
& Fas)         !IN, mixture phase, composition normalized
!--
!-- not complete, check the SITE models !!!
!--
  USE M_Dtb_Const,ONLY: R_jk
  USE M_T_Species,ONLY: T_Species
  !
  REAL(dp),        INTENT(IN):: TdgK,Pbar
  TYPE(T_Species), INTENT(IN):: vSpc(:)
  TYPE(T_MixModel),INTENT(IN):: MM
  TYPE(T_MixPhase),INTENT(IN):: Fas !phase
  !
  !~ REAL(dp),ALLOCATABLE:: vXAtom(:)
  REAL(dp):: Gmix,Gmeca
  INTEGER :: I
  !
  Gmix= MixModel_GibbsMixRT( & !
  & TdgK,Pbar,  & !IN
  & MM,         & !IN, mixing MM
  & Fas%vLPole, & !IN
  & Fas%vXPole)   !IN
  !
  !---------------------------------------------- "mechanical" mixing --
  Gmeca= Zero
  DO I=1,MM%NPole
    IF(Fas%vLPole(I)) &
    & Gmeca= Gmeca &
    &      + Fas%vXPole(I) *vSpc(MM%vIPole(I))%G0rt
  ENDDO
  !----------------------------------------------/"mechanical" mixing --
  !
  !~ IF(ALLOCATED(vXAtom)) DEALLOCATE(vXAtom)
  !
  MixPhase_GibbsRT= Gmix +Gmeca
  !
  RETURN
ENDFUNCTION MixPhase_GibbsRT

SUBROUTINE MixPhase_CalcActivs( & !
& TdgK,Pbar, & ! in
& MM,        & ! in:    mixing model
& F)           ! inout: mixture phase
!--
!-- calculate activities of end-members in phase F at given T,P
!--
  USE M_Dtb_Const,ONLY: R_jk
  !
  REAL(dp),        INTENT(IN)   :: Pbar, TdgK
  TYPE(T_MixModel),INTENT(IN)   :: MM
  TYPE(T_MixPhase),INTENT(INOUT):: F
  !
  LOGICAL :: Ok
  INTEGER :: N
  CHARACTER(LEN=80):: Msg
  !
  REAL(dp),ALLOCATABLE:: vLGam(:),vLIdeal(:),vLnAct(:) !!,vLMarg(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)' ) "< MixPhase_CalcActivs"
  !
  N= SIZE(F%vLnAct)
  ALLOCATE(vLGam(N),vLIdeal(N),vLnAct(N)) !!,vLMarg(N)
  !
  CALL MixModel_Activities( & !
  & TdgK,Pbar, & ! in
  & MM,        & ! in: mixing model
  & F%vXPole,  & ! in
  & F%vLPole,  & ! in
  & Ok, Msg,   & ! out
  & vLGam,     & ! out
  & vLIdeal,   & ! out
  & vLnAct)      ! out
  !
  F%vLIdeal(:)= vLIdeal(:)
  !! F%vLMarg(:)=  vLMarg(:)
  F%vLGam(:)=   vLGam(:)
  F%vLnAct(:)=  vLnAct(:)
  !
  DEALLOCATE(vLGam,vLIdeal,vLnAct) !!,vLMarg
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)' ) "</ MixPhase_CalcActivs"
  !
ENDSUBROUTINE MixPhase_CalcActivs

SUBROUTINE MixPhase_Constraint( &
& S,         & ! in: mixing model
& F,         & ! in: mixture phase
& Ok, Msg,   & ! out
& vC)          !
  TYPE(T_MixModel),INTENT(IN) :: S
  TYPE(T_MixPhase),INTENT(IN) :: F
  LOGICAL,         INTENT(OUT):: Ok
  CHARACTER(*),    INTENT(OUT):: Msg
  REAL(dp),        INTENT(OUT):: vC(:)
  !
  INTEGER:: I
  !
  SELECT CASE(TRIM(S%Model))

  CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR" or "FELSPAR"
    vC(1)= One
    DO I=1,S%NPole
      IF(F%vLPole(I)) &
      & vC(1)= vC(1) - F%vXPole(I)
    ENDDO
    !
  CASE("SITE")
    vC(1)= One
    DO I=1,S%NPole
      IF(F%vLPole(I)) &
      & vC(1)= vC(1) - F%vXPole(I)
    ENDDO
    !
  CASE DEFAULT
    Ok= .FALSE.
    Msg= TRIM(S%Model)//"= invalid S%Model in MixPhase_Constraint"

  END SELECT
  !
ENDSUBROUTINE MixPhase_Constraint

SUBROUTINE MixPhase_ConstraintGrad( &
& S,         & ! in: mixing model
& F,         & ! in: mixture phase
& Ok, Msg,   & ! out
& vGC)         !
  TYPE(T_MixModel),INTENT(IN) :: S
  TYPE(T_MixPhase),INTENT(IN) :: F
  LOGICAL,         INTENT(OUT):: Ok
  CHARACTER(*),    INTENT(OUT):: Msg
  REAL(dp),        INTENT(OUT):: vGC(:,:)
  !
  INTEGER:: I
  !
  vGC= Zero
  !
  SELECT CASE(TRIM(S%Model))

  CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
  !! should consider "IDEAL","POLE" as obsolete,
  !! and recommand using "MOLECULAR" or "FELSPAR"
    DO I=1,S%NPole
      IF(F%vLPole(I)) vGC(1,I)= -One
    ENDDO
    !
  CASE("SITE")
    DO I=1,S%NPole
      IF(F%vLPole(I)) vGC(1,I)= -One
    ENDDO
    !
  CASE DEFAULT
    Ok= .FALSE.
    Msg= TRIM(S%Model)//"= invalid S%Model in MixPhase_ConstraintGrad"

  END SELECT
  !
ENDSUBROUTINE MixPhase_ConstraintGrad

ENDMODULE M_T_MixPhase

!~ SUBROUTINE MixPhase_CalcActivities( & !
!~ & TdgK,Pbar, & ! in
!~ & MM,         & ! in: mixing model
!~ & F,         & ! in: mixture phase
!~ & Ok, Msg,   & ! out
!~ & vLGam,     & !
!~ !! & vLMarg,    & !
!~ & vLIdeal,   & !
!~ & vLnAct)      !
!~ !--
!~ !-------- calculate activities of end-members in phase F at given T,P --
!~ !--
  !~ USE M_Dtb_Const,ONLY: R_jk
  !~ USE M_MixModel_Special
  !~ !
  !~ REAL(dp),        INTENT(IN) :: Pbar, TdgK
  !~ TYPE(T_MixModel),INTENT(IN) :: MM
  !~ TYPE(T_MixPhase),INTENT(IN) :: F
  !~ LOGICAL,         INTENT(OUT):: Ok
  !~ CHARACTER(*),    INTENT(OUT):: Msg
  !~ REAL(dp),        INTENT(OUT):: vLGam(:)
  !~ !! REAL(dp),        INTENT(OUT):: vLMarg(:)
  !~ REAL(dp),        INTENT(OUT):: vLIdeal(:)
  !~ REAL(dp),        INTENT(OUT):: vLnAct(:)
  !~ !
  !~ INTEGER :: iP,iM
  !~ REAL(dp):: P
  !~ REAL(dp),ALLOCATABLE:: vMonome(:)
  !~ !
  !~ P=Pbar !for future use ??
  !~ !
  !~ Ok= .TRUE.
  !~ Msg= "Ok"
  !~ !
  !~ !F%vLPole(1:MM%NPole)= vX(1:MM%NPole)>Zero
  !~ !
  !~ vLGam(:)=   Zero !default
  !~ vLIdeal(:)= Zero !default
  !~ !
  !~ IF(TRIM(MM%Model)=="SPECIAL") THEN
    !~ !
    !~ CALL MixModel_Special_Activities( &
    !~ & MM%Name,      &
    !~ & TdgK, Pbar,   &
    !~ & F%vXpole,     &
    !~ & F%vLPole,     &
    !~ & vLIdeal,      &
    !~ & vLGam         )
    !~ !
  !~ ELSE

  !~ SELECT CASE(TRIM(MM%Model))
    !~ !
    !~ CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
      !~ CALL MixModel_Pole_LnActivsIdeal(MM,F%vXpole,F%vLpole,vLIdeal)

    !~ CASE("SITE")
      !~ CALL MixModel_Site_LnActivsIdeal(MM,F%vXatom,F%vLpole,Ok,vLIdeal)

    !~ CASE DEFAULT
      !~ Ok= .FALSE.
      !~ Msg= TRIM(MM%Model)//"= invalid MM%Model in MixPhase_CalcActivities"

    !~ END SELECT
    !~ !
    !~ vLGam(1:MM%NPole)=Zero
    !~ !
    !~ !-------------------------- activ coeff related to Margules Terms --
    !~ IF(MM%NMarg>0) THEN
      !~ ALLOCATE(vMonome(MM%NMarg))
      !~ !
      !~ SELECT CASE(TRIM(MM%Model))
      !~ !
      !~ CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
        !~ DO iM=1,MM%NMarg
          !~ vMonome(iM)= MixModel_Margules_Monome(MM%vMarg(iM),F%vXPole)
        !~ ENDDO
        !~ DO iP=1,MM%NPole
          !~ IF(F%vLPole(iP)) THEN
            !~ vLGam(iP)= MixModel_Pole_LnGammaMargules(MM,iP,vMonome,F%vXPole) /R_jk/TdgK
          !~ ENDIF
        !~ ENDDO
      !~ !
      !~ CASE("SITE")
        !~ DO iM=1,MM%NMarg
          !~ vMonome(iM)= MixModel_Margules_Monome(MM%vMarg(iM),F%vXatom)
        !~ ENDDO
        !~ DO iP=1,MM%NPole
          !~ IF(F%vLPole(iP)) THEN
            !~ vLGam(iP)= MixModel_Site_LnGammaMargules(MM,iP,vMonome,F%vXAtom) /R_jk/TdgK
          !~ ENDIF
        !~ ENDDO
      !~ !
      !~ CASE DEFAULT
        !~ Ok= .FALSE.
        !~ Msg= TRIM(MM%Model)//"= invalid MM%Model in MixPhase_CalcActivities"
      !~ !
      !~ END SELECT
      !~ !
      !~ DEALLOCATE(vMonome)
      !~ !
    !~ ENDIF
    !~ !-------------------------/ activ coeff related to Margules Terms --
    !~ !
  !~ ENDIF
  !~ !
  !~ vLnAct(1:MM%NPole)= vLIdeal(1:MM%NPole) + vLGam(1:MM%NPole)
  !~ !
  !~ IF(iDebug>0) THEN !------------------------------------------ trace --
    !~ WRITE(fTrc,'(A)') "MixPhase_CalcActivs -> X,ActIdeal,Gamma,Activ"
    !~ WRITE(fTrc,'(4A)') "Phase=",F%Name, "MixModel=",MM%Name
    !~ DO iP=1,MM%NPole
      !~ IF(F%vLPole(iP)) THEN
        !~ WRITE(fTrc,'()')
        !~ WRITE(fTrc,'(A,I2,A1,A15,A1,4(A4,G11.6,A1))') &
        !~ & "POLE",iP,           T_,&
        !~ & TRIM(MM%vNamPole(iP)),T_,&
        !~ & "Frc=",F%vXPole(iP), T_,&
        !~ & "XId=",EXP(vLIdeal(iP)),T_,&
        !~ & "Gam=",EXP(vLGam(iP)),  T_,&
        !~ & "Act=",EXP(vLnAct(iP)), T_
      !~ ENDIF
    !~ ENDDO
  !~ ENDIF !-----------------------------------------------------/ trace --
  !~ !
!~ ENDSUBROUTINE MixPhase_CalcActivities

!~ SUBROUTINE MixPhase_CalcMixing_( & !
!~ & TdgK,Pbar, & !IN
!~ & MM,        & !IN
!~ & Fas,       & !IN
!~ & GMix,G_IdMix,G_XsMix)  !OUT
!~ !--
!~ !-- Output -
!~ !--   G_IdMix: ideal part of free energy of mixing of mixture Phase for a given P,T
!~ !--   GMix: free energy of mixing of phase Phase for a given P,T
!~ !--   Gibbs free energy of the solution at P,T - GMeca + GMix
!~ !--
  !~ USE M_Dtb_Const,ONLY: R_jk
  !~ !
  !~ REAL(dp),        INTENT(IN) :: Pbar, TdgK
  !~ TYPE(T_MixModel),INTENT(IN) :: MM
  !~ TYPE(T_MixPhase),INTENT(IN) :: Fas
  !~ !REAL(dp),        INTENT(IN) :: vLPole(:),vXPole(:)
  !~ REAL(dp),        INTENT(OUT):: GMix,G_IdMix,G_XsMix
  !~ !
  !~ REAL(dp),ALLOCATABLE:: vXAtom(:)
  !~ !
  !~ G_XsMix= Zero
  !~ G_XsMix= Zero
  !~ !
  !~ IF(TRIM(MM%Model)=="SPECIAL") THEN

  !~ ELSE
    !~ !
    !~ SELECT CASE(TRIM(MM%Model))
    !~ !
    !~ CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
      !~ G_IdMix= -TdgK *MixModel_Pole_SConf(MM,Fas%vLPole,Fas%vXPole)
      !~ !
      !~ IF(MM%NMarg>0) G_XsMix= &
      !~ & MixModel_Pole_XsMargules(MM, Fas%vLPole, Fas%vXpole)
    !~ !
    !~ CASE("SITE")
      !~ ALLOCATE(vXAtom(MM%NAtom))
      !~ CALL MixModel_XPoleToXSite(MM,Fas%vXPole,vXAtom)
      !~ !
      !~ G_IdMix= -TdgK *MixModel_Site_SConf(MM,vXAtom)
      !~ !
      !~ IF(MM%NMarg>0) G_XsMix= &
      !~ & MixModel_Site_XsMargules(MM,vXAtom)
    !~ !
    !~ END SELECT
    !~ !
  !~ ENDIF
  !~ !
  !~ IF(ALLOCATED(vXAtom)) DEALLOCATE(vXAtom)
  !~ !
  !~ GMix= G_IdMix +G_XsMix
  !~ !
  !~ RETURN
!~ ENDSUBROUTINE MixPhase_CalcMixing_

