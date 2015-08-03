MODULE M_System_Tools
  USE M_Kinds
  USE M_Trace,ONLY: fTrc,T_,Stop_,iDebug,Pause_,fHtm
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: System_Build
  PUBLIC:: System_Build_Custom
  PUBLIC:: System_TP_Update
  !
CONTAINS

SUBROUTINE System_Build
!--
!-- build "master" system --
!-- system with all components (elements) involved in any system in the run
!--
  !!--tools--!!
  USE M_T_Species,    ONLY: T_Species,Species_Stoikio_Calc
  USE M_Global_Alloc, ONLY: MixModels_Alloc,MixPhases_Alloc,Phases_Alloc !_New
  USE M_SolModel_Alloc
  USE M_T_SolModel,   ONLY: SolModel_Spc_Init
  USE M_T_SolPhase,   ONLY: SolPhase_Init
  USE M_Solmodel_Read, ONLY: Solmodel_Solvent_Read
  USE M_SolModel_Tools,ONLY: SolModel_TP_Update
  USE M_Solmodel_Pitzer_Dtb,ONLY: Solmodel_Pitzer_Dtb_Init,Solmodel_Pitzer_Dtb_TPtest
  USE M_Global_Alloc, ONLY: DiscretParam_Alloc
  USE M_DiscretModel_Read
  USE M_DiscretModel_Tools
  !
  !--global variables--
  USE M_Global_Vars,  ONLY: vEle,vSpc,vMixFas,tFormula
  USE M_Global_Vars,  ONLY: vMixModel,vDiscretModel,vDiscretParam
  USE M_Global_Vars,  ONLY: SolModel
  USE M_Global_Vars,  ONLY: vSolModel,vSolFas
  USE M_Global_Vars,  ONLY: nAq,nMn,nGs
  !
  !--database variables--
  USE M_Dtb_Vars,     ONLY: DtbFormat,DtbLogK_Dim,DtbLogK_vTPCond
  USE M_Solmodel_Vars, ONLY: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  USE M_Solmodel_Vars, ONLY: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl,T_Spline
  !
  !--system variables--
  USE M_System_Vars,  ONLY: vCpn,TdgK,Pbar
  USE M_System_Vars,  ONLY: CpnIsSpc
  USE M_System_Vars,  ONLY: System_Zero,System_Type
  !
  INTEGER:: I,N
  LOGICAL:: Ok,fOk
  CHARACTER(LEN=80):: Msg
  REAL(dp),ALLOCATABLE:: vTdgC(:)
  TYPE(T_Species),ALLOCATABLE:: vSpcTmp(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< System_Build"
  !
  CALL System_Zero
  !
  !-- read local element list (-current system) -> build vCpn
  CALL Components_Alloc( & !
  & "SYSTEM",            & !
  & vEle,vSpc,vMixFas,   & !in
  & TdgK,Pbar,           & !
  & System_Type)           !out
  !
  !-- rebuild vEle, with only elements pointed to by vCpn(1:N)
  CALL Elements_Alloc_forSystem(vCpn)
  !
  CALL System_TP_Check(TdgK,Pbar,Ok,Msg)
  IF(.not. Ok) CALL Stop_(TRIM(Msg))
  !
  IF(iDebug>0) THEN
    DO i=1,SIZE(vCpn)
      WRITE(fTrc,'(2A,2(G15.6,1X))') &
      & vCpn(i)%NamCp,vCpn(i)%Statut,vCpn(i)%Mole,vCpn(i)%LnAct
    ENDDO
  ENDIF
  !
  !-- rebuild vSpc, updating vCpn%iSpc
  CALL Species_Alloc_forSystem(vEle,vCpn)
  !
  !--------------------------------------------- update vSpc(:)%vStoikio
  CALL Species_Stoikio_Calc(vEle,vSpc,fOk)
  !--------------------------------------------/ update vSpc(:)%vStoikio
  !
  !-------------------------- read solution models -> allocate vMixModel
  CALL MixModels_Alloc(vSpc) !
  !
  !------------------------------------------------ add discrete species
  !--- allocate & read vDiscretModel--
  CALL DiscretModel_Read(vMixModel)
  !
  IF(SIZE(vDiscretModel)>0) THEN
    !
    CALL DiscretParam_Alloc(vDiscretModel) !-> allocate vDiscretParam
    !
    ALLOCATE(vSpcTmp(SIZE(vDiscretParam)))
    !
    CALL DiscretParam_Init( &
    & vEle,vSpc,vMixModel,vDiscretModel, &
    & vDiscretParam,vSpcTmp) !-> build vSpcDiscret
    !
    CALL Species_Append(vSpcTmp) !-> new vSpc !!!
    !
    DEALLOCATE(vSpcTmp)
    !
    CALL DiscretSpecies_Stoikio_Calc( & !
    & vEle,          & !IN
    & vMixModel,     & !IN
    & vDiscretModel, & !IN
    & vDiscretParam, & !IN
    & vSpc)            !INOUT
    !
  ENDIF
  !-----------------------------------------------/ add discrete species
  !
  nAq= COUNT(vSpc%Typ=="AQU")
  nMn= COUNT(vSpc%Typ=="MIN")
  nGs= COUNT(vSpc%Typ=="GAS")
  nMn= nMn + nGs
  !
  !! !--- NEW2010-10-29
  !! N= SIZE(vCpn)
  !! IF(CpnIsSpc) THEN
  !!   DO i=1,N
  !!     vCpn(i)%iEle= I
  !!     J= vCpn(i)%iSpc
  !!     vCpn(i)%vStoikCp(0:N+1)= vSpc(J)%vStoikio(0:N+1)
  !!   ENDDO
  !! ELSE
  !!   DO i=1,N
  !!     vCpn(i)%iEle= I
  !!     vCpn(i)%vStoikCp(0)= 1 !formula divider !!
  !!     vCpn(i)%vStoikCp(:)= 0
  !!     vCpn(i)%vStoikCp(i)= 1
  !!   ENDDO
  !! ENDIF
  !! !---/ NEW2010-10-29
  !
  !------------------------------------------------ update formula table
  DEALLOCATE(tFormula)
  ALLOCATE(tFormula(1:SIZE(vEle),1:SIZE(vSpc)))
  CALL FormulaTable_Calc(vSpc,tFormula)
  !
  IF(iDebug>0) CALL FormulaTable_Sho(vEle,vSpc,tFormula)
  !-----------------------------------------------/ update formula table
  !
  !--- read phase compositions, allocate vMixFas --
  CALL MixPhases_Alloc(vSpc,vMixModel)
  CALL MixPhase_CheckFound(vEle,vSpc,vMixModel,vMixFas,vCpn)
  !
  !------------------------------------------------- initialize SolModel
  IF(System_Type=="AQUEOUS") THEN
    
    N= 0
    IF(DtbFormat=="LOGKTBL") N= DtbLogK_Dim
    
    ALLOCATE(vTdgC(N))
    IF(N>0) vTdgC(1:N)= DtbLogK_vTPCond(1:N)%TdgC
    
    ! even for HSV base, Solmodel_Read must be called,
    ! for reading activity model
    CALL Solmodel_Solvent_Read( &
    & N,vTdgC,SolModel, &
    & Ok_Rho, Ok_Eps, Ok_DHA, Ok_DHB, Ok_BDot, &
    & Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl)

    DEALLOCATE(vTdgC)

    IF(SolModel%iActModel== 8  .OR. &    !! "PITZER"
    &  SolModel%iActModel== 11 ) THEN    !! "SIT"
      CALL Solmodel_Pitzer_Dtb_Init(vSpc)
    ENDIF

    CALL SolModel_TP_Update(TdgK,Pbar,SolModel)

    IF(iDebug>2) THEN
      IF(SolModel%iActModel== 8  .OR. &    !! "PITZER"
      &  SolModel%iActModel== 11 ) &    !! "SIT"
      & CALL Solmodel_Pitzer_Dtb_TPtest !!!"WRK"!!!
    ENDIF

    !----------------------------------------------------------------NEW
    !--- initialize indexes of the sol'model
    CALL SolModel_Spc_Init("H2O",vSpc,SolModel,Ok,Msg)

    CALL SolModel_Alloc ! allocate vSolModel
    vSolModel(1)= SolModel

    CALL SolPhase_Alloc ! allocate vSolFas
    CALL SolPhase_Init( & !
    & 1,          & !IN, model index
    & vSolModel,  & !IN, base of solution models
    & vSolFas(1), & !INOUT, solution phase
    & Ok,Msg)       !OUT

    !! CALL SolPhase_Model_Init(SolModel%Model,vSolModel,vSolFas(1),Ok,Msg)
    !! IF(.NOT. Ok) CALL Stop_(TRIM(Msg))

    IF(.NOT. Ok) &
    & CALL Stop_("SolModel_Spc_Init "//TRIM(Msg))
    !---------------------------------------------------------------/NEW

  ENDIF
  !-------------------------------------------------/initialize SolModel
  !
  !--- allocate vFas >--
  ! CALL Phases_Alloc_New(vSpc,vMixFas,vSolFas)
  CALL Phases_Alloc(vSpc,vMixFas)
  !
  !------------------------------------- update global thermo'parameters
  !--------------------- (vSpc,vMixModel,vMixFas,vFas) at (TdgK,Pbar) --
  CALL System_TP_Update(TdgK,Pbar)
  !------------------------------------/ update global thermo'parameters
  !
  IF(System_Type=="AQUEOUS") THEN
    !-- check presence of species H2O, H+, OH- in vSpc --
    CALL System_Species_Check(vSpc,Ok,Msg)
    IF(.NOT. Ok) CALL Stop_(Msg)
  END IF
  !
  !! CALL Basis_Init(vCpn,TdgK,Pbar)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ System_Build"
  !
ENDSUBROUTINE System_Build

SUBROUTINE Species_Append(vSpcAdd)
!--
!-- Append vSpcAdd to current vSpc
!-- -> produce a new vSpc
!--
  USE M_T_Species, ONLY: T_Species
  !
  USE M_Global_Vars, ONLY: vSpc
  !
  TYPE(T_Species),INTENT(IN):: vSpcAdd(:)
  !
  TYPE(T_Species),ALLOCATABLE:: vSpcAll(:)
  INTEGER:: M, N
  !
  M= SIZE(vSpc)
  N= SIZE(vSpcAdd)
  !
  ALLOCATE(vSpcAll(M+N))
  vSpcAll(1   :M )= vSpc(1:M)
  vSpcAll(M+1:M+N)= vSpcAdd(1:N)
  !
  DEALLOCATE(vSpc)
  ALLOCATE(vSpc(N+M)) ; vSpc= vSpcAll
  !
  DEALLOCATE(vSpcAll)
  !
ENDSUBROUTINE Species_Append
  !
SUBROUTINE Species_Alloc_forSystem(vEle,vCpn)
!--
!-- from general vSpc built from database,
!-- reduce vSpc to species involved in the run
!--
  USE M_Numeric_Const,       ONLY: Ln10
  USE M_T_Component, ONLY: T_Component
  USE M_T_Element,   ONLY: T_Element,Element_Index,Formula_Read
  USE M_T_Species,   ONLY: T_Species,Species_Index
  !
  USE M_Global_Vars,ONLY: vSpc,tFormula
  USE M_System_Vars,ONLY: System_Type
  !
  TYPE(T_Element),   INTENT(INOUT):: vEle(:)
  TYPE(T_Component), INTENT(INOUT):: vCpn(:)
  !
  CHARACTER(LEN=23),DIMENSION(:),ALLOCATABLE:: vNamSpc
  !
  INTEGER,ALLOCATABLE:: vStoik(:)
  !! LOGICAL,ALLOCATABLE:: vExclude(:)
  !! LOGICAL,ALLOCATABLE:: vInclude(:)
  !
  TYPE(T_Species),ALLOCATABLE:: vSpcNew(:)
  !
  TYPE(T_Species):: S_
  LOGICAL        :: fOk
  INTEGER        :: nCp,iCp,nSp,iSp,iOx_,ZSp,Z_,i,nAq,nMn,nGs
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "< Species_Alloc"
  !
  nCp=SIZE(vCpn)
  !
  !before restructuring vSpc, save vCpn-vSpc links
  ALLOCATE(vNamSpc(1:nCp))
  DO iCp=1,nCp
    vNamSpc(iCp)= TRIM(vSpc(vCpn(iCp)%iSpc)%NamSp)
  ENDDO
  !
  CALL Redox_Calc(vSpc, vEle,vCpn)
  !
  !! !------------------ read SPECIES.INCLUDE and SPECIES.EXCLUDE blocks --
  !! ALLOCATE(vExclude(1:SIZE(vSpc)))  ;  vExclude(:)= .FALSE.
  !! ALLOCATE(vInclude(1:SIZE(vSpc)))  ;  vInclude(:)= .TRUE.
  !! CALL Species_Read_Excluded(vSpc,vExclude,vInclude)
  !! !-----------------/ read SPECIES.INCLUDE and SPECIES.EXCLUDE blocks --
  !
  !------------------------------------------------ species selection --
  !-- from the vSpc built from database,
  !-- retrieve species consistent
  !-- with vEle & redox
  IF(iDebug>0) WRITE(fTrc,'(A)') &
  & "< Restrict database to species consistent with element list and redox state"
  !
  ALLOCATE(vSpcNew(SIZE(vSpc)))
  !
  ALLOCATE(vStoik(1:nCp))
  iOx_= Element_Index("OX_",vEle)
  nSp=0
  DO iSp=1,SIZE(vSpc)
    !
    !! IF( vExclude(iSp) .OR. (.NOT. vInclude(iSp)) ) CYCLE
    !
    S_=vSpc(iSp)
    IF(S_%NamSp(1:1)=="!") CYCLE !!!§§§ bricollage !!!
    !
    CALL Formula_Read(S_%Formula,vEle,ZSp,S_%Div,fOk,vStoik)
    !
    IF(fOk) THEN
    !fOk=FormulaIsOk (as regards elements in run), build new LnkSpc, called Lnk
      !
      Z_=DOT_PRODUCT(vStoik(1:nCp),vEle(1:nCp)%Z) !check charge of species
      !
      IF(iOx_/=0 .OR. ZSp==Z_) THEN !IF charge is OK, then save in vSpcNew
        !
        nSp=  nSp+1
        S_%Z= ZSp
        !
        vSpcNew(nSp)= S_
        !
        IF(iDebug>0) WRITE(fTrc,'(A,A1,A15,A1,A39,A1,2I3)') &
        & "ACCEPT",T_,S_%NamSp,T_,S_%Formula,T_,S_%Z,nSp
      ENDIF
      !
    ENDIF
    !
  ENDDO
  !
  !! DEALLOCATE(vExclude)
  !! DEALLOCATE(vInclude)
  !! DEALLOCATE(vStoik)
  !-----------------------------------------------/ species selection --
  !
  !-------------------------------------------- build new sorted vSpc --
  DEALLOCATE(vSpc)
  ALLOCATE(vSpc(1:nSp))
  !
  nAq= 0; nMn= 0; nGs= 0
  iSp= 0
  DO i=1,nSp
    IF(vSpcNew(i)%Typ=="AQU") THEN
      iSp= iSp +1
      vSpc(iSp)= vSpcNew(i)
    ENDIF
  ENDDO
  DO i=1,nSp
    IF(vSpcNew(i)%Typ=="MIN") THEN
      iSp= iSp +1
      vSpc(iSp)= vSpcNew(i)
    ENDIF
  ENDDO
  DO i=1,nSp
    IF(vSpcNew(i)%Typ=="GAS") THEN
      iSp= iSp +1
      vSpc(iSp)= vSpcNew(i)
    ENDIF
  ENDDO
  DEALLOCATE(vSpcNew)
  !-------------------------------------------/ build new sorted vSpc --
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Restrict database"
  !
  IF(COUNT(vSpc%Typ=="AQU")==0 .AND. System_Type=="AQUEOUS") &
  & CALL Stop_("Found NO Aqueous Species ...") !=================<STOP==
  !
  !-----------------------------------------then relate vCpn with vSpc--
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') &
  & "Species_Alloc_forSystem, Element/ Component/ Status/ iSpc"
  !
  !vSpc may have changed
  !-> need to check whether the species pointed to by Cpn's are in the new vSpc
  !
  DO iCp=1,nCp
  !---------------------find which vSpc(iSp) matches with vNamSpc(iCp)--
    !
    iSp= Species_Index(vNamSpc(iCp),vSpc) !-> index of Prim'species in new vSpc
    !
    IF(iSp==0) CALL Stop_(TRIM(vNamSpc(iCp))//" Not Found as Species")
    !
    IF(System_Type=="AQUEOUS") THEN
      IF(vSpc(iSp)%Typ /= "AQU") THEN
        IF(TRIM(vCpn(iCp)%Statut)=="INERT") &
        & CALL Stop_( & ! !===================================== STOP ==
        & TRIM(vNamSpc(iCp))//" : ALL INERT SPECIES SHOULD BE AQUEOUS")
      ENDIF
    ENDIF
    !
    !update %iSpc (link between vSpc and vCpn)
    vCpn(iCp)%iSpc= iSp
    !
    IF(iDebug>0) WRITE(fTrc,'(A3,A23,A7)') &
    & vEle(iCp)%NamEl,vNamSpc(iCp),vCpn(iCp)%Statut
    !
  ENDDO
  !
  DEALLOCATE(vNamSpc)
  !
  IF(iDebug>1) THEN
    WRITE(fTrc,'(/,A,/)') "Species list, final"
    DO iSp=1,nSp
      WRITE(fTrc,'(I7,A1,A23,A1,A3)') iSp,T_,vSpc(iSp)%NamSp,T_,vSpc(iSp)%Typ
    ENDDO
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "</ Species_Alloc"
  !
ENDSUBROUTINE Species_Alloc_forSystem

SUBROUTINE System_Species_Check(vSpc,Ok,Msg)
!--
!-- check presence of species H2O, H+, OH- in vSpc --
!-- (for an aqueous system)
!--
  USE M_T_Species,ONLY: T_Species,Species_Index
  !
  TYPE(T_Species),  INTENT(IN)  :: vSpc(:)
  LOGICAL,          INTENT(OUT) :: Ok
  CHARACTER(LEN=80),INTENT(OUT):: Msg
  !
  Ok=  .true.
  Msg= "Ok"
  !
  IF(Species_Index("H2O",vSpc)==0) THEN
    Msg= "species H2O not found !!!"
    Ok=  .false.
    RETURN
  ENDIF
  IF(Species_Index("OH-",vSpc)==0 .and. Species_Index("OH[-]",vSpc)==0) THEN
    Msg= "species OH- not found !!!"
    Ok=  .false.
    RETURN
  ENDIF
  IF(Species_Index("H+", vSpc)==0 .and. Species_Index("H[+]",vSpc)==0) THEN
    Msg= "species H+ not found !!!"
    Ok=  .false.
    RETURN
  ENDIF
  !
ENDSUBROUTINE System_Species_Check

SUBROUTINE System_TP_Check(TdgK,Pbar,Ok,Msg)
  USE M_T_Tpcond,ONLY: T_TPCond
  USE M_Dtb_Vars,ONLY: DtbFormat,DtbLogK_vTPCond,Psat_Auto
  USE M_Dtb_Calc,ONLY: Dtb_TP_Check
  !
  REAL(dp),        INTENT(INOUT):: TdgK,Pbar
  LOGICAL,         INTENT(OUT)  :: Ok
  CHARACTER(LEN=*),INTENT(OUT)  :: Msg
  !
  Ok=  .true.
  Msg= "OK"
  CALL Dtb_TP_Check(DtbFormat,DtbLogK_vTPCond,Psat_Auto,TdgK,Pbar,Ok,Msg)
  !
ENDSUBROUTINE System_TP_Check

SUBROUTINE System_TP_Update(TdgK,Pbar)
!--
!-- update thermo parameters --
!--
  USE M_Global_Tools,  ONLY: Global_TP_Update
  USE M_SolModel_Tools,ONLY: SolModel_TP_Update
  !
  USE M_Global_Vars, ONLY: vSpcDtb,vSpc,vMixModel,vDiscretModel,vDiscretParam
  USE M_Global_Vars, ONLY: vMixFas,vFas,SolModel
  USE M_System_Vars, ONLY: System_Type
  !
  REAL(dp),INTENT(IN):: TdgK,Pbar

  CALL Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)

  IF(System_Type=="AQUEOUS") CALL SolModel_TP_Update(TdgK,Pbar,SolModel)

ENDSUBROUTINE System_TP_Update

SUBROUTINE MixPhase_CheckFound( &
& vEle,vSpc,vMixModel, &
& vMixFas,vCpn)
!--
!-- check whether phases have been found
!-- for all components with %namSol /- "AQU"
!--
  USE M_T_Component,ONLY: T_Component
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Species,  ONLY: T_Species
  USE M_T_MixModel, ONLY: T_MixModel,MixModel_Site_ActivIdeal,MixModel_XPoleToXSite
  USE M_T_MixPhase, ONLY: T_MixPhase, MixPhase_Index
  !
  TYPE(T_Element),  INTENT(IN)   :: vEle(:)
  TYPE(T_Species),  INTENT(IN)   :: vSpc(:)
  TYPE(T_MixModel), INTENT(IN)   :: vMixModel(:)
  TYPE(T_MixPhase), INTENT(INOUT):: vMixFas(:) !-> update vMixFas%vLPole
  TYPE(T_Component),INTENT(INOUT):: vCpn(:) !-> update vCpn%iPol, vCpn%iMix
  !
  TYPE(T_MixModel):: SM
  INTEGER :: nCp,iSp,iCp,iMix,J,K
  REAL(dp):: Act
  LOGICAl, ALLOCATABLE:: Phase_Found(:)
  REAL(dp),ALLOCATABLE:: vXAtom(:)
  !-> check that phase for each component is described in input
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "< MixPhases_CheckFound"
  !
  nCp= SIZE(vCpn)
  !
  ALLOCATE(Phase_Found(1:nCp))  ;  Phase_Found=.FALSE.
  !
  DO iCp=1,nCp !check components
    !
    iSp=vCpn(iCp)%iSpc
    !
    IF(iDebug>0) WRITE(fTrc,'(I3,6A)') &
    & iCp, &
    & " >> ELE=",  vEle(vCpn(iCp)%iEle)%NamEl, &
    & " >> SPC=",  vSpc(vCpn(iCp)%iSpc)%NamSp, &
    & " >> PHASE=",TRIM(vCpn(iCp)%namSol)
    !
    IF(iSp==0) CALL Stop_ &
    & ("NO SPECIES FOR COMPONENT "//TRIM(vEle(vCpn(iCp)%iEle)%NamEl))
    !
    IF(TRIM(vCpn(iCp)%namSol)=="Z") THEN
      !
      Phase_Found(iCp)=.TRUE.
      !
    ELSE
      !
      iMix= MixPhase_Index(vMixFas,vCpn(iCp)%namSol)
      !
      IF(iMix==0) CALL Stop_ & !------------------------------------STOP
      & ("PHASE NOT FOUND FOR COMPONENT "//TRIM(vEle(vCpn(iCp)%iEle)%NamEl))
      !
      Phase_Found(iCp)=.TRUE.
      vCpn(iCp)%iMix= iMix !!vSpc(iSp)%iMix=iMix
      !
      IF(iDebug>0) WRITE(fTrc,'(4(A,1X))') &
      & "COMPONENT", vEle(vCpn(iCp)%iEle)%NamEl, &
      & "PHASE", TRIM(vMixFas(vCpn(iCp)%iMix)%Name)
      !
      !----------------- find index of "species" in "solution"%vIPole --
      SM= vMixModel(vMixFas(iMix)%iModel)
      J=0
      DO K=1,SM%NPole
        IF(SM%vIPole(K)==iSp) J=K
      ENDDO
      !-----------------/
      !
      IF(iDebug>0) WRITE(fTrc,'(2A,/)') &
      & "END-MEMB  ", TRIM(vSpc(iSp)%NamSp)
      !
      IF(J==0) CALL Stop_& !----------------------------------------STOP
      & ("SPECIES NOT IN MIXTURE:"//TRIM(vSpc(iSp)%NamSp))
      !
      vCpn(iCp)%iPol=J
      !--- iPol- index of species in end-member list of phase
      IF(TRIM(SM%Model)=="SITE") THEN
        ALLOCATE(vXAtom(SM%NAtom))
        CALL MixModel_XPoleToXSite( &
        & SM, vMixFas(iMix)%vXPole, &
        & vXAtom)
        !
        Act= MixModel_Site_ActivIdeal(SM,J,vXAtom)
        !
        IF(iDebug>0) &
        & WRITE(fTrc,'(A,A,G15.6)') TRIM(vSpc(iSp)%NamSp),">> ACT=",Act
        !
        IF(Act/=0.0D0) vMixFas(iMix)%vLPole=.TRUE.
        !
        DEALLOCATE(vXAtom)
      ENDIF
      !
    ENDIF
    !
  ENDDO
  !
  DO J=1,nCp
    IF(.NOT. Phase_Found(J)) & !------------------------------------STOP
    & CALL Stop_( "DESCRIPTION NOT FOUND FOR "//TRIM(vCpn(J)%namSol) )
  ENDDO
  !
  DEALLOCATE(Phase_Found)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "</ MixPhases_CheckFound"
ENDSUBROUTINE MixPhase_CheckFound

SUBROUTINE Components_Alloc(  & !
!--
!-- read SYSTEM block, build vCpn
!-- (using temporary vCpn)
!--
& sKeyWord,vEle,vSpc,vMixFas, & !IN
& TdgK,Pbar,                  & !INOUT
& SysType)                      !OUT
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Component,ONLY: T_Component,Component_Zero
  USE M_T_Tpcond,   ONLY: T_TPcond
  USE M_T_MixPhase, ONLY: T_MixPhase
  USE M_Component_Read,ONLY: Components_Read
  !
  USE M_System_Vars,ONLY: vCpn
  !
  CHARACTER(LEN=*),INTENT(IN):: sKeyWord
  TYPE(T_Element), INTENT(IN):: vEle(:)
  TYPE(T_Species), INTENT(IN):: vSpc(:)
  TYPE(T_MixPhase),INTENT(IN):: vMixFas(:)
  !
  REAL(dp),        INTENT(INOUT):: TdgK,Pbar
  CHARACTER(LEN=7),INTENT(OUT)  :: SysType
  !
  INTEGER,PARAMETER:: MaxCpn= 60  ! dimension of temporary vCpn
  !
  TYPE(T_Component):: vTmp(MaxCpn)
  CHARACTER(LEN=80):: Msg
  LOGICAL:: Ok
  INTEGER:: N,I
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Components_Alloc"
  !
  DO I=1,SIZE(vTmp)
    CALL Component_Zero(vTmp(I))
  ENDDO

  CALL Components_Read( & !
  & sKeyWord,   & !IN
  & vEle,vSpc,vMixFas, & !IN
  & TdgK,Pbar,  & !INOUT
  & N,          & !OUT
  & SysType,    & !OUT
  & Ok,         & !OUT
  & Msg,        & !OUT
  & vTmp)         !OUT
  !
  IF(.NOT.Ok) CALL Stop_(TRIM(Msg))
  IF(N==0) CALL Stop_("NO COMPONENTS FOUND")

  IF(ALLOCATED(vCpn)) DEALLOCATE(vCpn)
  ALLOCATE(vCpn(1:N))
  vCpn(1:N)= vTmp(1:N)

  IF(iDebug>0) THEN
    DO i=1,N
      WRITE(fTrc,'(A)') vCpn(I)%NamCp
    ENDDO
  ENDIF

  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Components_Alloc"

  RETURN
ENDSUBROUTINE Components_Alloc

SUBROUTINE System_Build_Custom( &
!--
!-- from block sKeyWord, build a (sub)system vCpn of master system vCpnMaster
!--
& sKeyWord,          & !IN
& vEle,vSpc,vMixFas, & !IN
& vCpnMaster,        & !IN
& TdgK,Pbar,         & !INOUT
& vCpn,              & !INOUT
& Ok)                  !OUT
!--
!-- CAVEAT=
!-- the input master system must have all components INERT (or BUFFER) !!WHY ??
!--
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Species,  ONLY: T_Species
  USE M_T_Tpcond,   ONLY: T_TPcond
  USE M_T_MixPhase, ONLY: T_MixPhase
  USE M_T_Component,ONLY: T_Component,Component_Index,Component_Print,CpnMolMinim
  USE M_Component_Read,ONLY: Components_Read
  !
  CHARACTER(LEN=*), INTENT(IN):: sKeyWord
  TYPE(T_Element),  INTENT(IN):: vEle(:)
  TYPE(T_Species),  INTENT(IN):: vSpc(:)
  TYPE(T_MixPhase), INTENT(IN):: vMixFas(:)
  TYPE(T_Component),INTENT(IN):: vCpnMaster(:)
  !
  REAL(dp),         INTENT(INOUT):: TdgK,Pbar
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  LOGICAL,          INTENT(OUT)  :: Ok
  !
  CHARACTER(LEN=7):: SysType
  !
  LOGICAL,          ALLOCATABLE:: vCpnFound(:)
  TYPE(T_Component),ALLOCATABLE:: vCpnNew(:)
  !
  TYPE(T_Component):: C
  INTEGER:: iCp,I,N,nCp
  !!LOGICAL:: Ok
  CHARACTER(LEN=80):: Msg
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< System_Build_Custom"
  !
  !--- default values for component not found in BuildLnk
  vCpn= vCpnMaster  !! <-must be done "outside" ???
  vCpn(:)%Mole=   CpnMolMinim
  WHERE(vCpn(:)%Statut/= "BUFFER") vCpn(:)%Statut= "INERT"
  !---/
  !
  nCp= SIZE(vCpn)
  ALLOCATE(vCpnFound(1:nCp))  ;  vCpnFound=  .FALSE.
  ALLOCATE(vCpnNew(1:nCp))    ;  vCpnNew(:)= vCpn(:)
  !
  CALL Components_Read( & !
  & sKeyWord,   & !IN
  & vEle,vSpc,vMixFas, & !IN
  & TdgK,Pbar,  & !INOUT
  & N,          & !OUT
  & SysType,    & !OUT
  & Ok,         & !OUT
  & Msg,        & !OUT
  & vCpnNew)      !OUT
  !
  !~ IF(.NOT. Ok) CALL Stop_(TRIM(Msg))
  IF(.NOT. Ok) RETURN
  !
  CALL System_TP_Check(TdgK,Pbar,Ok,Msg)
  !~ IF(.NOT. Ok) CALL Stop_(TRIM(Msg))
  IF(.NOT. Ok) RETURN
  !
  !IF(N==0) CALL Stop_("NO COMPONENTS FOUND")
  Ok= N>0
  IF(.NOT. Ok) RETURN !------------------------------------N=0 -> RETURN

  !--------------------------------------components listed in new system
  DO I=1,N
    C= vCpnNew(I)
    iCp= Component_Index(C%NamCp,vCpnMaster)
    IF(iCp==0) THEN
      Ok= .FALSE.
      CALL Stop_("Error in System_Build_Custom: new Cpn not in vCpnMaster")
    ENDIF
    !!print *, "Ele,iCp=",vEle(C%iEle)%Name, iCp  ;  pause
    !-> will have same element order in vCpn as in vCpnMaster
    vCpnFound(iCp)=.TRUE.
    !! IF(C%Statut=="BUFFER") THEN
    !!   IF(vCpnMaster(iCp)%Statut/="BUFFER") &
    !!   & CALL Stop_ &
    !!   & ("A Cpn cannot be BUFFER in subsystem IF it's not also BUFFER in master system !!!")
    !! ENDIF
    !IF(vCpnMaster(iCp)%Statut/="BUFFER")
    vCpn(iCp)= C
    !-> when a Cpn is BUFFER in master system, it is BUFFER in any system
  ENDDO

  !--------------------------------- components not listed in new system
  !old- when an element is not listed, it is "absent": %Mole-CpnMolMinim
  !new- when an element is not listed, it is like in the master system
  DO iCp=1,nCp
    IF(.NOT. vCpnFound(iCp)) THEN
      !
      vCpn(iCp)= vCpnMaster(iCp)
      !
      IF(vCpnMaster(iCp)%Statut/="BUFFER") THEN
        !
        vCpn(iCp)%Statut= "INERT"
        !
        IF(vSpc(vCpn(iCp)%iSpc)%Typ /= "AQU") THEN
          PRINT '(A)', "unlisted components in subsystem are INERT !!"
          PRINT '(2A)',vSpc(vCpn(iCp)%iSpc)%NamSp,"-> species should be aqueous !!!"
          Ok= .FALSE.
        ENDIF
      ENDIF
      !
    ENDIF
  ENDDO
  !
  DEALLOCATE(vCpnNew)
  DEALLOCATE(vCpnFound)
  !
  !--------update global thermo'parameters (vSpc,vMixModel,vMixFas,vFas)
  !--------at (TdgK,Pbar)
  CALL System_TP_Update(TdgK,Pbar)
  !-------------------------------------/update global thermo'parameters
  
  !----------------------------------------------------------------trace
  IF(iDebug>2) THEN
    PRINT '(A)',"< System_Build_Custom "//TRIM(sKeyWord)
    DO I=1,SIZE(vCpn)
      CALL Component_Print(6,vEle,vSpc,vCpn(I))
    ENDDO
    PRINT '(A)',"</ System_Build_Custom"
    CALL Pause_
  ENDIF
  !---------------------------------------------------------------/trace
  
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ System_Build_Custom"
  
ENDSUBROUTINE System_Build_Custom

SUBROUTINE Elements_Alloc_forSystem(vCpn)
!--
!-- rebuild vEle, with only elements pointed to by vCpn(1:N)
!--
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Component,ONLY: T_Component
  USE M_Global_Vars,ONLY: vEle
  !
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  TYPE(T_Element),ALLOCATABLE:: vEleB(:)
  INTEGER:: I,N
  !
  N= SIZE(vCpn)
  !
  ALLOCATE(vEleB(1:SIZE(vEle)))
  vEleB= vEle !-------------------------------save current vEle to vEleB
  !
  DEALLOCATE(vEle) !----------------------------------- destroy old vEle
  ALLOCATE(vEle(1:N))
  vEle(1:N)= vEleB(vCpn(1:N)%iEle) !----------------------build new vEle
  !
  !----------------------------- update fields in component's constructs
  DO i=1,N
    vCpn(i)%iEle= I
    vCpn(i)%vStoikCp(0)= 1 !formula divider !!
    vCpn(i)%vStoikCp(:)= 0
    vCpn(i)%vStoikCp(i)= 1
  ENDDO
  !--------------------------------------------------------------------/
  !
  DEALLOCATE(vEleB)
  !
ENDSUBROUTINE Elements_Alloc_forSystem

SUBROUTINE Redox_Calc(vSpc,vEle,vCpn)
!--
!-- first, from the names of the primary species,
!-- compute the valencies of "potentially redox elements" (Fe,Cl,...)
!--
  USE M_T_Component, ONLY: T_Component
  USE M_T_Element,   ONLY: T_Element,Element_Index,Formula_Read
  USE M_T_Species,   ONLY: T_Species,Species_Stoikio
  !
  TYPE(T_Species),   INTENT(IN)   :: vSpc(:)
  TYPE(T_Element),   INTENT(INOUT):: vEle(:)
  TYPE(T_Component), INTENT(INOUT):: vCpn(:)
  !
  TYPE(T_Species):: S
  INTEGER:: ieOx, iCp, nCp, Z !, Zsp
  LOGICAL:: fOk
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Redox_Calc"
  !
  ieOx= Element_Index("OX_",vEle)
  !
  nCp= SIZE(vEle)
  !
  DO iCp=1,nCp
    !
    ! restrict this procedure to elements with no valency
    ! assigned in dtb (Fe,S,Cr,As)
    IF(vEle(iCp)%Z==0) THEN !or Ele%ReDox=="VAR" ???

      S= vSpc(vCpn(iCp)%iSpc)

      CALL Species_Stoikio(vEle,ieOx,S,fOk)
      ! nota: true stoikiometry is vStoik/S%Div, true charge is Zsp/S%Div
      ! IF(abs(vStoik(iCp))/=1) &
      ! & CALL Stop_("coeff' should be +/-1")

      ! IF(fOk) THEN !fOk=FormulaIsOk (as regards elements in run)
      ! calculate reDox state of element in its Prim'Species -> basis valency

        Z= DOT_PRODUCT( S%vStoikio(1:nCp),  vEle(1:nCp)%Z ) &
        &             - S%vStoikio(iCp  ) * vEle(iCp  )%Z
        !
        vEle(iCp)%Z= (S%Z-Z) / S%vStoikio(iCp)
        !
        IF(iDebug>0) WRITE(fTrc,'(2I3,A,A3,A1,3I3)') &
        & S%Z-Z, S%vStoikio(iCp), &
        & " -> Default Valency of ", vEle(iCp)%NamEl, "=", vEle(iCp)%Z
        !
      !ENDIF

    ENDIF
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Redox_Calc"
  !
ENDSUBROUTINE Redox_Calc
  !
SUBROUTINE FormulaTable_Calc(vSpc,T)
!--
!-- compute the formula table elements / species
!--
  USE M_T_Species,  ONLY: T_Species
  !
  TYPE(T_Species),INTENT(IN) :: vSpc(:)
  REAL(dp),       INTENT(OUT):: T(:,:)
  !
  INTEGER:: I, N
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "< FormulaTable_Calc"
  !
  N= SIZE(T,1)
  DO I=1,SIZE(vSpc)
    T(1:N,I)= vSpc(i)%vStoikio(1:N) / REAL(vSpc(I)%vStoikio(0))
    !column I of T = Formula of species I
    !-> T is the "Formula Matrix" a la SmithMissen
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "</ FormulaTable_Calc"
  !
ENDSUBROUTINE FormulaTable_Calc

SUBROUTINE FormulaTable_Sho(vEle,vSpc,T)
  USE M_IoTools,  ONLY: GetUnit
  USE M_Files,    ONLY: DirOut,Files_Index_Write
  USE M_T_Element,ONLY: T_Element
  USE M_T_Species,ONLY: T_Species
  !
  TYPE(T_Element),  DIMENSION(:),  INTENT(IN):: vEle
  TYPE(T_Species),  DIMENSION(:),  INTENT(IN):: vSpc
  REAL(dp),         DIMENSION(:,:),INTENT(IN):: T
  !
  INTEGER:: F,I,iEl
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_stoik.log")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_stoik.log",&
  & "STOIKIO: formula matrix tStoikio, transposed")
  !
  DO iEl=1,SIZE(vEle)
    WRITE(F,'(A3,A1)',ADVANCE='NO') vEle(iEl)%NamEl, T_
  ENDDO
  WRITE(F,'(3(A,A1))') "Chg",T_,"species",T_
  !
  DO I=1,SIZE(vSpc)
    DO iEl=1,SIZE(vEle)
      WRITE(F,'(F7.2,A1)',ADVANCE='NO') T(iEl,I), T_
    ENDDO
    WRITE(F,'(I3,A1,A12,A1)') vSpc(I)%Z, T_, TRIM(vSpc(I)%NamSp),T_
  ENDDO
  !
  CLOSE(F)
ENDSUBROUTINE FormulaTable_Sho

ENDMODULE M_System_Tools
