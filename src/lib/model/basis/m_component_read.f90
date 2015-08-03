MODULE M_Component_Read
  USE M_Kinds
  USE M_T_Component,ONLY: T_Component
  USE M_Trace
  !
  IMPLICIT NONE

  PRIVATE

  PUBLIC:: Components_Read
  !
CONTAINS

SUBROUTINE Components_Read( & !
!--
!-- scan the SYSTEM block -> build vCpn, read (T,P) conditions --
!--
& sKeyWord,  & !IN
& vEle,      & !IN
& vSpc,      & !IN
& vMixFas,   & !IN
& TdgK,Pbar, & !INOUT
& N,         & !OUT
& SysType,   & !OUT
& Ok,        & !OUT
& Msg,       & !OUT
& vCpn)        !INOUT
  !---------------------------------------------------------------------
  USE M_IoTools
  USE M_Files,      ONLY: NamFInn
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_Numeric_Const,ONLY: Ln10
  USE M_T_Element,  ONLY: T_Element,Element_Index
  USE M_T_Species,  ONLY: T_Species,Species_Index,Species_Rename
  USE M_T_MixPhase, ONLY: T_MixPhase,MixPhase_Index
  USE M_T_Component,ONLY: Component_Zero,CpnMolMinim
  !
  USE M_Dtb_Vars,ONLY: DtbLogK_vTPCond
  !
  USE M_System_Vars,ONLY: BufferIsExtern
  !---------------------------------------------------------------------
  CHARACTER(LEN=*), INTENT(IN):: sKeyWord
  TYPE(T_Element),  INTENT(IN):: vEle(:)
  TYPE(T_Species),  INTENT(IN):: vSpc(:)
  TYPE(T_MixPhase), INTENT(IN):: vMixFas(:)
  !
  REAL(dp),         INTENT(INOUT):: TdgK,Pbar
  INTEGER,          INTENT(OUT)  :: N
  CHARACTER(LEN=7), INTENT(OUT)  :: SysType
  LOGICAL,          INTENT(OUT)  :: Ok
  CHARACTER(*),     INTENT(OUT)  :: Msg
  TYPE(T_Component),INTENT(OUT)  :: vCpn(:)
  !---------------------------------------------------------------------
  !
  !---------------------------------------------------------development
  LOGICAL,PARAMETER:: MoleIsTrue= .FALSE.  !.TRUE. !
  != MoleIsTrue => when Statut is MOLE or GRAM
  != species stoichiometry is considered for
  != computation of element mole numbers
  !---------------------------------------------------------/development
  TYPE(T_Component) :: Cpn
  CHARACTER(LEN=255):: L,sListElem
  CHARACTER(LEN=80) :: W0,W,V1,V2,Wnum,WType
  CHARACTER(LEN=7)  :: Statut
  LOGICAL :: EoL
  LOGICAL :: bForceH2O =.false.
  INTEGER :: f,ios,iTP
  INTEGER :: I,iEl,iSp,iW,iH_,iO_,nEl
  REAL(dp):: X1,Eps
  REAL(dp),ALLOCATABLE:: vX(:)
  INTEGER, ALLOCATABLE:: tStoikCp(:,:)
  !---------------------------------------------------------------------
  
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Read_Input_System"
  !
  Ok= .TRUE.
  Msg= "Ok"
  !
  CALL GetUnit(f)
  CALL OPENFILE(f,FILE=TRIM(NamFInn))
  !
  iTP= 0
  SysType="Z"
  sListElem="" !used to prevent reading same element more than once
  !
  N=   0
  iH_= 0
  iO_= 0
  nEl= SIZE(vEle)
  !
  ALLOCATE(tStoikCp(SIZE(vCpn),0:nEl+1))
  tStoikCp(:,:)= 0 ! stoikiometry
  tStoikCp(:,0)= 1 ! divisor
  !
  DoFile: DO
    !
    READ(f,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W0,EoL)
    IF(W0(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W0,EoL)
    !
    IF(W0=="ENDINPUT") EXIT DoFile
    !!SELECT CASE(W0)
    !!  CASE("SYSTEM","SYSTEM.AQUEOUS"); SysType="AQUEOUS"
    !!  CASE("SYSTEM.MOMAS");            SysType="MOMAS"
    !!END SELECT
    IF(TRIM(W0)==TRIM(sKeyWord)) THEN
      !
      !-- possibility to have several SYSTEM blocks, each with a Name,
      !-- -> to define several systems, make them react with another, etc.
      IF ( TRIM(sKeyWord)=="SYSTEM"         .OR. &
      &    TRIM(sKeyWord)=="SYSTEM.AQUEOUS" .OR. &
      &    TRIM(sKeyWord)=="SYSTEM.BOX"     .OR. &
      &    TRIM(sKeyWord)=="SYSTEM.MIX"     .OR. &
      &    TRIM(sKeyWord)=="SYSTEM.INJECT")  &
      & SysType="AQUEOUS"
      !
      CALL Component_Zero(Cpn)
      !
      Statut= "INERT" !default value
      !
      IF(SysType=="AQUEOUS") THEN
        !
        IF(Element_Index("O__",vEle)==0) THEN
          Ok= .FALSE.
          Msg= "O__ not found among elements !!!"
          RETURN
        ENDIF
        !
        IF(Element_Index("H__",vEle)==0) THEN
          Ok= .FALSE.
          Msg= "H__ not found among elements !!!"
          RETURN
        ENDIF
        !
        IF(Species_Index("H2O",vSpc)==0) THEN
          Ok= .FALSE.
          Msg= "H2O not found among species  !!!"
          RETURN
        ENDIF
        !
        iW= Species_Index("H2O",vSpc)
        !
        IF(bForceH2O) THEN
          !---------------make water (solvent) as Cpn Nr1
          !---------------------------ele-"O__",spc-'H2O"
          iO_= 1
          Cpn%iEle=     Element_Index("O__",vEle)
          Cpn%NamCp=    vEle(Cpn%iEle)%NamEl
          Cpn%iSpc=     Species_Index("H2O",vSpc)
          Cpn%Mole=     One /vSpc(iW)%WeitKg
          Cpn%LnAct=    Zero
          Cpn%Statut=   "INERT"
          Cpn%namSol=   "Z"
          !
          sListElem=    "O__"
          N=N+1
          !
          vCpn(N)= Cpn
          !vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
          tStoikCp(N,Cpn%iEle)= 1
          !
        ENDIF
        !
      ENDIF
      !
      DoBlock: DO
        !
        READ(f,'(A)',IOSTAT=ios) L
        IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL) !-> READ 1st word -> W
        IF(W(1:1)=='!') CYCLE DoBlock !-> skip comment lines
        CALL AppendToEnd(L,W,EoL)
        !
        IF(    W=="ENDINPUT" &
        & .OR. W=="END"      &
        & .OR. W=="END"//TRIM(sKeyWord)) EXIT DoFile
        !
        !------------------------------------read Temperature / Pressure
        !IF(LEN_TRIM(W)>3) THEN
        IF(IsKeyword("_TP.POINT_TPPOINT_TDGK_TDGC_PBAR_PMPA_TDEGK_TDEGC_",W)) &
        & THEN
          !
          CALL LinToWrd(L,Wnum,EoL)
          !-> !!!§§§ add something to check input !!
          SELECT CASE(TRIM(W))
          !
          CASE("TP.POINT","TPPOINT")
            ! CALL Warning_(TRIM(W)//" is Obsolete, better give TdgC, Pbar !!!")
            !
            CALL WrdToInt(Wnum,iTP)
            !
            IF(iTP>0) THEN
              !! CALL Warning_ &
              !! & ("keyword TP.POINT soon Obsolete !! better give TdgC, Pbar in SYSTEM")
              !
              IF(iTP>SIZE(DtbLogK_vTPCond)) THEN
                Ok= .FALSE.
                Msg= "TP.POINT outside range"
                RETURN
              ENDIF
              !
              TdgK= DtbLogK_vTPCond(iTP)%TdgC +T_CK
              Pbar= DtbLogK_vTPCond(iTP)%Pbar
              !
            ENDIF
          !
          CASE("TDGK","TDGC","TDEGK","TDEGC","PBAR","PMPA")
            CALL WrdToReal(Wnum,X1)
            !
            SELECT CASE(TRIM(W))
            !
            CASE("TDGK");  TdgK= X1
            CASE("TDGC");  TdgK= X1 + T_CK
            !
            CASE("PBAR");  Pbar= X1
            CASE("PMPA");  Pbar= X1 / 10.0D0
            !
            CASE("TDEGK")
              TdgK= X1
              CALL Warning_("Depreciated Keyword TDEGK")
            CASE("TDEGC")
              TdgK= X1 + T_CK
              CALL Warning_("Depreciated Keyword TDEGC")
            !
            ENDSELECT
          !
          CASE DEFAULT
            Ok= .FALSE.
            Msg= "Unknown (or Obsolete) Keyword: "//TRIM(W)
            RETURN
          !
          END SELECT
          !
          CYCLE DoBlock
          !
        ENDIF
        !-----------------------------------/read Temperature / Pressure
        !
        !----------------------------------------------read element name
        IF(.NOT. IsElement(vEle,W)) THEN
          Ok= .FALSE.
          Msg=           "Element "//TRIM(W)//" Not in Element base !!!"
          RETURN !------------------------------------------------RETURN
        ENDIF
        !
        CALL Str_Append(W,3)
        !
        !--position of current element in vEle
        iEl= Element_Index(TRIM(W),vEle)
        !
        IF(INDEX(sListElem,TRIM(W))>0) THEN
          Ok= .FALSE.
          Msg=              "Component "//TRIM(W)//" Already listed !!!"
          RETURN !------------------------------------------------RETURN
        ENDIF
        !
        N=N+1
        IF(N>SIZE(vCpn)) THEN
          Ok= .FALSE.
          Msg=                                     "Too many components"
          RETURN !------------------------------------------------RETURN
        ENDIF
        !
        Cpn%NamCp= TRIM(W)
        Cpn%iEle= iEl
        !Cpn%NamCp= TRIM(vEle(Cpn%iEle)%NamEl)
        !Cpn%vStoikCp(iEl)= 1
        tStoikCp(N,iEl)= 1
        !
        IF(W=="H__") iH_= N
        IF(W=="O__") iO_= N
        !
        sListElem=TRIM(sListElem)//TRIM(W)
        !-> prevent reading same component in subsequent loops
        !
        !---------------------------------------------/read element name
        !
        !------------------------------------------read component status
        CALL LinToWrd(L,WType,EoL)
        !
        IF(.NOT. IsKeyword &
        & ("_INERT_MOLE_GRAM_MOBILE_ACTIVITY_PK_BUFFER_BALANCE_",WType)) THEN
          Ok= .FALSE.
          Msg=         TRIM(WType)//" -> unrecognized species status !!"
          RETURN !------------------------------------------------RETURN
        ENDIF
        !
        SELECT CASE(TRIM(WType))
        CASE("INERT","MOLE","GRAM")              ;  Statut="INERT"
        CASE("MOBILE","ACTIVITY","PK")           ;  Statut="MOBILE"
        CASE("BUFFER")                           ;  Statut="BUFFER"
        CASE("BALANCE")                          ;  Statut="BALANCE"
        !! CASE("STOIKIO")                ;  Statut="STOIKIO"
        !! CASE DEFAULT
        !!   Ok= .FALSE.
        !!   Msg=      TRIM(WType)//" -> unrecognized species status !!"
        !!   RETURN !---------------------------------------------RETURN
        END SELECT
        !
        !-----------------------------------------------------for Oxygen
        IF(N==iO_ .AND. Statut/="INERT") THEN
          Ok= .FALSE.
          Msg= "Oxygen should be INERT !!!"
          RETURN
        ENDIF
        !----------------------------------------------------/for Oxygen
        !
        !-- assign component mobility status
        Cpn%Statut= TRIM(Statut)
        !-----------------------------------------/read component status
        !
        !----------------------------------------------read species name
        IF(EoL) THEN
          Ok= .FALSE.
          Msg=                   "Species lacking for "//TRIM(Cpn%NamCp)
          RETURN
        ENDIF
        CALL LinToWrd(L,V1,EoL)!= species name
        !
        CALL Species_Rename(V1)
        !
        iSp=Species_Index(V1,vSpc)
        !
        IF(iSp==0) THEN
          Ok= .FALSE.
          Msg=              TRIM(V1)//" <-Species NOT in database ???!!"
          RETURN
        ELSE
          IF(Statut=="INERT"  .AND. vSpc(iSp)%Typ/="AQU") THEN
            Msg=  TRIM(V1)//" <- an INERT SPECIES should be AQUEOUS !!!"
            RETURN
          ENDIF
          !<BUFFER MODIFIED
          IF(BufferIsExtern   .AND. &
          &  Statut=="BUFFER" .AND. &
          &  vSpc(iSp)%Typ=="AQU") THEN
            Msg= TRIM(V1)//" <-a BUFFER SPECIES should NOT be AQUEOUS !!!"
            RETURN
          ENDIF
          !</BUFFER MODIFIED
        ENDIF
        !
        !------------------------------------------------- for Oxygen --
        IF(N==iO_ .AND. TRIM(V1)/="H2O") THEN
          Ok= .FALSE.
          Msg=                   "for Oxygen, species should be H2O !!!"
          RETURN
        ENDIF
        !------------------------------------------------/ for Oxygen --
        !
        Cpn%iSpc=iSp != index of "associated" species
        !
        IF(Cpn%Statut=="INERT") THEN
          IF(vSpc(iSp)%vStoikio(iEl)==0) THEN
            Ok= .FALSE.
            Msg= &
            & "Element "//TRIM(vEle(iEl)%NamEl)// &
            & " is not contained in species "//TRIM(vSpc(iSp)%NamSp)
            RETURN
          ENDIF
        ENDIF
        !
        !<new, 200911>
        IF(TRIM(WType)=="MOLE" .OR. TRIM(WType)=="GRAM") THEN
        !--------------- write species'stoichio to component'stoichio --
          !Cpn%vStoikCp(0:nEl+1)= vSpc(iSp)%vStoikio(0:nEl+1)
          tStoikCp(N,0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
        ENDIF
        !</new>
        !
        !---------------------------------------------/read species name
        !
        !------------------------- read numeric data (mole, pk, etc.) --
        X1=Zero
        Cpn%namSol= "Z" !TRIM(vSpc(iSp)%Typ) !default value
        !
        IF(Statut/="BALANCE") THEN
          !
          IF(EoL) THEN
          !----------------------if no numeric data, take default values
            !
            !------------ for mobile pure species, default is saturation
            IF(Statut=="MOBILE" .OR. Statut=="BUFFER") X1=Zero
            !
            IF(Statut=="INERT") X1= CpnMolMinim
            !
          ELSE
            !
            CALL LinToWrd(L,V2,EoL)
            !
            IF(TRIM(V2)/="SOLUTION" .AND. TRIM(V2)/="MIXTURE") THEN
              ! read real value -> total amount, colog10(activity), etc
              CALL WrdToReal(V2,X1)
            ELSE
              !
              IF(TRIM(V2)=="SOLUTION") &
              & CALL Warning_(TRIM(V2)//" soon Obsolete, better use MIXTURE !!!")
              !
              IF(Cpn%Statut=="INERT") THEN
                Ok= .FALSE.
                Msg= "SOLUTION/MIXTURE is NOT VALID KEYWORD for INERT components !!!"
                RETURN
              ENDIF
              !-------------------------read name of "controlling phase"
              CALL LinToWrd(L,V2,EoL)
              IF(MixPhase_Index(vMixFas,V2)==0) THEN
                Ok= .FALSE.
                Msg=          TRIM(V2)//" <-this phase is unknown !!!!!"
                RETURN
              ENDIF
              !
              Cpn%namSol=TRIM(V2)
              !!print *, TRIM(V2)  ;  pause
              !
            ENDIF
            !
          ENDIF
          !
        ENDIF !!IF(Statut/="BALANCE")
        !---------------------------------------------/read numeric data
        !
        !-------------------------------------------process numeric data
        SELECT CASE(TRIM(Statut))
        !
        CASE("INERT") !---------------------------------------case INERT
          !
          Cpn%Mole=X1  != tot.amount of element
          !
          Cpn%Factor= 1.D0
          !---------------------------------------------unit conversions
          IF (TRIM(WType)=="GRAM") THEN
            Cpn%Factor= vSpc(iSp)%WeitKg *1.0D3
          ELSE
            !---------------------------------------------------OBSOLETE
            !-------------------------------maintained for compatibility
            IF(.NOT. EoL) THEN
              !------------------------- read unit: MG/KG__, UG/KG__,...
              CALL LinToWrd(L,V2,EoL)
              SELECT CASE(TRIM(V2))
              CASE("MG/KG","PPM")
                Cpn%Factor= vSpc(iSp)%WeitKg *1.0D6
              CASE("UG/KG","PPB")
                Cpn%Factor= vSpc(iSp)%WeitKg *1.0D9
              CASE DEFAULT
              !------------------------------------unknown keyword, STOP
                Ok= .FALSE.
                Msg= TRIM(V2)//" = unknown unit..."
                RETURN
              END SELECT
            ENDIF
            !--------------------------------------------------/OBSOLETE
          ENDIF
          !--------------------------------------------/unit conversions
          Cpn%Mole= Cpn%Mole /Cpn%Factor
          !
        CASE("MOBILE","BUFFER") !----------------------------case MOBILE
          !
          Cpn%LnAct= -X1 *Ln10 !from colog10 to ln
          !
          !---------------------------------------------unit conversions
          SELECT CASE(TRIM(WType))
          !
          CASE("PK")
            Cpn%LnAct= -X1 *Ln10
          !
          CASE("ACTIVITY")
            Eps=EPSILON(X1)
            IF(ABS(X1)<Eps) THEN
              Ok= .FALSE.
              Msg= "Activity cannot be Zero !!!"
              RETURN
            ENDIF
            Cpn%LnAct= LOG(X1)
          !
          END SELECT
          !--------------------------------------------/unit conversions
          !
          !<OBSOLETE>
          !!<keyword BUFFER is now placed differently>
          !!IF(.NOT. EoL) THEN
          !!  CALL LinToWrd(L,V2,EoL)
          !!  IF(TRIM(V2)=="BUFFER") Statut="BUFFER"
          !!  !!! BUFFER MODIFIED
          !!  IF(BufferIsExtern .AND. Statut=="BUFFER" .AND. vSpc(iSp)%Typ=="AQU") &
          !!  & CALL Stop_(TRIM(V1)//" <-a BUFFER SPECIES should NOT be AQUEOUS !!!!!")
          !!  !!! BUFFER MODIFIED END
          !!ENDIF
          !</OBSOLETE>
          !
        END SELECT
        !--------------------------------------/ process numeric data --
        !
        IF (iDebug>2) WRITE(fTrc,'(3(A,A1),2(G15.6,A1))') &
        & Cpn%Statut,          T_, &
        & vEle(Cpn%iEle)%NamEl,T_, &
        & vSpc(Cpn%iSpc)%NamSp,T_, &
        & Cpn%Mole,            T_, &
        & Cpn%LnAct,           T_
        !
        vCpn(N)= Cpn
        !
      ENDDO DoBlock
    ENDIF !
  ENDDO DoFile
  CALL CLOSEFILE(f)
  !
  IF (iDebug>2) WRITE(fTrc,'(A)') &
  & "========================================================================="
  !
  IF(N==0) RETURN
  !
  !~ IF(Check_TP) THEN
    !~ IF(DtbFormat=="LOGK") THEN
      !~ iTP_= TPcond_IndexT(TdgK,vTPCond)
      !~ IF(iTP_==0) THEN
        !~ PRINT '(A,G15.6)',"TdgC=", TdgK -T_CK
        !~ CALL Stop_("Temperature not found in TP-series")
      !~ ELSE
        !~ TdgK= vTPCond(iTP_)%TdgC +T_CK
        !~ Pbar= vTPCond(iTP_)%Pbar
      !~ ENDIF
    !~ ELSE
      !~ iTP_= TPcond_IndexTP(TdgK,Pbar,vTPCond)
      !~ IF(iDebug==4) CALL Pause_("TPcond_IndexTP= "//TRIM(FIntToStr3(iTP_)))
      !~ !
      !~ IF(TdgK <Zero) TdgK= Zero
      !~ !
      !~ ! prevent pressure to be below vapor saturation for pure H2O
      !~ IF (TdgK<=647.25D0) THEN; CALL fluid_h2o_psat(TdgK,PSatBar)
      !~ ELSE                    ; PSatBar= 220.55D0
      !~ ENDIF
      !~ IF(Pbar<PSatBar) Pbar=PSatBar
    !~ ENDIF
    !~ IF(iDebug==4) PRINT '(/,A,2G15.6,/)',"TdgC,Pbar=",TdgK-T_CK,Pbar
  !~ ENDIF
  !
  IF(SysType=="AQUEOUS") THEN
    !----------------------- add O__/H2O IF not found in species list --
    IF(iO_==0) THEN
      !
      N= N+1
      IF(N>SIZE(vCpn)) THEN
        Ok= .FALSE.
        Msg= "Too many components"
        RETURN
      ENDIF
      !
      iO_= N
      !
      Cpn%iEle=   Element_Index("O__",vEle)
      Cpn%NamCp=  vEle(Cpn%iEle)%NamEl
      Cpn%iSpc=   Species_Index("H2O",vSpc)
      Cpn%Mole=   One /vSpc(iW)%WeitKg
      Cpn%LnAct=  Zero
      Cpn%Statut= "INERT"
      Cpn%namSol= "Z"
      !
      vCpn(N)= Cpn
      !Cpn%vStoikCp(0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
      !tStoikCp(N,0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
      tStoikCp(N,Cpn%iEle)= 1
      !
      sListElem= TRIM(sListElem)//"O__"
      !
    ENDIF
    !----------------------------------------------------/ add O__/H2O--
    !
    !------------------------ add H__/H+ IF not found in species list --
    IF(iH_==0) THEN !
      N=   N+1
      IF(N>SIZE(vCpn)) THEN
        Ok= .FALSE.
        Msg= "Too many components"
        RETURN
      ENDIF
      !
      iH_= N
      !
      Cpn%iEle=   Element_Index("H__",vEle)
      Cpn%NamCp=  vEle(Cpn%iEle)%NamEl
      Cpn%iSpc=   Species_Index("H+", vSpc)
      !!Cpn%Mole=   2.0D0 /vSpc(iW)%WeitKg
      Cpn%LnAct=  Zero
      Cpn%Statut= "INERT"
      Cpn%namSol= "Z"
      !
      vCpn(N)= Cpn
      tStoikCp(N,Cpn%iEle)= 1
      !Cpn%vStoikCp(0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
      !tStoikCp(N,0:nEl+1)= vSpc(Cpn%iSpc)%vStoikio(0:nEl+1)
      !
      sListElem=TRIM(sListElem)//"H__"
      !
    ENDIF
    !----------------------------------------------------/ add H__/H+--
  ENDIF !IF(SysType=="AQUEOUS")
  !
  !--------------------------------- compute mole numbers of elements --
  IF(MoleIsTrue) THEN
    ALLOCATE(vX(N))
    vX(:)= Zero
    !
    DO iEl=1,N
      DO I=1,N
        IF(vCpn(I)%Statut=="INERT") &
        & vX(iEl)= vX(iEl) +vCpn(I)%Mole *tStoikCp(I,vCpn(iEl)%iEle)
      ENDDO
    ENDDO
    !
    !--- trace
    IF(iDebug==4) THEN
      PRINT *,"== Components_Read =="
      DO iEl=1,N
        PRINT *,vCpn(iEl)%NamCp,vCpn(iEl)%iEle
        DO I=1,N
          PRINT *, &
          & "        ",      &
          & vCpn(I)%NamCp,    &
          & TRIM(vSpc(vCpn(I)%iSpc)%NamSp), &
          & tStoikCp(I,vCpn(iEl)%iEle)
          !!& vCpn(I)%vStoikCp(vCpn(iEl)%iEle)
        ENDDO
        PRINT '(A,G15.6)',"=======================================tot=",vX(iEl)
      ENDDO
      CALL Pause_
    ENDIF !==</ trace
    !
    vCpn(1:N)%Mole= vX(1:N)
    !
    DEALLOCATE(vX)
  ENDIF
  !--------------------------------/ compute mole numbers of elements --
  !
  DEALLOCATE(tStoikCp)
  !
  IF(iH_/=0) vCpn(iH_)%Mole= 2.0D0*vCpn(iO_)%Mole
  !
  !----------------------------------------- place Oxygen as species 1--
  IF(iO_/=0) THEN
    IF(iO_/=1) THEN
      Cpn=       vCpn(1)
      vCpn(1)=   vCpn(iO_)
      vCpn(iO_)= Cpn
    ENDIF
  ENDIF
  !----------------------------------------/ place Oxygen as species 1--
  !
  IF (iDebug>0) THEN
    DO I=1,N
      Cpn= vCpn(I)
      WRITE(fTrc,'(3(A,A1),2(G15.6,A1))') &
      & Cpn%Statut,          T_, &
      & vEle(Cpn%iEle)%NamEl,T_, &
      & vSpc(Cpn%iSpc)%NamSp,T_, &
      & Cpn%Mole,            T_, &
      & Cpn%LnAct,           T_
    ENDDO
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A10,A)') "LISTELEM= ",TRIM(sListElem)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ ReadInput_System"
  !
ENDSUBROUTINE Components_Read

LOGICAL FUNCTION IsSpecies(vSpc,W)
  USE M_T_Species,ONLY: T_Species,Species_Index,Species_Rename
  TYPE(T_Species), INTENT(IN):: vSpc(:)
  CHARACTER(LEN=*),INTENT(IN):: W
  !
  IsSpecies= (Species_Index(TRIM(W),vSpc)/=0)
  !
END FUNCTION IsSpecies

LOGICAL FUNCTION IsElement(vEle,W)
  USE M_T_Element,ONLY: T_Element,Element_Index
  USE M_IoTools,  ONLY: Str_Append
  TYPE(T_Element), INTENT(IN):: vEle(:)
  CHARACTER(LEN=*),INTENT(IN):: W
  !
  CHARACTER(LEN=3):: W2
  !
  W2= W(1:3)
  CALL Str_Append(W2,3)
  IsElement= (Element_Index(TRIM(W2),vEle)>0)
  !
END FUNCTION IsElement

LOGICAL FUNCTION IsKeyword(sList,W)
  CHARACTER(LEN=*),INTENT(IN):: sList
  CHARACTER(LEN=*),INTENT(IN):: W
  !
  IsKeyword= (INDEX(TRIM(sList),"_"//TRIM(W)//"_")>0)
  !
END FUNCTION IsKeyword

ENDMODULE M_Component_Read
