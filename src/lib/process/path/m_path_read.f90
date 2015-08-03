MODULE M_Path_Read
!--
!-- read parameters for path calculations --
!--
!
  USE M_Kinds
  USE M_Trace
  USE M_T_MixPhase,ONLY: T_MixPhase
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Path_ReadParam
  PUBLIC:: Path_ReadParam_new
  PUBLIC:: Path_ReadMode
  ! 
CONTAINS

SUBROUTINE Path_ReadMode( &
& NamFInn, &
& PathMode,Ok,Msg)
!--
!-- read path mode --
!--
  USE M_IOTools
  !
  CHARACTER(LEN=*),INTENT(IN) :: NamFInn
  CHARACTER(LEN=*),INTENT(OUT):: PathMode
  LOGICAL,         INTENT(OUT):: Ok
  CHARACTER(LEN=*),INTENT(OUT):: Msg
  !
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W
  LOGICAL:: EoL
  INTEGER:: f, ios
  !
  Ok=.FALSE.
  Msg= ""
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  !
  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L  ; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!')   CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    IF(W=="PATH") THEN
      !
      IF(EoL) THEN
        Ok=  .false.
        Msg= "Error in PATH: Path type undefined"
        RETURN
      ENDIF
      ! 
      !------------------------------------------- read the path mode --
      CALL LinToWrd(L,W,EoL)
      !
      Ok=.TRUE.
      !
      SELECT CASE(TRIM(W))
      
      CASE("MIX")         ; PathMode= "MIX"
      CASE("ADD")         ; PathMode= "ADD"
      CASE("ADDALL")      ; PathMode= "ADA"
      CASE("CHANGE")      ; PathMode= "CHG"
      CASE("GRID")        ; PathMode= "GRD"
      CASE("EVAPORATE")   ; PathMode= "EVP"
      CASE("LOGK")        ; PathMode= "LGK"
      
      CASE DEFAULT
        Ok=  .false.
        Msg= "Error in PATH: Path mode should be CHANGE | ADD | MIX | LOGK"
        RETURN
      
      ENDSELECT
      !------------------------------------------/ read the path mode --
      !
    ENDIF !Cod=="PATH"
    !
  ENDDO DoFile
  !
  CLOSE(f)
  !
  IF(.not. Ok) Msg= "PATH block not found !!"
  !
  RETURN
ENDSUBROUTINE Path_ReadMode

SUBROUTINE Path_ReadParam( &
& NamFInn,   &  !IN
& PathMode,  &  !IN
& vCpn,      &  !IN
& TdgK,Pbar, &  !IN
& Ok,Msg)       !OUT
!--
!-- read the path parameters --
!--
  USE M_IOTools
  USE M_Dtb_Const,   ONLY: T_CK
  USE M_System_Tools,ONLY: System_TP_Update
  USE M_T_Element,   ONLY: Formula_Read,Element_Index
  USE M_T_Species,   ONLY: Species_Index
  USE M_T_Component, ONLY: T_Component,Component_Index
  USE M_Basis,       ONLY: Basis_Change
  USE M_T_MixPhase,  ONLY: MixPhase_Index
  !
  USE M_Global_Vars, ONLY: vEle,vSpc,vFas,vMixFas,vMixModel
  USE M_Basis_Vars,  ONLY: tAlfFs
  USE M_Path_Vars,   ONLY: vLPath,vTPpath,tPathData,DimPath
  USE M_Path_Vars,   ONLY: vPhasBegin,vPhasFinal
  USE M_Path_Vars,   ONLY: iLogK,vPathLogK
  !
  CHARACTER(LEN=*), INTENT(IN)   :: NamFInn
  CHARACTER(LEN=*), INTENT(IN)   :: PathMode
  TYPE(T_Component),INTENT(IN)   :: vCpn(:)
  REAL(dp),         INTENT(IN)   :: TdgK,Pbar
  LOGICAL,          INTENT(OUT)  :: Ok
  CHARACTER(LEN=*), INTENT(OUT)  :: Msg
  !
  INTEGER,PARAMETER:: DimMax= 255
  !
  CHARACTER(LEN=1024):: L
  CHARACTER(LEN=80) :: W,W1
  LOGICAL :: EoL, fOk, OkMix, PathAdd, BuildTable, TestMixture
  INTEGER :: ios, ieOx, Z, Zxs, nDiv, iMix
  INTEGER :: I, J, K, f, iCp, nCp, nStep
  REAL(dp):: rBegin,rFinal,rRatio,rDelta,AddedTot,X
  INTEGER,DIMENSION(1:SIZE(vCpn))::vStoik !for PathAdd
  !
  REAL(dp):: vTmp(DimMax)
  REAL(dp),ALLOCATABLE:: tTmp(:,:)
  INTEGER, ALLOCATABLE:: vNstep(:)
  INTEGER, ALLOCATABLE:: vStAdd(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Path_Read"
  !
  !IF(iDebug==4) PRINT '(A)',"Path_Read"
  !
  Ok=.FALSE.
  Msg= ""
  TestMixture= .FALSE.
  !
  nCp=SIZE(vCpn)
  !
  ALLOCATE(tTmp(1:nCp+2,1:DimMax))
  ! 1:nCp = components
  ! nCp+1 = TdgC
  ! nCp+2 = Pbar
  tTmp(:,:)= Zero
  tTmp(nCp+1,:)= TdgK -T_CK
  tTmp(nCp+2,:)= Pbar
  !
  ALLOCATE(vNstep(nCp+2))  ;  vNstep(:)= 0
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  !
  IF(iDebug>2) THEN
    DO I=1,SIZE(vCpn)
      PRINT '(I3,1X,3A)', &
      & I,vCpn(I)%NamCp," -> STATUT=",vCpn(I)%Statut
    ENDDO
    CALL Pause_
  ENDIF
  !
  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile !---------------------------------end of file
    !
    CALL LinToWrd(L,W,EoL)
    !
    IF(W(1:1)=='!') CYCLE DoFile !--------------------skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    BuildTable= .false.
    !
    IF(W=="PATH") THEN 
      !
      Ok=.TRUE.
      !
      !--------------------allocate variables according to the Path mode
      SELECT CASE(PathMode) 
      CASE("ADD","ADA")
        ALLOCATE(vStAdd(0:nCp))  ; vStAdd=0
        !
      CASE("CHG","EVP")
        ALLOCATE(vLPath(1:nCp+2))  ;  vLPath=.FALSE.
        !
        ALLOCATE(vPhasBegin(1:nCp))  ;  vPhasBegin=0
        ALLOCATE(vPhasFinal(1:nCp))  ;  vPhasFinal=0
        !
      END SELECT
      !-------------------------------------------------------/ allocate
      !
      !---------------------------------------- read the path parameters
      DoBlock: DO
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!')   CYCLE DoBlock !skip comment lines
        CALL AppendToEnd(L,W,EoL)
        !
        IF(W=="ENDINPUT") EXIT  DoFile
        IF(W=="END" .OR. W=="ENDPATH") EXIT  DoFile
        !-> READ ONLY the first PATH block found !!
        !
        SELECT CASE(TRIM(PathMode))
        !-------------------------------------------- CASE "LGK" (-LOGK)
        CASE("LGK")
          
          iLogK= Species_Index(TRIM(W),vSpc)
          IF(iLogK==0) THEN
            Ok=  .false.
            Msg=      "In PATH LOGK, "//TRIM(W)//"= unknown species !!!"
            RETURN !----------------------------------------------------
          ENDIF
          
          !--------------------------------------------------read values
          rBegin=Zero; rFinal=Zero; rDelta=Zero
          DO
            IF(EoL) EXIT
            
            CALL LinToWrd(L,W1,EoL)
            
            SELECT CASE(W1)
            CASE("INITIAL")
              CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rBegin)
            CASE("FINAL")
              CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rFinal)
            CASE("DELTA")
              CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rDelta)
            CASE DEFAULT
              Ok=  .false.
              Msg= "In PATH LOGK, "//TRIM(W1)//"= invalid keyword !!!"
              RETURN !--------------------------------------------RETURN
            END SELECT
          ENDDO
          
          IF(ABS(rDelta)<Epsilon(rDelta)) THEN
            Ok=  .false.
            Msg=                       "In PATH LOGK, Delta not defined"
            RETURN !----------------------------------------------RETURN
          ENDIF
          
          ! IF(ABS(rFinal-rBegin) < 0.1) rFinal= rBegin +One
          
          rDelta= SIGN(rDelta,rFinal-rBegin)
          !-- function SIGN(a,b) returns abs(a)*sign(b) --
          !
          !------------------------------------------------/ read values
          !
          !--------------------------------------------- build vPathLogK
          I= 1
          vTmp(I)= rBegin
          DO
            IF(I>DimMax) EXIT
            IF(ABS(vTmp(I)-rBegin)>ABS(rFinal-rBegin)) EXIT
            I=I+1
            vTmp(I)= vTmp(I-1) + rDelta
            print '(A,I3,G15.6)',"logK=",I,vTmp(i)
          ENDDO
          ! call pause_
          !
          DimPath= I
          ALLOCATE(vPathLogK(1:I))
          vPathLogK(1:I)= vTmp(1:I)
          !-----------------------------------------/ build vPathLogK --
          !
        !----------------------------------------/ CASE "LGK" (-LOGK) --
        !
        !------------------------------------------------- CASE "ADD" --
        CASE("ADD","ADA")
          !
          !------------------------- read stoichio of added material --
          CALL Formula_Read(W,vEle,Z,nDiv,fOk,vStoik)
          !
          IF(Z/=0) THEN
            Ok=  .false.
            Msg=     "In PATH ADD, "//TRIM(W)//"= should be neutral !!!"
            RETURN !----------------------------------------------RETURN
          ENDIF
          !
          !---------------------compute stoikio coeff of redox component
          ieOx= Element_Index("OX_",vEle)
          Zxs= DOT_PRODUCT(vStoik(:),vEle(:)%Z)
          !-- -> oxidation state
          IF(Zxs/=Z) THEN
            IF(ieOx/=0) THEN
              vStoik(ieOx)= Z - Zxs
              fOk= .TRUE.
            ELSE
              fOk= .FALSE.
            ENDIF
          ENDIF
          !-------------------/ compute stoikio coeff of redox component
          !
          IF(.NOT.fOk) THEN
            Ok=  .false.
            Msg= "In PATH ADD, "//TRIM(W)//"= problem in stoikiometry ?"
            RETURN !----------------------------------------------RETURN
          ENDIF
          !
          vStAdd(1:nCp)= vStoik(vCpn(1:nCp)%iEle)
          vStAdd(0)=     nDiv
          !--------------- vStAdd(1:nCp)/vStAdd(0) is the stoikio'vector
          !
          !----------------------------/ read stoichio of added material
          !
          rBegin=  1.D-6
          rFinal=  1.0D0
          rRatio=  1.5D0
          rDelta=  Zero
          !
          nStep= 1
          AddedTot= Zero
          tTmp(1:nCp,1)= vCpn(:)%Mole
          !
          !---------------------------------------------- read values --
          DO
            IF(EoL) EXIT
            !
            CALL LinToWrd(L,W,EoL)
            SELECT CASE(W)
              
              CASE("INITIAL","FINAL","RATIO","DELTA")
                BuildTable= .true.
                SELECT CASE(W)
                  CASE("INITIAL"); CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rBegin)
                  CASE("FINAL");   CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rFinal)
                  CASE("RATIO");   CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rRatio)
                  CASE("DELTA");   CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rDelta)
                END SELECT
              
              CASE DEFAULT
                IF(BuildTable) THEN
                  Ok=  .false.
                  Msg=       TRIM(W)//"= unknown keyword in PATH CHANGE"
                  RETURN !==============================================
                ELSE
                  DO
                    nStep= nStep+1
                    CALL WrdToReal(W,X)
                    !
                    DO iCp=1,SIZE(vCpn)
                      tTmp(iCp,nStep)= tTmp(iCp,1) &
                      + vStAdd(iCp) /REAL(vStAdd(0)) *X
                    ENDDO
                    !
                    IF(EoL) EXIT
                    IF(nStep==DimMax) EXIT
                    !
                    CALL LinToWrd(L,W,EoL)
                  ENDDO
                  !vNstep(I)= nStep
                ENDIF
                
            END SELECT
          ENDDO
          !---------------------------------------------/ read values --
          !
          !------------------------------------------ build tPathData --
          PathAdd= (rDelta>Zero) !.FALSE.
          !
          IF(iDebug>2) THEN
            PRINT '(A)',"conditions for PATH ADD"
            PRINT '(A)',"stoikiometry vStAdd" ; CALL Pause_
            DO I=1,nCp
              IF(vStAdd(I)/=0) PRINT '(I3,G9.2)',I,vStAdd(I)  
            ENDDO
            !~ PRINT '(A,3G15.6,I3)',"rAddBegin,rAddFinal,rAddRatio,rAddDelta", &
            !~ & rAddBegin,rAddFinal,rAddRatio,rAddDelta
            CALL Pause_
          ENDIF
          !
          IF(BuildTable) THEN
            DO
              nStep= nStep+1
              !
              DO iCp=1,SIZE(vCpn)
                !
                IF(vStAdd(iCp)/=0) THEN
                  !
                  IF(PathAdd) THEN
                    tTmp(iCp,nStep)= tTmp(iCp,nStep-1) &
                    + vStAdd(iCp)/REAL(vStAdd(0)) *rDelta !*(nStep-1)
                  ELSE
                    tTmp(iCp,nStep)= tTmp(iCp,nStep-1) &
                    + vStAdd(iCp)/REAL(vStAdd(0)) *rBegin *rRatio**(nStep-1)
                    !AddedTot= AddedTot +rBegin *rRatio**(nStep-1)
                  ENDIF
                  !
                ELSE
                  tTmp(iCp,nStep)= tTmp(iCp,nStep-1)
                ENDIF
                !
                IF(iDebug>2) WRITE(61,'(G15.6,1X)',ADVANCE="no") tTmp(iCp,nStep)
                !
              ENDDO
              !
              IF(iDebug>2) WRITE(61,*)
              !
              IF(PathAdd) THEN
                AddedTot= AddedTot +rDelta
              ELSE
                AddedTot= AddedTot +rBegin *rRatio**(nStep-1)
              ENDIF
              !
              IF(AddedTot>=rFinal) EXIT
              !
              IF(nStep>=DimMax) EXIT
              !
            ENDDO
          ENDIF
          !
          DimPath= nStep
          !
          ALLOCATE(tPathData(nCp,DimPath))
          tPathData(1:nCp,1:DimPath)= tTmp(1:nCp,1:DimPath)
          !
          ALLOCATE(vTPpath(DimPath))
          vTPpath(:)%TdgC= tTmp(nCp+1,:)
          vTPpath(:)%Pbar= tTmp(nCp+2,:)
          !-----------------------------------------/ build tPathData --
          !
        !_ADD
        !------------------------------------------------/ CASE "ADD" --
        !!! CASE("MIX")
        !!!   !
        !!!   CALL Str_Append(W,3)
        !!!   I=Element_Index(W,vEle)
        !!!   IF(I<1) THEN
        !!!     CALL Stop_(TRIM(W)//"<-NOT in the current system !!!!")
        !!!   ELSE 
        !!!     CALL LinToWrd(L,W,EoL)
        !!!     CALL WrdToReal(W,X1)
        !!!     vSys2(I)=X1
        !!!   ENDIF
        !!! !_MIX
        !
        !------------------------------------------------ CASE "CHG" --
        CASE("CHG","EVP")
          !
          SELECT CASE(TRIM(W))
          !
          CASE("TDGC")  ;  I= nCp+1
          !
          CASE("PBAR")  ;  I= nCp+2
          !
          CASE DEFAULT
            CALL Str_Append(W,3)
            !-- from component name (max 3 chars), find index in vCpn
            I= Component_Index(W,vCpn)
            IF(I<=0) THEN
              Ok=  .false.
              Msg= "PATH CHANGE: Component "//TRIM(W)//" Not in the system .."
              RETURN !==================================================
            ENDIF
          !
          END SELECT
          !
          vLPath(I)=.TRUE.
          !
          IF(I<=nCp) THEN
            J=    vCpn(I)%iSpc !-> find related species
            iMix= vCpn(I)%iMix !-> component is an end-member of mixture phase iMix
          ELSE
            J= 0
            iMix= 0
          ENDIF
          !
          !----------- prim' mobile species is aqueous or pure phase --
          IF(iMix==0) THEN
            !
            rBegin= Zero; rFinal= Zero
            rRatio= One;  rDelta= Zero
            !-------------------------------------------- read values --
            DO
              !
              IF(EoL) EXIT
              CALL LinToWrd(L,W,EoL)
              !
              SELECT CASE(W)
              !
              CASE("INITIAL","FINAL","RATIO","DELTA")
                BuildTable= .true.
                SELECT CASE(W)
                CASE("INITIAL")
                  CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rBegin)
                CASE("FINAL")
                  CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rFinal)
                CASE("RATIO")
                  CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rRatio)
                CASE("DELTA")
                  CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rDelta)
                END SELECT
              !
              CASE DEFAULT
                IF(BuildTable) THEN
                  Ok=  .false.
                  Msg=     TRIM(W)//"= unknown keyword in PATH CHANGE"
                  RETURN !==============================================
                ELSE
                  nStep= 0
                  DO
                    nStep= nStep+1
                    CALL WrdToReal(W,tTmp(I,nStep))
                    IF(EoL) EXIT
                    IF(nStep==DimMax) EXIT
                    CALL LinToWrd(L,W,EoL)
                  ENDDO
                  vNstep(I)= nStep
                ENDIF
                  
              END SELECT
              !
            ENDDO
            !-------------------------------------------/ read values --
            PathAdd= (rDelta>Zero)
            !
            !--------------------------------------------- BuildTable --
            IF(BuildTable) THEN
              IF(I<=nCp) THEN
                !---------------------------- case of chemical change --
                IF(vCpn(I)%Statut=="INERT" .AND.(.NOT. PathAdd)) THEN
                  !->  default CASE (geometric)
                  IF(rRatio < 1.10D0 ) rRatio= 1.10D0
                  IF(rRatio > 3.0D0  ) rRatio= 3.00D0
                  IF(rFinal < rBegin ) rRatio= One/rRatio
                ELSE !MOBILE or BUFFER or DELTA
                  IF(rDelta==Zero)  rDelta= ABS(rFinal-rBegin) /20.0D0
                  !-> default value, in CASE not in input
                  IF(rFinal>rBegin) rDelta= ABS(rDelta)
                  IF(rFinal<rBegin) rDelta=-ABS(rDelta)
                ENDIF
                !---/
              ELSE
                !------------------------------ case of (T,P) changes --
                IF(.NOT. PathAdd) THEN
                  Ok=  .false.
                  Msg=           "PATH CHANGE: for (T,P) only DELTA ..."
                  RETURN !==============================================
                ENDIF
                IF(rDelta==Zero)  rDelta= ABS(rFinal-rBegin) /20.0D0
                !-> default value, in CASE not in input
                IF(rFinal>rBegin) rDelta= ABS(rDelta)
                IF(rFinal<rBegin) rDelta=-ABS(rDelta)
                !---/
              ENDIF
              !
              IF(iDebug>0) WRITE(fTrc,'(I3,A20,4G12.3)') &
              & I," =Begin,Final,Step= ",rBegin,rFinal,rDelta,rRatio
              !
              nStep= 0
              DO
                !
                nStep= nStep+1
                IF(PathAdd) THEN
                  tTmp(I,nStep)= rBegin + (nStep-1) *rDelta
                ELSE
                  tTmp(I,nStep)= rBegin *rRatio**(nStep-1)
                ENDIF
                !
                IF(rFinal>rBegin) THEN
                  IF(tTmp(I,nStep) > rFinal) EXIT
                ELSE
                  IF(tTmp(I,nStep) < rFinal) EXIT
                ENDIF
                !
                IF(nStep==DimMax) EXIT
                !
              ENDDO
              vNstep(I)= nStep-1
              !
            ENDIF
            !--------------------------------------------/ BuildTable --
            !
            ! print *,"unit change"
            IF(I<=nCp) THEN
              IF(vCpn(I)%Statut=="INERT") &
              & tTmp(I,1:vNstep(I))= tTmp(I,1:vNstep(I)) /vCpn(I)%Factor
            ENDIF
            !
            !----------------------------- add amounts in SYSTEM.ROCK --
            IF(I<=nCp) THEN
              IF(vCpn(I)%Statut=="INERT") THEN
                DO J=1,vNstep(I)
                  !! print *,"size(tTmp)=",size(tTmp,1),size(tTmp,2)
                  !! print *,"size(vFas)=",size(vFas)  ;  pause
                  DO K=1,SIZE(vFas)
                    tTmp(I,J)= tTmp(I,J) + &
                    & tAlfFs(I,K) &
                    & *vFas(K)%Mole
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
            !-----------------------------/add amounts in SYSTEM.ROCK --
            !
            IF(iDebug>0) THEN
              WRITE(fTrc,'(I3,A1)',ADVANCE="no") I,t_
              DO J= 1,vNstep(I)
                WRITE(fTrc,'(G12.3,A1)',ADVANCE="no") tTmp(I,J),t_
              ENDDO
              WRITE(fTrc,*)
            ENDIF
            !
          ELSE
          !------ iMix/-0 -> species activity controlled by non-aqueous phase --
            TestMixture= .TRUE.
            !-------------------------------------------- read values --
            DO
              !
              IF(EoL) EXIT
              CALL LinToWrd(L,W,EoL)
              !
              SELECT CASE(W)
                !
              CASE("INITIAL","FINAL")
                !
                CALL LinToWrd(L,W1,EoL) !-> name of initial phase
                K= MixPhase_Index(vMixFas,W1)
                !
                IF(K==0) THEN
                  Ok=  .false.
                  Msg= TRIM(W1)//" = UNKNOWN PHASE in PATH CHANGE!!!"
                  RETURN
                ENDIF
                !
                IF(vMixFas(K)%iModel /= vMixFas(iMix)%iModel)  THEN
                  Ok=  .false.
                  Msg=                                  "Phase "//TRIM(W1)// &
                  &   " should be of same model as "//TRIM(vMixFas(iMix)%Name)
                  RETURN !============================================
                ENDIF
                !
                SELECT CASE(W)
                CASE("INITIAL")
                  vPhasBegin(I)=K
                  
                  IF(iDebug>0) WRITE(fTrc,'(4A)') &
                  & "PhasBegin=",vMixFas(K)%Name, &
                  & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
                  
                CASE("FINAL")
                  vPhasFinal(I)=K
                  
                  IF(iDebug>0) WRITE(fTrc,'(4A)') &
                  & "PhasFinal=",vMixFas(K)%Name, &
                  & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
                  
                ENDSELECT
                !
                !! CASE("RATIO")
                !!   Ok=  .false.
                !!   Msg=                             "ONLY DELTA IN THIS CASE !!!"
                !!   RETURN !----------------------------------------------
                !!   !CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rRatio)
                !! !
                !! CASE("DELTA")
                !!   CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rDelta)
                !!   !currently not used ...!!! (always on 100 steps)
                !
                CASE DEFAULT
                  Ok=  .false.
                  Msg=           TRIM(W)//"= unknown keyword in PATH CHANGE !!!"
                  RETURN !==============================================
                !  
              END SELECT
              !
            ENDDO
            !
            IF(iDebug>1) THEN
              PRINT '(4A)', &
              & "PhasBegin=",vMixFas(vPhasBegin(I))%Name, &
              & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
              PRINT '(4A)', &
              & "PhasFinal=",vMixFas(vPhasFinal(I))%Name, &
              & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
            ENDIF
            !-------------------------------------------/ read values --
            !
            IF(vPhasBegin(I)*vPhasFinal(I)==0) THEN
              Ok= .false.
              Msg=    "Problem in reading PATH conditions on solid solutions ??"
              RETURN !==================================================
            ENDIF
            !
          ENDIF
          !-----/ iMix/-0 -> species activity controlled by non-aqueous phase --
        ! END CASE CHANGE
        !------------------------------------------------/ CASE "CHG" --
        !
        END SELECT
      ENDDO DoBlock
      !------------------------------------/ read the path parameters --
      !! IF(iDebug>0) CALL Pause_
    ENDIF !Cod=="PATH"
    !
  ENDDO DoFile
  CLOSE(f)
  !
  !------------------------------------------------ !! DEVELOPMENT !! --
  IF(TestMixture) THEN
    DimPath= 99
    ALLOCATE(vTPpath(DimPath))
    vTPpath(:)%TdgC= TdgK -T_CK
    vTPpath(:)%Pbar= Pbar
  ENDIF
  !------------------------------------------------/!! DEVELOPMENT !! --
  !
  IF(COUNT(vNstep>0)>0) THEN
    !
    DimPath= MINVAL(vNstep(:),MASK=vNstep(:)>0)
    ALLOCATE(tPathData(nCp,DimPath))
    tPathData(1:nCp,1:DimPath)= tTmp(1:nCp,1:DimPath)
    !! PRINT *,"minval(vNstep,mask=vNstep>0)", minval(vNstep,mask=vNstep>0)
    !! PAUSE
    !
    ALLOCATE(vTPpath(DimPath))
    vTPpath(:)%TdgC= tTmp(nCp+1,1:DimPath)
    vTPpath(:)%Pbar= tTmp(nCp+2,1:DimPath)
    !
  ENDIF
  !
  DEALLOCATE(tTmp)
  DEALLOCATE(vNstep)
  !
  IF(ALLOCATED(vStAdd))  DEALLOCATE(vStAdd)
  !
  IF(iDebug>1) CALL Pause_
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Path_Read"
  
  ! do i= 1, size(tPathData,2)
  !   do j=1,nCp
  !     write(11,'(G11.3,1X)',advance="no") tPathData(j,i)
  !   end do
  !   write(11,*)
  ! end do

END SUBROUTINE Path_ReadParam

SUBROUTINE Path_ReadParam_new( &
& NamFInn,   &  !IN
& PathMode,  &  !IN
& vCpn,      &  !IN
& TdgK,Pbar, &  !IN
& Ok,Msg)       !OUT
!--
!-- read the path parameters --
!-- -> initialize variables of M_Path_Vars
!--   vLPath,vTPpath
!--   vPhasBegin,vPhasFinal
!--   iLogK,vPathLogK
!--   tPathData,DimPath
!--
  USE M_Dtb_Const,   ONLY: T_CK
  USE M_IOTools
  USE M_System_Tools,ONLY: System_TP_Update
  USE M_T_Element,   ONLY: Formula_Read,Element_Index
  USE M_T_Species,   ONLY: Species_Index
  USE M_T_Component, ONLY: T_Component,Component_Index
  USE M_Basis,       ONLY: Basis_Change
  USE M_T_MixPhase,  ONLY: MixPhase_Index
  !
  !--- vars
  USE M_Global_Vars, ONLY: vEle,vSpc,vFas,vMixFas,vMixModel
  !
  USE M_Path_Vars,   ONLY: vLPath,vTPpath
  USE M_Path_Vars,   ONLY: vPhasBegin,vPhasFinal
  USE M_Path_Vars,   ONLY: tPathData,DimPath
  USE M_Path_Vars,   ONLY: iLogK,vPathLogK
  !---/vars

  CHARACTER(LEN=*), INTENT(IN)   :: NamFInn
  CHARACTER(LEN=*), INTENT(IN)   :: PathMode
  TYPE(T_Component),INTENT(IN)   :: vCpn(:)
  REAL(dp),         INTENT(IN)   :: TdgK,Pbar
  LOGICAL,          INTENT(OUT)  :: Ok
  CHARACTER(LEN=*), INTENT(OUT)  :: Msg
  !
  INTEGER,PARAMETER:: DimMax= 255
  !
  CHARACTER(LEN=1024):: L
  CHARACTER(LEN=80) :: W,W1
  !
  LOGICAL :: EoL, fOk, PathAdd, BuildTable, TestMixture
  INTEGER :: ios, ieOx, Z, Zxs, nDiv, iMix
  INTEGER :: I, J, K, f, iCp, nCp, nStep
  REAL(dp):: rBegin,rFinal,rRatio,rDelta,AddedTot,X
  !
  INTEGER,DIMENSION(1:SIZE(vCpn))::vStoik !for PathAdd
  ! REAL(dp):: TdgKAdd,PbarAdd
  !
  ! TYPE(T_Component),ALLOCATABLE:: vCpnAdd(:)
  !
  REAL(dp):: vTmp(1:DimMax)
  REAL(dp):: tTmp(1:SIZE(vCpn)+2,1:DimMax)
  
  INTEGER, ALLOCATABLE:: vNstep(:)
  REAL(dp),ALLOCATABLE:: vStAdd(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Path_Read"
  !
  !IF(iDebug==4) PRINT '(A)',"Path_Read"
  !
  Ok=.FALSE.
  Msg= ""
  TestMixture= .FALSE.
  !
  nCp=SIZE(vCpn)
  !
  ! 1:nCp = components
  ! nCp+1 = TdgC
  ! nCp+2 = Pbar
  tTmp(:,:)= Zero
  tTmp(nCp+1,:)= TdgK -T_CK
  tTmp(nCp+2,:)= Pbar
  !
  ALLOCATE(vNstep(nCp+2))  ;  vNstep(:)= 0
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  !
  IF(iDebug>2) THEN
    DO I=1,SIZE(vCpn)
      !PRINT '(I3,1X,3A)',I,vEle(vCpn(I)%iEle)%NamEl," -> STATUT=",vCpn(I)%Statut
      PRINT '(I3,1X,3A)', &
      & I,vCpn(I)%NamCp," -> STATUT=",vCpn(I)%Statut
    ENDDO
    CALL Pause_
  ENDIF
  !
  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile !============================< end of file ==
    !
    CALL LinToWrd(L,W,EoL)
    !
    IF(W(1:1)=='!') CYCLE DoFile !===============< skip comment lines ==
    CALL AppendToEnd(L,W,EoL)
    !
    BuildTable= .false.
    !
    IF(W=="PATH") THEN
      !
      Ok=.TRUE.
      !
      !---------------- allocate variables according to the Path mode --
      SELECT CASE(PathMode)
      
      CASE("ADD","ADA")
        ALLOCATE(vStAdd(1:nCp))  ; vStAdd= 0.D0
        !
      CASE("CHG")
        ALLOCATE(vLPath(1:nCp+2))  ;  vLPath=.FALSE.
        !
        ALLOCATE(vPhasBegin(1:nCp))  ;  vPhasBegin=0
        ALLOCATE(vPhasFinal(1:nCp))  ;  vPhasFinal=0
        !
      END SELECT
      !----------------------------------------------------/ allocate --
      !
      if(iDebug>2) then
        print *,'Path_ReadParam'
        do i=1,SIZE(vCpn)
          print *,i,vCpn(i)%Mole
        end do
        call pause_
      endif

      !~ DO I=1,SIZE(vCpn)
        !~ IF(vCpn(I)%Statut=="INERT") tTmp(I,:)= vCpn(:)%Mole
        !~ !IF() tTmp(I,:)= vCpn(:)%Mole
      !~ END DO
      !
      !------------------------------------- read the path parameters --
      DoBlock: DO
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!')   CYCLE DoBlock !skip comment lines
        CALL AppendToEnd(L,W,EoL)
        !
        IF(W=="ENDINPUT") EXIT  DoFile
        IF(W=="END" .OR. W=="ENDPATH") EXIT  DoFile
        !-> READ ONLY the first PATH block found !!
        !
        SELECT CASE(TRIM(PathMode))
        !----------------------------------------- CASE "LGK" (-LOGK) --
        CASE("LGK")
          !
          iLogK= Species_Index(TRIM(W),vSpc)
          IF(iLogK==0) THEN
            Ok=  .false.
            Msg=      "In PATH LOGK, "//TRIM(W)//"= unknown species !!!"
            IF(iDebug>0) WRITE(fTrc,'(A)') Msg
            RETURN !====================================================
          ENDIF
          !---------------------------------------------- read values --
          rBegin=Zero; rFinal=Zero; rDelta=Zero
          DO
            IF(EoL) EXIT
            !
            CALL LinToWrd(L,W1,EoL)
            SELECT CASE(W1)
              !
              CASE("INITIAL")
                CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rBegin)
              CASE("FINAL")
                CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rFinal)
              CASE("DELTA")
                CALL LinToWrd(L,W1,EoL); CALL WrdToReal(W1,rDelta)
              CASE DEFAULT
                Ok=  .false.
                Msg= "In PATH LOGK, "//TRIM(W1)//"= invalid keyword !!!"
                IF(iDebug>0) WRITE(fTrc,'(A)') Msg
                RETURN !================================================
              !
            END SELECT
            !
          ENDDO
          !
          IF(ABS(rDelta)<Epsilon(rDelta)) THEN
            Ok=  .false.
            Msg=                       "In PATH LOGK, Delta not defined"
            IF(iDebug>0) WRITE(fTrc,'(A)') Msg
            RETURN !====================================================
          ENDIF
          IF(ABS(rFinal-rBegin) < 0.1) rFinal= rBegin +One
          rDelta= SIGN(rDelta,rFinal-rBegin)
          !-- function SIGN(a,b) returns abs(a)*sign(b) --
          !
          !---------------------------------------------/ read values --
          !
          !------------------------------------------ build vPathLogK --
          I= 1
          vTmp(I)= rBegin
          DO
            IF(I>DimMax) EXIT
            IF(ABS(vTmp(I)-rBegin)>ABS(rFinal-rBegin)) EXIT
            I=I+1
            vTmp(I)= vTmp(I-1) + rDelta
            !!print '(I3,G15.6)',I,vTmp(i)
          ENDDO
          !!pause
          !
          DimPath= I
          ALLOCATE(vPathLogK(1:I))
          vPathLogK(1:I)= vTmp(1:I)
          !-----------------------------------------/ build vPathLogK --
          !
        !----------------------------------------/ CASE "LGK" (-LOGK) --
        !
        !------------------------------------------------- CASE "ADD" --
        CASE("ADD","ADA")
          !
          ! IF(TRIM(W)=="MIX") THEN
          !
          !   ALLOCATE(vCpnAdd(1:SIZE(vCpn)))
          !   vCpnAdd= vCpn
          !
          !   CALL Path_Calc_Fluid( &
          !   & "SYSTEM.MIX",vCpnAdd, &
          !   & TdgKAdd,PbarAdd,Ok)
          !
          !   IF(.NOT. Ok) THEN
          !     !~ CALL Warning_("SYSTEM.MIX NOT FOUND")
          !     DEALLOCATE(vCpnAdd)
          !     Msg= "SYSTEM.MIX NOT FOUND"
          !     RETURN
          !   ENDIF
          !
          !   vStAdd(1:nCp)= vCpnAdd(1:nCp)%Mole
          !
          !   DEALLOCATE(vCpnAdd)
          !
          ! ELSE

            !------------------------- read stoichio of added material --
            CALL Formula_Read(W,vEle,Z,nDiv,fOk,vStoik)
            !
            IF(Z/=0) THEN
              Ok=  .false.
              Msg=     "In PATH ADD, "//TRIM(W)//"= should be neutral !!!"
              IF(iDebug>0) WRITE(fTrc,'(A)') Msg
              RETURN !====================================================
            ENDIF
            !
            !----------------- compute stoikio coeff of redox component --
            ieOx= Element_Index("OX_",vEle)
            Zxs= DOT_PRODUCT(vStoik(:),vEle(:)%Z)
            !-- -> oxidation state
            IF(Zxs/=Z) THEN
              IF(ieOx/=0) THEN
                vStoik(ieOx)= Z - Zxs
                fOk= .TRUE.
              ELSE
                fOk= .FALSE.
              ENDIF
            ENDIF
            !----------------/ compute stoikio coeff of redox component --
            !
            IF(.NOT.fOk) THEN
              Ok=  .false.
              Msg= "In PATH ADD, "//TRIM(W)//"= problem in stoikiometry ?"
              IF(iDebug>0) WRITE(fTrc,'(A)') Msg
              RETURN !====================================================
            ENDIF
            !
            vStAdd(1:nCp)= vStoik(vCpn(1:nCp)%iEle) /REAL(nDiv)
            !------------ vStAdd(1:nCp)/vStAdd(0) is the stoikio'vector --
          ! ENDIF
          !-------------------------/ read stoichio of added material --
          !
          rBegin=  1.D-6
          rFinal=  1.0D0
          rRatio=  1.5D0
          rDelta=  Zero
          !
          nStep= 1
          AddedTot= Zero
          tTmp(1:nCp,1)= vCpn(:)%Mole
          !
          !---------------------------------------------- read values --
          DO
            IF(EoL) EXIT
            !
            CALL LinToWrd(L,W,EoL)
            SELECT CASE(W)

            CASE("INITIAL","FINAL","RATIO","DELTA")
              BuildTable= .true.
              SELECT CASE(W)
                CASE("INITIAL"); CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rBegin)
                CASE("FINAL");   CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rFinal)
                CASE("RATIO");   CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rRatio)
                CASE("DELTA");   CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rDelta)
              END SELECT

            CASE DEFAULT
              IF(BuildTable) THEN
                Ok=  .false.
                Msg=          TRIM(W)//"= unknown keyword in PATH ADD"
                IF(iDebug>0) WRITE(fTrc,'(A)') Msg
                RETURN !==============================================
              ELSE
                DO
                  nStep= nStep+1
                  CALL WrdToReal(W,X)
                  !
                  DO iCp=1,SIZE(vCpn)
                    tTmp(iCp,nStep)= tTmp(iCp,1) + vStAdd(iCp) *X
                  ENDDO
                  !
                  IF(EoL) EXIT
                  IF(nStep==DimMax) EXIT
                  CALL LinToWrd(L,W,EoL)
                ENDDO
                !vNstep(I)= nStep
              ENDIF

            END SELECT
          ENDDO
          !---------------------------------------------/ read values --
          !
          !------------------------------------------ build tPathData --
          PathAdd= (rDelta>Zero) !.FALSE.
          !
          !--- debug
          IF(iDebug>2) THEN
            PRINT '(A)',"conditions for PATH ADD"
            PRINT '(A)',"stoikiometry vStAdd" ; CALL Pause_
            DO I=1,nCp
              IF(vStAdd(I)/=0.D0) PRINT '(I3,G15.6)',I,vStAdd(I)
            ENDDO
            !~ PRINT '(A,3G15.6,I3)',"rAddBegin,rAddFinal,rAddRatio,rAddDelta", &
            !~ & rAddBegin,rAddFinal,rAddRatio,rAddDelta
            CALL Pause_
          ENDIF
          !---/debug
          !
          IF(BuildTable) THEN
            DO
              nStep= nStep+1
              !
              DO iCp=1,SIZE(vCpn)
                !
                IF(vStAdd(iCp)/=0.D0) THEN
                  !
                  IF(PathAdd) THEN
                    tTmp(iCp,nStep)= tTmp(iCp,nStep-1) &
                    + vStAdd(iCp) *rDelta !*(nStep-1)
                  ELSE
                    tTmp(iCp,nStep)= tTmp(iCp,nStep-1) &
                    + vStAdd(iCp) *rBegin *rRatio**(nStep-1)
                    !AddedTot= AddedTot +rBegin *rRatio**(nStep-1)
                  ENDIF
                  !
                ELSE
                  tTmp(iCp,nStep)= tTmp(iCp,nStep-1)
                ENDIF
                !
                IF(iDebug>2) WRITE(61,'(G15.6,1X)',ADVANCE="no") tTmp(iCp,nStep)
                !
              ENDDO
              !
              IF(iDebug>2) WRITE(61,*)
              !
              IF(PathAdd) THEN ;  AddedTot= AddedTot +rDelta
              ELSE             ;  AddedTot= AddedTot +rBegin *rRatio**(nStep-1)
              ENDIF
              !
              IF(AddedTot>=rFinal) EXIT
              !
              IF(nStep>=DimMax) EXIT
              !
            ENDDO
          ENDIF
          !
          DimPath= nStep
          !
          IF(ALLOCATED(tPathData)) DEALLOCATE(tPathData)
          ALLOCATE(tPathData(nCp,DimPath))
          tPathData(1:nCp,1:DimPath)= tTmp(1:nCp,1:DimPath)
          !
          IF(ALLOCATED(vTPpath)) DEALLOCATE(vTPpath)
          ALLOCATE(vTPpath(DimPath))
          vTPpath(:)%TdgC= tTmp(nCp+1,:)
          vTPpath(:)%Pbar= tTmp(nCp+2,:)
          !-----------------------------------------/ build tPathData --
          !
        !_ADD
        !------------------------------------------------/ CASE "ADD" --
        !!! CASE("MIX")
        !!!   !
        !!!   CALL Str_Append(W,3)
        !!!   I=Element_Index(W,vEle)
        !!!   IF(I<1) THEN
        !!!     CALL Stop_(TRIM(W)//"<-NOT in the current system !!!!")
        !!!   ELSE
        !!!     CALL LinToWrd(L,W,EoL)
        !!!     CALL WrdToReal(W,X1)
        !!!     vSys2(I)=X1
        !!!   ENDIF
        !!! !_MIX
        !
        !------------------------------------------------- CASE "CHG" --
        CASE("CHG","EVP")
          !
          SELECT CASE(TRIM(W))
          !
          CASE("TDGC")  ;  I= nCp+1
          !
          CASE("PBAR")  ;  I= nCp+2
          !
          CASE DEFAULT
            CALL Str_Append(W,3)
            !-- from component name (max 3 chars), find index in vCpn
            I= Component_Index(W,vCpn)
            IF(I<=0) THEN
              Ok=  .false.
              Msg= "PATH CHANGE: Component "//TRIM(W)//" Not in the system .."
              IF(iDebug>0) WRITE(fTrc,'(A)') Msg
              RETURN !--------------------------------------------return
            ENDIF
          !
          END SELECT
          !
          vLPath(I)=.TRUE.
          !
          IF(I<=nCp) THEN
            J=    vCpn(I)%iSpc !-> find related species
            iMix= vCpn(I)%iMix !-> component is an end-member of mixture phase iMix
          ELSE
            J= 0
            iMix= 0
          ENDIF
          !
          ! if(iDebug>2) then
          !   if(I<=nCp) print *,vCpn(I)%NamCp
          !   print *,'Path_ReadParam, iMix=',iMix
          ! endif
          !----------- prim' mobile species is aqueous or pure phase --
          IF(iMix==0) THEN
          
            rBegin= Zero; rFinal= Zero
            rRatio= One;  rDelta= Zero
            
            !-------------------------------------------- read values --
            DO
              IF(EoL) EXIT
              CALL LinToWrd(L,W,EoL)
              SELECT CASE(W)

                CASE("INITIAL","FINAL","RATIO","DELTA")
                  BuildTable= .true.
                  SELECT CASE(W)
                    CASE("INITIAL"); CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rBegin)
                    CASE("FINAL")  ; CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rFinal)
                    CASE("RATIO")  ; CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rRatio)
                    CASE("DELTA")  ; CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rDelta)
                  END SELECT

                CASE DEFAULT
                  IF(BuildTable) THEN
                    Ok=  .false.
                    Msg=     TRIM(W)//"= unknown keyword in PATH CHANGE"
                    IF(iDebug>0) WRITE(fTrc,'(A)') Msg
                    RETURN !============================================
                  ELSE
                    nStep= 0
                    DO
                      nStep= nStep+1
                      CALL WrdToReal(W,tTmp(I,nStep))
                      IF(EoL) EXIT
                      IF(nStep==DimMax) EXIT
                      CALL LinToWrd(L,W,EoL)
                    ENDDO
                    vNstep(I)= nStep
                  ENDIF

              END SELECT
            ENDDO
            !-------------------------------------------/ read values --
            !
            PathAdd= (rDelta>Zero)
            !
            !--------------------------------------------- BuildTable --
            IF(BuildTable) THEN
              IF(I<=nCp) THEN
                !---------------------------- case of chemical change --
                IF(vCpn(I)%Statut=="INERT" .AND.(.NOT. PathAdd)) THEN
                  !->  default CASE (geometric)
                  IF(rRatio < 1.10D0 ) rRatio= 1.10D0
                  IF(rRatio > 3.0D0  ) rRatio= 3.00D0
                  IF(rFinal < rBegin ) rRatio= One/rRatio
                ELSE !MOBILE or BUFFER or DELTA
                  IF(rDelta==Zero)  rDelta= ABS(rFinal-rBegin) /20.0D0
                  !-> default value, in CASE not in input
                  IF(rFinal>rBegin) rDelta= ABS(rDelta)
                  IF(rFinal<rBegin) rDelta=-ABS(rDelta)
                ENDIF
                !---/
              ELSE
                !------------------------------ case of (T,P) changes --
                IF(.NOT. PathAdd) THEN
                  Ok=  .false.
                  Msg=           "PATH CHANGE: for (T,P) only DELTA ..."
                  RETURN !==============================================
                ENDIF
                IF(rDelta==Zero)  rDelta= ABS(rFinal-rBegin) /20.0D0
                !-> default value, in CASE not in input
                IF(rFinal>rBegin) rDelta= ABS(rDelta)
                IF(rFinal<rBegin) rDelta=-ABS(rDelta)
                !---/
              ENDIF
              !
              IF(iDebug>0) WRITE(fTrc,'(I3,A20,4G12.3)') &
              & I," =Begin,Final,Step= ",rBegin,rFinal,rDelta,rRatio
              !
              nStep= 0
              DO
                !
                nStep= nStep+1
                IF(PathAdd) THEN
                  tTmp(I,nStep)= rBegin + (nStep-1) *rDelta
                ELSE
                  tTmp(I,nStep)= rBegin *rRatio**(nStep-1)
                ENDIF
                !
                IF(rFinal>rBegin) THEN
                  IF(tTmp(I,nStep) > rFinal) EXIT
                ELSE
                  IF(tTmp(I,nStep) < rFinal) EXIT
                ENDIF
                !
                IF(nStep==DimMax) EXIT
                !
              ENDDO
              vNstep(I)= nStep-1
            ENDIF
            !--------------------------------------------/ BuildTable --
            !
            IF(I<=nCp) THEN
              !IF(vCpn(I)%Statut=="INERT") THEN
                DO J=1,vNstep(I)
                  tTmp(I,J)= tTmp(I,J) /vCpn(I)%Factor
                  !if (iDebug>2) then
                  !  print *,"unit change"
                  !  print *,tTmp(I,J)
                  !  call pause_
                  !endif
                ENDDO
              !ENDIF
            ENDIF
            !
            !----------------------------- add amounts in SYSTEM.ROCK --
            ! IF(I<=nCp .AND. vCpn(I)%Statut=="INERT") THEN
            !   DO J=1,vNstep(I)
            !     print *,"size(tTmp)=",size(tTmp,1),size(tTmp,2)
            !     print *,"size(vFas)=",size(vFas)  ;  pause
            !     DO K=1,SIZE(vFas)
            !       tTmp(I,J)= tTmp(I,J) + &
            !       & tAlfFs(I,K) &
            !       & *vFas(K)%Mole
            !     ENDDO
            !   ENDDO
            ! ENDIF
            !-----------------------------/add amounts in SYSTEM.ROCK --
            !
            IF(iDebug>0) THEN
              WRITE(fTrc,'(I3,A1)',ADVANCE="no") I,t_
              DO J= 1,vNstep(I)
                WRITE(fTrc,'(G12.3,A1)',ADVANCE="no") tTmp(I,J),t_
              ENDDO
              WRITE(fTrc,*)
            ENDIF
            !
          ELSE
          !------ iMix/-0 -> species activity controlled by non-aqueous phase --
            TestMixture= .TRUE.
            !-------------------------------------------- read values --
            DO
              !
              IF(EoL) EXIT
              CALL LinToWrd(L,W,EoL)
              !
              SELECT CASE(W)
                !
                CASE("INITIAL","FINAL")
                  !
                  CALL LinToWrd(L,W1,EoL) !-> name of initial phase
                  K= MixPhase_Index(vMixFas,W1)
                  !
                  IF(K==0) THEN
                    Ok=  .false.
                    Msg= TRIM(W1)//" = UNKNOWN PHASE in PATH CHANGE!!!"
                    RETURN
                  ENDIF
                  !
                  IF(vMixFas(K)%iModel /= vMixFas(iMix)%iModel)  THEN
                    Ok=  .false.
                    Msg=                                  "Phase "//TRIM(W1)// &
                    &   " should be of same model as "//TRIM(vMixFas(iMix)%Name)
                    RETURN !============================================
                  ENDIF
                  !
                SELECT CASE(W)
                  CASE("INITIAL")
                    vPhasBegin(I)=K

                    IF(iDebug>0) WRITE(fTrc,'(4A)') &
                    & "PhasBegin=",vMixFas(K)%Name, &
                    & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name

                  CASE("FINAL")
                    vPhasFinal(I)=K

                    IF(iDebug>0) WRITE(fTrc,'(4A)') &
                    & "PhasFinal=",vMixFas(K)%Name, &
                    & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name

                ENDSELECT
                !
                ! CASE("RATIO")
                !   Ok=  .false.
                !   Msg=                             "ONLY DELTA IN THIS CASE !!!"
                !   RETURN !----------------------------------------------
                !   !CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rRatio)
                ! !
                ! CASE("DELTA")
                !   CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,rDelta)
                !   !currently not used ...!!! (always on 100 steps)
                ! !
                ! CASE DEFAULT
                !   Ok=  .false.
                !   Msg=           TRIM(W)//"= unknown keyword in PATH CHANGE !!!"
                !   RETURN !----------------------------------------------
                !
              END SELECT
            ENDDO
            !
            IF(iDebug>1) THEN
              PRINT '(4A)', &
              & "PhasBegin=",vMixFas(vPhasBegin(I))%Name, &
              & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
              PRINT '(4A)', &
              & "PhasFinal=",vMixFas(vPhasFinal(I))%Name, &
              & "-> Model=",vMixModel(vMixFas(K)%iModel)%Name
            ENDIF
            !-------------------------------------------/ read values --
            !
            IF(vPhasBegin(I)*vPhasFinal(I)==0) THEN
              Ok= .false.
              Msg=    "Problem in reading PATH conditions on solid solutions ??"
              RETURN !==================================================
            ENDIF
            !
          ENDIF
          !-----/ iMix/-0 -> species activity controlled by non-aqueous phase --
        !_CHANGE
        !------------------------------------------------/ CASE "CHG" --
        !
        END SELECT
      ENDDO DoBlock
      !------------------------------------/ read the path parameters --
      !! IF(iDebug>0) CALL Pause_
    ENDIF !Cod=="PATH"
    !
  ENDDO DoFile
  CLOSE(f)
  !
  !------------------------------------------------ !! DEVELOPMENT !! --
  IF(TestMixture) THEN
    DimPath= 99
    ALLOCATE(vTPpath(DimPath))
    vTPpath(:)%TdgC= TdgK -T_CK
    vTPpath(:)%Pbar= Pbar
  ENDIF
  !------------------------------------------------/!! DEVELOPMENT !! --
  !
  IF(COUNT(vNstep>0)>0) THEN
    !
    DimPath= MINVAL(vNstep,MASK=vNstep>0)
    ALLOCATE(tPathData(nCp,DimPath))
    tPathData(1:nCp,1:DimPath)= tTmp(1:nCp,1:DimPath)
    !! PRINT *,"minval(vNstep,mask=vNstep>0)", minval(vNstep,mask=vNstep>0)
    !! PAUSE
    !
    ALLOCATE(vTPpath(DimPath))
    vTPpath(:)%TdgC= tTmp(nCp+1,1:DimPath)
    vTPpath(:)%Pbar= tTmp(nCp+2,1:DimPath)
    !
  ENDIF
  !
  DEALLOCATE(vNstep)
  !
  IF(ALLOCATED(vStAdd))  DEALLOCATE(vStAdd)
  !
  IF(iDebug>1) CALL Pause_
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Path_Read"
  
END SUBROUTINE Path_ReadParam_new

END MODULE M_Path_Read
