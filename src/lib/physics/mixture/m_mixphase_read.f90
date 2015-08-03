MODULE M_MixPhase_Read
!--
!-- tools for reading mixture phases (- phase of variable composition) from text files,
!-- and build an array of T_MixPhase (- container for a mixture phase)
!--
  USE M_Trace,ONLY: iDebug,fTrc,T_,Warning_
  USE M_Kinds
  USE M_T_MixPhase,ONLY: T_MixPhase
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: MixPhase_BuildLnk
  PUBLIC:: MixPhase_LnkToVec
  !
  PUBLIC:: T_LnkFas
  !
  TYPE:: T_LnkFas
    TYPE(T_MixPhase):: Value
    TYPE(T_LnkFas),POINTER ::Next
  ENDTYPE T_LnkFas

CONTAINS

SUBROUTINE LnkFas_Build(First,Input_,L,P)
  LOGICAL               :: First
  TYPE(T_MixPhase)      :: Input_
  TYPE(T_LnkFas),POINTER:: L
  TYPE(T_LnkFas),POINTER:: P
  !
  IF(First) NULLIFY(L)
  !
  IF(First) THEN
    ALLOCATE(L)     ; NULLIFY(L%next)     ; L%Value=Input_     ; P=>L
  ELSE
    ALLOCATE(P%next); NULLIFY(P%next%next); P%next%Value=Input_; P=>P%next
  ENDIF
  !
ENDSUBROUTINE LnkFas_Build

SUBROUTINE MixPhase_BuildLnk(vSpc,vMixModel,N,LnkFas,Ok,MsgError)
!--
!-- scan the PHASE block(s)
!--
  USE M_IoTools
  USE M_Files,     ONLY: NamFInn
  USE M_T_Species, ONLY: T_Species,Species_Index
  USE M_T_MixModel !,ONLY: T_MixModel,MixModel_Index,MixModel_XPoleToXSite
  USE M_T_MixPhase,ONLY: T_MixPhase !,MixPhase_ActPole
  !
  TYPE(T_Species), INTENT(IN) :: vSpc(:)
  TYPE(T_MixModel),INTENT(IN) :: vMixModel(:)
  INTEGER,         INTENT(OUT):: N
  TYPE(T_LnkFas),  POINTER    :: LnkFas
  LOGICAL,         INTENT(OUT):: Ok
  CHARACTER(*),    INTENT(OUT):: MsgError
  !
  TYPE(T_LnkFas),POINTER:: pFas
  !
  TYPE(T_MixPhase)  :: MixFas
  TYPE(T_MixModel)  :: MixModel
  CHARACTER(LEN=255):: L !,sList0 !=elements in Database
  CHARACTER(LEN=80) :: W !,V1,V2
  LOGICAL :: EoL !,Ok
  INTEGER :: iS,iEl,I,J,K !,M,iPhas,iCp,iSp,I,J,K,M,nCp,fInn !,iEl !,J
  REAL(dp):: X1 !,Act
  INTEGER :: ios,fInn
  REAL(dp),ALLOCATABLE:: vXAtom(:)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< MixPhase_BuildLnk"
  !
  Ok= .TRUE.
  MsgError= "Ok"
  !
  CALL GetUnit(fInn)
  OPEN(fInn,FILE=TRIM(NamFInn))
  !Ok=.FALSE.
  !
  N=0
  !
  DoFile: DO
    !
    READ(fInn,'(A)',IOSTAT=ios) L ; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)
    !
    CASE("ENDINPUT"); EXIT DoFile
    !
    CASE("SOLUTION","MIXTURE") !mixture name and composition
      !
      IF(TRIM(W)=="SOLUTION") &
      & CALL Warning_(TRIM(W)//" soon Obsolete, better use MIXTURE !!!")
      !
      !SOLUTION SOL1
      !  MODEL FELDSPAR_SS
      !  COMPOSITION
      !    POLE2 0.5
      !    POLE3 0.5
      !  ENDCOMPOSITION
      !ENDSOLUTION
      !../..
      CALL LinToWrd(L,W,EoL) !-> phase name, SOL1

      IF(Species_Index(W,vSpc)>0) THEN
        Ok= .FALSE.
        MsgError="Mixture name should not match with any species name !!!"
        RETURN !-------------------------------------------------=RETURN
      ENDIF

      MixFas%Name=TRIM(W)
      !
      DoMixFas: DO
        !
        READ(fInn,'(A)',IOSTAT=ios) L  ;  IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE DoMixFas
        CALL AppendToEnd(L,W,EoL)
        !
        SELECT CASE(W)
        !
        CASE("ENDINPUT")                       ; EXIT DoFile
        !
        CASE("END","ENDSOLUTION","ENDMIXTURE") ; EXIT DoMixFas
        !
        !------------------------------------------------read model name
        CASE("MODEL")
        ! read model name and check it is among the model base vMixModel
        ! -> return its index, iS, in array vMixModel_
          CALL LinToWrd(L,W,EoL)
          !
          !find the mixture name in vMixModel%Name
          iS= MixModel_Index(vMixModel,W)

          IF(iS==0) THEN
            Ok= .FALSE.
            MsgError=                 TRIM(W)//"= MIXTURE.MODEL UNKNOWN"
            RETURN !----------------------------------------------RETURN
          ENDIF

          IF(iDebug>0) WRITE(fTrc,'(2A)') "MODEL=", TRIM(W)
          MixFas%iModel=   iS
          !
          !=MixModel
        !ENDCASE("MODEL")
        !-----------------------------------------------/read model name
        !
        !-----------------------------------------/read phase compositon
        CASE("COMPOSITION")
          !
          IF(iDebug>0) WRITE(fTrc,'(A)') "!!! !!!Read_COMPOSITION"
          !
          IF(MixFas%iModel==0) THEN
            Ok= .FALSE.
            MsgError=             "must define the model beforehand ..."
            RETURN !----------------------------------------------RETURN
          ENDIF
          !
          MixModel=vMixModel(MixFas%iModel)
          !
          !! SELECT CASE(TRIM(MixModel%Model))

          !! CASE("IDEAL","POLE","MOLECULAR","FELSPAR")
          !-------------------input on several lines, one per end-member
            DO
              !
              READ(fInn,'(A)') L; CALL LinToWrd(L,W,EoL)
              !-> end-member name, or "!", or "END", or "ENDCOMPOSITION"
              IF(W=="!") CYCLE !-> comment line
              CALL AppendToEnd(L,W,EoL)
              IF(W=="END") EXIT
              IF(W=="ENDCOMPOSITION") EXIT
              !
              I= Species_Index(TRIM(W),vSpc)
              !
              IF(I==0) THEN
                Ok= .FALSE.
                MsgError=   TRIM(W)//" <-end-member NOT FOUND in species list"
                RETURN !------------------------------------------RETURN
              ENDIF
              !
              !--- search for the end-member with index I
              !--- in the end-member list, vIPole(:), of mixing model MixModel
              !--- and read mole fraction
              K=0
              DO
                K=K+1
                IF(K>MixModel%NPole) EXIT
                !
                IF(I==MixModel%vIPole(K)) THEN
                  !
                  CALL LinToWrd(L,W,EoL) !-> rest of line -> composition
                  CALL WrdToReal(W,X1)
                  MixFas%vLPole(K)= (X1>Zero)
                  MixFas%vXPole(K)=  X1 !save to MixFas%composition
                  !
                  IF(iDebug>0) WRITE(fTrc,'(2(A,I3))') &
                  & "K=",K,", (MixFas%iModel)%vIPole(K)=",MixModel%vIPole(K)
                  !
                  EXIT
                ENDIF
              ENDDO
              !---/ search for the end-member named W
              !
              IF(K>MixModel%NPole) THEN
                Ok= .FALSE.
                MsgError=      TRIM(W)//" is NOT AN END-MEMBER of this model"
                RETURN !------------------------------------------RETURN
              ENDIF
              !
            ENDDO
          !ENDCASE("IDEAL","POLE","MOLECULAR","FELSPAR")
          !
          !-------------------------for SITE models: read site fractions
          !!! CASE("SITE")
          !!! ! input should be in number of atoms (not fractions ...) (???)
          !!! ! one line per site, reproducing the structure
          !!! ! of the SITE block of the MIXTURE.MODEL
          !!!   M=0
          !!!   K=0
          !!!   DoSite: DO !M=1,MixModel%NSite
          !!!     READ(fInn,'(A)') L
          !!!     CALL LinToWrd(L,W,EoL)
          !!!     CALL AppendToEnd(L,W,EoL)
          !!!     !
          !!!     IF(W=="!")              CYCLE DoSite !-> comment line
          !!!     IF(W=="END")            EXIT DoSite
          !!!     IF(W=="ENDCOMPOSITION") EXIT DoSite
          !!!     !
          !!!     M= M+1 !site index
          !!!     IF(M > MixModel%NSite) EXIT DoSite
          !!!     !
          !!!     DO
          !!!     !-- scan the rest of the line
          !!!     !-- -> retrieve composition data
          !!!       CALL LinToWrd(L,W,EoL)
          !!!       IF(W=="!") CYCLE DoSite !END of line -> READ a new line
          !!!       CALL WrdToReal(W,X1) !SAVE to MixFas%composition
          !!!       !MixFas%vLPole(K)=(X1>Zero)
          !!!       !
          !!!       K=K+1
          !!!       IF(K>MixModel%NAtom) EXIT DoSite
          !!!       ! save to vMixFas(J)%composition
          !!!       MixFas%vXAtom(K)= X1 / MixModel%vMulti(M)
          !!!       !
          !!!       IF(iDebug>0) &
          !!!       & WRITE(fTrc,'(I3,1X,2A,G15.6)') K,MixModel%vNamAtom(K),"=",X1
          !!!       !
          !!!     ENDDO
          !!!   ENDDO DoSite
          !!! !ENDCASE("SITE")
          !------------------------/for SITE models: read site fractions
          !
          !! END SELECT !CASE (TRIM(MixModel%Model))
          !
          IF(TRIM(MixModel%Model)=="SITE") THEN
          ! to initialize vLPole, calculate the ideal activities of the endmembers
            !
            ALLOCATE(vXAtom(MixModel%NAtom))
            CALL MixModel_XPoleToXSite( &
            & MixModel, &
            & MixFas%vXPole(1:MixModel%NPole), &
            & vXAtom(1:MixModel%NAtom))
            !
            !------------------------------------------------------trace
            IF(iDebug>0) THEN
              !
              WRITE(fTrc,'(/,A,/)') "POLE FRACTIONS"
              DO i=1,MixModel%NPole
                WRITE(fTrc,'(G12.4,1X)',ADVANCE="NO") MixFas%vXPole(I)
              ENDDO
              WRITE(fTrc,*)
              !
              WRITE(fTrc,'(/,A,/)') "SITE FRACTIONS"
              I= MixModel%vIAtomSite(1)
              DO iEl=1,MixModel%NAtom
                J= MixModel%vIAtomSite(iEl)
                IF(J/=I) THEN
                  I= J
                  WRITE(fTrc,*)
                ENDIF
                WRITE(fTrc,'(A,1X,G12.4,1X)',ADVANCE="NO") &
                & MixModel%vNamAtom(iEl), vXAtom(iEl)
              ENDDO
              WRITE(fTrc,*)
              !
            ENDIF
            !-----------------------------------------------------/trace

            IF(iDebug>0) &
            & WRITE(fTrc,'(/,A,/)') "for site models, initialize vLPole"
            !
            DO K=1,MixModel%NPole
              !IF(MixModel%vHasPole(K)) THEN
                X1= MixModel_Site_ActivIdeal(MixModel,K,vXAtom)
                IF(iDebug>0) WRITE(fTrc,'(A,1X,G15.6)') MixModel%vNamPole(K),X1
                MixFas%vLPole(K)= (X1>Zero)
              !ELSE
              !  MixFas%vLPole(K)= .FALSE.
              !ENDIF
            ENDDO
            !
            IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "end vLPole"
            !
          ENDIF
          !
          !-------------------------------------------------------------
          !-- for the model of the phase, MixModel%Model, to be accepted,
          !-- must check whether all 'active' endmembers of MixModel%Model
          !-- are found in vSpc
          !-------------------------------------------------------------
          DO K=1,MixModel%NPole
            IF(MixFas%vLPole(K)) THEN !ONLY for "active" endmembers
              !
              J=Species_Index(MixModel%vNamPole(K),vSpc)
              !
              IF(J==0) THEN
                Ok= .FALSE.
                MsgError= TRIM(MixFas%Name) &
                & //" SPECIES NOT IN BASE:"//TRIM(MixModel%vNamPole(K))
                RETURN !------------------------------------------RETURN
              ENDIF
              !
              MixModel%vIPole(K)=J
              !
            ENDIF
          ENDDO
          N=N+1
          !
          IF(iDebug>0) &
          & WRITE(fTrc,'(I3,1X,A24,A24)') N, MixFas%Name, MixModel%Name
          CALL LnkFas_Build(N==1,MixFas,LnkFas,pFas)
          !
        !ENDCASE("COMPOSITION")
        !----------------------------------------/read phase composition
        !
        END SELECT !CASE(W)!
        !
        IF(ALLOCATED(vXAtom)) DEALLOCATE(vXAtom)
        !
      ENDDO DoMixFas
      !ENDCASE("MIXTURE")
    END SELECT !CASE(W)
  ENDDO DoFile
  !
  CLOSE(fInn)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ MixPhase_BuildLnk"
  !
  RETURN
ENDSUBROUTINE MixPhase_BuildLnk
  !
SUBROUTINE MixPhase_LnkToVec(Lnk,V)
  USE M_T_MixPhase,ONLY: T_MixPhase
  !
  TYPE(T_LnkFas),               POINTER    :: Lnk
  TYPE(T_MixPhase),DIMENSION(:),INTENT(OUT):: V
  !
  TYPE(T_LnkFas),POINTER:: pCur,pPrev
  INTEGER::I
  !
  pCur=>Lnk
  I=0
  DO WHILE (ASSOCIATED(pCur))
    I=I+1
    V(I)=pCur%Value
    pPrev=>pCur; pCur=> pCur%next; DEALLOCATE(pPrev)
  ENDDO
  !
  RETURN
ENDSUBROUTINE MixPhase_LnkToVec

ENDMODULE M_MixPhase_Read
