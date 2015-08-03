MODULE M_SolPhase_Read
!--
!-- tools for reading solution phases
!-- (- phase of variable composition, with assymetric model)
!-- from text files,
!-- and build an array of T_SolPhase (- container for a mixture phase)
!--
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Warning_
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: SolPhase_Read

CONTAINS

SUBROUTINE SolPhase_Read(vSpc,vSolModel,Ok,MsgError)
!--
!-- scan the ELECTROLYTE block(s)
!--
  USE M_IoTools
  USE M_Files,     ONLY: NamFInn
  USE M_T_Species, ONLY: T_Species,Species_Index
  USE M_T_SolModel,ONLY: T_SolModel
  USE M_T_SolPhase,ONLY: T_SolPhase !,SolPhase_ActPole
  !
  USE M_Global_Vars,ONLY: vSolFas
  !
  TYPE(T_Species), INTENT(IN) :: vSpc(:)
  TYPE(T_SolModel),INTENT(IN) :: vSolModel(:)
  LOGICAL,         INTENT(OUT):: Ok
  CHARACTER(*),    INTENT(OUT):: MsgError
  !
  CHARACTER(LEN=255):: L !,sList0 !=elements in Database
  CHARACTER(LEN=80) :: W !,V1,V2
  LOGICAL :: EoL
  INTEGER :: ios,fInn
  !
  TYPE(T_SolModel)  :: SolModel
  TYPE(T_SolPhase)  :: SolFas
  TYPE(T_SolPhase)  :: vSolTmp(10)
  !
  INTEGER :: I,J,K,nS,N
  REAL(dp):: X1

  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< SolPhase_Read"

  Ok= .TRUE.
  MsgError= "Ok"

  CALL GetUnit(fInn)
  OPEN(fInn,FILE=TRIM(NamFInn))

  N=0
  !
  DoFile: DO
    READ(fInn,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)

    CASE("ENDINPUT"); EXIT DoFile

    CASE("ELECTROLYTE")

      ! ELECTROLYTE WATER1
      !   MODEL MYAQUEOUS-MODEL
      !   COMPOSITION
      !     H2O   55.55
      !     CA+2  0.001
      !     NA+   0.003
      !     NACL(AQ) 0.002
      !     H+   1.e-7
      !     OH-  1.e-7
      !   END
      ! END

      CALL LinToWrd(L,W,EoL) !-> phase name

      IF(Species_Index(W,vSpc)>0) THEN
        Ok= .FALSE.
        MsgError= "Solution name should not match with any species name !!!"
        RETURN !--------------------------------------------------RETURN
      ENDIF

      SolFas%Name=TRIM(W)

      DoSolFas: DO
        !
        READ(fInn,'(A)',IOSTAT=ios) L  ;  IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE DoSolFas
        CALL AppendToEnd(L,W,EoL)
        !
        SELECT CASE(W)

        CASE("ENDINPUT")             ; EXIT DoFile

        CASE("END","ENDELECTROLYTE") ; EXIT DoSolFas

        !----------------------------------------------- read model name 
        CASE("MODEL")
        ! read model name and check it is among the model base vSolModel
        ! -> return its index in array vSolModel
          CALL LinToWrd(L,W,EoL)
          !
          !find the model name in vSolModel%Name
          K=0
          DO I=1,SIZE(vSolModel)
            IF(TRIM(W)==TRIM(vSolModel(I)%Name)) THEN
              K= I
              EXIT
            ENDIF
          ENDDO
          IF(K==0) THEN
            Ok= .FALSE.
            MsgError=             TRIM(W)//"= ELECTROLYTE.MODEL UNKNOWN"
            RETURN !----------------------------------------------RETURN
          ENDIF

          IF(iDebug>0) WRITE(fTrc,'(2A)') "MODEL=", TRIM(W)
          SolFas%iModel= K

          nS= vSolModel(K)%nSpecies
          ALLOCATE(SolFas%vXSpecies(nS))
          SolFas%vXSpecies(:)= Zero !!1.0D-32

        !/CASE("MODEL")
        !-----------------------------------------------/read model name 

        !-----------------------------------------/read phase compositon
        CASE("COMPOSITION")
          !
          IF(iDebug>0) WRITE(fTrc,'(A)') "!!! !!!Read_COMPOSITION"
          !
          IF(SolFas%iModel==0) THEN
            Ok= .FALSE.
            MsgError=             "must define MODEL before COMPOSITION"
            RETURN !----------------------------------------------RETURN
          ENDIF
          !
          SolModel=vSolModel(SolFas%iModel)
          !-------------------input on several lines, one per end-member
          DO

            READ(fInn,'(A)') L
            CALL LinToWrd(L,W,EoL)
            IF(W=="!") CYCLE !-> comment line
            CALL AppendToEnd(L,W,EoL)

            IF(W=="END") EXIT
            IF(W=="ENDCOMPOSITION") EXIT

            I= Species_Index(TRIM(W),vSpc)
            !--- search for the species with global index I
            !--- in the species list, vISpecies(:), of solution model SolModel
            !--- and read mole fraction
            J= 0
            !~ PRINT *,"I=",I
            DO K=1,SolModel%nSpecies
              !~ PRINT *,"  vISpecies(K)=",SolModel%vISpecies(K)
              IF(SolModel%vISpecies(K)==I) THEN
                J= I
                EXIT
              ENDIF
            ENDDO

            IF(J==0) THEN
              Ok= .FALSE.
              MsgError=   TRIM(W)//" NOT FOUND in species list of MODEL"
              RETURN !--------------------------------------------RETURN
            ENDIF

            CALL LinToWrd(L,W,EoL) !-> rest of line -> composition
            CALL WrdToReal(W,X1)

            SolFas%vXSpecies(K)= X1

          ENDDO
          !
          N=N+1
          !
          IF(iDebug>0) &
          & WRITE(fTrc,'(I3,1X,A24,A24)') N, SolFas%Name, SolModel%Name

          vSolTmp(N)%Name=   SolFas%Name
          vSolTmp(N)%iModel= SolFas%iModel
          nS= vSolModel(SolFas%iModel)%nSpecies
          ALLOCATE(vSolTmp(N)%vXSpecies(nS))
          vSolTmp(N)%vXSpecies(:)= SolFas%vXSpecies(:)

          DEALLOCATE(SolFas%vXSpecies)

        !/CASE("COMPOSITION")
        !-----------------------------------------/read phase compositon
        !
        END SELECT !CASE(W)
        !
      ENDDO DoSolFas
    !/CASE(ELECTROLYTE)

    END SELECT !CASE(W)

  ENDDO DoFile
  !
  CLOSE(fInn)

  IF(N>0) THEN
    IF(ALLOCATED(vSolFas)) DEALLOCATE(vSolFas)
    ALLOCATE(vSolFas(N))
    DO I=1,N
      vSolFas(I)%Name= vSolTmp(I)%Name
      vSolFas(I)%iModel= vSolTmp(I)%iModel
      nS= vSolModel(vSolFas(I)%iModel)%nSpecies
      ALLOCATE(vSolFas(I)%vXSpecies(nS))
      vSolFas(I)%vXSpecies(:)= vSolTmp(I)%vXSpecies(:)
    ENDDO
  ENDIF

  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ SolPhase_BuildLnk"

  RETURN
ENDSUBROUTINE SolPhase_Read

ENDMODULE M_SolPhase_Read
