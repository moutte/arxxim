MODULE M_Solmodel_Read
!--
!-- reading solvent parameters from database
!--
  USE M_Trace,ONLY: iDebug,fTrc,Stop_,T_
  USE M_Kinds

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: Solmodel_Solvent_Read
  PUBLIC:: SolModel_Read

CONTAINS

!--SUBROUTINE Solmodel_Solvent_Read-------------------------------------
!< reads Solvent parameters from database
!! density, dielectric, Debye-Hueckel param's
!-----------------------------------------------------------------------
SUBROUTINE Solmodel_Solvent_Read( &
& nTP,              & ! IN
& vTdgC,            & ! IN
& SolModel,         & ! INOUT
& Rho_Ok, Eps_Ok, DHA_Ok, DHB_Ok, BDot_Ok, & !OUT
& Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl)  !OUT

  USE M_IOTools !, ONLY:dimV,LinToWrd
  USE M_Dtb_Const, ONLY: T_CK
  USE M_Files,     ONLY: NamFLogK,NamFPtz,NamFInn
  USE M_T_SolModel,ONLY: T_SolModel,T_SolModelDat
  USE M_T_SolModel,ONLY: vSolModelAct
  USE M_Solmodel_Vars,ONLY: T_Spline

  INTEGER,           INTENT(IN)   :: nTP       !< dim'n of vTdgC
  REAL(dp),          INTENT(IN)   :: vTdgC(:)  !< temperature array for spline
  TYPE(T_SolModel),  INTENT(INOUT):: SolModel  !< solution model
  LOGICAL,           INTENT(OUT)  :: Rho_Ok,Eps_Ok,DHA_Ok,DHB_Ok,BDot_Ok
  TYPE(T_Spline),    INTENT(INOUT):: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl

  CHARACTER(LEN=512):: L,W
  REAL(dp)::vX(dimV)
  LOGICAL :: Ok,EoL
  !LOGICAL:: Rho_Ok,Eps_Ok,DHA_Ok,DHB_Ok,BDot_Ok
  INTEGER:: ios,mDum,f,I

  TYPE(T_SolModelDat):: vSolvDat(nTP)
  ! REAL(dp),ALLOCATABLE:: B(:), C(:), D(:)

  !--

  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Solmodel_Solvent_Read"

  Ok=.FALSE.
  !
  Rho_Ok= .FALSE.
  Eps_Ok= .FALSE.
  DHA_Ok= .FALSE.
  DHB_Ok= .FALSE.
  BDot_Ok=.FALSE.

  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFLogK),STATUS='OLD')
  !
  DoFile: DO

    READ(f,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile

    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile

    CALL AppendToEnd(L,W,EoL)

    IF(W=="ENDINPUT") EXIT DoFile

    IF(W=="SOLVENT") THEN

      Ok=.TRUE.

      DoSolvent: DO

        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)

        IF(W(1:1)=='!')               CYCLE DoSolvent
        CALL AppendToEnd(L,W,EoL)
        IF(W=="ENDINPUT")             EXIT  DoFile
        IF(W=="END") EXIT DoSolvent
        IF(W=="ENDSOLVENT")           EXIT DoSolvent

        IF(INDEX("RHO_EPS_BDOT_DHA_DHB",TRIM(W))>0) THEN
          CALL ReadRValsV(L,mDum,vX)
          IF(mDum<nTP) CALL Stop_("!!! in SolModel block, DATA nr < dim' TP.table !!!")
        ENDIF

        SELECT CASE(TRIM(W))

        CASE("NAME")
          CALL LinToWrd(L,W,EoL)
          SolModel%Name=TRIM(W)

        CASE("MODEL")
          CALL LinToWrd(L,W,EoL)
          !IF(LEN_TRIM(W)>15) CALL Stop_(TRIM(W)//" = SolModel model name too long !!!")
          !
          SolModel%iActModel= 0
          DO i=1,SIZE(vSolModelAct)
            IF(TRIM(W)==TRIM(vSolModelAct(i))) THEN
              SolModel%iActModel= i
              EXIT
            ENDIF
          ENDDO
          !~ IF(INDEX(Solmodel_ModelList,TRIM(W)//"_")==0) &
          !~ & CALL Stop_(TRIM(W)//" = UNKNOWN ACTIVITY MODEL !!!")

          IF(SolModel%iActModel==0) &
          & CALL Stop_(TRIM(W)//" = UNKNOWN ACTIVITY MODEL !!!")

          !~ SolModel%ActModel=TRIM(W)
          !~ IF(SolModel%ActModel=="PITZER") NamFPtz= TRIM(NamFInn)
          IF(SolModel%iActModel==8) NamFPtz= TRIM(NamFInn)

        CASE("RHO")
          IF(nTP>0) THEN
            vSolvDat(1:nTP)%Rho=  vX(1:nTP)  ;   Rho_Ok= .TRUE.
          ENDIF

        CASE("EPS")
          IF(nTP>0) THEN
            vSolvDat(1:nTP)%Eps=  vX(1:nTP)  ;   Eps_Ok= .TRUE.
          ENDIF

        CASE("DHA")
          IF(nTP>0) THEN
            vSolvDat(1:nTP)%DHA=  vX(1:nTP)  ;   DHA_Ok= .TRUE.
          ENDIF

        CASE("DHB")
          IF(nTP>0) THEN
            vSolvDat(1:nTP)%DHB=  vX(1:nTP)  ;   DHB_Ok= .TRUE.
          ENDIF

        CASE("BDOT")
          IF(nTP>0) THEN
            vSolvDat(1:nTP)%BDot= vX(1:nTP)  ;   BDot_Ok=.TRUE.
          ENDIF

        END SELECT

      ENDDO DoSolvent

    ENDIF !W=="SolModel"

  ENDDO DoFile

  CLOSE(f)

  !~ IF(nTP>0) THEN

    IF(Rho_Ok) &
    & CALL Spline_Init(nTP,vTdgC(:),vSolvDat(:)%Rho,Rho_spl)
    IF(Eps_Ok) &
    & CALL Spline_Init(nTP,vTdgC(:),vSolvDat(:)%Eps,Eps_spl)
    IF(DHA_Ok) &
    & CALL Spline_Init(nTP,vTdgC(:),vSolvDat(:)%DHA,DHA_spl)
    IF(DHB_Ok) &
    & CALL Spline_Init(nTP,vTdgC(:),vSolvDat(:)%DHB,DHB_spl)
    IF(BDot_Ok) &
    & CALL Spline_Init(nTP,vTdgC(:),vSolvDat(:)%BDot,BDot_spl)

  !~ ENDIF

  ! A.M Error Detected Here
  !!IF(iDebug>0) WRITE(fTrc,'(A,E15.8)') "DHA",vSolvDat(1)%DHA
  !!IF(iDebug>0) WRITE(fTrc,'(A,E15.8)') "DHB",vSolvDat(1)%DHB

  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Solmodel_Solvent_Read"

ENDSUBROUTINE Solmodel_Solvent_Read

SUBROUTINE SolModel_Read(vSpc,Ok,MsgError)
!--
!-- initialize vSolModel --
!--

  USE M_IOTools

  USE M_Files,     ONLY: NamFInn
  USE M_T_Species, ONLY: T_Species,Species_Index
  USE M_T_SolModel,ONLY: T_SolModel,T_SolModelDat
  USE M_T_SolModel,ONLY: vSolModelAct
  !
  USE M_Global_Vars,ONLY: vSolModel
  !---------------------------------------------------------------------
  TYPE(T_Species), INTENT(IN) :: vSpc(:)
  LOGICAL,         INTENT(OUT):: Ok
  CHARACTER(LEN=*),INTENT(OUT):: MsgError
  !---------------------------------------------------------------------
  TYPE(T_SolModel):: S
  TYPE(T_SolModel):: vSolTmp(10)
  INTEGER:: vIsolTmp(300)
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL :: sEoL
  INTEGER:: ios
  INTEGER:: F
  !LOGICAL:: Rho_Ok,Eps_Ok,DHA_Ok,DHB_Ok,BDot_Ok
  INTEGER:: I,K,N,nS,M
  !---------------------------------------------------------------------
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "< SolModel_Read"

  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFInn))

  Ok= .TRUE.
  MsgError= "OK"

  N= 0

  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,sEol)
    IF(W(1:1)=='!')   CYCLE DoFile
    CALL AppendToEnd(L,W,sEoL)
    !
    IF(W=="ENDINPUT") EXIT DoFile
    !
    !-----------------------------------------build a new solution model
    IF(W=="ELECTROLYTE.MODEL") THEN

      !~ ELECTROLYTE.MODEL MYAQUEOUS-MODEL
        !~ MODEL IDEAL
        !~ SOLVENT H2O
        !~ POLE
          !~ H2O
          !~ OH-
          !~ H+
          !~ CA+2
          !~ NA+
          !~ NACL(AQ)
        !~ END
      !~ END

      !~ CALL MixModel_Zero(S)
      !
      CALL LinToWrd(L,W,sEol) 
      S%Name=TRIM(W)
      !
      IF(iDebug>1) PRINT '(/,A,/)',"Reading electrolyte model "//TRIM(S%Name)
      !
      DoSol: DO
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        !
        CALL LinToWrd(L,W,sEol) !; IF(iDebug>2) WRITE(fTrc,'(A2,A)') "W=",TRIM(W)
        !
        IF(W(1:1)=='!') CYCLE DoSol !skip comment lines
        CALL AppendToEnd(L,W,sEoL)
        !
        SELECT CASE(W)
        !
        CASE("ENDINPUT"); EXIT DoFile
        !
        CASE("END","ENDELECTROLYTE.MODEL")
        !---------------------------------- end of one ELECTROLYTE block
        !-------------------------- save the T_SolModel value in vSolTmp
          IF(S%nSpecies>0 .AND. S%iSolvent>0) THEN
            N=N+1
            vSolTmp(N)%Name=      S%Name
            vSolTmp(N)%iActModel= S%iActModel
            vSolTmp(N)%iSolvent=  S%iSolvent
            vSolTmp(N)%nSpecies=  S%nSpecies
            ALLOCATE(vSolTmp(N)%vISpecies(S%nSpecies))
            DO I=1,S%nSpecies
              vSolTmp(N)%vISpecies(I)= vIsolTmp(I)
            ENDDO
          ENDIF
          !
          EXIT DoSol
        !----------------------/save the T_MixModel value in linked list
        !
        CASE("MODEL") !MODEL IDEAL / POLE / SITE / SPECIAL
        !---------------------------------------------------- read MODEL
          CALL LinToWrd(L,W,sEol)

          M= 0
          DO i=1,SIZE(vSolModelAct)
            IF(TRIM(W)==TRIM(vSolModelAct(i))) THEN
              !SolModel%iActModel= i
              M= I
              EXIT
            ENDIF
          ENDDO
          !~ IF(INDEX(Solmodel_ModelList,TRIM(W)//"_")==0) &
          !~ & CALL Stop_(TRIM(W)//" = UNKNOWN ACTIVITY MODEL !!!")

          IF(M==0) THEN
            Ok= .FALSE.
            MsgError=         TRIM(W)//" = UNKNOWN ACTIVITY MODEL !!!"
            RETURN !----------------------------------------------return
          ENDIF
          !
          !S%ActModel=TRIM(W)
          S%iActModel=M
          !
          IF(iDebug>0) WRITE(fTrc,'(2A)') "MODEL=",TRIM(W)
          !
          CYCLE DoSol
          !not really necessary,
          !... unless the current value of W has changed to POLE or MARGULES ...
        !----------------------------------------------------/read MODEL 
        !
        CASE("SOLVENT")
        !----------------------------------------------read SOLVENT name
          CALL LinToWrd(L,W,sEol)
          K= Species_Index(W,vSpc)
          IF(K==0) THEN
            Ok= .FALSE.
            MsgError=                            "error in SOLVENT name"
            RETURN !----------------------------------------------------
          ENDIF
          S%iSolvent= K
          !
          IF(iDebug>0) WRITE(fTrc,'(2A)') "SOLVENT=",vSpc(S%iSolvent)%NamSp
          !
        !-------------------------------------------/read SOLVENT name--
        !
        !------------------------------------------------read POLE block
        CASE("POLE")
          !
          nS= 0
          !
          IF(iDebug>0) WRITE(fTrc,'(A)') "<< Read_POLE_block"
          !
          DoPole: DO
            !
            READ(F,'(A)') L; CALL LinToWrd(L,W,sEol)
            IF(W(1:1)=='!') CYCLE DoPole
            CALL AppendToEnd(L,W,sEoL)
            !
            IF(TRIM(W)=="END") EXIT DoPole
            IF(TRIM(W)=="ENDPOLE") EXIT DoPole
            !
            K= Species_Index(W,vSpc)
            !
            IF(K<1) THEN
              Ok= .FALSE.
              MsgError=                TRIM(W)//"species not in SPECIES"
              RETURN !--------------------------------------------RETURN
            ENDIF
            IF(vSpc(K)%Typ/="AQU") THEN
              Ok= .FALSE.
              MsgError=                       TRIM(W)//" is not aqueous"
              RETURN !--------------------------------------------RETURN
            ENDIF
            !
            nS= nS+1
            vIsolTmp(nS)= K
            !~ print *,"vIsolTmp(nS)=",vIsolTmp(nS)
            !
          ENDDO DoPole

          S%nSpecies= nS

        END SELECT
        !
      END DO DoSol
      !
    ENDIF !IF(W=="MIXTURE.MODEL" .OR. W=="SOLUTION.MODEL")
  END DO DoFile
  !
  CLOSE(F)

  IF(N>0) THEN

    IF(ALLOCATED(vSolModel)) DEALLOCATE(vSolModel)
    ALLOCATE(vSolModel(N))

    DO I=1,N
      vSolModel(I)%Name=      vSolTmp(I)%Name
      vSolModel(I)%iActModel= vSolTmp(I)%iActModel
      vSolModel(I)%iSolvent=  vSolTmp(I)%iSolvent
      vSolModel(I)%nSpecies=  vSolTmp(I)%nSpecies
      nS= vSolModel(I)%nSpecies
      ALLOCATE(vSolModel(I)%vISpecies(nS))
      vSolModel(I)%vISpecies(1:nS)= vSolTmp(I)%vISpecies(1:nS)
    ENDDO

    IF(iDebug>2) THEN
      PRINT *,"ELECTROLYTE MODELS"
      DO I=1,N
        S= vSolModel(I)
        PRINT *,S%Name,vSolModelAct(S%iActModel),S%nSpecies
        PRINT *,S%vISpecies(1:nS)
      ENDDO
    ENDIF

  ENDIF

  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ SolModel_Read"

END SUBROUTINE SolModel_Read

SUBROUTINE Spline_Init( &
& N, vX, vY, &
& S)
  USE M_Solmodel_Vars,ONLY: T_Spline
  USE M_CMM_Spline

  INTEGER,       INTENT(IN) :: N
  REAL(dp),      INTENT(IN) :: vX(N)
  REAL(dp),      INTENT(IN) :: vY(N)
  !
  TYPE(T_Spline),INTENT(OUT):: S

  REAL(dp):: B(N),C(N),D(N)

  CALL CMM_Spline_Compute(N,vX(1:N),vY(1:N),B,C,D)

  S%Dimm= N
  S%vX(1:N)= vX(1:N)
  S%vY(1:N)= vY(1:N)
  S%vSplineB(1:N)= B(1:N)
  S%vSplineC(1:N)= C(1:N)
  S%vSplineD(1:N)= D(1:N)
  
  RETURN
END SUBROUTINE Spline_Init

END MODULE M_Solmodel_Read

