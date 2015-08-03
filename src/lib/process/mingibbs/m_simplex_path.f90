MODULE M_Simplex_Path
  !
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_,Pause_
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Simplex_Path
  !-- series of equilibrium calculations ON PURE PHASES using simplex method
  !-- under variable T,P conditions, or under variable bulk composition
  !
  REAL(dp),ALLOCATABLE,PUBLIC:: tSimplexResult(:,:)
  LOGICAL, ALLOCATABLE,PUBLIC:: vSimplex_Ok(:)
  LOGICAL, ALLOCATABLE,PUBLIC:: vFasIsPresent(:)
  !
CONTAINS

SUBROUTINE Simplex_Path(Cod)
!--
!-- series of equilibrium calculations on pure phases
!-- using simplex method
!-- under variable T,P conditions,
!-- or under variable bulk composition
!--
  USE M_IoTools,     ONLY: OutStrVec
  USE M_Dtb_Const,   ONLY: T_CK, R_JK
  USE M_Global_Tools,ONLY: Global_TP_Update
  USE M_TPcond_Read, ONLY: TPpath_Read
  USE M_Path_Read,   ONLY: Path_ReadMode, Path_ReadParam_new
  USE M_Path_Vars,   ONLY: Path_Vars_Clean
  !
  USE M_Simplex_Build,ONLY: Simplex_CloseFiles
  USE M_Simplex_Calc, ONLY: Simplex_Calc
  USE M_Simplex_Vars, ONLY: Simplex_Vars_Alloc,Simplex_Vars_Clean
  !
  USE M_Simplex_Vars, ONLY: tSimplex
  !
  USE M_Files,       ONLY: NamFInn
  USE M_Global_Vars, ONLY: vFas,vSpcDtb,vSpc,vMixModel
  USE M_Global_Vars, ONLY: vDiscretModel,vDiscretParam,vMixFas
  !
  USE M_Path_Vars,   ONLY: vTPpath,vLPath,tPathData,DimPath
  !
  USE M_GEM_Vars,    ONLY: vCpnGEM,TdgK,Pbar,tStoikioGEM
  !---------------------------------------------------------------------
  CHARACTER(LEN=*),INTENT(IN) :: Cod
  !---------------------------------------------------------------------
  INTEGER :: I,J,nC,nF
  REAL(dp):: TdgK0,Pbar0
  LOGICAL :: Ok
  CHARACTER(LEN=3) :: PathMod3
  CHARACTER(LEN=80):: Msg
  !
  ! REAL(dp),ALLOCATABLE:: tSimplexResult(:,:)
  ! LOGICAL, ALLOCATABLE:: vSimplex_Ok(:)
  ! LOGICAL, ALLOCATABLE:: vFasIsPresent(:)
  ! 
  ! REAL(dp):: Delta
  !---------------------------------------------------------------------
  Ok= .TRUE.
  
  nC= SIZE(vCpnGEM)
  nF= SIZE(vFas)
  
  CALL Simplex_Vars_Alloc(nC,nF)
  
  ! if(iDebug>2) then
  ! print *,'done Simplex_Vars_Alloc' !  ;  call pause_
  ! endif
  
  ! tSimplex=
  ! row 0 =    Gibbs energy of phases 1:nF at T,P
  ! column 0 = bulk compos'n
  ! main =     stoikiometry matrix
  tSimplex(1:nC,1:nF)= -TRANSPOSE(tStoikioGEM(1:nF,1:nC)) !stoikio
  tSimplex(0,   1:nF)= -vFas(1:nF)%Grt  !Gibbs energy
  tSimplex(1:nC,0   )=  vCpnGEM(1:nC)%Mole !bulk compos'n
  
  IF(ALLOCATED(vFasIsPresent)) DEALLOCATE(vFasIsPresent)
  ALLOCATE(vFasIsPresent(1:nF))  ; vFasIsPresent=.FALSE.
  
  CALL Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)
  
  SELECT CASE(TRIM(Cod))
  
  CASE DEFAULT
    Msg= TRIM(Cod)//"= invalid code in SIMPLEX PATH"
    Ok= .FALSE.
    PRINT *,TRIM(Msg)
    RETURN
  
  CASE("PATH")
  !------------------------------------------------------------CASE PATH
  !--- calculation at fixed T,P,
  !--- changing step by step the amount of components' mole numbers
    !
    !--------------------------- read path parameters from PATH block --
    CALL Path_ReadMode(NamFInn,PathMod3,Ok,Msg)
    !
    IF(PathMod3 /= "CHG") THEN
      Msg= "Global equilibrium paths: only in CHANGE mode"
      Ok= .FALSE.
      PRINT *,TRIM(Msg)
      RETURN
    ENDIF
    !
    CALL Path_ReadParam_new( &
    & NamFInn,  &
    & PathMod3, &
    & vCpnGEM, &
    & TdgK,Pbar, &
    & Ok,Msg)
    !
    if(iDebug>2) then
    print *,'done Path_ReadParam'  ;  call pause_
    endif
    !
    ALLOCATE(vSimplex_Ok(1:dimPath))  ; vSimplex_Ok=.FALSE.
    !
    IF(ALLOCATED(tSimplexResult)) DEALLOCATE(tSimplexResult)
    ALLOCATE(tSimplexResult(1:nC+nF+2,1:dimPath)); tSimplexResult=Zero
    !
    !~ CALL Global_TP_Update( &
    !~ & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
    !~ & vSpc,vMixModel,vMixFas,vFas)
    !
    DO i=1,dimPath
      !
      TdgK0= TdgK
      Pbar0= Pbar
      !--------------------if change T,P update T,P-dependent properties
      TdgK= vTPpath(I)%TdgC +T_CK
      Pbar= vTPpath(I)%Pbar
      !
      IF(TdgK0 /= TdgK .OR. Pbar0 /= Pbar) &
      & CALL Global_TP_Update( &
      & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
      & vSpc,vMixModel,vMixFas,vFas)
      !
      DO J=1,SIZE(vSpc)
        WRITE(31,'(A,2F7.1,G15.7)') &
        & TRIM(vSpc(J)%NamSp),TdgK-T_CK,Pbar,vSpc(J)%G0rt*R_JK*TdgK
      END DO
      !----------------------------------------------------------------/
      !
      !--- system composition --
      DO J=1,nC
        IF(vLPath(J)) THEN
          vCpnGEM(J)%Mole= tPathData(J,I)
        ENDIF
      ENDDO
      !---/system composition --
      !
      !--------------------------------------------rebuild tSimplex(:,:)
      ! main = stoikiometry matrix
      tSimplex(1:nC,1:nF)= -TRANSPOSE(tStoikioGEM(1:nF,1:nC))
      ! first row= Gibbs energy of phases 1:nF at 'point' iTP
      tSimplex(0,1:nF)= -vFas(1:nF)%Grt
      ! first column= bulk compos'n
      tSimplex(1:nC,0)=  vCpnGEM(1:nC)%Mole
      !----------------------------------------------------------------/
      
      if(iDebug==1) print *,"Simplex_Run_Step=",i
      if(iDebug>2) then
      print *,'call Simplex_Run'  ;  call pause_
      endif
      
      !--------------------------------------------------------------run
      CALL Simplex_Run( &
      & i,TdgK,Pbar,vCpnGEM,tStoikioGEM, &
      & vSimplex_Ok(i))
      !old! CALL Simplex_Run( &
      !old! & i,TdgK,Pbar, &
      !old! & vFasIsPresent,tSimplexResult, &
      !old! & vSimplex_Ok(i))
      !----------------------------------------------------------------/
      !
    ENDDO
  !-----------------------------------------------------------/CASE PATH
  !
  CASE("TP")
  !--------------------------------------------------------------CASE TP
  !--- calculation at fixed composition,
  !--- changing the T,P condition along the TP.TABLE
    !
    !! CALL Dtb_Tabulate("GIBBS")
    !
    IF(iDebug>1) CALL Trace_1
    !
    CALL TPpath_Read
    !
    dimPath= SIZE(vTPpath)
    !
    IF(dimPath<2) THEN
      PRINT '(A)',"NO TP.TABLE available for Simplex TP.run"
      RETURN
    ENDIF
    !
    ALLOCATE(vSimplex_Ok(1:dimPath))
    vSimplex_Ok=.FALSE.
    !
    IF(ALLOCATED(tSimplexResult)) DEALLOCATE(tSimplexResult)
    ALLOCATE(tSimplexResult(1:nC+nF+2,1:dimPath))
    tSimplexResult=Zero
    !
    DO i=1,dimPath
      !
      !---------------------change T,P ; update T,P-dependent properties
      TdgK= vTPpath(I)%TdgC +T_CK
      Pbar= vTPpath(I)%Pbar
      !
      CALL Global_TP_Update( &
      & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
      & vSpc,vMixModel,vMixFas,vFas)
      !----------------------------------------------------------------/
      !
      !--------------------------------------------rebuild tSimplex(:,:)
      !-------------------------------------- main - stoikiometry matrix
      tSimplex(1:nC,1:nF)= -TRANSPOSE(tStoikioGEM(1:nF,1:nC))
      !------------------------------------- first column- bulk compos'n
      tSimplex(1:nC,0)=  vCpnGEM(1:nC)%Mole
      !----------- first row- Gibbs energy of phases 1:nF at 'point' iTP
      tSimplex(0,1:nF)= -vFas(1:nF)%Grt
      !
      if(iDebug==1) print *,"Simplex_Run_Step=",i
      IF(iDebug==4) THEN
        WRITE(fTrc,'(I3,A1,2(G15.6,A1))',ADVANCE="NO") &
        & I,T_, TdgK-T_CK,T_, Pbar,T_
        CALL OutStrVec(fTrc,vFas%Grt)
      ENDIF
      !-------------------------------------------/rebuild tSimplex(:,:)
      !
      !--------------------------------------------------------------run
      CALL Simplex_Run( &
      & i,TdgK,Pbar,vCpnGEM,tStoikioGEM, &
      & vSimplex_Ok(i))
      !old! CALL Simplex_Run( &
      !old! & i,TdgK,Pbar, &
      !old! & vFasIsPresent,tSimplexResult, &
      !old! & vSimplex_Ok(i))
      !-------------------------------------------------------------/run
      !
    ENDDO
    !
    DEALLOCATE(vTPpath)
  !-------------------------------------------------------------/CASE TP
  ENDSELECT
  !
  !CALL Files_Close
  !
  CALL WriteSysComp(dimPath,vCpnGEM) !,vSimplex_Ok,vFasIsPresent,tSimplexResult)
  !
  IF(ALLOCATED(vSimplex_Ok)) DEALLOCATE(vSimplex_Ok)
  DEALLOCATE(tSimplexResult)
  DEALLOCATE(vFasIsPresent)
  !
  IF(ALLOCATED(tStoikioGEM)) DEALLOCATE(tStoikioGEM)
  !
  IF(ALLOCATED(vTPpath)) DEALLOCATE(vTPpath)
  !
  CALL Simplex_Vars_Clean
  CALL Simplex_CloseFiles
  !
  CALL Path_Vars_Clean
  !
CONTAINS

SUBROUTINE Trace_1

  WRITE(fTrc,'(/,A)') "SPL2"
  WRITE(fTrc,'(3(A,A1))',ADVANCE="NO") "Count",T_,"TdgC",T_,"Pbar",T_
  DO I=1,SIZE(vFas)
    WRITE(fTrc,'(A15,A1)',ADVANCE="NO") vFas(I)%NamFs,T_
  ENDDO
  WRITE(fTrc,*)

  DO I=1,SIZE(vFas) !! check phase / species matching ...!!
    WRITE(fTrc,'(A15,A1)',ADVANCE="NO") vSpc(vFas(I)%iSpc)%NamSp,T_
  ENDDO
  WRITE(fTrc,*)

ENDSUBROUTINE Trace_1

ENDSUBROUTINE Simplex_Path

SUBROUTINE Simplex_Run( &
& iCount,TdgK,Pbar,vCpn,&
& tStoikio,Ok)

  USE M_T_Component, ONLY: T_Component
  USE M_Simplex_Vars,ONLY: tSimplex,IPOSV,iZRov
  USE M_Simplex_Calc,ONLY: Simplex_Calc
  !
  INTEGER,          INTENT(IN) :: iCount
  REAL(dp),         INTENT(IN) :: TdgK,Pbar
  TYPE(T_Component),INTENT(IN) :: vCpn(:)
  REAL(dp),         INTENT(IN) :: tStoikio(:,:)
  LOGICAL,          INTENT(OUT):: Ok
  !
  INTEGER:: iError !,n1,n2

  CALL Simplex_Calc(iError) ! &
  !! & tSimplex,IZROV,IPOSV,  &
  !! & iError,n1,n2)
  !
  Ok= iError==0

  SELECT CASE(iError)

  CASE(0)

    CALL Simplex_WriteResult(TdgK,Pbar,vCpn,tStoikio)
    CALL Simplex_StoreResult(iCount,vCpn(:)%Mole,TdgK,Pbar)

    ! CALL Simplex_WriteResult( &
    ! & TdgK,Pbar, &
    ! & tSimplex,iPosV,iZRov, &
    ! & vFasIsPresent)

    ! ! save results in table tSimplexResult
    ! CALL Simplex_StoreResult( &
    ! & iCount,TdgK,Pbar, &
    ! & IPOSV,tSimplex,   &
    ! & vFasIsPresent,tSimplexResult)

    IF(iDebug>2) CALL Simplex_Print(SIZE(vCpn))

  CASE(1)

    WRITE(fTrc,'(A)') "Unbounded Objective Function"
    IF(iDebug>0) &
    & PRINT *,"SIMPLEX, Error: Unbounded Objective Function"

  CASE(-1)

    WRITE(fTrc,'(A)') "No Solutions Satisfy Constraints Given"
    IF(iDebug>0) &
    & PRINT *,"SIMPLEX, Error: No Solutions Satisfy Constraints Given"

  CASE(-2)

    WRITE(fTrc,'(A)') "BAD INPUT TABLEAU IN SIMPLX"
    IF(iDebug>0) &
    & PRINT *,"SIMPLEX, Error: BAD INPUT TABLEAU IN SIMPLX"

  ENDSELECT

ENDSUBROUTINE Simplex_Run

SUBROUTINE Simplex_Print(nC)
  USE M_Simplex_Vars,ONLY: IPOSV,tSimplex
  USE M_Global_Vars, ONLY: vFas
  !
  INTEGER,INTENT(IN):: nC
  !
  INTEGER:: i
  !
  DO i=1,nC
    PRINT '(G15.6,2X,A)',tSimplex(i,0),TRIM(vFas(IPOSV(i))%NamFs)
  ENDDO
  PRINT *,""
  !
ENDSUBROUTINE Simplex_Print

!-----------------------------------------------------------------------
!-- write tabulated results for all phases found stable at any step ----
!-----------------------------------------------------------------------
SUBROUTINE WriteSysComp(DimPath,vCpn)
  USE M_Files,      ONLY: DirOut
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_IoTools,    ONLY: GetUnit
  USE M_T_MixModel, ONLY: MaxPole
  USE M_T_Component,ONLY: T_Component
  !
  USE M_Global_Vars,ONLY: vFas,vSpc,tFormula
  USE M_Global_Vars,ONLY: vMixModel,vDiscretModel,vDiscretParam
  !
  INTEGER,          INTENT(IN):: DimPath
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  !
  INTEGER :: I,J,K
  INTEGER :: iPath,iFs,iCnt
  INTEGER :: nC,nF
  INTEGER :: iDis,iDisMod
  INTEGER :: F
  LOGICAL,ALLOCATABLE:: vDisIsPresent(:)
  INTEGER,ALLOCATABLE:: tDisParam(:,:,:)
  !TYPE(T_Species):: Spc
  !
  nC=SIZE(vCpn)
  nF=SIZE(vFas)
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_moles.restab")
  !
  IF(SIZE(vDiscretModel)>0) THEN
    ALLOCATE(vDisIsPresent(SIZE(vDiscretModel)))
    vDisIsPresent(:)= .FALSE.
    DO iFs=1,nF
     IF(vFasIsPresent(iFs)) THEN
        !Spc= vSpc(vFas(iFs)%iSpc)
        iDis= vSpc(vFas(iFs)%iSpc)%iDiscret
        IF(iDis/=0) &
        & vDisIsPresent(vDiscretParam(iDis)%iModel)= .TRUE.
      ENDIF
    ENDDO
    !
    IF(COUNT(vDisIsPresent(:))>0) THEN
      ! N=0
      ! DO J=1,SIZE(vDiscretModel)
      !   IF(vDisIsPresent(J)) &
      !   & N= N +vMixModel(vDiscretModel(J)%iMix)%NPole
      ! ENDDO
      !print *,"dimension N=", N   ;  pause
      ALLOCATE(tDisParam(DimPath,3*SIZE(vDiscretModel),3)) !MaxPole))
      tDisParam(:,:,:)= 0
    ENDIF
  ENDIF
  !
  !-----------------------------------------------------------title line
  WRITE(F,'(3(A,A1))',ADVANCE='NO') "count",T_, "TdgC",T_, "Pbar",T_
  DO iFs=1,nC
    WRITE(F,'(A,A1)',ADVANCE='NO') TRIM(vCpn(iFs)%NamCp)//"_cpn", T_
  ENDDO
  ! DO iFs=1,nC
  !   WRITE(F,'(A,A1)',ADVANCE='NO') TRIM(vCpn(iFs)%NamCp)//"_frac", T_
  ! ENDDO
  DO iFs=1,nF
    IF(vFasIsPresent(iFs)) &
    & WRITE(F,'(A,A1)',ADVANCE='NO') TRIM(vFas(iFs)%NamFs), T_
  ENDDO
  WRITE(F,*)
  !----------------------------------------------------------/title line
  !
  iCnt=0
  DO iPath=1,DimPath

    IF(vSimplex_Ok(iPath)) THEN

      iCnt=iCnt+1
      WRITE(F,'(I3,A1)',     ADVANCE='NO') iCnt,T_


      WRITE(F,'(2(F12.3,A1))',ADVANCE='NO') &
      & tSimplexResult(nC+nF+1,iPath)-T_CK, T_, &
      & tSimplexResult(nC+nF+2,iPath),      T_

      !--- mole numbers of components
      DO iFs=1,nC
        WRITE(F,'(G15.6,A1)',ADVANCE='NO') tSimplexResult(iFs,iPath),T_
      ENDDO
      !---/mole numbers of components

      !~ Sum_=SUM(tSimplexResult(1:nC,iPath))
      !~ DO iFs=1,nC
        !~ WRITE(F,'(G15.6,A1)',ADVANCE='NO') tSimplexResult(iFs,iPath)/Sum_,T_
      !~ ENDDO

      DO iFs=1,nF

        IF(vFasIsPresent(iFs)) THEN

          WRITE(F,'(G15.3,A1)',ADVANCE='NO') &
          ! nr'moles phase          *nr'oxygen in formula
          & tSimplexResult(nC+iFs,iPath)*tFormula(iFs,1),T_

          iDis= vSpc(vFas(iFs)%iSpc)%iDiscret
          IF(iDis/=0) THEN
            IF(tSimplexResult(nC+iFs,iPath)>Zero) THEN
              iDisMod= vDiscretParam(iDis)%iModel
              !
              I= vDiscretParam(iDis)%I
              J= vDiscretParam(iDis)%J
              K= vDiscretParam(iDis)%K
              !
              IF(I>=J .AND. I>=K)      THEN  ;  iDisMod= iDisMod*3 -2
              ELSEIF(J>=I .AND. J>=K)  THEN  ;  iDisMod= iDisMod*3 -1
              ELSE                           ;  iDisMod= iDisMod*3
              ENDIF
              tDisParam(iCnt,iDisMod,1)= I
              tDisParam(iCnt,iDisMod,2)= J
              tDisParam(iCnt,iDisMod,3)= K
              !
              !WRITE(15,'(3(I3,1X))',ADVANCE="NO") I,J,K
            ENDIF
          ENDIF

        ENDIF

      ENDDO
      WRITE(F,*)

    ENDIF

  ENDDO
  !
  WRITE(F,'(A1)') "_"
  !
  CLOSE(F)
  !
  IF(ALLOCATED(tDisParam)) THEN
    CALL GetUnit(F)
    OPEN(F,FILE=TRIM(DirOut)//"_mixtures.restab")
    !
    DO I=1,iCnt
      DO J=1,SIZE(tDisParam,2)
        DO K=1,3
          WRITE(F,'(I3,1X)',ADVANCE="NO") tDisParam(I,J,K)
        ENDDO
      ENDDO
      WRITE(F,*)
    ENDDO
    DEALLOCATE(tDisParam)
    CLOSE(F)
  ENDIF
  !
  IF(ALLOCATED(vDisIsPresent)) DEALLOCATE(vDisIsPresent)
  !
  WRITE(*,'(A)') "Results in ",TRIM(DirOut)//"_moles.restab"
  !
ENDSUBROUTINE WriteSysComp

SUBROUTINE Simplex_StoreResult(iCount,vMolCpn,TdgK,Pbar)
!--
!-- save results in table tSimplexResult
!--
  USE M_Simplex_Vars,ONLY: IPOSV,tSimplex
  !
  INTEGER, INTENT(IN):: iCount
  REAL(dp),INTENT(IN):: vMolCpn(:)
  REAL(dp),INTENT(IN):: TdgK,Pbar
  !
  INTEGER::iC,nC,iFs,N
  !
  nC= SIZE(vMolCpn)
  N= SIZE(tSimplexResult,1)
  !
  tSimplexResult(1:nC,iCount)= vMolCpn(1:nC)
  !
  DO iC=1,nC
    iFs= IPOSV(iC)
    vFasIsPresent(iFs)= .TRUE.
    tSimplexResult(nC+iFs,iCount)= tSimplex(iC,0)
  ENDDO
  !
  tSimplexResult(N-1,iCount)= TdgK
  tSimplexResult(N,  iCount)= Pbar
  !
ENDSUBROUTINE Simplex_StoreResult

SUBROUTINE Simplex_WriteResult(TdgK,Pbar,vCpn,tStoikioCpn)
!--
!-- write system, results, ...
!-- to files xxx_check1.log and xxx_check2.log
!--
  USE M_IoTools
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_Files,      ONLY: DirOut
  USE M_T_Component,ONLY: T_Component
  !
  USE M_Global_Vars, ONLY: vFas
  USE M_Simplex_Vars,ONLY: iZRov,tSimplex,iPosV
  USE M_Simplex_Build,ONLY: fSpl1,fSpl2
  !
  REAL(dp),         INTENT(IN):: TdgK,Pbar
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  REAL(dp),         INTENT(IN):: tStoikioCpn(:,:)
  !
  INTEGER,PARAMETER:: maxOut1=100
  INTEGER :: nCpn,nFs,iCpn,iFas,I
  REAL(dp):: X
  ! CHARACTER(LEN=30):: sFormat
  !
  nCpn= SIZE(vCpn)
  nFs=  SIZE(vFas)
  !
  IF(fSpl2==0) THEN
    CALL GetUnit(fSpl2)
    OPEN(fSpl2,FILE=TRIM(DirOut)//"_check2.log")
  ENDIF
  !
  IF(nFs<=maxOut1) THEN
    !
    ! file xxx_check2.log is written
    ! only when number of phases is not too big !!
    !
    IF(fSpl1==0) THEN
      CALL GetUnit(fSpl1)
      OPEN(fSpl1,FILE=TRIM(DirOut)//"_check1.log")
      WRITE(fSpl1,'(4(A,/),/)') &
      & "line1= list of unstable phases", &
      & "line2= values for unstable phases", &
      & "lines 1..nC= stable phase name, amount, ...", &
      & "... stoikio of formation of unstable phase from stable assemblage"
    ENDIF

    ! ! first line=
    ! ! list of unstable phases
    ! WRITE(fSpl1,'(3(A,A1))',ADVANCE='NO') "_",T_,"_",T_,"_",T_
    ! DO iFas=1,nFs
    !   IF(izrov(iFas)<=nFs) &
    !   & WRITE(fSpl1,'(A,A1)',ADVANCE='NO') TRIM(vFas(IZROV(iFas))%NamFs),T_
    ! ENDDO
    ! WRITE(fSpl1,*)
    !
    WRITE(fSpl1,'(2(A,A1))',ADVANCE='NO')   "_",T_,"_",T_
    WRITE(fSpl1,'(G15.6,A1)',ADVANCE='NO') tSimplex(0,0),T_

    ! second line=
    ! values for unstable phases
    DO iFas=1,nFs
      IF(izrov(iFas)<=nFs) &
      & WRITE(fSpl1,'(G15.6,A1)',ADVANCE='NO') tSimplex(0,iFas),T_
    ENDDO
    WRITE(fSpl1,*)
    !
    DO iCpn=1,nCpn
      ! names of stable phases
      WRITE(fSpl1,'(A,A1)',  ADVANCE='NO') TRIM(vFas(IPOSV(iCpn))%NamFs),T_
      ! amount of stable phase
      WRITE(fSpl1,'(G15.6,A1)', ADVANCE='NO') tSimplex(iCpn,0),T_
      ! value for stable phase
      WRITE(fSpl1,'(G15.6,A1)', ADVANCE='NO') tSimplex(0,IPOSV(iCpn)),T_
      ! stoikio of formation of unstable phase from stable assemblage
      DO iFas=1,nFs
        IF (izrov(iFas)<=nFs) &
        & WRITE(fSpl1,'(G15.6,A1)',ADVANCE='NO') tSimplex(iCpn,iFas),T_
      ENDDO !iFas
      WRITE(fSpl1,*)
    ENDDO !iCpn
    WRITE(fSpl1,*)
    !
  ENDIF
  !
  WRITE(fSpl2,'(A)') "!"
  WRITE(fSpl2,'(A,2G15.3)') "TdgC/PBar",TdgK-T_CK,Pbar
  WRITE(fSpl2,'(A)') "!"
  !
  !----------------------------------------------------Write Composition
  WRITE(fSpl2,'(A,A1)',ADVANCE='NO') "vNamCpn",T_
  DO iCpn=1,nCpn
    WRITE(fSpl2,'(A,A1)',  ADVANCE='NO') TRIM(vCpn(iCpn)%NamCp), T_
  ENDDO
  WRITE(fSpl2,*)
  !
  WRITE(fSpl2,'(A,A1)',ADVANCE='NO') "vMolCpn",T_
  DO iCpn=1,nCpn
    WRITE(fSpl2,'(G15.6,A1)',ADVANCE='NO') vCpn(iCpn)%Mole, T_
  ENDDO
  WRITE(fSpl2,*)
  !
  WRITE(fSpl2,'(A,A1)',  ADVANCE='NO') "_",T_
  DO iCpn=1,nCpn
    WRITE(fSpl2,'(A,A1)',  ADVANCE='NO') TRIM(vCpn(iCpn)%NamCp), T_
  ENDDO
  WRITE(fSpl2,*)
  !---------------------------------------------------/Write Composition

  ! !--------------------------------------------------Write Composition
  ! WRITE(sFormat,'(A,I3,A)') '((A15,A1),',nCpn,'(A15,A1))'
  ! WRITE(fSpl2,sFormat) "vNamCpn",T_,(vCpn(iCpn)%NamCp, T_),iCpn=1,nCpn)
 
  ! WRITE(sFormat,'(A,I3,A)') '((A15,A1),',nCpn,'(G15.6,A1))'
  ! WRITE(fSpl2,sFormat) "vMolCpn",T_,((vCpn(iCpn)%Mole,  T_),iCpn=1,nCpn)
 
  ! WRITE(sFormat,'(A,I3,A)') '((A15,A1),',nCpn,'(A15,A1))'
  ! WRITE(fSpl2,sFormat) "_",T_,      ((vCpn(iCpn)%NamCp,  T_),iCpn=1,nCpn)
  ! !-------------------------------------------------/Write Composition

  !---------------------------------------------------Check Mass Balance
  DO iCpn=1,nCpn

    iFas=IPOSV(iCpn)
    WRITE(fSpl2,'(A,A1)',  ADVANCE='NO') TRIM(vFas(iFas)%NamFs),T_
    !
    vFasIsPresent(iFas)=.TRUE.
    !
    DO I=1,nCpn
      IF( tStoikioCpn(iFas,I) /=Zero ) THEN
        WRITE(fSpl2,'(G15.6,A1)',ADVANCE='NO') &
        & tSimplex(iCpn,0)*tStoikioCpn(iFas,I), T_
      ELSE
        WRITE(fSpl2,'(A,A1)',ADVANCE='NO') "0",T_
      ENDIF
    ENDDO
    WRITE(fSpl2,*)

  ENDDO
  !
  WRITE(fSpl2,'(A)') "___"
  !
  WRITE(fSpl2,'(A15,A1)',ADVANCE='NO') "BALANCE",T_
  DO I=1,nCpn
    !X=Zero
    !DO iCpn=1,nCpn
    !  X=X + tSimplex(iCpn,0)*tStoikioCpn(IPOSV(iCpn),I)
    !ENDDO
    X=SUM(tSimplex(1:nCpn,0)*tStoikioCpn(IPOSV(1:nCpn),I))
    WRITE(fSpl2,'(G15.6,A1)',ADVANCE='NO') X, T_
  ENDDO
  WRITE(fSpl2,*)
  !--------------------------------------------------/Check Mass Balance
  !
  RETURN
ENDSUBROUTINE Simplex_WriteResult

ENDMODULE M_Simplex_Path
