MODULE M_Simplex_Theriak
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!--
  USE M_Kinds
  USE M_Trace,       ONLY: iDebug,fTrc,T_,Stop_,Pause_
  USE M_T_Phase,     ONLY: T_Phase
  USE M_T_MixPhase,  ONLY: T_MixPhase
  USE M_T_MixModel,  ONLY: T_MixModel,MaxPole
  USE M_GEM_Vars,    ONLY: T_SavPhase
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Simplex_Theriak
  PUBLIC:: Simplex_Theriak_Path
  !
  TYPE:: T_GemPhase
    CHARACTER(LEN=23):: NamFs
    INTEGER :: nC
    REAL(dp):: vStoik(1:12) ! vStoik(1:MaxEle)
    INTEGER :: iModel
    REAL(dp):: vXPole(1:MaxPole)
    REAL(dp):: Grt
    REAL(dp):: Mole
  END TYPE T_GemPhase
  !
  ! TYPE:: T_Phase
  !   CHARACTER(LEN=23):: NamFs
  !   INTEGER :: iSpc != index of species in vSpc, in case of pure phase 
  !   INTEGER :: iMix != index of mixture in vMixFas, in case of mixture phase
  !   INTEGER :: iSol != index of solution in vSolFas, in case of solution phase
  !   REAL(dp):: Grt,VolM3,WeitKg
  !   REAL(dp):: Mole
  ! END TYPE T_Phase
  !
  !----------------------------------------------------private variables
  INTEGER:: F1,F2,FF
  !
  TYPE(T_SavPhase),ALLOCATABLE:: vSavPhase(:)
  TYPE(T_SavPhase),ALLOCATABLE:: tResultMix(:,:)
  !
  REAL(dp),ALLOCATABLE:: tResult(:,:)
  REAL(dp),ALLOCATABLE:: vFasMolPur(:)
  REAL(dp),ALLOCATABLE:: vFasMolMix(:)
  LOGICAL, ALLOCATABLE:: vFasIsPresent(:)
  LOGICAL, ALLOCATABLE:: vModelConvex(:)
  !--------------------------------------------------------------------/
  !
  !---------------------------------------------------private parameters
  REAL(dp),PARAMETER:: MixMinim_TolX= 1.D-4
  REAL(dp),PARAMETER:: GEM_G_Iota= 1.0D-9 !1.0D-9
  INTEGER, PARAMETER:: GEM_IterMax= 500
  !
  LOGICAL:: WarmRestart= .FALSE.  !!.TRUE.   !!
  LOGICAL:: WriteTrace=  .TRUE.   !!.FALSE.  !!
  !--------------------------------------------------------------------/
  !
CONTAINS

SUBROUTINE Simplex_Theriak
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!--
  USE M_IoTools,     ONLY: GetUnit
  USE M_Files,       ONLY: DirOut
  USE M_Simplex_Vars,ONLY: tSimplex
  USE M_Simplex_Vars,ONLY: Simplex_Vars_Alloc,Simplex_Vars_Clean
  !
  USE M_Global_Vars, ONLY: vFas
  USE M_GEM_Vars,    ONLY: vCpnGEM,tStoikioGEM
  !
  !--- for Global_TP_Update
  USE M_Global_Tools,ONLY: Global_TP_Update
  USE M_Global_Vars, ONLY: vSpcDtb,vSpc,vMixModel
  USE M_Global_Vars, ONLY: vDiscretModel,vDiscretParam,vMixFas
  USE M_GEM_Vars,    ONLY: TdgK,Pbar
  !---/
  !
  TYPE(T_SavPhase):: SavPhaseZero
  INTEGER:: iError
  INTEGER:: I,nC,nFpur,nMix
  !---------------------------------------------------------------------
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Simplex_Theriak"
  !
  F1= 0
  F2= 0
  FF= 0
  !
  SavPhaseZero%iModel=      0
  SavPhaseZero%nFas=        0
  SavPhaseZero%tXPole(:,:)= Zero
  SavPhaseZero%vMole(:)=    Zero
  SavPhaseZero%vGrt0(:)=    Zero
  SavPhaseZero%vVol0(:)=    Zero
  !
  nC= SIZE(vCpnGEM)
  nFpur= SIZE(vFas)
  nMix= SIZE(vMixModel)
  !
  !----------------------------------------------------------allocations
  ALLOCATE(vFasMolPur(1:nFpur))          ;  vFasMolPur(:)= Zero
  ALLOCATE(vFasIsPresent(1:nFpur+nMix))  ;  vFasIsPresent=.FALSE.
  !
  if(iDebug>2) then
    print *,"nMix=",nMix
    call pause_
  end if
  !
  ALLOCATE(vFasMolMix(nMix))
  ALLOCATE(vSavPhase(nMix))
  ALLOCATE(vModelConvex(nMix))
  !
  IF(nMix>0) THEN
    !
    vFasMolMix(:)= Zero
    vSavPhase(:)= SavPhaseZero
    !
    vModelConvex(:)=.FALSE.
    !
  END IF
  !---------------------------------------------------------/allocations
  
  !-----------------------------------------------------------open files
  IF(iDebug>2) THEN
    IF(nMix>0) THEN
      CALL Write_Log_Entete(vFas)
      IF(iDebug>3) THEN
        CALL GetUnit(ff)
        OPEN(ff,FILE='log_theriak_minimize.log')
      END IF
    END IF
  END IF
  !----------------------------------------------------------/open files
  
  CALL Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)
  
  !--------------------------------------------------------------SIMPLEX
  CALL Simplex_GEM(iError)
  !-------------------------------------------------------------/SIMPLEX
  !
  IF(iError/=0) THEN
    CALL Error_Show(iError)
  ELSE
    DO I=1,nFpur
      IF(vFasMolPur(I)>Zero) WRITE(11,'(G15.6,2X,A)') &
      & vFasMolPur(I),vFas(I)%NamFs
    END DO
    IF(nMix>0) THEN
      DO I=1,nMix
        IF(vFasMolMix(I)>Zero) WRITE(11,'(G15.6,2X,A)') &
        & vFasMolMix(I),vMixModel(I)%Name
      END DO
    ENDIF
  END IF

  IF(F1>0) CLOSE(F1)
  IF(F2>0) CLOSE(F2)
  IF(FF>0) CLOSE(ff)

  DEALLOCATE(vFasMolPur)
  DEALLOCATE(vFasIsPresent)
  !
  DEALLOCATE(vFasMolMix)
  DEALLOCATE(vSavPhase)
  DEALLOCATE(vModelConvex)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Simplex_Theriak"

  RETURN
END SUBROUTINE Simplex_Theriak

SUBROUTINE Simplex_Theriak_Path
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!-- changing step by step the amount of components' mole numbers
!--
  USE M_IoTools,     ONLY: GetUnit
  USE M_Dtb_Const,   ONLY: T_CK
  USE M_Path_Read,   ONLY: Path_ReadMode, Path_ReadParam_new
  USE M_Simplex_Vars,ONLY: Simplex_Vars_Alloc,Simplex_Vars_Clean
  !
  USE M_Files,       ONLY: NamFInn
  USE M_Files,       ONLY: DirOut
  USE M_IoTools,     ONLY: GetUnit
  !
  USE M_Global_Vars, ONLY: vFas
  USE M_GEM_Vars,    ONLY: vCpnGEM,tStoikioGEM
  USE M_GEM_Vars,    ONLY: TdgK,Pbar
  !
  USE M_Path_Vars,   ONLY: vTPpath,vLPath,tPathData,DimPath
  USE M_Path_Vars,   ONLY: Path_Vars_Clean
  USE M_Simplex_Vars,ONLY: tSimplex
  USE M_GEM_Write,   ONLY: GEM_Write_Phases,GEM_Write_Mixtures
  !
  !--for Global_TP_Update
  USE M_Global_Tools,ONLY: Global_TP_Update
  USE M_Global_Vars, ONLY: vSpcDtb,vSpc,vMixModel
  USE M_Global_Vars, ONLY: vDiscretModel,vDiscretParam,vMixFas
  !---------------------------------------------------------------------
  !
  TYPE(T_SavPhase):: SavPhaseZero
  INTEGER :: iPath,I,J
  INTEGER :: iError
  INTEGER :: nC,nFpur,nMix
  REAL(dp):: TdgK0,Pbar0
  !
  LOGICAL :: Ok
  CHARACTER(LEN=3) :: PathMod3
  CHARACTER(LEN=80):: Msg
  !
  LOGICAL, ALLOCATABLE:: vSimplex_Ok(:)
  !---------------------------------------------------------------------
  
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Simplex_Theriak_Path"
  !
  nC= SIZE(vCpnGEM)
  nFpur= SIZE(vFas) ! + SIZE(vMixModel)
  nMix= SIZE(vMixModel)
  !
  F1= 0
  F2= 0
  FF= 0
  !
  SavPhaseZero%iModel=      0
  SavPhaseZero%nFas=        0
  SavPhaseZero%tXPole(:,:)= Zero
  SavPhaseZero%vMole(:)=    Zero
  SavPhaseZero%vGrt0(:)=    Zero
  SavPhaseZero%vVol0(:)=    Zero
  !
  !---------------------------------read path parameters from PATH block
  CALL Path_ReadMode(NamFInn,PathMod3,Ok,Msg)  !out
  !
  IF(PathMod3 /= "CHG" .AND. &
  &  PathMod3 /= "GRD" ) THEN
    Msg= "Global equilibrium paths: only in CHANGE or GRID modes"
    Ok= .FALSE.
    PRINT *,TRIM(Msg)
    RETURN
  ENDIF
  !
  IF(PathMod3 == "CHG") CALL Path_ReadParam_new( &
  & NamFInn,  &
  & PathMod3, &
  & vCpnGEM, &
  & TdgK,Pbar, &
  & Ok,Msg)
  !
  IF(.NOT. Ok) THEN
    PRINT *,TRIM(Msg)
    RETURN
  ENDIF
  !--------------------------------/read path parameters from PATH block
  
  !----------------------------------------------------------allocations
  ALLOCATE(vFasMolPur(1:nFpur))          ;  vFasMolPur(:)= Zero
  ALLOCATE(vFasIsPresent(1:nFpur+nMix))  ;  vFasIsPresent=.FALSE.
  !
  ALLOCATE(vFasMolMix(nMix))
  ALLOCATE(vSavPhase(nMix))
  ALLOCATE(vModelConvex(nMix))
  !
  IF(ALLOCATED(tResult)) DEALLOCATE(tResult)
  ALLOCATE(tResult(1:nC+nFpur+nMix+2,1:dimPath))  ; tResult=Zero
  IF(nMix>0) ALLOCATE(tResultMix(1:nMix,1:dimPath))
  !
  ALLOCATE(vSimplex_Ok(1:dimPath))  ;  vSimplex_Ok=.FALSE.
  !
  IF(nMix>0) THEN
    !
    vFasMolMix(:)= Zero
    vSavPhase(:)= SavPhaseZero
    !
    vModelConvex(:)=.FALSE.
    !
  ENDIF
  !---------------------------------------------------------/allocations
  !
  !-----------------------------------------------------------open files
  IF(nMix>0) THEN
    IF(iDebug>2) CALL Write_Log_Entete(vFas)
    !
    IF(iDebug>3) THEN
      CALL GetUnit(ff)
      OPEN(ff,FILE='log_theriak_minimize.log')
    ENDIF
  ENDIF
  !----------------------------------------------------------/open files
  !
  CALL Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
  & vSpc,vMixModel,vMixFas,vFas)
  !
  !do i=1,size(vFas)
  !  print *, vFas(i)%VolM3*1.D6," = ",TRIM(vFas(i)%NamFs)
  !end do
  !pause

  DO iPath=1,dimPath
    !
    TdgK0= TdgK
    Pbar0= Pbar
    !
    !----------------------if change T,P update T,P-dependent properties
    TdgK= vTPpath(iPath)%TdgC +T_CK
    Pbar= vTPpath(iPath)%Pbar
    !
    IF(TdgK0 /= TdgK .OR. Pbar0 /= Pbar) &
    & CALL Global_TP_Update( &
    & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, &
    & vSpc,vMixModel,vMixFas,vFas)
    !------------------------------------------------------------------/
    !
    !-------------------------------------------------system composition
    ! print *,"system composition:"
    DO J=1,nC
      IF(vLPath(J)) THEN
        vCpnGEM(J)%Mole= tPathData(J,iPath)
      ENDIF
    ENDDO
    ! call pause_
    !------------------------------------------------------------------/
    !
    tResult(1,     iPath)= TdgK -T_CK
    tResult(2,     iPath)= Pbar
    tResult(3:nC+2,iPath)= vCpnGEM(1:nC)%Mole
    !
    !--------------------------------------------------------SIMPLEX_GEM
    CALL Simplex_GEM(iError)
    !-------------------------------------------------------/SIMPLEX_GEM
    !
    !------------------------store results in tables tResult, tResultMix
    IF(iError/=0) THEN
      !
      vSimplex_Ok(iPath)= .FALSE.
      CALL Error_Show(iError)
      !
    ELSE
      !
      vSimplex_Ok(iPath)= .TRUE.
      !
      tResult(nC+3:nC+nFpur+2,iPath)= vFasMolPur(1:nFpur)
      IF(nMix>0) THEN
        tResult(nC+nFpur+3:nC+nFpur+nMix+2,iPath)= vFasMolMix(1:nMix)
        tResultMix(1:nMix,iPath)= vSavPhase(1:nMix)
        DO i=1,nMix
          ! print *,vSavPhase(i)%vVol0(1)
          IF(vSavPhase(i)%nFas > 1) CALL Check_vXMean(i)
        ENDDO
        ! call pause_ !!
      ENDIF
      !
      ! CALL Path_StoreResult(iPath,vCpnGEM(:)%Mole,TdgK,Pbar)
      !
    ENDIF
    !-----------------------/store results in tables tResult, tResultMix
    !
    IF(iDebug==1) PRINT *,iPath
    IF(iDebug>1) THEN
      PRINT *, &
      & "============================================================", &
      & iPath
      if(iDebug>2) call pause_
    END IF
    !
  ENDDO
  !
  !--------------------------------------------------write result tables
  CALL GEM_Write_Phases( &
  & DimPath, &
  & vFasIsPresent, &
  & vSimplex_Ok, &
  & tResult, &
  & tResultMix)
  !
  IF(SIZE(vMixModel)>0) CALL GEM_Write_Mixtures( &
  & DimPath, &
  & vFasIsPresent, &
  & vSimplex_Ok,&
  & MixMinim_TolX, &
  & tResult, &
  & tResultMix)
  !-------------------------------------------------/write result tables
  !
  DEALLOCATE(vSimplex_Ok)
  DEALLOCATE(tResult)
  DEALLOCATE(vFasMolPur)
  DEALLOCATE(vFasIsPresent)
  DEALLOCATE(vFasMolMix)
  DEALLOCATE(vSavPhase)
  DEALLOCATE(vModelConvex)
  !
  IF(nMix>0) DEALLOCATE(tResultMix)
  !
  CALL Path_Vars_Clean
  !
  IF(F1>0) CLOSE(F1)
  IF(F2>0) CLOSE(F2)
  IF(ff>0) CLOSE(ff)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Simplex_Theriak_Path"
  !
CONTAINS

SUBROUTINE Check_vXMean(i)
!--
!-- for a mixing model with more than one phase in stable assemblage
!-- compute the Gibbs for the mean composition of the different phases.
!-- if the Gibbs is not close to zero, then there is unmixing
!--
  INTEGER,INTENT(IN):: i

  INTEGER :: J,nP,ff
  REAL(dp):: G
  REAL(dp),ALLOCATABLE:: vX(:)

  J=  vSavPhase(i)%iModel
  nP= vMixModel(J)%NPole
  ALLOCATE(vX(1:nP))
  vX(:)= Zero
  !
  DO ff=1,vSavPhase(i)%nFas
    vX(:)= vX(1:nP) &
    &    + vSavPhase(i)%vMole(ff) *vSavPhase(i)%tXPole(ff,1:nP)
  ENDDO
  vX(:)= vX(:) /SUM(vX(:))
  !
  G= MixPhase_Grt( &
  & TdgK,Pbar,    &
  & nP, &
  & vMixModel(J), &
  & vSavPhase(i)%vGrt0(1:nP), &
  & vX(1:nP))
  !
  DEALLOCATE(vX)

  !~ print *,"G(vXMean)..",vMixModel(J)%Name, G  ;  pause

  RETURN
END SUBROUTINE Check_vXMean

END SUBROUTINE Simplex_Theriak_Path

SUBROUTINE Simplex_GEM(iError)
!--
!-- equilibrium calculations on assemblage of pure phases and mixtures,
!-- with THERIAK method (Capitani-Brown,1987)
!--
  USE M_Simplex_Calc,ONLY: Simplex_Calc
  USE M_Simplex_Vars,ONLY: Simplex_Vars_Alloc,Simplex_Vars_Clean
  !
  USE M_Global_Vars, ONLY: vFas,vMixModel
  USE M_GEM_Vars,    ONLY: vCpnGEM,tStoikioGEM
  USE M_GEM_Vars,    ONLY: TdgK,Pbar
  !
  USE M_Simplex_Vars,ONLY: iPosV,tSimplex
  !
  INTEGER, INTENT(OUT):: iError
  !
  REAL(dp),        ALLOCATABLE:: tStoikio(:,:)
  TYPE(T_Phase),   ALLOCATABLE:: vFasPur(:)
  TYPE(T_GemPhase),ALLOCATABLE:: vFas0(:)
  TYPE(T_MixPhase),ALLOCATABLE:: vMixFas0(:)
  TYPE(T_SavPhase),ALLOCATABLE:: vMixFas_Xpole_Init(:)
  TYPE(T_SavPhase):: SavPhaseZero
  !
  ! TYPE(T_MixModel):: MM
  !
  INTEGER :: nC,nF
  INTEGER :: nFpur,nFmix,nMix
  INTEGER :: nFmixSpl
  INTEGER :: Iter
  INTEGER :: I,K
  REAL(dp):: G_Iota
  ! CHARACTER(LEN=30):: sFMT
  !
  G_Iota= GEM_G_Iota
  !
  SavPhaseZero%iModel=      0
  SavPhaseZero%nFas=        0
  SavPhaseZero%tXPole(:,:)= Zero
  SavPhaseZero%vMole(:)=    Zero
  SavPhaseZero%vGrt0(:)=    Zero
  SavPhaseZero%vVol0(:)=    Zero
  !
  nC=    SIZE(vCpnGEM)
  nFpur= SIZE(vFas) !== at this point, only pure phases are in vFas ==
  nMix=  SIZE(vMixModel)
  !
  nFmix= nMix
  !
  !-----------------------------------------------------------Allocation
  ALLOCATE(tStoikio(1:nFpur+nFmix+nC,1:nC))
  ALLOCATE(vFasPur(1:nFpur))
  ALLOCATE(vFas0(1:nFpur+nFmix+nC))
  !-- vMixFas0: One mixture phase for Each mixture model
  ALLOCATE(vMixFas0(nFmix))
  ALLOCATE(vMixFas_Xpole_Init(nFmix))
  !----------------------------------------------------------/Allocation
  !
  tStoikio(1:nFpur,:)= tStoikioGEM(1:nFpur,:)
  vFasPur(1:nFpur)= vFas(1:nFpur)
  !
  DO I=1,nFpur
    vFas0(I)%NamFs=     vFas(I)%NamFs
    vFas0(I)%nC=        nC
    vFas0(I)%vStoik(1:nC)= tStoikioGEM(I,1:nC)
    vFas0(I)%iModel=    0 != no mixing model, because phase is pure
    vFas0(I)%Grt=       vFas(I)%Grt
    vFas0(I)%Mole=      Zero
  ENDDO
  !
  IF(nMix>0) THEN
    DO K=1,nMix
      vMixFas0(K)%Name= vMixModel(K)%Name
      vMixFas0(K)%iModel= K
      vMixFas0(K)%vXPole(:)= Zero
    ENDDO
    CALL MixPhase_XPole_Init(vMixModel,vMixFas_Xpole_Init)
  ENDIF

  !------------------------------------------------------initial simplex
  nFmixSpl= 0
  nF= nFpur
  !
  CALL Simplex_Vars_Alloc(nC,nF)
  !
  tSimplex(1:nC, 0)=  vCpnGEM(1:nC)%Mole
  tSimplex(0,    1:nF)= -vFas0(1:nF)%Grt
  tSimplex(1:nC, 1:nF)= -TRANSPOSE(tStoikio(1:nF,1:nC))
  !
  CALL Simplex_Calc(iError)
  !
  IF(iError/=0) THEN
    DEALLOCATE(tStoikio)
    DEALLOCATE(vFasPur)
    DEALLOCATE(vFas0)
    ! IF(nMix>0) DEALLOCATE(vMixFas0)
    ! IF(nMix>0) DEALLOCATE(vMixFas_Xpole_Init)
    DEALLOCATE(vMixFas0)
    DEALLOCATE(vMixFas_Xpole_Init)
    RETURN
  ENDIF
  !
  !-----------------compute new free energies of formation of all phases
  !---------------------------------------from current stable assemblage
  CALL Gibbs_Change( &
  & nC,nF,IPOSV, &
  & vFas0)
  !--------------------------------------------------------------------/
  !
  ! IF(iDebug>2) &
  ! & CALL GEM_ShowResults(nC,vMixModel,vFas0,IPOSV,tSimplex(:,0))
  ! IF(iDebug>2) pause
  !
  !-----------------------------------------------------/initial simplex

  IF(nMix==0) THEN
    !
    IF(iDebug>2) CALL ShowResults
    CALL StoreResults
    !
  ELSE
    !
    Iter= 0
    DO
      !
      Iter= Iter+1
      !
      !-------------------------------calc' minimals G's of all mixtures
      vFasPur(1:nFpur)%Grt= vFas0(1:nFpur)%Grt
      !
      CALL Mixture_Minimize( &
      & TdgK,Pbar, &
      & vFasPur,vMixModel,  &
      & vMixFas_Xpole_Init, &
      & vMixFas0)
      !----------------------------------------------------------------/
      IF(ALL(vMixFas0(:)%Grt >=-G_Iota)) THEN
        IF(iDebug>1) CALL ShowResults
        CALL StoreResults
        EXIT
      END IF
      !
      !-------- to the current table of G and stoichiometry for simplex,
      !----------------------- append all mixture phases with negative G
      !---------------------------------- as phases of fixed composition
      !------------- (with the composition computed in Mixture_Minimize)
      nF= nFpur + nFmixSpl
      CALL GEM_AddMixtures( &
      & nC, G_Iota, vMixModel, vMixFas0, &
      & vFas0, nF)
      !
      DO I=nFpur+1,nF
        tStoikio(I,1:nC)= vFas0(I)%vStoik(1:nC)
      ENDDO
      !
      CALL Simplex_Vars_Clean
      CALL Simplex_Vars_Alloc(nC,nF)
      !
      tSimplex(1:nC,0   )=  vCpnGEM(1:nC)%Mole
      tSimplex(0,   1:nF)= -vFas0(1:nF)%Grt
      tSimplex(1:nC,1:nF)= -TRANSPOSE(tStoikio(1:nF,1:nC))
      !
      CALL Simplex_Calc(iError)
      !
      IF(iError/=0) THEN
        DEALLOCATE(tStoikio)
        DEALLOCATE(vFasPur)
        DEALLOCATE(vFas0)
        ! IF(nMix>0) DEALLOCATE(vMixFas0)
        ! IF(nMix>0) DEALLOCATE(vMixFas_Xpole_Init)
        DEALLOCATE(vMixFas0)
        DEALLOCATE(vMixFas_Xpole_Init)
        RETURN
      ENDIF
      !
      IF(iDebug>1) CALL ShowResults
      CALL StoreResults
      !
      !-------------compute new free energies of formation of all phases
      !-----------------------------------from current stable assemblage
      CALL Gibbs_Change( &
      & nC,nF,IPOSV, &
      & vFas0)
      !--/
      !
      IF(iDebug>2) THEN
        WRITE(92,'(/,A,I3)') "vFas0%Grt, iter= ", Iter
        DO I=1,nF
          WRITE(92,'(4X,I3,2X,G15.6,2X,A)') &
          & I,vFas0(I)%Grt,TRIM(vFas0(I)%NamFs)
        ENDDO
      ENDIF
      !
      !--------------------------------to vFas(1:nFpur) (= pure phases),
      !---append mixture phases that belong to current stable assemblage
      !-------------------------------------found from last simplex call
      CALL GEM_AppendSPLMixtures(nC,nFpur,iPosV,vFas0,nFmixSpl)
      !
      !CALL Simplex_Vars_Clean
      !
      IF(Iter > GEM_IterMax) EXIT
      !
    ENDDO
    !
    IF(iDebug>2) PRINT '(A,I3)',"Iter= ",Iter
    !
  END IF
  !
  CALL GEM_SaveResults(nFpur,Iter,vMixModel)
  !
  CALL Simplex_Vars_Clean
  !--------------------------------------------------------DesAllocation
  DEALLOCATE(tStoikio)
  DEALLOCATE(vFasPur)
  DEALLOCATE(vFas0)
  ! IF(nMix>0) DEALLOCATE(vMixFas0)
  ! IF(nMix>0) DEALLOCATE(vMixFas_Xpole_Init)
  DEALLOCATE(vMixFas0)
  DEALLOCATE(vMixFas_Xpole_Init)
  !-------------------------------------------------------/DesAllocation
  !
CONTAINS

  SUBROUTINE ShowResults
    INTEGER:: I,J,K,P,M
    
    WRITE(6,'(A)') "PURE="
    DO K=1,nC
      I= IPOSV(K)
      IF(vFas0(I)%iModel==0) &
      & WRITE(6,'(4X,G15.6,2X,A)') tSimplex(K,0), TRIM(vFas0(I)%NamFs)
    ENDDO
    
    WRITE(6,'(A)') "MIX="
    DO K=1,nC
      I= IPOSV(K)
      IF(vFas0(I)%iModel/=0) THEN
        WRITE(6,'(4X,2I3,2X,G15.6,A)',ADVANCE="NO") &
        & I,K,tSimplex(K,0)," X="
        J= vFas0(I)%iModel
        DO P=1,vMixModel(J)%NPole
          WRITE(6,'(F7.3,1X)',ADVANCE="NO") vFas0(I)%vXPole(P)
        ENDDO
        WRITE(6,'(2X,A)') TRIM(vFas0(I)%NamFs)
      ENDIF
    ENDDO
    WRITE(6,*)
    
    RETURN
  END SUBROUTINE ShowResults

  SUBROUTINE StoreResults
    INTEGER:: I,J,K,P
    INTEGER:: nP

    vFasMolPur(:)= Zero
    !
    vSavPhase(:)= SavPhaseZero
    DO I=1,SIZE(vSavPhase)
      vSavPhase(I)%iModel= I
    ENDDO

    DO K=1,nC
      !
      I= IPOSV(K)
      J= vFas0(I)%iModel
      !
      IF(vFas0(I)%iModel==0) THEN
        !
        vFasMolPur(I)= tSimplex(K,0)
        !
      ELSE
        !
        vSavPhase(J)%iModel= J
        vSavPhase(J)%nFas= vSavPhase(J)%nFas +1
        vSavPhase(J)%vMole(vSavPhase(J)%nFas)= tSimplex(K,0)
        nP= vMixModel(J)%NPole
        DO P=1,nP
          vSavPhase(J)%tXPole(vSavPhase(J)%nFas,P)= vFas0(I)%vXPole(P)
          vSavPhase(J)%vGrt0(P)= vFasPur(vMixModel(J)%vIPole(P))%Grt
          vSavPhase(J)%vVol0(P)= vFasPur(vMixModel(J)%vIPole(P))%VolM3
          !print *,"StoreResults",vSavPhase(J)%vVol0(P)
        ENDDO
        !
      ENDIF
      !
    ENDDO

    !--------------------------------------------------------------trace
    IF(iDebug>2) THEN
      DO K=1,nC
        I= IPOSV(K)
        IF(vFas0(I)%iModel/=0) THEN
          WRITE(93,'(A,A1,G15.6,A1)',ADVANCE="NO") &
          & TRIM(vFas0(I)%NamFs),T_, &
          & tSimplex(K,0),       T_
          J= vFas0(I)%iModel
          DO P=1,vMixModel(J)%NPole
            WRITE(93,'(G15.6,A1)',ADVANCE="NO") vFas0(I)%vXPole(P),T_
          ENDDO
        ENDIF
      ENDDO
      WRITE(93,*)
    ENDIF
    !-------------------------------------------------------------/trace

    RETURN
  END SUBROUTINE StoreResults

END SUBROUTINE Simplex_GEM

REAL(dp) FUNCTION MixPhase_Grt(TdgK,Pbar,nP,MM,vMu0rt,vX)
  USE M_T_MixModel,ONLY: MixModel_Activities
  !---------------------------------------------------------------------
  REAL(dp),        INTENT(IN):: TdgK,Pbar
  INTEGER,         INTENT(IN):: nP
  TYPE(T_MixModel),INTENT(IN):: MM        ! mixing model
  REAL(dp),        INTENT(IN):: vMu0rt(:) !
  REAL(dp),        INTENT(IN):: vX(:)     ! phase composition
  !---------------------------------------------------------------------
  REAL(dp):: vLGam(nP),vLIdeal(nP),vLnAct(nP)
  LOGICAL :: vLPole(nP)
  REAL(dp):: G
  INTEGER :: i
  LOGICAL :: Ok
  CHARACTER(LEN=30):: Msg
  !---------------------------------------------------------------------
  vLPole(:)= (vX(:)>Zero) ! .AND. MM%vHasPole(:)

  CALL MixModel_Activities( & !
  & TdgK,Pbar, & ! in
  & MM,        & ! in
  & vX,        & ! in
  & vLPole,    & ! in
  & Ok, Msg,   & ! out
  & vLGam,     & !
  & vLIdeal,   & !
  & vLnAct)      !
  !
  G= Zero
  DO i=1,nP
    IF(vLPole(i)) &
    ! vMu0rt(i)= vFasPur(MM%vIPole(i))%Grt
    & G= G &
    &  + vX(i) *(vMu0rt(i) + vLnAct(i))
  END DO
  !
  MixPhase_Grt= G

  RETURN
END FUNCTION MixPhase_Grt

SUBROUTINE GEM_SaveResults(nFpur,Iter,vMixModel)
  INTEGER,         INTENT(IN):: nFpur
  INTEGER,         INTENT(IN):: Iter
  TYPE(T_MixModel),INTENT(IN):: vMixModel(:)
  !
  INTEGER :: I,J
  
  !--- pure phases --
  DO I=1,nFpur
    IF(vFasMolPur(I)>Zero) vFasIsPresent(I)= .TRUE.
  ENDDO
  
  !--- mixtures --
  vFasMolMix(:)= Zero
  DO I=1,SIZE(vSavPhase)
    IF(vSavPhase(I)%nFas >0) THEN
      vFasIsPresent(nFpur +I)= .TRUE.
      DO J=1,vSavPhase(I)%nFas
        vFasMolMix(I)= vFasMolMix(I) +vSavPhase(I)%vMole(J)
      ENDDO
    ENDIF
  ENDDO
  
  RETURN
END SUBROUTINE GEM_SaveResults

SUBROUTINE GEM_AppendSPLMixtures(nC,nFpur,iPosV,vFas,nF)
!--
!-- to vFas(1:nFpur) (- pure phases),
!-- append mixture phases that belong to current stable assemblage
!-- found from last simplex call
!--
  INTEGER,         INTENT(IN)   :: nC, nFpur
  INTEGER,         INTENT(IN)   :: IPOSV(:)
  TYPE(T_GemPhase),INTENT(INOUT):: vFas(:)
  INTEGER,         INTENT(INOUT):: nF
  !
  TYPE(T_GemPhase),ALLOCATABLE:: vFas1(:)
  INTEGER:: I,K
  !
  ALLOCATE(vFas1(SIZE(vFas)))
  vFas1= vFas
  nF= 0
  DO K=1,nC
    I= IPOSV(K)
    IF(vFas1(I)%iModel/=0) THEN
      nF= nF+1
      vFas(nFpur+nF)= vFas1(I)
    ENDIF
  ENDDO
  DEALLOCATE(vFas1)
  !
END SUBROUTINE GEM_AppendSPLMixtures

SUBROUTINE GEM_AddMixtures( &
& nC, G_Iota, vMixModel, vMixFas, &
& vFas, nFas)
!--
!-- to the current table of G and stoichiometry for simplex calculations,
!-- append all mixture phases with negative G
!-- as phases of fixed composition
!-- (with the composition that makes G minimal, computed in Mixture_Minimize)
!--
  INTEGER,         INTENT(IN)   :: nC
  REAL(dp),        INTENT(IN)   :: G_Iota
  TYPE(T_MixModel),INTENT(IN)   :: vMixModel(:)
  TYPE(T_MixPhase),INTENT(IN)   :: vMixFas(:)
  TYPE(T_GemPhase),INTENT(INOUT):: vFas(:)
  INTEGER,         INTENT(INOUT):: nFas
  !
  TYPE(T_MixModel):: MM
  INTEGER :: iMF,iMdl,iMix,iCp,N

  iMF= nFas
  !
  DO iMix=1,SIZE(vMixFas)
  
    IF(vMixFas(iMix)%Grt <-G_Iota) THEN
      !
      iMdl= vMixFas(iMix)%iModel
      iMF=  iMF +1
      MM=   vMixModel(iMdl)
      N=    MM%NPole
      !
      vFas(iMF)%NamFs=       vMixFas(iMix)%Name
      vFas(iMF)%vXPole(1:N)= vMixFas(iMix)%vXPole(1:N)
      vFas(iMF)%iModel=      iMdl
      vFas(iMF)%Grt=         vMixFas(iMix)%Grt
      !
      DO iCp=1,nC ! stoichiometry of mixture vs components
        vFas(iMF)%vStoik(iCp)= &
        & SUM( vFas(MM%vIPole(1:N))%vStoik(iCp) &
        &    * vMixFas(iMix)%vXPole(1:N) )
      ENDDO
      !
      IF(iDebug>2) PRINT *,"ADDED= ",TRIM(vFas(iMF)%NamFs)
      !
    ENDIF
  
  ENDDO
  !
  nFas= iMF
  !
END SUBROUTINE GEM_AddMixtures

SUBROUTINE Error_Show(iErr)
  INTEGER,INTENT(IN):: iErr
  SELECT CASE(iErr)
    CASE(1)
      WRITE(*,'(A)') "Unbounded Objective Function"
    CASE(-1)
      WRITE(*,'(A)') "No Solutions Satisfy Constraints Given"
    CASE(2)
      WRITE(*,'(A)') "Max Iteration Reached !!!"
  ENDSELECT
END SUBROUTINE Error_Show

SUBROUTINE Gibbs_Change( &
& nC,nF,IPOSV, &
& vFas0)
!--
!-- use stable phase assemblage as the new component set
!-- and recalculate free energies of formation of all species,
!-- as free enegies of formation from these components
!--
  USE M_Numeric_Mat,ONLY: LU_Decomp, LU_BakSub
  !---------------------------------------------------------------------
  INTEGER, INTENT(IN):: nC,nF
  INTEGER, INTENT(IN):: IPOSV(:)
  TYPE(T_GemPhase),INTENT(INOUT):: vFas0(:)
  !---------------------------------------------------------------------
  REAL(dp),ALLOCATABLE:: tTransform(:,:)
  INTEGER, ALLOCATABLE:: vIndx(:)
  REAL(dp),ALLOCATABLE:: vY(:)
  REAL(dp),ALLOCATABLE:: vGrt(:),vGrt0(:)
  !
  INTEGER :: I,K
  LOGICAL :: bSingul
  REAL(dp):: D
  !---------------------------------------------------------------------
  ALLOCATE(tTransform(1:nC,1:nC))
  ALLOCATE(vIndx(1:nC),vY(1:nC))
  ALLOCATE(vGrt(1:nF))
  ALLOCATE(vGrt0(1:nC))
  !
  !---------------------------------build stoikio table of stable phases
  DO K=1,nC
    I= IPOSV(K) ! index of stable phase in phase list
    tTransform(1:nC,K)= vFas0(I)%vStoik(1:nC)
    vGrt0(K)= vFas0(I)%Grt !-> free energy of stable phase
  ENDDO
  !--------------------------------------------------------------------/
  !
  CALL LU_Decomp(tTransform,vIndx,D,bSingul)
  IF(bSingul) CALL Stop_("SINGUL IN Gibbs_Change")
  !
  DO I=1,nF
    vY(1:nC)= vFas0(I)%vStoik(1:nC)
    !-----------------------Solve(A,Row,B,X) Yout= inv(tTransform) * Yin
    CALL LU_BakSub(tTransform,vIndx,vY)
    vGrt(I)= vFas0(I)%Grt - SUM(vY(:)*vGrt0(:))
  ENDDO
  vFas0(1:nF)%Grt= vGrt(1:nF)
  !
  DEALLOCATE(tTransform)
  DEALLOCATE(vGrt,vGrt0)
  DEALLOCATE(vIndx,vY)
  !
END SUBROUTINE Gibbs_Change

SUBROUTINE MixPhase_XPole_Init(vMixModel,vSavFas)
  TYPE(T_MixModel),INTENT(IN):: vMixModel(:)
  TYPE(T_SavPhase),INTENT(INOUT):: vSavFas(:)
  !
  INTEGER :: I,P,Q,nP
  
  DO I=1,SIZE(vSavFas)
    !nP= vMixModel(vSavFas(I)%iModel)%nPole
    nP= vMixModel(I)%nPole
    DO P=1,nP
      !--start from compos'n close to end-member P
      vSavFas(i)%tXPole(P,P)= One - 1.0D-3
      DO Q=1,nP
        IF(Q/=P) vSavFas(i)%tXPole(P,Q)= 1.0D-3/REAL(nP-1)
      ENDDO
    ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE MixPhase_XPole_Init

SUBROUTINE Mixture_Minimize( &
& TdgK,Pbar, &
& vFasPur,vMixModel,  &
& vMixFas_Xpole_Init, &
& vMixFas0)
!--
!-- for each mixture phase,
!-- compute the composition X that minimizes G
!-- and compute G at X, to see whether it is -0
!-- (the value istself is useful only for output to trace file)
!--
!-- if there is no phase with negative G,
!-- then the current phase assemblage is stable
!--
  USE M_Numeric_Const, ONLY: Ln10
  USE M_Safe_Functions,ONLY: FSafe_Exp
  !~ USE M_GEM_Vars,   ONLY: TdgK,Pbar
  !
  USE M_Optimsolver_Theriak
  USE M_MixModel_Optim
  !
  REAL(dp),        INTENT(IN)   :: TdgK,Pbar
  TYPE(T_Phase),   INTENT(IN)   :: vFasPur(:)
  TYPE(T_MixModel),INTENT(IN)   :: vMixModel(:)
  TYPE(T_SavPhase),INTENT(INOUT):: vMixFas_Xpole_Init(:)
  TYPE(T_MixPhase),INTENT(INOUT):: vMixFas0(:)
  !
  TYPE(T_MixModel):: MM
  REAL(dp),ALLOCATABLE:: vX(:),vXmin(:),vMu(:),vMuMin(:)
  REAL(dp):: G,Gmin
  INTEGER :: N,I,J,K,P !,iFas
  INTEGER :: Multi
  REAL(dp):: S
  !
  REAL(dp):: TolX,DeltaInit
  INTEGER :: its,nCallG
  LOGICAL :: OkConverge
  !
  TolX= MixMinim_TolX
  DeltaInit= 0.05D0
  !
  DO I=1,SIZE(vMixFas0)
    !
    MM= vMixModel(vMixFas0(I)%iModel)
    N= MM%NPole
    Multi= MM%vMulti(1)
    !
    !print '(A,2I3,1X,A,1X,A)', &
    !& "I,N,FAS,MODEL=", &
    !& I,N,TRIM(vMixFas0(I)%Name),TRIM(MM%Name)
    !pause
    !
    ALLOCATE(vXmin(N))
    ALLOCATE(vMuMin(N))
    !
    IF( TRIM(MM%Model)=="MOLECULAR" .AND. MM%NMarg==0 ) THEN
      !-------------------------------------------------analytic minimum
      Multi=  MM%vMulti(1)
      !
      S= Zero
      DO P=1,N
        vXmin(P)= FSafe_Exp(-vFasPur(MM%vIPole(P))%Grt /REAL(Multi))
        S= S + vXmin(P)
      ENDDO
      vXmin(1:N)=  vXmin(1:N) /S
      vMuMin(1:N)= vFasPur(MM%vIPole(1:N))%Grt + Multi*LOG(vXmin(1:N))
      Gmin= SUM(vXmin(1:N) * vMuMin(1:N))
      !------------------------------------------------/analytic minimum
    ELSE
      !---------------------------------------------numerical minimum(s)
      CALL MixModel_Optim_SetParams(TdgK,Pbar,MM)
      ALLOCATE(Mixmodel_Optim_vMu0rt(N))
      ALLOCATE(Mixmodel_Optim_vLPole(N))
      Mixmodel_Optim_vMu0rt(1:N)= vFasPur(MM%vIPole(1:N))%Grt
      Mixmodel_Optim_vLPole(1:N)= .TRUE.
      !
      ALLOCATE(vX(N))
      ALLOCATE(vMu(N))
      Gmin= 1.0D30
      !
      DO J=1,N
        !
        !!! !--- start from compos'n close to end-member J
        !!! vX(J)= One - 1.0D-3
        !!! DO K=1,N
        !!!   IF(K/=J) vX(K)= 1.0D-3/REAL(N-1)
        !!! ENDDO
        !
        vX(1:N)= vMixFas_Xpole_Init(I)%tXPole(J,1:N)
        !
        CALL Optimsolver_Theriak( & !
        !& MixModel_Optim_GetGibbs,      & !
        !& MixModel_Optim_GetPotentials, & !
        & MixModel_Optim_GetMu, & ! INTERFACE
        & FF,                   & ! IN
        & DeltaInit,            & ! IN
        & TolX,                 & ! IN
        & vX,                   & ! INOUT
        & vMu,                  & ! OUT
        & G,                    & ! OUT
        & OkConverge,           & ! OUT
        & its,nCallG)             ! OUT
        !
        IF(iDebug>2) WRITE(91,'(I3,2X,A)') its, TRIM(MM%Name)
        !
        IF(WarmRestart) vMixFas_Xpole_Init(I)%tXPole(J,1:N)= vX(1:N)
        !
        IF(G<Gmin) THEN
          Gmin= G
          vXmin(:)= vX(:)
          vMuMin(:)= vMu(:)
        ENDIF
        !
      END DO
      !PAUSE
      !
      DEALLOCATE(vX)
      DEALLOCATE(vMu)
      !
      DEALLOCATE(Mixmodel_Optim_vMu0rt)
      DEALLOCATE(Mixmodel_Optim_vLPole)
      !--------------------------------------------/numerical minimum(s)
      !
    ENDIF
    !
    IF(iDebug>2) WRITE(91,'(A)') "================================="
    !
    vMixFas0(i)%vLPole(1:N)= .TRUE.
    vMixFas0(i)%vXPole(1:N)= vXmin(1:N)
    vMixFas0(i)%Grt= Gmin
    vMixFas0(i)%vLnAct(1:N)= vMuMin(1:N) -vFasPur(MM%vIPole(1:N))%Grt
    !
    DEALLOCATE(vXmin)
    DEALLOCATE(vMuMin)
    !
  ENDDO
  
  !------------------------------------------------------------log files
  IF(F1>0) THEN
    WRITE(F1,'(/,A,/)') "Mixture_Minimize,EndMember,Xi,Gi"
    DO I=1,SIZE(vMixModel)
      !IF(vMixFas0(I)%Grt > 1.D-3) CYCLE
      MM= vMixModel(I)
      WRITE(F1,'(2A)') "MODEL= ",MM%Name
      DO K=1,MM%NPole
        WRITE(F1,'(A,1X,G15.6)') &
        & MM%vNamPole(K),        &
        & vMixFas0(I)%vXPole(K)
      ENDDO
      WRITE(F1,'(A15,G15.6)') "G Mimim= ",vMixFas0(I)%Grt/Ln10
      WRITE(F1,*)
    ENDDO
    WRITE(F1,*)
  ENDIF
  !
  IF(F2>0) THEN
    DO I=1,SIZE(vMixModel)
    !IF(vMixFas0(I)%Grt < Zero) THEN
      MM= vMixModel(I)
      !
      WRITE(F2,'(A,A1)',ADVANCE="NO") MM%Name,T_
      DO K=1,MM%NPole
        WRITE(F2,'(G15.6,A1)',ADVANCE="NO") vMixFas0(I)%vXPole(K),T_
      ENDDO
      !DO K=1,MM%NPole
      !  WRITE(F2,'(G15.6,A1)',ADVANCE="NO") vFasPur(MM%vIPole(K))%Grt,T_
      !ENDDO
      WRITE(F2,'(G15.6,A1)',ADVANCE="NO") vMixFas0(I)%Grt,T_
      !
    !ENDIF
    ENDDO
    WRITE(F2,*)
  ENDIF
  !-----------------------------------------------------------/log files
  
  RETURN
END SUBROUTINE Mixture_Minimize

SUBROUTINE Speciation
END SUBROUTINE Speciation

SUBROUTINE Write_Log_Entete(vFas)
  USE M_IoTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirOut
  !
  TYPE(T_Phase),INTENT(IN):: vFas(:)
  !
  CALL GetUnit(F1)
  OPEN(F1,FILE=TRIM(DirOut)//"_gibbs.log")
  WRITE(F1,'(A)') "g of all phases, relative to current assemblage"
  !
  CALL GetUnit(F2)
  OPEN(F2,FILE=TRIM(DirOut)//"_mixcomp.log")
  WRITE(F2,'(A)')   "properties of mixtures with min(g)<0 at each iteration"
  WRITE(F2,'(A,/)') "mix.model.name, Xi, Gmin"
  !
END SUBROUTINE Write_Log_Entete

ENDMODULE M_Simplex_Theriak
