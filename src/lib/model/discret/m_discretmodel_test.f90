MODULE M_DiscretModel_Test
  !.object describing "discretization" of a mixture phase
  !.to an array of pure phases
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_,Pause_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: DiscretModel_Test
  !
CONTAINS

SUBROUTINE Stoikio_Test
  USE M_IoTools,ONLY: GetUnit
  USE M_T_Species,ONLY: T_Species
  USE M_Global_Vars,ONLY: vEle,vFas,vSpc,vMixModel
  USE M_Global_Vars,ONLY: vDiscretModel,vDiscretParam
  !
  INTEGER:: f
  INTEGER:: iSp,iEl,iDis
  TYPE(T_Species):: S
  !
  CALL GetUnit(f)
  OPEN(f,file="debug_discretmodel_stoik.log")
  !
  WRITE(f,'(2(A,A1))',ADVANCE="no") "name",T_,"discri",T_
  WRITE(f,'(3(A,A1))',ADVANCE="no") "disc%I",T_,"disc%iEl",T_,"disc%iDis",T_
  DO iEl=1,SIZE(vEle)
    WRITE(f,'(A,A1)',ADVANCE="no") vEle(iEl)%NamEl,T_
  ENDDO
  WRITE(f,'(A)') "div"
  !
  DO iSp=1,SIZE(vSpc)
    !
    S= vSpc(iSp)
    iDis= S%iDiscret
    !
    IF(iDis==0) CYCLE
    !
    WRITE(f,'(A,A1)',ADVANCE="no") S%NamSp,T_
    WRITE(f,'(A,A1)',ADVANCE="no") &
    & vDiscretModel(vDiscretParam(iDis)%iModel)%Name,T_
    !
    WRITE(f,'(3(I3,A1))',ADVANCE="no") &
    & vDiscretParam(iDis)%i,T_, &
    & vDiscretParam(iDis)%j,T_, &
    & vDiscretParam(iDis)%k,T_
    !
    DO iEl=1,SIZE(vEle)
      WRITE(f,'(I3,A1)',ADVANCE="no")  S%vStoikio(iEl),T_
    ENDDO
    WRITE(f,'(I3,A1)') S%vStoikio(0)
    !
  ENDDO
  !
  CLOSE(f)
  !
  IF(iDebug>1) THEN
    PRINT *,">> results in debug_discretmodel_stoik.log !!!"
    CALL Pause_
  ENDIF
ENDSUBROUTINE Stoikio_Test

SUBROUTINE DiscretModel_Test(TdgK,Pbar)
  USE M_Global_Vars,ONLY: vEle,vFas,vSpc,vMixModel
  USE M_Global_Vars,ONLY: vDiscretModel,vDiscretParam
  USE M_DiscretModel_Tools
  !
  REAL(dp),INTENT(IN):: TdgK,Pbar
  !
  INTEGER:: i
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DiscretModel_Test"
  !
  CALL DiscretModel_Test_Init(TdgK,Pbar)
  !
  IF(iDebug>1) THEN
    !
    PRINT *,"<============================================== SPECIES =="
    DO i=1,SIZE(vSpc)
      PRINT *,vSpc(i)%NamSp
    ENDDO
    CALL Pause_
    !
    PRINT *,"<=========================================== MIX.MODELS =="
    DO i=1,SIZE(vMixModel)
      PRINT *,vMixModel(i)%Name
    ENDDO
    CALL Pause_
    !
    PRINT *,"<======================================= DISCRET.MODELS =="
    DO i=1,SIZE(vDiscretModel)
      PRINT *,vDiscretModel(i)%Name
    ENDDO
    CALL Pause_
    !
    PRINT *,"<=============================================== PHASES =="
    DO i=1,SIZE(vFas)
      PRINT *,vFas(i)%NamFs
    ENDDO
    CALL Pause_
    !
    PRINT *,"<======================================= DISCRET.PARAMS =="
    DO i=1,SIZE(vDiscretParam)
      PRINT *, &
      & vDiscretModel(vDiscretParam(i)%iModel)%Name, &
      & vDiscretParam(i)%I, &
      & vDiscretParam(i)%J, &
      & vDiscretParam(i)%K
    ENDDO
    CALL Pause_
    !
  ENDIF
  !
  CALL DiscretSpecies_Stoikio_Calc( &
  & vEle,          & !IN
  & vMixModel,     & !IN
  & vDiscretModel, & !IN
  & vDiscretParam, & !IN
  & vSpc)            !INOUT with updated stoichio's for discretized species
  !
  CALL Stoikio_Test
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DiscretModel_Test"
  !
ENDSUBROUTINE DiscretModel_Test

SUBROUTINE DiscretModel_Test_Init(TdgK,Pbar)
  USE M_T_Species,   ONLY: T_Species,Species_Stoikio_Calc
  USE M_Dtb_Const,   ONLY: T_CK,Tref,Pref
  USE M_Dtb_Read,    ONLY: Dtb_Read_Species
  USE M_Dtb_Calc,    ONLY: SpeciesDtb_ToSpecies
  USE M_Element_Read,ONLY: Element_Read_Redox,Element_Read_Entropy
  USE M_Global_Tools,ONLY: Global_TP_Update
  USE M_DiscretModel_Read
  USE M_DiscretModel_Tools
  USE M_Global_Alloc
  !
  USE M_Global_Vars, ONLY: vEle,vSpc,vSpcDtb,vFas
  USE M_Global_Vars, ONLY: vMixModel,vMixFas
  USE M_Global_Vars, ONLY: vDiscretModel,vDiscretParam
  USE M_Dtb_Vars,    ONLY: DtbLogK_vTPCond
  !
  REAL(dp),INTENT(IN):: TdgK,Pbar
  !
  INTEGER:: N
  LOGICAL:: fOk
  TYPE(T_Species),ALLOCATABLE:: vSpcTmp(:)
  !
  CALL Elements_Alloc_forDtb !-> builds vEle
  !
  CALL Element_Read_Redox(vEle)
  CALL Element_Read_Entropy(vEle)
  !
  CALL Dtb_Read_Species(vEle) !-> read databases -> build vDtbXxx
  !
  CALL SpeciesDtb_Alloc       !-> from vDtb*, build vSpcDtb
  !
  N= SIZE(vSpcDtb)
  IF(N<1) RETURN
  !
  DEALLOCATE(vSpc)  ;  ALLOCATE(vSpc(1:N))
  !
  !! CALL Species_Alloc_forDtb(vEle) !-> read databases -> vSpc,tFormula
  CALL SpeciesDtb_ToSpecies(vSpcDtb,vSpc)
  !
  CALL Species_Stoikio_Calc(vEle,vSpc,fOk)
  !
  CALL MixModels_Alloc(vSpc) !-> build MixModel DATAbase vMixModel
  !
  CALL DiscretModel_Read(vMixModel) !-> allocate vDiscretModel
  !
  N= SIZE(vDiscretModel)
  IF(N>0) THEN
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
  ENDIF
  !
  CALL MixPhases_Alloc(vSpc,vMixModel) !read phase compositions, build vMixFas
  !
  CALL Phases_Alloc(vSpc,vMixFas) !-> build vFas
  !
  CALL Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel,vMixFas,vFas) !inout
  !
ENDSUBROUTINE DiscretModel_Test_Init

SUBROUTINE Species_Append(vSpcAdd)
!--
!-- append vSpcAdd to current vSpc -> produce a new vSpc
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
ENDMODULE M_DiscretModel_Test
