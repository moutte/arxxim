MODULE M_Global_Build
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Global_Build
  PUBLIC:: Global_Build_New

CONTAINS

SUBROUTINE Global_Build
!--
!-- build "global" vEle and vSpc
!-- (all available elements and species, not restricted to those used in the run)
!--
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc  !,Stop_,Warning_
  !
  USE M_T_Species,   ONLY: T_Species,Species_Stoikio_Calc
  USE M_Dtb_Const,   ONLY: T_CK,Tref,Pref
  USE M_Dtb_Read,    ONLY: Dtb_Read_Species
  USE M_Dtb_Calc,    ONLY: SpeciesDtb_ToSpecies
  USE M_Element_Read,ONLY: Element_Read_Redox,Element_Read_Entropy
  USE M_SpeciesDtb_Read,ONLY: Species_Read_AquSize,Species_Write_AquSize
  USE M_Global_Tools,ONLY: Global_TP_Update,Global_Species_Select
  USE M_Global_Alloc
  !
  !--global variables--
  USE M_Global_Vars, ONLY: vEle,vSpc,vSpcDtb
  USE M_Global_Vars, ONLY: vMixModel,vMixFas,vFas
  USE M_Global_Vars, ONLY: vDiscretModel,vDiscretParam
  !
  INTEGER :: iTP, N
  LOGICAL :: fOk
  REAL(dp):: TdgK,Pbar
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Global_Build"
  !
  CALL Elements_Alloc_forDtb !-> builds vEle
  !
  CALL Element_Read_Redox(vEle)
  !CALL Element_Read_Entropy(vEle)
  !
  CALL Dtb_Read_Species(vEle) !-> read databases -> build vDtbXxx
  !
  CALL SpeciesDtb_Alloc       !-> from vDtb*, build vSpcDtb
  !
  ! write(11,'(A)') "<=======vSpcDtb=1==="
  ! do n=1,size(vSpcDtb)
  !   write(11,'(I3,1X,I3,1X,A)') n, vSpcDtb(n)%Indx, vSpcDtb(n)%DtbModel
  ! end do
  ! write(11,'(A)') "===================>"
  !
  N= SIZE(vSpcDtb)
  IF(N<1) RETURN
  !
  DEALLOCATE(vSpc)  ;  ALLOCATE(vSpc(N))
  !
  CALL SpeciesDtb_ToSpecies(vSpcDtb,vSpc)
  !
  ! write(11,'(A)') "<=======species=2==="
  ! do n=1,size(vSpc)
  !   write(11,'(I3,1X,A)') n, vSpc(n)%NamSp
  ! end do
  ! write(11,'(A)') "===================>"
  !
  CALL Global_Species_Select
  !
  CALL Species_Stoikio_Calc(vEle,vSpc,fOk)
  !
  CALL Species_Read_AquSize(vSpc)
  ! CALL Species_Write_AquSize(vSpc)
  !
  CALL vSpecies_Rename(vSpc)
  !
  CALL MixModels_Alloc(vSpc) !-> build MixModel database vMixModel
  !
  CALL MixPhases_Alloc(vSpc,vMixModel) !READ phase compositions, build vMixFas
  !
  CALL Phases_Alloc(vSpc,vMixFas) !-> build vFas
  !
  CALL Condition_Read(iTP)
  !
  TdgK= Tref
  Pbar= Pref
  !
  CALL Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel,vMixFas,vFas) !inout
  !
  CALL Species_Write_AquSize(vSpc)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Global_Build"
  !
END SUBROUTINE Global_Build

SUBROUTINE Global_Build_New
!--
!-- build "global" vEle and vSpc
!-- (all available elements and species,
!-- not restricted to those used in the run)
!--
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,Stop_  !,Warning_
  !
  USE M_T_Species,   ONLY: T_Species,Species_Stoikio_Calc
  USE M_T_SolModel,  ONLY: T_SolModel
  !
  USE M_Dtb_Const,   ONLY: T_CK,Tref,Pref
  USE M_Dtb_Read,    ONLY: Dtb_Read_Species
  USE M_Dtb_Calc,    ONLY: SpeciesDtb_ToSpecies
  !
  USE M_Element_Read,ONLY: Element_Read_Redox,Element_Read_Entropy
  USE M_SpeciesDtb_Read,ONLY: Species_Read_AquSize,Species_Write_AquSize
  USE M_Global_Tools,ONLY: Global_TP_Update,Global_Species_Select
  USE M_Global_Alloc
  !
  !~ USE M_SolModel_Alloc
  USE M_SolModel_Read
  USE M_SolPhase_Read
  USE M_T_SolModel,   ONLY: SolModel_Spc_Init
  USE M_T_SolPhase,   ONLY: SolPhase_Init
  !
  USE M_DiscretModel_Read
  USE M_DiscretModel_Tools
  !
  ! !--database variables--
  ! USE M_Solmodel_Vars, ONLY: Ok_Rho,Ok_Eps,Ok_DHA,Ok_DHB,Ok_BDot
  ! USE M_Solmodel_Vars, ONLY: Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl,T_Spline
  !
  !--global variables--
  USE M_Global_Vars, ONLY: vEle,vSpc,vSpcDtb
  USE M_Global_Vars, ONLY: vMixModel,vMixFas,vFas
  USE M_Global_Vars, ONLY: vSolModel,vSolFas
  USE M_Global_Vars, ONLY: vDiscretModel,vDiscretParam
  !
  INTEGER :: iTP, N
  LOGICAL :: fOk
  REAL(dp):: TdgK,Pbar
  LOGICAL:: Ok
  CHARACTER(LEN=80):: Msg
  TYPE(T_Species),ALLOCATABLE:: vSpcTmp(:)
  ! TYPE(T_SolModel):: SolModel
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Global_Build"
  !
  CALL Elements_Alloc_forDtb !-> builds vEle
  !
  CALL Element_Read_Redox(vEle)
  !CALL Element_Read_Entropy(vEle)
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
  CALL SpeciesDtb_ToSpecies(vSpcDtb,vSpc)
  !
  CALL Global_Species_Select
  !
  CALL Species_Stoikio_Calc(vEle,vSpc,fOk)
  !
  CALL Species_Read_AquSize(vSpc)
  !CALL Species_Write_AquSize(vSpc)
  !
  CALL vSpecies_Rename(vSpc)
  !
  CALL MixModels_Alloc(vSpc) !-> build MixModel database vMixModel
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
  !----------------------------------------------------- build vSolModel
  IF(COUNT(vSpc(:)%Typ=="AQU") >0) THEN

    ! N= 0
    ! IF(DtbFormat=="LOGKTBL") N= DtbLogK_Dim

    ! ALLOCATE(vTdgC(N))
    ! IF(N>0) vTdgC(1:N)= DtbLogK_vTPCond(1:N)%TdgC

    ! ! even for HSV base, Solmodel_Read must be called,
    ! ! for reading activity model
    ! CALL Solmodel_Solvent_Read( &
    ! & N,vTdgC,SolModel, &
    ! & Ok_Rho, Ok_Eps, Ok_DHA, Ok_DHB, Ok_BDot, &
    ! & Rho_Spl,Eps_Spl,DHA_Spl,DHB_Spl,BDot_Spl)

    ! IF(SolModel%ActModel=="PITZER" .OR. &
    ! &  SolModel%ActModel=="SIT"  ) THEN
    !   CALL Pitzer_Dtb_Init(vSpc)
    !   CALL Pitzer_Dtb_TPUpdate(TdgK,Pbar)
    !   IF(iDebug>2) CALL Pitzer_Dtb_TPtest !!!"WRK"!!!
    ! ENDIF

    !!--- initialize indexes of the sol'model
    ! CALL SolModel_Spc_Init("H2O",vSpc,SolModel,Ok,Msg)

    ! CALL SolModel_Alloc ! allocate vSolModel
    ! vSolModel(1)= SolModel

    CALL SolModel_Read(vSpc,Ok,Msg)
    IF(.NOT. Ok) &
    & CALL Stop_("SolModel_Read==> "//TRIM(Msg))

    CALL SolPhase_Read(vSpc,vSolModel,Ok,Msg)
    IF(.NOT. Ok) &
    & CALL Stop_("SolPhase_Read==> "//TRIM(Msg))
    ! CALL SolPhase_Alloc ! allocate vSolFas

    ! CALL SolPhase_Init( & !
    ! & 1,          & !IN, model index
    ! & vSolModel,  & !IN, base of solution models
    ! & vSolFas(1), & !INOUT, solution phase
    ! & Ok,Msg)       !OUT
 
    ! CALL SolPhase_Model_Init(SolModel%Model,vSolModel,vSolFas(1),Ok,Msg)
    ! IF(.NOT. Ok) CALL Stop_(TRIM(Msg))


  ENDIF
  !-----------------------------------------------------/build vSolModel
  !
  CALL MixPhases_Alloc(vSpc,vMixModel) !read phase compositions, build vMixFas
  !
  CALL Phases_Alloc_New(vSpc,vMixFas,vSolFas)
  !
  CALL Condition_Read(iTP)
  !
  TdgK= Tref
  Pbar= Pref
  !
  CALL Global_TP_Update( &
  & TdgK,Pbar,vSpcDtb,vDiscretModel,vDiscretParam, & !in
  & vSpc,vMixModel,vMixFas,vFas) !inout
  !
  CALL Species_Write_AquSize(vSpc)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Global_Build"
  !
END SUBROUTINE Global_Build_New

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

SUBROUTINE vSpecies_Rename(vSpc)
  USE M_T_Species,ONLY: T_Species,Species_Rename
  TYPE(T_Species),INTENT(INOUT):: vSpc(:)
  !
  INTEGER:: I
  !
  DO I=1,SIZE(vSpc)
    CALL Species_Rename(vSpc(I)%NamSp)
  END DO
  !
END SUBROUTINE vSpecies_Rename

SUBROUTINE Condition_Read(iTP)
!--
!-- read the block CONDITIONS
!--
  USE M_Trace,ONLY: fHtm,iDebug,fTrc
  USE M_Files,ONLY: NamFInn,cTitle,File_Path
  USE M_Files,ONLY: DirOut,DirLog,DirDtbLog,DirDtbOut
  USE M_Files,ONLY: Files_Index_Close,Files_Index_Open
  USE M_IoTools
  !
  INTEGER,INTENT(OUT):: iTP
  !
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W
  LOGICAL:: EoL,Ok
  INTEGER:: f,ios,N
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Condition_Read"
  !
  iTP= 0
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))

  Ok=.FALSE.

  DoFile: DO
    !
    READ(f,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)

    CASE("ENDINPUT")
      EXIT  DoFile
    !
    CASE("CONDITIONS")
      Ok=.TRUE.
      !
      LoopRead: DO
        !
        READ(f,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE LoopRead !skip comment lines
        CALL AppendToEnd(L,W,EoL)
        !
        SELECT CASE(W)
        !
        CASE("ENDINPUT")            ; EXIT DoFile
        CASE("END","ENDCONDITIONS") ; EXIT LoopRead
        !
        CASE("TP.POINT","TPPOINT")
          !! CALL Stop_(TRIM(W)//" is Obsolete !! Give TdgC, Pbar in SYSTEM")
          CALL LinToWrd(L,W,EoL); CALL WrdToInt(W,iTP)
        CASE("TITLE")
          cTitle=TRIM(L)
        CASE("OUTPUT")
          CALL LinToWrd(L,W,EoL,"NO"); DirOut=File_Path(TRIM(W))
        CASE("DIRLOG")
          CALL LinToWrd(L,W,EoL,"NO"); DirLog=File_Path(TRIM(W))
        CASE("DTBDIROUT")
          CALL LinToWrd(L,W,Eol,"NO"); DirDtbOut=File_Path(TRIM(W))
        CASE("DTBDIRLOG")
          CALL LinToWrd(L,W,Eol,"NO"); DirDtbLog=File_Path(TRIM(W))
        CASE("DEBUG")
          CALL LinToWrd(L,W,EoL); CALL WrdToInt(W,N)
          IF(N>0 .AND. iDebug==0) iDebug=N
        !
        END SELECT
        !
      ENDDO LoopRead
    !
    END SELECT
    !
  ENDDO DoFile

  CLOSE(f)

  IF(.NOT.Ok) THEN
    IF(iDebug>0) WRITE(fTrc,'(A)') "block CONDITIONS not found ???!!!"
    IF(iDebug>2) PRINT '(A)',"block CONDITIONS not found ???!!!"
  ENDIF

  IF(iDebug>0) THEN
    IF(fHtm>0) CALL Files_Index_Close(fHtm)
    CALL Files_Index_Open(fHtm)
  ENDIF

  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Condition_Read"

END SUBROUTINE Condition_Read

ENDMODULE M_Global_Build
