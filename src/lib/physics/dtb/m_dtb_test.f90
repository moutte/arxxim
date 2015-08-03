MODULE M_Dtb_Test
!--
!-- routines or testing databases, solution models, etc
!--
  USE M_Kinds
  USE M_T_Tpcond,ONLY: T_TPCond
  USE M_Trace,   ONLY: iDebug,fTrc,fHtm,Stop_,T_,Pause_,Warning_

  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Dtb_Tabulate
  PUBLIC:: Dtb_Tabulate_ForSystem
  PUBLIC:: Dtb_Test_Fluid
  PUBLIC:: Dtb_Test_H2OHkf
  PUBLIC:: Dtb_Test_AquHkf
  PUBLIC:: Dtb_Test_Species
  PUBLIC:: Dtb_Tabulate_PSatH2O
  PUBLIC:: Dtb_Test_EQ36
  !
CONTAINS

SUBROUTINE Dtb_Tabulate(sCode)
!-----------------------------------------------------------------------
!--tabulate the species logK's along the vTPpath
!--retrieved from the TP.PATH block
!-----------------------------------------------------------------------
  USE M_Dtb_Const,   ONLY: T_CK
  USE M_T_SolModel,  ONLY: T_SolModel,T_SolModelDat
  USE M_SolModel_Tools,ONLY: SolModel_TP_Update
  USE M_Dtb_Calc,    ONLY: DtbSpc_GrtTable_Build
  USE M_TPcond_Read, ONLY: TPpath_Read
  USE M_Path_Vars,   ONLY: vTPpath
  !
  USE M_Global_Vars,ONLY: vSpcDtb,vSpc
  !
  CHARACTER(LEN=*),INTENT(IN):: sCode
  !
  REAL(dp),ALLOCATABLE:: tGrt(:,:)
  TYPE(T_SolModelDat),ALLOCATABLE:: vSolvDat(:)
  TYPE(T_SolModel):: Slv
  INTEGER :: i,N,M
  REAL(dp):: TdgK,Pbar

  !--read the TP path
  CALL TPpath_Read
  !
  N= SIZE(vTPpath)
  M= SIZE(vSpc)
  !
  ALLOCATE(vSolvDat(1:N))
  ALLOCATE(tGrt(1:SIZE(vSpc),1:N))
  tGrt= Zero
  !
  !----------compute the solvent properties for the series of T,P points
  DO i=1,N

    TdgK= vTPpath(i)%TdgC +T_CK
    Pbar= vTPpath(i)%Pbar

    CALL SolModel_TP_Update(TdgK,Pbar,Slv)

    vSolvDat(i)= Slv%Dat

  ENDDO
  !-----------/
  !
  !---------compute the species' thermodyn' properties- build table tGrt
  CALL DtbSpc_GrtTable_Build(vTPpath,vSpcDtb,vSpc,tGrt)
  !--------/
  !
  !--- write tGrt to a file --
  CALL Dtb_GrtTable_Tabulate(TRIM(sCode),vTPpath,vSolvDat,vSpcDtb,vSpc,tGrt)
  !
  DEALLOCATE(tGrt)
  DEALLOCATE(vSolvDat)
  DEALLOCATE(vTPpath)
  !
ENDSUBROUTINE Dtb_Tabulate

!-----------------------------------------------------------------------
!> write a logK database for current run,
!! setting logK=0 for primary species
SUBROUTINE Dtb_Tabulate_ForSystem( &
& vTPpath,  &
& vCpn,vSpc,vSpcDtb, &
& vLCi,vLAx, &
& vOrdPr, &
& tNuSp)
  USE M_IoTools,    ONLY: GetUnit, OutStrVec
  USE M_Files,      ONLY: DirOut,Files_Index_Write
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_T_Tpcond,   ONLY: T_TPCond
  USE M_T_SolModel,  ONLY: T_SolModel,T_SolModelDat
  USE M_SolModel_Tools,ONLY: SolModel_TP_Update
  USE M_Numeric_Const,ONLY: Ln10
  USE M_T_Component,ONLY: T_Component
  USE M_T_Species,  ONLY: T_Species,T_SpeciesDtb
  !
  USE M_Dtb_Calc,   ONLY: DtbSpc_GrtTable_Build
  !
  !< NEW 2011-03, for discretized species
  USE M_Global_Vars,ONLY: vMixModel,vDiscretModel,vDiscretParam
  USE M_T_MixModel, ONLY: MixModel_Margul_Wg_Calc
  USE M_DiscretModel_Tools,ONLY: DiscretSpecies_TP_Update
  !</
  !
  TYPE(T_TPCond),    INTENT(IN):: vTPpath(:)
  TYPE(T_Component), INTENT(IN):: vCpn(:)
  TYPE(T_Species),   INTENT(IN):: vSpc(:)
  TYPE(T_SpeciesDtb),INTENT(IN):: vSpcDtb(:)
  LOGICAL,           INTENT(IN):: vLCi(:),vLAx(:)
  INTEGER,           INTENT(IN):: vOrdPr(:)
  REAL(dp),          INTENT(IN):: tNuSp(:,:)
  !
  INTEGER :: iTP,iSp,I,J,N,nCp,F,iMM
  REAL(dp):: X,Tk,Pb !, T, P
  CHARACTER(LEN=15):: Cc
  TYPE(T_SolModel):: Slv
  TYPE(T_Species):: S
  TYPE(T_SolModelDat),ALLOCATABLE:: vSolvDat(:)
  TYPE(T_Species),   ALLOCATABLE:: vSpcTmp(:)
  REAL(dp),          ALLOCATABLE:: tGrt(:,:)
  !
  N= SIZE(vTPpath)
  !
  ALLOCATE(vSolvDat(1:N))
  ALLOCATE(tGrt(1:SIZE(vSpc),1:N))
  tGrt= Zero
  !
  DO iTP=1,N
    !! CALL SolventDat_Default(vSolvDat(iTP))
    Tk= vTPpath(iTP)%TdgC +T_CK
    Pb= vTPpath(iTP)%Pbar
    CALL SolModel_TP_Update(Tk,Pb,Slv)
    vSolvDat(iTP)= Slv%Dat
  ENDDO
  !
  CALL DtbSpc_GrtTable_Build(vTPpath,vSpcDtb,vSpc,tGrt)

  !---NEW 2011-03
  !-------------------------------update species built by discretization
  IF(SIZE(vDiscretModel)>0) THEN

    ALLOCATE(vSpcTmp(SIZE(vSpc)))  ;  vSpcTmp(:)= vSpc(:)

    DO iTP=1,N

      Tk= vTPpath(iTP)%TdgC +T_CK
      Pb= vTPpath(iTP)%Pbar

      !---------------update T,P dependent parameters in mixing models--
      DO iMM=1,SIZE(vMixModel)
        IF(vMixModel(iMM)%NMarg>0) &
        & CALL MixModel_Margul_Wg_Calc( &
        & Tk,Pb,       &  !in
        & vMixModel(iMM)) !out
      ENDDO
      !--------------/

      DO iSp= 1,SIZE(vSpcTmp)
        IF(vSpc(iSp)%iDtb>0) &
        & vSpcTmp(iSp)%G0rt= tGrt(iSp,iTP)
      ENDDO

      CALL DiscretSpecies_TP_Update( &
      & vMixModel,     & !IN
      & Tk,Pb,         & !IN
      & vDiscretModel, & !IN
      & vDiscretParam, & !IN
      & vSpcTmp)         !INOUT

      DO iSp= 1,SIZE(vSpcTmp)
        IF(vSpc(iSp)%iDiscret/=0) &
        & tGrt(iSp,iTP)= vSpcTmp(iSp)%G0rt
      ENDDO

    ENDDO

    DEALLOCATE(vSpcTmp)

  ENDIF
  !------------------------------/update species built by discretization
  !
  nCp=SIZE(vCpn)
  !
  !--------------------------------logK from current prim'species set--/
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_logk_prim.dtb")
  !
  !--------------------------------------------------------index in html
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_logk_prim.dtb",&
  & "table of logK from prim'species, for current run")
  !----------------------------------------------------/
  !
  !-------------------------------------------------------TP.TABLE block
  WRITE(F,'(/,A)') "TP.TABLE"
  CALL OutStrVec(F,vTPpath(1:N)%TdgC,Opt_S="TdgC")
  CALL OutStrVec(F,vTPpath(1:N)%Pbar,Opt_S="Pbar")
  WRITE(F,'(A,/)') "END TP.TABLE"
  !---------------------------------------------------/
  !
  !--------------------------------------------------------SOLVENT block
  WRITE(F,'(/,A)') "SOLVENT"
  WRITE(F,'(A)') "NAME H2O"
  CALL OutStrVec(F,vSolvDat(1:N)%Rho, Opt_S="Rho")
  CALL OutStrVec(F,vSolvDat(1:N)%Eps, Opt_S="Eps")
  CALL OutStrVec(F,vSolvDat(1:N)%dhA, Opt_S="DHA")
  CALL OutStrVec(F,vSolvDat(1:N)%dhB, Opt_S="DHB")
  CALL OutStrVec(F,vSolvDat(1:N)%BDOt,Opt_S="BDOT")
  WRITE(F,'(A,/)') "END SOLVENT"
  !----------------------------------------------------/
  !
  !------------------------------------------------------- SPECIES block
  WRITE(F,'(/,A)') "SPECIES" !
  !-------------------------------------first, write the primary species
  DoPrim: DO I=1,SIZE(vSpc)
    !
    IF(.NOT. (vLCi(I) .OR. vLAx(I))) CYCLE DoPrim
    !
    S=vSpc(I)
    IF(S%iDtb>0) THEN
      Cc= TRIM(vSpcDtb(S%iDtb)%DtbTrace)
    ELSE
      Cc= "_"
    ENDIF
    !
    WRITE(F, '(A3,A1, A15,A1, A23,A1, A63,A1)',ADVANCE='NO') &
    & S%Typ,    T_, &
    & Cc,       T_, &
    & S%NamSp,  T_, &
    & S%Formula,T_
    
    !WRITE(F,'(4(A,A1))',ADVANCE='NO') &
    !& TRIM(S%Typ),    T_, &
    !& TRIM(Cc),       T_, &
    !& TRIM(S%NamSp),   T_, &
    !& TRIM(S%Formula),T_
    
    IF(S%Typ=="AQU") THEN
      !for aqu'species, write size parameter
      WRITE(F, '(F15.1,A1)',ADVANCE='NO') S%AquSize, T_
    ELSE
      !for non aqueous species, write density
      WRITE(F, '(G15.8,A1)',ADVANCE='NO') S%WeitKg / S%V0, T_
    ENDIF
    
    WRITE(F,'(*(A15,A1))') ("0.000", T_,J=1,N)
    
  ENDDO DoPrim
  !------------------------------------/first, write the primary species
  !
  !------------------------------------then, write the secondary species
  DoSecond: DO I=1,SIZE(vSpc)
    !
    IF(vLCi(I) .OR. vLAx(I)) CYCLE DoSecond
    !
    S=vSpc(I)
    IF(S%iDtb>0) THEN
      Cc= TRIM(vSpcDtb(S%iDtb)%DtbTrace)
    ELSE
      Cc= "_"
    ENDIF
    !
    WRITE(F, '(A3,A1, A15,A1, A23,A1, A63,A1)',ADVANCE='NO') &
    & S%Typ,    T_, &
    & Cc,       T_, &
    & S%NamSp,   T_, &
    & S%Formula,T_
    !
    IF(S%Typ=="AQU") THEN
      !-- for aqu'species, write size parameter
      WRITE(F, '(F15.1,A1)',ADVANCE='NO') S%AquSize, T_
    ELSE
      !-- for non aqueous species, WRITE density
      IF(S%V0<1.D-9) THEN
        WRITE(F,'(A,A1)',ADVANCE='NO') "_?_",T_
      ELSE
        WRITE(F, '(G15.8,A1)',ADVANCE='NO') S%WeitKg / S%V0, T_
      ENDIF
    ENDIF
    !
    DO J=1,N
      ! X=  vSpc(vOrdAs(I))%G0rt &
      ! & - DOT_PRODUCT(tNuAs(I,1:nCp),vSpc(vOrdPr(1:nCp))%G0rt)
      X=  tGrt(I,J) &
      & - DOT_PRODUCT(tNuSp(I,1:nCp),tGrt(vOrdPr(1:nCp),J))
      WRITE(F,'(G15.8,A1)',ADVANCE='NO') -X /Ln10, T_
    ENDDO
    !
    WRITE(F,*)
    !
  ENDDO DoSecond
  !-----------------------------------/then, write the secondary species
  !
  WRITE(F,'(A,/)') "END SPECIES"
  !-------------------------------------------------------/SPECIES block
  !
  CLOSE(F)
  !
  IF(iDebug>0) PRINT '(A)',"Table of LogK in "//TRIM(DirOut)//"_logk.dtb"
  !------------------------------/ logK from current prim'species set --
  !
  !---------------------------------------------------------------------
  !-------------------------------------------------- logK from database
  !-------------------------without refering to current prim'species set
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_logk_elem.dtb")
  !
  !-------------------------------------------------------index in html
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_logk_elem.dtb",&
  & "table of logK from elements, for current run")
  !----------------------------------------------------/
  !
  !-------------------------------------------------------TP.TABLE block
  WRITE(F,'(/,A)') "TP.TABLE"
  CALL OutStrVec(F,vTPpath(1:N)%TdgC,Opt_S="TdgC")
  CALL OutStrVec(F,vTPpath(1:N)%Pbar,Opt_S="Pbar")
  WRITE(F,'(A,/)') "END TP.TABLE"
  !------------------------------------------------------/
  !
  !--------------------------------------------------------SOLVENT block
  WRITE(F,'(/,A)') "SOLVENT"
  WRITE(F,'(A)') "NAME H2O"
  CALL OutStrVec(F,vSolvDat(1:N)%Rho, Opt_S="Rho")
  CALL OutStrVec(F,vSolvDat(1:N)%Eps, Opt_S="Eps")
  CALL OutStrVec(F,vSolvDat(1:N)%dhA, Opt_S="DHA")
  CALL OutStrVec(F,vSolvDat(1:N)%dhB, Opt_S="DHB")
  CALL OutStrVec(F,vSolvDat(1:N)%BDOt,Opt_S="BDOT")
  WRITE(F,'(A,/)') "END SOLVENT"
  !-------------------------------------------------------/SOLVENT block
  !
  !--------------------------------------------------------SPECIES block
  WRITE(F,'(/,A)') "SPECIES" !
  !-------------------------------------first, write the primary species
  DoPrim2: DO I=1,SIZE(vSpc)
    !
    IF(.NOT. (vLCi(I) .OR. vLAx(I))) CYCLE DoPrim2
    !
    S=vSpc(I)
    IF(S%iDtb>0) THEN
      Cc= TRIM(vSpcDtb(S%iDtb)%DtbTrace)
    ELSE
      Cc= "_"
    ENDIF
    !
    WRITE(F, '(A3,A1, A15,A1, A23,A1, A63,A1)',ADVANCE='NO') &
    & S%Typ,    T_, &
    & Cc,       T_, &
    & S%NamSp,   T_, &
    & S%Formula,T_
    !WRITE(F,'(4(A,A1))',ADVANCE='NO') &
    !& TRIM(S%Typ),    T_, &
    !& TRIM(Cc),       T_, &
    !& TRIM(S%NamSp),   T_, &
    !& TRIM(S%Formula),T_
    IF(S%Typ=="AQU") THEN
      !for aqu'species, write size parameter
      WRITE(F, '(F15.1,A1)',ADVANCE='NO') S%AquSize, T_
    ELSE
      !for non aqueous species, write density
      WRITE(F, '(G15.8,A1)',ADVANCE='NO') S%WeitKg / S%V0, T_
    ENDIF
    DO J=1,N
      X=  tGrt(I,J)
      WRITE(F,'(G15.8,A1)',ADVANCE='NO') -X /Ln10, T_
    ENDDO
    !
    WRITE(F,*)
  ENDDO DoPrim2
  !------------------------------------/first, write the primary species
  !
  !------------------------------------then, write the secondary species
  DoSecond2: DO I=1,SIZE(vSpc)
    !
    IF(vLCi(I) .OR. vLAx(I)) CYCLE DoSecond2
    !
    S=vSpc(I)
    IF(S%iDtb>0) THEN
      Cc= TRIM(vSpcDtb(S%iDtb)%DtbTrace)
    ELSE
      Cc= "_"
    ENDIF
    !
    WRITE(F, '(A3,A1, A15,A1, A23,A1, A63,A1)',ADVANCE='NO') &
    & S%Typ,    T_, &
    & Cc,       T_, &
    & S%NamSp,  T_, &
    & S%Formula,T_
    !
    IF(S%Typ=="AQU") THEN
      !-- for aqu'species, WRITE SIZE PARAMETER
      WRITE(F, '(F15.1,A1)',ADVANCE='NO') S%AquSize, T_
    ELSE
      !-- for non aqueous species, WRITE density
      IF(S%V0<1.D-9) THEN
        WRITE(F,'(A,A1)',ADVANCE='NO') "_?_",T_
      ELSE
        WRITE(F, '(G15.8,A1)',ADVANCE='NO') S%WeitKg / S%V0, T_
      ENDIF
    ENDIF
    !
    DO J=1,N
      X=  tGrt(I,J)
      WRITE(F,'(G15.8,A1)',ADVANCE='NO') -X /Ln10, T_
    ENDDO
    !
    WRITE(F,*)
    !
  ENDDO DoSecond2
  !-----------------------------------/then, write the secondary species
  !
  WRITE(F,'(A,/)') "END SPECIES"
  !-------------------------------------------------------/SPECIES block
  CLOSE(F)
  !
  IF(iDebug>0) PRINT '(A)',"Table of LogK in "//TRIM(DirOut)//"_logk_elem.dtb"
  !--------------------------------------------------/logK from database
  !
  DEALLOCATE(tGrt)
  DEALLOCATE(vSolvDat)
  !
END SUBROUTINE Dtb_Tabulate_ForSystem

SUBROUTINE Dtb_Tabulate_EntropyZero(vS0Ele)
  USE M_T_Species,  ONLY: Species_EntropyZero
  USE M_Global_Vars,ONLY: vEle,vSpc

  REAL(dp):: vS0Ele(:)
  INTEGER :: I,N

  N= SIZE(vEle)

  DO I=1,SIZE(vSpc)
    vS0Ele(I)= Species_EntropyZero(vEle,vSpc(I))
    IF(iDebug>2) WRITE(71,'(A,A1,G15.6)') vSpc(I)%NamSp,T_,vS0Ele(I)
  ENDDO
  
  RETURN
ENDSUBROUTINE Dtb_Tabulate_EntropyZero

SUBROUTINE Dtb_GrtTable_Tabulate(sCode,vTPCond,vSolModlDat,vSpcDtb,vSpc,tGrt)
  USE M_Numeric_Const,ONLY: Ln10
  USE M_T_Species, ONLY: T_Species, T_SpeciesDtb
  USE M_Dtb_Const, ONLY: T_CK,R_jK,Tref
  USE M_IoTools,   ONLY: GetUnit, OutStrVec
  USE M_Files,     ONLY: DirOut,Files_Index_Write
  USE M_T_Tpcond,  ONLY: T_TPCond
  USE M_T_SolModel,ONLY: T_SolModelDat
  !
  CHARACTER(LEN=*),   INTENT(IN):: sCode
  TYPE(T_SolModelDat),INTENT(IN):: vSolModlDat(:)
  TYPE(T_TPCond),     INTENT(IN):: vTPCond(:)
  TYPE(T_SpeciesDtb), INTENT(IN):: vSpcDtb(:)
  TYPE(T_Species),    INTENT(IN):: vSpc(:)
  REAL(dp),           INTENT(IN):: tGrt(:,:)
  !
  TYPE(T_Species):: S
  INTEGER :: f,iSp,jTP,N,M
  REAL(dp):: X
  CHARACTER(LEN=15):: Cc
  REAL(dp),ALLOCATABLE:: vS0Ele(:)

  N= SIZE(vTPCond)
  M= SIZE(vSpc)

  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(DirOut)//"_"//TRIM(sCode)//".dtb")
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_"//TRIM(sCode)//".dtb",&
  & TRIM(sCode)//" table for current run")

  IF(TRIM(sCode)=="GIBBS2") THEN
    ALLOCATE(vS0Ele(1:M))
    CALL Dtb_Tabulate_EntropyZero(vS0Ele)
  ENDIF
  !
  !----------------------------------------------------- write header --
  SELECT CASE(TRIM(sCode))
  !
  CASE("LOGK","LOGK_EQ36")
    !
    WRITE(F,'(/,A)') "TP.TABLE"
    !
    CALL OutStrVec(F,vTPCond(:)%TdgC,Opt_S="TdgC")
    CALL OutStrVec(F,vTPCond(:)%Pbar,Opt_S="Pbar")
    !
    WRITE(F,'(A,/)') "END TP.TABLE"
    !
    WRITE(F,'(/,A)') "SOLVENT"
    !
    WRITE(F,'(A)') "NAME H2O"
    CALL OutStrVec(F,vSolModlDat(:)%Rho, Opt_S="Rho")
    CALL OutStrVec(F,vSolModlDat(:)%Eps, Opt_S="Eps")
    CALL OutStrVec(F,vSolModlDat(:)%dhA, Opt_S="DHA")
    CALL OutStrVec(F,vSolModlDat(:)%dhB, Opt_S="DHB")
    CALL OutStrVec(F,vSolModlDat(:)%BDOt,Opt_S="BDOT")
    !
    WRITE(F,'(A,/)') "END SOLVENT"
    !
    WRITE(F,'(/,A)') "SPECIES" !
    !
  !
  CASE("GIBBS","GIBBS2")
    !
    SELECT CASE(TRIM(sCode))
    CASE("GIBBS")   ;  WRITE(F,'(A)') "!Gibbs_FreeEnergy_alaBenson"
    CASE("GIBBS2")  ;  WRITE(F,'(A)') "!Gibbs_FreeEnergy_alaBerman"
    END SELECT

    !-- write array of Temperatures --
    WRITE(F,'(4(A,A1))',ADVANCE='NO') &
    & ".Type",T_,".Trace",T_,".Species",T_,".ECFORM",T_
    ! CALL OutStrVec(F,vTPCond(:)%TdgC,Opt_S="TdgC")
    CALL OutStrVec(F,vTPCond(:)%TdgC)
    
    !-- write array of Pressures --
    WRITE(F,'(4(A,A1))',ADVANCE='NO') &
    & ".Type",T_,".Trace",T_,".Species",T_,".ECFORM",T_
    ! CALL OutStrVec(F,vTPCond(:)%Pbar,Opt_S="Pbar")
    CALL OutStrVec(F,vTPCond(:)%Pbar)
  !
  END SELECT
  !----------------------------------------------------/ write header --
  !
  !---------------------------------------------------- write results --
  DO iSp=1,M
    !
    S=vSpc(iSp)
    !
    IF(S%iDtb>0) THEN
      Cc= TRIM(vSpcDtb(S%iDtb)%DtbTrace)
    ELSE
      Cc= "_"
    ENDIF
    !
    WRITE(F,'(A3,A1, A15,A1, A23,A1, A63,A1)',ADVANCE='NO') &
    & S%Typ,    T_, &
    & Cc,       T_, &
    & S%NamSp,   T_, &
    & S%Formula,T_
    !
    SELECT CASE(sCode)
    !
    CASE("LOGK")
      IF(S%Typ=="AQU") THEN
        !-- for aqu'species, write size parameter --
        WRITE(F, '(F15.1,A1)',ADVANCE='NO') S%AquSize, T_
      ELSE
        !-- for non aqueous species, write density --
        WRITE(F, '(G15.8,A1)',ADVANCE='NO') S%WeitKg /S%V0, T_
      ENDIF
      !-- write array of log10(-G/RT), i.e. "logK" --
      DO jTP=1,N
        WRITE(F,'(G15.8,A1)',ADVANCE='NO') -tGrt(iSp,jTP) /Ln10, T_
      ENDDO
    !
    CASE("GIBBS")
      !--- write Gibbs free energy --
      DO jTP= 1,N
        X= tGrt(iSp,jTP) *(vTPCond(jTP)%TdgC +T_CK) *R_jK
        ! X= X - Tref *vS0Ele(iSp) !-> to Berman-Brown
        WRITE(F,'(G15.8,A1)',ADVANCE='NO') X, T_
      ENDDO
    !
    CASE("GIBBS2")
      !--- write Gibbs free energy, in Berman-Brown convention --
      DO jTP= 1,N
        X= tGrt(iSp,jTP) *(vTPCond(jTP)%TdgC +T_CK) *R_jK
        X= X - Tref *vS0Ele(iSp) !-> to Berman-Brown
        WRITE(F,'(G15.8,A1)',ADVANCE='NO') X, T_
      ENDDO
    !
    END SELECT
    !
    WRITE(F,*)
    !
  ENDDO
  !---------------------------------------------------/ write results --
  !
  IF (TRIM(sCode)=="LOGK") WRITE(F,'(A,/)') "END SPECIES"
  !
  IF(TRIM(sCode)=="GIBBS2") DEALLOCATE(vS0Ele)
  !
  CLOSE(F)
  !
  IF(iDebug>0) THEN
    PRINT '(A)',"Table of LogK in "//TRIM(DirOut)//"_"//TRIM(sCode)//".dtb"
    IF(iDebug>2) CALL Pause_
  ENDIF
  !
ENDSUBROUTINE Dtb_GrtTable_Tabulate

SUBROUTINE Dtb_Test_EQ36
!--
!-- tabulate logK of formation of all species of vSpc
!-- along the EQ36 series of T,P points
!--
  USE M_Dtb_Const, ONLY: T_CK
  USE M_Dtb_Calc,  ONLY: DtbSpc_GrtTable_Build
  USE M_T_SolModel,ONLY: T_SolModel,T_SolModelDat !,SolventDat_Default
  USE M_SolModel_Tools,ONLY: SolModel_TP_Update
  !
  USE M_Global_Vars,ONLY: vSpcDtb,vSpc
  !
  TYPE(T_TPCond)      :: vTPcond(8)
  TYPE(T_SolModel)     :: Slv
  TYPE(T_SolModelDat)  :: vSolvDat(8)
  REAL(dp),ALLOCATABLE:: tGrt(:,:)
  INTEGER:: i
  !
  ALLOCATE(tGrt(1:SIZE(vSpc),1:8))
  tGrt= Zero
  !
  vTPcond(1)%TdgC= 0.01D0  ;  vTPcond(1)%Pbar= 1.013D0
  vTPcond(2)%TdgC= 25.0D0  ;  vTPcond(2)%Pbar= 1.013D0
  vTPcond(3)%TdgC= 60.0D0  ;  vTPcond(3)%Pbar= 1.013D0
  vTPcond(4)%TdgC= 100.D0  ;  vTPcond(4)%Pbar= 1.013D0
  vTPcond(5)%TdgC= 150.D0  ;  vTPcond(5)%Pbar= 4.757D0
  vTPcond(6)%TdgC= 200.D0  ;  vTPcond(6)%Pbar= 15.34D0
  vTPcond(7)%TdgC= 250.D0  ;  vTPcond(7)%Pbar= 39.74D0
  vTPcond(8)%TdgC= 300.D0  ;  vTPcond(8)%Pbar= 85.84D0
  !
  DO i=1,8
    !~ CALL SolventDat_Default(vSolvDat(i))
    !! provisional !! should CALL calculation of vSolvDat !!
    CALL SolModel_TP_Update(vTPcond(i)%TdgC +T_CK,vTPcond(i)%Pbar,Slv)
    vSolvDat(i)= Slv%Dat
  ENDDO
  !
  CALL DtbSpc_GrtTable_Build(vTPCond,vSpcDtb,vSpc,tGrt)
  !
  CALL Dtb_GrtTable_Tabulate("LOGK",vTPCond,vSolvDat,vSpcDtb,vSpc,tGrt)
  !
  DEALLOCATE(tGrt)
  !
ENDSUBROUTINE Dtb_Test_EQ36

SUBROUTINE Dtb_Tabulate_PSatH2O
!--
!-- tabulate (T,P) series a la EQ36, up to critical temp' of water
!--
  USE M_Dtb_Const,ONLY: T_CK
  USE M_IOTools,  ONLY: GetUnit
  USE M_Fluid_Calc
  !
  REAL(dp):: TdgC,TdgK,Pbar,Psat
  INTEGER :: f
  !
  CALL GetUnit(f)
  OPEN(f,FILE="psath2o.restab")
  !CALL Files_Index_Write(fHtm,&
  !& "psath2o.restab",&
  !& "T(deg)C-P(bar) along gas saturation")

  TdgC=Zero
  DO

    TdgK= TdgC +T_CK

    IF (TdgK<=647.25D0) CALL Eos_H2O_Psat(TdgK,Psat)

    IF(TdgC<=100._dp) THEN  ;  Pbar= 1.013D0
    ELSE                    ;  Pbar= Psat
    ENDIF

    WRITE(f,'(F7.2,A1,F7.2)') TdgC,t_,Pbar

    TdgC= TdgC +5.D0
    IF(TdgK>647.25D0) EXIT

  ENDDO

  CLOSE(f)
ENDSUBROUTINE Dtb_Tabulate_PSatH2O

SUBROUTINE Dtb_Test_Fluid
  USE M_Numeric_Const,ONLY: ln10
  USE M_Dtb_Const,  ONLY: T_CK, R_JK,Tref,Pref
  USE M_Files,      ONLY: DirOut
  USE M_IOTools,    ONLY: GetUnit
  USE M_TPcond_Read,ONLY: TPpath_Read,TPgrid_Build
  USE M_Path_Vars,  ONLY: vTPpath
  !
  USE M_Fluid_Calc
  !~ USE M_FluidMix_KerJac
  !
  LOGICAL :: Ok
  REAL(dp):: TdgC,TdgK,Pbar,Grt,H,S,V_m3
  INTEGER :: i,f1,f2,f3
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !~ REAL(dp):: XCO2
  !~ REAL(dp):: FugC_H2O,FugC_CO2,ActH2O,ActCO2

  !~ TdgC= 400.0D0
  !~ Pbar= 10.0D3
  !~ TdgK= TdgC + T_CK

  !~ XCO2= 0.01D0
  !~ DO

    !~ CALL EoS_KerJac_CalcMix( &
    !~ & TdgK,Pbar,XCO2, & !in
    !~ & FugC_H2O,FugC_CO2,ActH2O,ActCO2) !out

    !~ WRITE(11,'(5(G15.6,1X))') XCO2,FugC_H2O,FugC_CO2,ActH2O,ActCO2

    !~ XCO2= XCO2 + 0.01D0

    !~ IF(XCO2>0.995D0) EXIT

  !~ ENDDO

  !~ RETURN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF(iDebug>0) &
  & PRINT *,"Compute the thermodynamic properties of H2O"

  TdgK= Tref
  Pbar= Pref
  !
  CALL TPgrid_Build(Ok)
  IF(.NOT. Ok) CALL TPpath_Read(TdgK,Pbar)
  !
  CALL GetUnit(f1)
  OPEN(f1,FILE=TRIM(DirOut)//"_dtbfluid_ghio.restab")
  WRITE(f1,'(10(A,A1))') &
  & "Model",t_,"TdgC",t_,"log(Pbar)",t_,"Pbar",t_,"logK",t_, &
  & "G",    t_,"H",   t_,"S",        t_,"H-TS",t_,"Rho", t_
  !
  CALL GetUnit(f2)
  OPEN(f2,FILE=TRIM(DirOut)//"_dtbfluid_haar.restab")
  WRITE(f2,'(7(A,A1))') &
  & "Model",t_,"TdgC",t_,"log(Pbar)",t_,"Pbar",t_,"logK",t_,"G",t_,"Rho",t_
  !
  CALL GetUnit(f3)
  OPEN(f3,FILE=TRIM(DirOut)//"_dtbfluid_supcrt.restab")
  WRITE(f3,'(10(A,A1))') &
  & "Model",t_,"TdgC",t_,"log(Pbar)",t_,"Pbar",t_,"logK",t_, &
  & "G",    t_,"H",   t_,"S",        t_,"H-TS",t_,"Rho", t_
  !
  DO i=1,SIZE(vTPpath)

    TdgC= vTPpath(i)%TdgC
    Pbar= vTPpath(i)%Pbar
    TdgK= TdgC +T_CK

    CALL CalcFluid( &
    & "H2O","HAAR.GHIORSO",TdgK,Pbar, &
    & Grt,H,S,V_m3,Ok)
    IF(Ok) THEN
      WRITE(f1,'(A,A1,3(G15.6,A1))',ADVANCE='no') &
      & "HAARGHIO",t_,TdgC,t_,log10(Pbar),t_,Pbar,t_
      WRITE(f1,'(6(G15.6,A1))') &
      & -Grt/ln10,t_,  Grt*R_JK*TdgK, t_,  &
      &  H,       t_,  S,             t_,  &
      &  H-TdgK*S,t_,  0.01852D0/V_m3,t_
    ENDIF

    CALL CalcFluid( &
    & "H2O","HAAR",TdgK,Pbar, &
    & Grt,H,S,V_m3,Ok)
    IF(Ok) THEN
      WRITE(f2,'(A,A1,3(G15.6,A1))',ADVANCE='no') &
      & "HAAR",t_,TdgC,t_,log10(Pbar),t_,Pbar,t_
      WRITE(f2,'(3(G15.6,A1))') &
      &  -Grt/ln10,t_,  Grt*R_JK*TdgK,t_, 0.01852D0/V_m3,t_
    ENDIF

    CALL CalcFluid( &
    & "H2O","SUPCRT",TdgK,Pbar, &
    & Grt,H,S,V_m3,Ok)
    IF(Ok) THEN
      WRITE(f3,'(A,A1,3(G15.6,A1))',ADVANCE='no') &
      & "SUPCRT",t_,TdgC,t_,log10(Pbar),t_,Pbar,t_
      WRITE(f3,'(6(G15.6,A1))') &
      & -Grt/ln10,t_,  Grt*R_JK*TdgK,  t_,  &
      &  H,       t_,  S,              t_,  &
      &  H-TdgK*S,t_,  0.01852D0/V_m3, t_
    ENDIF

  ENDDO

  CLOSE(f1)
  CLOSE(f2)
  CLOSE(f3)

  IF(iDebug>0) THEN
    PRINT *,">> results in files dtbfluid_xxx.restab"
    IF(iDebug>2) CALL Pause_
  ENDIF

  DEALLOCATE(vTPpath)
  !
ENDSUBROUTINE Dtb_Test_Fluid

SUBROUTINE Dtb_Test_Species
!--
!-- tabulate the properties of all species
!--

  USE M_TPcond_Read,ONLY: TPpath_Read
  USE M_Dtb_Write
  !
  USE M_Dtb_Vars
  USE M_Path_Vars, ONLY: vTPpath
  !
  REAL(dp):: TdgK,Pbar
  
  TdgK= 298.15D0
  Pbar= 1.0D0
  !
  CALL TPpath_Read(TdgK,Pbar)
  !
  IF(SIZE(vDtbAquHkf)>0) CALL DtbAquHkf_Tabulate(vTPpath)
  ! IF(SIZE(vDtbAquThr)>0) CALL DtbAquThr_Tabulate(vTPpath)
  ! IF(SIZE(vDtbMinHkf)>0) CALL DtbMinHkf_Tabulate(vTPpath)
  IF(SIZE(vDtbMinThr)>0) CALL DtbMin_Tabulate("THR",vTPpath)
  IF(SIZE(vDtbMinHkf)>0) CALL DtbMin_Tabulate("HKF",vTPpath)
  
  IF(SIZE(vDtbAquHkf)>0) THEN
    IF(SIZE(vDtbMinThr)>0) CALL DtbSys_Tabulate("THR",vTPpath)
    IF(SIZE(vDtbMinHkf)>0) CALL DtbSys_Tabulate("HKF",vTPpath)
  END IF
  
  RETURN
END SUBROUTINE Dtb_Test_Species

SUBROUTINE Dtb_Test_H2OHkf
!--
!-- tabulate the properties of H2O used by the HKF model
!--
  USE M_Dtb_Const,  ONLY: T_CK
  USE M_IOTools,    ONLY: GetUnit
  USE M_TPcond_Read,ONLY: TPpath_Read,TPgrid_Build
  USE M_Path_Vars,  ONLY: vTPpath
  USE M_Fluid_Calc
  USE M_T_DtbH2OHkf
  !
  LOGICAL :: Ok
  REAL(dp):: TdgC,TdgK,Pbar,Psat
  INTEGER :: f1,f2
  INTEGER :: i
  TYPE(T_DtbH2OHkf):: pW

  TdgK= 298.15D0
  Pbar= 1.01325D0
  !
  CALL TPgrid_Build(Ok)
  IF(.NOT. Ok) CALL TPpath_Read(TdgK,Pbar)
  !
  CALL GetUnit(f1)
  OPEN(f1,FILE="dtbh2o.restab")
  WRITE(f1,'(11(A,A1))') &
  & "TdgC",t_,"log(Pbar)",t_,"Pbar",t_,&
  & "Rho",t_,"Eps",t_,"Q",t_,"X",t_,"Y",t_,"Z",t_,"dhA",t_,"dhB",t_
  !
  CALL GetUnit(f2)
  OPEN(f2,FILE="dtbh2o_alfa.restab")
  WRITE(f2,'(9(A,A1))') &
  & "TdgC",t_,"log(Pbar)",t_,"Pbar",t_,&
  & "Rho",t_,"Alfa",t_,"dAlfdT",t_,"Beta",t_,"dBetdT",t_,"dBetdP",t_
  !
  DO i=1,SIZE(vTPpath)

    TdgC= vTPpath(i)%TdgC
    TdgK= TdgC +T_CK
    Pbar= vTPpath(i)%Pbar

    Psat= Zero
    IF (TdgK<=647.25D0) CALL Eos_H2O_Psat(TdgK,Psat)

    IF(TdgK>647.25D0 .OR. Pbar>Psat) THEN
    !-> restrict to liquid or supercritic'domain
      
      !!CALL CalcRhoH2O( &
      !!& TdgK,Pbar, & !IN
      !!& RH2O,Alfa,Beta,dAlfdT,dBetdT,dBetdP) !OUT
      !!WRITE(f1,'(8(G15.6,A1))') &
      !!& TdgC,t_,log10(Pbar),t_, &
      !!& RH2O,t_,Alfa,t_,Beta,t_,dAlfdT,t_,dBetdT,t_,dBetdP,t_
      
      CALL DtbH2OHkf_Calc(TdgK,Pbar,pW)
      !
      WRITE(f1,'(11(G15.6,A1))') &
      & TdgC,t_,        log10(Pbar),t_, Pbar,t_, &
      & pW%Rho*1.D3,t_, pW%Eps,t_,      pW%Q,t_, &
      & pW%X,t_,        pW%Y,t_,        pW%Z,t_, &
      & pW%dhA,t_,      pW%dhB,t_

      WRITE(f2,'(9(G15.6,A1))') &
      & TdgC,t_,        log10(Pbar),t_, Pbar,t_, &
      & pW%Rho*1.D3,t_, &
      & pW%Alfa,t_,     pW%dAlfdT,t_, &             !thermal expansion.
      & pW%Beta,t_,     pW%dBetdT,t_, pW%dBetdP,t_  !compressibility
      
    ENDIF

  ENDDO
  CLOSE(f1)
  CLOSE(f2)
  !
  DEALLOCATE(vTPpath)
  
  RETURN
ENDSUBROUTINE Dtb_Test_H2OHkf

SUBROUTINE Dtb_Test_AquHkf
!--
!-- test computations on aqueous species (HKF model, HSV data)
!-- tabulate values of logK and partial molar volume
!-- if vSpc contains also minerals,
!-- their logK and density (kg/m3) are tabulated
!-- along points of a TP.GRID
!--
  USE M_Dtb_Const,    ONLY: T_CK,Tref,Pref
  USE M_Files,        ONLY: DirOut
  USE M_IOTools,      ONLY: GetUnit
  USE M_Numeric_Const,ONLY: Ln10
  USE M_T_Species,    ONLY: T_Species
  USE M_Dtb_Calc,     ONLY: Species_TP_Update_fromDtb
  USE M_T_DtbH2OHkf,  ONLY: T_DtbH2OHkf,DtbH2OHkf_Calc
  USE M_Fluid_Calc,   ONLY: Eos_H2O_Psat
  USE M_TPcond_Read,  ONLY: TPgrid_Build,TPpath_Read
  !
  USE M_Global_Vars,  ONLY: vSpcDtb,vSpc
  USE M_Dtb_Vars,     ONLY: vDtbAquHkf
  USE M_Path_Vars,    ONLY: vTPpath
  !
  REAL(dp):: TdgC,TdgK,Pbar,X
  !
  LOGICAL:: Ok
  INTEGER:: N, I, J
  INTEGER:: f1,f2
  TYPE(T_DtbH2OHkf):: PropsH2O
  !
  CALL TPgrid_Build(Ok)
  IF(.not. Ok) CALL TPpath_Read(Tref,Pref)
  !
  N= SIZE(vSpc)
  !
  IF(N<=0) RETURN
  !
  IF(SIZE(vDtbAquHkf)<1) THEN
    PRINT '(A)',"no HSV DATA on species"
    RETURN
  ENDIF
  !
  CALL GetUnit(f1)
  OPEN(f1,FILE=TRIM(DirOut)//"_tpgrid_logk.restab")
  WRITE(f1,'(3(A,A1))',ADVANCE='no') "TdgC",t_,"log(Pbar)",t_,"Pbar",t_
  !
  CALL GetUnit(f2)
  OPEN(f2,FILE=TRIM(DirOut)//"_tpgrid_dens.restab")
  WRITE(f2,'(3(A,A1))',ADVANCE='no') "TdgC",t_,"log(Pbar)",t_,"Pbar",t_
  !
  DO I=1,SIZE(vSpc) !N
    !!IF(vSpc(I)%Typ=="AQU") WRITE(f1,'(A,A1)',ADVANCE='no') TRIM(vSpc(I)%NamSp),t_
    WRITE(f1,'(A,A1)',ADVANCE='no') TRIM(vSpc(I)%NamSp),t_
  ENDDO
  WRITE(f1,*)
  DO I=1,SIZE(vSpc) !N
    !!IF(vSpc(I)%Typ=="AQU") WRITE(f1,'(A,A1)',ADVANCE='no') TRIM(vSpc(I)%NamSp),t_
    WRITE(f2,'(A,A1)',ADVANCE='no') TRIM(vSpc(I)%NamSp),t_
  ENDDO
  WRITE(f2,*)

  DO N=1,SIZE(vTPpath)

    TdgC= vTPpath(N)%TdgC
    TdgK= TdgC +T_CK
    Pbar= vTPpath(N)%Pbar

    !! PRINT '(A,F7.2,1X,F7.2)',"Dtb_Test_AquHkf ",TdgC,Pbar
    WRITE(f1,'(3(G15.6,A1))',ADVANCE='no') TdgC,t_,LOG10(Pbar),t_,Pbar,t_
    WRITE(f2,'(3(G15.6,A1))',ADVANCE='no') TdgC,t_,LOG10(Pbar),t_,Pbar,t_
    !
    CALL DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O)
    !-> solvent properties, for aqu'species
    !
    DO J=1,SIZE(vSpc)
      !
      CALL Species_TP_Update_fromDtb(TdgK,Pbar,PropsH2O,vSpcDtb,vSpc(J))
      !
      WRITE(f1,'(G15.6,A1)',ADVANCE='no') -vSpc(J)%G0rt/Ln10,t_
      !
      IF(vSpc(J)%Typ=="AQU") THEN
        WRITE(f2,'(G15.6,A1)',ADVANCE='no')  vSpc(J)%V0,t_
      ELSE
        X= vSpc(J)%WeitKg
        WRITE(f2,'(G15.6,A1)',ADVANCE='no')  X/vSpc(J)%V0,t_
      ENDIF
      !
    ENDDO

    WRITE(f1,*)
    WRITE(f2,*)

  ENDDO

  CLOSE(f1)
  CLOSE(f2)
  !
  DEALLOCATE(vTPpath)
  !
  RETURN
ENDSUBROUTINE Dtb_Test_AquHkf

ENDMODULE M_Dtb_Test

!~ SUBROUTINE Dtb_Test_AquHkf_
  !~ USE M_Dtb_Const,  ONLY: T_CK
  !~ USE M_IOTools,    ONLY: GetUnit
  !~ USE M_Files,      ONLY: DirOut
  !~ USE M_Numeric_Const,ONLY: Ln10
  !~ USE M_T_Species,  ONLY: T_Species
  !~ USE M_Dtb_Vars,   ONLY: vDtbAquHkf
  !~ USE M_Global_Vars,ONLY: vSpc
  !~ USE M_Dtb_Calc,   ONLY: DtbSpc_TP_Update
  !~ USE M_Dtb_Read,   ONLY: Dtb_Read_TPgrid
  !~ USE M_T_DtbH2OHkf,ONLY: T_DtbH2OHkf,DtbH2OHkf_Calc
  !~ !USE M_T_DtbAquHkf,ONLY: DtbAquHkf_Calc,DtbAquHkf_CalcThr
  !~ USE M_Fluid_Calc ,ONLY: Eos_H2O_Psat
  !~ !
  !~ !REAL(dp),INTENT(IN):: T_Min, T_Max, T_delta
  !~ !REAL(dp),INTENT(IN):: P_Min, P_Max, P_delta, P_ratio
  !~ !
  !~ LOGICAL :: Ok
  !~ REAL(dp):: T_Min, T_Max, T_delta, T_ratio
  !~ REAL(dp):: P_Min, P_Max, P_delta, P_ratio
  !~ REAL(dp):: TdgC,TdgK,Pbar,Psat,X
  !~ !
  !~ INTEGER:: N, I, J
  !~ INTEGER:: f1,f2

  !~ TYPE(T_DtbH2OHkf):: PropsH2O
  !~ !
  !~ T_Min=   0._dp     !5._dp     !
  !~ T_Max=   200._dp   !500._dp   !
  !~ T_delta= 5._dp     !5._dp     !
  !~ T_ratio= One
  !~ P_Min=   1._dp     !1._dp     !
  !~ P_Max=   1000._dp  !1000._dp  !
  !~ P_ratio= 1.2_dp    !1.1_dp    !
  !~ P_delta= 5._dp
  !~ !
  !~ CALL Dtb_Read_TPgrid( &
  !~ & Ok,T_Min,T_Max,T_ratio,T_delta,P_Min,P_Max,P_ratio,P_delta)
  !~ !
  !~ IF(iDebug>3) THEN
    !~ PRINT '(A,4G15.6)',"T_Min,T_Max,T_ratio,T_delta",T_Min,T_Max,T_ratio,T_delta
    !~ PRINT '(A,4G15.6)',"P_Min,P_Max,P_ratio,P_delta",P_Min,P_Max,P_ratio,P_delta
  !~ CALL Pause_
  !~ ENDIF
  !~ !
  !~ N= SIZE(vSpc)
  !~ IF(N<=0) RETURN
  !~ IF(SIZE(vDtbAquHkf)<1) THEN
    !~ PRINT '(A)',"no HSV DATA on species"
    !~ RETURN
  !~ ENDIF
  !~ !
  !~ CALL GetUnit(f1) ; OPEN(f1,FILE=TRIM(DirOut)//"_tpgrid_logk.restab")
  !~ CALL GetUnit(f2) ; OPEN(f2,FILE=TRIM(DirOut)//"_tpgrid_dens.restab")
  !~ WRITE(f1,'(3(A,A1))',ADVANCE='no') "TdgC",t_,"log(Pbar)",t_,"Pbar",t_
  !~ WRITE(f2,'(3(A,A1))',ADVANCE='no') "TdgC",t_,"log(Pbar)",t_,"Pbar",t_
  !~ DO I=1,SIZE(vSpc) !N
    !~ !!IF(vSpc(I)%Typ=="AQU") WRITE(f1,'(A,A1)',ADVANCE='no') TRIM(vSpc(I)%NamSp),t_
    !~ WRITE(f1,'(A,A1)',ADVANCE='no') TRIM(vSpc(I)%NamSp),t_
    !~ WRITE(f2,'(A,A1)',ADVANCE='no') TRIM(vSpc(I)%NamSp),t_
  !~ ENDDO
  !~ WRITE(f1,*)
  !~ WRITE(f2,*)
  !~ !
  !~ TdgC= T_Min
  !~ Do_T: DO WHILE(TdgC<T_Max)
    !~ Pbar= P_Min
    !~ Do_P: DO WHILE(Pbar<=P_Max)
      !~ TdgK= TdgC +T_CK
      !~ !PRINT '(F7.2,1X,F7.2)',TdgC,Pbar
      !~ !WRITE(f1,'(2(G15.6,A1))',ADVANCE='no') TdgC,t_,Pbar,t_
      !~ !
      !~ Psat= Zero
      !~ IF (TdgK<=647.25D0) CALL Eos_H2O_Psat(TdgK,Psat)
      !~ IF(TdgK>647.25D0 .OR. Pbar>Psat) THEN
        !~ WRITE(f1,'(3(G15.6,A1))',ADVANCE='no') TdgC,t_,log10(Pbar),t_,Pbar,t_
        !~ WRITE(f2,'(3(G15.6,A1))',ADVANCE='no') TdgC,t_,log10(Pbar),t_,Pbar,t_
        !~ !
        !~ CALL DtbH2OHkf_Calc(TdgK,Pbar,PropsH2O) !-> solvent properties, for aqu'species
        !~ !
        !~ DO J=1,SIZE(vSpc)
          !~ CALL DtbSpc_TP_Update(TdgK,Pbar,PropsH2O,vSpc(J))
          !~ WRITE(f1,'(G15.6,A1)',ADVANCE='no') -vSpc(J)%G0rt/Ln10,t_
          !~ IF(vSpc(J)%Typ=="GAS") THEN
            !~ X= EXP(vSpc(J)%LnPhi -LOG(Pbar))
            !~ WRITE(f1,'(G15.6,A1)',ADVANCE='no') X,t_
          !~ ENDIF
          !~ IF(vSpc(J)%Typ=="AQU") THEN
            !~ WRITE(f2,'(G15.6,A1)',ADVANCE='no')  vSpc(J)%V0,t_
          !~ ELSE
            !~ X= vSpc(J)%WeitKg
            !~ WRITE(f2,'(G15.6,A1)',ADVANCE='no')  X/vSpc(J)%V0,t_
          !~ ENDIF
        !~ ENDDO
        !~ WRITE(f1,*)
        !~ WRITE(f2,*)
      !~ ENDIF
      !~ !WRITE(f1,'(2(G15.6,A1))') Grt,t_,H,t_ !,S,t_,18.052D0/V_m3,t_
      !~ !
      !~ IF(P_ratio/=One) THEN; Pbar= Pbar *P_ratio
      !~ ELSE                 ; Pbar= Pbar + P_delta
      !~ ENDIF
      !~ !
    !~ ENDDO Do_P
    !~ TdgC= TdgC + T_delta
  !~ ENDDO Do_T
  !~ CLOSE(f1)
  !~ CLOSE(f2)
  !~ !
  !~ RETURN
!~ ENDSUBROUTINE Dtb_Test_AquHkf_

