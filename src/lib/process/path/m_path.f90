MODULE M_Path
!--
!-- implements various path calculations
!--
  USE M_Kinds
  USE M_Trace,ONLY: Stop_,fTrc,T_,iDebug,Pause_,Warning_
  USE M_T_Component,ONLY: T_Component

  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Path_Execute
  !
  INTEGER,PARAMETER:: SafeStep= 255
  !
  TYPE(T_Component),ALLOCATABLE:: vCpnMix(:)
  REAL(dp),         ALLOCATABLE:: vSys1(:), vSys2(:)
  !
  REAL(dp):: TdgKMix, PbarMix
  REAL(dp):: pH_,pO2
  !
CONTAINS

SUBROUTINE Path_Execute(EquMod5)
  !---------------------------------------------------------------------
  USE M_Files,       ONLY: NamFInn
  USE M_System_Tools,ONLY: System_TP_Update
  USE M_TPcond_Read, ONLY: TPpath_Read, TPgrid_Build
  USE M_Path_Read,   ONLY: Path_ReadMode, Path_ReadParam
  USE M_Basis,       ONLY: Basis_Change
  USE M_Equil,       ONLY: Equil_Calc
  USE M_Equil_Read,  ONLY: Equil_Read_PhaseAdd
  !
  USE M_Global_Vars, ONLY: vFas
  USE M_System_Vars, ONLY: TdgK,Pbar,vCpn
  USE M_Path_Vars,   ONLY: Path_Vars_Clean,vTPpath,vFasFound
  !---------------------------------------------------------------------
  CHARACTER(LEN=5), INTENT(IN):: EquMod5
  !
  CHARACTER(LEN=3) :: PathMod3
  CHARACTER(LEN=3) :: EquMod3
  CHARACTER(LEN=80):: Msg
  LOGICAL:: Ok,OkGrid,OkMix,Singular
  LOGICAL:: TPpath
  INTEGER:: iErr,I
  !---------------------------------------------------------------------
  
  ALLOCATE(vFasFound(1:SIZE(vFas)))
  vFasFound= .FALSE.
  !
  TPpath= (EquMod5(4:5)=="TP")
  EquMod3= EquMod5(1:3)
  !
  !------------------------------------------------------- initialize --
  IF(TPpath) THEN
    !
    CALL TPgrid_Build(OkGrid)
    IF(.NOT. OkGrid) CALL TPpath_Read(TdgK,Pbar)
    !
    PathMod3= "TP_"
    !
    IF(EquMod3(1:2)=="EQ") THEN
      IF(COUNT(vCpn(:)%Statut=="MOBILE")>0 .OR. &
      &  COUNT(vCpn(:)%Statut=="BUFFER")>0) THEN
        CALL Equil_Calc("SPC")
        CALL Basis_Change("DYN",vCpn)
      ENDIF
    ENDIF
    !
  ELSE
    !
    CALL Path_ReadMode( & 
    & NamFInn, &
    & PathMod3,Ok,Msg)
    !
    IF(.NOT.Ok) THEN
      IF(iDebug>0) CALL Warning_(TRIM(Msg))
      RETURN !===================================< i/o error, return ==
    ENDIF
    !
    !--------------------------- compute speciation of current system --
    CALL Equil_Calc("SPC")
    ! CALL Equil_Calc(EquMod3)
    !
    SELECT CASE(PathMod3)
    !
    CASE("ADD")
      !
      !--- change basis, MOBILE > INERT
      CALL Basis_Change("EQU",vCpn)
      CALL System_TP_Update(TdgK,Pbar)
      !
    CASE("ADA")
      !
      !--- change basis, MOBILE > INERT
      CALL Basis_Change("EQU",vCpn)
      CALL System_TP_Update(TdgK,Pbar)
      !
      IF(EquMod3(1:2)=="EQ") CALL Equil_Read_PhaseAdd(vFas)
      !
    CASE("MIX")
      !
      !--- save composition of system 1 to vSys1
      CALL System_TP_Update(TdgK,Pbar)
      !
      CALL Basis_Change("EQU",vCpn)
      !
      ALLOCATE(vSys1(1:SIZE(vCpn)))
      vSys1(:)= vCpn(:)%Mole
      !
      !--- compute composition of system 2 and save to vSys2
      ALLOCATE(vCpnMix(1:SIZE(vCpn)))
      !
      CALL Path_Calc_FluidMix(TdgKMix,PbarMix,OkMix)
      !
      ALLOCATE(vSys2(1:SIZE(vCpn)))
      vSys2(:)= vCpnMix(:)%Mole
      !
    ENDSELECT
    !---------------------------/compute speciation of current system --
    !
    IF(EquMod3(1:2)=="EQ") CALL Equil_Read_PhaseAdd(vFas)
    !
    CALL Path_ReadParam( &
    & NamFInn,   &  !IN
    & PathMod3,  &  !IN
    & vCpn,      &  !IN
    & TdgK,Pbar, &  !IN
    & Ok,Msg)       !OUT
    !
    IF(PathMod3=="PHP") THEN
      pH_= 1.0D0
      pO2= 50.0D0
    ENDIF
    !
    IF(.NOT.Ok) THEN
      IF(iDebug>0) CALL Warning_(TRIM(Msg))
      RETURN !====================================< i/o error, return ==
    ENDIF
    !
  ENDIF
  !------------------------------------------------------/ initialize --
  !
  !----------------------------------------------------- compute path --
  CALL Path_Single(PathMod3,EquMod3)
  !----------------------------------------------------/ compute path --
  !
  !------------------------------------------------------------ clean --
  !
  IF(PathMod3=="MIX") DEALLOCATE(vSys1,vSys2,vCpnMix)
  !
  CALL Path_Vars_Clean
  !-----------------------------------------------------------/ clean --
  !
  RETURN
ENDSUBROUTINE Path_Execute

SUBROUTINE Path_Double
  
END SUBROUTINE Path_Double

SUBROUTINE Path_Single(PathMod3,EquMod3)
  !---------------------------------------------------------------------
  USE M_Dtb_Const,   ONLY: T_CK
  USE M_System_Tools,ONLY: System_TP_Update
  USE M_Basis
  USE M_Equil_Tools, ONLY: Equil_Zero,Equil_Restart,Equil_Save
  USE M_Equil_Tools, ONLY: Equil_Trace_Init,Equil_Trace_Close
  USE M_Equil_Tools, ONLY: Equil_Errors,Equil_Errors_Warning
  USE M_Equil_Specia,ONLY: Equil_Specia
  USE M_Equil_1,     ONLY: Equil_Eq1
  USE M_Equil_2,     ONLY: Equil_Eq2
  USE M_Equil_Write
  USE M_Path_Write
  !
  USE M_Global_Vars, ONLY: vFas
  USE M_System_Vars, ONLY: TdgK,Pbar,vCpn
  USE M_Equil_Vars,  ONLY: Equil_Vars_Clean,vYesList
  USE M_Path_Vars,   ONLY: vTPpath,vFasFound,DimPath,TotalMixStep,vLPath
  !---------------------------------------------------------------------
  CHARACTER(LEN=3),INTENT(IN):: PathMod3
  CHARACTER(LEN=3),INTENT(IN):: EquMod3
  !---------------------------------------------------------------------
  REAL(dp):: CpuBegin,CpuEnd
  INTEGER :: I,nCp,nStep,iErr
  LOGICAL :: PathExit,PathRedox
  !---------------------------------------------------------------------
  
  CALL System_TP_Update(TdgK,Pbar)
  !
  CALL Basis_Change("SPC",vCpn)
  !
  CALL Equil_Zero(EquMod3)
  CALL Equil_Trace_Init
  !
  !-------------------- initial speciation in case PATH is ADD or MIX --
  IF(PathMod3=="ADD" .OR. &
  &  PathMod3=="ADA" .OR. &
  &  PathMod3=="MIX") THEN
    !
    CALL Equil_Restart
    !
    SELECT CASE(EquMod3)
    
    CASE("SPC")
      CALL Equil_Specia(iErr)
    
    CASE("EQ1")
      CALL Equil_Eq1(.TRUE.,iErr)
    
    CASE("EQ2")
      CALL Equil_Eq2(.TRUE.,iErr)
      IF(iErr<0) CALL Equil_Eq1(.TRUE.,iErr)
      
    CASE("EQM")
      CALL Equil_Eq2(.FALSE.,iErr)
      
    END SELECT
    !
    IF(iErr<0) CALL Equil_Errors(iErr)
    !
    CALL Equil_Save
    !
  ENDIF
  !-------------------/ initial speciation in case PATH is ADD or MIX --
  !
  IF (TRIM(PathMod3)=="MIX") THEN
    CALL Basis_Change("EQU",vCpn)
    vSys1(:)= vCpn(:)%Mole
    !
    DimPath= TotalMixStep
    !
    ALLOCATE(vTPpath(DimPath))
    DO I=1,TotalMixStep
      vTPpath(I)%TdgC= TdgK + (TdgKMix - TdgK)*I/TotalMixStep -T_CK
      vTPpath(I)%Pbar= Pbar
    ENDDO
  ENDIF
  !
  CALL Path_Write_FasEnTete
  !
  IF(iDebug>2) CALL CPU_TIME(CpuBegin)
  !
  !-------------------------------------------------------- path loop --
  nStep=1
  DO
    
    IF(nStep>SafeStep) THEN
      IF(iDebug>1) PRINT '(A)',"Reached max' number of steps !!"
      EXIT !to prevent unterminated loops !!!
    ENDIF
    !
    !IF(nStep>SIZE(vTPpath)) EXIT !------------------- EXIT path loop --
    !
    IF(PathMod3 /= "LGK") THEN
      IF( ABS(vTPpath(nStep)%TdgC +T_CK -TdgK) > 1.D-2 .OR. &
      &   ABS(vTPpath(nStep)%Pbar       -Pbar) > 1.D-2) THEN
        !
        TdgK= vTPpath(nStep)%TdgC +T_CK
        Pbar= vTPpath(nStep)%Pbar
        IF(iDebug>1) PRINT &
        & '(A,I3,2G15.6)',"TP changed, at Step ",nStep,TdgK-T_CK,Pbar
        !
        CALL System_TP_Update(TdgK,Pbar)
        !
      ENDIF
    ENDIF
    !
    !----------------------------------------- conditions of new step --
    SELECT CASE(PathMod3)
    !
    CASE("TP_")
    !--- path along the TP table --
      !
      IF(nStep>SIZE(vTPpath)) EXIT !=================< EXIT path loop ==
      !
      !! !--- conditions of new step --
      !! TdgK= vTPpath(nStep)%TdgC +T_CK
      !! Pbar= vTPpath(nStep)%Pbar
      !! CALL System_TP_Update(TdgK,Pbar)
      !! !---/
      !
    CASE("PHP")
      CALL Path_NewStep_pHpE(PathExit)
      IF(PathExit) pH_= pH_ +1.0D0
      IF(pH_>13.0D0) EXIT
      !
    CASE("ADD","ADA","MIX","CHG","EVP","LGK")
    !--- chemical path defined by the PATH block --
      !
      !--- conditions of new step --
      CALL Path_NewStep(PathMod3,nStep,PathExit)
      !---/
      !
    END SELECT
    !----------------------------------------/ conditions of new step --
    !
    !------------------------------------------------ compute new step--
    CALL Equil_Restart
    !
    SELECT CASE(EquMod3)
    !
    CASE("SPC")
      CALL Equil_Specia(iErr)
    CASE("EQM")
      CALL Equil_Eq2(.FALSE.,iErr)
      CALL FasFound_Record(vFas,vFasFound)
    CASE("EQ1")
      CALL Equil_Eq1(.TRUE.,iErr)
      CALL FasFound_Record(vFas,vFasFound)
    CASE("EQ2")
      CALL Equil_Eq2(.TRUE.,iErr)
      IF(iErr<0) CALL Equil_Eq1(.TRUE.,iErr)
      CALL FasFound_Record(vFas,vFasFound)
    !~ CASE("EQ3")
      !~ CALL Equil_Eq3(iErr)
      !~ CALL FasFound_Record(vFas,vFasFound)
    !
    END SELECT
    !
    IF(iErr<0) THEN
      CALL Equil_Errors_Warning(iErr)
      iErr= 0
      IF(PathExit) EXIT !============================< EXIT path loop ==
      nStep= nStep+1
      CYCLE
    END IF
    !
    !~ IF(iDebug==4) CALL Basis_Change_Wrk(vCpn)
    !
    CALL Equil_Save
    !-----------------------------------------------/ compute new step--
    !
    SELECT CASE(EquMod3)
    CASE("SPC"); CALL Path_Write_Line(WrCount=nStep,WrCod="SPC")
    CASE("EQM"); CALL Path_Write_Line(WrCount=nStep,WrCod="EQM",vYes=vYesList)  
    CASE("EQ1"); CALL Path_Write_Line(WrCount=nStep,WrCod="EQ1",vYes=vYesList)  
    CASE("EQ2"); CALL Path_Write_Line(WrCount=nStep,WrCod="EQ2",vYes=vYesList)   
    END SELECT
    !
    IF(iDebug==4) CALL Path_Write_Distrib(nStep)
    !
    CALL Path_Write_FasAff(nStep)
    !IF(OkKinetics) CALL KinRateActiv_Test(nStep)
    !
    IF(PathExit) EXIT !==============================< EXIT path loop ==
    !
    nStep= nStep+1
    !
  ENDDO
  !-------------------------------------------------------/ path loop --
  IF(iDebug>2) CALL CPU_TIME(CpuEnd)
  IF(iDebug>2) PRINT '(A,G15.6)',"CPU=", CpuEnd -CpuBegin
  !
  IF(iDebug>0 .AND. EquMod3(1:2)=="EQ") CALL FasFound_Sho(vFas,vFasFound)
  !
  CALL Equil_Trace_Close
  CALL Equil_Vars_Clean
  CALL Path_Files_Close
  !
  RETURN
END SUBROUTINE Path_Single

SUBROUTINE Path_NewStep_pHpE(PathExit)
!
! 2.H2O=  O2 + 4.H+ + 4.e-
! 2.G(H2O)+lnA_H2O=  G(O2)+lnA_O2  +4.lnA_H+ +4.lnA_e-
! 4.(pH + pE)=  [G(O2) +lnA_O2 -2.G(H2O) -2.lnA_H2O]/ln10
! log10(Act(O2))=4*(pH_+pE_)+ 2[G(H2O) + lnA_H2O] - [G(O2) + lnA_O2]
!
! IF redox species is H2aq:
! 2H+ + 2E-=  H2 -> pH + pE=  -( vSpc(iOx)%G0 + vLnAct(iOx)) /2
!
!  pE= & 
!  & (  vSpc(isO2)%G0rt + vSpc(isO2)%Dat%LAct          &
!  &  -(vSpc(isW )%G0rt + vSpc(isW )%Dat%LAct) *Two  ) &
!  & /LN10 /4.0D0 &
!  & - pH
!  
!  pE= vSpc(isO2)%G0rt     -vSpc(isW )%G0rt     *Two &
!  & + vSpc(isO2)%Dat%LAct -vSpc(isW )%Dat%LAct *Two
!  pE= pE /LN10 /4.0D0 - pH
!
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Global_Vars,  ONLY: vSpc
  USE M_Basis_Vars,   ONLY: isH_,isO2
  !
  LOGICAL,INTENT(OUT):: PathExit
  !
  REAL(dp):: pHmin,pHmax !,pEmin,pEmax
  REAL(dp):: pH,pK_O2 !,pE
  !
  vSpc(isH_)%Dat%LAct= -pH_ *LN10
  vSpc(isO2)%Dat%LAct= -pO2 *LN10
  !
  pO2= pO2 + 1.0D0
  PathExit= (pO2>1.0D0)
  !
END SUBROUTINE Path_NewStep_pHpE

SUBROUTINE FasFound_Sho(vFas,vFasFound)
  USE M_T_Phase
  TYPE(T_Phase),INTENT(IN):: vFas(:)
  LOGICAL,      INTENT(IN):: vFasFound(:)
  !
  INTEGER:: I
  !
  IF(COUNT(vFasFound)>0) THEN
    IF(iDebug>0) THEN
      PRINT '(/,A,/)',"phases encountered in PATHEQU:"
      DO I=1,SIZE(vFas)
        IF(vFasFound(I)) PRINT '(A)',vFas(I)%NamFs
      ENDDO
    ENDIF
  ENDIF
  !
END SUBROUTINE FasFound_Sho

SUBROUTINE FasFound_Record(vFas,vFasFound)
!--
!-- keep record of phases found during whole path --
!--
  USE M_T_Phase
  TYPE(T_Phase),INTENT(IN):: vFas(:)
  LOGICAL,   INTENT(INOUT):: vFasFound(:)
  !
  INTEGER:: iFs
  !
  DO iFs=1,SIZE(vFas) 
    IF(vFas(iFs)%Mole>Zero) vFasFound(iFs)=.TRUE.
  ENDDO
  !
  RETURN
END SUBROUTINE FasFound_Record

SUBROUTINE Path_NewStep(PathMod3,nStep,PathExit)
  USE M_Files
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Dtb_Const,    ONLY: T_CK
  USE M_T_Species,    ONLY: Species_Index
  USE M_T_MixPhase,   ONLY: T_MixPhase, MixPhase_CalcActivs
  !
  USE M_Global_Vars,ONLY: vEle,vSpc,vFas,vMixFas,vMixModel,nAq
  USE M_System_Vars,ONLY: vCpn,TdgK,Pbar
  USE M_Basis_Vars, ONLY: isW,iH_,MWSv,tAlfFs,iBal
  USE M_Path_Vars,  ONLY: &
  & vPhasBegin,vPhasFinal,vFasFound, &
  & vLPath,tPathData,DimPath,TotalMixStep, &
  & iLogK,vPathLogK
  !
  CHARACTER(LEN=3),INTENT(IN) :: PathMod3
  INTEGER,         INTENT(IN) :: nStep
  LOGICAL,         INTENT(OUT):: PathExit
  !
  INTEGER,PARAMETER:: AdjustMode= 1
  !
  INTEGER :: iCp,J,iFs,iP,NPole,nCp
  REAL(dp):: X
  !
  TYPE(T_MixPhase):: Fas, PhasBegin, PhasFinal
  
  PathExit=.FALSE.
  !
  nCp= SIZE(vCpn)
  !
  IF (iDebug>0) WRITE(fTrc,'(/,A,I3)') "< Path_NewStep ==== step= ", nStep
  
  SELECT CASE(PathMod3)
  
  !------------------------------------------------------- CASE "LGK" --
  CASE("LGK")
    IF(iDebug>1) PRINT '(A,I3,G15.6)',"nStep",nStep,vPathLogK(nStep)
    !
    vSpc(iLogK)%G0rt= -vPathLogK(nStep) *Ln10
    DO J=1,SIZE(vFas)
      IF(vFas(J)%iSpc == iLogK) vFas(J)%Grt= vSpc(iLogK)%G0rt
    ENDDO
    PathExit= nStep==SIZE(vPathLogK)
    !print *,"nStep          = ",nStep
    !print *,"SIZE(vPathLogK)= ",SIZE(vPathLogK)
    !pause
  !------------------------------------------------------/ CASE "LGK" --
  
  !------------------------------------------ CASE "EVP" (-EVAPORATE) --
  CASE("EVP")
    X= vCpn(isW)%Mole *0.10D0
    vCpn(isW)%Mole= vCpn(isW)%Mole - X
    vCpn(iH_)%Mole= vCpn(iH_)%Mole - X*2.0D0
    PathExit= nStep==255
    
  !--------------------------------------------- CASE "CHG" (-CHANGE) --
  CASE("CHG")
  ! CASE("CHG","EVP")
    !
    IF(iDebug>0) PRINT '(A,I3)',"nStep",nStep
    !
    DoComponent: DO iCp=1,nCp
      !
      IF(vLPath(iCp)) THEN
      !-------------------------------------- change in component iCp --
        J= vCpn(iCp)%iSpc !-> related species
        !
        SELECT CASE(TRIM(vCpn(iCp)%Statut))
        
        CASE("MOBILE","BUFFER")
        !--------------- case of CHANGE on mobile or buffer component --
        !----------------- -> change ACTIVITY -> change vSpc(J)%LnAct --
          
          IF(vCpn(iCp)%iMix==0) THEN
          !--------------- mobile species is aqueous or in pure phase --
            vCpn(iCp)%LnAct= - tPathData(iCp,nStep) *Ln10
            vSpc(J)%Dat%LAct= vCpn(iCp)%LnAct
            !
            PathExit= (nStep==DimPath)
            !
          ELSE
          !-------------- mobile species in non-aqueous mixture phase --
            PhasBegin= vMixFas(vPhasBegin(iCp))
            PhasFinal= vMixFas(vPhasFinal(iCp))
            Fas= PhasBegin
            !
            NPole= vMixModel(Fas%iModel)%NPole
            DO iP=1,NPole
              IF(Fas%vLPole(iP)) THEN
                Fas%vXPole(iP)= &
                & (PhasBegin%vXPole(iP)*(100 - nStep)  &
                +  PhasFinal%vXPole(iP)*(nStep)      ) &
                & /100.0D0
              ENDIF
            ENDDO
            !
            CALL MixPhase_CalcActivs( &
            & TdgK,Pbar,&
            & vMixModel(Fas%iModel),&
            & Fas)
            !
            !-------------------------------------------------- debug --
            IF(iDebug>2) THEN 
              WRITE(fTrc,'(A15,A1)',ADVANCE='NO') "PoleFraction=",T_
              DO iP=1,NPole
                IF(Fas%vLPole(iP)) &
                & WRITE(fTrc,'(G15.6,A1)',ADVANCE='NO') Fas%vXPole(iP),T_
              ENDDO
              WRITE(fTrc,*)
              WRITE(fTrc,'(A15,A1)',ADVANCE='NO') "Activity=",T_
              DO iP=1,NPole
                IF(Fas%vLPole(iP)) &
                & WRITE(fTrc,'(G15.6,A1)',ADVANCE='NO') EXP(Fas%vLnAct(iP)),T_
              ENDDO
              WRITE(fTrc,*)
            ENDIF
            !--------------------------------------------------/debug --
            !
            vSpc(J)%Dat%LAct= Fas%vLnAct(vCpn(iCp)%iPol) !not necessary ???
            vCpn(iCp)%LnAct=  vSpc(J)%Dat%LAct
            !
            IF(iDebug>0) WRITE(fTrc,'(A4,A3,A1, A4,A23,A1, A7,G15.6)') &
            & "Cpn=",vEle(vCpn(iCp)%iEle)%NamEl,T_, &
            & "Spc=",vSpc(J)%NamSp,T_, &
            & "LogAct=", vSpc(J)%Dat%LAct/Ln10
            !
            PathExit=(nStep==99)
            !
          ENDIF
          !------------/ mobile species in non-aqueous sol'phase -------
        !--------------/ CASE of CHANGE on mobile or buffer component --

        CASE("INERT") 
        !-------------------------- CASE of CHANGE on INERT component --
        !---- for inert species, change MOLE -> change vCpn(iCp)%Mole --
          vCpn(iCp)%Mole= tPathData(iCp,nStep)
          PathExit= (nStep==DimPath)
          !
          !------------------------------ bri-collage !!!§§§->rework !!!
          !-- vCpn will be used for the FLUID system !!!
          !-- it coincides with total system in SPC mode,
          !-- but, in EQU mode, we control here the total system ...
          !-- substract mole nrs from non-fluid phases
          !-- (case of EQn calculation)
          !-------------------------------------------------------------
          DO iFs=1,SIZE(vFas)
            IF(vFas(iFs)%Mole>Zero) &
            & vCpn(iCp)%Mole= vCpn(iCp)%Mole - tAlfFs(iCp,iFs) *vFas(iFs)%Mole
          ENDDO
          !
          DO iFs=1,SIZE(vFas)
            IF(vFas(iFs)%Mole>Zero) vFasFound(iFs)=.TRUE.
          ENDDO
          !
        !-------------------------/ CASE of CHANGE on INERT component --
        END SELECT
      ENDIF ! IF(vLPath(iCp))
    ENDDO DoComponent ! 
    !
    IF(PathMod3=="EVP") THEN
      DO iFs=1,SIZE(vFas)
        IF(vFas(iFs)%Mole>Zero) vFas(iFs)%Mole= 1.D-12
      ENDDO
    ENDIF
    !
    CALL NewStep_Adjust_TotO_TotH( & !
    & AdjustMode,        &           !
    & iBal,isW,iH_,vEle, &           !IN
    & vCpn)                          !INOUT
    !
  !-------------------------------------------------------/CASE "CHG" --
  
  !------------------------------------------------------- CASE "ADD" --
  CASE("ADD","ADA")
    !
    IF(iDebug>0) PRINT '(A,I3)',"nStep",nStep
    !
    vCpn(:)%Mole= tPathData(:,nStep)
    !
    !------------------------------------ bri-collage !!!§§§->rework !!!
    !-- vCpn is the FLUID system !!!
    !-- it coincides with total system in SPC mode,
    !-- but, in EQU mode, we control here the total system ...
    !----- substract mole nrs from non-fluid phases (case of EQn calculation) --
    DO iCp=1,SIZE(vCpn)
      DO iFs=1,SIZE(vFas)
        IF(vFas(iFs)%Mole>Zero) &
        & vCpn(iCp)%Mole= vCpn(iCp)%Mole - tAlfFs(iCp,iFs) *vFas(iFs)%Mole
      ENDDO
    ENDDO
    !-----/
    !
    IF(iDebug>2) THEN
      WRITE(61,'(A,1X)',ADVANCE="no") "Path_NewStep"
      DO J=1,SIZE(vCpn)
        WRITE(61,'(G15.6,1X)',ADVANCE="no") vCpn(J)%Mole
      ENDDO
      WRITE(61,*)
    ENDIF
    !
    PathExit= nStep==DimPath
  !-------------------------------------------------------/CASE "ADD" --
  
  !------------------------------------------------------- CASE "MIX" --
  CASE("MIX")
    !
    IF(nStep==TotalMixStep) THEN
      PathExit=.TRUE.
    ELSE
      DO iCp=1,SIZE(vCpn)
        !
        vCpn(iCp)%Mole= vSys1(iCp) &
        &             + (vSys2(iCp) -vSys1(iCp)) *nStep/REAL(TotalMixStep)
        !
        IF(iDebug>0) WRITE(fTrc,'(A3,1X,3G15.6,1X,A)') &
        & vEle(vCpn(iCp)%iEle)%NamEl, &
        & vSys1(iCp),vSys2(iCp),vCpn(iCp)%Mole,&
        & vCpn(iCp)%Statut
        !
      ENDDO
      !
    ENDIF
  !-------------------------------------------------------/CASE "MIX" --
  
  END SELECT ! CASE(PathMod3)
  
  !--------------------- keep record of phase found during whole path --
  DO iFs=1,SIZE(vFas) !
    IF(vFas(iFs)%Mole>Zero) vFasFound(iFs)=.TRUE.
  ENDDO
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Path_NewStep"
  !
END SUBROUTINE Path_NewStep

SUBROUTINE NewStep_Adjust_TotO_TotH( & !
& AdjustMode, &                        !
& iBal,ic_W,ic_H,vEle, &               !
& vCpn)
!--
!-- adjust total molar amounts of O and H according to molar amounts of other elements
!-- -> values of vTotF
!--
  USE M_T_Element,  ONLY: T_Element
  USE M_T_Component,ONLY: T_Component
  !
  INTEGER,          INTENT(IN)   :: AdjustMode
  INTEGER,          INTENT(IN)   :: iBal
  INTEGER,          INTENT(IN)   :: ic_W,ic_H
  TYPE(T_Element),  INTENT(IN)   :: vEle(:)
  TYPE(T_Component),INTENT(INOUT):: vCpn(:)
  !
  REAL(dp):: Zbalance,TotO_
  INTEGER :: I,iBal_
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< NewStep_Adjust_TotO_TotH"
  !
  TotO_= vCpn(ic_W)%Mole
  !
  IF(iBal==0) THEN ;  iBal_= ic_H
  ELSE             ;  iBal_= iBal
  ENDIF

  IF(iBal_ /= ic_H .AND. vCpn(ic_H)%Statut=="INERT") THEN
    vCpn(ic_W)%Mole= TotO_
    vCpn(ic_H)%Mole= TotO_*Two
    !RETURN
  ENDIF
  
  Zbalance= Zero
  DO I=1,SIZE(vCpn)
    IF(vCpn(I)%Statut=="INERT" .AND. I/=ic_W .AND. I/=ic_H) THEN
      !print *,vEle(vCpn(I)%iEle)%NamEl, vEle(vCpn(I)%iEle)%Z
      Zbalance= Zbalance + vEle(vCpn(I)%iEle)%Z *vCpn(I)%Mole
    ENDIF
  ENDDO
  !pause

  SELECT CASE(AdjustMode)
  
  CASE(1)
    !--- IF excess cation, equiv'Na > equiv'Cl, adjust with OH
    IF(Zbalance > Zero) THEN
      vCpn(ic_W)%Mole=  TotO_     +Zbalance
      vCpn(iBal_)%Mole= TotO_*Two +Zbalance
    !--- IF excess anion, Cl>Na, adjust with H
    ELSE
      vCpn(ic_W)%Mole=  TotO_
      vCpn(iBal_)%Mole= (TotO_*Two -Zbalance) /REAL(vEle(iBal_)%Z)
    ENDIF

  CASE(2)
    !--- Oxygen number unchanged, equilibrate using hydrogen only
    vCpn(ic_W)%Mole= TotO_
    vCpn(iBal_)%Mole= TotO_*Two -Zbalance
    
  END SELECT
  
  !write(12,'(3G15.6)') Zbalance,vCpn(ic_W)%Mole,vCpn(iBal_)%Mole

  !~ IF(iBal /= ic_H) THEN
    !~ vCpn(ic_W)%Mole= TotO_
    !~ vCpn(ic_H)%Mole= TotO_*Two
    !~ !RETURN
  !~ ENDIF
  
  !~ Zbalance= Zero
  !~ DO I=1,SIZE(vCpn)
    !~ IF(vCpn(I)%Statut=="INERT" .AND. I/=ic_W .AND. I/=ic_H) &
    !~ & Zbalance= Zbalance + vEle(vCpn(I)%iEle)%Z *vCpn(I)%Mole
  !~ ENDDO
  !~ !
  !~ SELECT CASE(AdjustMode)

  !~ CASE(1)
    !~ !--- IF excess cation, Na>Cl, adjust with OH
    !~ IF(Zbalance > Zero) THEN
      !~ vCpn(ic_W)%Mole= TotO_     +Zbalance
      !~ vCpn(iBal)%Mole= TotO_*Two +Zbalance
    !~ !--- IF excess anion, Cl>Na, adjust with H
    !~ ELSE
      !~ vCpn(ic_W)%Mole=  TotO_
      !~ vCpn(iBal)%Mole= (TotO_*Two -Zbalance) /REAL(vEle(iBal)%Z)
    !~ ENDIF

  !~ CASE(2)
    !~ !--- Oxygen number unchanged, equilibrate using hydrogen only
    !~ vCpn(ic_W)%Mole= TotO_
    !~ vCpn(iBal)%Mole= TotO_*Two -Zbalance

  !~ ENDSELECT
  !
  IF(iDebug>2) THEN
    WRITE(fTrc,'(A,G15.6)')  "Zbalance       =", Zbalance
    WRITE(fTrc,'(A,G15.6)')  "vCpn(ic_W)%Mole=", vCpn(ic_W)%Mole
    WRITE(fTrc,'(A,G15.6)')  "vCpn(ic_H)%Mole=", vCpn(ic_H)%Mole
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ NewStep_Adjust_TotO_TotH"
  !
ENDSUBROUTINE NewStep_Adjust_TotO_TotH
  !
SUBROUTINE Path_Calc_FluidMix(TdgKMix,PbarMix,Ok)
!--
!-- read constraints on mix fluid,
!-- then calc' speciation, then make all cpn inert,
!--
  USE M_T_Component, ONLY: T_Component
  USE M_Basis,       ONLY: Basis_Change
  USE M_System_Tools,ONLY: System_Build_Custom,System_TP_Update
  USE M_Equil
  !
  USE M_Global_Vars, ONLY: vEle,vSpc,vMixFas
  USE M_System_Vars, ONLY: vCpn,TdgK,Pbar
  !
  REAL(dp),INTENT(OUT):: TdgKMix,PbarMix
  LOGICAL, INTENT(OUT):: Ok
  !
  TYPE(T_Component),DIMENSION(:),ALLOCATABLE:: vC
  INTEGER :: I,J,N
  REAL(dp):: T0,P0
  !
  ALLOCATE(vC(1:SIZE(vCpn)))
  vC= vCpn
  T0= TdgK
  P0= Pbar
  !
  !IF(ALLOCATED(vCpnMix)) DEALLOCATE(vCpnMix)
  vCpnMix= vCpn
  !
  CALL System_Build_Custom( &
  & "SYSTEM.MIX",vEle,vSpc,vMixFas,vCpn, & !IN
  & TdgK,Pbar, &
  & vCpnMix,Ok)
  !
  IF(Ok) THEN
    vCpn= vCpnMix
    !
    CALL System_TP_Update(TdgK,Pbar)
    
    CALL Equil_Calc("SPC")
    !-> calc' speciation for mixing end-member system
    CALL Basis_Change("EQU",vCpn)
    !-> make all components inert
    vCpnMix= vCpn
  ENDIF
  !
  TdgKMix= TdgK
  PbarMix= Pbar
  !
  vCpn= vC
  TdgK= T0
  Pbar= P0
  !
  !! IF(PbarMix/=Pbar) THEN
  !!   CALL Warning_("pressure should be the same as master system !!")
  !!   PbarMix= Pbar
  !! ENDIF
  !
  IF(Ok) THEN
    !permute vCpnMix -> its element order must be same as vCpnBox
    !-> USE vC for swapping
    vC= vCpnMix
    N= SIZE(vCpn)
    DO I=1,N
      DO J=1,N; IF(vC(J)%iEle==vCpn(I)%iEle) EXIT; ENDDO
      vCpnMix(I)= vC(J)
    ENDDO
  ENDIF
  DEALLOCATE(vC)
  !
  RETURN
END SUBROUTINE Path_Calc_FluidMix

SUBROUTINE Path_Record
  USE M_Global_Vars,ONLY: vFas
  USE M_Path_Vars,  ONLY: vFasFound
  !  
  INTEGER:: iFs
  DO iFs=1,SIZE(vFas) !keep record of phase found during whole path
    IF(vFas(iFs)%Mole>Zero) vFasFound(iFs)=.TRUE.
  ENDDO
  !
END SUBROUTINE Path_Record

ENDMODULE M_Path

