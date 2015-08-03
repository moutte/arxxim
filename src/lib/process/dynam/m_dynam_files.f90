MODULE M_Dynam_Files
  USE M_Kinds
  USE M_Trace,ONLY: fTrc,fHtm,T_,iDebug
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Dynam_Files_Init
  PUBLIC:: Dynam_Files_Close
  PUBLIC:: Dynam_Files_OpenLogs
  PUBLIC:: Dynam_Files_WriteFasAff
  !
  INTEGER,PUBLIC:: &
  & fDynAct= 0, & !results: log10(activities) 
  & fDynMol= 0, & !mole numbers
  & fDynEle= 0, & !component mole numbers (different from element, in orthogonal option)
  & fDynElement= 0, & !element mole numbers
  & fDynQsK= 0, & !minerals
  & fDynMnK= 0, & !minerals
  & fDynGam= 0    !gammas

CONTAINS

SUBROUTINE Dynam_Files_Init
  USE M_Files,      ONLY: DirOut,cTitle,Files_Index_Write
  USE M_Numeric_Tools,ONLY: fNewtF,fNewtR,fNewt_I
  USE M_IoTools,    ONLY: GetUnit
  USE M_Equil_Write,ONLY: Equil_Write_EnTete
  !
  USE M_Global_Vars,ONLY: vEle,vSpc
  USE M_System_Vars,ONLY: vCpn
  USE M_Basis_Vars, ONLY: vOrdAq,vPrmFw
  !
  USE M_Dynam_Vars, ONLY: TUnit,dTSav,fSavTime,fSavRate,DebNewt,DebJacob
  !
  INTEGER:: I,nAq
  !
  nAq= COUNT(vSpc%Typ(1:3)=="AQU")
  !
  IF(fDynGam==0) THEN
    !
    CALL GetUnit(fDynGam)
    OPEN(fDynGam,FILE=TRIM(DirOut)//"_gamma.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_gamma.restab", &
    & "DYNAMIC: activity coeff's of aqueous species")
    !
    WRITE(fDynGam,'(2(A,A1))',ADVANCE='NO') "step",T_,"time",T_
    CALL Equil_Write_EnTete(fDynGam,vSpc,vPrmFw,vSpc%Typ=="AQU") 
    !
  ENDIF
  !
  IF(fDynAct==0) THEN
    !
    CALL GetUnit(fDynAct)
    OPEN(fDynAct,FILE=TRIM(DirOut)//"_activ.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_activ.restab", &
    & "DYNAMIC: at each time step, species activities")
    !
    WRITE(fDynAct,'(2(A,A1))',ADVANCE='NO') "STEP",T_,"TIME/"//TUnit,T_
    DO I=1,nAq; WRITE(fDynAct,'(A12,A1)',ADVANCE='NO') vSpc(vOrdAq(I))%NamSp,T_; ENDDO
    WRITE(fDynAct,*)
    !
  ENDIF
  !
  IF(fDynMol==0) THEN
    !
    CALL GetUnit(fDynMol)
    OPEN(fDynMol,FILE=TRIM(DirOut)//"_moleinbox.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_moleinbox.restab", &
    & "DYNAMIC: mole numbers of aqu.species in box")
    !
    WRITE(fDynMol,'(2(A,A1))',ADVANCE='NO') "STEP",T_,"TIME/"//TUnit,T_
    DO I=1,nAq; WRITE(fDynMol,'(A,A1)',ADVANCE='NO') TRIM(vSpc(vOrdAq(I))%NamSp),T_; ENDDO
    WRITE(fDynMol,*)
    !
  ENDIF
  !
  IF(iDebug>2 .AND. fDynElement==0) THEN
    !
    CALL GetUnit(fDynElement)
    OPEN(fDynElement,FILE=TRIM(DirOut)//"_elements.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_elements.restab", &
    & "DYNAMIC: at each time step,"//  &
    & " mole numbers of elements within Box, in Fluid and in Minerals")
    !
    !-- element order will be the order in tStoikioAqu / tStoikioKin,
    !-- which are retrieved from vSpc(:)%vStoikio(:)
    !-- it is thus the element order in original vEle(:)
    !
    WRITE(fDynElement,'(2(A,A1))',ADVANCE='NO') "STEP",T_,"TIME/"//TUnit,T_
    DO I=1,SIZE(vEle); WRITE(fDynElement,'(A,A1)',ADVANCE='NO') &
    & vEle(I)%NamEl//"inFluid", T_; ENDDO
    DO I=1,SIZE(vEle); WRITE(fDynElement,'(A,A1)',ADVANCE='NO') &
    & vEle(I)%NamEl//"inMin", T_; ENDDO
    DO I=1,SIZE(vEle); WRITE(fDynElement,'(A,A1)',ADVANCE='NO') &
    & vEle(I)%NamEl//"Total", T_; ENDDO
    WRITE(fDynElement,*)
    !
  ENDIF
  !
  IF(fDynEle==0) THEN
    !
    CALL GetUnit(fDynEle)
    OPEN(fDynEle,FILE=TRIM(DirOut)//"_elem.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_elem.restab", &
    & "DYNAMIC: at each time step,"//  &
    & " mole numbers of elements within Box, in Fluid and in Minerals")
    !
    !-- this header will be valid only when component--element !
    !-- -> to be modified in other cases,
    !-- e.g. when using orthogonal components !!
    !
    WRITE(fDynEle,'(4(A,A1))',ADVANCE='NO') "STEP",T_,"TIME/"//TUnit,T_,"DeltaDarcy",T_
    !~ DO I=1,SIZE(vEle); WRITE(fDynEle,'(A,A1)',ADVANCE='NO') vEle(iCpnEle(I))%NamEl//"totF", T_; ENDDO
    !~ DO I=1,SIZE(vEle); WRITE(fDynEle,'(A,A1)',ADVANCE='NO') vEle(iCpnEle(I))%NamEl//"molal",T_; ENDDO
    !~ DO I=1,SIZE(vEle); WRITE(fDynEle,'(A,A1)',ADVANCE='NO') vEle(iCpnEle(I))%NamEl//"inMin",T_; ENDDO
    DO I=1,SIZE(vCpn); WRITE(fDynEle,'(A,A1)',ADVANCE='NO') TRIM(vCpn(I)%NamCp)//"totF", T_; ENDDO
    DO I=1,SIZE(vCpn); WRITE(fDynEle,'(A,A1)',ADVANCE='NO') TRIM(vCpn(I)%NamCp)//"molal",T_; ENDDO
    DO I=1,SIZE(vCpn); WRITE(fDynEle,'(A,A1)',ADVANCE='NO') TRIM(vCpn(I)%NamCp)//"inMin",T_; ENDDO
    WRITE(fDynEle,*)
    !
  ENDIF
  !
  IF(fDynMnK==0) THEN
    !
    CALL GetUnit(fDynMnK)
    OPEN(fDynMnK,FILE=TRIM(DirOut)//"_minmol.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_minmol.restab", &
    & "DYNAMIC: at each time step, pH, vol'fractions fluid and minerals, logQsK, etc")
    !
    CALL EnTeteFMnk(fDynMnK) !,bCell)
    !
  ENDIF
  !
  IF(dTSav>Zero) THEN
    CALL GetUnit(fSavTime)
    OPEN(fSavTime,FILE=TRIM(DirOut)//"_time.restab")
    !
    CALL EnTeteFMnk(fSavTime) !,bCell)
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_time.restab", &
    & "DYNAMIC: at regular time laps, pH, vol'fractions fluid and minerals, logQsK, etc")
  ENDIF
  !
  IF(iDebug>2) THEN
    IF(fSavRate==0) THEN
      !
      CALL GetUnit(fSavRate)
      OPEN(fSavRate,FILE=TRIM(DirOut)//"_rate.restab")
      !
      CALL EnTeteFRate(fSavRate)
      !
      CALL Files_Index_Write(fHtm, &
      & TRIM(DirOut)//"_rate.restab", &
      & "DYNAMIC: at each time step, details on mineral dissol/precip rates (surface,...)")
     ENDIF
  ENDIF
  !
  IF(DebNewt) THEN
    !
    !fNewtF:  trace for Function values on aqu'species
    !fNewtR:  trace for Residual values on aqu'species
    !fNewtFm: trace for Function values on minerals
    !fNewtRm: trace for Residual values on minerals
    !
    CALL GetUnit(fNewtF);  OPEN(fNewtF, FILE=TRIM(DirOut)//"_newtaqu1.log")
    CALL GetUnit(fNewtR);  OPEN(fNewtR, FILE=TRIM(DirOut)//"_newtaqu2.log") !trace for Residual values
    !CALL GetUnit(fNewtFm); OPEN(fNewtFm,FILE=TRIM(DirOut)//"_newtmin1.log")
    !CALL GetUnit(fNewtRm); OPEN(fNewtRm,FILE=TRIM(DirOut)//"_newtmin2.log")
    !
    fNewt_I= 0
    !
    WRITE(fNewtF,'(A,A1)',ADVANCE='NO') "Index",T_
    CALL Equil_Write_EnTete(fNewtF,vSpc,vPrmFw,vSpc%Typ=="AQU",Str1="indx") 
    !
    WRITE(fNewtR,'(A,A1)',ADVANCE='NO') "Index",T_
    CALL Equil_Write_EnTete(fNewtR,vSpc,vPrmFw,vSpc%Typ=="AQU",Str1="indx") 
    !
    DebJacob= .TRUE.
    !
  ENDIF
  !
ENDSUBROUTINE Dynam_Files_Init

SUBROUTINE EnTeteFMnk(f) !,bCell)
  USE M_Files, ONLY: cTitle
  USE M_Global_Vars,ONLY: vKinFas
  USE M_Dynam_Vars, ONLY: TUnit
  !
  INTEGER,INTENT(IN):: f
  !!LOGICAL,INTENT(IN):: bCell
  !
  INTEGER::i,N
  !
  N=SIZE(vKinFas)
  !
  !WRITE(f,'(4(A,A1))',ADVANCE='NO') "!iStep", T_,"Time",T_,"pH",T_,"PhiFluid",T_
  !DO J=1,3
  !  DO i=1,N; WRITE(f,'(A,A1)', ADVANCE='NO') TRIM(vKinFas(i)%Name),T_; ENDDO
  !ENDDO
  !WRITE(f,*)
  !
  !WRITE(f,'(4(A,A1))',ADVANCE='NO') "!iStep", T_,"Time",T_,"pH",T_,"PhiFluid",T_
  !DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "PhiM",T_;  ENDDO
  !DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "LogQsK",T_;  ENDDO
  !DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "MolNr",T_;  ENDDO
  !WRITE(f,*)
  !
  !!IF(bCell) WRITE(f,'(A,A1)',ADVANCE='NO') "iCell",T_
  WRITE(f,'(4(A,A1))',ADVANCE='NO') "iStep", T_,"Time/"//TUNit,T_,"pH",T_,"PhiFluid",T_
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "PhiM_"//TRIM(vKinFas(i)%NamKF),T_;   ENDDO
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "LogQsK_"//TRIM(vKinFas(i)%NamKF),T_; ENDDO
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "MolNr"//TRIM(vKinFas(i)%NamKF),T_;   ENDDO
  WRITE(f,*)
ENDSUBROUTINE EnTeteFMnk

SUBROUTINE EnTeteFRate(f)
  USE M_Files, ONLY: cTitle
  USE M_Global_Vars,ONLY: vKinFas
  !
  INTEGER,INTENT(IN)::f
  !
  INTEGER::i,N
  !
  N=SIZE(vKinFas)
  !
  WRITE(f,'(4(A,A1))',ADVANCE='NO') "iStep", T_,"Time",T_,"pH",T_,"PhiFluid",T_
  !
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "Phi_"    //TRIM(vKinFas(i)%NamKF),T_; ENDDO
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "SurfR_"  //TRIM(vKinFas(i)%NamKF),T_; ENDDO
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "SurfKg_" //TRIM(vKinFas(i)%NamKF),T_; ENDDO
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "lgQsK_"  //TRIM(vKinFas(i)%NamKF),T_; ENDDO
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "ratQsK_" //TRIM(vKinFas(i)%NamKF),T_; ENDDO
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "ratAct_" //TRIM(vKinFas(i)%NamKF),T_; ENDDO
  DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "Radius_" //TRIM(vKinFas(i)%NamKF),T_; ENDDO
  !DO i=1,N; WRITE(f,'(A,A1)',ADVANCE='NO') "SferNum_"//TRIM(vKinFas(i)%NamKF),T_; ENDDO
  WRITE(f,*)
  !
ENDSUBROUTINE EnTeteFRate

SUBROUTINE Dynam_Files_Close
  USE M_Trace,ONLY: DebugCoores
  USE M_Numeric_Tools,ONLY: fNewtF,fNewtR
  USE M_Dynam_Vars, ONLY: fSavTime,fSavRate
  !! USE M_Stockvar_Kinxim,ONLY: LSTOCK,DEL_STOCKVAR,PRINT_STOCKVAR !,SET_STOCKVAR,INIT_STOCKVAR
  !
  IF(fSavTime>0) THEN; CLOSE(fSavTime); fSavTime= 0; ENDIF
  IF(fSavRate>0) THEN; CLOSE(fSavRate); fSavRate= 0; ENDIF
  !
  IF(fDynMol>0)  THEN; WRITE(fDynMol,*)  ; CLOSE(fDynMol);  fDynMol= 0; ENDIF 
  IF(fDynAct>0)  THEN; WRITE(fDynAct,*)  ; CLOSE(fDynAct);  fDynAct= 0; ENDIF
  IF(fDynGam>0)  THEN; WRITE(fDynAct,*)  ; CLOSE(fDynGam);  fDynGam= 0; ENDIF
  !
  IF(fNewtF>0)   THEN; CLOSE(fNewtF);  fNewtF=  0; ENDIF
  IF(fNewtR>0)   THEN; CLOSE(fNewtR);  fNewtR=  0; ENDIF
  !
  !IF(fNewtFm>0)  THEN; CLOSE(fNewtFm); fNewtFm= 0; ENDIF
  !IF(fNewtRm>0)  THEN; CLOSE(fNewtRm); fNewtRm= 0; ENDIF
  !
  IF(fDynEle>0)  THEN; CLOSE(fDynEle); fDynEle=   0; ENDIF 
  IF(fDynMnK>0)  THEN; CLOSE(fDynMnK); fDynMnK=    0; ENDIF
  !
  IF(fDynElement>0)  THEN; CLOSE(fDynElement); fDynElement=   0; ENDIF 
  
  !---del table to store the time evolution of species 
  !--- stockage solution
  !! IF (LSTOCK) THEN
  !!    IF ( DebugCoores ) CALL Print_Stockvar('output.var');
  !!    CALL Del_Stockvar
  !! END IF
  
ENDSUBROUTINE Dynam_Files_Close

SUBROUTINE Dynam_Files_OpenLogs(f1,f2,f3,f4)
  USE M_IOTools,ONLY: GetUnit
  USE M_Files,  ONLY: DirOut,DirLog,Files_Index_Write
  USE M_Global_Vars,ONLY: vEle,vKinFas
  !~ USE M_Basis_Vars, ONLY: iCpnEle
  USE M_Dynam_Vars,ONLY: vCpnBox
  USE M_Dynam_Vars, ONLY: TUnit
  !
  INTEGER,INTENT(OUT):: f1,f2,f3,f4
  !
  INTEGER:: I
  !
  !----------------------------------------------------------------------------!
  CALL GetUnit(f1)
  !!OPEN(f1,FILE=TRIM(DirLog)//"calcdyn_totout.log")
  OPEN(f1,FILE=TRIM(DirOut)//"_bilans.restab")
  !
  CALL Files_Index_Write(fHtm,&
  !! & TRIM(DirLog)//"calcdyn_totout.log",& !
  & TRIM(DirOut)//"_bilans.restab", &
  & "DYNAMIC/LOG: balances on elements")
  !
  WRITE(f1,'(2(A,A1))',ADVANCE='NO') "STEP",T_,"TIME/"//TUnit,T_
  DO I=1,SIZE(vCpnBox)
    WRITE(f1,'(A,A1)',ADVANCE='NO') &
    & TRIM(vEle(vCpnBox(I)%iEle)%NamEl)//'inj',T_; ENDDO 
  DO I=1,SIZE(vCpnBox)
    WRITE(f1,'(A,A1)',ADVANCE='NO') &
    & TRIM(vEle(vCpnBox(I)%iEle)%NamEl)//'bilan',T_; ENDDO 
  DO I=1,SIZE(vCpnBox)
    WRITE(f1,'(A,A1)',ADVANCE='NO') &
    & TRIM(vEle(vCpnBox(I)%iEle)%NamEl)//'diff',T_; ENDDO 
  WRITE(f1,*)
  !
  !--------------------------------------------------------------------------!
  CALL GetUnit(f2)
  !OPEN(f2,FILE=TRIM(DirLog)//"calcdyn_steps.log")
  OPEN(f2,FILE=TRIM(DirOut)//"_calcdyn_steps.log")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_calcdyn_steps.log",&
  & "DYNAMIC/LOG: iMaxDelta,VarMax,iDo1,iDo2,etc.")
  !
  WRITE(f2,'(14(A,A1))') &
  & ".MaxVal",T_,"iStep",T_,"Time",T_,"dTime",T_,"PhiF",T_,&
  & "iMaxDelta",T_,"VarMAx",T_,&
  & "iDo1",T_,"iDo2",T_,"Newt_iDo",T_,&
  & "Newt_iErr",T_,"NewtErrF",T_,"NewtErrX",T_,"NewtErrG",T_
  !
  !--------------------------------------------------------------------------!
  CALL GetUnit(f3)
  !OPEN(f3,FILE=TRIM(DirLog)//"calcdyn_maxvars.log")
  OPEN(f3,FILE=TRIM(DirOut)//"_calcdyn_maxvars.log")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//"_calcdyn_maxvars.log",&
  & "DYNAMIC/LOG: MaxAqu,MaxMin,MaxAquRel,MaxMinRel,...")
  !
  WRITE(f3,'(9(A,A1))') &
  & ".MaxAqu",T_,".MaxMin",T_,&
  & "iStep",T_,"Time",T_,"dTime",T_,&
  & "MaxAqu",T_,"MaxMin",T_,&
  & "MaxAquRel",T_,"MaxMinRel",T_
  !
  !--------------------------------------------------------------------------!
  CALL GetUnit(f4)
  OPEN(f4,FILE=TRIM(DirLog)//"calcdyn_satur.log")
  !
  RETURN
ENDSUBROUTINE Dynam_Files_OpenLogs

SUBROUTINE Dynam_Files_WriteFasAff(iStep,Time,vLnAct,vLnBuf)
!.WRITE affinity of pure phases (i.e. non-aqu'species)
  USE M_IOTools
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Files,        ONLY: cTitle,DirOut,Files_Index_Write
  USE M_Basis_Vars,   ONLY: vOrdPr,tNuFas,nCi
  USE M_Global_Vars,  ONLY: vFas,vSpc
  !
  INTEGER, INTENT(IN):: iStep
  REAL(dp),INTENT(IN):: Time
  REAL(dp),INTENT(IN):: vLnAct(:)
  REAL(dp),INTENT(IN):: vLnBuf(:)
  !
  INTEGER :: iFs,nFs,nCp,nCx
  REAL(dp):: X
  !
  nFs= SIZE(vFas)
  nCp= SIZE(vOrdPr)
  nCx= nCp -nCi
  !
  IF(iStep>9999) RETURN
  !
  IF (nFs>0) THEN
    !
    !------------------------------------------- open fDynQsK, write en tete --!
    IF(fDynQsK==0) THEN
      CALL GetUnit(fDynQsK)
      OPEN(fDynQsK,FILE=TRIM(DirOut)//"_minqsk.restab")
      CALL Files_Index_Write(fHtm, &
      & TRIM(DirOut)//"_minqsk.restab", &
      & "DYNAMIC: Q/K of all phases from DATAbase")
      WRITE(fDynQsK,'(2(A,A1))',ADVANCE= 'NO') "Step",T_,"Time",T_ !,"pH",T_
      DO iFs=1,nFs
        WRITE(fDynQsK,'(A,A1)',ADVANCE= 'NO') TRIM(vFas(iFs)%NamFs),T_
      ENDDO
      WRITE(fDynQsK,*)
    ENDIF
    !--------------------------------------------------------------------------!
    !
    WRITE(fDynQsK,'(I5,A1,G15.6,A1)',ADVANCE= 'NO') iStep,T_,Time,T_
    DO iFs=1,nFs
      X= &
      & vFas(iFs)%Grt &
      - DOT_PRODUCT( tNuFas(iFs,1:nCi), &
      &              vLnAct(1:nCi) +vSpc(vOrdPr(1:nCi))%G0rt )
      IF(nCx>0) &
      & X= X &
      &  - DOT_PRODUCT( tNuFas(iFs,nCi+1:nCp), &
      &              vLnBuf(1:nCx) +vSpc(vOrdPr(nCi+1:nCp))%G0rt )
      !
      WRITE(fDynQsK,'(F12.3,A1)',ADVANCE='NO') -X /Ln10,T_
      !                                       log10(QsK)= - Affinity /Ln10
    ENDDO
    WRITE(fDynQsK,*)
  ENDIF
ENDSUBROUTINE Dynam_Files_WriteFasAff    

ENDMODULE M_Dynam_Files
