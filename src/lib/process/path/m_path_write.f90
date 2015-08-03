MODULE M_Path_Write
  USE M_Kinds
  USE M_Trace,ONLY: T_,fHtm !Stop_,fTrc,T_,iDebug,Pause_
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC:: Path_Write_Line
  PUBLIC:: Path_Write_FasAff
  PUBLIC:: Path_Write_FasEnTete
  PUBLIC:: Path_Write_Distrib
  PUBLIC:: Path_Files_Close
  
  INTEGER::     & !
  & fQsK=    0, & !saturation states
  & fActiv=  0, & !log10(activities)
  & fPoten=  0, & !potentials /RT, Mu/RT(:)= Mu°(:)/RT + ln(act(:))
  & fMoles=  0, & !mole numbers
  & fMolal=  0, & !molality solutes
  & fGamma=  0    !gammas
  !
  LOGICAL:: bOpenDist=.FALSE.
  !
  CHARACTER(LEN=15),ALLOCATABLE:: vNamFasYes(:)

CONTAINS

SUBROUTINE Path_Write_FasEnTete
!--
!-- write entete for the file of affinities of pure phases
!-- (i.e. non-aqu'species)
!--
  USE M_IOTools,   ONLY: GetUnit
  USE M_Files,     ONLY: cTitle,DirOut, Files_Index_Write
  USE M_Global_Vars,ONLY: vFas
  !
  INTEGER :: iFs
  !
  CALL GetUnit(fQsK)
  OPEN(fQsK,FILE=TRIM(DirOut)//"_minqsk.restab")
  CALL Files_Index_Write(fHtm, &
  & TRIM(DirOut)//"_minqsk.restab", &
  & "Q/K of all phases along path ")
  !
  WRITE(fQsK,'(A,A1)',ADVANCE= 'NO') "Count",T_ !,"pH",T_
  DO iFs=1,SIZE(vFas)
    WRITE(fQsK,'(A,A1)',ADVANCE= 'NO') TRIM(vFas(iFs)%NamFs),T_
  ENDDO
  WRITE(fQsK,*)
  !
END SUBROUTINE Path_Write_FasEnTete

SUBROUTINE Path_Write_FasAff(iStep)
!--
!-- write affinity of pure phases (i.e. non-aqu'species)
!--
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Global_Vars,  ONLY: vFas,vSpc
  USE M_Basis_Vars,   ONLY: tNuFas,vOrdPr
  !
  INTEGER,INTENT(IN):: iStep
  !
  INTEGER :: iFs
  REAL(dp):: X
  !
  WRITE(fQsK,'(I4,A1)',ADVANCE= 'NO') iStep,T_
  DO iFs=1,SIZE(vFas)
    X= &
    & vFas(iFs)%Grt &
    - DOT_PRODUCT( tNuFas(iFs,:), &
    &              vSpc(vOrdPr(:))%Dat%LAct+vSpc(vOrdPr(:))%G0rt )
    X= - X /SUM(ABS(tNuFas(iFs,:))) /Ln10
    !!WRITE(fQsK,'(G15.6,A1)',ADVANCE='NO') -vFasAff(iFs) /Ln10,T_ !-> log10(QsK)= - Affinity /Ln10
    WRITE(fQsK,'(G15.6,A1)',ADVANCE='NO') X,T_ !-> log10(QsK)= - Affinity /Ln10
  ENDDO
  WRITE(fQsK,*)
  !
END SUBROUTINE Path_Write_FasAff

SUBROUTINE Path_Write_Distrib(N)
!--- distribution of all elements among their different species --
!--- to be used after Equil_Save (to update vSpc%Mole) --
!--
!Element iEl is distributed among all Species iAq that have tFormula(iEl,iAq)/=0
!IF Speciation succeeds, 
!THEN we should have
!  vCpn(iEl)%Mole=DOt_product(tFormula(iEl,1:nAq),vSpc(1:nAq)%Dat%Mole)),
!except for element iBal or for a mobile element
!--
!Output format is designed for easy data processing w/ spreadsheet
!(sort according to column vNamEl, -> produce tables for the different elements)
!currently no special treatment for ReDOx !!! should DO something...
!--
  USE M_IoTools,ONLY: GetUnit
  USE M_Numeric_Const,  ONLY: Ln10
  USE M_SolModel_Tools,ONLY: Solmodel_pHpE
  !
  USE M_Global_Vars, ONLY: vEle,nAq,vSpc,tFormula
  USE M_System_Vars, ONLY: vCpn
  USE M_Files,       ONLY: DirOut
  USE M_Basis_Vars,  ONLY: isW,iOx,isH_,isO2
  !
  INTEGER ::fDist,N,iAq,iEl,nCp
  REAL(dp)::Tot, X, pH_,pE_
  !
  CALL Solmodel_pHpE(isW,iOx,isH_,isO2,vSpc,pH_,pE_)
  nCp=SIZE(vCpn)
  IF(.NOT.bOpenDist) THEN
    CALL EnTeteDistrib
    bOpenDist=.TRUE.
  ENDIF
  DO iEl=2,nCp
    CALL GetUnit(fDist)
    OPEN(fDist, &
    & FILE=TRIM(DirOut)//"_distrib_"//TRIM(vEle(vCpn(iEl)%iEle)%NamEl)//".restab", &
    & ACCESS='APPEND')
    Tot=DOT_PRODUCT(tFormula(iEl,1:nAq),vSpc(1:nAq)%Dat%Mole)
    IF(TRIM(vCpn(iEl)%Statut)/="INERT") vCpn(iEl)%Mole=Tot
    WRITE(fDist,'(I3,A1,A,A1,F7.3,A1,F7.3,A1)',ADVANCE='NO') &
    &             N,T_,TRIM(vEle(vCpn(iEl)%iEle)%NamEl),T_,pH_,T_,pE_,T_
    DO iAq=2,nAq
      !IF(iAq/=iH_) THEN
        IF(tFormula(iEl,iAq)/=0) THEN
          X=vSpc(iAq)%Dat%Mole/tFormula(iEl,iAq)/Tot*1000.0D0
          !divides vMolF of species iAq by number of moles of iEl in one mole of iAq
          !IF(X>1.0D-3) THEN; WRITE(fDist,'(F12.3,A1)',ADVANCE='NO') X,T_ 
          !             ELSE; WRITE(fDist,'(A12,A1)',  ADVANCE='NO') "0",T_; ENDIF 
          WRITE(fDist,'(G15.6,A1)',ADVANCE='NO') X,T_
        ENDIF
      !ENDIF
    ENDDO
    DO iAq=2,nAq
      IF(tFormula(iEl,iAq)/=0) &
      & WRITE(fDist,'(G15.6,A1)',ADVANCE='NO') 1.0E6*vSpc(iAq)%Dat%Mole,T_
      !*1.0E6 -> output in MICROMOLES !!!
    ENDDO
    DO iAq=2,nAq
      IF(tFormula(iEl,iAq)/=0) &
      & WRITE(fDist,'(G15.6,A1)',ADVANCE='NO') EXP(vSpc(iAq)%Dat%LGam),T_ !vSpc(iAq)%LnGam/Ln10
    ENDDO
    WRITE(fDist,'(2(G15.6,A1))') 1.0E6*Tot,T_,1.0E6*vCpn(iEl)%Mole
    CLOSE(fDist)
  ENDDO
END SUBROUTINE Path_Write_Distrib

SUBROUTINE EnteteDistrib
!--
!-- write head line for distribution of elements among species --
!--
  USE M_IoTools,     ONLY: GetUnit
  USE M_Global_Vars, ONLY: nAq,vEle,vSpc,tFormula
  USE M_System_Vars, ONLY: vCpn
  USE M_Files,       ONLY: cTitle,DirOut
  !
  INTEGER::fDist,iAq,iEl
  !
  DO iEl=2,SIZE(vCpn)
    !
    CALL GetUnit(fDist)
    OPEN(fDist,FILE=TRIM(DirOut)//"_distrib_"//TRIM(vEle(vCpn(iEl)%iEle)%NamEl//".restab"))
    !
    !IF(iEl/=iH_) THEN
    WRITE(fDist,'(A,/A)') &
    & "!."//TRIM(cTitle), &
    & "!.species distribution for a given element, relative (permil) and absolute (micromoles/kgH2O)"
    !
    WRITE(fDist,'(4(A,A1))',ADVANCE='NO') &
    &            "count",T_,TRIM(vEle(vCpn(iEl)%iEle)%NamEl),T_,"pH",T_,"pE",T_
    !
    DO iAq=2,nAq
      IF(tFormula(iEl,iAq)/=0) &
      &  WRITE(fDist,'(A,A1)',  ADVANCE='NO') TRIM(vSpc(iAq)%NamSp)//"_rel",T_
    ENDDO
    DO iAq=2,nAq
      IF(tFormula(iEl,iAq)/=0) &
      &  WRITE(fDist,'(A,A1)',  ADVANCE='NO') TRIM(vSpc(iAq)%NamSp),T_
    ENDDO
    DO iAq=2,nAq
      IF(tFormula(iEl,iAq)/=0) &
      &  WRITE(fDist,'(A,A1)',  ADVANCE='NO') TRIM(vSpc(iAq)%NamSp)//"_Gam",T_
    ENDDO
    !
    WRITE(fDist,'(A,A1,A)',ADVANCE='NO') "TOT,umoles",T_,"Check"
    WRITE(fDist,*)
    !
    !ENDIF
    CLOSE(fDist)
  ENDDO
  !
  RETURN
END SUBROUTINE EnTeteDistrib

SUBROUTINE Path_Write_Line(WrCount,WrCod,vYes,WrTitl)
!--
!-- one result on one line
!-- writes on several files:
!-- species log(activity) on _activ,
!-- species molalities on _molal,
!-- element abundance on _moles
!--
  USE M_IOTools
  USE M_Files,        ONLY: DirOut,cTitle,Files_Index_Write
  USE M_Dtb_Const,    ONLY: T_CK
  USE M_Dtb_Const,    ONLY: R_JK,Tref
  USE M_Numeric_Const,ONLY: Ln10
  USE M_SolModel_Tools, ONLY: Solmodel_pHpE,Solmodel_CalcMolal
  USE M_T_Species,    ONLY: Species_EntropyZero
  !
  USE M_Global_Vars,ONLY: vEle,vSpc,vFas,nAq,SolModel !!,Solvent
  USE M_System_Vars,ONLY: TdgK,Pbar,vCpn
  USE M_Basis_Vars, ONLY: isW,MWsv,iOx,isH_,isO2,tStoikio,tAlfFs
  USE M_Basis_Vars, ONLY: vLCi,vLAs,vLAx,vLMx !,nAx,nMx,vYesList
  USE M_Path_Vars,  ONLY: DimPath,tPathResults
  !
  INTEGER,              INTENT(IN):: WrCount
  CHARACTER(LEN=3),     INTENT(IN):: WrCod
  LOGICAL,     OPTIONAL,INTENT(IN):: vYes(:)
  CHARACTER(*),OPTIONAL,INTENT(IN):: WrTitl
  !
  LOGICAL,PARAMETER:: OutputInGrams = .FALSE.  !! .TRUE.
  !
  REAL(dp):: vMolal(1:nAq)
  REAL(dp):: ZBal,ZRdx,IonStr,xTotF
  REAL(dp):: pH_,pE_,X,nOx
  INTEGER :: iPr,I,N,K,nCp
  !
  REAL(dp),ALLOCATABLE:: vS0Ele(:)
  REAL(dp),ALLOCATABLE:: vPot(:)
  
  ALLOCATE(vS0Ele(SIZE(vSpc)))
  ALLOCATE(vPot(SIZE(vSpc)))
  !
  nCp= SIZE(vCpn)
  !
  DO I=1,SIZE(vSpc)
    !vS0Ele(I)= DOT_PRODUCT(vSpc(I)%vStoikio(1:N),vEle(1:N)%S0) /REAL(vSpc(I)%vStoikio(0))
    vS0Ele(I)= Species_EntropyZero(vEle,vSpc(I))
  ENDDO
  !
  CALL Solmodel_pHpE(isW,iOx,isH_,isO2,vSpc,pH_,pE_)
  !
  !------------------------------------------------ACTIV'COEFFS.HEADER--
  IF(fGamma==0) THEN
    !
    CALL GetUnit(fGamma)
    OPEN(fGamma,FILE=TRIM(DirOut)//"_gamma.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_gamma.restab", &
    & "PATH: log10(activity coeff's)")
    !
    WRITE(fGamma,'(A,A1)',   ADVANCE='NO') "Count",T_
    WRITE(fGamma,'(3(A,A1))',ADVANCE='NO') "TempC",T_,"Pbar",T_,"RhoW",T_
    WRITE(fGamma,'(3(A,A1))',ADVANCE='NO') "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !
    IF(iOx/=0) &
    & WRITE(fGamma,'(A2,A1,A5,A1)',ADVANCE='NO') "pE",T_,"RedOx",T_
    !
    !WRITE component names
    DO iPr=1,SIZE(vCpn) 
      WRITE(fGamma,'(A,A1)',ADVANCE='NO') TRIM(vEle(vCpn(iPr)%iEle)%NamEl),T_
    ENDDO
    !
    !WRITE species names (aqu'species ONLY)
    DO I=1,SIZE(vSpc)
      IF(vSpc(I)%Typ(1:3)=="AQU") & !!(vLCi(I).OR.vLAs(I).OR.vLAx(I))
      & WRITE(fGamma,'(A,A1)',ADVANCE='NO') TRIM(vSpc(I)%NamSp),T_
    ENDDO
    !
    WRITE(fGamma,*)
    !
  ENDIF
  !-----------------------------------------------/ACTIV'COEFFS.HEADER--
  !
  !--------------------------------------------------MOLALITIES.HEADER--
  IF(fMolal==0) THEN
    !
    CALL GetUnit(fMolal)
    OPEN(fMolal,FILE=TRIM(DirOut)//"_molal.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_molal.restab", &
    & "PATH: molalities")
    !
    WRITE(fMolal,'(A,A1)',   ADVANCE='NO') "Count",T_
    WRITE(fMolal,'(3(A,A1))',ADVANCE='NO') "TempC",T_,"Pbar",T_,"RhoW",T_
    WRITE(fMolal,'(3(A,A1))',ADVANCE='NO') "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !
    IF(iOx/=0) &
    & WRITE(fMolal,'(A2,A1,A5,A1)',ADVANCE='NO') "pE",T_,"RedOx",T_
    !
    !--- write component names
    DO iPr=1,SIZE(vCpn) 
      WRITE(fMolal,'(A,A1)',ADVANCE='NO') TRIM(vEle(vCpn(iPr)%iEle)%NamEl),T_
    ENDDO
    !
    !--- write species names (aqu'species ONLY)
    DO I=1,SIZE(vSpc)
      IF(vSpc(I)%Typ(1:3)=="AQU") & !!(vLCi(I).OR.vLAs(I).OR.vLAx(I))
      & WRITE(fMolal,'(A,A1)',ADVANCE='NO') TRIM(vSpc(I)%NamSp),T_
    ENDDO
    !
    WRITE(fMolal,*)
    !
  ENDIF
  !-------------------------------------------------/MOLALITIES.HEADER--
  !
  !--------------------------------------------------POTENTIALS.HEADER--
  IF(fPoten==0) THEN
    !
    CALL GetUnit(fPoten)
    OPEN(fPoten,FILE=TRIM(DirOut)//"_potent.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_potent.restab", &
    & "PATH result: log10(species activities)")
    !
    IF(PRESENT(WrTitl)) WRITE(fPoten,'(A,A1)',ADVANCE='NO') ".Title",T_
    WRITE(fPoten,'(A,A1)',   ADVANCE='NO') "Count",T_
    WRITE(fPoten,'(3(A,A1))',ADVANCE='NO') "TempC",T_,"Pbar",T_,"RhoW",T_
    WRITE(fPoten,'(3(A,A1))',ADVANCE='NO') "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !
    IF(iOx/=0) &
    & WRITE(fPoten,'(A2,A1,A5,A1)',ADVANCE='NO') "pE",T_,"RedOx",T_
    !
    !--- write component names
    DO iPr=1,SIZE(vCpn) 
      WRITE(fPoten,'(A,A1)',ADVANCE='NO') TRIM(vEle(vCpn(iPr)%iEle)%NamEl),T_
    ENDDO
    !
    !--- write species names (all aqu' + mobile non'aqu')
    DO I=1,SIZE(vSpc)
      IF(vLCi(I).OR.vLAs(I).OR.vLAx(I).OR.vLMx(I)) &
      & WRITE(fPoten,'(A,A1)',ADVANCE='NO') TRIM(vSpc(I)%NamSp),T_
    ENDDO
    WRITE(fPoten,*)
    !
    !~ !-------------------------------------------------------------------
    !~ !
    !~ IF(PRESENT(WrTitl)) WRITE(fPoten,'(A,A1)',ADVANCE='NO') ".Title",T_
    !~ WRITE(fPoten,'(A,A1)',   ADVANCE='NO') "Count",T_
    !~ WRITE(fPoten,'(3(A,A1))',ADVANCE='NO') "TempC",T_,"Pbar",T_,"RhoW",T_
    !~ WRITE(fPoten,'(3(A,A1))',ADVANCE='NO') "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !~ !
    !~ IF(iOx/=0) &
    !~ & WRITE(fPoten,'(A2,A1,A5,A1)',ADVANCE='NO') "pE",T_,"RedOx",T_
    !~ !
    !~ !--- write component names
    !~ DO iPr=1,SIZE(vCpn) 
      !~ WRITE(fPoten,'(A,A1)',ADVANCE='NO') TRIM(vEle(vCpn(iPr)%iEle)%NamEl),T_
    !~ ENDDO
    !~ !
    !~ DO I=1,SIZE(vSpc)
      !~ IF(vLCi(I).OR.vLAs(I).OR.vLAx(I).OR.vLMx(I)) &
      !~ & WRITE(fPoten,'(G15.6,A1)',ADVANCE='NO') vSpc(I)%G0rt,T_
    !~ ENDDO
    !~ WRITE(fPoten,*)
    !
  ENDIF
  !-----------------------------------------------/ POTENTIALS.HEADER --
  !
  !------------------------------------------------ ACTIVITIES.HEADER --
  IF(fActiv==0) THEN
    !
    CALL GetUnit(fActiv)
    OPEN(fActiv,FILE=TRIM(DirOut)//"_activ.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_activ.restab", &
    & "PATH result: log10(species activities)")
    !
    IF(PRESENT(WrTitl)) WRITE(fActiv,'(A,A1)',ADVANCE='NO') ".Title",T_
    WRITE(fActiv,'(A,A1)',   ADVANCE='NO') "Count",T_
    WRITE(fActiv,'(3(A,A1))',ADVANCE='NO') "TempC",T_,"Pbar",T_,"RhoW",T_
    WRITE(fActiv,'(3(A,A1))',ADVANCE='NO') "IonStr",T_,"pH",T_,"ChargeDiff",T_
    !
    IF(iOx/=0) &
    & WRITE(fActiv,'(A2,A1,A5,A1)',ADVANCE='NO') "pE",T_,"RedOx",T_
    !
    !--- write component names
    DO iPr=1,SIZE(vCpn) 
      WRITE(fActiv,'(A,A1)',ADVANCE='NO') TRIM(vEle(vCpn(iPr)%iEle)%NamEl),T_
    ENDDO
    !
    !--- write species names (all aqu' + mobile non'aqu')
    DO I=1,SIZE(vSpc)
      IF(vLCi(I).OR.vLAs(I).OR.vLAx(I).OR.vLMx(I)) &
      & WRITE(fActiv,'(A,A1)',ADVANCE='NO') TRIM(vSpc(I)%NamSp),T_
    ENDDO
    WRITE(fActiv,*)
    !
  ENDIF
  !-------------------------------------------------/ACTIVITIES.HEADER--
  !
  !CALL OutStrVec(fMoles,vCpn(1:nCp)%Mole,I=WrCount,S="TotEle=",C="F")
  !_____________________________________Output
  !CALL OutStrVec(fMoles,vSpc(1:nAq)%Dat%Mole,I=WrCount,S="MolNum=",C="F")
  !
  ZBal=DOT_PRODUCT(vSpc(1:nAq)%Dat%Mole,vSpc(1:nAq)%Z)
  !
  CALL Solmodel_CalcMolal(vSpc,isW,vMolal,IonStr)
  !
  !------------------------------------------------left part of tables--
  !
  !----------------------------------------------ACTIV'COEFF'.RESULTS--
  WRITE(fGamma,'(I3,A1)',      ADVANCE='NO') WrCount,T_
  WRITE(fGamma,'(3(G15.3,A1))',ADVANCE='NO') TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  WRITE(fGamma,'(3(G15.6,A1))',ADVANCE='NO') IonStr,T_,pH_,T_,ZBal,T_
  !
  !-------------------------------------------------ACTIVITIES.RESULTS--
  IF(PRESENT(WrTitl)) WRITE(fActiv,'(A,A1)',ADVANCE='NO') TRIM(WrTitl),T_
  WRITE(fActiv,'(I3,A1)',      ADVANCE='NO') WrCount,T_
  WRITE(fActiv,'(3(G15.3,A1))',ADVANCE='NO') TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  WRITE(fActiv,'(3(G15.6,A1))',ADVANCE='NO') IonStr,T_,pH_,T_,ZBal,T_
  !
  !-------------------------------------------------POTENTIALS.RESULTS--
  IF(PRESENT(WrTitl)) WRITE(fActiv,'(A,A1)',ADVANCE='NO') TRIM(WrTitl),T_
  WRITE(fPoten,'(I3,A1)',      ADVANCE='NO') WrCount,T_
  WRITE(fPoten,'(3(G15.3,A1))',ADVANCE='NO') TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  WRITE(fPoten,'(3(G15.6,A1))',ADVANCE='NO') IonStr,T_,pH_,T_,ZBal,T_
  !
  !-------------------------------------------------MOLALITIES.RESULTS--
  WRITE(fMolal,'(I3,A1)',      ADVANCE='NO') WrCount,T_
  WRITE(fMolal,'(3(G15.3,A1))',ADVANCE='NO') TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  WRITE(fMolal,'(3(G15.6,A1))',ADVANCE='NO') IonStr,T_,pH_,T_,ZBal,T_
  !
  IF(iOx/=0) THEN
    ZRdx= DOT_PRODUCT(vSpc(1:nAq)%Dat%Mole,tStoikio(iOx,1:nAq))
    WRITE(fActiv,'(2(G15.6,A1))',ADVANCE='NO') pE_,T_,ZRdx,T_
    WRITE(fPoten,'(2(G15.6,A1))',ADVANCE='NO') pE_,T_,ZRdx,T_
    WRITE(fGamma,'(2(G15.6,A1))',ADVANCE='NO') pE_,T_,ZRdx,T_
    WRITE(fMolal,'(2(G15.6,A1))',ADVANCE='NO') pE_,T_,ZRdx,T_
  ENDIF
  !
  DO iPr=1,SIZE(vCpn)
    !
    xTotF= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole)
    !
    WRITE(fActiv,'(G15.8,A1)',ADVANCE='NO') xTotF,T_
    WRITE(fPoten,'(G15.8,A1)',ADVANCE='NO') xTotF,T_
    WRITE(fGamma,'(G15.8,A1)',ADVANCE='NO') xTotF,T_
    WRITE(fMolal,'(G15.8,A1)',ADVANCE='NO') xTotF,T_
    !
  ENDDO
  !-----------------------------------------------/left part of tables--
  !
  !!RdxBal=Zero
  !!IF(iOx/=0) THEN !redox balance
  !!  RdxBal=                    dot_product(tAlfPr(iOx,1:nCi),        vMolF(1:        nCi))
  !!  IF(nAs>0) RdxBal= RdxBal + dot_product(tAlfAs(iOx,1:    nAs),    vMolF(nCi+1:    nCi+nAs))
  !!  IF(nAx>0) RdxBal= RdxBal + dot_product(tAlfPr(iOx,nCi+1:nCi+nAx),vMolF(nCi+nAs+1:nCi+nAs+nAx))
  !!ENDIF
  !
  !
  !!! DO I=1,SIZE(vSpc)
  !!!   IF(vLCi(I).OR.vLAs(I).OR.vLAx(I).OR.vLMx(I)) &
  !!!   & WRITE(fActiv,'(G15.6,A1)',ADVANCE='NO') vSpc(I)%LnAct/Ln10,T_
  !!! ENDDO
  !!! WRITE(fActiv,*)
  !
  !-----------------------------------------------RIGHT PART OF TABLES--
  !-------------------------------------------------ACTIVITIES.RESULTS--
  CALL OutStrVec( &
  & fOut= fActiv,&
  & Vec=  vSpc(:)%Dat%LAct/Ln10, &
  & vL=   vLCi.OR.vLAs.OR.vLAx.OR.vLMx)
  !
  !-------------------------------------------------POTENTIALS.RESULTS--
  vPot(:)= vSpc(:)%Dat%LAct +vSpc(:)%G0rt !-> "reduced"potential, dimensionless
  vPot(:)= vPot(:) *R_JK *TdgK            !-> potential in Joule
  ! vPot(:)= vPot(:) - Tref*vS0Ele(:)       !-> in Berman-Brown convention
  !
  CALL OutStrVec( &
  & fOut= fPoten,&
  & Vec=  vPot(:), &
  & vL=   vLCi.OR.vLAs.OR.vLAx.OR.vLMx)
  !
  !------------------------------------------------ACTIV'COEFF'.RESULTS--
  CALL OutStrVec( &
  & fOut= fGamma,&
  & Vec=  vSpc(:)%Dat%LGam/Ln10, &
  & vL=   vLCi.OR.vLAs.OR.vLAx)
  !
  !-------------------------------------------------MOLALITIES.RESULTS--
  CALL OutStrVec( &
  & fOut= fMolal,&
  & Vec=  vMolal(:))
  !----------------------------------------------/RIGHT PART OF TABLES--
  !
  !---------------------------------------------------------------------
  !------------------------------------------------MOLE.NUMBERS.HEADER--
  IF(fMoles==0) THEN
    !
    !IF(ALLOCATED(tPathResults)) DEALLOCATE(tPathResults)
    !IF(ALLOCATED(vNamFasYes)) DEALLOCATE(vNamFasYes)
    !
    CALL GetUnit(fMoles)
    OPEN(fMoles,FILE=TRIM(DirOut)//"_moles.restab")
    !
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_moles.restab", &
    & "PATH: mole nrs of elements, influid, in system, mole nrs of other phases")
    !
    IF(PRESENT(WrTitl)) WRITE(fMoles,'(A,A1)',ADVANCE='NO') ".Title",T_
    WRITE(fMoles,'(4(A,A1))',ADVANCE='NO') "COUNT",T_,"TdgC",T_,"Pbar",T_,"RhoW",T_
    WRITE(fMoles,'(2(A,A1))',ADVANCE='NO') "pH",T_,"pE",T_
    !
    !--------------------------------------------write component names--
    DO iPr=1,SIZE(vCpn)
      WRITE(fMoles,'(A,A1)',ADVANCE='NO') TRIM(vEle(vCpn(iPr)%iEle)%NamEl),T_
    ENDDO
    !
    !---------------------------------------write non'aqu'phases names--
    IF(PRESENT(vYes) .AND. WrCod(1:2)=="EQ") THEN
      !
      !------------------------------------------ALLOCATE BUFFER TABLE--
      !-------will be used to build a result file without zero columns--
      !------used especially for EQUPATH with large number of minerals--
      !
      IF(COUNT(vYes)>0) THEN
        ALLOCATE(vNamFasYes(COUNT(vYes)))
        ALLOCATE(tPathResults(DimPath,2+2*nCp+COUNT(vYes)))
        tPathResults(:,:)= Zero
      ENDIF
      !-----------------------------------------/ALLOCATE BUFFER TABLE--
      !
      DO iPr=1,SIZE(vCpn)
        WRITE(fMoles,'(A,A1)',ADVANCE='NO') &
        & TRIM(vEle(vCpn(iPr)%iEle)%NamEl)//"Tot",T_
      ENDDO
      !
      K=0
      DO I=1,SIZE(vFas)
        IF(vYes(I)) THEN
          WRITE(fMoles,'(A,A1)',ADVANCE='NO') TRIM(vFas(I)%NamFs),T_
          K=K+1
          IF(ALLOCATED(vNamFasYes)) vNamFasYes(K)= TRIM(vFas(I)%NamFs)
        ENDIF
      ENDDO
    ENDIF
    !
    WRITE(fMoles,*)
    !
  ENDIF
  !-----------------------------------------------/MOLE.NUMBERS.HEADER--
  !
  !-----------------------------------------------MOLE.NUMBERS.RESULTS--
  !WRITE(fMoles,'(A4,A1,I3,A1)',ADVANCE='NO') "ELEM",T_,WrCount,T_
  IF(PRESENT(WrTitl)) WRITE(fMoles,'(A,A1)',ADVANCE='NO') TRIM(WrTitl),T_
  !
  WRITE(fMoles,'(I4,A1,3(G15.3,A1))',ADVANCE='NO')  &
  & WrCount,T_,TdgK -T_CK,T_,Pbar,T_,SolModel%Dat%Rho,T_
  !
  WRITE(fMoles,'(2(G15.3,A1))',ADVANCE='NO')  pH_,T_,pE_,T_
  !
  IF(ALLOCATED(tPathResults)) THEN !-------------WRITE IN BUFFER TABLE--
    tPathResults(WrCount,1)= TdgK -T_CK
    tPathResults(WrCount,2)= Pbar
  ENDIF
  !
  !---------------------------------------write total element in fluid--
  DO iPr=1,SIZE(vCpn)
    X= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole) !-> mole nrs in fluid
    WRITE(fMoles,'(G15.8,A1)',ADVANCE='NO') X,T_
    IF(ALLOCATED(tPathResults)) &
    & tPathResults(WrCount,iPr+2)= X !-----------WRITE IN BUFFER TABLE--
  ENDDO
  !
  !-------------------------------write total element in fluid + min's--
  IF(PRESENT(vYes) .AND. WrCod(1:2)=="EQ") THEN
    !
    DO iPr=1,SIZE(vCpn)
      !
      X= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole) !-> mole nrs in fluid
      !
      DO I=1,SIZE(vFas)
        IF(vYes(I) .AND. vFas(I)%Mole>Zero) &
        & X= X + tAlfFs(iPr,I) *vFas(I)%Mole !-> mole nrs in fluid + min's
      ENDDO
      !
      WRITE(fMoles,'(G15.8,A1)',ADVANCE='NO') X,T_
      !
      IF(ALLOCATED(tPathResults)) &
      & tPathResults(WrCount,iPr+nCp+2)= X !-----WRITE IN BUFFER TABLE--
      !
    ENDDO
    !
    K= nCp*2 +2
    DO I=1,SIZE(vFas) !mineral mole numbers from Equil_N
      !IF(.NOT.vFasYes(I))  vFas(I)%Mole=-999.D0 
      IF(vYes(I)) THEN
        !
        K=K+1
        IF(vFas(I)%Mole>Zero) THEN
          !
          !! nOx= One !*tFormula(iFs,1)
          !! IF(tAlfFs(1,I)/=0) nOx= tAlfFs(1,I) !-> nr of oxygens in formula
          !! WRITE(fMoles,'(G15.6,A1)',ADVANCE='NO') vFas(I)%Mole*nOx,T_ !  &
          !
          !! WRITE(fMoles,'(G15.6,A1)',ADVANCE='NO') 
          !! & vFas(I)%Mole /vCpn(1)%Mole/MWsv,T_
          !!-> mole nrs of phase / kgH2O
          !
          nOx= One !*tFormula(iFs,1)
          !
          IF(OutputInGrams) THEN
            !
            WRITE(fMoles,'(G15.6,A1)',ADVANCE='NO') &
            & vFas(I)%Mole*vFas(I)%WeitKg*1.0D3,T_
            !
            !-- WRITE IN BUFFER TABLE --
            IF(ALLOCATED(tPathResults)) &
            & tPathResults(WrCount,K)= vFas(I)%Mole*vFas(I)%WeitKg*1.0D3
            !
          ELSE
            !
            IF(tAlfFs(1,I)/=0) nOx= tAlfFs(1,I) !-> nr of oxygens in formula
            WRITE(fMoles,'(G15.6,A1)',ADVANCE='NO') vFas(I)%Mole*nOx,T_ !  &
            !
            !-- WRITE IN BUFFER TABLE --
            IF(ALLOCATED(tPathResults)) &
            & tPathResults(WrCount,K)= vFas(I)%Mole*nOx
            !
          END IF
          !
        ELSE
          !
          WRITE(fMoles,'(A,A1)',ADVANCE='NO') "0.0",T_
          IF(ALLOCATED(tPathResults)) &
          & tPathResults(WrCount,K)= Zero
          !
        ENDIF
        !
        !
      ENDIF
    ENDDO
    !
  ENDIF
  !------------------------------/write total element in fluid + min's--
  WRITE(fMoles,*)
  !----------------------------------------------/MOLE.NUMBERS.RESULTS--
  !---------------------------------------------------------------------
  !
  DEALLOCATE(vS0Ele)
  DEALLOCATE(vPot)
  !
  RETURN
END SUBROUTINE Path_Write_Line

SUBROUTINE WriteBuffer
  USE M_IOTools
  USE M_Files,    ONLY: DirOut
  !
  USE M_Global_Vars,ONLY: vFas,vEle
  USE M_System_Vars,ONLY: vCpn
  USE M_Path_Vars,ONLY: tPathResults
  !
  INTEGER:: FF,I,J,K,nYes,nCp
  LOGICAL,ALLOCATABLE:: vIsZero(:)
  !
  nCp= SIZE(vCpn)
  nYes= SIZE(tPathResults,2) -(2+2*nCp)
  ALLOCATE(vIsZero(nYes))
  !
  DO J=1,nYes
    vIsZero(J)= (ABS(MAXVAL(tPathResults(:,2+2*nCp+J)))<1.D-16)
  ENDDO
  !
  !PRINT *,"=====nYes-COUNT(vIsZero)===",nYes-COUNT(vIsZero)
  !RETURN
  !
  CALL GetUnit(FF)
  OPEN(FF,FILE=TRIM(DirOut)//"_moles_clean.restab")
  !
  !-------------------------------------------------------------HEADER--
  WRITE(FF,'(A,A1)',ADVANCE='NO') "COUNT",T_
  WRITE(FF,'(2(A,A1))',ADVANCE='NO') "TdgC",T_,"Pbar",T_
  !~ WRITE(FF,'(3(A,A1))',ADVANCE='NO') "RhoW",T_,"pH",T_,"pE",T_
  !~ !
  !-----------------------------------------------write component names--
  DO I=1,SIZE(vCpn)
    WRITE(FF,'(A,A1)',ADVANCE='NO') &
    & TRIM(vEle(vCpn(I)%iEle)%NamEl)//"Fluid",T_
  ENDDO
  !
  DO I=1,SIZE(vCpn)
    WRITE(FF,'(A,A1)',ADVANCE='NO') &
    & TRIM(vEle(vCpn(I)%iEle)%NamEl)//"Tot",T_
  ENDDO
  !-----------------------------------------write non'aqu'phases names--
  K=0
  DO I=1,nYes
    K=K+1
    IF(.NOT. vIsZero(I)) THEN
      WRITE(FF,'(A,A1)',ADVANCE='NO') TRIM(vNamFasYes(K)),T_
    ENDIF
  ENDDO
  WRITE(FF,*)
  !------------------------------------------------------------/HEADER--
  !
  DO I=1,SIZE(tPathResults,1)
    !
    WRITE(FF,'(I3,A1)',ADVANCE='NO') I,T_
    WRITE(FF,'(G15.8,A1)',ADVANCE='NO') tPathResults(I,1),T_ !! TdgC
    WRITE(FF,'(G15.8,A1)',ADVANCE='NO') tPathResults(I,2),T_ !! Pbar
    !
    DO J=1,2*nCp
      WRITE(FF,'(G15.8,A1)',ADVANCE='NO') tPathResults(I,J+2),T_
    ENDDO
    !
    DO J=1,nYes
      IF(.NOT. vIsZero(J)) &
      & WRITE(FF,'(G15.8,A1)',ADVANCE='NO') tPathResults(I,J+2+2*nCp),T_
    ENDDO
    !
    WRITE(FF,*)
    !
  ENDDO
   !
  DEALLOCATE(vIsZero)
  !
  CLOSE(FF)
  !
END SUBROUTINE WriteBuffer

SUBROUTINE Path_Files_Close
  USE M_Path_Vars,ONLY: tPathResults
  !
  IF(fActiv>0) THEN  ;  WRITE(fActiv,*)  ;  CLOSE(fActiv) ;  fActiv=0  ;  ENDIF
  IF(fPoten>0) THEN  ;  WRITE(fPoten,*)  ;  CLOSE(fPoten) ;  fPoten=0  ;  ENDIF
  IF(fMolal>0) THEN  ;  WRITE(fMolal,*)  ;  CLOSE(fMolal) ;  fMolal=0  ;  ENDIF
  IF(fGamma>0) THEN  ;                      CLOSE(fGamma) ;  fGamma=0  ;  ENDIF
  IF(fQsk>0)   THEN  ;                      CLOSE(fQsk)   ;  fQsk=0    ;  ENDIF
  !
  IF(fMoles>0) THEN
  
    WRITE(fMoles,*)
    CLOSE(fMoles)
    fMoles=0
    
    IF(ALLOCATED(tPathResults)) THEN
      CALL WriteBuffer
      DEALLOCATE(tPathResults)
      DEALLOCATE(vNamFasYes)
    ENDIF
    
  ENDIF
  !
END SUBROUTINE Path_Files_Close

ENDMODULE M_Path_Write

