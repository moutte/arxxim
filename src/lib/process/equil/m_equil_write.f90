MODULE M_Equil_Write
!--
!-- routines writing results of equilibrium or speciation calculations
!--
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,fHtm,T_,Pause_
  IMPLICIT NONE

  PRIVATE

  PUBLIC:: Equil_Write_Detail
  PUBLIC:: Equil_Write_ShoDetail
  PUBLIC:: Equil_Write_EnTete
  PUBLIC:: Equil_Write_LogK

  LOGICAL:: Done_WriteLogK= .false.

CONTAINS

SUBROUTINE Equil_Write_LogK(vCpn,TdgK,Pbar)
!--
!-- write a table of logK of species
!-- calculated with logK-0 for the SYSTEM's primary species
!--
  USE M_T_Component,ONLY: T_Component
  USE M_Dtb_Test,   ONLY: Dtb_Tabulate_ForSystem
  USE M_Basis,      ONLY: Basis_Build
  USE M_Basis_Vars, ONLY: vLCi,vLAx,vOrdPr,tNuSp
  USE M_TPcond_Read,ONLY: TPpath_Read
  !
  USE M_Global_Vars,ONLY: vSpc,vSpcDtb
  USE M_Path_Vars,  ONLY: vTPpath
  !
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  REAL(dp),         INTENT(IN):: TdgK,Pbar
  !
  TYPE(T_Component):: vCpnTmp(SIZE(vCpn))
  !
  IF(Done_WriteLogK) RETURN !-------------------------------------return
  !
  vCpnTmp= vCpn
  CALL Basis_Build(.FALSE.,vCpnTmp)
  !
  CALL TPpath_Read(TdgK,Pbar)
  !
  CALL Dtb_Tabulate_ForSystem( &
  & vTPpath, &
  & vCpnTmp,vSpc,vSpcDtb, &
  & vLCi,vLAx,vOrdPr,tNuSp)
  !
  DEALLOCATE(vTPpath)
  !
  Done_WriteLogK= .true.
  !
  RETURN
ENDSUBROUTINE Equil_Write_LogK

SUBROUTINE Equil_Write_EnTete(iFile,vSpc,vPrm,vBool,Str1,Str2,WriteCR)
!--
!-- write title line: a name (s) followed by species names
!--
  USE M_T_Species,ONLY: T_Species
  !
  INTEGER,        INTENT(IN):: iFile
  TYPE(T_Species),INTENT(IN):: vSpc(:)
  INTEGER,        INTENT(IN):: vPrm(:)
  LOGICAL,        INTENT(IN):: vBool(:)
  !
  CHARACTER(*),   INTENT(IN),OPTIONAL:: Str1
  CHARACTER(*),   INTENT(IN),OPTIONAL:: Str2
  LOGICAL,        INTENT(IN),OPTIONAL:: WriteCR
  !
  INTEGER::I,J
  !
  IF(PRESENT(Str1)) WRITE(iFile,'(A,A1)',ADVANCE='NO') TRIM(Str1),T_
  IF(PRESENT(Str2)) WRITE(iFile,'(A,A1)',ADVANCE='NO') TRIM(Str2),T_

  DO J=1,SIZE(vSpc)
    I=vPrm(J)
    IF(vBool(I)) WRITE(iFile,'(A,A1)',ADVANCE='NO') TRIM(vSpc(I)%NamSp),T_
  ENDDO
  IF(.NOT. PRESENT(WriteCR)) THEN
    WRITE(iFile,*)
  ELSE
    IF(WriteCR) WRITE(iFile,*)
  ENDIF

  RETURN
ENDSUBROUTINE Equil_Write_EnTete

SUBROUTINE Equil_Write_Detail(cSelec,TdgK,Pbar,vCpn)
!--
!-- write detailed results for a single speciation
!--
  USE M_Numeric_Tools !->CalcPermut
  USE M_IoTools !GetUnit
  USE M_Files,       ONLY: DirOut,cTitle,File_Write_Date_Time,Files_Index_Write
  USE M_Dtb_Const,   ONLY: T_CK, R_JK, F_JV
  USE M_Numeric_Const,ONLY: MinExpDP,MaxExpDP,Ln10
  USE M_T_Species,   ONLY: T_Species
  USE M_T_Component, ONLY: T_Component
  USE M_SolModel_Tools,ONLY: Solmodel_CalcMolal,Solmodel_pHpE
  USE M_Basis,       ONLY: Basis_CpnInert,Basis_Change
  !
  USE M_Global_Vars, ONLY: vSpc,vSpcDtb,vEle,nAq
  USE M_Global_Vars, ONLY: vFas,vMixFas,vMixModel
  USE M_Basis_Vars,  ONLY: iH_,iBal,isW,MWsv,iOx,isH_,isO2,tStoikio
  USE M_Basis_Vars,  ONLY: vOrdPr,tAlfFs,tNuFas
  !
  CHARACTER(LEN=3), INTENT(IN):: cSelec
  REAL(dp),         INTENT(IN):: TdgK,Pbar
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  !
  TYPE(T_Species)  :: S
  CHARACTER(LEN=80):: Str1,Str2,Str
  CHARACTER(LEN=80):: strFile
  INTEGER :: nCp,nFs,nSp,nPur
  INTEGER :: iPr,iAq,iFs,fo,J,K,iModel
  REAL(dp):: IonStr
  REAL(dp):: X,Y,Affin,pH_,pE_,Eh_
  REAL(dp):: Z_Plus,Z_Minus
  LOGICAL,PARAMETER:: ReportMajorOnly= .FALSE.
  !
  INTEGER, ALLOCATABLE::vPrm(:),vPrmMs(:)
  REAL(dp),ALLOCATABLE::vMolal(:),vQsK(:)
  !
  TYPE(T_Component),DIMENSION(:),ALLOCATABLE:: vCpnInert
  ! TYPE(T_Component),DIMENSION(:),ALLOCATABLE:: vC
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Write_Detail"
  !
  nCp=SIZE(vCpn)
  nSp=SIZE(vSpc)
  nFs=SIZE(vFas)
  !
  ALLOCATE(vPrm(nAq))
  ALLOCATE(vMolal(nAq))
  !
  Z_Plus= SUM(vSpc(1:nAq)%Dat%Mole*vSpc(1:nAq)%Z, MASK=(vSpc(1:nAq)%Z >0))
  Z_Minus=SUM(vSpc(1:nAq)%Dat%Mole*vSpc(1:nAq)%Z, MASK=(vSpc(1:nAq)%Z <0))
  !
  CALL Solmodel_CalcMolal(vSpc,isW,vMolal,IonStr)
  CALL Solmodel_pHpE(isW,iOx,isH_,isO2,vSpc,pH_,pE_)
  !
  SELECT CASE(TRIM(cSelec))
  CASE("SPC")  ;  strFile="_specia.res"
  CASE("MIX")  ;  strFile="_specia_mix.res"
  CASE("BOX")  ;  strFile="_specia_box.res"
  CASE("INJ")  ;  strFile="_specia_inj.res"
  CASE("INI")  ;  strFile="_specia_ini.res"
  CASE("DYN")  ;  strFile="_specia_end.res"
  CASE("EQ1","EQ2","EQM")  ; strFile="_equil.res"
  ENDSELECT
  !
  CALL GetUNit(fo)
  OPEN(fo,FILE=TRIM(DirOut)//TRIM(strFile))

  CALL Files_Index_Write(fHtm,&
  & TRIM(DirOut)//TRIM(strFile),&
  & "SPC/EQU: result of speciation or equilibrium")
  !
  CALL File_Write_Date_Time(fo)
  !
  WRITE(fo,'(A)') "!."//TRIM(cTitle)
  !
  SELECT CASE(TRIM(cSelec))
    CASE("EQ1","EQ2","EQM")
      WRITE(fo,'(A)') "!.fluid in equilibrium with other phases"
    CASE("SPC")
      WRITE(fo,'(A)') "!.fluid speciation (before interaction)"
    CASE("DYN")
      WRITE(fo,'(A)') "!.fluid speciation (after interaction)"
    CASE("INI")
      WRITE(fo,'(A)') "!.fluid speciation (begin interaction)"
    CASE("MIX")
      WRITE(fo,'(A)') "!.fluid speciation of mixing fluid"
  ENDSELECT
  !
  !WRITE(fo,'(A)') "INPUT="//TRIM(NamFInn)
  !-------------------------------------------------write Charge Balance
  X= Zero
  IF(Z_Plus - Z_Minus > EPSILON(X)) X= (Z_Plus + Z_Minus)/(Z_Plus - Z_Minus)
  WRITE(fo,'(/,A)') &
  & "Charge balance, sum(+), sum(-), delta/sum"
  WRITE(fo,'(3G15.6,/)') Z_Plus,Z_Minus,X
  !
  IF(iBal==0 .AND. vCpn(iH_)%Statut/="INERT") &
  & WRITE(fo,'(A,/,A,/,A)') &
  & "!!!", &
  & "!CAVEAT! NO Element For Electron Balance -> ELECTRONEUTRALITY NOT ENFORCED", &
  & "!!!"
  !--/
  
  !---------------------------------------------------GLOBAL COMPOSITION
  ALLOCATE(vCpnInert(1:nCp))
  CALL Basis_CpnInert(isW,tStoikio,vCpn,vCpnInert)
  !
  WRITE(fo,'(A,/,A,/,A)') &
  & "!-----------------------------------------------------------------------", &
  & "!--fluid composition, mole nr / kgH2O-----------------------------------", &
  & "!--can be used as SYSTEM for a new all-inert run------------------------"
  WRITE(fo,'(2X,A,G15.6)') "TdgC  ",TdgK-T_CK
  WRITE(fo,'(2X,A,G15.6)') "Pbar  ",Pbar
  DO iPr=1,SIZE(vCpn)
    WRITE(fo,'(2X,3(A,1X),G15.6)') &
    & vEle(vCpnInert(iPr)%iEle)%NamEl, &
    & "INERT", &
    & vSpc(vCpnInert(iPr)%iSpc)%NamSp,&
    & vCpnInert(iPr)%Mole /MWsv /vCpnInert(isW)%Mole
  ENDDO
  WRITE(fo,'(A,/)') &
  & "!-----------------------------------------------------------------------"
  !
  DEALLOCATE(vCpnInert)
  !
  IF(cSelec(1:2)=="EQ") THEN

    WRITE(fo,'(A,/,A)') "!","!total composition"
    !WRITE(fo,'(A)') "!can be used as input for fluid mixing"
    DO iPr=1,SIZE(vCpn)
      X= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole) &
      +  SUM(tAlfFs(iPr,1:nFs) *vFas(1:nFs)%Mole)
      WRITE(fo,'(A,2(A1,G24.17))') vEle(vCpn(iPr)%iEle)%NamEl,T_,X !,T_,Y
    ENDDO

  ENDIF
  !--------------------------------------------------/GLOBAL COMPOSITION
  !
  !--------------------------------------------------ELEMENTS,COMPONENTS
  WRITE(fo,'(A)') &
  & "!-----------------------------------------------------------------------"
  WRITE(fo,'(A)') "!element, status, mole balance, mole number, mg/kg"
  !
  DO iPr=1,SIZE(vCpn)
    !
    X= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole) !-> abundance in fluid
    Y= X *vEle(vCpn(iPr)%iEle)%WeitKg *1.D6 !-> kg to mg

    !---COMMENT--
    !  value of vCpn(:)%Mole has not been updated !!
    !  seems it is the input value ....
    !---/

    WRITE(fo,'(I3,A1,2(A,A1),3(G15.8,A1))',ADVANCE="NO") &
    & iPr,                       T_, &
    & vEle(vCpn(iPr)%iEle)%NamEl,T_, &
    & vCpn(iPr)%Statut,          T_, &
    & vCpn(iPr)%Mole,            T_, &
    & X,T_,Y,T_
    !
    !IF(vCpn(iPr)%Statut=="MOBILE") & !
    IF(    TRIM(vCpn(iPr)%Statut)=="MOBILE" &
    & .OR. TRIM(vCpn(iPr)%Statut)=="BUFFER") &
    & WRITE(fo,'(G15.6,A1,A,A1)',ADVANCE="NO") &
    & EXP(vSpc(vCpn(iPr)%iSpc)%Dat%LAct),T_,vSpc(vCpn(iPr)%iSpc)%NamSp,T_
    !
    WRITE(fo,*)
    !
  ENDDO
  !
  !-------------------------------------------------/ELEMENTS,COMPONENTS
  
  !--------------------------------------------------------------pH, etc
  WRITE(fo,'(A)') &
  & "!-----------------------------------------------------------------------"
  WRITE(fo,'(A,G12.3)') "pH=         ",pH_ !vSpc(isH_)%LnAct/Ln10
  IF(iOx>0) THEN
    Eh_= pE_ *TdgK *Ln10 *R_JK /F_JV
    WRITE(fo,'(A,G12.3)') "pE=         ",pE_
    WRITE(fo,'(A,G12.3)') "Eh(Volts)=  ",Eh_
  ENDIF
  WRITE(fo,'(A12,G12.3)')            "IonStrength=",IonStr
  !------------------------------------------------------------/ pH, etc
  
  !---------------------------------------------------EQUILIBRIUM.PHASES
  IF(cSelec(1:2)=="EQ") THEN
    !
    WRITE(fo,'(A)') &
    & "!-----------------------------------------------------------------------"

    WRITE(fo,'(A)') "Equilibrium Species (result of Equil_n)"

    DO iFs=1,SIZE(vFas)

      IF(vFas(iFs)%Mole > Zero) THEN

        WRITE(fo,'(2A,G15.6)', ADVANCE="NO") &
        & vFas(iFs)%NamFs," MOLE=",vFas(iFs)%Mole

        IF(vFas(iFs)%iMix>0) THEN
          WRITE(fo,'(A)',ADVANCE="NO") " X(:)="
          !print *, "vFas(iFs)%iMix",vFas(iFs)%iMix
          !pause
          iModel= vMixFas(vFas(iFs)%iMix)%iModel
          !print *, "iModel, nPole=",iModel,vMixModel(iModel)%nPole
          !pause
          DO J=1,vMixModel(iModel)%nPole
            WRITE(fo,'(G15.6,1X)',ADVANCE="NO") &
            & vMixFas(vFas(iFs)%iMix)%vXPole(J)
          ENDDO
        ENDIF
        WRITE(fo,*)

      ENDIF

    ENDDO

    WRITE(fo,'(A)') &
    & "!-----------------------------------------------------------------------"

    !-------------------------------------------------------------header
    WRITE(fo,'(A7,A1)', ADVANCE="NO") "Element",T_
    WRITE(fo,'(A15,A1)',ADVANCE="NO") "Fluid          ",T_
    DO iFs=1,nFs
      IF(vFas(iFs)%Mole > Zero) &
      WRITE(fo,'(A15,A1)',ADVANCE="NO") vFas(iFs)%NamFs,T_
    ENDDO
    WRITE(fo,*)
    !------------------------------------------------------------------/

    !--------------results: distribution of element amounts among phases
    DO iPr=1,SIZE(vCpn)
      WRITE(fo,'(A7,A1)',ADVANCE="NO") vEle(vCpn(iPr)%iEle)%NamEl,T_
      X= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole)
      WRITE(fo,'(G15.6,A1)',ADVANCE="NO") X,T_
      DO iFs=1,nFs
        IF(vFas(iFs)%Mole > Zero .AND. vFas(iFs)%iSpc/=0) THEN
          !!!toDO!!! PRINT *,vFas(iFs)%NamFs
          WRITE(fo,'(G15.6,A1)',ADVANCE="NO") &
          & tStoikio(iPr,vFas(iFs)%iSpc)*vFas(iFs)%Mole,T_
        ENDIF
      ENDDO
      WRITE(fo,*)
    ENDDO
    !------------------------------------------------------------------/

    WRITE(fo,'(A)') &
    & "!-----------------------------------------------------------------------"

    WRITE(fo,'(A)') "Fluid Composition"
    DO iPr=1,SIZE(vCpn)
      X=SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole)
      WRITE(fo,'(I3,A1,2(A,A1),G15.6,A1)',ADVANCE="NO") &
      & iPr,T_, vEle(vCpn(iPr)%iEle)%NamEl,T_, vCpn(iPr)%Statut,T_,X,T_
      !IF(TRIM(vCpn(iPr)%Statut)=="MOBILE") &!
      IF(TRIM(vCpn(iPr)%Statut)=="MOBILE" .OR. &
      &  TRIM(vCpn(iPr)%Statut)=="BUFFER") &
      & WRITE(fo,'(G15.6,A1,A,A1)',ADVANCE="NO") &
      & EXP(vSpc(vCpn(iPr)%iSpc)%Dat%LAct),T_,vSpc(vCpn(iPr)%iSpc)%NamSp,T_
      WRITE(fo,*)
    ENDDO

    WRITE(fo,'(A)') &
    & "!-----------------------------------------------------------------------"

  ENDIF !(cSelec=="EQ*")
  !--------------------------------------------------/EQUILIBRIUM.PHASES
  !
  !--- order species --
  CALL CalcPermut(vSpc(1:nAq)%Dat%Mole,vPrm)
  !-> order species by increasing abundance
  !to apply permutation vPermut to array Arr: Arr=Arr(vPermut)
  !---/
  !
  !----------------------------------------------------------AQU.SPECIES
  WRITE(fo,'(A)') &
  & "!-----------------------------------------------------------------------"
  WRITE(fo,'(A)') "All Species, in order of decreasing Mole number"
  !
  WRITE(fo,'(A15,A1)',   ADVANCE='NO') "________SPECIES",T_
  WRITE(fo,'(2(A12,A1))',ADVANCE='NO') "_MOLE_NUMBER",T_,"_MOLALITY___",T_
  WRITE(fo,'(2(A12,A1))',ADVANCE='NO') "_LOG(GAMMA)_",T_,"_GAMMA______",T_
  WRITE(fo,'(2(A12,A1))',ADVANCE='NO') "LOG(ACTIVIT)",T_,"_ACTIVITY___",T_
  WRITE(fo,*)
  !
  !-----------------------------------------ordered in decreasing amount
  !DO iPr=1,nCp

    DO iAq=nAq,1,-1

      S= vSpc(vPrm(iAq))
      X= S%Dat%Mole /vSpc(isW)%Dat%Mole /vSpc(isW)%WeitKg
      !-> X is molality

      WRITE(fo,'(A15,A1)',ADVANCE='NO')  TRIM(S%NamSp),T_

      WRITE(fo,'(5(G12.5,A1))',ADVANCE='NO') &
      & S%Dat%Mole,      T_, &
      & X,               T_, &
      & S%Dat%LGam/Ln10, T_, &
      & EXP(S%Dat%LGam), T_, &
      & S%Dat%LAct/Ln10, T_

      IF(S%Dat%LAct>MinExpDP .AND. S%Dat%LAct<MaxExpDP) THEN
        WRITE(fo,'(G12.5,A1)') EXP(S%Dat%LAct),T_
      ELSE
        WRITE(fo,*)
      ENDIF

    ENDDO

  !ENDDO
  WRITE(fo,*)
  !---------------------------------------------------------/AQU.SPECIES

  !-----------------------------------------------------RELATIVE AMOUNTS
  WRITE(fo,'(/,2A,/)') &
  & "for each Element, Species in order of decreasing Mole number", &
  & "Relative Amounts in Permil of Total Amount"
  WRITE(fo,*)
  !
  WRITE(fo,'(A15,A1,A7,A1)',ADVANCE='NO')  &
  & "________ELEMENT",T_,"_STOKIO",T_
  WRITE(fo,'(7(A15,A1))',   ADVANCE='NO')  &
  & "________SPECIES",T_, &
  & "RELATIVE_AMOUNT",T_,"_______MOLALITY",T_,"___MOLE_NUMBER_",T_, &
  & "____ACTIV_COEFF",T_,"__LOG(ACTIVITY)",T_,"_______ACTIVITY",T_
  WRITE(fo,*)
  !
  DO iPr=2,nCp
    !
    WRITE(fo,'(A15,A1,A7,A1,A15,A1,F15.3)') &
    & vEle(vCpn(iPr)%iEle)%NamEl,T_, &
    & "___",                    T_, &
    & "        TOT=",           T_, &
    & vCpn(iPr)%Mole
    !
    DO iAq=nAq,1,-1

      IF(tStoikio(iPr,vPrm(iAq))/=0) THEN
        S= vSpc(vPrm(iAq))
        !
        X= 1000.0D0 *ABS(tStoikio(iPr,vPrm(iAq))) *S%Dat%Mole &
        &  /SUM(tStoikio(iPr,:)*vSpc(:)%Dat%Mole)
        !
        WRITE(fo,'(A15,A1,F7.2,A1)',ADVANCE='NO') &
        & vEle(vCpn(iPr)%iEle)%NamEl, T_, &
        & tStoikio(iPr,vPrm(iAq)),   T_
        WRITE(fo,'(A15,A1)',     ADVANCE='NO') TRIM(S%NamSp),T_
        WRITE(fo,'(5(G15.6,A1))',ADVANCE='NO') &
        & X,                T_, &
        & vMolal(vPrm(iAq)),T_, &
        & S%Dat%Mole,       T_, &
        & EXP(S%Dat%LGam),  T_, &
        & S%Dat%LAct/Ln10
        !
        IF(S%Dat%LAct>MinExpDP .AND. S%Dat%LAct<MaxExpDP) THEN
          WRITE(fo,'(G15.6,A1)') EXP(S%Dat%LAct),T_
        ELSE
          WRITE(fo,*)
        ENDIF
        !
      ENDIF

    ENDDO
    !
    WRITE(fo,'(A)') &
    & "!-----------------------------------------------------------------------"
    !
  ENDDO
  !----------------------------------------------------/RELATIVE AMOUNTS

  !--------------------------------------------RELATIVE AMOUNTS (MAJORS)
  IF(ReportMajorOnly) THEN
    !
    WRITE(fo,'(/,A,/)') & !the same, ONLY major species ....
    & "for each Element, Species in order of decreasing Mole number, ONLY Major"
    !
    DO iPr=2,nCp

      WRITE(fo,'(A12,A1,A3,A1,A12,A1,F15.3)') &
      & vEle(vCpn(iPr)%iEle)%NamEl,T_, &
      & "___",                     T_, &
      & "        TOT=",            T_, &
      & vCpn(iPr)%Mole

      DO iAq=nAq,1,-1

        IF(tStoikio(iPr,vPrm(iAq))/=0) THEN
          S=vSpc(vPrm(iAq))
          X=1.0D3 *ABS(tStoikio(iPr,vPrm(iAq))) &
          &       *S%Dat%Mole &
          &       /SUM(tStoikio(iPr,:)*vSpc(:)%Dat%Mole)
          IF(X>=10.0D0) THEN
            WRITE(fo,'(A12,A1,F7.2,A1)',    ADVANCE='NO') &
            & vEle(vCpn(iPr)%iEle)%NamEl, T_, &
            & tStoikio(iPr,vPrm(iAq)),   T_
            WRITE(fo,'(A12,A1)',            ADVANCE='NO') &
            & TRIM(S%NamSp),T_
            WRITE(fo,'(F15.3, A1,F15.4,A1)',ADVANCE='NO') &
            & X,                T_, &
            & vMolal(vPrm(iAq)),T_
            WRITE(fo,'(F15.8, A1,G15.4,A1)',ADVANCE='NO') &
            & S%Dat%Mole, T_, &
            & S%Dat%Mole, T_
            WRITE(fo,'(G15.4,A1,G15.4,A1,G15.4)') &
            & EXP(S%Dat%LGam),T_, &
            & EXP(S%Dat%LAct),T_, &
            & S%Dat%LAct/Ln10
          ENDIF
        ENDIF

      ENDDO

    ENDDO
    !
  ENDIF
  !-------------------------------------------/RELATIVE AMOUNTS (MAJORS)
  !
  !-------------------------------------------AFFINITIES PHASES vs FLUID
  IF(nFs>0) THEN
    !
    nPur= COUNT(vFas(:)%iSpc>0)
    !
    ALLOCATE(vQsK(1:nPur))
    ALLOCATE(vPrmMs(1:nPur))
    !
    WRITE(fo,'(/,A,/)') &
    & "logQsK, all non-aqueous phases, Minerals & Gases"
    !
    !--- calculate logQsk -> store in vQsk
    DO iFs=1,nPur
      Affin= vFas(iFs)%Grt &
      &    - DOT_PRODUCT( tNuFas(iFs,1:nCp),  &
      &        vSpc(vOrdPr(1:nCp))%Dat%LAct + vSpc(vOrdPr(1:nCp))%G0rt )
      vQsK(iFs)= - Affin /Ln10
    ENDDO
    !
    !-------------sort minerals in order of increasing saturation degree
    CALL CalcPermut(vQsk(:),vPrmMs(:))
    !
    !---------------------write sorted (in both F12.3 and G15.8 formats)
    DO iFs=nPur,1,-1
      J=vPrmMs(iFs)
      WRITE(fo,'(A24,A1,F12.3,A1,G15.8,A1)',ADVANCE='NO') &
      & TRIM(vFas(J)%NamFs),T_,&
      & vQsk(J),T_,&
      & vQsk(J),T_
      WRITE(fo,*)
    ENDDO
    !---/
    !
    WRITE(fo,'(/,2(A,/))') &
    & "logQsK, all non-aqueous phases, Minerals & Gases", &
    & "scaled to stoichiometric numbers of moles of prim'species in reaction"
    !--- compute logQsk -> store in vQsk --
    DO iFs=1,nPur
      Affin= vFas(iFs)%Grt &
      &    - DOT_PRODUCT( tNuFas(iFs,1:nCp), &
      &                   vSpc(vOrdPr(1:nCp))%Dat%LAct &
      &                   +vSpc(vOrdPr(1:nCp))%G0rt )
      vQsK(iFs)= -Affin /Ln10 /SUM(ABS(tNuFas(iFs,1:nCp)))
    ENDDO
    !---/
    !
    !---------------sort phases in order of increasing saturation degree
    CALL CalcPermut(vQsk(:),vPrmMs(:))
    !
    !---------------------write sorted (in both F12.3 and G15.8 formats)
    DO iFs=nPur,1,-1
      J= vPrmMs(iFs)
      !! IF(vFas(J)%iSpc/=0) THEN
      Str1="_"
      Str2="_"
      IF(vFas(J)%iSpc/=0) THEN
        K= vFas(J)%iSpc
        IF(vSpc(K)%iDtb>0) Str1= TRIM(vSpcDtb(vSpc(K)%iDtb)%DtbTrace)
        IF(vSpc(K)%iDtb>0) Str2= TRIM(vSpc(K)%Formula)
      ENDIF
      WRITE(fo,'(A24,A1,F12.3,A1,G15.8,A1,A7,A1,A,A1)',ADVANCE='NO') &
      & TRIM(vFas(J)%NamFs),T_, &
      & vQsk(J),            T_, &
      & vQsk(J),            T_, &
      & TRIM(Str1),         T_, &
      & TRIM(Str2),         T_
      WRITE(fo,*)
    ENDDO
    !---/
    !
    DEALLOCATE(vQsK)
    DEALLOCATE(vPrmMs)
    !
  ENDIF !IF(nFs>0)
  !------------------------------------------/AFFINITIES PHASES vs FLUID
  !
  DEALLOCATE(vPrm)
  DEALLOCATE(vMolal)
  !
  CLOSE(fo)
  !
  IF(iDebug>0) THEN
    PRINT '(/,A,/)',&
    & "Speciation detail saved in file "//TRIM(DirOut)//TRIM(strFile)
    IF(iDebug>1) CALL Pause_
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Write_Detail"
  !
ENDSUBROUTINE Equil_Write_Detail

SUBROUTINE Equil_Write_ShoDetail(cSelec,TdgK,Pbar,vCpn)
!--
!-- write detailed results for a single speciation
!--
  USE M_Numeric_Tools !->CalcPermut
  USE M_IoTools !GetUnit
  USE M_Files,       ONLY: DirOut,cTitle
  !
  USE M_Dtb_Const,   ONLY: T_CK
  USE M_Numeric_Const,ONLY: MinExpDP,MaxExpDP,Ln10
  !
  USE M_SolModel_Tools,ONLY: Solmodel_CalcMolal,Solmodel_pHpE
  !
  USE M_T_Component, ONLY: T_Component
  USE M_T_Species,   ONLY: T_Species
  !
  USE M_Global_Vars, ONLY: vSpc,vEle,vFas,nAq
  USE M_Basis_Vars,  ONLY: iH_,isW,iOx,isH_,isO2,MWsv
  USE M_Basis_Vars,  ONLY: vOrdPr,tNuFas,tAlfFs,tStoikio
  !
  CHARACTER(LEN=3), INTENT(IN):: cSelec
  REAL(dp),         INTENT(IN):: TdgK,Pbar
  TYPE(T_Component),INTENT(IN):: vCpn(:)
  !
  TYPE(T_Species):: S
  INTEGER :: nCp,nFs,nSp
  INTEGER :: iPr,iAq,iMs !,J
  REAL(dp):: IonStr
  REAL(dp):: X,pH_,pE_ !,Y,Affin
  REAL(dp):: Z_Plus,Z_Minus
  INTEGER :: vPrm(nAq)
  REAL(dp):: vMolal(nAq)
  !
  INTEGER,PARAMETER:: FF= 6 != write on screen

  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Write_ShoDetail"
  !
  nCp= SIZE(vCpn)
  nSp= SIZE(vSpc)
  nFs= SIZE(vFas)
  nFs= COUNT(vFas(:)%iSpc>0)
  !
  Z_Plus= SUM(vSpc(1:nAq)%Dat%Mole*vSpc(1:nAq)%Z, MASK=(vSpc(1:nAq)%Z >0))
  Z_Minus=SUM(vSpc(1:nAq)%Dat%Mole*vSpc(1:nAq)%Z, MASK=(vSpc(1:nAq)%Z <0))
  !
  CALL Solmodel_CalcMolal(vSpc,isW,vMolal,IonStr)
  CALL Solmodel_pHpE(isW,iOx,isH_,isO2,vSpc,pH_,pE_)
  !
  WRITE(FF,'(A)') "----------------------------------------------------"
  WRITE(FF,'(A)') "----------------------------------------- results --"
  WRITE(FF,'(A)') "----------------------------------------------------"
  WRITE(FF,'(2(A,G15.6,A1))') "TdgC=",TdgK-T_CK,T_,"Pbar=",Pbar,T_
  !
  !---------------------------------------------------global composition
  WRITE(FF,'(A)') "< --------------------- global composition (mole) --"
  DO iPr=1,SIZE(vCpn)
    X=  SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole) & !-> abundance in fluid
    & + SUM(tAlfFs(iPr,1:nFs) *vFas(1:nFs)%Mole)
    WRITE(FF,'(A,A1,G15.6)') vEle(vCpn(iPr)%iEle)%NamEl,T_,X
  ENDDO
  WRITE(FF,'(A)') "</--"
  WRITE(FF,*)

  !----------------------------------------------------fluid composition
  WRITE(FF,'(A)') "< ---------------- fluid composition (mole|activ) --"
  DO iPr=1,SIZE(vCpn)
    !
    WRITE(FF,'(I3,1X,2(A,A1))',ADVANCE="NO") &
    & iPr,vEle(vCpn(iPr)%iEle)%NamEl,T_,vCpn(iPr)%Statut,T_
    !
    IF(TRIM(vCpn(iPr)%Statut)=="MOBILE" .OR. &
    &  TRIM(vCpn(iPr)%Statut)=="BUFFER") THEN

      !------for MOBILE/BUFFER comp't print activity and related species
      WRITE(FF,'(G15.6,A1,A)') &
      & EXP(vSpc(vCpn(iPr)%iSpc)%Dat%LAct), T_, &
      & vSpc(vCpn(iPr)%iSpc)%NamSp

    ELSE

      !---------------------for INERT/BALANCE species, print mole number
      X= SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole)
      WRITE(FF,'(G15.6,A1,G15.6)') vCpn(iPr)%Mole,T_,X

    ENDIF
    !
  ENDDO
  WRITE(FF,'(A)') "</--"
  WRITE(FF,*)

  !---------------------------------------------------------------------
  WRITE(FF,'(A)') "< ----------------------- neutrality, pH, pE, ... --"
  !
  WRITE(FF,'(3(A,G15.6,/))') &
  & "  Sum(+)=     ", Z_Plus,    &
  & "  Sum(-)=     ", Z_Minus,   &
  & "  Balance=    ", ABS(Z_Plus + Z_Minus)
  !
  WRITE(FF,'(A,G12.3)')           "  pH=         ",pH_ !vSpc(isH_)%LnAct/Ln10
  IF(iOx>0) WRITE(FF,'(A,G12.3)') "  pE=         ",pE_
  WRITE(FF,'(A,G12.3)')           "  IonStrength=",IonStr
  !
  WRITE(FF,'(A)') "</-- " !neutrality, pH, pE, ... --"
  WRITE(FF,*)

  !-------------------------------------------------for equilibrium runs
  IF(cSelec(1:2)=="EQ") THEN

    WRITE(FF,'(A)') "< ------------------------- equilibrium results --"
    !
    WRITE(FF,'(A)') "------------------------------- minerals, gases --"
    DO iMs=1,SIZE(vFas)
      IF(vFas(iMs)%Mole > Zero) &
      WRITE(FF,'(G15.6,2A)') vFas(iMs)%Mole,"= ",TRIM(vFas(iMs)%NamFs)
    ENDDO
    WRITE(FF,'(A)') "--"
    !
    WRITE(FF,'(A)') "----------------------------------------- fluid --"
    DO iPr=1,SIZE(vCpn)
      !
      X=SUM(tStoikio(iPr,1:nAq)*vSpc(1:nAq)%Dat%Mole)
      WRITE(FF,'(I3,A1,2(A,A1),G15.8,A1)',ADVANCE="NO") &
      & iPr,                                        T_, &
      & vEle(vCpn(iPr)%iEle)%NamEl,                 T_, &
      & vCpn(iPr)%Statut,                           T_, &
      & X,                                          T_
      !IF(TRIM(vCpn(iPr)%Statut)=="MOBILE") & !
      !
      IF(TRIM(vCpn(iPr)%Statut)=="MOBILE" .OR. &
      &  TRIM(vCpn(iPr)%Statut)=="BUFFER") &
      & WRITE(FF,'(G15.6,A1,A,A1)',ADVANCE="NO") &
      & EXP(vSpc(vCpn(iPr)%iSpc)%Dat%LAct),  T_, &
      & vSpc(vCpn(iPr)%iSpc)%NamSp,          T_
      !
      WRITE(FF,*)
    
    ENDDO
    
    WRITE(FF,'(A)') "--"
    WRITE(FF,'(A)') "</------------------------- equilibrium results --"

  ENDIF
  !------------------------------------------------/for equilibrium runs

  WRITE(FF,'(A)') &
  & "!-----------------------------------------------------------------------"
  WRITE(FF,'(A)') "All Species, in order of decreasing Mole number"
  !
  ! order species by increasing abundance
  !,to apply permutation vPrm to array Arr: Arr=Arr(vPrm)
  CALL CalcPermut(vSpc(1:nAq)%Dat%Mole,vPrm)
  !
  WRITE(FF,'(4(A15,A1))') &
  & "________SPECIES",T_,"____MOLE_NUMBER",T_, &
  & "____ACTIV_COEFF",T_,"__LOG(ACTIVITY)",T_
  DO iAq=nAq,1,-1
    S=vSpc(vPrm(iAq))
    WRITE(FF,'(A15,A1,3(G15.6,A1))') &
    & TRIM(S%NamSp),  T_, &
    & S%Dat%Mole,     T_, &
    & EXP(S%Dat%LGam),T_, &
    & S%Dat%LAct/Ln10,T_
  ENDDO
  WRITE(FF,'(A)') &
  & "!-----------------------------------------------------------------------"
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Write_ShoDetail"

  RETURN
ENDSUBROUTINE Equil_Write_ShoDetail

ENDMODULE M_Equil_Write
