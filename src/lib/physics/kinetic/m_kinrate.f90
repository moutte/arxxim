MODULE M_KinRate
!--
!-- rate calculations on kinetic species
!--

  USE M_Kinds
  USE M_Trace,ONLY: T_,fHtm,Pause_
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC:: KinRate_CalcActivFactor
  PUBLIC:: KinRate_CalcQsK
  PUBLIC:: KinRate_SatState
  PUBLIC:: KinRate_CalcQsKFactor
  !
  PUBLIC:: KinRate_ActivTest
  
  REAL(dp),PARAMETER:: QsK_Max=  1.D+12 !
  LOGICAL, PARAMETER:: LimitQsK= .TRUE.
  
CONTAINS

SUBROUTINE KinRate_CalcActivFactor(&
!-----------------------------------------------------------------------
!-- compute VmAct, the factor depending on fluid chemistry only
!-----------------------------------------------------------------------
& cSatur, &  !IN:  "D" / "P" ! "DISSOLU"/"PRECIPI"
& M,      &  !IN:  kinetic model of mineral
& vLnAct, &  !IN:  fluid species activities
& VmAct)     !OUT: factor depending on fluid chemistry only
!! & dVmAdLnX_M) !OUT: its derivative vs ln(X)
  !
  USE M_T_KinFas,  ONLY: T_KinFas
  USE M_T_Kinmodel,ONLY: T_KinModel
  !
  CHARACTER,       INTENT(IN) :: cSatur
  TYPE(T_KinModel),INTENT(IN) :: M
  REAL(dp),        INTENT(IN) :: vLnAct(:)
  !
  REAL(dp), INTENT(OUT):: VmAct
  ! REAL(dp), INTENT(OUT):: dVmAdLnX_M(:) !(1:nAq)
  !
  INTEGER ::I
  REAL(dp)::X
  !
  ! dVmAdLnX_M=Zero
  !provisionally, dVmAdLnX_M is not computed in this new version  !!!!!!!!!!!!!
  !-> we assume that this factor will not be implicited           !!!!!!!!!!!!!
  !
  !! WRITE(51,*) -LnActH_/Ln10, CHAR(9), -LnActOH/Ln10
  VmAct=Zero
  SELECT CASE(cSatur)

  CASE("D") !"DISSOLU")
    DO I=1,M%NTermD
      X= EXP( M%N_d(I) *vLnAct(M%ISpcDiss(I)) )
      IF(M%JSpcDiss(I)>0) X= X * EXP( M%NJ_d(I)*vLnAct(M%JSpcDiss(I)) )
      VmAct= VmAct + M%Kd(I) * X
    ENDDO
    !IF(M%Special/=0) THEN !-> additional terms for carbonates, sulfides ...
    !  SELECT CASE(M%Special)
    !    CASE(1) !CALCITE
    !    CASE(2) !etcetera
    !    .......
    !  END SELECT
    !ENDIF

  CASE("P") !"PRECIPI")
    DO I=1,M%NTermP
      X= EXP( M%N_p(I) *vLnAct(M%ISpcPrec(I)) )
      IF(M%JSpcPrec(I)>0) X= X * EXP( M%NJ_p(I) *vLnAct(M%JSpcPrec(I)) )
      VmAct= VmAct + M%Kp(I) * X
    ENDDO

  END SELECT
  !
  RETURN
ENDSUBROUTINE KinRate_CalcActivFactor

!-----------------------------------------------------------------------
!-- calc. Q/K and d(Q/K)/dLnX of a phase of properties DG,vNu
!-- for a fluid composition vLnX
!-- (composition input in log(mole numbers))
!-----------------------------------------------------------------------
SUBROUTINE KinRate_CalcQsK( &
& nCi,       & !IN
& nCx,       & !IN
& DG0rt,     & !IN
& vNu,       & !IN
& vLnX,      & !IN
& vLnGam,    & !IN
& vLnActBuf, & !IN
& LnActW,    & !IN
& QsK,       & !OUT
& dQsKdLnXi)   !OUT

  USE M_Numeric_Const,     ONLY: MinExpDP,MaxExpDP
  USE M_Basis_Vars,ONLY: isW,MWSv
  !
  INTEGER, INTENT(IN):: nCi          !nr inert comp'nt
  INTEGER, INTENT(IN):: nCx          !nr mobile comp'nt
  REAL(dp),INTENT(IN):: DG0rt        !deltaG/RT of formation reaction
  REAL(dp),INTENT(IN):: vNu(:)       !stoikio of formation reaction
  REAL(dp),INTENT(IN):: vLnX(:)      !log(Mole Nrs. Prim.Species)  !1:nCi
  REAL(dp),INTENT(IN):: vLnGam(:)    !log(Act.Coeff. Prim.Species) !1:nCi
  REAL(dp),INTENT(IN):: vLnActBuf(:) !log(ActivityBufferSpecies)   !1:nCx
  REAL(dp),INTENT(IN):: LnActW       !log(Activity Solvent)
  !
  REAL(dp),INTENT(OUT):: QsK
  REAL(dp),INTENT(OUT):: dQsKdLnXi(:) !1:nCi
  !
  REAL(dp)::X
  !
  !lnAct= lnMolal + lnGamma=  lnMol + lnGam - lnMolH2O - lnMWH2O
  !QsK=EXP(DOT_PRODUCT(tNuMk(iMk,1:nCp),vLnAct(1:nCp)) - vDG_Mk(iMk))
  !or,directly from G0,
  !QsK=EXP(DOT_PRODUCT(tNuMk(iMk,1:nCp),vLnAct(1:nCp) +vSpc(1:nCp)%G0) - vKinFas(iMk)%G0)
  !
  !caveat:
  !Solvent should follow mole fraction scale, not molality ...
  !-> here, we assume explicited Solvent activity:
  X=            vNu(isW)  *LnActW &
  + DOT_PRODUCT(vNu(2:nCi),vLnX(2:nCi)) & !log(mole numbers)
  + DOT_PRODUCT(vNu(2:nCi),vLnGam(2:nCi)) & !log(gammas)
  -         SUM(vNu(2:nCi))*(vLnX(isW)+LOG(MWSv)) & !log(mass Solvent)
  - DG0rt !equiv -logK
  !
  IF(nCx>0) & != mobile species
  & X= X &
  &  + DOT_PRODUCT(vNu(nCi+1:nCi+nCx),vLnActBuf(1:nCx))
  !
  IF(X>MinExpDP .AND. X<MaxExpDP) THEN
    QsK=EXP(X)
  ELSE
    IF(X<=MinExpDP) QsK=EXP(MinExpDP)
    IF(X>=MaxExpDP) QsK=EXP(MaxExpDP)
  ENDIF
  !
  !<new 200911>
  IF(LimitQsK) QsK= MIN(QsK,QsK_Max)
  !</new>
  !
  dQsKdLnXi(:)=     Zero
  dQsKdLnXi(isW)= - SUM(vNu(2:nCi))*QsK
  dQsKdLnXi(2:nCi)=     vNu(2:nCi) *QsK
  !
  RETURN
ENDSUBROUTINE KinRate_CalcQsK

!-----------------------------------------------------------------------
!= from QsK deduce cSatur, for branching to dissol/precip
!-----------------------------------------------------------------------
SUBROUTINE KinRate_SatState( &
& M,         & !IN
& nMol,      & !IN
& nMolMinim, & !IN
& QsK,       & !IN
& Iota,      & !IN !mod 10/06/2008 17:02 added
& cSatur)      !OUT
  USE M_T_KinFas,ONLY: T_KinFas
  !
  TYPE(T_KinFas), INTENT(IN) :: M      !IN
  REAL(dp),       INTENT(IN) :: nMol   !IN, Nr Moles of mineral M
  REAL(dp),       INTENT(IN) :: nMolMinim !IN
  REAL(dp),       INTENT(IN) :: QsK    !IN
  REAL(dp),       INTENT(IN) :: Iota   !IN
  CHARACTER,      INTENT(OUT):: cSatur !OUT
  !
  cSatur=M%Dat%cSat != current value, cSatur will be the new value
  !
  cSatur="I" ! INERT, normal assumption by default
  IF(QsK < One - Iota) THEN
    IF(nMol<=nMolMinim) THEN  ;  cSatur="M" ! "MINIMAL"
    ELSE                      ;   cSatur="D"  ! "DISSOLU"
    ENDIF
  ELSEIF(QsK >= M%QsKSeuil + Iota) THEN
    cSatur="P" ! PRECIPI
  ENDIF
  !
ENDSUBROUTINE KinRate_SatState
  !
  ! IF ( QsK < One - Iota ) THEN
  ! IF ( nMol > nMolMinim ) THEN
  ! cSatur="DISSOLU"
  ! ELSE
  ! cSatur="MINIMAL"
  ! END IF
  ! ELSEIF ( QsK > QsKseuil + Iota ) THEN
  ! cSatur="PRECIPI"
  ! ELSE
  ! cSatur="INERT"
  ! END IF


  !one possible sequence:
  !!!  IF (QsK<One .AND. nMol>nMolMinim .AND. cSatur/="MINIMAL") THEN
  !!!    cSatur="DISSOLU"
  !!!  ELSEIF (QsK<=One.AND.nMol<=nMolMinim) THEN
  !!!    cSatur="MINIMAL" !QsK<1, but amount Too low to dissolve anymore
  !!!  ELSEIF (QsK>M%Dat%QsKSeuil .OR. (QsK>One .AND. M%Dat%cSat=="PRECIPI")) THEN
  !!!    cSatur="PRECIPI" !supersaturation or min. alREADy precipitating
  !!!  ENDIF
  !another possible sequence:
  !
  !!  !! IF(ABS(Qsk - One) < Iota) THEN
  !!  !!   cSatur="INERT"
  !!  !! ELSE!
  !!  !!   SELECT CASE(cSatur)
  !!  !!     CASE("PRECIPI"); IF(QsK<One)   cSatur="DISSOLU" !ELSE it continues as "PRECIPI"
  !!  !!     CASE("DISSOLU"); IF(QsK>One)   cSatur="PRECIPI" !ELSE it continues as "DISSOLU"
  !!  !!     CASE("MINIMAL")
  !!  !!       IF(QsK > One + Iota)        cSatur="PRECIPI"
  !!  !!     CASE("INERT")
  !!  !!       IF(QsK > M%Dat%QsKSeuil + Iota) cSatur="PRECIPI" !ELSE it continues as "INERT_"
  !!  !!       IF(QsK < One - Iota)        cSatur="DISSOLU" !
  !!  !!     !!CASE("PRIMARY")
  !!  !!     !!  IF(QsK > One + Iota)        cSatur="PRECIPI" !
  !!  !!     !!  IF(QsK < One - Iota)        cSatur="DISSOLU" !
  !!  !!     !!CASE("SECONDA")
  !!  !!     !!  IF(QsK > M%Dat%QsKSeuil + Iota) cSatur="PRECIPI"
  !!  !!   END SELECT
  !!  !! ENDIF
  !!  !! IF(cSatur=="DISSOLU" .AND. nMol<nMolMinim) cSatur="MINIMAL"
  !
  !SELECT CASE(cSatur)
  !CASE("PRECIPI"); IF(QsK<One)   cSatur="DISSOLU" !ELSE it continues as "PRECIPI"
  !CASE("DISSOLU"); IF(QsK>=One)  cSatur="PRECIPI" !ELSE it continues as "DISSOLU"
  !CASE("MINIMAL")
  !  IF(QsK >= One + Iota)        cSatur="PRECIPI"
  !CASE("INERT")
    !IF(QsK >= M%QsKSeuil + Iota) cSatur="PRECIPI" !ELSE it continues as "INERT_"
    !IF(QsK < One - Iota)         cSatur="DISSOLU" !
  !CASE("PRIMARY")
  !  IF(QsK >= One + Iota)        cSatur="PRECIPI" !
  !  IF(QsK < One - Iota)        cSatur="DISSOLU" !
  !CASE("SECONDA")
  !  IF(QsK >= M%QsKSeuil + Iota) cSatur="PRECIPI"
  !END SELECT
  !! IF(cSatur=="DISSOLU" .AND. nMol<=nMolMinim) cSatur="MINIMAL"
  !! IF(cSatur=="DISSOLU" .AND. nMol<=1.E-6) cSatur="MINIMAL"

SUBROUTINE KinRate_CalcQsKFactor(&
& cSatur,        & !IN
& M,             & !IN: Kinetic model
& QsK,           & !IN
& dQsKdLnXi,     & !IN
& VmQsK,         & !OUT
& dVmQdLnX_M)      !OUT
  !
  USE M_T_Kinmodel,ONLY: T_KinModel
  !
  CHARACTER,       INTENT(IN):: cSatur
  TYPE(T_KinModel),INTENT(IN):: M
  REAL(dp),        INTENT(IN):: QsK
  REAL(dp),        INTENT(IN):: dQsKdLnXi(:) !(1:nCi)
  !
  REAL(dp), INTENT(OUT):: VmQsK
  REAL(dp), INTENT(OUT):: dVmQdLnX_M(:) !(1:nAq)
  !
  REAL(dp):: dVmQdQsK
  INTEGER :: nCi_
  !
  nCi_= SIZE(dQsKdLnXi)
  VmQsK=     Zero !-> values for cSatur=="MINIMAL"
  dVmQdLnX_M=Zero !
  dVmQdQsK=  Zero
  !
  SELECT CASE(cSatur)

  CASE("D") !("DISSOLU")
  ! QsK<1 -> (One - QsK**M%AlfaD)>0 -> VmQsK<0
  ! QsK<1 -> -LOG(QsK)>0 -> VmQsK<0
    
    IF(M%BetaD /=Zero) THEN
      VmQsK=    - (One - QsK**M%AlfaD)**M%BetaD
      dVmQdQsK= M%BetaD &
      &         *(One - QsK**M%AlfaD)**(M%BetaD-One) &
      &         *M%AlfaD *QsK**(M%AlfaD-One)
    ELSE
    !cf Soler_Lasaga
      VmQsK=    (-LOG(QsK))**M%AlfaD
      dVmQdQsK= M%AlfaD *VmQsK /QsK /LOG(QsK)
    ENDIF

  CASE("P") !("PRECIPI")
  ! QsK>1 -> (QsK**M%AlfaD - One)>0 -> VmQsK>0
  ! QsK>1 -> LOG(QsK)>0 -> VmQsK>0
    
    IF(M%BetaP /=Zero) THEN
      VmQsK=(QsK**M%AlfaP - One)**M%BetaP
      !IF(VmF<0.0) THEN
      !  cSatur="INERT" !"I"nert mineral
      !  VmF=0.0
      !  !test!!!dVSdQsK=0.0
      !  RETURN
      !ENDIF
      dVmQdQsK= M%BetaP &
      &       *(QsK**M%AlfaP - One)**(M%BetaP-One) &
      &       *M%AlfaP *QsK**(M%AlfaP-One)
    ELSE
    !cf Soler_Lasaga
      VmQsK=    (LOG(QsK))**M%AlfaD
      dVmQdQsK= M%AlfaD *VmQsK /QsK /LOG(QsK)
    ENDIF

  CASE("M") !("MINIMAL")
    VmQsK=    Zero
    dVmQdQsK= Zero
    
  END SELECT

  dVmQdLnX_M(1:nCi_)= dVmQdQsK*dQsKdLnXi(1:nCi_)

  !! WRITE(51,'(A,3A1,3(G15.9,A1))') &
  !! & TRIM(M%Name),T_,cSatur,T_,LOG(QsK)/Ln10,T_,VmQsK,T_,QsK**M%AlfaD,T_
  RETURN
END SUBROUTINE KinRate_CalcQsKFactor

SUBROUTINE KinRate_ActivTest(fMnK,iCount) !,LnActH_,LnActOH,LnActCO2)
!-----------------------------------------------------------------------
!= test the "validity" of KinRate_CalcActivFactor
!= of module M_T_KinFas for rate calculations
!= (= the factor depending on species'activities)
!-----------------------------------------------------------------------
  USE M_IoTools,  ONLY: GetUnit
  USE M_Files,    ONLY: DirOut,cTitle,Files_Index_Write !,bOpenMnk
  !
  USE M_Numeric_Const,     ONLY: Ln10,TinyDP
  USE M_Dtb_Const, ONLY: T_CK
  !
  USE M_T_Species,  ONLY: Species_Index
  USE M_Basis_Vars, ONLY: nMx, isH_
  USE M_T_Kinmodel, ONLY: T_KinModel,KinModel_CalcCoeffs,KinModel_PrmIndex
  !
  USE M_System_Vars,ONLY: TdgK
  USE M_Global_Vars,ONLY: nAq,vSpc,vKinModel
  
  ! USE M_KinRate,    ONLY: KinRate_CalcActivFactor
  
  INTEGER,INTENT(INOUT):: fMnK
  INTEGER,INTENT(IN)   :: iCount
  !
  TYPE(T_KinModel),ALLOCATABLE::vKinMod(:)
  TYPE(T_KinModel)::M
  INTEGER :: I,iKm,nKm
  REAL(dp):: VmAct,X
  ! REAL(dp):: dum(1:nAq)
  REAL(dp):: vLnAct(1:nAq) !for dissolution rates
  INTEGER :: vPrm(1:nAq)
  REAL(dp):: pH
  
  !! print *,'KinRate_ActivTest'  ;  CALL Pause_
  
  nKm=SIZE(vKinModel)
  !
  IF(nKm==0) RETURN
  !
  DO iKm=1,nKm
    CALL KinModel_CalcCoeffs(vKinModel(iKm),One,TdgK)
  ENDDO
  !
  vLnAct(1:nAq)=vSpc(1:nAq)%Dat%LAct
  DO I=1,nAq; vPrm(i)=i; ENDDO
  !
  ALLOCATE(vKinMod(1:SIZE(vKinModel))); vKinMod=vKinModel
  !
  DO iKm=1,nKm
    CALL KinModel_PrmIndex(vPrm,vKinModel(iKm),vKinMod(iKm))
  ENDDO
  !
  pH=-vLnAct(isH_)/Ln10
  !
  IF(fMnk==0) THEN
  
    CALL GetUnit(fMnk)
    
    OPEN(fMnk,FILE=TRIM(DirOut)//"_vmact.restab")
    
    CALL Files_Index_Write(fHtm, &
    & TRIM(DirOut)//"_vmact.restab", &
    & "dependence of mineral dissolution rates on solution chemistry" )
    WRITE(fMnk,'(3(A,A1))',ADVANCE= 'NO') "count",T_,"TdgC",T_,"pH",T_ !,"pOH",T_
    
    DO iKm=1,nKm
      WRITE(fMnk,'(A12,A1)',ADVANCE= 'NO') vKinMod(iKm)%Name,T_
    ENDDO
    WRITE(fMnk,*)
  
  ENDIF
  !
  WRITE(fMnk,'(I3,A1,2(G15.3,A1))',ADVANCE= 'NO') iCount,T_,TdgK-T_CK,T_,pH,T_
  !
  DO iKm=1,nKm
    
    M=vKinMod(iKm)
    !iMn=Species_Index(M%Name,vSpc)-nAq-nMx
    
    !IF(iMn>0) THEN !mineral is not "mobile"
      CALL KinRate_CalcActivFactor(&
      & "D",      & !IN: "D"/"P"
      & M,        & !IN: kinetic model of mineral
      & vLnAct,   & !LnActH_,LnActOH,LnActCO2, & !"activators"
      & VmAct)      !OUT
      !! & dum) !dVmAdLnX_M) !OUT
      !X=ABS(VmAct)
      X=VmAct
      IF(X>TinyDP) THEN; X=LOG(X)/Ln10
      ELSE             ; X=Zero
      ENDIF
      WRITE(fMnk,'(G15.6,A1)',ADVANCE='NO') X,T_ !/(QsK-One)
    !ENDIF
    
  ENDDO
  !
  WRITE(fMnk,*)
  !
  DEALLOCATE(vKinMod)
  
  RETURN
END SUBROUTINE KinRate_ActivTest

END MODULE M_KinRate
