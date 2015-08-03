MODULE M_Equil_2
!--
!-- derived from M_Equil_3 --
!--
!-- calculate equilibrium speciation --
!-- considering saturation with minerals or gases --
!-- this version considers pure phases and mixtures --
!--
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_,Pause_
  IMPLICIT NONE

  PRIVATE

  PUBLIC:: Equil_Eq2

  REAL(dp),PARAMETER:: FasMinim= 1.D-6
  ! FasMinim= mimim' amount of phase
  REAL(dp),PARAMETER:: AffinIota= 1.D-3
  ! AffinIota= tolerance for equilibrium condition
  ! ( equilibrium IFF ABS(Affinity)<AffinIota )
  REAL(dp),PARAMETER:: MixMinim_TolX= 1.D-3
  ! for mixture minimisation (non ideal case)

  INTEGER, ALLOCATABLE:: tIPole(:,:)

  INTEGER,PARAMETER:: MaxIter= 30
  CHARACTER(LEN=30):: sFMT

  CHARACTER:: cTest_

CONTAINS

SUBROUTINE Equil_Eq2(NoMixture,iErr,cTest)
!--
!-- this version considers pure phases and mixtures
!-- (NoMixture - pure phase only )
!--
  USE M_IOTools,      ONLY: GetUnit
  USE M_Files,        ONLY: DirOut
  !
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Numeric_Tools,ONLY: iMinLoc_R
  !
  USE M_Dtb_Const,    ONLY: T_CK
  USE M_IoTools,      ONLY: OutStrVec
  !
  USE M_T_MixModel,   ONLY: MaxPole
  USE M_T_MixPhase,   ONLY: T_MixPhase
  USE M_T_Phase,      ONLY: T_Phase
  !
  USE M_Equil_Vars,   ONLY: T_EquPhase
  USE M_Equil_Solve,  ONLY: Equil_Solve
  !
  USE M_Global_Vars,  ONLY: vSpc,vFas !,tFormula
  USE M_Global_Vars,  ONLY: vMixModel
  USE M_Basis_Vars,   ONLY: tNuFas,vOrdPr
  USE M_Equil_Vars,   ONLY: cEquMode
  USE M_Equil_Vars,   ONLY: vYesList
  USE M_Equil_Vars,   ONLY: vDeltaG_Fs,vAffScale,vLnAct,vFasMole
  USE M_Equil_Vars,   ONLY: vEquFas,nEquFas
  USE M_Equil_Vars,   ONLY: vDeltaG_Eq,tNuEq
  !---------------------------------------------------------------------
  LOGICAL,INTENT(IN) :: NoMixture
  INTEGER,INTENT(OUT):: iErr
  CHARACTER,INTENT(IN),OPTIONAL:: cTest
  !---------------------------------------------------------------------
  INTEGER :: nCp,nFs,nPur,nSol
  INTEGER :: iFs,iMix,iPur
  INTEGER :: FF
  INTEGER :: I,J,K,C
  INTEGER :: nMix,nP !,nEquMix
  INTEGER :: iPrecip,iElimin
  INTEGER :: nIts,iDo0,iDo1
  ! REAL(dp):: X
  ! LOGICAL :: MixIsStable
  !
  INTEGER, ALLOCATABLE:: vIPole(:)
  REAL(dp),ALLOCATABLE:: vNuFas(:)
  !
  REAL(dp),ALLOCATABLE:: vFasAff(:)
  REAL(dp),ALLOCATABLE:: vFasAffPur(:)
  !
  TYPE(T_EquPhase),ALLOCATABLE:: vEquFas0(:)
  TYPE(T_MixPhase),ALLOCATABLE:: vMixFas0(:)
  TYPE(T_Phase),   ALLOCATABLE:: vFas0(:) !,vFasPur(:)
  TYPE(T_EquPhase):: EquFasNew
  !---------------------------------------------------------------------
  !! TYPE:: T_EquPhase
  !!   CHARACTER(LEN=23):: NamEq
  !!   INTEGER :: iPur
  !!   INTEGER :: iMix
  !!   REAL(dp):: vXPole(1:MaxPole)
  !!   !REAL(dp):: Grt
  !!   REAL(dp):: Mole
  !! END TYPE T_EquPhase

  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Eq2_Mix"

  IF(PRESENT(cTest)) THEN  ;  cTest_= cTest
  ELSE                     ;  cTest_= "1"
  ENDIF

  FF=0
  IF(iDebug>1) THEN
    !~ iCountEq=0
    !CALL Equil_FasTrace_EnTete(vFas,vYesList,FF)
    CALL GetUnit(FF)
    OPEN(FF,FILE=TRIM(DirOut)//"_equil.log")
  ENDIF

  cEquMode="EQ2"
  !-> flag in Equil_Solve, Equil_Residual, Equil_Jacobian

  nCp= SIZE(vOrdPr)
  ! nPur= COUNT(vFas(:)%Typ=="PURE")
  nPur= COUNT(vFas(:)%iSpc/=0)
  nMix= SIZE(vMixModel)
  nSol= 1 !! SIZE(vSolFas)
  nFs= nPur +nMix !+ nSol

  !! DEBUGG !!
  IF(iDebug>1) THEN
    WRITE(6,'(A)') "==Equil_Eq2=="
    PRINT *,"nPur=",nPur
    DO J=1,SIZE(vFas)
      !WRITE(6,'(A)')  TRIM(vFas(J)%NamFs)
      PRINT *,TRIM(vFas(J)%NamFs),vYesList(J)
    END DO
  END IF
  !PAUSE

  IF(NoMixture) nMix= 0

  ALLOCATE(vFasAff(nFs+nCp))
  ALLOCATE(vFasAffPur(nPur))
  ALLOCATE(vEquFas0(nCp))
  
  ALLOCATE(vFas0(nFs))
  vFas0(1:nPur)= vFas(1:nPur)
  
  !print *,"nSol, nPur, size(vFas), size(vFas0)=", &
  !&        nSol, nPur, size(vFas), size(vFas0)  ;  pause

  !IF(nSol>0) vFas0(nPur+1)= vFas(nPur+1)
  
  !--------------------------------------------------------- mixtures --
  IF(nMix>0) THEN

    DO K=1,nMix
      vFas0(nPur+nSol+K)%NamFs= vMixModel(K)%Name
    ENDDO

    ! for each mixing model,
    ! we shall compute ONE composition of minimal G
    ALLOCATE(vMixFas0(nMix))
    DO K=1,nMix
      vMixFas0(K)%Name= vMixModel(K)%Name
      vMixFas0(K)%iModel= K
      vMixFas0(K)%vXPole(:)= Zero
    ENDDO

    ALLOCATE(tIPole(nMix,MaxPole))
    tIPole(:,:)= 0
    CALL MixModel_Init(nPur,vFas0,vMixModel)

  ENDIF
  !--------------------------------------------------------/ mixtures --

  I=0
  DO iPur=1,nPur
    IF(vFas0(iPur)%Mole>FasMinim) THEN
      I=I+1
      vEquFas(I)%iPur= iPur
      vEquFas(I)%iMix= 0
      vEquFas(I)%NamEq= TRIM(vFas0(iPur)%NamFs)
      vEquFas(I)%Mole= vFas0(iPur)%Mole
    ENDIF
  ENDDO
  nEquFas= I

  !-------------------------------------- loop on mineral saturation, --
  !------------------------ EXIT when all minerals have vPurYes false --
  iDo0= 0
  Do0: DO
    iDo0= iDo0 +1

    IF(iDo0> MaxIter) THEN
      IF(iDebug>2) THEN
        WRITE(6,*) "iDo0>MaxIter in Equil_2"
        CALL Pause_
      ENDIF
      iErr= -7
      EXIT Do0
    ENDIF

    iDo1= 0
    !---------------------- solve the system until no negative phases --
    Do1: DO
      iDo1= iDo1 +1

      IF(iDebug>2) THEN
        WRITE(6,'(A)') "<-- vEquFas --"
        DO I=1,nEquFas
          CALL Trace_EquPhase(6,vEquFas(I))
        ENDDO
        WRITE(6,'(A)') "</- vEquFas --"
      ENDIF

      IF(nEquFas>0) THEN
        CALL EquFas_Clean
        CALL EquFas_Alloc(nCp,nEquFas)
        CALL EquFas_Update(nCp,nEquFas,vEquFas)
      ENDIF

      !------------------------------------------- solve the system --
      CALL Equil_Solve(nIts,iErr)
      !------------------------------------------/ solve the system --

      IF(nEquFas==0) EXIT Do1

      !----------------------------------- exclude "near zero" phases --
      vEquFas0(:)= vEquFas(:)
      J= 0
      DO I=1,nEquFas
        IF(ABS(vEquFas0(I)%Mole)>FasMinim) THEN
          J=J+1
          vEquFas(J)= vEquFas0(I)
        ENDIF
      ENDDO
      nEquFas= J
      !----------------------------------/ exclude "near zero" phases --
      
      !------------------------ check for phases with negative amount --
      IF(nEquFas>0) THEN
        !
        !------------------------------------------------------ trace --
        IF(iDebug>2) THEN
          WRITE(6,'(A)') "<-- CURRENT --"
          DO I=1,nEquFas
            CALL Trace_EquPhase(6,vEquFas(I))
          ENDDO
          WRITE(6,'(A)') "</- CURRENT --"
        ENDIF
        !------------------------------------------------------/trace --
        !
        !--------------------- remove most "negative" phase and CYCLE --
        IF(MINVAL(vEquFas(1:nEquFas)%Mole) < -FasMinim) THEN

          iElimin= iMinLoc_R(vEquFas(1:nEquFas)%Mole)

          IF(iDebug>2) &
          & WRITE(6,'(16X,2A)') "NEG= ", TRIM(vEquFas(iElimin)%NamEq)

          CALL EquFas_Remove(iElimin,nEquFas)

          CYCLE Do1 

        ENDIF
        !---/
        !
      ENDIF
      !
      EXIT Do1
      !-----------------------/ check for phases with negative amount --
      !
    ENDDO Do1
    !---------------------/ solve the system until no negative phases --

    IF(iErr<0) EXIT Do0 != error in Newton -> EXIT

    vFasAff= Zero
    vFasAffPur= Zero

    DO iPur=1,nPur

      IF(vYesList(iPur)) &
      & vFasAffPur(iPur)= &
      & (  vDeltaG_Fs(iPur) &
      &  - DOT_PRODUCT(tNuFas(iPur,1:nCp),vLnAct(vOrdPr(1:nCp))) )

      ! scaling ala EQ3/6
      IF(vYesList(iPur)) &
      & vFasAff(iPur)= vFasAffPur(iPur) /SUM(ABS(tNuFas(iPur,:))) !vAffScale(iPur)

    ENDDO

    !-----------------------------------------------------------mixtures
    IF(nMix>0) THEN

      DO iMix=1,nMix

        ! for each mixing model,
        ! compute the composition that minimize
        ! the affinity of formation from end-members
        ! at current affinity w.r.t. the aqu'phase

        !vIPole(:)= tIPole(iMix,:)
        CALL Mixture_Minimize( &
        & vFasAffPur(:),tIPole(iMix,:),vMixModel(iMix), &
        & vMixFas0(iMix))

        vFasAff(nPur+iMix)= vMixFas0(iMix)%Grt

        !--------------------------------------------------------scaling
        nP= vMixModel(iMix)%NPole
        ALLOCATE(vNuFas(nCp))
        ALLOCATE(vIPole(nP))
        vIPole(1:nP)= tIPole(iMix,1:nP) !vMixModel(iMix)%vIPole(1:nP)
        DO C=1,nCp
          vNuFas(C)= Zero
          DO K=1,nP
            vNuFas(C)= vNuFas(C) &
            &        + vMixFas0(iMix)%vXPole(K) *tNuFas(vIPole(K),C)
          ENDDO
          !~ ! vNuFas(C)= SUM( vMixFas0(iMix)%vXPole(1:nP) &
          !~ ! &              * tNuFas(vIPole(1:nP),C) )
        ENDDO
        DEALLOCATE(vIPole)
        vFasAff(nPur+iMix)= vMixFas0(iMix)%Grt /SUM(ABS(vNuFas(:)))
        DEALLOCATE(vNuFas)
        !-------------------------------------------------------/scaling

      ENDDO

    ENDIF
    !----------------------------------------------------------/mixtures

    !--------------------------------------------------------------trace
    IF(iDebug>2) THEN
      WRITE(6,'(A)') "<-- OVERSAT --"
      ! print *,size(vFasAff),nPur  ;  pause
      DO I=1,nPur
        IF (vYesList(I)) &
        ! .AND. vFasAff(I)<AffinIota) &
        & WRITE(6,'(A,G15.6,1X,A15,1X)') &
        & "Aff=",vFasAff(I),vFas0(I)%NamFs
      ENDDO
      DO I=1,nMix
        IF (vFasAff(nPur+I)<AffinIota) THEN
          WRITE(6,'(A,G15.6)',ADVANCE="NO") "Aff=",vFasAff(nPur+I)
          WRITE(6,'(A)',ADVANCE="NO") " X(:)="
          DO J=1,vMixModel(I)%nPole
            WRITE(6,'(G15.6,1X)',ADVANCE="NO") vMixFas0(I)%vXPole(J)
          ENDDO
          WRITE(6,'(A)') vMixModel(I)%Name
        ENDIF
      ENDDO
      WRITE(6,'(A)') "</- OVERSAT --"
      CALL Pause_
    ENDIF
    !-------------------------------------------------------------/trace

    !---------------------add phase with highest supersaturation, if any
    !---------if no phase found with affinity below -AffinIota, EXIT Do0
    IF(MINVAL(vFasAff)<-AffinIota) THEN

      iPrecip= iMinLoc_R(vFasAff)
      !-> the saturated phase with lowest vFasAff
      !~ print *,"vFas(iPp)%NamFs= ",vFas(iPp)%NamFs  ;  pause
      !~ vNewFas(iPrecip)= .TRUE.

    ELSE

      iPrecip= 0 != there is no saturated phase ==
      EXIT Do0 ! no phase found with affinity below -AffinIota > EXIT ==

    ENDIF

    !-------------------------------------------iPrecip>0, add new phase
    
    IF(iPrecip<=nPur) THEN
    !--------------------------------------------------new phase is pure

      EquFasNew%iPur= iPrecip
      EquFasNew%iMix= 0
      EquFasNew%NamEq= TRIM(vFas0(iPrecip)%NamFs)
      EquFasNew%Mole= Zero

    ELSE
    !-----------------------------------------------new phase is mixture

      iMix= iPrecip -nPur

      EquFasNew%iPur= 0
      EquFasNew%iMix= iMix
      EquFasNew%NamEq= TRIM(vMixModel(iMix)%Name)
      nP= vMixModel(iMix)%NPole
      EquFasNew%nPole= nP
      EquFasNew%vIPole(1:nP)= tIPole(iMix,1:nP)

      EquFasNew%vXPole(1:nP)= vMixFas0(iMix)%vXPole(1:nP)
      EquFasNew%vLnAct(1:nP)= vMixFas0(iMix)%vLnAct(1:nP)
      EquFasNew%Mole= Zero

    ENDIF

    CALL CheckAssemblage(FF,nCp,nEquFas,EquFasNew,vEquFas,iElimin)

    !----------------------if necessary, remove one of the active phases
    IF(iElimin>0) THEN

      CALL EquFas_Remove(iElimin,nEquFas)

      IF(iDebug>2) &
      & WRITE(6,'(16X,2A)') "OUT= ",TRIM(vEquFas(iElimin)%NamEq)

    ENDIF
    !------------------------------------/

    vEquFas(nEquFas+1)= EquFasNew
    nEquFas= nEquFas +1

    IF(iDebug>2) &
    & WRITE(6,'(16X,2A)') "ADD= ",TRIM(EquFasNew%NamEq)
    !-----------------------------------------------------/add new phase
    
    !--------------------------------------------------------------trace
    IF(iDebug>3) THEN
      DO I=1,nEquFas
        iFs= vEquFas(I)%iPur
        IF(iFs>0) THEN
          WRITE(6,'(A,I4,1X,A15,2(A,G15.6))') &
          & "Do0",iDo0, &
          & vFas0(iFs)%NamFs, &
          & "/ lQsK=", -vFasAff(iFs)/Ln10 ,&
          & "/ Mole=",vEquFas(I)%Mole
        ENDIF
      ENDDO
      !~ IF(iDebug>3)
      CALL Pause_
    ENDIF
    !-------------------------------------------------------------/trace
    !
    !----------------/ add phase with highest supersaturation, if any --
    !
    IF(nEquFas==0) EXIT Do0
    !
    IF(iDebug>2) PRINT *,"============================================="
    !~ IF(iDebug>2) PAUSE

  ENDDO Do0
  !-------------------------------------/ loop on mineral saturation, --

  !~ IF (ALLOCATED(vFasMole)) DEALLOCATE(vFasMole)
  !~ ALLOCATE(vFasMole(nPur+nMix))
  !~ vFasMole(:)= Zero
  !~ DO I=1,nEquFas
    !~ IF(vEquFas(I)%iPur>0) vFasMole(vEquFas(I)%iPur)= vEquFas(I)%Mole
    !~ IF(vEquFas(I)%iMix>0) vFasMole(vEquFas(I)%iMix+nPur)= vEquFas(I)%Mole
  !~ ENDDO

  IF(nMix>0) CALL Check_vXMean

  CALL EquFas_Clean

  DEALLOCATE(vFas0)
  DEALLOCATE(vFasAff)
  DEALLOCATE(vFasAffPur)
  DEALLOCATE(vEquFas0)

  IF(nMix>0) THEN
    DEALLOCATE(vMixFas0)
    DEALLOCATE(tIPole)
  ENDIF

  IF(FF>0) CLOSE(FF)

  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Eq2_Mix"
  !
CONTAINS

SUBROUTINE Check_vXMean
!--
!-- for each mixing model with nFas>1
!-- compute the free energy of the mixture
!-- with the mean composition
!--
  USE M_System_Vars,ONLY: TdgK,Pbar
  !
  !--- T_TabPhase: same TYPE as used for GEM routine:
  ! for a given mixing model,
  ! there can be, in the stable assemblage,
  ! up to NPole phases of different compositions.
  ! T_TabPhase is design to contain the different compositions
  ! that may be active for a same mixing model.
  ! nFas is the number of active phases
  TYPE:: T_TabPhase
    CHARACTER(LEN=23):: NamTbl
    INTEGER :: nFas                        ! -> number of phases
    REAL(dp):: tXPole(1:Maxpole,1:MaxPole) ! table of phase compositions
    REAL(dp):: vMole(1:MaxPole)            ! mole nr of each phase
  END TYPE T_TabPhase
  TYPE(T_TabPhase),ALLOCATABLE:: vTabPhase(:)
  !---/(same as for GEM routine)
  !
  TYPE(T_EquPhase):: EquFas
  INTEGER :: I,J,K,P,nP,N,nF
  INTEGER :: nTabFas
  REAL(dp):: G
  LOGICAL :: MustUpdate
  REAL(dp),ALLOCATABLE:: vX(:)
  REAL(dp),ALLOCATABLE:: vG0(:)

  MustUpdate= .FALSE.
  ALLOCATE(vTabPhase(nMix))
  vTabPhase(:)%nFas= 0
  nTabFas= 0

  DO I=1,nEquFas
    EquFas= vEquFas(I)
    IF(EquFas%iMix>0) THEN
      J= EquFas%iMix
      vTabPhase(J)%NamTbl= TRIM(EquFas%NamEq)
      vTabPhase(J)%nFas= vTabPhase(J)%nFas +1
      vTabPhase(J)%vMole(vTabPhase(J)%nFas)= EquFas%Mole
      nP= vMixModel(J)%NPole
      DO P=1,nP
        vTabPhase(J)%tXPole(vTabPhase(J)%nFas,P)= EquFas%vXPole(P)
      ENDDO
    ENDIF
  END DO

  IF(COUNT(vTabPhase(:)%nFas>0)==0) RETURN

  DO J=1,nMix

    nF= vTabPhase(J)%nFas
    IF (nF>1) THEN

      nP= vMixModel(J)%NPole

      !-- compute the mean composition of the mixture
      ALLOCATE(vX(nP))
      vX(:)= Zero
      DO K=1,nF
        vX(1:nP)= vX(1:nP) &
        &       + vTabPhase(J)%tXPole(K,1:nP) &
        &         *vTabPhase(J)%vMole(K)
      END DO
      vX(:)= vX(:) /SUM(vX(:))

      !-- compute its free energy,
      !-- using vFasAffPur(:),tIPole(J,:),vMixModel(J)
      ALLOCATE(vG0(nP))
      DO P=1,nP
        vG0(P)= vFasAffPur(tIPole(J,P))
      END DO
      !
      G= MixPhase_Grt( &
      & TdgK,Pbar, &
      & nP, &
      & vMixModel(J), &
      & vG0(1:nP), &
      & vX(1:nP))
      !
      print *,"G=",G
      ! if G is near zero,
      ! then there is only one phase, with the composition vX
      IF(ABS(G)<AffinIota *1.D1) THEN
        MustUpdate= .TRUE.
        vTabPhase(J)%nFas= 1
        vTabPhase(J)%tXPole(1,1:nP)= vX(1:nP)
        vTabPhase(J)%vMole(1)= SUM(vTabPhase(J)%vMole(1:nF))
      END IF

      DEALLOCATE(vG0)
      DEALLOCATE(vX)

    END IF

  END DO

  IF(MustUpdate) THEN

    K= 0
    vEquFas0(:)= vEquFas(:)
    vEquFas(:)%Mole= Zero

    DO I=1,nEquFas
      IF(vEquFas0(I)%iMix==0) THEN
        K=K+1
        vEquFas(K)= vEquFas0(I)
      ENDIF
    ENDDO

    DO J=1,nMix
      IF(vTabPhase(J)%nFas>0) THEN
        !~ ! SUM( vTabPhase(J)%vMole(1:vTabPhase(J)%nFas) )
        !~ J= vEquFas0(I)%iMix
        !~ print *, "J=", J  ;  pause
        DO N= 1,vTabPhase(J)%nFas
          K=K+1
          vEquFas(K)%NamEq= TRIM(vTabPhase(J)%NamTbl)
          vEquFas(K)%iPur= 0
          vEquFas(K)%iMix= J
          vEquFas(K)%vXPole(1:nP)= vTabPhase(J)%tXPole(N,1:nP)
          vEquFas(K)%Mole= vTabPhase(J)%vMole(N)
        END DO
      ENDIF
    END DO

    nEquFas= K

  ENDIF

  DEALLOCATE(vTabPhase)

  RETURN
END SUBROUTINE Check_vXMean

SUBROUTINE EquFas_Remove(N,M)
  INTEGER,INTENT(IN)   :: N
  INTEGER,INTENT(INOUT):: M

  INTEGER:: I,J

  vEquFas0(:)= vEquFas(:)
  J= 0
  DO I=1,M
    IF(I/=N) THEN
      J=J+1
      vEquFas(J)= vEquFas0(I)
    ENDIF
  ENDDO
  M= J

  RETURN
END SUBROUTINE EquFas_Remove

ENDSUBROUTINE Equil_Eq2

SUBROUTINE EquFas_Alloc(nC,nEq)
  USE M_Equil_Vars,ONLY: vDeltaG_Eq,tAlfEq,tNuEq
  !
  INTEGER,INTENT(IN):: nC,nEq
  !
  ALLOCATE(vDeltaG_Eq(nEq))
  ALLOCATE(tAlfEq(nC,nEq))
  ALLOCATE(tNuEq(nEq,nC))
  !
  RETURN
ENDSUBROUTINE EquFas_Alloc

SUBROUTINE EquFas_Update(nCp,nEquFas,vEquFas)
!--
!-- update vDeltaG_Eq, tAlfEq, tNuEq,
!-- which are used in Equ_Residual
!--
  USE M_T_MixModel, ONLY: MaxPole
  USE M_Global_Vars,ONLY: vMixModel
  USE M_Equil_Vars, ONLY: T_EquPhase
  USE M_Basis_Vars, ONLY: tAlfFs,tNuFas
  USE M_Equil_Vars, ONLY: vDeltaG_Eq,tAlfEq,tNuEq
  USE M_Equil_Vars, ONLY: vDeltaG_Fs
  !
  INTEGER,INTENT(IN):: nCp
  INTEGER,INTENT(IN):: nEquFas
  TYPE(T_EquPhase),INTENT(IN):: vEquFas(:)

  INTEGER:: iPur,iMix,I,C

  DO I=1,nEquFas

    iPur= vEquFas(I)%iPur
    iMix= vEquFas(I)%iMix

    IF(iPur>0) THEN
      vDeltaG_Eq(I)= vDeltaG_Fs(iPur)
      tAlfEq(:,I)= tAlfFs(:,iPur)
      tNuEq(I,:)= tNuFas(iPur,:)
    ENDIF

    IF(iMix>0) CALL EquFasMix_Update(I,nCp,vEquFas(I))

  ENDDO

  WRITE(sFMT,'(a,i3,a)') '(A,1X,',nCp,'(G12.3,1X))'

  IF(iDebug>2) THEN
    WRITE(75,'(A)') "===========================EquFas_Update=========="
    DO I=1,nEquFas
      WRITE(75,'(A,G15.6,1X,A)') &
      & "DeltaG_Eq= ",vDeltaG_Eq(I),TRIM(vEquFas(I)%NamEq)
      WRITE(75,sFMT) "tAlfEq=    ",(tAlfEq(C,I),C=1,nCp)
      WRITE(75,sFMT) "tNuEq=     ",(tNuEq(I,C),C=1,nCp)
    ENDDO
  ENDIF
  !
  RETURN
ENDSUBROUTINE EquFas_Update

SUBROUTINE EquFasMix_Update(I,nCp,EquFas)
!--
!-- compute vDeltaG_Eq(I),tAlfEq(:,I),tNuEq(I,:)
!-- using EquFas%vXPole(:), and vDeltaG_Fs,tAlfFs,tNuFas
!--
  USE M_Equil_Vars, ONLY: T_EquPhase
  USE M_Global_Vars,ONLY: vMixModel
  USE M_Basis_Vars, ONLY: tAlfFs,tNuFas
  USE M_Equil_Vars, ONLY: vDeltaG_Eq,tAlfEq,tNuEq,vDeltaG_Fs
  !
  INTEGER,INTENT(IN):: I
  INTEGER,INTENT(IN):: nCp
  TYPE(T_EquPhase),INTENT(IN):: EquFas
  !
  INTEGER:: C,nP
  INTEGER,ALLOCATABLE:: vIPole(:) ! indexes of end-members in vFas

  nP= EquFas%NPole
  ALLOCATE(vIPole(nP))
  vIPole(1:nP)= EquFas%vIPole(1:nP)

  vDeltaG_Eq(I)= SUM( EquFas%vXPole(1:nP) &
  &                 * vDeltaG_Fs(vIPole(1:nP)) ) &
  &            + SUM( EquFas%vXPole(1:nP) &
  &                 * EquFas%vLnAct(1:nP) )

  DO C=1,nCp
    tAlfEq(C,I)= SUM( EquFas%vXPole(1:nP) &
    &               * tAlfFs(C,vIPole(1:nP)) )
    tNuEq(I,C)= SUM( EquFas%vXPole(1:nP) &
    &              * tNuFas(vIPole(1:nP),C) )
  ENDDO

  DEALLOCATE(vIPole)

  RETURN
ENDSUBROUTINE EquFasMix_Update

SUBROUTINE Trace_EquPhase(iFile,Fas)
  USE M_Equil_Vars,ONLY: T_EquPhase
  
  INTEGER,         INTENT(IN):: iFile
  TYPE(T_EquPhase),INTENT(IN):: Fas

  INTEGER:: J
  
  WRITE(iFile,'(A,G15.6,1X)',ADVANCE="NO") "N=",Fas%Mole
  IF(Fas%iMix>0) THEN
    WRITE(iFile,'(A)',ADVANCE="NO") " X(:)="
    DO J=1,Fas%nPole
      WRITE(iFile,'(G15.6,1X)',ADVANCE="NO") Fas%vXPole(J)
    ENDDO
  ENDIF
  WRITE(iFile,'(A)') TRIM(Fas%NamEq)

ENDSUBROUTINE Trace_EquPhase

SUBROUTINE EquFas_Clean
  USE M_Equil_Vars,ONLY: vDeltaG_Eq,tAlfEq,tNuEq

  IF(ALLOCATED(vDeltaG_Eq)) DEALLOCATE(vDeltaG_Eq)
  IF(ALLOCATED(tAlfEq))     DEALLOCATE(tAlfEq)
  IF(ALLOCATED(tNuEq))      DEALLOCATE(tNuEq)

  RETURN
ENDSUBROUTINE EquFas_Clean

SUBROUTINE CheckAssemblage( & ! IN
& FF,                       & ! IN
& nCp,                      & ! IN
& nEquFas,                  & ! IN
& EquFasNew,                & ! IN, the new phase
& vEquFas,                  & ! IN, the current phase set
& iElimin)                    ! OUT
!--
!-- check whether new species, indexed vFas(iPp)%iPur
!-- is independent of current phase assemblage
!--
  USE M_Numeric_Mat,  ONLY: LU_BakSub
  USE M_Numeric_Tools,ONLY: iMinLoc_R,iFirstLoc,iMaxLoc_R
  USE M_Basis_Tools,  ONLY: Basis_FreeSet_Select
  USE M_T_MixModel,   ONLY: MaxPole
  !
  USE M_Global_Vars,  ONLY: vFas,vSpc,vMixModel
  USE M_Basis_Vars,   ONLY: tAlfPr,tAlfFs,vOrdPr
  USE M_Equil_Vars,   ONLY: T_EquPhase
  !---------------------------------------------------------------------
  INTEGER, INTENT(IN)   :: FF
  INTEGER, INTENT(IN)   :: nCp
  INTEGER, INTENT(IN)   :: nEquFas
  TYPE(T_EquPhase),INTENT(IN):: EquFasNew
  TYPE(T_EquPhase),INTENT(IN):: vEquFas(:)
  INTEGER, INTENT(OUT)  :: iElimin
  !---------------------------------------------------------------------
  INTEGER :: iPur,iMix,I,J,C,nP
  LOGICAL :: IsNotFree
  LOGICAL :: bSingul
  INTEGER :: vIPole(MaxPole)
  !
  ! INTEGER :: iElimin1, iElimin2
  ! REAL(dp):: xMin, xMax
  !
  REAL(dp),ALLOCATABLE:: tStoikCpn(:,:),tStoikCpn2(:,:)
  REAL(dp),ALLOCATABLE:: tBase(:,:)
  INTEGER, ALLOCATABLE:: vIndex(:)
  ! REAL(dp),ALLOCATABLE:: tTransform(:,:)
  INTEGER, ALLOCATABLE:: vIndx(:)
  REAL(dp),ALLOCATABLE:: vY(:)
  REAL(dp),ALLOCATABLE:: vTest(:),vTest2(:)
  !---------------------------------------------------------------------
  !~ IF(iDebug>2) WRITE(6,'(A)') "< CheckAssemblage"
  !
  IF(iDebug>2) WRITE(FF,'(2A)') "NEW= ", EquFasNew%NamEq
  !
  iElimin= 0
  !
  ALLOCATE(vIndex(nCp))
  vIndex(:)=0
  DO i=1,nEquFas+2
    vIndex(i)= i
  ENDDO
  ALLOCATE(tStoikCpn(nCp,nEquFas+2))
  !
  !----------------------------------------------------compute tStoikCpn
  !----------------------= stoikio table of current set of stable phases
  tStoikCpn(:,1)= tAlfPr(:,1) ! water
  vIndex(1)= 1
  !
  DO I=1,nEquFas

    iPur= vEquFas(I)%iPur
    iMix= vEquFas(I)%iMix

    IF(iPur>0) tStoikCpn(:,1+I)= tAlfFs(:,iPur)

    IF(iMix>0) THEN
      nP= vEquFas(I)%nPole
      vIPole(1:nP)= vEquFas(I)%vIPole(1:nP)
      DO C=1,nCp
        tStoikCpn(C,1+I)= SUM( vEquFas(I)%vXPole(1:nP) &
        &                    * tAlfFs(C,vIPole(1:nP)) )
      ENDDO
    ENDIF

  ENDDO

  iPur= EquFasNew%iPur
  iMix= EquFasNew%iMix

  IF(iPur>0) tStoikCpn(:,nEquFas+2)= tAlfFs(:,EquFasNew%iPur)

  IF(iMix>0) THEN
    nP= EquFasNew%nPole
    vIPole(1:nP)= EquFasNew%vIPole(1:nP)
    DO C=1,nCp
      tStoikCpn(C,nEquFas+2)= SUM( EquFasNew%vXPole(1:nP) &
      &                          * tAlfFs(C,vIPole(1:nP)) )
    ENDDO
  ENDIF
  !---------------------------------------------------/compute tStoikCpn

  !----------------------------------------------------------------trace
  IF(iDebug>2) THEN
    do i=1,nEquFas
      do C=1,nCp
        write(71,'(F7.3,1X)',advance="no") tStoikCpn(C,i+1)
      enddo
      write(71,'(A)') TRIM(vEquFas(i)%NamEq)
    enddo
    do C=1,nCp
      write(71,'(F7.3,1X)',advance="no") tStoikCpn(C,nEquFas+2)
    enddo
    write(71,'(A)') TRIM(EquFasNew%NamEq)
    write(71,'(A)') "=================================================="
  ENDIF
  !---------------------------------------------------------------/trace

  !-- using Gaussian elimination (without permutations),
  !-- find a set of independent species (of stoikios tAlfPr)
  !-- independent also of those in tStoikioCpn
  CALL Basis_FreeSet_Select( & !
  & tAlfPr,      & !IN
  & vIndex,      & !OUT
  & IsNotFree,   & !OUT
  & tStoikioCpn= tStoikCpn)     !IN
  !
  DEALLOCATE(vIndex)

  IF(IsNotFree) THEN
    !-- new phase not independent from current assemblage
    !-- -> build a basis without the new phase,
    !-- and compute stoichio of new phase versus current assemblage
    !
    ALLOCATE(tStoikCpn2(nCp,nEquFas+1))
    tStoikCpn2(:,1:nEquFas+1)= tStoikCpn(:,1:nEquFas+1)

    ALLOCATE(vIndex(nCp))
    vIndex(:)=0
    DO i=1,nEquFas+1
      vIndex(i)= i
    ENDDO

    !-- find a set of independent species (of stoikios tAlfPr)
    !-- independent also of those in tStoikioCpn2
    !-- vIndex(:) gives the indexes in tAlfPr of these independent species
    CALL Basis_FreeSet_Select( & !
    & tAlfPr,                  & !IN
    & vIndex,                  & !OUT
    & IsNotFree,               & !OUT
    & tStoikioCpn= tStoikCpn2)   !IN

    ALLOCATE(tBase(nCp,nCp))
    tBase(:,1:nEquFas+1)= tStoikCpn(:,1:nEquFas+1)
    DO I=nEquFas+2,nCp
      tBase(:,I)= tAlfPr(:,vIndex(I))
    ENDDO

    !--------------------------------------------------------------trace
    IF(iDebug>2) THEN
      WRITE(73,'(A)') "==CheckAssemblage=="
      DO I=1,nCp
        DO J=1,nCp
          WRITE(73,'(G11.2,1X)',ADVANCE="NO") tBase(J,I)
        ENDDO
        WRITE(73,*)
      ENDDO
      WRITE(73,'(A)') "==**=="
    ENDIF
    !-------------------------------------------------------------/trace

    ALLOCATE(vIndx(nCp))
    ALLOCATE(vY(nCp))

    !------------------------------------------LU decomposition of tBase
    CALL Compute_Transform(tBase,vIndx,bSingul)
    IF(bSingul) CALL Stop_("SINGUL IN CheckAssemblage")

    iPur= EquFasNew%iPur
    iMix= EquFasNew%iMix
    IF(iPur>0) vY(:)= tAlfFs(:,iPur)
    IF(iMix>0) THEN
      nP= EquFasNew%nPole
      vIPole(1:nP)= EquFasNew%vIPole(1:nP)
      DO C=1,nCp
        vY(C)= SUM( EquFasNew%vXPole(1:nP) &
        &         * tAlfFs(C,vIPole(1:nP)) )
      ENDDO
    ENDIF
    CALL LU_BakSub(tBase,vIndx,vY)

    !-----------------------------------find which phase must be removed
    ALLOCATE(vTest(nEquFas))   ;  vTest= Zero
    ALLOCATE(vTest2(nEquFas))  ;  vTest2= 1.D9
    
    DO I=1,nEquFas
      IF(ABS(vY(I+1))>1.D-9) vTest2(I)= vEquFas(I)%Mole/vY(I+1)
      vTest(I)= ABS(vY(I+1)) /vEquFas(I)%Mole
    ENDDO

    !--------------------------------------------------------------trace
    IF(iDebug>2) THEN
      WRITE(73,'(2A)') "EquFasNew=",TRIM(EquFasNew%NamEq)
      DO I=1,nEquFas
        WRITE(73,'(3(G15.6,1X),A)') &
        & vY(I+1),vEquFas(I)%Mole,vTest(I),TRIM(vEquFas(I)%NamEq)
        !~ WRITE(73,'(2(G15.6,1X),A)') vY(I+1),vTest(I),TRIM(vEquFas(I)%NamEq)
      ENDDO
    ENDIF
    !-------------------------------------------------------------/trace

    !!X= MAXVAL(vTest)
    SELECT CASE(cTest_)
    CASE("2")  ;  iElimin= iMinLoc_R(vTest2)
    CASE("1")  ;  iElimin= iMaxLoc_R(vTest)
    END SELECT

    DEALLOCATE(vTest)
    DEALLOCATE(vTest2)

    !----------------------------------/find which phase must be removed
    !
    DEALLOCATE(vIndx)
    DEALLOCATE(vY)
    DEALLOCATE(tBase)
    !
    DEALLOCATE(vIndex)
    DEALLOCATE(tStoikCpn2)
    !
  ENDIF
  !
  DEALLOCATE(tStoikCpn)
  !
  RETURN
ENDSUBROUTINE CheckAssemblage

SUBROUTINE MixModel_Init(nPur,vFas,vMixModel)
!--
!-- initialize tIPole(iMix,iP),
!-- - address of e-m iP of mixture iMix in vSpc or vFas
!--
  USE M_T_Phase,    ONLY: T_Phase
  USE M_T_MixModel, ONLY: T_MixModel

  INTEGER,         INTENT(IN):: nPur
  TYPE(T_Phase),   INTENT(IN):: vFas(:)
  TYPE(T_MixModel),INTENT(IN):: vMixModel(:)

  INTEGER:: iModel,P,P1,P2,iPur,nP
  !
  DO iModel=1,SIZE(vMixModel)

    nP= vMixModel(iModel)%NPole

    ! find the indexes of end-members of mixture vMixFas(iMix)
    ! in the pure phase list vFas(1:nPur)
    ! -> indices in tIPole(iMix,:)
    !~ print *,TRIM(vMixModel(iModel)%Name)
    DO P=1,nP
      P1= vMixModel(iModel)%vIPole(P)
      !-> index in vSpc
      !-> find the pure phases that points to this species
      P2= 0
      DO iPur=1,nPur
        IF(vFas(iPur)%iSpc==P1) THEN
          P2= iPur
          !~ vPurIsPole(iPur)= .TRUE.
          EXIT
        ENDIF
      ENDDO
      tIPole(iModel,P)= P2
      !~ print *,TRIM(vFas(P2)%NamFs)
      !
    ENDDO

  ENDDO
  !
  RETURN
ENDSUBROUTINE MixModel_Init

SUBROUTINE Compute_Transform(tTransform,vIndx,Error)
  USE M_Numeric_Mat,ONLY: LU_Decomp
  !
  REAL(dp),INTENT(INOUT):: tTransform(:,:)
  INTEGER, INTENT(OUT)  :: vIndx(:)
  LOGICAL, INTENT(OUT)  :: Error
  !
  REAL(dp):: D
  !
  ! the transformation matrix, tTransform, must be invertible !!
  CALL LU_Decomp(tTransform,vIndx,D,Error)
  !
END SUBROUTINE Compute_Transform

REAL(dp) FUNCTION MixPhase_Grt(TdgK,Pbar,nP,MM,vMu0rt,vX)
  USE M_T_MixModel,ONLY: T_MixModel,MixModel_Activities
  !----------------------------------------------------------------inout
  REAL(dp),        INTENT(IN):: TdgK,Pbar
  INTEGER,         INTENT(IN):: nP
  TYPE(T_MixModel),INTENT(IN):: MM        ! mixing model
  REAL(dp),        INTENT(IN):: vMu0rt(:) !
  REAL(dp),        INTENT(IN):: vX(:)     ! phase composition
  !-----------------------------------------------------------------vars
  REAL(dp):: vLGam(nP),vLIdeal(nP),vLnAct(nP)
  LOGICAL :: vLPole(nP)
  REAL(dp):: G
  INTEGER :: i
  LOGICAL :: Ok
  CHARACTER(LEN=30):: Msg
  !---------------------------------------------------------------------
  vLPole(:)= (vX(:)>Zero)
  !
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
  !
  RETURN
END FUNCTION MixPhase_Grt

SUBROUTINE Mixture_Minimize( &
& vAffPole,vIPole,MixModel, &
& MixFas) !,nFmix)
!--
!-- for a given mixture model, MixModel,
!-- and given affinities of its end-members, vAffPole,
!-- compute the composition X that minimizes G
!-- and compute G at X
!--
  USE M_Numeric_Const, ONLY: Ln10
  USE M_Safe_Functions,ONLY: FSafe_Exp
  USE M_System_Vars,   ONLY: TdgK,Pbar
  !
  USE M_T_Phase,   ONLY: T_Phase
  USE M_T_MixPhase,ONLY: T_MixPhase
  USE M_T_MixModel,ONLY: T_MixModel,MaxPole
  !
  USE M_Optimsolver_Theriak
  USE M_MixModel_Optim
  !
  USE M_Global_Vars,ONLY: vFas
  !
  REAL(dp),        INTENT(IN)   :: vAffPole(:)
  INTEGER,         INTENT(IN)   :: vIPole(:)
  TYPE(T_MixModel),INTENT(IN)   :: MixModel
  TYPE(T_MixPhase),INTENT(INOUT):: MixFas
  !
  REAL(dp):: vX(MaxPole),vMu(MaxPole),vLnAct(MaxPole)
  REAL(dp):: vXmin(MaxPole),vMuMin(MaxPole)
  REAL(dp):: G,Gmin
  INTEGER :: nP,K,P
  INTEGER :: FF
  INTEGER :: Multi
  !
  REAL(dp):: TolX,DeltaInit
  INTEGER :: its,nCallG
  LOGICAL :: OkConverge

  FF= 0 ! log file
  TolX= MixMinim_TolX
  DeltaInit= 0.05D0

  nP= MixModel%NPole
  Multi= MixModel%vMulti(1)
  !
  !! !---- compute affin' of all reactions between end-memb' and water --
  !! DO P= 1,nP
  !!   ! index of pole P in pure phase list
  !!   ! K= tIPole(iModel,P)
  !!   K= MixModel%vIPole(P)
  !!   ! compute affinity of reaction (pole P <> water )
  !!   !vAffPole(P)= vFas(K)%Grt &
  !!   !&          - DOT_PRODUCT(tNuFas(K,:), vSpc(vOrdPr(:))%G0rt) &
  !!   !&          - DOT_PRODUCT(tNuFas(K,:), vLnAct(vOrdPr(:)))
  !!   vAffPole(P)= vDeltaGPole(K) &
  !!   &          - DOT_PRODUCT(tNuFas(K,:), vLnAct(vOrdPr(:)))
  !! ENDDO
  !! !---

  IF( TRIM(MixModel%Model)=="IDEAL" .AND. MixModel%NMarg==0 ) THEN

    !---------------------------------------------------analytic minimum

    Multi= MixModel%vMulti(1)
    Multi= 1

    DO P=1,nP
      ! K= MixModel%vIPole(P)
      ! vXmin(P)= FSafe_Exp(-vFasPur(MixModel%vIPole(P))%Grt /REAL(Multi))
      vXmin(P)= FSafe_Exp(-vAffPole(vIPole(P)) /REAL(Multi))
    ENDDO
    vXmin(1:nP)=  vXmin(1:nP) /SUM(vXmin(1:nP))
    vLnAct(1:nP)= Multi*LOG(vXmin(1:nP))
    vMuMin(1:nP)= vAffPole(vIPole(1:nP)) + vLnAct(1:nP)

    Gmin= SUM(vXmin(1:nP) * vMuMin(1:nP))

    !--------------------------------------------------------------trace
    IF(iDebug>2) THEN
    WRITE(74,'(A)') "--Mixture_Minimize--"
      DO P=1,nP
        WRITE(74,'(A,2G15.6,1X,A)') &
        & "vAffPole(P),vXmin(P)= ",&
        & vAffPole(vIPole(P)), &
        & vXmin(P), &
        & TRIM(vFas(vIPole(P))%NamFs)
      ENDDO
    ENDIF
    !-------------------------------------------------------------/trace

  ELSE

    !-----------------------------------------------numerical minimum(s)

    CALL MixModel_Optim_SetParams(TdgK,Pbar,MixModel)
    ALLOCATE(Mixmodel_Optim_vMu0rt(nP))
    ALLOCATE(Mixmodel_Optim_vLPole(nP))
    ! Mixmodel_Optim_vMu0rt(1:nP)= vFasPur(MixModel%vIPole(1:nP))%Grt
    Mixmodel_Optim_vMu0rt(1:nP)= vAffPole(vIPole(1:nP))
    Mixmodel_Optim_vLPole(1:nP)= .TRUE.

    Gmin= 1.0D30

    DO P=1,nP

      !--- start from compos'nP close to end-member P
      vX(P)= One - 1.0D-3
      DO K=1,nP
        IF(K/=P) vX(K)= 1.0D-3/REAL(nP-1)
      ENDDO

      !~ vX(1:nP)= vMixFas_Xpole_Init(I)%tXPole(P,1:nP)

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

      !~ IF(.NOT. OkConverge) THEN
        !~ PRINT *,"NoConvergence in Mixture_Minimize"
        !~ CYCLE
      !~ ENDIF
      !vMixFas_Xpole_Init(I)%tXPole(P,1:nP)= vX(1:nP)

      IF(G<Gmin) THEN
        Gmin= G
        vXmin(1:nP)= vX(1:nP)
        vMuMin(1:nP)= vMu(1:nP)
      ENDIF

    END DO
    !PAUSE

    DEALLOCATE(Mixmodel_Optim_vMu0rt)
    DEALLOCATE(Mixmodel_Optim_vLPole)

  ENDIF

  MixFas%vLPole(:)= .TRUE.
  write(76,'(G15.6)') SUM(vXmin(:))
  vXmin(:)= vXmin(:) /SUM(vXmin(:))
  MixFas%vXPole(:)= vXmin(:)
  MixFas%Grt= Gmin
  MixFas%vLnAct(1:nP)= vMuMin(1:nP) -vAffPole(vIPole(1:nP))

  !------------------------------------------------------------log files
  !! IF(F1>0) THEN
  !!   WRITE(F1,'(/,A,/)') "Mixture_Minimize,EndMember,Xi,Gi"
  !!   DO I=1,SIZE(vMixModel)
  !!     !IF(vMixFas0(I)%Grt > 1.D-3) CYCLE
  !!     MixModel= vMixModel(I)
  !!     WRITE(F1,'(2A)') "MODEL= ",MixModel%Name
  !!     DO K=1,MixModel%NPole
  !!       WRITE(F1,'(A,1X,G15.6)') &
  !!       & MixModel%vNamPole(K),        &
  !!       & vMixFas0(I)%vXPole(K)
  !!     ENDDO
  !!     WRITE(F1,'(A15,G15.6)') "G Mimim= ",vMixFas0(I)%Grt/Ln10
  !!     WRITE(F1,*)
  !!   ENDDO
  !!   WRITE(F1,*)
  !! ENDIF
  !! !
  !! IF(F2>0) THEN
  !!   DO I=1,SIZE(vMixModel)
  !!   !IF(vMixFas0(I)%Grt < Zero) THEN
  !!     MixModel= vMixModel(I)
  !!     !
  !!     WRITE(F2,'(A,A1)',ADVANCE="NO") MixModel%Name,T_
  !!     DO K=1,MixModel%NPole
  !!       WRITE(F2,'(G15.6,A1)',ADVANCE="NO") vMixFas0(I)%vXPole(K),T_
  !!     ENDDO
  !!     !DO K=1,MixModel%NPole
  !!     !  WRITE(F2,'(G15.6,A1)',ADVANCE="NO") vFasPur(MixModel%vIPole(K))%Grt,T_
  !!     !ENDDO
  !!     WRITE(F2,'(G15.6,A1)',ADVANCE="NO") vMixFas0(I)%Grt,T_
  !!     !
  !!   !ENDIF
  !!   ENDDO
  !!   WRITE(F2,*)
  !! ENDIF
  !-----------------------------------------------------------/log files
  !
END SUBROUTINE Mixture_Minimize

ENDMODULE M_Equil_2
