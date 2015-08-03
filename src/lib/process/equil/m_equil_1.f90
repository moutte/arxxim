MODULE M_Equil_1
!--
!-- calculate equilibrium speciation --
!-- considering saturation with minerals or gases --
!-- this version considers pure phases and mixtures --
!--
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_,Pause_
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Equil_Eq1
  !
  REAL(dp),PARAMETER:: MixMinim_TolX= 1.D-4
  != parameter for Mixture_Minimize
  !
  REAL(dp),PARAMETER:: FasMinim= 1.D-9
  ! FasMinim= mimim' amount of phase
  !
  REAL(dp),PARAMETER:: AffinIota= 1.D-6
  ! AffinIota= tolerance for equilibrium condition
  ! ( equilibrium == ABS(Affinity)<AffinIota )
  !
  REAL(dp),PARAMETER:: dXi_Minim= 1.D-16
  REAL(dp),PARAMETER:: dXi_Maxim= 1.D+12
  !
  INTEGER, ALLOCATABLE:: tIPole(:,:)
  CHARACTER(LEN=30):: sFMT

CONTAINS

SUBROUTINE Equil_Eq1(NoMixture,iErr)
!--
!-------------------- this version considers pure phases and mixtures --
!--
  USE M_IOTools,      ONLY: GetUnit
  USE M_Files,        ONLY: DirOut
  USE M_Numeric_Const,ONLY: Ln10
  USE M_Numeric_Tools,ONLY: iMinLoc_R
  USE M_Dtb_Const,    ONLY: T_CK
  USE M_IoTools,      ONLY: OutStrVec
  USE M_T_MixModel,   ONLY: MaxPole
  USE M_T_MixPhase,   ONLY: T_MixPhase
  USE M_T_Phase,      ONLY: T_Phase
  !
  USE M_Equil_Solve,  ONLY: Equil_Solve
  USE M_Equil_Tools,  ONLY: Equil_FasTrace_EnTete,Equil_FasTrace_Write
  !
  USE M_Global_Vars,  ONLY: vSpc,vFas
  USE M_Global_Vars,  ONLY: vMixModel
  USE M_Basis_Vars,   ONLY: nCi,nAs,isW
  USE M_Basis_Vars,   ONLY: tNuFas,vOrdPr
  USE M_Equil_Vars,   ONLY: cEquMode,dXi
  USE M_Equil_Vars,   ONLY: vYesList
  USE M_Equil_Vars,   ONLY: vDeltaG_Fs,vAffScale,vLnAct,vFasMole
  USE M_Equil_Vars,   ONLY: vEquFas,nEquFas,T_EquPhase
  USE M_Equil_Vars,   ONLY: iCountEq,fTrcEq
  !
  LOGICAL,INTENT(IN) :: NoMixture
  INTEGER,INTENT(OUT):: iErr
  !
  INTEGER :: nCp,nFs,nPur,nSol
  INTEGER :: iFs,iMix,iPur
  INTEGER :: FF
  INTEGER :: nCheck,I,J,K
  INTEGER :: nMix,nP
  INTEGER :: iPrecip,iElimin
  INTEGER :: nIts,iDo0,iDo1
  ! REAL(dp):: rDum= One
  ! INTEGER :: vIPole(MaxPole)
  LOGICAL :: dXiTooLow
  LOGICAL :: Use_nIts= .TRUE.  ! .FALSE. !

  LOGICAL, ALLOCATABLE:: vYesFas(:) ! phase is in current equil.phase list
  REAL(dp),ALLOCATABLE:: vFasAff0(:)

  TYPE(T_EquPhase),ALLOCATABLE:: vEquFas0(:)
  TYPE(T_MixPhase),ALLOCATABLE:: vMixFas0(:)
  TYPE(T_Phase),   ALLOCATABLE:: vFas0(:)
  TYPE(T_EquPhase):: EquFasNew
  ! TYPE(T_MixPhase):: MixFas

  !~ TYPE:: T_EquPhase
    !~ CHARACTER(LEN=23):: NamEq
    !~ INTEGER :: iSpc
    !~ INTEGER :: iMix
    !~ REAL(dp):: vXPole(1:MaxPole)
    !~ !REAL(dp):: Grt
    !~ REAL(dp):: Mole
  !~ END TYPE T_EquPhase

  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Eq1_Mix"

  FF=0
  IF(iDebug>1 .AND. (FF==0)) THEN
    !~ iCountEq=0
    !CALL Equil_FasTrace_EnTete(vFas,vYesList,FF)
    CALL GetUnit(FF)
    OPEN(FF,FILE=TRIM(DirOut)//"_equil.log")
  ENDIF

  cEquMode="EQ1" !-> flag in Equil_Solve
  dXi=       One !1.D3 !1.D-6 !reaction increment
  dXiTooLow= .FALSE.

  nCheck= 0

  nCp= SIZE(vOrdPr)
  !! nPur= COUNT(vFas(:)%Typ=="PURE")
  nPur= COUNT(vFas(:)%iSpc/=0)
  nMix= SIZE(vMixModel)
  nSol= 1 !! SIZE(vSolFas)

  IF(NoMixture) nMix= 0

  nFs= nPur + nMix + nSol

  ALLOCATE(vYesFas(nFs)) ; vYesFas(:)= .FALSE.
  ALLOCATE(vFasAff0(nFs))
  ALLOCATE(vEquFas0(nCp))

  IF(nMix>0) THEN

    ALLOCATE(vFas0(nPur))
    vFas0(:)= vFas(1:nPur)
    DEALLOCATE(vFas)
    ALLOCATE(vFas(nFs))
    vFas(1:nPur)= vFas0(1:nPur)
    DO K=1,nMix
      vFas(nPur+K)%NamFs= vMixModel(K)%Name
    ENDDO
    DEALLOCATE(vFas0)

    ALLOCATE(vMixFas0(nMix))
    DO K=1,nMix
      vMixFas0(K)%Name= vMixModel(K)%Name
      vMixFas0(K)%iModel= K
      vMixFas0(K)%vXPole(:)= Zero
    ENDDO
    !
    ALLOCATE(tIPole(nMix,MaxPole))
    CALL MixModel_Init(nPur,vFas,vMixModel)

    !~ CALL MixPhase_XPole_Init(vMixModel,vMixFas_Xpole_Init)
    !
    !~ ALLOCATE(vDeltaGPole(nPur))
    !~ DO iFs= 1,nPur
      !~ vDeltaGPole(iFs)= vFas(iFs)%Grt &
      !~ &               - DOT_PRODUCT(tNuFas(iFs,:),vSpc(vOrdPr(:))%G0rt)
    !~ ENDDO

  ENDIF
  !
  I=0
  DO iPur=1,nPur
    IF(vFas(iPur)%Mole>FasMinim) THEN
      vYesFas(iPur)= .TRUE.
      I=I+1
      vEquFas(I)%iPur= iPur
      vEquFas(I)%NamEq= TRIM(vFas(iPur)%NamFs)
      vEquFas(I)%iMix= 0
      vEquFas(I)%Mole= vFas(iPur)%Mole
    ENDIF
  ENDDO
  nEquFas= I
  !
  !-------------------------------------- loop on mineral saturation, --
  !------------------------ EXIT when all minerals have vPurYes false --
  iDo0= 0
  Do0: DO
    iDo0= iDo0 +1

    iDo1= 0
    !---------------------- solve the system until no negative phases --
    Do1: DO
      iDo1= iDo1 +1

      !~ IF(iDebug>2) THEN
        !~ WRITE(6,'(A)') "<-- CURRENT --"
        !~ DO I=1,nEquFas
          !~ WRITE(6,'(2X,A)') vEquFas(I)%NamEq
        !~ ENDDO
        !~ WRITE(6,'(A)') "</- CURRENT --"
      !~ ENDIF

      IF(nEquFas>0) THEN
        CALL EquFas_Clean
        CALL EquFas_Alloc(nCp,nEquFas)
        CALL EquFas_Update(nCp,nEquFas,vEquFas)
      ENDIF

      !-------------------------------------------------- loop on dXi --
      Do2: DO

        !------------------------------------------- solve the system --
        CALL Equil_Solve(nIts,iErr)
        !------------------------------------------/ solve the system --

        IF(iErr==0) EXIT Do2 !when system solved-> EXIT

        !IF necessary, try smaller dXi
        dXi= dXi /2.0_dp  !ELSE DO again with smaller step
        IF(dXi<dXi_Minim) dXiTooLow=.TRUE.
        IF(dXiTooLow) EXIT Do0 !-> zannennagara,  found no solution ...

      ENDDO Do2
      !--------------------------------------------------/loop on dXi --
      !
      IF(nEquFas==0) EXIT Do1
      !
      IF(Use_Nits) THEN
        IF(nIts>12 .AND. dXi<dXi_Minim) dXi= dXi /2.0_dp
        IF(nIts<6  .AND. dXi<dXi_Maxim) dXi= dXi *2.0_dp
      ENDIF
      !
      IF(iDebug>3) CALL Trace_1
      !
      !----------------------------------- exclude "near zero" phases --
      vEquFas0= vEquFas
      J= 0
      DO I=1,nEquFas
        IF(ABS(vEquFas0(I)%Mole)>FasMinim) THEN
          J=J+1
          vEquFas(J)= vEquFas0(I)
          IF(vEquFas0(I)%iPur>0) vYesFas(vEquFas0(I)%iPur)= .TRUE.
          IF(vEquFas0(I)%iMix>0) vYesFas(nPur+vEquFas0(I)%iPur)= .TRUE.
        ENDIF
      ENDDO
      nEquFas= J

      vYesFas(:)= .FALSE.
      vYesFas(vEquFas(1:nEquFas)%iPur)= .TRUE.
      !----------------------------------/ exclude "near zero" phases --
      !
      !------------------------ check for phases with negative amount --
      IF(nEquFas>0) THEN
        !
        !~ !------------------------------------------------------ trace --
        !~ IF(iDebug>2) THEN
          !~ WRITE(6,'(A)') "<-- CURRENT --"
          !~ DO I=1,nEquFas
            !~ WRITE(6,'(G15.6,1X,A15,1X)',ADVANCE="NO") &
            !~ & vEquFas(I)%Mole,TRIM(vEquFas(I)%NamEq)
            !~ IF(vEquFas(I)%iMix>0) THEN
              !~ DO J=1,vEquFas(I)%nPole
                !~ WRITE(6,'(G15.6,1X)',ADVANCE="NO") vEquFas(I)%vXPole(J)
              !~ ENDDO
            !~ ENDIF
            !~ WRITE(6,*)
          !~ ENDDO
          !~ WRITE(6,'(A)') "</- CURRENT --"
        !~ ENDIF
        !~ !------------------------------------------------------/trace --
        !
        IF(MINVAL(vEquFas(1:nEquFas)%Mole)<-FasMinim) THEN
          iElimin= iMinLoc_R(vEquFas(1:nEquFas)%Mole)
          IF(iDebug>2) THEN
            DO I=1,nEquFas
              WRITE(6,'(2X,A)') vEquFas(I)%NamEq
            ENDDO
            WRITE(6,'(16X,2A)') "NEG= ", TRIM(vEquFas(iElimin)%NamEq)
          ENDIF
          CALL EquFas_Remove(iElimin)

          IF(vEquFas(iElimin)%iPur>0) &
          & vYesFas(vEquFas(iElimin)%iPur)= .FALSE.
          IF(vEquFas(iElimin)%iMix>0) &
          & vYesFas(nPur+vEquFas(iElimin)%iMix)= .FALSE.

          CYCLE Do1 !========< remove most "negative" phase and CYCLE ==
        ENDIF
        !
      ENDIF
      !
      EXIT Do1
      !-----------------------/ check for phases with negative amount --
      !
    ENDDO Do1
    !---------------------/ solve the system until no negative phases --
    !
    IF(iErr<0) EXIT Do0 != error in Newton -> EXIT
    !
    vFasAff0= Zero
    DO iPur=1,nPur
      IF(vYesList(iPur)) &
      & vFasAff0(iPur)= &
      & (  vDeltaG_Fs(iPur) &
      &  - DOT_PRODUCT(tNuFas(iPur,1:nCp),vLnAct(vOrdPr(1:nCp))) ) &
      & /vAffScale(iPur) !scaling ala EQ3/6
    ENDDO
    !
    IF(nMix>0) THEN
      !!CALL Mixture_Minimize(vFasAff0,vMixModel,vMixFas0)
      DO iMix=1,nMix
        !vIPole(:)= tIPole(iMix,:)
        CALL Mixture_Minimize( &
        & vFasAff0,tIPole(iMix,:),vMixModel(iMix), &
        & vMixFas0(iMix))
      ENDDO
      vFasAff0(nPur+1:nPur+nMix)= vMixFas0(1:nMix)%Grt
    ENDIF
    !
    IF(iDebug>3) CALL Trace_2
    !
    !---------------------------------------------------------- trace --
    !~ IF(iDebug>2) THEN
      !~ WRITE(6,'(A)') "<-- OVERSAT --"
      !~ DO I=1,SIZE(vFasAff0)
        !~ IF (vFasAff0(I)<AffinIota) THEN
          !~ IF(I<nPur) THEN
            !~ WRITE(6,'(G15.6,1X,A15,1X)') vFasAff0(I),vFas(I)%NamFs
          !~ ELSE
            !~ WRITE(6,'(G15.6,1X,A15,1X)') vFasAff0(I),vFas(I)%NamFs
          !~ ENDIF
        !~ ENDIF
      !~ ENDDO
      !~ WRITE(6,'(A)') "</- OVERSAT --"
    !~ ENDIF
    !----------------------------------------------------------/trace --
    !
    !----------------- add phase with highest supersaturation, if any --
    IF(MINVAL(vFasAff0)<-AffinIota) THEN
      iPrecip= iMinLoc_R(vFasAff0)
      !-> the saturated phase with lowest vFasAff
      IF(vYesFas(iPrecip)) iPrecip= 0
    ELSE
      iPrecip= 0 != there is no saturated phase ==
      EXIT Do0 ! no phase found with affinity below -AffinIota > EXIT ==
    ENDIF

    IF(iPrecip>0) THEN

      IF(iPrecip<=nPur) THEN

        vYesFas(iPrecip)= .TRUE.

        EquFasNew%iPur= iPrecip
        EquFasNew%NamEq= TRIM(vFas(iPrecip)%NamFs)
        EquFasNew%iMix= 0
        EquFasNew%Mole= Zero

      ELSE

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

      !~ CALL CheckAssemblage(FF,nEquFas,EquFasNew,vEquFas,iElimin)
      !
      !~ IF(iElimin>0) THEN
        !~ iPur= vEquFas(iElimin)%iPur
        !~ iMix= vEquFas(iElimin)%iMix
        !~ IF(iPur > 0) vFas(iPur)%Mole= Zero
        !~ IF(iMix > 0) vFas(nPur+iMix)%Mole= Zero
        !~ IF(iDebug>2) &
        !~ & WRITE(6,'(16X,2A)') "OUT= ",TRIM(vEquFas(iElimin)%NamEq)
        !~ CALL EquFas_Remove(iElimin)
      !~ ENDIF

      IF(iDebug>2) THEN
        DO I=1,nEquFas
          WRITE(6,'(2X,A)') vEquFas(I)%NamEq
        ENDDO
        WRITE(6,'(16X,2A)') "ADD= ",TRIM(EquFasNew%NamEq)
      ENDIF

      IF(nEquFas < SIZE(vEquFas)) THEN
        vEquFas(nEquFas+1)= EquFasNew
        nEquFas= nEquFas +1
      ENDIF

    ENDIF

    IF(iDebug>3) THEN !=======================================< trace ==
      DO I=1,nEquFas
        iFs= vEquFas(I)%iPur
        IF(iFs>0) THEN
          WRITE(6,'(A,I4,1X,A15,2(A,G15.6))') &
          & "Do0",iDo0, &
          & vFas(iFs)%NamFs, &
          & "/ lQsK=", -vFasAff0(iFs)/Ln10 ,&
          & "/ Mole=",vEquFas(I)%Mole
        ENDIF
      ENDDO
      IF(iDebug==4) CALL Pause_
    ENDIF !==================================================</ trace ==

    !----------------/ add phase with highest supersaturation, if any --

    IF(nEquFas==0) EXIT Do0

    !~ IF(iDebug>2) PAUSE
  ENDDO Do0
  !-------------------------------------/ loop on mineral saturation, --
  !
  !--- save mineral mole numbers in vFas
  !~ CALL Save_EquPhase
  !
  DEALLOCATE(vFasMole)
  ALLOCATE(vFasMole(nPur+nMix))
  vFasMole(:)= Zero
  DO I=1,nEquFas
    IF(vEquFas(I)%iPur>0) vFasMole(vEquFas(I)%iPur)= vEquFas(I)%Mole
    IF(vEquFas(I)%iMix>0) vFasMole(vEquFas(I)%iMix+nPur)= vEquFas(I)%Mole
  ENDDO
  !~ pause
  !---/
  !
  CALL EquFas_Clean
  !
  !~ DEALLOCATE(vNewFas)
  DEALLOCATE(vFasAff0)
  DEALLOCATE(vEquFas0)
  IF(nMix>0) THEN
    DEALLOCATE(vMixFas0)
    DEALLOCATE(tIPole)
  ENDIF
  !
  IF(FF>0) CLOSE(FF)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Eq1_Mix"
  !
CONTAINS

SUBROUTINE EquFas_Remove(N)
  INTEGER,INTENT(IN):: N
  INTEGER:: I,J
  vEquFas0= vEquFas
  J= 0
  DO I=1,nEquFas
    IF(I/=N) THEN
      J=J+1
      vEquFas(J)= vEquFas0(I)
    ENDIF
  ENDDO
  nEquFas= J
END SUBROUTINE EquFas_Remove

SUBROUTINE Trace_1
  WRITE(6,'(A)') "_"
  DO I=1,nEquFas
    IF(vEquFas(I)%iPur>0) &
    WRITE(6,'(A,2(A,G15.6))') &
    & vFas(vEquFas(I)%iPur)%NamFs, &
    & "/ lQsK=",-vFasAff0(vEquFas(I)%iPur)/Ln10, &
    & "/ Mole=",vEquFas(I)%Mole
  ENDDO
ENDSUBROUTINE Trace_1

SUBROUTINE Trace_2
  WRITE(6,'(6X,A)') "Do0"
  DO I=1,nEquFas
    IF(vEquFas(I)%iPur>0) &
    WRITE(6,'(A,2(A,G15.6))') &
    & vFas(vEquFas(I)%iPur)%NamFs, &
    & "/ lQsK=",-vFasAff0(vEquFas(I)%iPur)/Ln10, &
    & "/ Mole=",vEquFas(I)%Mole
  ENDDO
ENDSUBROUTINE Trace_2

ENDSUBROUTINE Equil_Eq1

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

    ! find the indices of end-members of mixture vMixFas(iMix)
    ! in the pure phase list vFas(1:nPur)
    ! -> indices in tIPole(iMix,:)
    DO P=1,nP
      P1= vMixModel(iModel)%vIPole(P)
      !-> index in vSpc
      !-> find the pure phases that points to this species
      P2= 0
      DO iPur=1,nPur
        IF(vFas(iPur)%iSpc==P1) THEN
          P2= iPur
          EXIT
        ENDIF
      ENDDO
      tIPole(iModel,P)= P2
      !
    ENDDO

  ENDDO

  RETURN
ENDSUBROUTINE MixModel_Init

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
  USE M_T_MixModel, ONLY: MaxPole
  USE M_Global_Vars,ONLY: vMixModel
  USE M_Basis_Vars, ONLY: tAlfFs,tNuFas
  USE M_Equil_Vars, ONLY: vDeltaG_Eq,tAlfEq,tNuEq
  USE M_Equil_Vars, ONLY: vDeltaG_Fs,T_EquPhase
  !
  INTEGER,INTENT(IN):: nCp
  INTEGER,INTENT(IN):: nEquFas
  TYPE(T_EquPhase),INTENT(IN):: vEquFas(:)
  INTEGER:: iPur,iMix,I,C
  INTEGER:: nP
  INTEGER:: vIPole(MaxPole)
  !
  DO I=1,nEquFas
    !
    iPur= vEquFas(I)%iPur
    iMix= vEquFas(I)%iMix
    !
    IF(iPur>0) THEN
      vDeltaG_Eq(I)= vDeltaG_Fs(iPur)
      tAlfEq(:,I)= tAlfFs(:,iPur)
      tNuEq(I,:)= tNuFas(iPur,:)
    ENDIF
    !
    IF(iMix>0) THEN
      nP= vEquFas(I)%NPole
      vIPole(1:nP)= vEquFas(I)%vIPole(1:nP)
      vDeltaG_Eq(I)= SUM( vEquFas(I)%vXPole(1:nP) &
      &                 * vDeltaG_Fs(vIPole(1:nP)) ) &
      &            + SUM( vEquFas(I)%vXPole(1:nP) &
      &                 * vEquFas(I)%vLnAct(1:nP) )
      DO C=1,nCp
        tAlfEq(C,I)= SUM( vEquFas(I)%vXPole(1:nP) &
        &               * tAlfFs(C,vIPole(1:nP)) )
        tNuEq(I,C)= SUM( vEquFas(I)%vXPole(1:nP) &
        &              * tNuFas(vIPole(1:nP),C) )
      ENDDO
    ENDIF
    !
  ENDDO
  !
  WRITE(sFMT,'(a,i3,a)') '(A,1X,',nCp,'(G12.3,1X))'
  DO I=1,nEquFas
    WRITE(73,'(A,G15.6,1X,A)') "DeltaG_Eq= ",vDeltaG_Eq(I),TRIM(vEquFas(I)%NamEq)
    WRITE(73,sFMT)        "tAlfEq=    ",(tAlfEq(C,I),C=1,nCp)
    WRITE(73,sFMT)        "tNuEq=     ",(tNuEq(I,C),C=1,nCp)
  ENDDO
  WRITE(73,'(A)') "===================================================="
  !
  RETURN
ENDSUBROUTINE EquFas_Update

SUBROUTINE EquFasMix_Update(nCp,I,EquFas)
  USE M_T_MixModel, ONLY: MaxPole
  USE M_Global_Vars,ONLY: vMixModel
  USE M_Basis_Vars, ONLY: tAlfFs,tNuFas
  USE M_Equil_Vars, ONLY: vDeltaG_Eq,tAlfEq,tNuEq
  USE M_Equil_Vars, ONLY: vDeltaG_Fs,T_EquPhase
  !
  INTEGER,INTENT(IN):: nCp
  INTEGER,INTENT(IN):: I
  TYPE(T_EquPhase),INTENT(IN):: EquFas
  !
  INTEGER:: C
  INTEGER:: nP
  INTEGER:: vIPole(MaxPole)
  !
  nP= EquFas%NPole
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
  !
  RETURN
ENDSUBROUTINE EquFasMix_Update

SUBROUTINE EquFas_Clean
  USE M_Equil_Vars,ONLY: vDeltaG_Eq,tAlfEq,tNuEq
  
  IF(ALLOCATED(vDeltaG_Eq)) DEALLOCATE(vDeltaG_Eq)
  IF(ALLOCATED(tAlfEq))     DEALLOCATE(tAlfEq)
  IF(ALLOCATED(tNuEq))      DEALLOCATE(tNuEq)
  
  RETURN
ENDSUBROUTINE EquFas_Clean

SUBROUTINE Save_EquPhase
  USE M_Equil_Vars, ONLY: T_EquPhase
  USE M_Global_Vars,ONLY: vMixModel,vMixFas,vFas
  USE M_Equil_Vars, ONLY: nEquFas,vEquFas
  !
  !~ INTEGER,INTENT(IN):: nEquFas
  !~ TYPE(T_EquPhase),INTENT(IN):: vEquFas(:)
  !
  INTEGER:: nPur,I,K,nP
  
  nPur= COUNT(vFas(:)%iSpc>0)
  !
  vFas(:)%Mole= Zero
  DO I=1,nEquFas
    IF(vEquFas(I)%iPur>0) vFas(vEquFas(I)%iPur)%Mole= vEquFas(I)%Mole
    IF(vEquFas(I)%iMix>0) vFas(vEquFas(I)%iMix+nPur)%Mole= vEquFas(I)%Mole
  ENDDO
  !
  K=0
  DO I=1,nEquFas
    IF(vEquFas(I)%iMix>0) K= K+1
  ENDDO
  DEALLOCATE(vMixFas)
  ALLOCATE(vMixFas(K))
  !
  K=0
  IF(SIZE(vMixFas)>0) THEN
    DO I=1,nEquFas
      IF(vEquFas(I)%iMix>0) THEN
        K= K+1
        !
        nP= vMixModel(vEquFas(I)%iMix)%nPole
        vMixFas(K)%Name= TRIM(vEquFas(I)%NamEq)
        vMixFas(K)%iModel= vEquFas(I)%iMix
        vMixFas(K)%vXPole(1:nP)= vEquFas(I)%vXPole(1:nP)
        !
        vFas(nPur+K)%NamFs= TRIM(vEquFas(I)%NamEq)
        !! vFas(nPur+K)%Typ=   "MIXT"
        vFas(nPur+K)%iSpc=  0
        vFas(nPur+K)%iSol=  0
        vFas(nPur+K)%iMix=  K
        !
      ENDIF
    ENDDO
  ENDIF
  
  RETURN
ENDSUBROUTINE Save_EquPhase

SUBROUTINE Mixture_Minimize( &
& vAffPole,vIPole,MixModel, &
& MixFas) !,nFmix)
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
  TYPE(T_MixModel):: MM
  REAL(dp):: vX(MaxPole),vMu(MaxPole),vLnAct(MaxPole)
  REAL(dp):: vXmin(MaxPole),vMuMin(MaxPole)
  REAL(dp):: G,Gmin
  INTEGER :: nP,K,P
  INTEGER :: FF
  INTEGER :: Multi
  REAL(dp):: S
  !
  REAL(dp):: TolX,DeltaInit
  INTEGER :: its,nCallG
  LOGICAL :: OkConverge
  !
  FF= 0 ! log file
  TolX= MixMinim_TolX
  DeltaInit= 0.05D0
  !
  !! DO I=1,SIZE(vMixFas0)
    !
    !! iModel= vMixFas0(I)%iModel
    !! iModel= MixFas%iModel
    !! MM= vMixModel(iModel)
    MM= MixModel
    nP= MM%NPole
    !! vIPole(1:nP)= tIPole(iModel,1:nP)
    Multi= MM%vMulti(1)
    !
    !~ !---- compute affin' of all reactions between end-memb' and water --
    !~ DO P= 1,nP
      !~ ! index of pole P in pure phase list
      !~ ! K= tIPole(iModel,P)
      !~ K= MM%vIPole(P)
      !~ ! compute affinity of reaction (pole P <> water )
      !~ !vAffPole(P)= vFas(K)%Grt &
      !~ !&          - DOT_PRODUCT(tNuFas(K,:), vSpc(vOrdPr(:))%G0rt) &
      !~ !&          - DOT_PRODUCT(tNuFas(K,:), vLnAct(vOrdPr(:)))
      !~ vAffPole(P)= vDeltaGPole(K) &
      !~ &          - DOT_PRODUCT(tNuFas(K,:), vLnAct(vOrdPr(:)))
    !~ ENDDO
    !~ !---
    !
    IF( TRIM(MM%Model)=="IDEAL" .AND. MM%NMarg==0 ) THEN
      !--------------------------------------------- analytic minimum --
      !
      Multi=  MM%vMulti(1)
      !
      S= Zero
      DO P=1,nP
        ! K= MM%vIPole(P)
        ! vXmin(P)= FSafe_Exp(-vFasPur(MM%vIPole(P))%Grt /REAL(Multi))
        vXmin(P)= FSafe_Exp(-vAffPole(vIPole(P)) /REAL(Multi))
        S= S + vXmin(P)
      ENDDO
      vXmin(1:nP)=  vXmin(1:nP) /S
      vLnAct(1:nP)= Multi*LOG(vXmin(1:nP))
      vMuMin(1:nP)= vAffPole(vIPole(1:nP)) + vLnAct(1:nP)
      DO P=1,nP
        WRITE(74,'(A,2G15.6,1X,A)') &
        & "vAffPole(P),vXmin(P)= ",&
        & vAffPole(vIPole(P)), &
        & vXmin(P), &
        & TRIM(vFas(vIPole(P))%NamFs)
      ENDDO
      !pause
      Gmin= SUM(vXmin(1:nP) * vMuMin(1:nP))
      !
    ELSE
      !----------------------------------------- numerical minimum(s) --
      !
      CALL MixModel_Optim_SetParams(TdgK,Pbar,MM)
      ALLOCATE(Mixmodel_Optim_vMu0rt(nP))
      ALLOCATE(Mixmodel_Optim_vLPole(nP))
      ! Mixmodel_Optim_vMu0rt(1:nP)= vFasPur(MM%vIPole(1:nP))%Grt
      Mixmodel_Optim_vMu0rt(1:nP)= vAffPole(vIPole(1:nP))
      Mixmodel_Optim_vLPole(1:nP)= .TRUE.
      !
      Gmin= 1.0D30
      !
      DO P=1,nP
        !
        !--- start from compos'nP close to end-member P
        vX(P)= One - 1.0D-3
        DO K=1,nP
          IF(K/=P) vX(K)= 1.0D-3/REAL(nP-1)
        ENDDO
        !
        !~ vX(1:nP)= vMixFas_Xpole_Init(I)%tXPole(P,1:nP)
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
        !vMixFas_Xpole_Init(I)%tXPole(P,1:nP)= vX(1:nP)
        !
        IF(G<Gmin) THEN
          Gmin= G
          vXmin(1:nP)= vX(1:nP)
          vMuMin(1:nP)= vMu(1:nP)
        ENDIF
        !
      END DO
      !PAUSE
      !
      DEALLOCATE(Mixmodel_Optim_vMu0rt)
      DEALLOCATE(Mixmodel_Optim_vLPole)
      !
    ENDIF
    !
    !!vMixFas0(i)%vLPole(:)= .TRUE.
    !!vMixFas0(i)%vXPole(:)= vXmin(:)
    !!vMixFas0(i)%Grt= Gmin
    !!vMixFas0(i)%vLnAct(1:nP)= vMuMin(1:nP) -vAffPole(vIPole(1:nP))
    !
    MixFas%vLPole(:)= .TRUE.
    MixFas%vXPole(:)= vXmin(:)
    MixFas%Grt= Gmin
    MixFas%vLnAct(1:nP)= vMuMin(1:nP) -vAffPole(vIPole(1:nP))
    !
  !! ENDDO
  !
  !~ !-------------------------------------------------------- log files --
  !~ IF(F1>0) THEN
    !~ WRITE(F1,'(/,A,/)') "Mixture_Minimize,EndMember,Xi,Gi"
    !~ DO I=1,SIZE(vMixModel)
      !~ !IF(vMixFas0(I)%Grt > 1.D-3) CYCLE
      !~ MM= vMixModel(I)
      !~ WRITE(F1,'(2A)') "MODEL= ",MM%Name
      !~ DO K=1,MM%NPole
        !~ WRITE(F1,'(A,1X,G15.6)') &
        !~ & MM%vNamPole(K),        &
        !~ & vMixFas0(I)%vXPole(K)
      !~ ENDDO
      !~ WRITE(F1,'(A15,G15.6)') "G Mimim= ",vMixFas0(I)%Grt/Ln10
      !~ WRITE(F1,*)
    !~ ENDDO
    !~ WRITE(F1,*)
  !~ ENDIF
  !~ !
  !~ IF(F2>0) THEN
    !~ DO I=1,SIZE(vMixModel)
    !~ !IF(vMixFas0(I)%Grt < Zero) THEN
      !~ MM= vMixModel(I)
      !~ !
      !~ WRITE(F2,'(A,A1)',ADVANCE="NO") MM%Name,T_
      !~ DO K=1,MM%NPole
        !~ WRITE(F2,'(G15.6,A1)',ADVANCE="NO") vMixFas0(I)%vXPole(K),T_
      !~ ENDDO
      !~ !DO K=1,MM%NPole
      !~ !  WRITE(F2,'(G15.6,A1)',ADVANCE="NO") vFasPur(MM%vIPole(K))%Grt,T_
      !~ !ENDDO
      !~ WRITE(F2,'(G15.6,A1)',ADVANCE="NO") vMixFas0(I)%Grt,T_
      !~ !
    !~ !ENDIF
    !~ ENDDO
    !~ WRITE(F2,*)
  !~ ENDIF
  !~ !-------------------------------------------------------/ log files --
  !
END SUBROUTINE Mixture_Minimize

ENDMODULE M_Equil_1

