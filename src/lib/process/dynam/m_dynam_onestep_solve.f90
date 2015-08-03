MODULE M_Dynam_OneStep_Solve
  USE M_Kinds
  USE M_Trace,ONLY: iDebug,fTrc,T_
  IMPLICIT NONE

  PRIVATE

  PUBLIC:: Dynam_OneStep_Solve

  REAL(dp),DIMENSION(:),ALLOCATABLE:: vTolF

CONTAINS

!-----------------------------------------------------------------------
!-- solve the dynamic system for one time step--------------------------
!-----------------------------------------------------------------------
SUBROUTINE Dynam_OneStep_Solve( & !
& F,iStep,                      & !IN
& vLKinActiv,vLEquActiv,        & !IN
& vMolM_Slope,                  & !IN
& NewtTolF,NewtTolF_Equil,NewtTolX,NewtTolMin, & !IN
& Time_Decrease,NewtMaxIts,     & !IN
& Time,dTmin,dTime_in,          & !IN
!! & vMolF,vMolK,                  & !INOUT
& dTime_out,Newt_iDo,Newt_iErr, & !OUT
& DTimeTooSmall,VarMax)           !OUT
  !
  USE M_Numeric_Tools,ONLY: iMaxLoc_R,iMinLoc_R
  USE M_IoTools
  !
  USE M_System_Vars,ONLY: TdgK,Pbar
  USE M_Global_Vars,ONLY: vSpc,nAq,vKinFas
  USE M_Basis_Vars, ONLY: vOrdAq,nCi,nAs
  !
  USE M_Dynam_OneStep_Tools,ONLY: Dynam_OneStep_Secondspecies
  USE M_Dynam_OneStep_Tools,ONLY: Dynam_OneStep_KinFasSurf_Update
  USE M_Dynam_OneStep_Tools,ONLY: Dynam_OneStep_DLGammaUpd
  USE M_Dynam_OneStep_Tools,ONLY: Dynam_OneStep_GammaUpd
  !
  USE M_Dynam_Vars, ONLY: &
  & TimeScale, &
  & Implicit_Surface, &
  & LogForAqu,LogForMin,cMethod, &
  & DirectSub,vMolSec,  &
  & VBox,PhiF, &
  & vMolF,vMolK,&
  & vKinMinim, &
  & vLnAct,vLnGam,pH_,vDLGam_As, &
  & dTime, &
  & Extrapole
  !
  USE M_Dynam_Vars, ONLY: PhiF0,vStatusK,vSurfK,vSurfK0,vMolK0
  !
  INTEGER, INTENT(IN) :: F,iStep
  LOGICAL, INTENT(IN) :: vLKinActiv(:),vLEquActiv(:)
  REAL(dp),INTENT(IN) :: vMolM_Slope(:)
  REAL(dp),INTENT(IN) :: NewtTolF,NewtTolF_Equil
  REAL(dp),INTENT(IN) :: NewtTolX,NewtTolMin
  REAL(dp),INTENT(IN) :: Time_Decrease
  INTEGER, INTENT(IN) :: NewtMaxIts
  REAL(dp),INTENT(IN) :: Time,dTmin
  REAL(dp),INTENT(IN) :: dTime_in
  !
  REAL(dp),INTENT(OUT):: dTime_out
  INTEGER, INTENT(OUT):: Newt_iDo,Newt_iErr
  LOGICAL, INTENT(OUT):: DTimeTooSmall
  REAL(dp),INTENT(OUT):: VarMax
  !---------------------------------------------------------------------
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vX,vX0,vDelta,vMinim
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vMolF1,vMolK1
  REAL(dp),DIMENSION(:),ALLOCATABLE:: vLnGamOld
  LOGICAL,ALLOCATABLE:: vXisPlus(:)
  CHARACTER(LEN=7),ALLOCATABLE:: vName(:)
  !
  INTEGER :: iDo1,iDo2,iMaxDelta,iMk,I
  INTEGER :: nF,nMkA,nEqA,nDim
  REAL(dp):: NewtErrF,NewtErrX,NewtErrG
  REAL(dp):: dTime_Min
  REAL(dp):: MolK
  LOGICAL :: Newt_Check
  REAL(dp):: Z_Plus,Z_Minus
  !
  REAL(dp),PARAMETER:: DynTolGam= 1.D-3
  !---------------------------------------------------------------------
  dTime= dTime_in
  dTimeTooSmall=.FALSE.
  !
  nMkA= COUNT(vLKinActiv) !number of "active" kinetic phases
  nEqA= COUNT(vLEquActiv) !number of "active" equil' phases
  !
  IF(DirectSub) THEN  ;  nF= nCi
  ELSE                ;  nF= nAq
  ENDIF
  !
  nDim= nF +nMkA +nEqA
  !
  ALLOCATE(vX(nDim))
  ALLOCATE(vX0(nDim))
  ALLOCATE(vName (nDim))
  ALLOCATE(vDelta(nDim))
  ALLOCATE(vMinim(nDim))
  ALLOCATE(vXisPlus(nDim))
  ALLOCATE(vTolF(nDim))
  !
  ALLOCATE(vMolF1(SIZE(vMolF)))  ;  vMolF1=vMolF
  ALLOCATE(vMolK1(SIZE(vMolK)))  ;  vMolK1=vMolK
  !
  ALLOCATE(vMolSec(nAs))
  !
  !CALL Dynam_OneStep_Alloc(nF +nMkA +nEqA, nAq -nF)
  !
  vXisPlus(:)= .FALSE.
  IF(.NOT. LogForAqu) vXisPlus(1:nF)= .TRUE.
  !
  ALLOCATE(vLnGamOld(SIZE(vLnGam)))
  !
  !---------------------------------------------loop on activity coeff's
  DoGamma: DO
  !
  CALL Dynam_OneStep_DLGammaUpd(vLnGam, vDLGam_As)
  !
  iDo1= 0
  !--------------------------------------------Solve for time laps dTime
  DO1: DO

    iDo1=iDo1+1
    !
    !--------------------------------------------------------assemble vX
    !
    !-------------------------------------------assemble vX: aqu'species
    IF(LogForAqu) THEN
      !-- switch to Log's for aqu'species --
      vX(1:nF)= LOG(vMolF1(1:nF))
    ELSE
      vX(1:nF)= vMolF1(1:nF)
    ENDIF
    !
    IF(F>0) vName(1:nF)= vSpc(vOrdAq(1:nF))%NamSp(1:7)
    !-------------------------------------------------------/aqu'species
    !
    !----------------------------------------assemble vX: minerals et al
    I=0
    DO iMk=1,SIZE(vLEquActiv)
      IF(vLEquActiv(iMk)) THEN
        I= I + 1
        vMinim(I)= vKinMinim(iMk) /1.5D0
        IF(F>0) vName(nF+I)= vKinFas(I)%NamKF(1:7)
      ENDIF
    ENDDO
    DO iMk=1,SIZE(vLKinActiv)
      IF(vLKinActiv(iMk)) THEN
        I= I + 1
        vMinim(I)= vKinMinim(iMk) /1.5D0
        IF(F>0) vName(nF+I)= vKinFas(I)%NamKF(1:7)
      ENDIF
    ENDDO
    !
    !-----------------------extrapolate mineral mole nr, or adjust dTime
    dTime_Min= dTime
    !
    I=0
    DO iMk=1,SIZE(vLEquActiv)
      
      IF(vLEquActiv(iMk)) THEN

        I= I + 1
        !
        IF(Extrapole) THEN
          MolK= vMolK1(iMk) +vMolM_Slope(iMk) *dTime
          !--
          IF(MolK<Zero) THEN
            dTime_Min= MIN( &
            &  dTime_Min, &
            &  ABS((vMolK1(iMk) -vKinMinim(iMk))/ vMolM_Slope(iMk)))
            MolK= vKinMinim(iMk)
          ENDIF
        ELSE
          MolK= vMolK1(iMk)
        ENDIF
        !--when LogForMin, switch to Log's for minerals
        IF(LogForMin) THEN  ;  vX(nF+I)= LOG(MolK)
        ELSE                ;  vX(nF+I)= MolK
        ENDIF
        !
      ENDIF
    
    ENDDO
    !
    DO iMk=1,SIZE(vLKinActiv)
      IF(vLKinActiv(iMk)) THEN

        I= I + 1
        !
        IF(Implicit_Surface) THEN

          MolK= vMolK1(iMk)

        ELSE
          !-- use as initial guess for vMolK1(iMk)
          IF(Extrapole) THEN
            !-- use as initial guess for vMolK1(iMk)
            !-- the value extrapolated using the mineral rate
            !-- observed at previous step
            !-- vMolM_Slope(:)- (vMolK1(:) -vMolK_0(:)) /dTime
            MolK= vMolK1(iMk) +vMolM_Slope(iMk) *dTime
            !
            IF(MolK<Zero) THEN
              !dTime_Min= MIN(dTime_Min,ABS(vMolK1(iMk) / vMolM_Slope(iMk)))
              dTime_Min= MIN( &
              &  dTime_Min, &
              &  ABS((vMolK1(iMk) -vKinMinim(iMk))/ vMolM_Slope(iMk)))
              MolK= vKinMinim(iMk)
            ENDIF
          ELSE
            MolK= vMolK1(iMk)
          ENDIF
          !
        ENDIF
        !
        !--when LogForMin, switch to Log's for minerals
        IF(LogForMin) THEN ; vX(nF+I)= LOG(MolK)
        ELSE               ; vX(nF+I)= MolK
        ENDIF
        !
      ENDIF
    ENDDO
    !
    dTime= MAX(dTime_Min,dTmin*Two)
    !----------------------/extrapolate mineral mole nr, or adjust dTime
    !
    IF(.NOT. Implicit_Surface) & !
    & CALL Dynam_OneStep_KinFasSurf_Update( & !
    & vLKinActiv,   & ! IN
    & vStatusK,     & ! IN
    & vMolK0,vMolK1,& ! IN
    & PhiF0,PhiF,   & ! IN
    & vSurfK0,      & ! IN
    & vSurfK)         ! OUT
    !
    !----------------------------------------------------/minerals et al
    !
    !IF(iDebug>2) THEN
    !  WRITE(72,'(A,2G15.6)') "vX,Min/Max",MAXVAL(vX(:)),MINVAL(vX(:))
    !ENDIF
    !-------------------------------------------------------/assemble vX
    !
    ! print *,"debug Dynam_OneStep_Solve"
    ! do i=1,size(vX)
    !   print *,vName(i),vX(i)
    ! end do
    ! pause
    !iMethod= 2
    !
    vTolF(:)= NewtTolF
    !--- convergence criterion for aqu'species equilibrium conditions --
    IF(.NOT. DirectSub) vTolF(1:nCi+1:nCi+nAs)= NewtTolF_Equil
    !
    vX0(:)= vX(:)
    !
    iDo2=1
    DO2: DO

      iDo2=iDo2+1
      !
      vX(:)= vX0(:)
      !------------------------------------------------------CALL SOLVER
      CALL Dynam_Solve( & !
      & cMethod,        & !
      & vX,             & !INOUT
      & vXisPlus,       & !index of items that should stay >0
      & NewtTolF,NewtTolX,NewtMaxIts, & !IN,NewtTolMin
      & NewtErrF,NewtErrX,NewtErrG,Newt_iDo,Newt_Check,Newt_iErr) !OUT
      !-----------------------------------------------------/CALL SOLVER
      !
      vDelta=   ABS( vX0(1:nF+nMkA) - vX(1:nF+nMkA) )
      VarMax=   MAXVAL(vDelta)
      iMaxDelta=iMaxLoc_R(vDelta) !index of species with largest variation
      !
      IF(nMkA+nEqA>0) THEN
        IF(.NOT.LogForMin .AND. (Newt_iErr>=0)) THEN
          !
          !IF Newton is OK, and log is not used for minerals
          !check that no mineral has "too negative" amount
          !
          !IF, for any mineral, vol.fraction is significantly <0
          ! (it may be because precipitation was "too quick"),
          !THEN set error code iErr=-6
          !-> induces new calculation with lower timestep
          !
          !! if(iDebug>2) write(51,'(G15.6)') vX(nF+1)
          !
          IF( ANY(vX(nF+1:nF+nEqA+nMkA) < -vMinim(1:nEqA+nMkA)) ) THEN
            Newt_iErr=-6
            IF(iDebug>2) print *, "negative phase" !! ; pause
            !find the most negative phase
            !! J= iMinLoc_R(vMinim(1:nEqA+nMkA))
            !! -> adjust dTime
          ENDIF
          !
          WHERE( ABS(vX(nF+1:nF+nEqA+nMkA)) < vMinim(1:nEqA+nMkA) ) &
          & vX(nF+1:nF+nEqA+nMkA)= vMinim(1:nEqA+nMkA)
          !
        END IF
      END IF
      !
      IF(iDebug>1) PRINT '(I7,3(G12.3,1X),6I4,1X,3G12.3)', &
      & iStep, &
      & Time/TimeScale,dTime/TimeScale,PhiF, &
      & nEqA,nMkA, &
      & iDo2,Newt_iDo,Newt_iErr,iMaxDelta, &
      & VarMax,NewtErrF,NewtErrX
      !
      IF(F>0) &
      & WRITE(F,'(A,A1,I7,A1,3(G12.3,A1), I7,A1,G12.3,A1, 4(I7,A1),3(G15.6,A1))') &
      & TRIM(vName(iMaxDelta)),T_, &
      & iStep,T_,Time/TimeScale,T_,dTime/TimeScale,T_,PhiF,T_,iMaxDelta,T_,VarMax,T_, &
      & iDo1,T_,iDo2,T_,Newt_iDo,T_,&
      & Newt_iErr,T_,NewtErrF,T_,NewtErrX,T_,NewtErrG,T_
      !
      dTime_out= dTime
      !
      IF(Newt_iErr==0) EXIT DO2 !-------------------------------EXIT_DO2
      !
      dTime= dTime *Time_Decrease
      !
      dTime_out= dTime
      IF(dTime<dTmin) THEN
        DTimeTooSmall=.TRUE.
        !RETURN !----------------------- RETURN with TimeTooSmall .true.
        EXIT DO1 !---------------------EXIT_Do1 with TimeTooSmall .true.
      ENDIF
      !
    ENDDO DO2
    !
    IF (PhiF>Zero) THEN
      ! actually, should have always PhiF>Zero,
      ! because PhiF is implicited in FuncPsi ???
      !
      IF(DirectSub) CALL Dynam_OneStep_SecondSpecies(nF,vX,vMolSec)
      !
      IF(LogForAqu) THEN
        !--- Back to Mole Quantities
        IF(DirectSub) THEN
          vMolF1(1:nCi)=         EXP(vX(1:nF))
          vMolF1(nCi+1:nCi+nAs)= vMolSec(1:nAs)
        ELSE
          vMolF1(1:nF)= EXP(vX(1:nF))
        ENDIF
      ELSE
        IF(DirectSub) THEN
          vMolF1(1:nCi)=         vX(1:nF)
          vMolF1(nCi+1:nCi+nAs)= vMolSec(1:nAs)
        ELSE
          vMolF1(1:nF)= vX(1:nF)
        ENDIF
      ENDIF
      !
      !! IF(LogForAqu) THEN
      !!   vMolF1(1:nAq)= EXP(vX(1:nAq))
      !! ELSE
      !!   vMolF1(1:nAq-nAs)= vX(1:nAq-nAs)
      !!   CALL Compute_SecondSpecies(vMolF1)
      !! ENDIF
      !
      I=0
      DO iMk=1,SIZE(vLEquActiv) !-> mole nrs equil'phases
        IF(vLEquActiv(iMk)) THEN
          I= I + 1
          IF(LogForMin) THEN ;  vMolK1(iMk)= EXP(vX(nF+I))
          ELSE               ;  vMolK1(iMk)= vX(nF+I)
          ENDIF
        ENDIF
      ENDDO
      DO iMk=1,SIZE(vLKinActiv)
        IF(vLKinActiv(iMk)) THEN !-> mole nrs kin'phases
          I= I + 1
          IF(LogForMin) THEN ; vMolK1(iMk)= EXP(vX(nF+I))
          ELSE               ; vMolK1(iMk)= vX(nF+I)
          ENDIF
        ENDIF
      ENDDO
      !
      EXIT DO1 !------------------------------------------------EXIT_DO1
    ENDIF
    !
    dTime= dTime *Time_Decrease
    dTime_out= dTime
    !
    IF(dTime<dTmin) THEN
      DTimeTooSmall=.TRUE.
      EXIT DO1 !---------------------------------EXIT_Do1 + TimeTooSmall
    ENDIF
    !
  ENDDO DO1
  !-------------------------------------------/Solve for time laps dTime
  !
  vLnGamOld= vLnGam
  CALL Dynam_OneStep_GammaUpd( &
  & TdgK,Pbar, &
  & vMolF1,vLnAct,vLnGam, &
  & Z_Plus,Z_Minus,pH_)
  !
  !! IF(Gamma_NoLoop) EXIT
  IF( MAXVAL(ABS(vLnGamOld - vLnGam)) < DynTolGam ) EXIT DoGamma
  !
  ENDDO DoGamma
  !--------------------------------------------/loop on activity coeff's
  !
  vMolF= vMolF1  ;  DEALLOCATE(vMolF1)
  vMolK= vMolK1  ;  DEALLOCATE(vMolK1)
  !
  DEALLOCATE(vX)
  DEALLOCATE(vX0)
  DEALLOCATE(vDelta)
  DEALLOCATE(vMinim)
  DEALLOCATE(vName)
  DEALLOCATE(vXisPlus)
  DEALLOCATE(vTolF)
  DEALLOCATE(vMolSec)
  !
  DEALLOCATE(vLnGamOld)
  !CALL Dynam_OneStep_Clean
  !
ENDSUBROUTINE Dynam_OneStep_Solve

!-----------------------------------------------------------------------
!-- Dynam_Solve - interface with solvers -------------------------------
!-----------------------------------------------------------------------
SUBROUTINE Dynam_Solve( & !
& cMethod, & !
& vX,      & !inout= the initial guess, and the root returned
& vXIsPlus,& !in=    enforce vX(vIsPos=TRUE)>0
& TolF,    & !in=    convergence criterion on function values (scalar)
& TolX,    & !in=    convergence criterion on dx
& MaxIts,  & !in=    maximum number of iterations
& Error_F, & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Gradient,& !out=
& Nits,    & !out=   number of iterations
& Check,   & !out=   if Check, should check convergence
& iErr)      !out=   error code
  !
  USE M_Numeric_Newton
  USE M_Numeric_Broyden
  ! USE M_Numeric_Tensolve
  !
  USE M_Dynam_Vars,ONLY: bFinDIF
  !
  USE M_Dynam_Residual
  USE M_Dynam_Jacobian
  !
  CHARACTER(LEN=*),INTENT(IN):: cMethod
  REAL(dp), INTENT(INOUT):: vX(:)
  LOGICAL,  INTENT(IN)   :: vXisPlus(:)
  REAL(dp), INTENT(IN)   :: TolF,TolX
  INTEGER,  INTENT(IN)   :: MaxIts
  REAL(dp), INTENT(OUT)  :: Error_F,Delta_X,Gradient
  LOGICAL,  INTENT(OUT)  :: Check
  INTEGER,  INTENT(OUT)  :: nIts,iErr
  !
  Error_F= Zero
  Delta_X= Zero
  Gradient= Zero
  Check= .FALSE.
  !
  ALLOCATE(vNewtTolF(SIZE(vX)))
  ALLOCATE(vTolCoef(SIZE(vX)))
  vTolCoef(:)=  One
  vNewtTolF(:)= vTolF(:)
  !
  SELECT CASE(TRIM(cMethod))

  CASE("NEWTONPRESS")
    CALL Newton_Press( & !
    & vX,       & !initial guess x for a root in n dimensions
    & vXisPlus, & !in=

    & Dynam_Residual, & !
    & Dynam_Jacobian, & !
    & Dynam_Converge, & !

    & TolF,     & !in=    convergence criterion on function values
    & TolX,     & !in=    convergence criterion on dx
    & bFinDif,  & !in=    use numeric Jacobian
    & MaxIts,   & !in=    maximum number of iterations

    & Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,     & !out=   number of iterations
    & iErr)       !out=   error code

  CASE("NEWTONKELLEY")
    CALL Newton_Kelley( & !
    & vX,       & !initial guess x for a root in n dimensions
    & vXisPlus, & !in=

    & Dynam_Residual, & !
    & Dynam_Jacobian, & !
    & Dynam_Converge, & !

    & TolF,     & !in=    convergence criterion on function values
    & TolX,     & !in=    convergence criterion on dx
    & bFinDif,  & !in=    use numeric Jacobian
    & MaxIts,   & !in=    maximum number of iterations

    & Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,     & !out=   number of iterations
    & iErr      & !out=   error code
    & )

  CASE("NEWTONWALKER")
    CALL Newton_Walker( & !
    & vX,       & !initial guess x for a root in n dimensions
    & vXisPlus, & !in=

    & Dynam_Residual, & !
    & Dynam_Jacobian, & !
    & Dynam_Converge, & !

    & TolF,     & !in=    convergence criterion on function values
    & TolX,     & !in=    convergence criterion on dx
    & bFinDIF,  & !in=    use numeric Jacobian
    & MaxIts,   & !in=    maximum number of iterations

    & Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,     & !out=   number of iterations
    & iErr)       !out=   error code

  CASE("NEWTLNSRCH")
    CALL NewtLnsrch( & !
    & vX,      & !inout= the initial guess, and the root RETURNed
    !!& vLPos,   & !in=    enforce vX(vLPos=TRUE)>0

    & Dynam_Residual, & !
    & Dynam_Jacobian, & !
    & Dynam_Converge, & !

    & TolF,    & !in=    convergence criterion on FUNCTION values
    & TolX,    & !in=    convergence criterion on dx
    !!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
    & bFinDIF, & !in=    USE numeric Jacobian
    & MaxIts,  & !in=    maximum number of iterations

    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Gradient,& !out=
    & Nits,    & !out=   number of iterations
    & Check,   & !out=   IF Check, should check convergence
    & iErr)      !out=   error code

    ! CALL Newton( &
    ! & vX,      & !inout= the initial guess, and the root RETURNed
    ! !& vLPos,   & !in=    enforce vX(vLPos=TRUE)>0
    ! & Dynam_Residual, Dynam_Jacobian, &
    ! & TolF,    & !in=    convergence criterion on FUNCTION values
    ! & TolX,    & !in=    convergence criterion on dx
    ! !!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
    ! !& bFinDIF, & !in=
    ! & MaxIts,  & !in=    maximum number of iterations
    ! & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    ! & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    ! !& Gradient, & !out=
    ! & Nits,    & !out=   number of iterations
    ! !& Check,   & !out=   IF Check, should check convergence
    ! & iErr)      !out=   error code

  CASE("NEWTONCHESS")
    CALL NewtonChess( & !
    & vX,       & !inout= the initial guess, and the root returned

    & Dynam_Residual, & !
    & Dynam_Jacobian, & !

    & TolF,     & !in=    convergence criterion on FUNCTION values
    & TolX,     & !in=    convergence criterion on dx
    & MaxIts,   & !in=    maximum number of iterations

    & Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,     & !out=   number of iterations
    & iErr)       !out=   error code

  CASE("BROYDEN")
    CALL Broyden( &
    & vX,      & !inout= the initial guess, and the root RETURNed

    & Dynam_Residual, &
    & Dynam_Jacobian, &

    & TolF,    & !in=    convergence criterion on FUNCTION values
    & TolX,    & !in=    convergence criterion on dx
    !!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
    & bFinDIF, & !in=    USE numeric Jacobian
    & MaxIts,  & !in=    maximum number of iterations

    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Gradient, & !out=
    & Nits,    & !out=   number of iterations
    & Check,   & !out=   IF Check, should check convergence
    & iErr)      !out=   error code

  ! CASE("TENSOLVE")
  !   bFinDIF= .FALSE. !.TRUE.
  !   CALL TenSolve( &
  !   & "NEWTON","TRUSTREGION", &
  !   & vX,       &
  !
  !   & Dynam_Residual, &
  !   & Dynam_Jacobian, &
  !
  !   & TolF,     & !in=    convergence criterion on FUNCTION values
  !   & bFinDIF,  & !in=    USE numeric Jacobian
  !   & MaxIts,   & !in=    maximum number of iterations
  !
  !   & Error_F,  & !out
  !   & Nits,     & !out=   number of iterations
  !   & iErr)       !out=   error code

  END SELECT
  !
  DEALLOCATE(vNewtTolF)
  DEALLOCATE(vTolCoef)
  !
  RETURN
ENDSUBROUTINE Dynam_Solve

ENDMODULE M_Dynam_OneStep_Solve


