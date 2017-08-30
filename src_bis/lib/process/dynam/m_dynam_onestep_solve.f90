module M_Dynam_OneStep_Solve
  use M_Kinds
  use M_Trace,only: iDebug,fTrc,T_
  implicit none

  private

  public:: Dynam_OneStep_Solve

  real(dp),dimension(:),allocatable:: vTolF

contains

!-----------------------------------------------------------------------
!-- solve the dynamic system for one time step--------------------------
!-----------------------------------------------------------------------
subroutine Dynam_OneStep_Solve( & !
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
  use M_Numeric_Tools,only: iMaxLoc_R,iMinLoc_R
  use M_IoTools
  !
  use M_System_Vars,only: TdgK,Pbar
  use M_Global_Vars,only: vSpc,nAq,vKinFas
  use M_Basis_Vars, only: vOrdAq,nCi,nAs
  !
  use M_Dynam_OneStep_Tools,only: Dynam_OneStep_Secondspecies
  use M_Dynam_OneStep_Tools,only: Dynam_OneStep_KinFasSurf_Update
  use M_Dynam_OneStep_Tools,only: Dynam_OneStep_DLGammaUpd
  use M_Dynam_OneStep_Tools,only: Dynam_OneStep_GammaUpd
  !
  use M_Dynam_Vars, only: &
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
  use M_Dynam_Vars, only: PhiF0,vStatusK,vSurfK,vSurfK0,vMolK0
  !
  integer, intent(in) :: F,iStep
  logical, intent(in) :: vLKinActiv(:),vLEquActiv(:)
  real(dp),intent(in) :: vMolM_Slope(:)
  real(dp),intent(in) :: NewtTolF,NewtTolF_Equil
  real(dp),intent(in) :: NewtTolX,NewtTolMin
  real(dp),intent(in) :: Time_Decrease
  integer, intent(in) :: NewtMaxIts
  real(dp),intent(in) :: Time,dTmin
  real(dp),intent(in) :: dTime_in
  !
  real(dp),intent(out):: dTime_out
  integer, intent(out):: Newt_iDo,Newt_iErr
  logical, intent(out):: DTimeTooSmall
  real(dp),intent(out):: VarMax
  !---------------------------------------------------------------------
  real(dp),dimension(:),allocatable:: vX,vX0,vDelta,vMinim
  real(dp),dimension(:),allocatable:: vMolF1,vMolK1
  real(dp),dimension(:),allocatable:: vLnGamOld
  logical,allocatable:: vXisPlus(:)
  character(len=7),allocatable:: vName(:)
  !
  integer :: iDo1,iDo2,iMaxDelta,iMk,I
  integer :: nF,nMkA,nEqA,nDim
  real(dp):: NewtErrF,NewtErrX,NewtErrG
  real(dp):: dTime_Min
  real(dp):: MolK
  logical :: Newt_Check
  real(dp):: Z_Plus,Z_Minus
  !
  real(dp),parameter:: DynTolGam= 1.D-3
  !---------------------------------------------------------------------
  dTime= dTime_in
  dTimeTooSmall=.false.
  !
  nMkA= count(vLKinActiv) !number of "active" kinetic phases
  nEqA= count(vLEquActiv) !number of "active" equil' phases
  !
  if(DirectSub) then  ;  nF= nCi
  else                ;  nF= nAq
  end if
  !
  nDim= nF +nMkA +nEqA
  !
  allocate(vX(nDim))
  allocate(vX0(nDim))
  allocate(vName (nDim))
  allocate(vDelta(nDim))
  allocate(vMinim(nDim))
  allocate(vXisPlus(nDim))
  allocate(vTolF(nDim))
  !
  allocate(vMolF1(size(vMolF)))  ;  vMolF1=vMolF
  allocate(vMolK1(size(vMolK)))  ;  vMolK1=vMolK
  !
  allocate(vMolSec(nAs))
  !
  !call Dynam_OneStep_Alloc(nF +nMkA +nEqA, nAq -nF)
  !
  vXisPlus(:)= .false.
  if(.not. LogForAqu) vXisPlus(1:nF)= .true.
  !
  allocate(vLnGamOld(size(vLnGam)))
  !
  !---------------------------------------------loop on activity coeff's
  DoGamma: do
  !
  call Dynam_OneStep_DLGammaUpd(vLnGam, vDLGam_As)
  !
  iDo1= 0
  !--------------------------------------------Solve for time laps dTime
  do1: do

    iDo1=iDo1+1
    !
    !--------------------------------------------------------assemble vX
    !
    !-------------------------------------------assemble vX: aqu'species
    if(LogForAqu) then
      !-- switch to Log's for aqu'species --
      vX(1:nF)= log(vMolF1(1:nF))
    else
      vX(1:nF)= vMolF1(1:nF)
    end if
    !
    if(F>0) vName(1:nF)= vSpc(vOrdAq(1:nF))%NamSp(1:7)
    !-------------------------------------------------------/aqu'species
    !
    !----------------------------------------assemble vX: minerals et al
    I=0
    do iMk=1,size(vLEquActiv)
      if(vLEquActiv(iMk)) then
        I= I + 1
        vMinim(I)= vKinMinim(iMk) /1.5D0
        if(F>0) vName(nF+I)= vKinFas(I)%NamKF(1:7)
      end if
    enddo
    do iMk=1,size(vLKinActiv)
      if(vLKinActiv(iMk)) then
        I= I + 1
        vMinim(I)= vKinMinim(iMk) /1.5D0
        if(F>0) vName(nF+I)= vKinFas(I)%NamKF(1:7)
      end if
    enddo
    !
    !-----------------------extrapolate mineral mole nr, or adjust dTime
    dTime_Min= dTime
    !
    I=0
    do iMk=1,size(vLEquActiv)
      
      if(vLEquActiv(iMk)) then

        I= I + 1
        !
        if(Extrapole) then
          MolK= vMolK1(iMk) +vMolM_Slope(iMk) *dTime
          !--
          if(MolK<Zero) then
            dTime_Min= MIN( &
            &  dTime_Min, &
            &  ABS((vMolK1(iMk) -vKinMinim(iMk))/ vMolM_Slope(iMk)))
            MolK= vKinMinim(iMk)
          end if
        else
          MolK= vMolK1(iMk)
        end if
        !--when LogForMin, switch to Log's for minerals
        if(LogForMin) then  ;  vX(nF+I)= log(MolK)
        else                ;  vX(nF+I)= MolK
        end if
        !
      end if
    
    enddo
    !
    do iMk=1,size(vLKinActiv)
      if(vLKinActiv(iMk)) then

        I= I + 1
        !
        if(Implicit_Surface) then

          MolK= vMolK1(iMk)

        else
          !-- use as initial guess for vMolK1(iMk)
          if(Extrapole) then
            !-- use as initial guess for vMolK1(iMk)
            !-- the value extrapolated using the mineral rate
            !-- observed at previous step
            !-- vMolM_Slope(:)- (vMolK1(:) -vMolK_0(:)) /dTime
            MolK= vMolK1(iMk) +vMolM_Slope(iMk) *dTime
            !
            if(MolK<Zero) then
              !dTime_Min= MIN(dTime_Min,ABS(vMolK1(iMk) / vMolM_Slope(iMk)))
              dTime_Min= MIN( &
              &  dTime_Min, &
              &  ABS((vMolK1(iMk) -vKinMinim(iMk))/ vMolM_Slope(iMk)))
              MolK= vKinMinim(iMk)
            end if
          else
            MolK= vMolK1(iMk)
          end if
          !
        end if
        !
        !--when LogForMin, switch to Log's for minerals
        if(LogForMin) then ; vX(nF+I)= log(MolK)
        else               ; vX(nF+I)= MolK
        end if
        !
      end if
    enddo
    !
    dTime= MAX(dTime_Min,dTmin*Two)
    !----------------------/extrapolate mineral mole nr, or adjust dTime
    !
    if(.not. Implicit_Surface) & !
    & call Dynam_OneStep_KinFasSurf_Update( & !
    & vLKinActiv,   & ! IN
    & vStatusK,     & ! IN
    & vMolK0,vMolK1,& ! IN
    & PhiF0,PhiF,   & ! IN
    & vSurfK0,      & ! IN
    & vSurfK)         ! OUT
    !
    !----------------------------------------------------/minerals et al
    !
    !if(iDebug>2) then
    !  write(72,'(A,2G15.6)') "vX,Min/Max",MAXVAL(vX(:)),MINVAL(vX(:))
    !end if
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
    if(.not. DirectSub) vTolF(1:nCi+1:nCi+nAs)= NewtTolF_Equil
    !
    vX0(:)= vX(:)
    !
    iDo2=1
    do2: do

      iDo2=iDo2+1
      !
      vX(:)= vX0(:)
      !------------------------------------------------------call SOLVER
      call Dynam_Solve( & !
      & cMethod,        & !
      & vX,             & !INOUT
      & vXisPlus,       & !index of items that should stay >0
      & NewtTolF,NewtTolX,NewtMaxIts, & !IN,NewtTolMin
      & NewtErrF,NewtErrX,NewtErrG,Newt_iDo,Newt_Check,Newt_iErr) !OUT
      !-----------------------------------------------------/call SOLVER
      !
      vDelta=   ABS( vX0(1:nF+nMkA) - vX(1:nF+nMkA) )
      VarMax=   MAXVAL(vDelta)
      iMaxDelta=iMaxLoc_R(vDelta) !index of species with largest variation
      !
      if(nMkA+nEqA>0) then
        if(.not.LogForMin .and. (Newt_iErr>=0)) then
          !
          !if Newton is OK, and log is not used for minerals
          !check that no mineral has "too negative" amount
          !
          !if, for any mineral, vol.fraction is significantly <0
          ! (it may be because precipitation was "too quick"),
          !then set error code iErr=-6
          !-> induces new calculation with lower timestep
          !
          !! if(iDebug>2) write(51,'(G15.6)') vX(nF+1)
          !
          if( ANY(vX(nF+1:nF+nEqA+nMkA) < -vMinim(1:nEqA+nMkA)) ) then
            Newt_iErr=-6
            if(iDebug>2) print *, "negative phase" !! ; pause
            !find the most negative phase
            !! J= iMinLoc_R(vMinim(1:nEqA+nMkA))
            !! -> adjust dTime
          end if
          !
          where( ABS(vX(nF+1:nF+nEqA+nMkA)) < vMinim(1:nEqA+nMkA) ) &
          & vX(nF+1:nF+nEqA+nMkA)= vMinim(1:nEqA+nMkA)
          !
        end if
      end if
      !
      if(iDebug>1) print '(I7,3(G12.3,1X),6I4,1X,3G12.3)', &
      & iStep, &
      & Time/TimeScale,dTime/TimeScale,PhiF, &
      & nEqA,nMkA, &
      & iDo2,Newt_iDo,Newt_iErr,iMaxDelta, &
      & VarMax,NewtErrF,NewtErrX
      !
      if(F>0) &
      & write(F,'(A,A1,I7,A1,3(G12.3,A1), I7,A1,G12.3,A1, 4(I7,A1),3(G15.6,A1))') &
      & trim(vName(iMaxDelta)),T_, &
      & iStep,T_,Time/TimeScale,T_,dTime/TimeScale,T_,PhiF,T_,iMaxDelta,T_,VarMax,T_, &
      & iDo1,T_,iDo2,T_,Newt_iDo,T_,&
      & Newt_iErr,T_,NewtErrF,T_,NewtErrX,T_,NewtErrG,T_
      !
      dTime_out= dTime
      !
      if(Newt_iErr==0) exit do2 !-------------------------------exit_do2
      !
      dTime= dTime *Time_Decrease
      !
      dTime_out= dTime
      if(dTime<dTmin) then
        DTimeTooSmall=.true.
        !return !----------------------- return with TimeTooSmall .true.
        exit do1 !---------------------exit_Do1 with TimeTooSmall .true.
      end if
      !
    enddo do2
    !
    if (PhiF>Zero) then
      ! actually, should have always PhiF>Zero,
      ! because PhiF is implicited in FuncPsi ???
      !
      if(DirectSub) call Dynam_OneStep_SecondSpecies(nF,vX,vMolSec)
      !
      if(LogForAqu) then
        !--- Back to Mole Quantities
        if(DirectSub) then
          vMolF1(1:nCi)=         exp(vX(1:nF))
          vMolF1(nCi+1:nCi+nAs)= vMolSec(1:nAs)
        else
          vMolF1(1:nF)= exp(vX(1:nF))
        end if
      else
        if(DirectSub) then
          vMolF1(1:nCi)=         vX(1:nF)
          vMolF1(nCi+1:nCi+nAs)= vMolSec(1:nAs)
        else
          vMolF1(1:nF)= vX(1:nF)
        end if
      end if
      !
      !! if(LogForAqu) then
      !!   vMolF1(1:nAq)= exp(vX(1:nAq))
      !! else
      !!   vMolF1(1:nAq-nAs)= vX(1:nAq-nAs)
      !!   call Compute_SecondSpecies(vMolF1)
      !! end if
      !
      I=0
      do iMk=1,size(vLEquActiv) !-> mole nrs equil'phases
        if(vLEquActiv(iMk)) then
          I= I + 1
          if(LogForMin) then ;  vMolK1(iMk)= exp(vX(nF+I))
          else               ;  vMolK1(iMk)= vX(nF+I)
          end if
        end if
      enddo
      do iMk=1,size(vLKinActiv)
        if(vLKinActiv(iMk)) then !-> mole nrs kin'phases
          I= I + 1
          if(LogForMin) then ; vMolK1(iMk)= exp(vX(nF+I))
          else               ; vMolK1(iMk)= vX(nF+I)
          end if
        end if
      enddo
      !
      exit do1 !------------------------------------------------exit_do1
    end if
    !
    dTime= dTime *Time_Decrease
    dTime_out= dTime
    !
    if(dTime<dTmin) then
      DTimeTooSmall=.true.
      exit do1 !---------------------------------exit_Do1 + TimeTooSmall
    end if
    !
  enddo do1
  !-------------------------------------------/Solve for time laps dTime
  !
  vLnGamOld= vLnGam
  call Dynam_OneStep_GammaUpd( &
  & TdgK,Pbar, &
  & vMolF1,vLnAct,vLnGam, &
  & Z_Plus,Z_Minus,pH_)
  !
  !! if(Gamma_NoLoop) exit
  if( MAXVAL(ABS(vLnGamOld - vLnGam)) < DynTolGam ) exit DoGamma
  !
  enddo DoGamma
  !--------------------------------------------/loop on activity coeff's
  !
  vMolF= vMolF1  ;  deallocate(vMolF1)
  vMolK= vMolK1  ;  deallocate(vMolK1)
  !
  deallocate(vX)
  deallocate(vX0)
  deallocate(vDelta)
  deallocate(vMinim)
  deallocate(vName)
  deallocate(vXisPlus)
  deallocate(vTolF)
  deallocate(vMolSec)
  !
  deallocate(vLnGamOld)
  !call Dynam_OneStep_Clean
  !
end subroutine Dynam_OneStep_Solve

!-----------------------------------------------------------------------
!-- Dynam_Solve - interface with solvers -------------------------------
!-----------------------------------------------------------------------
subroutine Dynam_Solve( & !
& cMethod, & !
& vX,      & !inout= the initial guess, and the root returned
& vXIsPlus,& !in=    enforce vX(vIsPos.true.>0
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
  use M_Numeric_Newton
  use M_Numeric_Broyden
  ! use M_Numeric_Tensolve
  !
  use M_Dynam_Vars,only: bFinDif
  !
  use M_Dynam_Residual
  use M_Dynam_Jacobian
  !
  character(len=*),intent(in):: cMethod
  real(dp), intent(inout):: vX(:)
  logical,  intent(in)   :: vXisPlus(:)
  real(dp), intent(in)   :: TolF,TolX
  integer,  intent(in)   :: MaxIts
  real(dp), intent(out)  :: Error_F,Delta_X,Gradient
  logical,  intent(out)  :: Check
  integer,  intent(out)  :: nIts,iErr
  !
  Error_F= Zero
  Delta_X= Zero
  Gradient= Zero
  Check= .false.
  !
  allocate(vNewtTolF(size(vX)))
  allocate(vTolCoef(size(vX)))
  vTolCoef(:)=  One
  vNewtTolF(:)= vTolF(:)
  !
  select case(trim(cMethod))

  case("NEWTONPRESS")
    call Newton_Press( & !
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

  case("NEWTONKELLEY")
    call Newton_Kelley( & !
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

  case("NEWTONWALKER")
    call Newton_Walker( & !
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

  case("NEWTLNSRCH")
    call NewtLnsrch( & !
    & vX,      & !inout= the initial guess, and the root returned
    !!& vLPos,   & !in=    enforce vX(vLPos.true.>0

    & Dynam_Residual, & !
    & Dynam_Jacobian, & !
    & Dynam_Converge, & !

    & TolF,    & !in=    convergence criterion on function values
    & TolX,    & !in=    convergence criterion on dx
    !!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
    & bFinDif, & !in=    use numeric Jacobian
    & MaxIts,  & !in=    maximum number of iterations

    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Gradient,& !out=
    & Nits,    & !out=   number of iterations
    & Check,   & !out=   if Check, should check convergence
    & iErr)      !out=   error code

    ! call Newton( &
    ! & vX,      & !inout= the initial guess, and the root returned
    ! !& vLPos,   & !in=    enforce vX(vLPos.true.>0
    ! & Dynam_Residual, Dynam_Jacobian, &
    ! & TolF,    & !in=    convergence criterion on function values
    ! & TolX,    & !in=    convergence criterion on dx
    ! !!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
    ! !& bFinDif, & !in=
    ! & MaxIts,  & !in=    maximum number of iterations
    ! & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    ! & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    ! !& Gradient, & !out=
    ! & Nits,    & !out=   number of iterations
    ! !& Check,   & !out=   if Check, should check convergence
    ! & iErr)      !out=   error code

  case("NEWTONCHESS")
    call NewtonChess( & !
    & vX,       & !inout= the initial guess, and the root returned

    & Dynam_Residual, & !
    & Dynam_Jacobian, & !

    & TolF,     & !in=    convergence criterion on function values
    & TolX,     & !in=    convergence criterion on dx
    & MaxIts,   & !in=    maximum number of iterations

    & Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,     & !out=   number of iterations
    & iErr)       !out=   error code

  case("BROYDEN")
    call Broyden( &
    & vX,      & !inout= the initial guess, and the root returned

    & Dynam_Residual, &
    & Dynam_Jacobian, &

    & TolF,    & !in=    convergence criterion on function values
    & TolX,    & !in=    convergence criterion on dx
    !!& TolMin,  & !in=    whether spurious convergence to a minimum of fmin has occurred
    & bFinDif, & !in=    use numeric Jacobian
    & MaxIts,  & !in=    maximum number of iterations

    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Gradient, & !out=
    & Nits,    & !out=   number of iterations
    & Check,   & !out=   if Check, should check convergence
    & iErr)      !out=   error code

  ! case("TENSOLVE")
  !   bFinDif= .false. !.true.
  !   call TenSolve( &
  !   & "NEWTON","TRUSTREGION", &
  !   & vX,       &
  !
  !   & Dynam_Residual, &
  !   & Dynam_Jacobian, &
  !
  !   & TolF,     & !in=    convergence criterion on function values
  !   & bFinDif,  & !in=    use numeric Jacobian
  !   & MaxIts,   & !in=    maximum number of iterations
  !
  !   & Error_F,  & !out
  !   & Nits,     & !out=   number of iterations
  !   & iErr)       !out=   error code

  end select
  !
  deallocate(vNewtTolF)
  deallocate(vTolCoef)
  !
  return
end subroutine Dynam_Solve

end module M_Dynam_OneStep_Solve


