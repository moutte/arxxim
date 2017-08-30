module M_TimeMng_Control

  !=====================================================
  ! TimeMng : High Level DTime Control Functions
  !=====================================================

  use M_Kinds
  use M_TimeMng_Vars
  implicit none
  private

  !// public functionS
  public:: TimeMng_Control_TimeStep

  private:: TimeMng_DTime_Increase
  private:: TimeMng_DTime_Reduce

contains

  subroutine TimeMng_Control_TimeStep(OKTimeStep)

    !========================================================
    ! Purpose : TimeStep Control by DTRatio
    !========================================================
    implicit none
    logical, intent(in) :: OKTimeStep
    !
    if (OKTimeStep) then
       call TimeMng_DTime_Increase
    else
       call TimeMng_DTime_Reduce
    end if

  end subroutine TimeMng_Control_TimeStep

  !---

  subroutine TimeMng_DTime_Increase

    !========================================================
    ! Purpose : TimeStep Control by DTRatio OKTimeStep
    !--------------------------------------------------------
    ! OKTimeStep => Increase TimeStep if OKTimeStep
    ! Mult by DTRatio_Increase and Check DTMax, Tfinal
    !========================================================
    implicit none
    !
    real(dp) :: TInit, TFinal
    real(dp) :: DTMax, IncreaseFactor
    real(dp) :: Time, DTime
    real(dp) :: dTime_new
    real(dp) :: epsilon = 1.d-14 
    !
    TInit=          TimeLine% TInit
    TFinal=         TimeLine% TFinal
    DTMax=          TimeStepMng% DTMax
    IncreaseFactor= TimeStepMng% DTRatio_Increase
    Time=           TimeStep% Time
    DTime=          TimeStep% DTime

    ERRORTimeStepTooSmall = .false.

    !// Increase the TimeStep if Time > Tinit    
    if (Time==TInit) then
       dTime_new= dTime 
    else
       dTime_new= dTime * IncreaseFactor
    end if

    dTime_new = Min(dTime_new, DTMax)

    !// Adjust the Time Step to Final Time
    OKFinalTimeStep= .false.
    if (Time + dTime_new >= Tfinal - epsilon) then
       dTime_new= ( Tfinal - Time )
       OKFinalTimeStep= .true.

    elseif (Time + 2*dTime_new > Tfinal) then
       dTime_new= ( Tfinal - Time )*0.5
       OKFinalTimeStep= .false.
    end if

    TimeStep%DTime = dTime_new

  end subroutine TimeMng_DTime_Increase

  !---

  subroutine TimeMng_DTime_Reduce

    !========================================================
    ! Purpose : TimeStep Control by DTRatio if ERRORTimeStep
    !           Update DTime + ERRORTimeStepTooSmall
    !--------------------------------------------------------
    !  NOT(OKTimeStep) => Reduce TimeStep
    !  Mult by DTRatio_Recuce and check DTMin
    !========================================================
    implicit none
    !
    real(kind=8) :: DTime
    real(kind=8) :: DTMin, ReduceFactor
    real(kind=8) :: dTime_new
    !
    DTMin=        TimeStepMng% DTMin
    ReduceFactor= TimeStepMng% DTRatio_Reduce
    DTime=        TimeStep% DTime

    OKFinalTimeStep= .false.

    !// Reduce the Time Step
    dTime_new = dTime* ReduceFactor

    ERRORTimeStepTooSmall = .false.
    if (dTime_new < DTMin ) ERRORTimeStepTooSmall=.true.

    TimeStep%DTime = dTime_new

  end subroutine TimeMng_DTime_Reduce


end module M_TimeMng_Control
