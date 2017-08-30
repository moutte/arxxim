module M_TimeMng_Base

  !=====================================================
  ! TimeMng : Basic Functions
  !=====================================================

  use M_Kinds
  use M_Trace
  use M_T_TimeLine
  use M_T_TimeStep
  use M_T_TimeStepMng
  use M_TimeMng_Vars
  
  implicit none
  private
  
  !// public functionS
  public :: TimeMng_Zero
  public :: TimeMng_Set_TimeLine
  public :: TimeMng_Set_TimeStepMng
  public :: TimeMng_Set_DTRatio
  public :: TimeMng_Reset_TimeStep
  public :: TimeMng_Check
  public :: TimeMng_Info
  
contains

  
  subroutine TimeMng_Zero
    implicit none
    
    !// TimeMng
    call TimeLine_Zero(TimeLine)
    call TimeStepMng_Zero(TimeStepMng)
    call TimeStep_Zero(TimeStep)
    
    !// TimeStep Status
    call TimeMng_Reset_TimeStep
    
  end subroutine TimeMng_Zero

  !---

  subroutine TimeMng_Set_TimeLine(T)
    
    implicit none
    type(T_TimeLine), intent(in) :: T
    !
    TimeLine% TInit=  T% TInit
    TimeLine% TFinal= T% TFinal
    
    !// Check TimeLine Validity
    call TimeLine_Check(T)
    
  end subroutine TimeMng_Set_TimeLine
  
  !---
    
  subroutine TimeMng_Set_TimeStepMng(T)
    implicit none
    !
    type(T_TimeStepMng), intent(in) :: T
    !
    TimeStepMng% DTInit= T% DTInit
    TimeStepMng% DTMin=  T% DTMin
    TimeStepMng% DTMax=  T% DTMax
    TimeStepMng% DTRatio_Increase= T% DTRatio_Increase
    TimeStepMng% DTRatio_Reduce=   T% DTRatio_Reduce
        
    !// Check TimeStepMng Validity
    call TimeStepMng_Check(T)

  end subroutine TimeMng_Set_TimeStepMng

  !--

  subroutine TimeMng_Set_DTRatio( DTRatio_Increase, DTRatio_Reduce ) 
    !===================================================
    ! Purpose : Set TimeStepMng Ratios Only
    !===================================================
    implicit none
    !
    real(dp), intent(in) :: DTRatio_Increase, DTRatio_Reduce 
    !-----
    call TimeStepMng_Set_DTRatio( TimeStepMng, DTRatio_Increase, DTRatio_Reduce ) 
    
  end subroutine TimeMng_Set_DTRatio

   !---

  subroutine TimeMng_Reset_TimeStep
    !===================================================
    ! Purpose : Reset the TimeStep to the Initial Values
    !===================================================
    implicit none
    !---
    TimeStep% Time  = TimeLine% TInit
    TimeStep% DTime = TimeStepMng% DTInit

    !// TimeStep Status
    ERRORTimeStepTooSmall =.false.
    OKFinalTimeStep =.false.
    
  end subroutine TimeMng_Reset_TimeStep

  !---
  
  subroutine TimeMng_Check
    !===================================================
    ! Purpose : Check TimeMng Minimal Validity
    !===================================================
    implicit none
    !---
    
    call TimeLine_Check(TimeLine)
    call TimeStep_Check(TimeStep)
    call TimeStepMng_Check(TimeStepMng)
    
    !-- TInit <= Time <= TFinal
    if ( TimeStep% Time <  TimeLine% TInit ) then
       call Fatal_("Time < TInit")
    end if
    
    if ( TimeStep% Time + TimeStep% DTime >  TimeLine% TFinal ) then
       call Fatal_("Time+DTime > TFinal")
    end if

    !-- DTMin <= DTime <= DTMax
    if ( TimeStep% DTime> TimeStepMng% DTMin) then
       call Fatal_("DTime < DTMin")
    end if

    if ( TimeStep% DTime< TimeStepMng% DTMax) then
       call Fatal_("DTime > DTMax")
    end if

  end subroutine TimeMng_Check

  !---

  subroutine TimeMng_Info
  
    implicit none
    !
    call Info_("TimeMng Info") 
    call TimeLine_Info(TimeLine)
    call TimeStep_Info(TimeStep)
    call TimeStepMng_Info(TimeStepMng)
    
  end subroutine TimeMng_Info
end module M_TimeMng_Base
