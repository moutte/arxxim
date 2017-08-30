module M_TimeMng_Time

  !=====================================================
  ! TimeMng : High Level Functions Concerning Time
  !=====================================================

  use M_Kinds
  use M_Trace
  use M_TimeMng_Vars
  implicit none
  private

  !// Public functionS
  public :: TimeMng_Update_Time
  public :: TimeMng_Get_Time
  public :: TimeMng_Info_TimeStep

contains

  !---

  subroutine TimeMng_Update_Time
    !===============================================
    ! Purpose : Update Time : Time = Time + DTime
    !===============================================
    use M_Trace
    implicit none
    real(DP) :: Time_new
    !
    Time_new = TimeStep%Time + TimeStep%DTime
    if (Time_New > TimeLine% TFinal ) then
       call Fatal_("Error Update : New Time > TFinal")
    else
       TimeStep%Time = Time_new
    end if
  end subroutine TimeMng_Update_Time

  !--

  subroutine TimeMng_Get_Time(Time, DTime)
    !===================================================
    ! Purpose : Get TimeStep State Values
    !===================================================
    implicit None
    real(dp), intent(out) :: Time, DTime
    !
    Time=  TimeStep% Time
    DTime= TimeStep% DTime

  end subroutine TimeMng_Get_Time

  !---

  subroutine TimeMng_Info_TimeStep
    !===================================================
    ! Purpose : Small Message About TimeStep Values
    !===================================================
    implicit none
    character(len=255) :: Message
    !
    write(Message,'(A,F12.3,F12.3)') "[Time, DTime] =", TimeStep% Time, TimeStep% DTime
    call Info_(Message)

  end subroutine TimeMng_Info_TimeStep


end module M_TimeMng_Time
