module M_T_TimeStep

  !======================================================
  ! T_TimeStep, Discretized Time Current State 
  ! Define a Time Interval [Time, Time+DTime]
  !======================================================
 
  use M_Kinds
  implicit none
  private
  
  public :: T_TimeStep
  type T_TimeStep
     real(dp) :: Time
     real(dp) :: DTime
  end type T_TimeStep

  !// Public Functions
  public :: TimeStep_Zero
  public :: TimeStep_Init
  public :: TimeStep_Check
  public :: TimeStep_Info

contains

  !---
  
 subroutine TimeStep_Zero( T )
    !===================================================
    ! Purpose : Init TimeLine by using given parameters
    !===================================================
    implicit none
    type(T_TimeStep) :: T
    !-----
    T% Time  = Zero
    T% DTime = Zero

  end subroutine TimeStep_Zero

  !---
  
  subroutine TimeStep_Init( T, Time, DTime)
    !===================================================
    ! Purpose : Init TimeLine by using given parameters
    !===================================================
    implicit none
    !
    type(T_TimeStep) :: T
    real(dp), intent(inout) :: Time, DTime
    
    !-----
    T% Time  = Time
    T% DTime = DTime

  end subroutine TimeStep_Init

  !---
  
  subroutine TimeStep_Check( T )
    !===================================================
    ! Purpose : Minimal Validity Check TimeStep
    !---
    ! DTime > 0
    !===================================================
    use M_Trace
    implicit none
    type(T_TimeStep) :: T
    !--
    if ( T% DTime < Zero ) then
       call Fatal_("Wrong TimeStep : DTime < Zero")
    end if
  end subroutine TimeStep_Check

  !---

  subroutine TimeStep_Info( T )
    !===================================================
    ! Information 
    !===================================================
    use M_Trace
    implicit none
    type(T_TimeStep) :: T
    !---
    call Info_("------ TimeStep Object")
    Write(*,*) "Time=", T% Time
    Write(*,*) "DTime=", T% DTime
    
  end subroutine TimeStep_Info
 
end module M_T_TimeStep
