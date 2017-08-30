module M_Box_TimeLoop_Vars

  use M_Kinds
  use M_Trace
  
  use M_TimeMng, only : &
       OKFinalTimeStep, &
       ERRORTimeStepTooSmall

  implicit none
  public 

  !// Loop Counter
  integer :: iTimeStep
  integer :: iTryStep

  !// TimeLoop Status
  logical :: StopTimeLoop
  logical :: TimeLoop_Succes
  logical :: TimeLoop_Error

  !// TimeStep Manager
  character(len=10) :: TUnit
  real(dp) :: TFinal
  real(dp) :: dTInit
  real(dp) :: dTRestart
  real(dp) :: dTMin
  real(dp) :: dTMax

  !// TimeStep Shared Logicals
  public :: OKFinalTimeStep
  public :: ERRORTimeStepTooSmall
 
contains

  subroutine Box_TimeLoop_Vars_New
    implicit none
    call Info_("Box_TimeLoop_Vars_New")

    call Box_TimeLoop_Vars_Zero

  end subroutine Box_TimeLoop_Vars_New

  !---

  subroutine Box_TimeLoop_Vars_Zero
    implicit none
    call Info_("Box_TimeLoop_Vars_Zero")

    TUnit = "none"
    TFinal = Zero
    dTInit = Zero
    dTRestart = Zero
    dTMin  = Zero
    dTMax  = Zero

  end subroutine Box_TimeLoop_Vars_Zero

  !---

  subroutine Box_TimeLoop_Vars_Delete
    implicit none
    call Info_("Box_TimeLoop_Vars_Delete")

    call Box_TimeLoop_Vars_Zero

  end subroutine Box_TimeLoop_Vars_Delete

end module M_Box_TimeLoop_Vars
