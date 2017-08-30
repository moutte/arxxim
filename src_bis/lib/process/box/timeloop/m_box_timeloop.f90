module M_Box_TimeLoop

  !===================================================
  ! Purpose : Box Dynamic Time Solver Algorihm 
  !===================================================

  use M_Trace
  use M_Box_Tools
  use M_Kinds
  use M_Box_TimeLoop_Vars
  use M_BOX_COUPLER_VARS
     
  implicit none
  private

  !// public functionS
  interface Box_TimeLoop_Solve
     module procedure TimeLoop_Solve
  end interface

  public :: Box_TimeLoop_Solve

  !// private functionS

  private :: TimeLoop_Solve

  private :: TimeLoop_Init
  private :: TimeLoop_Compute
  private :: TimeLoop_End 
  private :: TimeLoop_Info_Step

contains

  !---

  subroutine TimeLoop_Solve
    !===================================================
    ! General Facade
    !===================================================
    implicit none
    !//
    call TimeLoop_Init
    call TimeLoop_Compute
    call TimeLoop_End

  end subroutine TimeLoop_Solve

  !---

  subroutine TimeLoop_Init()
    !===================================================
    ! Init The TimeLoop : T = TInit
    !===================================================
    !---
    implicit none

    call Info_("Init TimeLoop")

    !// Init Time Loop
    call Init_Solution
    iTimeStep= 0
    iTryStep = 0
    call TimeLoop_Info_Step("Init")
    
    !// Init TimeLoop Status
    TimeLoop_Succes = .false.
    TimeLoop_Error  = .false.

  end subroutine TimeLoop_Init

  !---

  subroutine TimeLoop_Compute
    !===================================================
    ! Compute The TimeLoop : T = TInit -> TFinal
    ! --------------------------------------------------
    ! High Level Exception Handling and TimeStep Control
    !===================================================
    use M_TimeMng

    implicit none
    logical :: OKNewton
    logical :: OKSolution
    logical :: OKTimeStep
    logical :: OKLastTimeStep
    logical :: OKTryNewStep
    ! 
    call Info_("Compute TimeLoop")

    OKTimeStep=.true.
    !
    StopTimeLoop = .false.
    IerrorChemistry = 0
    TimeLoop: do while ( .not.(StopTimeLoop) ) 

       ! -------------- INIT TIMESTEP 
       call Info_("Init Time Step")
       
       if(OKTimeStep) then
          iTimeStep= iTimeStep + 1
          iTryStep = 0
       end if
       
       call Control_TimeStep(OKTimeStep)
       OKTryNewStep= (.not.ERRORTimeStepTooSmall) 
       
       ! -------------- COMPUTE TIMESTEP 
       OKTimeStep=.false.        
       if (OKTryNewStep) then

          call Info_("Try New Step")          
          iTryStep= iTryStep + 1
          call TimeLoop_Info_Step("Compute")  
          
          ! ------- SOLVE NEWTON
          call Prepare_Newton
          call Solve_Newton(OKNewton)
          
          ! ------- CHECK SOLUTION
          OKSolution= .false.
          if (OKNewton) then
             call Check_Solution(OKSolution)
          end if

          OKTimeStep = OKSolution

       end if

       ! -------------- FINISH TIMESTEP
       call Info_("Finish TimeStep")

       !// Update TimeStep Solution
       if (OKTimeStep) then 
          !call TimeLoop_Info_Step("Solution")
          call Update_Time
          call Update_Solution
       else
          call Recover_Solution
       end if

       !// Save TimeStep Results
       if(OKTimeStep) then 
          call Save_Result
       end if

       !// Compute TimeLoop Status 
       TimeLoop_Succes= (OKFinalTimeStep.and.OKTimeStep)
       TimeLoop_Error= ERRORTimeStepTooSmall

       StopTimeLoop= (TimeLoop_Succes .or. TimeLoop_Error .or. (IerrorChemistry>0))

    end do TimeLoop

  end subroutine TimeLoop_Compute

  !---

  subroutine TimeLoop_End
    !====================================================
    ! Finalize Time Loop : T = TFinal of Error
    !====================================================
    implicit none
    !
    call Info_("End TimeLoop")
    !
    call Control_Final_TimeStep
    !
    if (TimeLoop_Succes) then
       !
       ! End TimeLoop 
       ! Compute the final solution properties
       !
       iTimeStep= iTimeStep + 1
       iTryStep = 0
       !
       call Info_("Advance to Final Step")          
       call TimeLoop_Info_Step("Finish")  
       !
       call Prepare_Newton
       call Save_Result
       !
       call TimeLoop_Info_Step("End") 
       !! call Close_Solution
       !
       call Info_("OK, TimeLoop end without error") 
    end if

    if (TimeLoop_Error) then
       call TimeLoop_Info_Step("Error")
       call Close_Solution
       call Fatal_("ERROR, TimeLoop stop before TFinal")
    end if

    call Close_Solution
    
  end subroutine TimeLoop_End

  !---
  
  subroutine TimeLoop_Info_Step(Str)
    !==============================================
    ! print Info Message for Time Step
    !==============================================
    use M_TimeMng, only : TimeMng_Get_Time
    use M_Box_Debug_Vars
    implicit none
    character(len=*) :: Str
    character(len=255):: Message
    !
    real(dp) :: Time, DTime
    character(len=10):: TAG
    !
    call Info_("Info_TimeStep")
    
    if (LDebug_TimeStep) then
       TAG = "["//trim(Str)//"]"

       call TimeMng_Get_Time(Time, DTime)

       write(Message,'(A10, A,I8, A, I4,2(A, G12.4))') &
            TAG, &
            "  TIMESTEP= ", iTimeStep,&
            ", TRYSTEP= ", iTryStep, &
            ", T= ", Time,&
            ", DT= ", DTime

       call Message_(Message)
    end if

  end subroutine TimeLoop_Info_Step

  


end module M_Box_TimeLoop
