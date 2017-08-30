module M_Box_Tools

  !====================================================
  ! Tools to Operate on Box Model Variables
  !====================================================
  use M_Trace
  use M_Kinds
  use M_TimeMng
  use M_Box_TimeStorage
    
  implicit none
  private

  public ::Init_Solution
  public ::Solve_Newton
  public ::Prepare_Newton
  public ::Check_Solution
  public ::Control_TimeStep
  public ::Control_Final_TimeStep
  public ::Update_Time 
  public ::Update_Solution
  public ::Recover_Solution
  public ::Save_Result
  public ::Close_Solution

contains

  !---

  subroutine Init_Solution
    !==============================================
    ! Init Solution and Storage Variables and Files
    !==============================================
    use M_BOX_COUPLER_VARS
    use M_Box_Vars
    use M_Box_Info
    use M_Box_Entries
    use M_Box_Load_Vars
    use M_Box_Files_Vars
    use M_Box_Files
    implicit none
    call Info_("Init_Solution")
    call Box_Compute_Entries(vX)
    !call Box_Info(100)
        
    !// Clear TimeStorage
    call Box_TimeStorage_Zero
    
    !// Store the time step solution 
    call Box_TimeStorage_Store_StartSolution(Time, DTime)

    !// Init Output File Storage
    call Box_Files_Vars_New
     if (.not.(LCouplerActive)) call Box_Files_Init

  end subroutine Init_Solution

  !---
  
  subroutine Close_Solution
    !==============================================
    ! Close Storage Variables and Files
    !==============================================
    use M_Box_Vars
    use M_Box_Files_Vars
    use M_Box_Files
    !
    implicit none
    !// Store the time step solution 
    call Box_TimeStorage_Store_EndSolution(Time, DTime)
    
    !// Close the solution output file
    call Box_Files_Close
    call Box_Files_Vars_Delete
  
  end subroutine Close_Solution

  !--

  subroutine Control_TimeStep(OKTimeStep)
    !==============================================
    ! CONTROL The time Step
    !==============================================
    use M_Box_Vars
    implicit none
    logical, intent(in) :: OKTimeStep
    !
    call Info_("Control_TimeStep") 
    !write(*,*) "Control Time", Time, dTime
    !write(*,*) "TimeOK", OKTimeStep
    call TimeMng_Control_TimeStep(OKTimeStep)
    !write(*,*) "TimeOK", OKTimeStep
    call TimeMng_Get_Time(Time,dTime)

  end subroutine Control_TimeStep

  !---

  subroutine Control_Final_TimeStep
    !==============================================
    ! Control The time Step at the last time step
    !==============================================
    use M_Box_Vars
    !
    implicit none
    !
    call Info_("Control_Final_TimeStep") 
    call TimeMng_Get_Time(Time,dTime)
    
  end subroutine Control_Final_TimeStep

  !---  

  subroutine Prepare_Newton
    !==============================================
    ! Compute Static Explicit Variables 
    !==============================================
    use M_Box_Load_Vars, only : vX
    use M_Box_Prepare
    !---
    implicit none
    !
    call Info_("Prepare Newton_Solve")    
    call BOX_Prepare_TimeStep
    call BOX_Prepare_GammaStep(vX)
    !OKNewton = .true.    

  end subroutine Prepare_Newton
  
  !---

  subroutine Solve_Newton(OKNewton)
    !==============================================
    ! Solve Implicit Nonlinear Problem with Newton 
    !==============================================
    use M_Box_Load_Vars, only : vX
    use M_Box_Newton
    !---
    implicit none
    logical :: OKNewton
    !
    call Info_("Compute Newton_Solve")    
    call Box_Newton_Solve(vX,  OKNewton)
    !OKNewton = .true.    

  end subroutine Solve_Newton

  !---

  subroutine Check_Solution(OKSolution)
    !==============================================
    ! CONTROL The time Step
    !==============================================
    use M_Box_Entries
    use M_Box_Load_Vars, only : vX
    implicit none
    logical :: OKSolution  
    real(dp):: Time, DTime
    !
    call Info_("Check_Solution")    
    call BOX_Check_Variables(vX, OKSolution)
    
    call TimeMng_Get_Time(Time, DTime)
    !OKSolution = .true.
    
  end subroutine Check_Solution

  !---

  subroutine Recover_Solution

    use M_Box_Load_Vars, only : vX
    use M_Box_Entries
    implicit none

    call Box_Recover_Variables(vX)

  end subroutine Recover_Solution

  !---
  
  subroutine Update_Time
    !==============================================
    ! Update Time for the next Time Step
    !==============================================
    implicit none
    !
    call Info_("Update_Solution")
    call TimeMng_Update_Time
    
  end subroutine Update_Time

  !---
  
  subroutine Update_Solution
    !==============================================
    ! Update The Solution for the next Time Step
    !==============================================
    
    use M_Box_Solver_Vars,   only : LogForMin
    use M_Box_Load_Vars, only : vX
    use M_Box_Residual_Vars
    use M_Box_Vars
    use M_Box_Utils
    use M_Box_Entries
    use M_Box_Kinetics_Vars
    use M_Box_Kinetics
    use M_echomat
    use M_Box_Debug_Vars, only : LDebug_TimeLoop
    
    implicit none
    integer :: iMk
    
    call Box_Compute_Variables(vX)
    
    !// Aqueous Phase Update -------

    if (LDebug_TimeLoop) then 
       call Info_("=======================================================================")
       call echovec_real(vMolf, size(vMolf),     'Initial Value vMolf', 'T')
       call echovec_real(vXf,   size(vXf),       'Solution vXf',        'T')
       call echovec_real(vXf-vMolf, size(vMolf), 'Difference',          'T')
       call Info_("=======================================================================")
    end if
    vMolf(1:nAq) = vXf(1:nAq)

    !// Mineral Phase Update -------

    if (LDebug_TimeLoop) then 
       call echovec_real(vMolK, size(vMolK),       'Initial Value vMolK', 'T')
       call echovec_real(vXm,   size(vXm),         'Solution vXm',        'T')
       call echovec_real(vXm(1:nMk) - vMolK(1:nMk), nMk, 'Difference',          'T')
    end if
    
    !// Adjust Minimal Mineral Amounts
    do iMk = 1, nMk
       if (  vXm(iMk) < vKinMinim(iMk) ) vXm(iMk) = vKinMinim(iMk)
    end do
    
    vMolK(1:nMk) = vXm(1:nMk)

    !// Store Explicit Mineral Results -------
    
    vKinRateK(1:nMk) = vVm(1:nMk)    

    call Info_("=======================================================================")
    
  end subroutine Update_Solution

  !---

  subroutine Save_Result

    !==============================================
    ! Save The Solution of the TimeStep
    !==============================================
    use M_BOX_COUPLER_VARS
    use M_Box_Info    
    use M_Box_Files
    implicit none
    real(dp):: Time, DTime
    !
    call Info_("Save_Result")
    !call Box_Info(100)
    
    call TimeMng_Get_Time(Time, DTime)

    !// Store the time step solution 
    call Box_TimeStorage_Store_ComputeSolution(Time, DTime)
    
    !// Save Output File
    if (.not.(LCouplerActive)) then
       call Box_Files_ShoMin
       call Box_Files_ShoAqu
       call Box_Files_ShoEle
    end if
    
  end subroutine Save_Result

end module M_Box_Tools
