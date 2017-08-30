!! *******************************************************************
!! * File name:    M_Box_Calc
!! * Purpose:      Local Chemical Computation
!! * Author:       Anthony Michel (anthony.michel@ifp.fr)
!! * Created:      2008
!! ********************************************************************
module M_Box_Calc

  use M_Kinds
  use M_Trace
  implicit none
  private

  !// public members
  public :: BOX_Calc_Init
  public :: BOX_Calc_Start
  !--
  public :: BOX_Calc_Restart
  !--
  public :: Box_Calc_Prepare
  public :: Box_Calc_Compute
  public :: Box_Calc_Finalize
  !--
  public :: BOX_Calc_End

  !==================================================
  !  Init-Loop
  !        BOX_Calc_Init
  !        BOX_Calc_Start
  !
  ! Compute-Loop
  !        (( Box_Calc_Restart ))
  !        Box_Calc_Prepare
  !        Box_Calc_Compute
  !        Box_Calc_Finalize
  !
  !  End-Loop
  !        Box_Calc_End
  !==================================================

contains

  subroutine Box_Calc_Init
    !==================================================
    ! Purpose   : Init Static Using Input File 
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================
    use M_Box_Debug_Vars
    use M_Box_Thermo_Vars
    use M_Box_RefState_Vars
    use M_Box_System_Vars
    use M_Box_Kinetics_Vars
    use M_Box_Load_Vars
    use M_Box_Residual_Vars
    use M_Box_TimeLoop_Vars
    use M_Box_Newton_Vars
    use M_Box_Solver_Vars
    use M_Box_Init_Vars
    use M_Box_Vars
    use M_Box_Coupler_Vars
    use M_Box_TimeStorage
    use M_Box_Lumping_Vars
    use M_Box_Read
    use M_Box_Init
    use M_Box_TP_Update

    implicit none
    logical :: Ok

    call Info_("Box_Init")
    
    call Box_Debug_Vars_New
    call Box_Read_Debug(Ok)
    if (LDebug_Info) LInfo = .true.
    if (LDebug_Warning) LWarning = .true.

    !// Init System and Fluid Compositions
    call Box_Init_System    
    
    !// Init Box Solver Variables
    call Box_Init_FluidInject(Ok)
    call Box_Init_FluidBox(Ok)
   
    !// Reset Basis
    call Box_Init_Reset_Basis    
	
    !// Init Box Variables
    call Box_Vars_New
    call Box_RefState_Vars_New
    call Box_Thermo_Vars_New
    call Box_Kinetics_Vars_New
    call Box_Solver_Vars_New
    call Box_Newton_Vars_New
    call Box_Load_Vars_New
    call Box_Residual_Vars_New
    call Box_Coupler_Vars_New
    call Box_Lumping_Vars_New
    
    !// Init Time Storage
    call Box_TimeStorage_New
    
    !// Init Dynamic Resolution Parameters
    call Box_Read_Time(Ok)
    call Box_Read_Numeric(Ok)
    call Box_Read_Reactor(Ok)
    call Box_Read_Lumping(Ok)
   
    !// Init Reactor Properties
    call Box_Init_vMolF_FluidInject
    call Box_Init_vMolF_FluidBox

    !// Init Rock Properties
    call Box_Init_vMolK_Rock
    call Box_Init_Texture_Rock
    
    !// Update TP Model Properties if necessary
    call BOX_TP_Update
    
    !// Debug Info
    call Box_Init_Info    

  end subroutine Box_Calc_Init

 
  !----

  subroutine Box_Calc_Start
    !==================================================
    ! Purpose   : Start during the initial phase
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================
    use M_TimeMng
    use M_T_TimeLine
    use M_T_TimeStepMng
    use M_Box_Debug_Vars, only : LDebug_TimeLoop
    use M_Box_Info
    
    use M_Box_TimeLoop_Vars
    use M_Box_Load_Vars
    use M_Box_Vars
    use M_Box_Residual_Vars

    implicit none

    type(T_TimeLine)    :: TimeLine
    type(T_TimeStepMng) :: TimeStepMng

    call Info_("Box_Start")
    call Box_System_Info()

    !// Set TimeLoop
    call TimeLine_Init(TimeLine, Zero, TFinal)
    call TimeMng_Set_TimeLine(TimeLine)

    call TimeStepMng_Init(TimeStepMng, dTInit, dTMin, dTMax)
    call TimeMng_Set_TimeStepMng(TimeStepMng)

    call TimeMng_Reset_TimeStep    
    if (LDebug_TimeLoop) call TimeMng_Info  
       
    call TimeMng_Get_Time(Time,dTime)
    
  end subroutine Box_Calc_Start

  !----

  subroutine Box_Calc_Restart
    !==================================================
    ! Purpose   : ReStart during the compute phase
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================

    use M_TimeMng
    use M_T_TimeLine
    use M_T_TimeStepMng
    use M_Box_TimeLoop_Vars
    use M_Box_Debug_Vars, only : LDebug_TimeLoop
    use M_Box_Vars
    use M_Box_TP_Update
    implicit none
    logical :: LResetdTime = .true.

    type(T_TimeLine)    :: TimeLine
    type(T_TimeStepMng) :: TimeStepMng

    call Info_("Box_Restart")

    ! TOdo : implement restart if necessary
    call BOX_TP_Update

    !// Test continuation 
    !write(*,*) "Tfinal", Tfinal, Tfinal/86400
    call TimeLine_Init(TimeLine, Time, TFinal)
    call TimeMng_Set_TimeLine(TimeLine)
    
    !// Time Stepping    
    LResetdTime = .true.
    dTRestart = dTime
    if ( LResetdTime ) dTRestart = dTInit  
    
    call TimeStepMng_Init(TimeStepMng, dTRestart, dTMin, dTMax)
    call TimeMng_Set_TimeStepMng(TimeStepMng)

    call TimeMng_Set_TimeLine(TimeLine)

    call TimeMng_Reset_TimeStep
    if (LDebug_TimeLoop) call TimeMng_Info  

  end subroutine Box_Calc_Restart

  !----

  subroutine Box_Calc_Prepare
    !==================================================
    ! Purpose   : Prepare the compute phase
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================

    use M_Box_Prepare
    use M_Box_Test
    use M_Box_Newton

    implicit none

    call Info_("Box_Calc_Prepare")

    !// Init Box Solver Variables
    call Box_Prepare_Guess
    call Box_Newton_Init
    
    call Warning_("---------------------------------------------------")
    call Warning_("Check Residual At T = Tinit")
    call Warning_("---------------------------------------------------")
    call Box_Load_Check_Residual
    
  end subroutine Box_Calc_Prepare

  !---

  subroutine Box_Calc_Compute
    !==================================================
    ! Purpose   : Apply dynamic computation
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================
    use M_Box_TimeLoop
    implicit none

    !write(*,*) "********************************************************"
    !write(*,*) "            Box_Calc_Compute                            "
    !write(*,*) "********************************************************"
    call Info_("Box_Calc_Compute")
    call Box_TimeLoop_Solve
    
  end subroutine Box_Calc_Compute

  !---- 

  subroutine Box_Calc_Finalize
    !==================================================
    ! Purpose   : Finalize the compute phase
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================
    use M_Box_TimeLoop_Vars
    use M_Box_Load_Vars
    use M_Box_Residual_Vars

    !// Clean Box Solver Variables    
    call Box_Load_Vars_Zero
    call Box_Residual_Vars_Zero
    
  end subroutine Box_Calc_Finalize

 !----

  subroutine Box_Calc_End
    !==================================================
    ! Purpose   : Clear all at the end phase
    ! Author    : Anthony Michel (anthony.michel@ifp.fr)
    ! Update    : 2008
    !==================================================
    use M_Box_RefState_Vars
    use M_Box_Debug_Vars
    use M_Box_Thermo_Vars
    use M_Box_System_Vars
    use M_Box_Kinetics_Vars
    use M_Box_Load_Vars
    use M_Box_Residual_Vars
    use M_Box_TimeLoop_Vars
    use M_Box_Newton_Vars
    use M_Box_Solver_Vars
    use M_Box_Init_Vars
    use M_Box_Coupler_Vars
    use M_Box_TimeStorage
    use M_Box_Lumping_Vars
    use M_Box_Vars

    implicit none

    call Info_("Box_End")

    !// Delete Box Solver Variables
    call Box_Vars_Delete
    call Box_System_Vars_Delete
    call Box_RefState_Vars_Delete
    call Box_Thermo_Vars_Delete
    call Box_Kinetics_Vars_Delete
    call Box_Load_Vars_Delete
    call Box_Residual_Vars_Delete
    call Box_Init_Vars_Delete
    call Box_Solver_Vars_Delete
    call Box_Newton_Vars_Delete
    call Box_Coupler_Vars_Delete
    call Box_Lumping_Vars_Delete

    !// Delete Time Storage
    call Box_TimeStorage_Del

    call Box_Debug_Vars_Delete

  end subroutine Box_Calc_End


end module M_Box_Calc
