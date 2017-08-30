module M_Box_Solver_Vars

  !=====================================================
  ! BOX SOLVER parameterS AND VARIABLES
  !=====================================================
  
  use M_Kinds
  use M_Trace

  use M_Box_Param_Vars
  
  implicit none
  public
  
  !==// GENERAL OPTIONS
  
  logical:: LogForFluid= .true.  !.true. !-> works on log(MoleNr) for AqSpecies
  logical:: LogForMin= .false.   !.true. !-> works on log(MoleNr) for MkSpecies 
    
  logical:: Implicit_Surface= .true.      !.true. !-> use implicit surface 
  logical:: Implicit_ActivFactor= .false. !.true. !-> implicit formula for activities in catalyse factor
  
  real(dp) :: ThetaRate= 1.0_dp  ! 1. = implicit, 0.5 = trapeze, 0. = explicit

  integer:: fSavTime= 0 !file for saving results at (nearly) fixed time intervals (dTSav)
  integer:: fSavRate= 0 !file for details on the rate parameters of minerals

  !==// TIME MNG OPTIONS

  character(len=4):: TUnit      !"DAY","YEAR","HOUR"
  real(dp)        :: TimeFactor != the current time unit, in seconds
  !
  real(dp):: TInit     ! end of simulation
  real(dp):: TFinal    ! end of simulation
  !
  real(dp):: dTInit    ! minimal time step
  real(dp):: dTmin     ! minimal time step
  real(dp):: dTMax     ! maximal time step
  
  real(dp):: Time_Decrease= 0.5_dp  ! timestep decrease factor
  real(dp):: Time_Increase= 2.0_dp  ! timestep increase factor  

  !==// resultS

  integer:: Dynam_nStep
  integer:: Dynam_nNewtIter
  integer:: Dynam_nTotalNewtIter

  contains 
 
  subroutine Box_Solver_Vars_New
    implicit none
    call Info_("Box_Solver_Vars_New")
    call Box_Solver_Vars_Zero
  
  end subroutine Box_Solver_Vars_New
  
  !---
  
  subroutine Box_Solver_Vars_Zero
    implicit none
    call Info_("Box_Solver_Vars_Zero")
    
  end subroutine Box_Solver_Vars_Zero

  !---

  subroutine Box_Solver_Vars_Delete  
    implicit none
    call Info_("Box_Solver_Vars_Del")
    !---
    call Box_Solver_Vars_Zero
  end subroutine Box_Solver_Vars_Delete

end module M_Box_Solver_Vars
