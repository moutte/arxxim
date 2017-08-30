module M_TimeMng_Vars
  !===================================================
  ! Time Discretization Manager 
  ! Basic Manager for TimeLine and TimeStep 
  !===================================================
  use M_T_TimeLine
  use M_T_TimeStep
  use M_T_TimeStepMng

  use M_Kinds
  implicit none

  public
  
  !// private VARIABLES
  type (T_TimeLine):: TimeLine
  type (T_TimeStep):: TimeStep
  type (T_TimeStepMng) :: TimeStepMng
  
  logical:: ERRORTimeStepTooSmall
  logical:: OKFinalTimeStep

end module M_TimeMng_Vars
