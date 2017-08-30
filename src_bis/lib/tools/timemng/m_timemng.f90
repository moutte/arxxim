module M_TimeMng
  !===================================================
  ! Time Discretization Manager 
  ! Basic Manager for TimeLine and TimeStep 
  !===================================================
  use M_TimeMng_Vars
  use M_TimeMng_Base
  use M_TimeMng_Time
  use M_TimeMng_Control
  use M_Kinds
  use M_Kinds

  implicit none
  private

  !// Exported public typeS
  public:: T_TimeLine
  public:: T_TimeStep
  public:: T_TimeStepMng
  
  !// Exported public VARIABLES
  public:: ERRORTimeStepTooSmall
  public:: OKFinalTimeStep

  !// Exported Public Functions  :: Base
  public:: TimeMng_Zero
  public:: TimeMng_Set_TimeLine
  public:: TimeMng_Set_TimeStepMng
  public:: TimeMng_Set_DTRatio
  public:: TimeMng_Check
  public:: TimeMng_Reset_TimeStep
  public:: TimeMng_Info

  !// Exported Public Functions  :: Time
  public:: TimeMng_Get_Time
  public:: TimeMng_Update_Time
  public:: TimeMng_Info_TimeStep

  !// Exported Public Functions  :: Control
  public :: TimeMng_Control_TimeStep

end module M_TimeMng
