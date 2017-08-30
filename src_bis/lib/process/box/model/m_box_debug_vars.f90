module M_Box_Debug_Vars
  
  use M_Kinds
  use M_Trace

  implicit none 
  public
  
  logical :: LDebug_Info  
  logical :: LDebug_Warning 
  !--
  logical :: LDebug_System     
  logical :: LDebug_Init  
  logical :: LDebug_Check
  logical :: LDebug_Calc       
  !--
  logical :: LDebug_KinThermo  
  logical :: LDebug_KinRate   
  logical :: LDebug_Gamma      
  logical :: LDebug_RefState   
  logical :: LDebug_Inject 
  logical :: LDebug_Outflow
  logical :: LDebug_Volume
  logical :: LDebug_TPUpdate 
  !--
  logical :: LDebug_Lumping
  !-- 
  logical :: LDebug_Entries    
  logical :: LDebug_Variables  
  logical :: LDebug_Residual   
  logical :: LDebug_Jacobian   
  !--
  logical :: LDebug_Newton 
  logical :: LDebug_TimeLoop
  logical :: LDebug_TimeStep
  
  contains 
 
  subroutine Box_Debug_Vars_New
    implicit none
    call Info_("Box_Debug_Vars_New")
    call Box_Debug_Vars_Zero
  
  end subroutine Box_Debug_Vars_New
  
  !---
  
  subroutine Box_Debug_Vars_Zero
    implicit none
    
    call Info_("Box_Debug_Vars_Zero")

    LDebug_Info= .false.  
    LDebug_Warning= .false. 
    !--
    LDebug_System= .false.    
    LDebug_Init= .false.       
    LDebug_Calc= .false.
    LDebug_Check= .false. 
    !--
    LDebug_KinThermo= .false.   
    LDebug_KinRate= .false.   
    LDebug_Gamma= .false.     
    LDebug_RefState= .false.   
    LDebug_Inject= .false.   
    LDebug_Outflow= .false. 
    LDebug_Volume= .false.    
    LDebug_Outflow= .false. 
    LDebug_TPUpdate= .false.  
    !--
    LDebug_Lumping= .false.
    !-- 
    LDebug_Entries= .false.  
    LDebug_Variables= .false.  
    LDebug_Residual= .false.  
    LDebug_Jacobian= .false. 
    !--
    LDebug_Newton= .false.     
    LDebug_TimeLoop= .false. 
    LDebug_TimeStep= .false. 

  end subroutine Box_Debug_Vars_Zero

  !---

  subroutine Box_Debug_Vars_Delete  
    implicit none
    call Info_("Box_Debug_Vars_Del")
    !---
    call Box_Debug_Vars_Zero
  end subroutine Box_Debug_Vars_Delete

  
end module M_Box_Debug_Vars
