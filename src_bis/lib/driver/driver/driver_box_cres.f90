subroutine Driver_Box_Cres
  
  use M_Trace
  use M_API_KINXIM_BOX 
  implicit none
  
  !---
  
  call Init_API_KINXIM_INTERACT
  call Compute_API_KINXIM_noargs
  call Output_API_KINXIM
  call End_API_KINXIM
  
end subroutine Driver_Box_Cres
