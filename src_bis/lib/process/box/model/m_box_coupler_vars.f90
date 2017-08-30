module M_Box_Coupler_Vars

  use M_Kinds
  !ltuse M_Trace
  
  use M_Box_Param_Vars, only : nAq, nMk
  
  implicit none
  public
  
  !==// WORKING ARRAYS
  logical :: LCouplerInject 
  logical :: LCouplerOutflow
  logical :: LCouplerActive
  integer :: IERRORCHEMISTRY
  
contains
  
  subroutine Box_Coupler_Vars_New
    implicit none
    !lt call Info_("Box_Coupler_Vars_New")
    !---
    call Box_Coupler_Vars_Zero
  end subroutine Box_Coupler_Vars_New
  
  !---
  
  subroutine Box_Coupler_Vars_Zero
    implicit none
    !ltcall Info_("Box_Coupler_Vars_Zero")
    !---
    LCouplerInject  = .false.
    LCouplerOutflow = .false.
    LCouplerActive  = .false.
    IerrorChemistry = 0

  end subroutine Box_Coupler_Vars_Zero

  !---

  subroutine Box_Coupler_Vars_Delete 
    implicit none
    !ltcall Info_("Box_Coupler_Vars_Delete")
    !---
    call Box_Coupler_Vars_Zero
   
  end subroutine Box_Coupler_Vars_Delete

end module M_Box_Coupler_Vars
