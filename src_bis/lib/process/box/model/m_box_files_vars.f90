module M_Box_Files_Vars
  
  use M_Kinds
  use M_Trace

  implicit none 
  public
  
  !integer :: fBoxAct       ! aqueous species activities   
  !integer :: fBoxQsK       ! saturation index of minerals
  !integer :: fBoxGam       ! aqueous species activity coefficients
  
  integer :: fBoxEle       ! aqueous element mole numbers
  integer :: fBoxMnK       ! minerals amounts
  integer :: fBoxMol       ! aqueous species mole numbers
  
  contains 
 
  subroutine Box_Files_Vars_New
    implicit none
    call Info_("Box_Files_Vars_New")
    call Box_Files_Vars_Zero
  
  end subroutine Box_Files_Vars_New
  
  !---
  
  subroutine Box_Files_Vars_Zero
    implicit none
    
    call Info_("Box_Files_Vars_Zero")

    !fBoxAct = 0    
    !fBoxQsK = 0    
    !fBoxGam = 0
    
    fBoxEle = 0
    fBoxMnK = 0
    fBoxMol = 0
    
  end subroutine Box_Files_Vars_Zero

  !---

  subroutine Box_Files_Vars_Delete  
    implicit none
    call Info_("Box_Files_Vars_Del")
    !---
    call Box_Files_Vars_Zero
    
  end subroutine Box_Files_Vars_Delete
  
end module M_Box_Files_Vars

