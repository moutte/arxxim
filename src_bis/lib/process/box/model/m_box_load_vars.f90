module M_Box_Load_Vars

  use M_Kinds
  use M_Trace
  
  use M_Box_Param_Vars, only : nAq, nMk
  
  implicit none
  public
  
  !==// WORKING ARRAYS
  integer :: nX
  real(dp), allocatable :: vX(:)
  real(dp), allocatable :: vR(:)
  real(dp), allocatable :: tJac(:,:)
  
contains
  
  subroutine Box_Load_Vars_New
    implicit none
    call Info_("Box_Load_Vars_New")
    !---
    nX = nAq + nMk
    allocate (vX(nX))
    allocate (vR(nX))
    allocate (tJac(nX, nX))
    
    call Box_Load_Vars_Zero
  
  end subroutine Box_Load_Vars_New
  
  !---
  
  subroutine Box_Load_Vars_Zero
    implicit none
    call Info_("Box_Load_Vars_Zero")
    !---
    vX(:)= Zero
    vR(:)= Zero
    tJac(:,:)= Zero
    
  end subroutine Box_Load_Vars_Zero

  !---

  subroutine Box_Load_Vars_Delete 
    implicit none
    call Info_("Box_Load_Vars_Delete")
    !---
    call Box_Load_Vars_Zero
    deallocate (vX)
    deallocate (vR)
    deallocate (tJac)
    nX= 0
  end subroutine Box_Load_Vars_Delete

end module M_Box_Load_Vars
