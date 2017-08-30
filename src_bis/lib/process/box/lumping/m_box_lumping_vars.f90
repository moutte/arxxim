module M_Box_Lumping_Vars

  !=====================================================
  ! BOX LUMPING VARIABLES
  !=====================================================

  use M_Kinds
  use M_Trace

  use M_T_Lumping_Map
  use M_Box_Param_Vars, only : nCp, nAq, nMk, nSpc
 
  implicit none
  public 
  
  !// Lumping Map
  integer, parameter :: MaxSize_LumpMap = 100   
  type(T_Lumping_Map) , allocatable :: vLumpMap(:)
  integer :: Size_LumpMap
  character(len=20) :: CLump_OutputName
  
  !// C-LUMP Component  
  integer :: iCons_H2O
  integer :: iCons_CLUMP
  logical :: LOption_Lumping

contains

  !---

  subroutine Box_Lumping_Vars_New
    implicit none
    !--
    call Info_("Box_Lumping_Vars_New")

    allocate(vLumpMap(MaxSize_LumpMap))    

    call Box_Lumping_Vars_Zero

   
  end subroutine Box_Lumping_Vars_New

  !---

  subroutine Box_Lumping_Vars_Zero
    implicit none
    integer :: iMap
    !--
    call Info_("Box_Lumping_Vars_Zero")
    
    do iMap = 1, size(vLumpMap)
       call Lumping_Map_Zero(vLumpMap(iMap)) 
    end do
    Size_LumpMap = 0
    iCons_CLUMP = 0
    iCons_H2O = 0
    CLump_OutputName = 'none'  

  end subroutine Box_Lumping_Vars_Zero

  !---

  subroutine Box_Lumping_Vars_Delete
    implicit none
    !--
    call Info_("Box_Lumping_Vars_Delete")

    call Box_Lumping_Vars_Zero
    
    deallocate(vLumpMap)

  end subroutine Box_Lumping_Vars_Delete

end module M_Box_Lumping_Vars
