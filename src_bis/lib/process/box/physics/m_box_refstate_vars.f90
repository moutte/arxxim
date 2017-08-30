module M_Box_RefState_Vars

  !==================================================================================
  ! REFERENCE STATE VARIABLES
  !==================================================================================

  use M_Kinds
  use M_Trace
  use M_Box_Vars, only : nMk, nAq, nAs, nCp

  implicit none
  public

  real(dp), allocatable:: vG0rt_Aq(:)    !G0Rt of aqueous species in box
  real(dp), allocatable:: vG0rt_Mk(:)    !G0Rt of mineral species in box

  real(dp), allocatable:: vLnK_As(:)     !LnK of aqueous species in box
  real(dp), allocatable:: vLnK_Mk(:)     !LnK of mineral species in box

contains 

  subroutine Box_RefState_Vars_New
    use M_Box_Vars
    implicit none

    call Info_("Box_RefState_Vars_New")

    allocate(vG0rt_Aq(nAq))
    allocate(vG0rt_Mk(nMk))

    allocate(vLnK_As(nAs))
    allocate(vLnK_Mk(nMk))

    call Box_RefState_vars_Zero

  end subroutine Box_RefState_Vars_New

  !---

  subroutine Box_RefState_Vars_Zero
    implicit none

    call Info_("Box_RefState_Vars_Zero")
    
    vG0rt_Aq = Zero
    vG0rt_Mk = Zero

    vLnK_As = Zero
    vLnK_Mk = Zero

  end subroutine Box_RefState_Vars_Zero

  !---

  subroutine Box_RefState_Vars_Delete
    implicit none

    call Info_("Box_RefState_Vars_Del")
    
    call Box_RefState_vars_Zero

    deallocate(vG0rt_Aq)
    deallocate(vG0rt_Mk)

    deallocate(vLnK_As)
    deallocate(vLnK_Mk)

  end subroutine Box_RefState_Vars_Delete

end module M_Box_RefState_Vars
