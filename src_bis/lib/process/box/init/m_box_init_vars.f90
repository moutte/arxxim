module M_Box_Init_Vars

  !=====================================================
  ! BOX INIT VARIABLES
  !=====================================================

  use M_Kinds
  use M_Trace

  use M_Box_Param_Vars, only : nCp, nAq, nMk, nSpc

  implicit none
  public 

  public :: nCp
  public :: nAq
  public :: nMk
  public :: nSpc
  
  real(dp), allocatable :: vMolCpn_FluidSpc(:)
  real(dp), allocatable :: vMolCpn_FluidBox(:)
  real(dp), allocatable :: vMolCpn_FluidInject(:)

  real(dp), allocatable :: vMolSpc_FluidSpc(:)
  real(dp), allocatable :: vMolSpc_FluidBox(:)
  real(dp), allocatable :: vMolSpc_FluidInject(:)

contains

  !---

  subroutine Box_Init_Vars_New
    implicit none
    !--
    call Info_("Box_Init_Vars_New")

    allocate(vMolCpn_FluidSpc(nCp))
    allocate(vMolCpn_FluidBox(nCp))
    allocate(vMolCpn_FluidInject(nCp))

    allocate(vMolSpc_FluidSpc(nSpc))   
    allocate(vMolSpc_FluidBox(nSpc))    
    allocate(vMolSpc_FluidInject(nSpc))

    call Box_Init_Vars_Zero

   
  end subroutine Box_Init_Vars_New

  !---

  subroutine Box_Init_Vars_Zero
    implicit none
    !--
    call Info_("Box_Init_Vars_Zero")
    
    vMolCpn_FluidSpc= Zero
    vMolCpn_FluidBox= Zero
    vMolCpn_FluidInject= Zero
    
    vMolSpc_FluidSpc= Zero    
    vMolSpc_FluidBox= Zero    
    vMolSpc_FluidInject= Zero

  end subroutine Box_Init_Vars_Zero

  !---

  subroutine Box_Init_Vars_Delete
    implicit none
    !--
    call Info_("Box_Init_Vars_Delete")

    call Box_Init_Vars_Zero
    
    deallocate(vMolCpn_FluidSpc)
    deallocate(vMolCpn_FluidBox)
    deallocate(vMolCpn_FluidInject)

    deallocate(vMolSpc_FluidSpc)
    deallocate(vMolSpc_FluidBox)
    deallocate(vMolSpc_FluidInject)

  end subroutine Box_Init_Vars_Delete


end module M_Box_Init_Vars
