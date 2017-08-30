module M_Box_Outflow

  use M_Kinds
  use M_Trace

  use M_Box_Vars

  implicit none
  private

  public :: Box_Compute_Outflow

contains

  subroutine Box_Compute_Outflow

    use M_Box_Coupler_Vars
    implicit none

    call Info_('BOX_Compute_Outflow')

    if (LCouplerOutflow) then
       call Box_Compute_Outflow_Coupler
    else
       call Box_Compute_Outflow_System
    end if

  end subroutine Box_Compute_Outflow

  !---

  subroutine Box_Compute_Outflow_Coupler

    implicit none
    call Info_('BOX_Compute_Outflow_Coupler')

    !// Vout is given directly

  end subroutine Box_Compute_Outflow_Coupler

  !---

  subroutine Box_Compute_Outflow_System

    use M_Box_Debug_Vars, only : LDebug_Outflow
    implicit none

    call Info_('BOX_Compute_Outflow_System')

    Vout = Fout/VolF

    if (LDebug_Outflow) then
       write(*,*) "Darcy   Outflow Rate Fout = ", Fout
       write(*,*) "Reduced Outflow Rate Vout = ", Vout
    end if


  end subroutine Box_Compute_Outflow_System


end module M_Box_Outflow
