module M_Box_TimeStorage

  !=======================================
  ! Storage of Time Step Solutions
  !---------------------------------------
  ! Time integrals are used for
  ! local time stepping strategies
  !=======================================

  use M_stockvar_KINXIM_BOX
  use M_Kinds,only: dp
  implicit none
  private
  
  public :: Box_TimeStorage_New
  public :: Box_TimeStorage_Del
  public :: Box_TimeStorage_Zero
  !---
  public :: Box_TimeStorage_Store_StartSolution
  public :: Box_TimeStorage_Store_ComputeSolution
  public :: Box_TimeStorage_Store_EndSolution
  
  integer :: BOX_MAX_NBTIMESTEP = 50000
  integer :: BOX_TIMESTORAGE = .true.
  
contains
  
  subroutine Box_TimeStorage_New
    use M_Box_Param_Vars, only : nAq
    implicit none

    LSTOCK = BOX_TIMESTORAGE
    
    if (LSTOCK) then
       call Init_Stockvar ( BOX_MAX_NBTIMESTEP, nAq )
       call Box_TimeStorage_Zero
    end if
    
  end subroutine Box_TimeStorage_New

  !---

  subroutine Box_TimeStorage_Del
    implicit none
    
    if (LSTOCK) then
       call Del_Stockvar 
    end if
    
  end subroutine Box_TimeStorage_Del
  
  !---
  
  subroutine Box_TimeStorage_Zero
    implicit none
    
    if (LSTOCK) then
       call Reset_Stockvar
    end if
    
  end subroutine Box_TimeStorage_Zero

  !---
  
  subroutine Box_TimeStorage_Store_StartSolution(Time, DTime)
    use M_Box_Vars, only : Vout, vMolF, nAq
    use M_Box_Newton_Vars, only : Newton_Result
    implicit none
    real(dp), intent(in) :: Time, Dtime
    real(dp) :: ErrF
    !--
    ErrF =  Newton_Result%ErrF
    
    !--- store values of output flux     
    if (LSTOCK) then
       call Set_Stockvar( vMolF(1:nAq), Vout, Time, Dtime, ErrF, 'S') 
    end if
    
  end subroutine Box_TimeStorage_Store_StartSolution
  
  !--

  subroutine Box_TimeStorage_Store_ComputeSolution(Time, Dtime)
    use M_Box_Vars, only : Vout, vMolF, nAq
    use M_Box_Newton_Vars, only : Newton_Result
    implicit none
    real(dp), intent(in) :: Time, Dtime
    real(dp) :: ErrF
    !--
    ErrF =  Newton_Result%ErrF

    !--- store values of output flux 
    if (LSTOCK) then
       call Set_Stockvar( vMolF(1:nAq), Vout, Time, Dtime, ErrF, 'C') 
    end if
    
  end subroutine Box_TimeStorage_Store_ComputeSolution

  !---
  
  subroutine Box_TimeStorage_Store_EndSolution(Time, DTime)
    use M_Box_Vars, only : Vout, vMolF, nAq
    use M_Box_Newton_Vars, only : Newton_Result
    implicit none
    real(dp), intent(in) :: Time, Dtime
    real(dp) :: ErrF
    !--
    ErrF =  Newton_Result%ErrF
    
    !--- store values of output flux 
    if (LSTOCK) then
       call Set_StockVar(vMolF(1:nAq), Vout, Time, Dtime, ErrF, 'E') 
    end if

  end subroutine Box_TimeStorage_Store_EndSolution


end module M_Box_TimeStorage
