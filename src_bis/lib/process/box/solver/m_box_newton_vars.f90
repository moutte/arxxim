module M_Box_Newton_Vars

  !====================================================================
  ! NEWTON SOLVER parameterS 
  !====================================================================

  use M_Kinds
  use M_Trace
  
  use M_T_Nonlinear_Solver

  implicit none
  public

  character(len=30) :: cMethod

  !// Newton Numerical Parameters
  integer:: NewtMaxIts
  integer:: NewtIterMax ! if(iter>NewtIterMax) dT0=dT0*Time_Decrease
  integer:: NewtIterMin ! if(iter<NewtIterMin) dT0=dT0*Time_Increase

  real(dp):: NewtTOLF    ! convergence criterion on function values
  real(dp):: NewtTOLMIN  ! criterion for spurious convergence
  real(dp):: NewtTOLX    ! convergence criterion on dx

  logical::  bFinDif     ! finite difference jacobian compute 

  !// Debug Test Options
  logical:: DebNewt      ! Deb Newt
  logical:: DebJacob     ! Deb Jacob

  logical:: TestJacob    ! Test Jacob
  logical:: TestMax      ! Test Max

  type (T_Nonlinear_Solver_Param)  :: Newton_Param 
  type (T_Nonlinear_Solver_Result) :: Newton_Result

  integer:: Total_NewtonIter
  integer:: NewtonIter

contains
  
  subroutine Box_Newton_Vars_New    
    implicit none
    call Info_("Box_Newton_Vars_New")

    call Box_Newton_Vars_Zero

  end subroutine Box_Newton_Vars_New

  !---
  
  subroutine Box_Newton_Vars_Zero
    implicit none    
    call Info_("Box_Newton_Vars_Zero")
    
    cMethod="none"
    
    NewtMaxIts = 0
    NewtIterMax = 0
    NewtIterMin = 0
    !--
    NewtTOLF = Zero
    NewtTOLMIN = Zero
    NewtTOLX = Zero
    !
    bFinDif = .false.
    !
    DebNewt = .false.
    DebJacob = .false.
    TestJacob = .false.
    TestMax = .false.

    Total_NewtonIter = 0
    NewtonIter = 0
    
    call Nonlinear_Solver_Zero(Newton_Param, Newton_Result) 
    
  end subroutine Box_Newton_Vars_Zero

  !---
  
  subroutine Box_Newton_Vars_Delete
    implicit none
    call Info_("Box_Newton_Vars_Delete")

    call Box_Newton_Vars_Zero

  end subroutine Box_Newton_Vars_Delete

end module M_Box_Newton_Vars
