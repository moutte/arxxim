module M_T_Nonlinear_Solver

  use M_Kinds
  use M_Trace
  use M_Info_Value
  implicit none
  
  private
 
  public:: T_Nonlinear_Solver_Param
  public:: T_NonLinear_Solver_Result
 
  !// NEWTON SOLVER parameterS
  type T_Nonlinear_Solver_Param
    character(len=30) :: cMethod
    real(dp) :: TolF
    real(dp) :: TolX
    logical  :: bFinDif
    integer  :: MaxIts     
  end type T_Nonlinear_Solver_Param
  
  !// NEWTON SOLVER resultS
  type T_Nonlinear_Solver_Result
    real(dp):: ErrF
    real(dp):: ErrX
    real(dp):: ErrG
    integer:: Nits
    logical:: Check
    integer:: iErr     
  end type T_Nonlinear_Solver_Result

  public :: Nonlinear_Solver_Zero
  public :: Nonlinear_Solver_Result_Info
  public :: Nonlinear_Solver_Param_Info
  
contains
  
  subroutine Nonlinear_Solver_Zero(Solver_Param, Solver_Result) 
    
    implicit none
    type(T_Nonlinear_Solver_Param)  :: Solver_Param
    type(T_NonLinear_Solver_Result) :: Solver_Result
    !---
    Solver_Param%cMethod = "none"
    Solver_Param%TolF = Zero
    Solver_Param%TolX = Zero
    Solver_Param%bFinDif = .false.
    Solver_Param%MaxIts = 0
    
    Solver_Result%ErrF = Zero
    Solver_Result%ErrX = Zero
    Solver_Result%ErrG = Zero
    Solver_Result%Nits = 0
    Solver_Result%Check = .false.
    Solver_Result%iErr = 0
    
  end subroutine Nonlinear_Solver_Zero
    
  !---
   
  subroutine Nonlinear_Solver_Param_Info(Solver_Param)

    implicit none
    type(T_NonLinear_Solver_Param) :: Solver_Param
    !---
    call Info_("Nonlinear_Solver_Param_Info")
    call Info_Value("Solver_Param cMethod", Solver_Param% cMethod)
    call Info_Value("Solver_Param TolF",    Solver_Param% TolF)
    call Info_Value("Solver_Param TolF",    Solver_Param% TolX)
    call Info_Value("Solver_Param MaxIts",  Solver_Param% MaxIts)
    call Info_Value("Solver_Param TolF",    Solver_Param% bFinDif)
    
  end subroutine Nonlinear_Solver_Param_Info
   
  !---
   
  subroutine Nonlinear_Solver_Result_Info(Solver_Result)

    implicit none
    type(T_NonLinear_Solver_Result) :: Solver_Result
    !---
    call Info_("Nonlinear_Solver_Result_Info")
    call Info_Value("Solver_Result ErrF", Solver_Result% ErrF)
    call Info_Value("Solver_Result ErrX", Solver_Result% ErrX)
    call Info_Value("Solver_Result ErrG", Solver_Result% ErrG)
    call Info_Value("Solver_Result Nits", Solver_Result% Nits)
    call Info_Value("Solver_Result Check",Solver_Result% Check)
    call Info_Value("Solver_Result iErr", Solver_Result% iErr)
    
  end subroutine Nonlinear_Solver_Result_Info
 
end module M_T_NonLinear_Solver
