module M_Numeric_Nonlinear_Solver_Typed
  !--------------------------------------------------------------------------
  ! Purpose   : General Interface to Nonlinear Newton Type Solvers 
  ! With Types T_Newt_Param and T_Newt_Result
  !--------------------------------------------------------------------------

  use M_T_Nonlinear_Solver
  implicit none
  private
  !
  public:: Nonlinear_Solver_Typed
  public:: T_Nonlinear_Solver_Param
  public:: T_Nonlinear_Solver_Result
  !
contains

  !---

  subroutine Nonlinear_Solver_Typed( &
  & vX,        & !inout= the initial guess, and the root returned
  & vIsPlus,   & !in=    if(vIsPlus(i)) then vX(i) must be >Zero
  & Residual,  & !in= external residual function
  & Jacobian,  & !in= external jacobian function
  & Converge,  & !in= external converge test function
  & NLS_Param, & !in= input parameters
  & NLS_Result & !out= output results
  & )

    !------------------------------------------------
    ! Purpose : Nonlinear Solver Wrapper 
    !------------------------------------------------

    use M_Kinds
    use M_Numeric_Nonlinear_Solver
    implicit none

    real(dp),dimension(:), intent(inout):: vX
    logical, dimension(:), intent(in)   :: vIsPlus

    interface
      function Residual(v)
        use M_Kinds
        implicit none
        real(dp),dimension(:),intent(in):: v
        real(dp),dimension(size(v))     :: Residual
      end function Residual

      subroutine Jacobian(v,t)
        use M_Kinds
        implicit none
        real(dp),dimension(:),              intent(in) :: v
        real(dp),dimension(size(v),size(v)),intent(out):: t
      end subroutine Jacobian

      logical function Converge(v,tolv)
        use M_Kinds
        implicit none
        real(dp),intent(in):: v(:),tolv(:)
      end function Converge
    end interface

    type(T_Nonlinear_Solver_Param)::  NLS_Param
    type(T_Nonlinear_Solver_Result):: NLS_Result

    !//------------ call NONLINEAR SOLVER    
    call Nonlinear_Solver( &
    
    & NLS_Param% cMethod, & !input= name of the method
    
    & vX,           & !inout= the initial guess, and the root returned
    & vIsPlus,      & !in= implements positivity constraint
    & Residual,     & !in= external residual function
    & Jacobian,     & !in= external jacobian function
    & Converge,     & !in= external converge test function
    
    & NLS_Param% TolF,     & !in=    convergence criterion on function values 
    & NLS_Param% TolX,     & !in=    convergence criterion on dx
    & NLS_Param% bFinDif,  & !in=    use numeric jacobian
    & NLS_Param% MaxIts,   & !in=    maximum number of iterations
    !---
    & NLS_Result% ErrF,    & !out=   MAXVAL(ABS(fVec(:)))
    & NLS_Result% ErrX,    & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & NLS_Result% ErrG,    & !out=   value of the gradient ???
    & NLS_Result% Nits,    & !out=   number of iterations
    & NLS_Result% Check,   & !out=   if Check, should check convergence ???
    & NLS_Result% iErr     & !out=   error code
    & )

  end subroutine Nonlinear_Solver_Typed

end module M_Numeric_Nonlinear_Solver_Typed
