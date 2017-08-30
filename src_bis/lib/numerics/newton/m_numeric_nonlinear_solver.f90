module M_Numeric_Nonlinear_Solver
!-----------------------------------------------------------------------
! Purpose   : General Interface to Nonlinear Newton Type Solvers 
!
! optional METHODS :
! ----------------
!
!  - "NEWTON"       Newton Method                  ( M_Numeric_Newton )
!  - "NEWTONCHESS"  Newton Method like Chess       ( M_Numeric_Newton )
!  - "LINESEARCH"   Newton Method with Linesearch  ( M_Numeric_Newton )
!  - "BROYDEN"      Broyden Method                 ( M_Numeric_Broyden )
!  - "TENSOLVE"     Tensolve Method                ( M_Numeric_Tensolve )  
!
!-----------------------------------------------------------------------
  implicit none
  private
  !
  public:: Nonlinear_Solver
  !
contains

subroutine Nonlinear_Solver( &
& cMethod,  & !input= name of the method
& vX,       & !inout= the initial guess, and the root returned
& vIsPlus,  & !in=    if(vIsPlus(i)) then vX(i) must be >Zero

& Residual, & !in= external residual function
& Jacobian, & !in= external jacobian function
& Converge, & !in= external converge test function

& TolF,     & !in=    convergence criterion on function values 
& TolX,     & !in=    convergence criterion on dx
& bFinDif,  & !in=    use numeric jacobian
& MaxIts,   & !in=    maximum number of iterations
!---
& Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
& Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
& Gradient, & !out=   value of the gradient ???
& Nits,     & !out=   number of iterations
& Check,    & !out=   if Check, should check convergence ???
& iErr)       !out=   error code
!
!-----------------------------------------------------------------------
! Purpose : Nonlinear Solver Wrapper 
!-----------------------------------------------------------------------
!iErr= 0 : convergence reached -> OK
!iErr=-1 : MaxIts reached without convergence 
!iErr=-2 : singular jacobian
!iErr=-3 : roundoff problem in linesearch
!iErr=-4 : no convergence on vFunc, problem with gradient ?
!iErr=-5 : no convergence on vFunc, vX very close to previous vX -> stationary point ??
!-----------------------------------------------------------------------
  !
  use M_Numeric_Newton
  use M_Numeric_Broyden
  use M_Numeric_Tensolve
  use M_Kinds
  use M_Trace
  implicit none
  character(len=*) :: cMethod
  real(dp),dimension(:), intent(inout):: vX
  logical, dimension(:), intent(in)   :: vIsPlus
  real(dp),              intent(in)   :: TolF,TolX !!,TolMin
  logical,               intent(in)   :: bFinDif
  integer,               intent(in)   :: MaxIts
  real(dp),              intent(out)  :: Error_F,Delta_X,Gradient
  logical,               intent(out)  :: Check
  integer,               intent(out)  :: nIts,iErr
  interface
    function Residual(v)
      use M_Kinds
      implicit none
      real(dp),dimension(:),intent(in):: v
      real(dp),dimension(size(v))     :: Residual
    end function Residual
    !
    subroutine Jacobian(v,t)
      use M_Kinds
      implicit none
      real(dp),dimension(:),              intent(in) :: v
      real(dp),dimension(size(v),size(v)),intent(out):: t
    end subroutine Jacobian
    !
    logical function Converge(v,tolv)
      use M_Kinds
      implicit none
      real(dp),intent(in):: v(:),tolv(:)
    end function Converge
  end interface
  !
  select case(trim(cMethod))
  !
  case("NEWTLNSRCH")
    call NewtLnsrch( & 
    & vX,      & !inout= the initial guess, and the root returned
    
    & Residual, Jacobian, Converge, &
    
    & TolF,    & !in=    convergence criterion on function values 
    & TolX,    & !in=    convergence criterion on dx
    & bFinDif, & !in=    use numeric Jacobian
    & MaxIts,  & !in=    maximum number of iterations
    
    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Gradient,& !out=   
    & Nits,    & !out=   number of iterations
    & Check,   & !out=   if Check, should check convergence
    & iErr     & !out=   error code
    & )
    !
  case("NEWTONPRESS")
    call Newton_Press( & 
    & vX,      & !inout= the initial guess, and the root returned
    & vIsPlus, & !
    
    & Residual, Jacobian, Converge, &
    
    & TolF,    & !in=    convergence criterion on function values 
    & TolX,    & !in=    convergence criterion on dx
    & bFinDif, & !in=    use numeric Jacobian
    & MaxIts,  & !in=    maximum number of iterations
    
    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,    & !out=   number of iterations
    & iErr     & !out=   error code
    & )
    !
  case("NEWTONKELLEY")
    call Newton_Kelley( & 
    & vX,      & !inout= the initial guess, and the root returned
    & vIsPlus, & !
    
    & Residual, Jacobian, Converge, &
    
    & TolF,    & !in=    convergence criterion on function values 
    & TolX,    & !in=    convergence criterion on dx
    & bFinDif, & !in=    use numeric Jacobian
    & MaxIts,  & !in=    maximum number of iterations
    
    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,    & !out=   number of iterations
    & iErr)      !out=   error code

  case("NEWTONWALKER")
    call Newton_Walker( & 
    & vX,       & !inout= the initial guess, and the root returned
    & vIsPlus,  & !
    
    & Residual, &
    & Jacobian, &
    & Converge, &
    
    & TolF,     & !in=    convergence criterion on function values 
    & TolX,     & !in=    convergence criterion on dx
    & bFinDif,  & !in=    use numeric Jacobian
    & MaxIts,   & !in=    maximum number of iterations
    
    & Error_F,  & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X,  & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,     & !out=   number of iterations
    & iErr)      !out=   error code

  case("NEWTON")
    call Newton( & 
    & vX,      & !inout= the initial guess, and the root returned
    & Residual, Jacobian, &
    & TolF,    & !in=    convergence criterion on function values 
    & TolX,    & !in=    convergence criterion on dx
    & MaxIts,  & !in=    maximum number of iterations
    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,    & !out=   number of iterations
    & iErr)      !out=   error code

  case("NEWTONCHESS")
    call NewtonChess( & 
    & vX,      & !inout= the initial guess, and the root returned
    & Residual, Jacobian, &
    & TolF,    & !in=    convergence criterion on function values 
    & TolX,    & !in=    convergence criterion on dx
    & MaxIts,  & !in=    maximum number of iterations
    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Nits,    & !out=   number of iterations
    & iErr)      !out=   error code
    !
  case("BROYDEN")
    call Broyden( &
    & vX,      & !inout= the initial guess, and the root returned
    & Residual, Jacobian, &
    & TolF,    & !in=    convergence criterion on function values 
    & TolX,    & !in=    convergence criterion on dx
    & bFinDif, & !in=    use numeric Jacobian
    & MaxIts,  & !in=    maximum number of iterations
    & Error_F, & !out=   MAXVAL(ABS(fVec(:)))
    & Delta_X, & !out=   MAXVAL( ABS(vX(:)-vXOld(:)) / MAX(ABS(vX(:)),One) )
    & Gradient, & !out=   
    & Nits,    & !out=   number of iterations
    & Check,   & !out=   if Check, should check convergence
    & iErr)      !out=   error code
    !
  case("TENSOLVE")
    call TenSolve( &
    & "NEWTON","LINESEARCH", &
    & vX,       & 
    & Residual, &
    & Jacobian, &
    & TolF,     & !in=    convergence criterion on function values 
    & bFinDif,  & !in=    use numeric Jacobian
    & MaxIts,   & !in=    maximum number of iterations
    & Error_F,  & !out
    & Nits,     & !out=   number of iterations
    & iErr)       !out=   error code
    !
  case default
    call Fatal_("Nonlinear_Solver, wrong CMethod = "//trim(cMethod)) 
  end select
  !
end subroutine Nonlinear_Solver
  !
end module M_Numeric_Nonlinear_Solver
