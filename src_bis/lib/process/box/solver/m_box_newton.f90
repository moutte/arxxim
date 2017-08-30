module M_Box_Newton

  !=======================================================
  ! Purpose : Newton Solver Interface
  ! Provide Interface to NonLinear Solver Use
  ! Can be used to Add Numeric Specific Parameters
  !=======================================================

  use M_Kinds
  use M_Trace

  use M_Box_Newton_Vars
  use M_Box_Debug_Vars

  !---
  implicit none
  private

  !// public functionS
  public :: Box_Newton_Init
  public :: Box_Newton_Solve

  private :: Newton_Residual
  private :: Newton_Jacobian
  private :: Newton_Converge

  integer :: nResidual, nJacobian

contains

  !---

  function Newton_Residual(vX)
    !==================================================
    ! Purpose   : Call the model residual 
    !==================================================
    use M_Box_Residual
    implicit none
    !
    real(dp),dimension(:),intent(in) :: vX
    real(dp),dimension(size(vX))     :: Newton_Residual
    !--
    nResidual = nResidual + 1
    if (LDebug_Newton) call Warning_("Residual - NEval Residual = "//trim(Int2Str(nResidual)))
    Newton_Residual =  Box_Residual(vX)
    
  end function Newton_Residual

  !---

  subroutine Newton_Jacobian(vX,tJac) 
    !==================================================
    ! Purpose   : Call the model jacobian
    !==================================================
    use M_Box_Jacobian
    implicit none
    !
    real(dp),dimension(:),                 intent(in) :: vX
    real(dp),dimension(size(vX),size(vX)), intent(out):: tJac
    !--
    nJacobian = nJacobian + 1
    if (LDebug_Newton)  call Warning_("Jacobian - NEval Jacobian = "//trim(Int2Str(nJacobian)))
    call Box_Jacobian(vX,tJac) 
   

  end subroutine Newton_Jacobian

  !---

  logical function Newton_Converge(v, tolv)
    !==================================================
    ! Purpose   : makes possible adaptation of 
    ! the convergence criteria 
    !==================================================
    implicit none
    real(dp),intent(in):: v(:), tolv(:)
    !--
    call Info_("Test Convergence")
    call Stop_("Call Converge Function")
    Newton_Converge= MAXVAL(ABS(v)) < 1.E-6_dp

  end function Newton_Converge

  !--

  subroutine Box_Newton_InitParam

    call Nonlinear_Solver_Zero(Newton_Param, Newton_Result) 

    Newton_Param%cMethod = cMethod;
    Newton_Param%TolF    = NewtTolF; 
    Newton_Param%TolX    = NewtTolX; 
    Newton_Param%bFinDif = bFinDif
    Newton_Param%MaxIts  = NewtIterMax;
    
    nResidual = 0
    nJacobian = 0
    
  end subroutine Box_Newton_InitParam

  !---

  subroutine Box_Newton_Init

    implicit none   
    call Box_Newton_InitParam

  end subroutine Box_Newton_Init

  !---

  subroutine Box_Newton_Solve(vX_Sol,  OKNewton)

    use M_Numeric_Nonlinear_Solver_Typed

    implicit none  
    !---
    real(dp), intent(inout):: vX_Sol(:)
    logical, intent(out):: OKNewton
    !
    logical:: vIsPlus(size(vX_Sol))
    !---
    call Info_("Box_Newton_Solve")

    if (LDebug_Newton) call Nonlinear_Solver_Param_Info(Newton_Param)
    
    vIsPlus(:)= .false.
    
    call Nonlinear_Solver_Typed( &
    & vX_Sol,       & !inout= the initial guess, and the root returned
    & vIsPlus,      & !in=    if(vIsPlus(i)) then vX(i) must be >Zero
    & Newton_Residual, & !in= external residual function
    & Newton_Jacobian, & !in= external jacobian function
    & Newton_Converge, & !in= external converge test function
    & Newton_Param, & !in= input parameters
    & Newton_Result & !out= output results
    & )

    OKNewton = ( Newton_Result% iErr.eq.0 )
    NewtonIter = Newton_Result% Nits
    Total_NewtonIter = Total_NewtonIter + NewtonIter
    
    if (LDebug_Newton) call Nonlinear_Solver_Result_Info(Newton_Result)
    call Box_Newton_Solve_Info(OKNewton)
    nResidual = 0
    nJacobian = 0

  end subroutine Box_Newton_Solve

  !---

  subroutine Box_Newton_Solve_Info(OKNewton)

    use M_Numeric_Show
    use M_Info_Value
    implicit none
    logical, intent(in):: OKNewton
    
    if (LDebug_Newton) then
       write(*,*) '========================================================'
       if (OKNewton)        write(*,*) "OK ! NEWTON CONVERGED"
       if (.not.(OKNewton)) write(*,*) "ERRORS ! NEWTON NOT CONVERGED"
       write(*,*) '========================================================'
       call Info_Value("Number of Iterations ...  ", Newton_Result%Nits)
       call Info_Value("ErrF                 ...  ", Newton_Result%ErrF)
       call Info_Value("iErr                 ...  ", Newton_Result%iErr)
       write(*,*) '========================================================'
       call NewtonError_Sho(Newton_Result%iErr)
       write(*,*) '========================================================'
    end if
    
  end subroutine Box_Newton_Solve_Info

  !---
  
  function Int2Str(I) result (S)
    implicit none
    character(len=6) :: S
    character(len=6) :: Str
    integer, intent(in) :: I
    !---
    write(Str,'(I6)') I
    S = adjustl(Str)
    
  end function INT2STR
  
end module M_Box_Newton
