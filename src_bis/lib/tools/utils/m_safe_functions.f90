module M_Safe_Functions

  !=========================================
  ! Safe Functions Exp, Log, Log10
  ! Truncated to avoid Out of Range Errors
  !=========================================
  use M_Kinds

  implicit none
  private

  real(dp), parameter:: Ln10= 2.30258509299404D0
  ! real(dp), parameter:: Zero = 0.D0         ! is in M_Kinds
  ! real(dp), parameter:: One = 1.D0          ! is in M_Kinds

  real(dp) :: MinExpDP= -300.D0
  real(dp) :: MaxExpDP= +300.D0
  real(dp) :: TinyDP= 1.D-300

  !// Set Parameters
  public :: Set_MinExpDP
  public :: Set_MaxExpDP
  public :: Set_TinyDP

  !// Scalar values
  public :: Safe_Exp
  public :: Safe_Log
  public :: Safe_Log10

  !// Function Scalar values
  public :: FSafe_Exp
  public :: FSafe_Log
  public :: FSafe_Log10
  
  !// Vector values
  public :: Safe_vExp
  public :: Safe_vLog
  public :: Safe_vLog10

  !// Function Vector values
  public :: FSafe_vExp
  public :: FSafe_vLog
  public :: FSafe_vLog10

contains

  !==================== parameter SETTING =====================

  subroutine Set_MinExpDP(value)
    implicit none
    real(dp) :: value
    MinExpDP= value
  end subroutine Set_MinExpDP

  !---

  subroutine Set_MaxExpDP(value)
    implicit none
    real(dp) :: value
    MaxExpDP= value
  end subroutine Set_MaxExpDP

  !---
  
  subroutine Set_TinyDP(value)
    implicit none
    real(dp) :: value
    TinyDP= value
  end subroutine Set_TinyDP

  !====================>>  SCALAR valueS <<=====================

  subroutine Safe_Exp(X,Expx) 
    implicit none
    real(dp), intent(in)  :: X
    real(dp), intent(out) :: Expx
    !--
    Expx = 0.; 
    if( (X>MinExpDP) .and. (X<MaxExpDP) ) then
       Expx= exp(X)
    else
       if(X<=MinExpDP) Expx= exp(MinExpDP)
       if(X>=MaxExpDP) Expx= exp(MaxExpDP)
    end if

  end subroutine Safe_Exp

  !---
  
  function FSafe_Exp(X) result (Expx) 
    implicit none
    real(dp), intent(in)  :: X
    real(dp) :: Expx
    !--
    Expx = 0.; 
    if( (X>MinExpDP) .and. (X<MaxExpDP) ) then
       Expx= exp(X)
    else
       if(X<=MinExpDP) Expx= exp(MinExpDP)
       if(X>=MaxExpDP) Expx= exp(MaxExpDP)
    end if
  end function FSafe_Exp

  !---

  subroutine Safe_Log(X,Logx) 
    implicit none
    real(dp), intent(in)  :: X
    real(dp), intent(out) :: Logx
    !--
    if( X>TinyDP ) then
       Logx= log(X)
    else
       if ( X>= Zero ) Logx= log(TinyDP)
       if ( X< Zero )  Logx= log(X) ! Error Log Undefined !    
    end if

  end subroutine Safe_Log

  !---
  
  function FSafe_Log(X) result(LogX) 
    implicit none
    real(dp), intent(in)  :: X
    real(dp) :: Logx
    !--
    if( X>TinyDP ) then
      Logx= log(X)
    else
      if ( X>= Zero ) Logx= log(TinyDP)
      if ( X< Zero )  Logx= log(X) ! Error Log Undefined !    
    end if
    
  end function FSafe_Log

  !---

  subroutine Safe_Log10(X,Logx) 
    implicit none
    real(dp), intent(in)  :: X
    real(dp), intent(out) :: Logx
    !--
    call Safe_Log(X,Logx) 
    Logx = LogX + Ln10;

  end subroutine Safe_Log10

  !---
  
   function FSafe_Log10(X) result(Logx) 
    implicit none
    real(dp), intent(in)  :: X
    real(dp) :: Logx
    !--
    call Safe_Log(X,Logx) 
    Logx = LogX + Ln10;

  end function FSafe_Log10

  !====================>>  VECTOR valueS <<=====================


  subroutine Safe_vExp(vX,vExpx) 
    implicit none
    real(dp), intent(in)  :: vX(:)
    real(dp), intent(out) :: vExpx(:)
    !--
    vExpx = 0.
    where((vX>MinExpDP) .and. (vX<MaxExpDP))
      vExpx=exp(vX)
    elsewhere
      where(vX<=MinExpDP) vExpx= exp(MinExpDP)
      where(vX>=MaxExpDP) vExpx= exp(MaxExpDP)
    end where

  end subroutine Safe_vExp

  !---

  function FSafe_vExp(vX) result(vExpx) 
    implicit none
    real(dp), intent(in)  :: vX(:)
    real(dp) :: vExpx(size(vX))
    !--
    vExpx = 0.
    where((vX>MinExpDP) .and. (vX<MaxExpDP))
      vExpx=exp(vX)
    elsewhere
      where(vX<=MinExpDP) vExpx= exp(MinExpDP)
      where(vX>=MaxExpDP) vExpx= exp(MaxExpDP)
    end where
  end function FSafe_vExp
    
  !---
  
  subroutine Safe_vLog(vX,vLogx) 
    implicit none
    real(dp), intent(in)  :: vX(:)
    real(dp), intent(out) :: vLogx(:)
    !
    where( vX> TinyDP )
       vLogx= log(vX)
    elsewhere
       where(vX>=0.) vLogx= log(TinyDP)
       where(vX<0.)  vLogx= log(vX)     ! Error Log Undefined !    
    end where

  end subroutine Safe_vLog

  !---

  function FSafe_vLog(vX) result (vLogx)
    implicit none
    real(dp), intent(in)  :: vX(:)
    real(dp) :: vLogx(size(vX))
    !
    where( vX> TinyDP )
       vLogx= log(vX)
    elsewhere
       where(vX>=0.) vLogx= log(TinyDP)
       where(vX<0.)  vLogx= log(vX)     ! Error Log Undefined !    
    end where

  end function FSafe_vLog

  !---

  subroutine Safe_vLog10(vX,vLogx) 
    implicit none
    real(dp), intent(in)  :: vX(:)
    real(dp), intent(out) :: vLogx(:)
    !--
    call Safe_vLog(vX,vLogx) 
    vLogx = vLogx + Ln10

  end subroutine Safe_vLog10

  !---

  function FSafe_vLog10(vX) result (vLogx) 
    implicit none
    real(dp), intent(in)  :: vX(:)
    real(dp)  :: vLogx(size(vX))
    !--
    call Safe_vLog(vX,vLogx) 
    vLogx = vLogx + Ln10

  end function FSafe_vLog10

  !---

end module M_Safe_Functions
