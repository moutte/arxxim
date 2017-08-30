module M_Box_Utils

use M_Kinds
implicit none
private

public :: Limited_Exp
public :: Limited_Log

contains

!--------

  function Limited_Exp(vX)  result (vY)
    !==================================================
    ! Purpose   : truncated exponential to avoid overflow
    !==================================================
    use M_Numeric_Const,      only: MinExpDP,MaxExpDP,Ln10
    implicit none
    real(dp), intent(in)  :: vX(:)
    real(dp) :: vY(size(VX))

    !---------------------------------------------------
    where(vX>MinExpDP .and. vX<MaxExpDP) vY = exp(vX)
    where(vX<MinExpDP) vY =  exp(MinExpDP) 
    where(vX>MaxExpDP) vY =  exp(MaxExpDP) 

  end function Limited_Exp

  !---
  
  function Limited_Log(vX)  result (vLogx)
    use M_Numeric_Const,      only:TinyDP
    implicit none
    real(dp), intent(in)  :: vX(:)
    real(dp)  :: vLogx(size(VX))
    !
    where( vX> TinyDP )
       vLogx= log(vX)
    elsewhere
       where(vX>=0.) vLogx= log(TinyDP)
       where(vX<0.)  vLogx= log(vX)     ! Error Log Undefined !    
    end where

  end function Limited_Log


end module M_Box_Utils
