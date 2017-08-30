module M_Box_Residual

  use M_Kinds 
 
  implicit none
  private

  !-- Public Functions
  public:: Box_Residual

contains

  !-------

  function Box_Residual(vX)
    !==================================================
    ! Purpose : compute the residual for the box 
    ! problem used in reaction transport coupled mode
    !==================================================
    use M_Box_Residual_Base
    use M_Box_Residual_Dynam
    !--
    
    implicit none

    real(dp),dimension(:),intent(in) :: vX
    real(dp),dimension(size(vX))     :: Box_Residual
    
    !Box_Residual = Box_Residual_Base(vX)
    Box_Residual = Box_Residual_Dynam(vX)

  end function Box_Residual
  
end module M_Box_Residual
