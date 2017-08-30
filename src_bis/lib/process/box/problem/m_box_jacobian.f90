module M_Box_Jacobian
  
  use M_Kinds
  
  implicit none
  private
  
  !-- Public Functions
  public:: Box_Jacobian
  
contains
  
  subroutine Box_Jacobian(vX,tJac) 
    !==================================================
    ! Purpose   : compute the jacobian for the box 
    ! problem used in reaction transport coupled mode
    !==================================================
    use M_Kinds
    use M_Box_Jacobian_Base,  only : Box_Jacobian_Base
    use M_Box_Jacobian_Dynam, only : Box_Jacobian_Dynam
    implicit none

    real(dp),dimension(:),                 intent(in) :: vX
    real(dp),dimension(size(vX),size(vX)), intent(out):: tJac
    
    !call Box_Jacobian_Base(vX,tJac)
    call Box_Jacobian_Dynam(vX,tJac)    
    
  end subroutine Box_Jacobian

end module M_Box_Jacobian
