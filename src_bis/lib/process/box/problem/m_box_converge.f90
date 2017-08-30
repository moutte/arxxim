module M_Box_Converge

  use M_Kinds
  implicit none
  private

  public :: Box_Converge
  
contains

  logical function Box_Converge(v, tolv)
    !==================================================
    ! Purpose   : makes possible adaptation of 
    ! the convergence criteria to special type of 
    ! condition .e.g. material conservation 
    ! vs potential, kinetic, etc.
    !==================================================
    implicit none
    !
    real(dp),intent(in):: v(:), tolv(:)
    Box_Converge= MAXVAL(ABS(v)) < 1.D-6 !for example

  end function Box_Converge

end module M_Box_Converge
