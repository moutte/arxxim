module M_FArgC
  !=============================================================== 
  ! Fortran Command Arguments Functions
  ! Wrapper to System Functions : IARGC and GETARG
  !
  ! These Functions are not compatible with standard F90 
  ! Please Modify this file if you encounter difficulties
  ! concerning these functions with your compiler
  !
  !===============================================================
  implicit none
  private
  !---
  public :: F_IARGC
  public :: F_GETARG

contains

  function F_IARGC()
  
    integer :: F_IARGC
    integer IARGC
    
    F_IARGC = 0
    F_IARGC = IARGC() ! Comment this line if problems

  end function F_IARGC

  !---
  
  subroutine F_GETARG(I, STR)
  
    integer, intent(in) :: I
    character(len=*), intent(out):: STR
    
    STR = ""
    call GETARG(I, STR)       ! Comment this line if problems
    
  end subroutine F_GETARG


end module M_FArgC
