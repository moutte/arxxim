module M_Read_Value

  use M_Trace
  implicit none
  private
  
  interface Read_Value
    module procedure Read_Value_Real
    module procedure Read_Value_Integer
    module procedure Read_Value_Logical
    module procedure Read_Value_String
  end interface

  public :: Read_Value

  private :: Read_Value_Real
  private :: Read_Value_Integer
  private :: Read_Value_Logical
  private :: Read_Value_String

contains

    !---
  subroutine Read_Value_Real(Str, X) 
    implicit none
    character(len=*),intent(in)::Str
    real(kind=8),intent(out)::X
    !--
    call Message_("Input ",Str//" , value = ? ")
    read(*,*) X

  end subroutine Read_Value_Real
  
  !---
  
  subroutine Read_Value_Integer(Str, X) 
    implicit none
    character(len=*),intent(in)::Str
    integer,intent(out)::X
    !--
    call Message_("Input ",Str//" , value = ? ")
    read(*,*) X

  end subroutine Read_Value_Integer

  !---

  subroutine Read_Value_Logical(Str, X) 
    implicit none
    character(len=*),intent(in)::Str
    logical,intent(out)::X
    !--
    
    call Message_("Input ",Str//" , value = ? ")
    read(*,*) X   
 
  end subroutine Read_Value_Logical

  !---

  subroutine Read_Value_String(Str, X) 
    implicit none
    character(len=*),intent(in)::Str
    character(len=*),intent(out)::X
    !--
    call Message_("Input ",Str//" , value = ? ")
    read(*,*) X

  end subroutine Read_Value_String

end module M_Read_Value
