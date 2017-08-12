  module M_Info_Value

  use M_Trace
  implicit none
  private

  character(len=8), private :: RFormat = 'G16.6'
  character(len=8), private :: iformat = 'I8'

  interface Info_Value
    module procedure Info_Value_Real
    module procedure Info_Value_Integer
    module procedure Info_Value_Logical
    module procedure Info_Value_String
  end interface

  public :: Info_Value

  private :: Info_Value_Real
  private :: Info_Value_Integer
  private :: Info_Value_Logical
  private :: Info_Value_String

contains

    !---
  subroutine Info_Value_Real(Str, X) 
    use M_Kinds
    implicit none
    character(len=*),intent(in)::Str
    real(dp),intent(in)::X
    !--
    character(len=20)::W
    !--
    write(W, '('//RFormat//')' ) X
    call Message_("Info",Str//" , value =  "//W )

  end subroutine Info_Value_Real
  
  !---
  
  subroutine Info_Value_Integer(Str, X) 
    implicit none
    character(len=*),intent(in)::Str
    integer,intent(in)::X
    !--
    character(len=20)::W
    !--
    write(W, '('//iformat//')' ) X
    call Message_("Info",Str//" , value = "//W )

  end subroutine Info_Value_Integer

  !---

  subroutine Info_Value_Logical(Str, X) 
    implicit none
    character(len=*),intent(in)::Str
    logical,intent(in)::X
    !--
    character(len=20)::W
    !--
    W = ".false."
    if ( X ) W = ".true."
    
    call Message_("Info",Str//" , value = "//W )

  end subroutine Info_Value_Logical

  !---

  subroutine Info_Value_String(Str, X) 
    implicit none
    character(len=*),intent(in)::Str
    character(len=*),intent(in)::X
    !--
    call Message_("Info",Str//" , value = "//trim(X) )

  end subroutine Info_Value_String

end module M_Info_Value
