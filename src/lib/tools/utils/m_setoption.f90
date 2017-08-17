module M_SetOption

  use M_Kinds
  use M_Trace
  use M_IOTools
  implicit none
  private

  public :: SetOption
  public :: ComputeTimeFactor
  public :: SetOption_Time
  public :: SetOption_FlowRate

   !// SetOption Procedure Interface
  interface SetOption
     module procedure SetOption_Real
  end interface

  interface SetOption
     module procedure SetOption_Integer
  end interface

  interface SetOption
     module procedure SetOption_String
  end interface

  interface SetOption
     module procedure SetOption_Logical
  end interface

  !// Private Real Procedure
  private :: SetOption_Real
  private :: SetOption_Integer
  private :: SetOption_String
  private :: SetOption_Logical

contains

  !---

  subroutine ComputeTimeFactor( TimeFactor, Tunit)
    implicit none
    !---
    real(dp), intent(out) :: TimeFactor
    character(len=*), intent(in) :: Tunit
    !---
    select case(trim(TUNIT))
    case ("SECOND") ; TimeFactor = 1.D0
    case ("HOUR")   ; TimeFactor = 3600
    case ("DAY")    ; TimeFactor = 3600 * 24
    case ("YEAR")   ; TimeFactor = 3600 * 24 * 365
    case default
        call FATAL_("Wrong Time Unit "// TUnit )
    end select

  end subroutine ComputeTimeFactor

  !---

  subroutine SetOption_Time( Keyword , Str, Var, Tunit)
    implicit none
    character(len=*), intent(in) :: Keyword
    character(len=*), intent(in) :: Str
    character(len=*), intent(in) :: Tunit
    !---
    real(dp), intent(out) :: Var
    !---
    call SetOption_Info("SetOption_Time", Keyword, trim(Str)//" "//trim(Tunit) )
    call WrdToReal(Str//' ',Var)
    select case(trim(TUNIT))
    case ("SECOND") ; Var = Var
    case ("HOUR")    ; Var = Var * 3600
    case ("DAY")    ; Var = Var * 3600 * 24
    case ("YEAR")   ; Var = Var * 3600 * 24 * 365
    case default
        call FATAL_("Wrong Time Unit "// TUnit )
    end select

  end subroutine SetOption_Time

  !---

  subroutine SetOption_FlowRate( Keyword , Str, Var, Tunit)
    implicit none
    character(len=*), intent(in) :: Keyword
    character(len=*), intent(in) :: Str
    character(len=*), intent(in) :: Tunit
    !---
    real(dp), intent(out) :: Var
    !---
    call SetOption_Info("SetOption_FlowRate", Keyword, trim(Str)//" m3 by"//trim(Tunit) )
    call WrdToReal(Str//' ',Var)
    select case(trim(TUNIT))
    case ("SECOND")  ; Var = Var
    case ("HOUR")    ; Var = Var / 3600
    case ("DAY")     ; Var = Var / ( 3600 * 24 )
    case ("YEAR")    ; Var = Var / ( 3600 * 24 * 365 )
    case default
        call FATAL_("Wrong Time Unit "// TUnit )
    end select

  end subroutine SetOption_FlowRate

  !---

  subroutine SetOption_Real( Keyword , Str, Var)
    implicit none
    character(len=*), intent(in) :: Keyword
    character(len=*), intent(in) :: Str
    !---
    real(dp), intent(out) :: Var
    !---

    call SetOption_Info("SetOption_Real", Keyword,Str)
    call WrdToReal(Str//' ',Var)
  end subroutine SetOption_Real

  !---

  subroutine SetOption_Integer( Keyword , Str, Var)
    implicit none
    character(len=*), intent(in) :: Keyword
    character(len=*), intent(in) :: Str
    !---
    integer, intent(out) :: Var
    !---
    call SetOption_Info("SetOption_Integer", Keyword,Str)
    call WrdToInt(Str//' ',Var)
  end subroutine SetOption_Integer

  !---

  subroutine SetOption_String( Keyword , Str, Var)
    implicit none
    character(len=*), intent(in) :: Keyword
    character(len=*), intent(in) :: Str
    !---
    character(len=*), intent(inout) :: Var
    !---
    call SetOption_Info("SetOption_String", Keyword,Str)
    Var = trim(adjustl(Str))

  end subroutine SetOption_String

  !---

  subroutine SetOption_Logical( Keyword , Str, Var)
    implicit none
    character(len=*), intent(in) :: Keyword
    character(len=*), intent(in) :: Str
    !---
    logical, intent(out) :: Var
    !---
    call SetOption_Info("SetOption_Logical", Keyword,Str)
    call WrdToLogical(Str,Var)

  end subroutine SetOption_Logical

  !------------------------ Utils ------------------------------------------------------

  subroutine SetOption_Info(Tag, Keyword, Str)
    implicit none
    character(len=*), intent(in) :: Keyword
    character(len=*), intent(in) :: Str
    character(len=*), intent(in) :: Tag
    !---
    call Info_(Tag //" ["//trim(Keyword)//"|"//trim(Str)//"]");

  end subroutine SetOption_Info

  !---

  subroutine WrdToLogical(StrIn,LOut) !conversion of string StrIn to integer IOut
    implicit none
    character*(*),intent(in) :: StrIn
    logical,      intent(out):: LOut
    !---
    character(len=10):: StrVal
    !---
    StrVal = StrIn
    call Str_Upper(StrVal)

    LOut = .false.
    select case(StrVal)
    case("TRUE", "T", "YES")
       LOut = .true.
    case("FALSE", "F", "NO")
       LOut = .false.
    case default
       call Fatal_("WrdToLogical, error in translation of "//trim(StrIn))
    end select

  end subroutine WrdToLogical

end module M_SetOption
