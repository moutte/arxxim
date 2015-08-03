MODULE M_SetOption

  USE M_Kinds
  USE M_Trace
  USE M_IOTools
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SetOption
  PUBLIC :: ComputeTimeFactor
  PUBLIC :: SetOption_Time
  PUBLIC :: SetOption_FlowRate

   !// SetOption Procedure Interface
  INTERFACE SetOption
     MODULE PROCEDURE SetOption_Real
  END INTERFACE

  INTERFACE SetOption
     MODULE PROCEDURE SetOption_Integer
  END INTERFACE

  INTERFACE SetOption
     MODULE PROCEDURE SetOption_String
  END INTERFACE

  INTERFACE SetOption
     MODULE PROCEDURE SetOption_Logical
  END INTERFACE

  !// Private Real Procedure
  PRIVATE :: SetOption_Real
  PRIVATE :: SetOption_Integer
  PRIVATE :: SetOption_String
  PRIVATE :: SetOption_Logical

CONTAINS

  !---

  SUBROUTINE ComputeTimeFactor( TimeFactor, Tunit)
    IMPLICIT NONE
    !---
    REAL(dp), INTENT(out) :: TimeFactor
    CHARACTER(LEN=*), INTENT(in) :: Tunit
    !---
    SELECT CASE(trim(TUNIT))
    CASE ("SECOND") ; TimeFactor = 1.D0
    CASE ("HOUR")   ; TimeFactor = 3600
    CASE ("DAY")    ; TimeFactor = 3600 * 24
    CASE ("YEAR")   ; TimeFactor = 3600 * 24 * 365
    CASE DEFAULT
        CALL FATAL_("Wrong Time Unit "// TUnit )
    END SELECT

  END SUBROUTINE ComputeTimeFactor

  !---

  SUBROUTINE SetOption_Time( Keyword , Str, Var, Tunit)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: Keyword
    CHARACTER(LEN=*), INTENT(in) :: Str
    CHARACTER(LEN=*), INTENT(in) :: Tunit
    !---
    REAL(dp), INTENT(out) :: Var
    !---
    CALL SetOption_Info("SetOption_Time", Keyword, trim(Str)//" "//trim(Tunit) )
    CALL WrdToReal(Str//' ',Var)
    SELECT CASE(trim(TUNIT))
    CASE ("SECOND") ; Var = Var
    CASE ("HOUR")    ; Var = Var * 3600
    CASE ("DAY")    ; Var = Var * 3600 * 24
    CASE ("YEAR")   ; Var = Var * 3600 * 24 * 365
    CASE DEFAULT
        CALL FATAL_("Wrong Time Unit "// TUnit )
    END SELECT

  END SUBROUTINE SetOption_Time

  !---

  SUBROUTINE SetOption_FlowRate( Keyword , Str, Var, Tunit)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: Keyword
    CHARACTER(LEN=*), INTENT(in) :: Str
    CHARACTER(LEN=*), INTENT(in) :: Tunit
    !---
    REAL(dp), INTENT(out) :: Var
    !---
    CALL SetOption_Info("SetOption_FlowRate", Keyword, trim(Str)//" m3 by"//trim(Tunit) )
    CALL WrdToReal(Str//' ',Var)
    SELECT CASE(trim(TUNIT))
    CASE ("SECOND")  ; Var = Var
    CASE ("HOUR")    ; Var = Var / 3600
    CASE ("DAY")     ; Var = Var / ( 3600 * 24 )
    CASE ("YEAR")    ; Var = Var / ( 3600 * 24 * 365 )
    CASE DEFAULT
        CALL FATAL_("Wrong Time Unit "// TUnit )
    END SELECT

  END SUBROUTINE SetOption_FlowRate

  !---

  SUBROUTINE SetOption_Real( Keyword , Str, Var)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: Keyword
    CHARACTER(LEN=*), INTENT(in) :: Str
    !---
    REAL(dp), INTENT(out) :: Var
    !---

    CALL SetOption_Info("SetOption_Real", Keyword,Str)
    CALL WrdToReal(Str//' ',Var)
  END SUBROUTINE SetOption_Real

  !---

  SUBROUTINE SetOption_Integer( Keyword , Str, Var)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: Keyword
    CHARACTER(LEN=*), INTENT(in) :: Str
    !---
    INTEGER, INTENT(out) :: Var
    !---
    CALL SetOption_Info("SetOption_Integer", Keyword,Str)
    CALL WrdToInt(Str//' ',Var)
  END SUBROUTINE SetOption_Integer

  !---

  SUBROUTINE SetOption_String( Keyword , Str, Var)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: Keyword
    CHARACTER(LEN=*), INTENT(in) :: Str
    !---
    CHARACTER(LEN=*), INTENT(inout) :: Var
    !---
    CALL SetOption_Info("SetOption_String", Keyword,Str)
    Var = trim(adjustl(Str))

  END SUBROUTINE SetOption_String

  !---

  SUBROUTINE SetOption_Logical( Keyword , Str, Var)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: Keyword
    CHARACTER(LEN=*), INTENT(in) :: Str
    !---
    LOGICAL, INTENT(out) :: Var
    !---
    CALL SetOption_Info("SetOption_Logical", Keyword,Str)
    CALL WrdToLogical(Str,Var)

  END SUBROUTINE SetOption_Logical

  !------------------------ Utils ------------------------------------------------------

  SUBROUTINE SetOption_Info(Tag, Keyword, Str)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(in) :: Keyword
    CHARACTER(LEN=*), INTENT(in) :: Str
    CHARACTER(LEN=*), INTENT(IN) :: Tag
    !---
    CALL Info_(Tag //" ["//TRIM(Keyword)//"|"//TRIM(Str)//"]");

  END SUBROUTINE SetOption_Info

  !---

  SUBROUTINE WrdToLogical(StrIn,LOut) !conversion of string StrIn to integer IOut
    IMPLICIT NONE
    CHARACTER*(*),INTENT(IN) :: StrIn
    LOGICAL,      INTENT(OUT):: LOut
    !---
    CHARACTER(LEN=10):: StrVal
    !---
    StrVal = StrIn
    CALL Str_Upper(StrVal)

    LOut = .FALSE.
    SELECT CASE(StrVal)
    CASE("TRUE", "T", "YES")
       LOut = .TRUE.
    CASE("FALSE", "F", "NO")
       LOut = .FALSE.
    CASE DEFAULT
       CALL Fatal_("WrdToLogical, error in translation of "//trim(StrIn))
    END SELECT

  ENDSUBROUTINE WrdToLogical

END MODULE M_SetOption
