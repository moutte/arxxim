  MODULE M_Info_Value

  USE M_Trace
  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=8), PRIVATE :: RFormat = 'G16.6'
  CHARACTER(LEN=8), PRIVATE :: IFormat = 'I8'

  INTERFACE Info_Value
    MODULE PROCEDURE Info_Value_Real
    MODULE PROCEDURE Info_Value_Integer
    MODULE PROCEDURE Info_Value_Logical
    MODULE PROCEDURE Info_Value_String
  END INTERFACE

  PUBLIC :: Info_Value

  PRIVATE :: Info_Value_Real
  PRIVATE :: Info_Value_Integer
  PRIVATE :: Info_Value_Logical
  PRIVATE :: Info_Value_String

CONTAINS

    !---
  SUBROUTINE Info_Value_Real(Str, X) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    REAL(kind=8),INTENT(IN)::X
    !--
    CHARACTER(LEN=20)::W
    !--
    WRITE(W, '('//RFormat//')' ) X
    CALL Message_("Info",Str//" , VALUE =  "//W )

  ENDSUBROUTINE Info_Value_Real
  
  !---
  
  SUBROUTINE Info_Value_Integer(Str, X) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    INTEGER,INTENT(IN)::X
    !--
    CHARACTER(LEN=20)::W
    !--
    WRITE(W, '('//IFormat//')' ) X
    CALL Message_("Info",Str//" , VALUE = "//W )

  ENDSUBROUTINE Info_Value_Integer

  !---

  SUBROUTINE Info_Value_Logical(Str, X) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    LOGICAL,INTENT(IN)::X
    !--
    CHARACTER(LEN=20)::W
    !--
    W = ".FALSE."
    IF ( X ) W = ".TRUE."
    
    CALL Message_("Info",Str//" , VALUE = "//W )

  ENDSUBROUTINE Info_Value_Logical

  !---

  SUBROUTINE Info_Value_String(Str, X) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    CHARACTER(LEN=*),INTENT(IN)::X
    !--
    CALL Message_("Info",Str//" , VALUE = "//trim(X) )

  ENDSUBROUTINE Info_Value_String

END MODULE M_Info_Value
