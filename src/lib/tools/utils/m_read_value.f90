MODULE M_Read_Value

  USE M_Trace
  IMPLICIT NONE
  PRIVATE
  
  INTERFACE Read_Value
    MODULE PROCEDURE Read_Value_Real
    MODULE PROCEDURE Read_Value_Integer
    MODULE PROCEDURE Read_Value_Logical
    MODULE PROCEDURE Read_Value_String
  END INTERFACE

  PUBLIC :: Read_Value

  PRIVATE :: Read_Value_Real
  PRIVATE :: Read_Value_Integer
  PRIVATE :: Read_Value_Logical
  PRIVATE :: Read_Value_String

CONTAINS

    !---
  SUBROUTINE Read_Value_Real(Str, X) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    REAL(kind=8),INTENT(OUT)::X
    !--
    CALL Message_("Input ",Str//" , VALUE = ? ")
    READ(*,*) X

  ENDSUBROUTINE Read_Value_Real
  
  !---
  
  SUBROUTINE Read_Value_Integer(Str, X) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    INTEGER,INTENT(OUT)::X
    !--
    CALL Message_("Input ",Str//" , VALUE = ? ")
    READ(*,*) X

  ENDSUBROUTINE Read_Value_Integer

  !---

  SUBROUTINE Read_Value_Logical(Str, X) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    LOGICAL,INTENT(OUT)::X
    !--
    
    CALL Message_("Input ",Str//" , VALUE = ? ")
    READ(*,*) X   
 
  ENDSUBROUTINE Read_Value_Logical

  !---

  SUBROUTINE Read_Value_String(Str, X) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    CHARACTER(LEN=*),INTENT(OUT)::X
    !--
    CALL Message_("Input ",Str//" , VALUE = ? ")
    READ(*,*) X

  ENDSUBROUTINE Read_Value_String

END MODULE M_Read_Value
