MODULE M_Formula_Utils
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: Name_Index

  CONTAINS
  
  !---
  
  FUNCTION Name_Index(Str,V) RESULT(index)
    !=====================================================
    ! Position of Name with Name==Str in VName
    !=====================================================
    USE m_trace, ONLY : iDebug
    IMPLICIT NONE
    CHARACTER(LEN=*),              INTENT(IN)::Str
    CHARACTER(LEN=*), DIMENSION(:),INTENT(IN)::V
    !---
    INTEGER :: I
    INTEGER :: index
    !---
    index=0
    IF(SIZE(V)==0) RETURN
    I=0
    DO
      I=I+1 
      IF(iDebug==5) WRITE(*,*) "Testing" , TRIM(Str), TRIM(V(I))
      IF(TRIM(Str)==TRIM(V(I))) THEN
        index=I
        EXIT
      ENDIF
      IF(I==SIZE(V)) EXIT
    ENDDO !IF Str not found -> I=0
    RETURN
  END FUNCTION Name_Index
 
END MODULE M_Formula_Utils
