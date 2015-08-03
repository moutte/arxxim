MODULE M_Dtb_Read_Tools
!.tools for database reading routines
  USE M_Kinds
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: DtbRead_Build_vElement
  PUBLIC:: DtbRead_Build_ExtendedFormula
  !
CONTAINS

SUBROUTINE DtbRead_Build_vElement(vEle,vElement)
  USE M_IOTools,ONLY: cUpper,cLower
  USE M_T_Element,ONLY: T_Element
  !
  TYPE(T_Element), INTENT(IN) :: vEle(:)
  CHARACTER(LEN=2),INTENT(OUT):: vElement(:)
  !
  CHARACTER(LEN=2):: ss
  INTEGER:: i
  !
  vElement= "  "
  DO i= 1, SIZE(vEle)
    ss="  "
    ss(1:1)= cUpper(vEle(i)%NamEl(1:1))
    IF(vEle(i)%NamEl(2:2)/="_") ss(2:2)= cLower(vEle(i)%NamEl(2:2))
    vElement(i)= TRIM(ss)
  ENDDO
  !
  vElement(SIZE(vEle)+1)= '+'
  vElement(SIZE(vEle)+2)= '-'
  vElement(SIZE(vEle)+3)= '/'
  !
ENDSUBROUTINE DtbRead_Build_vElement

SUBROUTINE DtbRead_Build_ExtendedFormula(ff,vElement,W,isOk)
!--
!-- transform compact formula to extended formula
!--
  USE M_Trace,  ONLY: t_
  USE M_Formula_Parser
  !
  INTEGER,INTENT(IN):: ff
  CHARACTER(LEN=2),INTENT(IN):: vElement(:)
  CHARACTER(LEN=*),INTENT(INOUT):: W
  LOGICAL,INTENT(OUT):: isOk 
  !
  !! CHARACTER(LEN=10):: FilCode
  !
  REAL(dp) :: b(SIZE(vElement))
  REAL(dp) :: divREAL
  CHARACTER(LEN=71) :: s
  !
  IF(ff>0) WRITE(ff,'(A72,1X)',ADVANCE="NO") TRIM(W)
  !
  CALL Formula_Vector_Read ( W, vElement, s, b, isok, divREAL )
  !
  IF(.not. isok) THEN
    IF(ff>0) WRITE(ff,'(A)') " -> !!!CHECK!!!"
    RETURN
  ENDIF
  !
  !IF(ff>0 .and. nint(divREAL)/=1) &
  !& WRITE(ff,'(A,I3,A)',ADVANCE="NO") "divREAL=",nint(divREAL)," = "
  !
  CALL Formula_Arxim_Build ( vElement, nint(b*divREAL), W )
  !
  IF(nint(divREAL)/=1) THEN
    s=''
    WRITE(s,'(I6)') nint(divREAL)
    s= ADJUSTL(s)
    W= TRIM(W)//'/('//TRIM(s)//')'
  ENDIF
  !
  IF(ff>0) WRITE(ff,'(A)') TRIM(W)
ENDSUBROUTINE DtbRead_Build_ExtendedFormula

ENDMODULE M_Dtb_Read_Tools
