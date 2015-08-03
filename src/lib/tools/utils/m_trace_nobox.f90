MODULE M_Trace
  !===========================================================
  ! Trace File and Debug Status Manager
  !-----------------------------------------------------------
  ! > Trace Log File 
  ! > Debug Level
  ! > Trace Special Chars ( TAB_, SLASH_, BACKSLASH_ )
  ! > Special Trace Functions : 
  !  - Pause_, Debug_, Info_, Warning_ 
  !  - Fatal_, Stop_
  !=========================================================== 
  IMPLICIT NONE
  !
  PRIVATE
  !
  !// PARAMETERS
  CHARACTER(LEN=20), PARAMETER :: Default_TraceFileName =  "debug_all.log"

  !// PUBLIC DATA
  CHARACTER,PARAMETER,PUBLIC:: T_         = ACHAR(9)
  CHARACTER,PARAMETER,PUBLIC:: TAB_       = ACHAR(9)
  CHARACTER,PARAMETER,PUBLIC:: SLASH_     = ACHAR(47)
  CHARACTER,PARAMETER,PUBLIC:: BACKSLASH_ = ACHAR(92)

  INTEGER,PUBLIC:: iDebug= -1
  INTEGER,PUBLIC:: fTrc= 0
  INTEGER,PUBLIC:: fHtm= 0
  LOGICAL,PUBLIC:: DebugCoores=.FALSE.

  LOGICAL,PUBLIC:: LInfo= .FALSE.
  LOGICAL,PUBLIC:: LWarning= .FALSE.

  !// PUBLIC FUNCTIONS
  PUBLIC:: Trace_Init
  PUBLIC:: Trace_Close
  !
  PUBLIC:: Debug_
  PUBLIC:: Pause_
  PUBLIC:: Info_
  PUBLIC:: Fatal_
  PUBLIC:: Stop_
  PUBLIC:: Warning_
  PUBLIC:: Message_

  !// PRIVATE FUNCTIONS
  PRIVATE:: GetUnit_Trace
  !
  !// PRIVATE DATA
  CHARACTER(LEN=80):: TraceFileName = Default_TraceFileName
  !INTEGER:: fError= 0
  
CONTAINS

  !---

  SUBROUTINE Trace_Reset()
    IMPLICIT NONE
    CALL Trace_Close
    CALL Trace_Init(TraceFileName)
    
  ENDSUBROUTINE Trace_Reset

  !---

  SUBROUTINE Trace_Init(Str)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: Str
    !--
    CHARACTER(LEN=80) :: FileName
    !--
    FileName= Default_TraceFileName
    IF(PRESENT(Str)) FileName= Str 
    
    IF(fTrc==0) THEN
       !// open a new file
       CALL GetUnit_Trace(fTrc)
       TraceFileName= FileName
       OPEN(fTrc,FILE=TRIM(FileName))
    ELSE
       !// nothing to do
       CALL Warning_("TraceFile Already Initialized")
    END IF
   
  ENDSUBROUTINE Trace_Init

  !---

  SUBROUTINE Trace_Close
    IMPLICIT NONE
    IF(fTrc>0) THEN; CLOSE(fTrc)
       fTrc= 0
    ENDIF
  END SUBROUTINE Trace_Close

  !---

  SUBROUTINE Pause_(Str)
    IMPLICIT NONE
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: Str
    !--
    IF(PRESENT(Str)) THEN 
      CALL Message_("Pause",Str)
    END IF
    WRITE(*,'(A)') "... type return to continue." 
    READ(*,*)

  ENDSUBROUTINE Pause_

  !---

  SUBROUTINE Warning_(Str)
    IMPLICIT NONE
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: Str
    !--
    IF(PRESENT(Str)) THEN 
      IF (LWarning) CALL Message_("Warning",Str)
    ENDIF

  ENDSUBROUTINE Warning_

  !---

  SUBROUTINE Debug_(Str) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    !--
    IF (iDebug>0) CALL Message_("Debug",Str)

  ENDSUBROUTINE Debug_

  !---

  SUBROUTINE Info_(Str) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    !--
    IF (LInfo) CALL Message_("Info",Str)

  ENDSUBROUTINE Info_

  !---

  SUBROUTINE Fatal_(Str)
    !~ USE M_BOX_COUPLER_VARS
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    !--
    !~ IF (LCouplerActive) THEN
      !~ CALL Message_("Warning",Str)
      !~ IerrorChemistry = 1
    !~ ELSE
      CALL Message_("Fatal Error",Str)
      STOP "FATAL ERROR"
    !~ ENDIF

  ENDSUBROUTINE Fatal_

  !---

  SUBROUTINE Stop_(Str) 
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Str
    !--
    CALL Message_("Program Stop",Str)
    
    !IF(fError==0) THEN
    !  CALL GetUnit_Trace(fError)
    !  OPEN(fError,FILE="error.log")
    !ENDIF
    !WRITE(fError,'(A)') TRIM(Str)
    
    STOP "PROGRAM STOP" 
    
  ENDSUBROUTINE Stop_

  !---

  SUBROUTINE Message_(Header, Str)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::Header
    CHARACTER(LEN=*),INTENT(IN), optional :: Str

    !// Trace file
    IF(fTrc>0) THEN
       IF(PRESENT(Str)) THEN
          WRITE(fTrc,'(3A)') TRIM(Header) , ' : ', TRIM(Str)
       ELSE
          WRITE(fTrc,'(A)') TRIM(Header) 
       END IF
    END IF
    
    !// Screen
    IF(PRESENT(Str)) THEN
      WRITE(*,'(3A)')    TRIM(Header) , ' : ', TRIM(Str)
      IF (iDebug>1) CALL Pause_
    ELSE
       WRITE(*,'(A)')    TRIM(Header) 
    END IF

  END SUBROUTINE Message_

  !---

  SUBROUTINE GetUnit_Trace(F) 
    !===========================================================
    !returns a free unit number.
    !-------------------------------
    !Author:John Burkardt
    !A "free" FORTRAN unit number is an integer between 1 and 99
    !which is not currently associated with an I/O device.
    !A free FORTRAN unit number is needed in order to open a 
    !file with the OPEN command.
    !IUNIT=0
    !no free FORTRAN unit could be found, although all 99 units
    !were checked
    !===========================================================
    IMPLICIT NONE
    INTEGER,INTENT(OUT)::F
    INTEGER:: i, ios
    LOGICAL:: lOpen
    F= 0
    DO i=1,99
       IF (i/=5 .AND. i/=6 .AND. i/=9) THEN
          !units 5-6-9 are commonly reserved for console I/O
          INQUIRE(UNIT=i,OPENED=lopen,IOSTAT=ios)
          IF (ios==0) THEN
             IF (.NOT. lOpen) THEN
                F=i
                RETURN
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    RETURN
  ENDSUBROUTINE GetUnit_Trace


ENDMODULE M_Trace

