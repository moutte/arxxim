MODULE M_Files_Index
  !========================================================
  ! Index of output files produced during Arxim Simulation
  ! Implementation : single html flat file
  !========================================================
  IMPLICIT NONE
  PRIVATE

  !// Private Variables
  CHARACTER(LEN=30) :: Index_Filename = "output.html"
  
  PUBLIC :: Files_Index_Open
  PUBLIC :: Files_Index_Close
  PUBLIC :: Files_Index_Write
  PUBLIC :: Files_Index_Include
  PUBLIC :: Files_Index_InputFile

CONTAINS

  !---

  SUBROUTINE Files_Index_Open(f)
    !========================================================
    ! Open a new Index File
    !========================================================
    USE M_IOTools,ONLY: GetUnit
    USE M_FileUtils,ONLY: File_Write_Date_Time
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: f
    !--
    CALL GetUnit(f)
    OPEN(f,FILE=Index_Filename)
    WRITE(f,'(A)') "<html>"
    WRITE(f,'(A)') "<body>"
    WRITE(f,'(A)') "<P> <h2>Arxim Index Of Results</h2></P>"
    WRITE(f,'(A)') "<P>"
    CALL File_Write_Date_Time(f)
    WRITE(f,'(A)') "</P>"
    WRITE(f,'(A)') "<hr>"
  END SUBROUTINE Files_Index_Open

  !---

  SUBROUTINE Files_Index_Close(f)
    !========================================================
    ! Close an Index File
    !========================================================
    IMPLICIT NONE
    INTEGER,INTENT(INOUT):: f
    !--
    WRITE(f,'(A)') "</html>"
    WRITE(f,'(A)') "</body>"
    CLOSE(f)
    f= 0
  ENDSUBROUTINE Files_Index_Close

  
  !---

  SUBROUTINE Files_Index_Include(f,sTarget)
    !========================================================
    ! Add an index for an included File
    !========================================================
    IMPLICIT NONE
    INTEGER,INTENT(IN):: f
    CHARACTER(*),INTENT(IN):: sTarget
    !---
    IF(f>0) THEN
        WRITE(f,'(A)') '<P> File included: '
        WRITE(f,'(A)') '<a href="' //TRIM(sTarget) //'"><b>' //TRIM(sTarget) //'</b></a>'
        WRITE(f,'(A)') '</P>'
     ENDIF
  END SUBROUTINE Files_Index_Include

  !---

  SUBROUTINE Files_Index_InputFile(f,sTarget)
    !========================================================
    ! Add an index for an included File
    !========================================================
    USE M_FileUtils, only : Unix_Filename
    IMPLICIT NONE
    INTEGER,INTENT(IN):: f
    CHARACTER(*),INTENT(IN):: sTarget
    CHARACTER(LEN=255):: S
    !---
    S = sTarget
    CALL Unix_Filename(S)
    IF(f>0) THEN
        WRITE(f,'(A)') '<P> Input file: '
        WRITE(f,'(A)') '<a href="' //TRIM(S) //'"><b>' //TRIM(S) //'</b></a>'
        WRITE(f,'(A)') '</P>'
     ENDIF
  END SUBROUTINE Files_Index_InputFile

  !---

  SUBROUTINE Files_Index_Write(f,sTarget,sComment)
    !========================================================
    ! Add an index for an output file
    !========================================================
    USE M_FileUtils, only : Unix_Filename
    IMPLICIT NONE
    INTEGER,INTENT(IN):: f
    CHARACTER(*),INTENT(IN):: sTarget,sComment
    !--
    CHARACTER(LEN=255):: S
    !--
    S= TRIM(sTarget)
    IF(f>0) THEN
       CALL Unix_FileName(S)
       WRITE(f,'(A)') '<P> Output result:'
       WRITE(f,'(A)') '<a href="' //TRIM(S) //'"><b>' //TRIM(S) //'</b></a>,'
       WRITE(f,'(A)') TRIM(sComment)
       WRITE(f,'(A)') '</P>'
    ENDIF
  ENDSUBROUTINE Files_Index_Write

END MODULE M_Files_Index
