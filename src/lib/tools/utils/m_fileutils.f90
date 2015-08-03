MODULE M_FileUtils
  !========================================================
  ! M_Files_Utils
  !--------------------------------------------------------
  ! Library of utils for file edition  
  !========================================================
  IMPLICIT NONE
  PRIVATE

  !// Public Methods
  PUBLIC :: File_Exist
  PUBLIC :: File_Exist_System
  PUBLIC :: File_Write_Date_Time
  PUBLIC :: File_Path
  PUBLIC :: Unix_FileName
  
CONTAINS

  !--

   FUNCTION File_Exist_System(Str) 
    !=======================================
    ! reports whether a file exists 
    ! and modify the name ( / <-> \ ) 
    ! to adapt to systems if necessary
    !=======================================
    USE M_Trace, only : iDebug
    IMPLICIT NONE
    LOGICAL :: File_Exist_System
    CHARACTER(LEN=*),INTENT(INOUT)::Str
    !---
    CHARACTER(LEN=255) ::StrOut
    !! CHARACTER(LEN=255) ::StrUnix
    !! CHARACTER(LEN=8)   ::FS
    LOGICAL :: OK !!, OkUnixFS
    !---
    StrOut = File_Path(TRIM(Str))
    Ok = File_Exist(StrOut)
    
    !-- 
    File_Exist_System = Ok
    Str = StrOut
  
  ENDFUNCTION File_Exist_System

  !---
  
  FUNCTION File_Exist(Str) 
    !=======================================
    ! reports whether a file exists
    !=======================================
    IMPLICIT NONE
    LOGICAL :: File_Exist
    CHARACTER(LEN=*),INTENT(IN)::Str
    !---
    INQUIRE(FILE=Str,EXIST=File_Exist)
    RETURN  
  ENDFUNCTION File_Exist

  !--

  SUBROUTINE File_Write_Date_Time(F)
    !========================================
    ! write the date and time to a file unit
    !========================================
    IMPLICIT NONE
    INTEGER,INTENT(IN):: F
    CHARACTER(LEN=10):: Date_,Time_
    !
    CALL DATE_AND_TIME(DATE=Date_,TIME=Time_)
    WRITE(F,'(5(A1,A))') &
         & "!",Date_(1:4), &
         & "/",Date_(5:6), &
         & "/",Date_(7:8), &
         & "/",Time_(1:2), &
         & "H",Time_(3:4)
  END SUBROUTINE File_Write_Date_Time

  !--

  SUBROUTINE Windows_Filename(Str)
    USE M_Trace, only: SLASH_, BACKSLASH_
    IMPLICIT NONE
    CHARACTER(LEN=*), intent(inout) :: Str
    !--
    CALL Char_Replace(Str,SLASH_,BACKSLASH_)

  END SUBROUTINE Windows_Filename

  !--

  SUBROUTINE Unix_FileName(Str)
    USE M_Trace, only: SLASH_, BACKSLASH_
    IMPLICIT NONE
    CHARACTER(LEN=*), intent(inout) :: Str
    !--
    CALL Char_Replace(Str,BACKSLASH_,SLASH_)

  END SUBROUTINE Unix_FileName

  !--

  SUBROUTINE Char_Replace(Str, a, b)
    CHARACTER(LEN=*), intent(inout) :: Str
    CHARACTER, intent(in) :: a, b
    !---
    integer :: i
    !---
    do i=1, len(Str)
       if (Str(i:i) == a) Str(i:i) = b
    end do

  END SUBROUTINE Char_Replace
  
  !---

  FUNCTION File_Path(StrIn) 
    !=========================================================
    ! Return a FullFileName for the Current FileSystem
    ! Default Mode = UnixFS    ( Slash )
    ! Other Mode   = WindowsFS ( Backslash )
    !=========================================================
    USE M_Trace, only : iDebug
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: StrIn
    CHARACTER(LEN=LEN(StrIn))   :: File_Path
    !---
    CHARACTER(LEN=LEN(StrIn))   :: StrOut    
    LOGICAL :: OkUnixFS
    !---
    OkUnixFS = Is_UnixFS()
    StrOut=StrIn
    IF ( OkUnixFS ) THEN
       CALL Unix_Filename(StrOut)
    ELSE 
       CALL Windows_Filename(StrOut)
    END IF

    File_Path = StrOut

    !// debug
    !!if (iDebug>0) WRITE(*,*) "UnixFileSystem =", OkUnixFS, ", Input = ", StrIn, " Result = ", StrOut
    
  END FUNCTION File_Path

   !---

  FUNCTION Is_UnixFS() result(Ok)
    !======================================================
    ! Check if the File System is a UnixFS 
    ! Test if SLASH_ = "/" can be used in a FullFileName
    !======================================================
    USE M_IOTools, only : GetUnit
    IMPLICIT NONE
    CHARACTER,PARAMETER::SLASH_     = ACHAR(47)
    !---
    LOGICAL :: Ok
    !---
    INTEGER :: f
    CHARACTER(LEN=7) :: FileName
    CHARACTER(LEN=9) :: FullFileName
    !---
    CALL GetUnit(f)
    FileName = 'zzz_tmp'
    OPEN(unit = f, file = FileName)   
    CLOSE(f)
    
    FullFileName='.'//SLASH_//FileName
    INQUIRE(FILE=FullFileName,EXIST=Ok)
    
    OPEN(unit = f, file = FileName, status='unknown') 
    CLOSE(f, status='DELETE')
    
  END FUNCTION Is_UnixFS

  !--

  FUNCTION Is_WindowsFS() result(Ok)
    !======================================================
    ! Check if the File System is a WindowsFS 
    ! Test if BACKSLASH_ = "\" can be used in a FullFileName
    !======================================================
    USE M_IOTools, only : GetUnit
    IMPLICIT NONE
    CHARACTER,PARAMETER::BACKSLASH_ = ACHAR(92)
    !---
    LOGICAL :: Ok
    !---
    INTEGER :: f
    CHARACTER(LEN=7) :: FileName
    CHARACTER(LEN=9) :: FullFileName
    !---
    CALL GetUnit(f)
    FileName = 'zzz_tmp'
    OPEN(unit = f, file = FileName) 
    CLOSE(F)
        
    FullFileName='.'//BACKSLASH_//FileName
    INQUIRE(FILE=FullFileName,EXIST=Ok)
    
    OPEN(unit = f, file = FileName) !, status='SCRATCH') 
    CLOSE(f)
    
  END FUNCTION Is_WindowsFS


END MODULE M_FileUtils
