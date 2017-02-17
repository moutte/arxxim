module M_FileUtils
  !========================================================
  ! M_Files_Utils
  !--------------------------------------------------------
  ! Library of utils for file edition  
  !========================================================
  implicit none
  private

  !// Public Methods
  public :: File_Exist
  public :: File_Exist_System
  public :: File_Write_Date_Time
  public :: File_Path
  public :: Unix_FileName
  
contains

  !--

   function File_Exist_System(Str) 
    !--=====================================
    ! reports whether a file exists 
    ! and modify the name ( / <-> \ ) 
    ! to adapt to systems if necessary
    !--=====================================
    logical :: File_Exist_System
    character(len=*),intent(inout)::Str
    !---
    character(len=255) ::StrOut
    !! character(len=255) ::StrUnix
    !! character(len=8)   ::FS
    logical :: OK !!, OkUnixFS
    !---
    StrOut = File_Path(trim(Str))
    Ok = File_Exist(StrOut)
    
    !-- 
    File_Exist_System = Ok
    Str = StrOut
  
  end function File_Exist_System

  !---
  
  function File_Exist(Str) 
    !--=====================================
    ! reports whether a file exists
    !--=====================================
    implicit none
    logical :: File_Exist
    character(len=*),intent(in)::Str
    !---
    inquire(file=Str,EXIST=File_Exist)
    return  
  end function File_Exist

  !--

  subroutine File_Write_Date_Time(F)
    !--======================================
    ! write the date and time to a file unit
    !--======================================
    implicit none
    integer,intent(in):: F
    character(len=10):: Date_,Time_
    !
    call DATE_AND_TIME(DATE=Date_,TIME=Time_)
    write(F,'(5(A1,A))') &
         & "!",Date_(1:4), &
         & "/",Date_(5:6), &
         & "/",Date_(7:8), &
         & "/",Time_(1:2), &
         & "H",Time_(3:4)
  end subroutine File_Write_Date_Time

  !--

  subroutine Windows_Filename(Str)
    use M_Trace, only: SLASH_, BACKSLASH_
    implicit none
    character(len=*), intent(inout) :: Str
    !--
    call Char_Replace(Str,SLASH_,BACKSLASH_)
  end subroutine Windows_Filename

  !--

  subroutine Unix_FileName(Str)
    use M_Trace, only: SLASH_, BACKSLASH_
    implicit none
    character(len=*), intent(inout) :: Str
    !--
    call Char_Replace(Str,BACKSLASH_,SLASH_)
  end subroutine Unix_FileName

  !--

  subroutine Char_Replace(Str, a, b)
    character(len=*), intent(inout) :: Str
    character, intent(in) :: a, b
    !---
    integer :: i
    !---
    do i=1, len(Str)
       if (Str(i:i) == a) Str(i:i) = b
    end do

  end subroutine Char_Replace

  
  !---

  function File_Path(StrIn) 
    !=========================================================
    ! Return a FullFileName for the Current FileSystem
    ! Default Mode = UnixFS    ( Slash )
    ! Other Mode   = WindowsFS ( Backslash )
    !=========================================================
    implicit none
    character(len=*),intent(in) :: StrIn
    character(len=len(StrIn))   :: File_Path
    !---
    character(len=len(StrIn))   :: StrOut    
    logical :: OkUnixFS
    !---
    OkUnixFS = Is_UnixFS()
    StrOut=StrIn
    if ( OkUnixFS ) then
       call Unix_Filename(StrOut)
    else 
       call Windows_Filename(StrOut)
    end if

    File_Path = StrOut

    !// debug
    !!if (iDebug>0) write(*,*) "UnixFileSystem =", OkUnixFS, ", Input = ", StrIn, " Result = ", StrOut
    
  end function File_Path

   !---

  function Is_UnixFS() result(Ok)
    !======================================================
    ! Check if the File System is a UnixFS 
    ! Test if SLASH_ = "/" can be used in a FullFileName
    !======================================================
    use M_IOTools, only : GetUnit
    implicit none
    character,parameter::SLASH_     = Achar(47)
    !---
    logical :: Ok
    !---
    integer :: f
    character(len=7) :: FileName
    character(len=9) :: FullFileName
    !---
    call GetUnit(f)
    FileName = 'zzz_tmp'
    open(unit = f, file = FileName)   
    close(f)
    
    FullFileName='.'//SLASH_//FileName
    inquire(file=FullFileName,EXIST=Ok)
    
    open(unit = f, file = FileName, status='unknown') 
    close(f, status='DELETE')
    
  end function Is_UnixFS

  !--

  function Is_WindowsFS() result(Ok)
    !======================================================
    ! Check if the File System is a WindowsFS 
    ! Test if BACKSLASH_ = "\" can be used in a FullFileName
    !======================================================
    use M_IOTools, only : GetUnit
    implicit none
    character,parameter::BACKSLASH_ = Achar(92)
    !---
    logical :: Ok
    !---
    integer :: f
    character(len=7) :: FileName
    character(len=9) :: FullFileName
    !---
    call GetUnit(f)
    FileName = 'zzz_tmp'
    open(unit = f, file = FileName) 
    close(F)
        
    FullFileName='.'//BACKSLASH_//FileName
    inquire(file=FullFileName,EXIST=Ok)
    
    open(unit = f, file = FileName) !, status='SCRATCH') 
    close(f)
    
  end function Is_WindowsFS


end module M_FileUtils
