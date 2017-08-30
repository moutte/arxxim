module M_Files_Index
  !========================================================
  ! Index of output files produced during Arxim Simulation
  ! Implementation : single html flat file
  !========================================================
  implicit none
  private

  !// Private Variables
  character(len=30) :: Index_Filename = "output.html"
  
  public :: Files_Index_Open
  public :: Files_Index_Close
  public :: Files_Index_Write
  public :: Files_Index_Include
  public :: Files_Index_InputFile

contains

  !---

  subroutine Files_Index_Open(f)
    !========================================================
    ! Open a new Index File
    !========================================================
    use M_IOTools,only: GetUnit
    use M_FileUtils,only: File_Write_Date_Time
    implicit none
    integer,intent(out):: f
    !--
    call GetUnit(f)
    open(f,file=Index_Filename)
    write(f,'(A)') "<html>"
    write(f,'(A)') "<body>"
    write(f,'(A)') "<P> <h2>Arxim Index Of Results</h2></P>"
    write(f,'(A)') "<P>"
    call File_Write_Date_Time(f)
    write(f,'(A)') "</P>"
    write(f,'(A)') "<hr>"
  end subroutine Files_Index_Open

  !---

  subroutine Files_Index_Close(f)
    !========================================================
    ! Close an Index File
    !========================================================
    implicit none
    integer,intent(inout):: f
    !--
    write(f,'(A)') "</html>"
    write(f,'(A)') "</body>"
    close(f)
    f= 0
  end subroutine Files_Index_Close

  
  !---

  subroutine Files_Index_Include(f,sTarget)
    !========================================================
    ! Add an index for an included File
    !========================================================
    implicit none
    integer,intent(in):: f
    character(*),intent(in):: sTarget
    !---
    if(f>0) then
        write(f,'(A)') '<P> File included: '
        write(f,'(A)') '<a href="' //trim(sTarget) //'"><b>' //trim(sTarget) //'</b></a>'
        write(f,'(A)') '</P>'
     end if
  end subroutine Files_Index_Include

  !---

  subroutine Files_Index_InputFile(f,sTarget)
    !========================================================
    ! Add an index for an included File
    !========================================================
    use M_FileUtils, only : Unix_Filename
    implicit none
    integer,intent(in):: f
    character(*),intent(in):: sTarget
    character(len=255):: S
    !---
    S = sTarget
    call Unix_Filename(S)
    if(f>0) then
        write(f,'(A)') '<P> Input file: '
        write(f,'(A)') '<a href="' //trim(S) //'"><b>' //trim(S) //'</b></a>'
        write(f,'(A)') '</P>'
     end if
  end subroutine Files_Index_InputFile

  !---

  subroutine Files_Index_Write(f,sTarget,sComment)
    !========================================================
    ! Add an index for an output file
    !========================================================
    use M_FileUtils, only : Unix_Filename
    implicit none
    integer,intent(in):: f
    character(*),intent(in):: sTarget,sComment
    !--
    character(len=255):: S
    !--
    S= trim(sTarget)
    if(f>0) then
       call Unix_FileName(S)
       write(f,'(A)') '<P> Output result:'
       write(f,'(A)') '<a href="' //trim(S) //'"><b>' //trim(S) //'</b></a>,'
       write(f,'(A)') trim(sComment)
       write(f,'(A)') '</P>'
    end if
  end subroutine Files_Index_Write

end module M_Files_Index
