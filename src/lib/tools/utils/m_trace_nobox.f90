module M_Trace
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
  implicit none
  !
  private
  !
  !// parameterS
  character(len=20), parameter :: Default_TraceFileName =  "debug_all.log"
  character(len=9),  parameter :: Default_ErrorFileName =  "error.log"

  !// public data
  character,parameter,public:: T_         = Achar(9)
  character,parameter,public:: TAB_       = Achar(9)
  character,parameter,public:: SLASH_     = Achar(47)
  character,parameter,public:: BACKSLASH_ = Achar(92)

  integer,public:: iDebug= -1
  integer,public:: fTrc=   0
  integer,public:: fHtm=   0
  integer,public:: fError= 0
  logical,public:: DebugCoores=.false.

  logical,public:: LInfo= .false.
  logical,public:: LWarning= .false.

  !// public functionS
  public:: Trace_Init
  public:: Trace_Close
  !
  public:: Debug_
  public:: Pause_
  public:: Info_
  public:: Fatal_
  public:: Stop_
  public:: Warning_
  public:: Message_

  !// private functionS
  private:: GetUnit_Trace
  !
  !// private data
  character(len=80):: TraceFileName = Default_TraceFileName
  
contains

  !---

  subroutine Trace_Reset()
    implicit none
    call Trace_Close
    call Trace_Init(TraceFileName)
    
  end subroutine Trace_Reset

  !---

  subroutine Trace_Init(Str)
    implicit none
    character(len=*),intent(in),optional:: Str
    !--
    character(len=80) :: FileName
    !--
    FileName= Default_TraceFileName
    if(present(Str)) FileName= Str 
    
    if(fTrc==0) then
       !// open a new file
       call GetUnit_Trace(fTrc)
       TraceFileName= FileName
       open(fTrc,file=trim(FileName))
    else
       !// nothing to do
       call Warning_("TraceFile Already Initialized")
    end if
   
  end subroutine Trace_Init

  !---

  subroutine Trace_Close
    implicit none
    if(fTrc>0) then; close(fTrc)
       fTrc= 0
    end if
  end subroutine Trace_Close

  !---

  subroutine Pause_(Str)
    implicit none
    character(len=*),optional,intent(in):: Str
    !--
    if(present(Str)) then 
      call Message_("Pause",Str)
    end if
    write(*,'(A)') "... type return to continue." 
    read(*,*)

  end subroutine Pause_

  !---

  subroutine Warning_(Str)
    implicit none
    character(len=*),optional,intent(in):: Str
    !--
    if(present(Str)) then 
      if (LWarning) call Message_("Warning",Str)
    end if

  end subroutine Warning_

  !---

  subroutine Debug_(Str) 
    implicit none
    character(len=*),intent(in)::Str
    !--
    if (idebug>1) call Message_("Debug",Str)

  end subroutine Debug_

  !---

  subroutine Info_(Str) 
    implicit none
    character(len=*),intent(in)::Str
    !--
    if (LInfo) call Message_("Info",Str)

  end subroutine Info_

  !---

  subroutine Fatal_(Str)
    ! use M_BOX_COUPLER_VARS
    implicit none
    character(len=*),intent(in)::Str
    !--
    !if (LCouplerActive) then
    !  call Message_("Warning",Str)
    !  IerrorChemistry = 1
    !else
      call Message_("Fatal Error",Str)
      stop "FATAL ERROR"
    !end if
  end subroutine Fatal_

  !---

  subroutine Stop_(Str) 
    implicit none
    character(len=*),intent(in)::Str
    !--
    call Message_("Program Stop",Str)
    
    if(iDebug>0) then
      if(fError==0) then
        call GetUnit_Trace(fError)
        open(fError,file="error.log")
      end if
      write(fError,'(A)') trim(Str)
      close(fError)
    end if

    stop "program stop" 
    
  end subroutine Stop_

  !---

  subroutine Message_(Header, Str)
    implicit none
    character(len=*),intent(in)::Header
    character(len=*),intent(in), optional :: Str

    !// Trace file
    if(fTrc>0) then
       if(present(Str)) then
          write(fTrc,'(3A)') trim(Header) , ' : ', trim(Str)
       else
          write(fTrc,'(A)') trim(Header) 
       end if
    end if
    
    !// Screen
    if(present(Str)) then
      write(*,'(3A)')    trim(Header) , ' : ', trim(Str)
      if (iDebug>1) call Pause_
    else
       write(*,'(A)')    trim(Header) 
    end if

  end subroutine Message_

  !---

  subroutine GetUnit_Trace(F) 
    !===========================================================
    !returns a free unit number.
    !-------------------------------
    !Author:John Burkardt
    !A "free" FORTRAN unit number is an integer between 1 and 99
    !which is not currently associated with an I/O device.
    !A free FORTRAN unit number is needed in order to open a 
    !file with the open command.
    !IUNIT=0
    !no free FORTRAN unit could be found, although all 99 units
    !were checked
    !===========================================================
    implicit none
    integer,intent(out)::F
    integer:: i, ios
    logical:: lOpen
    F= 0
    do i=1,99
       if (i/=5 .and. i/=6 .and. i/=9) then
          !units 5-6-9 are commonly reserved for console I/O
          inquire(UNIT=i,opened=lopen,iostat=ios)
          if (ios==0) then
             if (.not. lOpen) then
                F=i
                return
             end if
          end if
       end if
    end do
    return
  end subroutine GetUnit_Trace


end module M_Trace

