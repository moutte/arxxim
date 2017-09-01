module M_Files

  use M_Files_Vars
  use M_Files_Index
  use M_FileUtils
  implicit none
  !
  private
  !
  !// Private variables
  character(len=80) :: NamFIn0

  !// Exported M_Files_Vars Variables
  public :: NamFInn
  public :: cTitle
  public :: NamFLogK,NamFEle,NamFKin,NamFSol,NamFPtz
  public :: NamDtbAqu,NamDtbMin,NamDtbMlt
  public :: DirOut,DirLog
  public :: DirDtbOut,DirDtbLog 

  !// Public Methods
  public :: Files_Vars_Init
  public :: Files_BuildInput
  public :: Files_PathRead
  public :: Files_PathRead_From_FileName
  public :: Info_PathRead
  
  ! Exported M_Files_Index methods
  public :: Files_Index_Write
  public :: Files_Index_Open
  public :: Files_Index_Close

  ! Exported M_Files_Utils methods
  public :: File_Exist
  public :: File_Write_Date_Time 
  public :: File_Path

  !// Private Methods
  private :: Files_PathRead_From_CommandArgs
  private :: Files_PathRead_From_UserInput
  
contains

  subroutine Info_PathRead
    !===================================================================
    ! Print the Input FileName 
    !===================================================================
    use M_Trace,  only: fHtm,iDebug
    implicit none
    !--
    if(iDebug>2) &
    & write(*,'(A)') "Arxim Script File : "// trim(NamFIn0)
    call Files_Index_InputFile(fHtm, NamFIn0)

  end subroutine Info_PathRead

  !---

  subroutine Files_PathRead(Ok)
    !===================================================================
    ! Read the Input FileName ( Default = "arxim.inn" )
    !===================================================================
    use M_Trace, only: Stop_
    implicit none
    logical, intent(out) :: Ok
    !--
    Ok = .false.
    if (.not.Ok) call Files_PathRead_From_CommandArgs(Ok)
    if (.not.Ok) call Files_PathRead_From_FileName("arxim.inn",Ok)
    if (.not.Ok) call Files_PathRead_From_UserInput(Ok)
    if (.not.Ok) call Stop_("Error Reading Input File")

    call Info_PathRead
    
  end subroutine Files_PathRead

  !---

  subroutine Files_BuildInput
    !===================================================================
    ! Build the input data file from main file + include
    ! Init Files Variables
    !===================================================================
    implicit none
    !---
    ! Merge include files
    call Files_BuildInput_IncludeFiles

    ! Init files names variable
    call Files_Vars_Init(NamFInn)

  end subroutine Files_BuildInput

  !---

  subroutine Files_PathRead_From_FileName(sFileName, Ok) 
    !===================================================================
    ! Read the Input FileName from the command line argumants
    ! -> Default Extension = .inn
    ! Output in NamFIn0 
    !===================================================================
    implicit none
    logical,intent(out):: Ok
    character(len=*), intent(in):: sFileName
    !--
    Ok=.false.
    NamFIn0="NOFILE"
    
    if(.not. File_Exist(trim(sFileName))) then
       Ok=.false.
       NamFIn0="NOFILE"
       print '(A)',trim(sFileName)//"-> File Not Found ...!!!"
    else
       Ok=.true.
       NamFIn0=trim(sFileName)
    end if
    
  end subroutine Files_PathRead_From_FileName

  !---

  subroutine Files_PathRead_From_CommandArgs(Ok) 
    !===================================================================
    ! Read the Input FileName from the command line argumants
    ! -> Default Extension = .inn
    ! Output in NamFIn0 
    !===================================================================
    use M_FArgC, only: F_GETARG, F_IARGC
    use M_VarStr,only: T_VarStr,VarStr_Get,VarStr_Char
    implicit none
    !--
    logical,intent(out):: Ok
    !--
    character(len=255):: Str
    !--
    Ok=.false.
    NamFIn0="NOFILE"

    if(F_IARGC()>0) then
    
      call F_GETARG(1,Str)
      
      if(INDEX(Str,".")<1) Str=trim(Str)//".inn"
      ! Default file name suffix is .inn
      
      if(.not. File_Exist(trim(Str))) then
        Ok=.false.
        print '(A)',trim(Str)//"-> File Not Found ...!!!"
      else
        Ok=.true.
        NamFIn0=trim(Str)
      end if
      
    end if
    
    return
  end subroutine Files_PathRead_From_CommandArgs

  !--

  subroutine Files_PathRead_From_UserInput(Ok) 
    !===================================================================
    ! Read the Input FileName by using interactive questions
    ! -> Default FileName  = arxim.inn
    ! -> Default Extension = .inn
    ! Output in NamFIn0 
    !===================================================================
    use M_VarStr,only: T_VarStr,VarStr_Get,VarStr_Char
    implicit none
    !--
    logical,intent(out):: Ok
    !--
    type(T_VarStr):: VStr
    character(len=255):: Str
    !--
    Ok=.false.
    NamFIn0="NOFILE"
    
    DoNam: do
       write(*,'(A)') "Input File Name [ Return = arxim.inn ] ?"
       call VarStr_Get(String=VStr)
       if(VarStr_Char(VStr)=="") then
          Str="arxim.inn"
       else
          Str=VarStr_Char(VStr)
          ! Default file name suffix is .inn 
          if(INDEX(Str,".")<1) Str=trim(Str)//".inn"        
       end if
       if(File_Exist(Str)) then
          Ok=.true.
          exit DoNam
       else                   
          Ok=.false.
          print '(A)',"File Not Found ...!!!"
       end if
    enddo DoNam
    NamFIn0=trim(Str)
    
  end subroutine Files_PathRead_From_UserInput

  !---

  subroutine Files_BuildInput_IncludeFiles
    !==========================================================
    ! Build the Input file NamFinn = "arxim_inn.tmp"
    !
    !  -> Merge nested include files 
    !  -> Skip comments
    ! 
    ! Input  File  = NamFIn0 ( + include Files )
    ! Output File  = NamFinn = "arxim_inn.tmp" 
    ! 
    !==========================================================
    use M_Trace,  only: iDebug,Pause_,fHtm,Stop_
    use M_IOTools,only: GetUnit,LinToWrd,AppendToEnd
    !
    integer:: fAll,fInn,fAdd,fAdd2
    character(len=512):: L, LL, W, W1, W2
    character(len=512):: sIncludeFile1, sIncludeFile2
    logical:: EoL
    integer:: ios     
    !
    NamFInn="arxim_inn.tmp"
    !
    call GetUnit(fAll); open(fAll,file=trim(NamFInn))
    call GetUnit(fInn); open(fInn,file=trim(NamFIn0))

    !=========== read ROOT FILE  ================================
    Do1: do
      !
      read(fInn,'(A)',iostat=ios) LL
      if(ios/=0) exit Do1 
      L=trim(LL)
      call LinToWrd(L,W,EoL)
      call AppendToEnd(L,W,EoL)
      !
      if(W(1:1)=='!') cycle Do1 ! skip comment lines
      if(trim(W)=='/*') then !skip comment block
        do
          read(fInn,'(A)',iostat=ios) L
          if(ios/=0) then
            print '(A)', "!!!WARNING!!! COMMENT block UNTERMINATED !!!WARNING!!!"
            exit Do1
          end if
          call LinToWrd(L,W,EoL)
          if(trim(W)=='*/') exit
        enddo
        cycle Do1
      end if
      !
      select case(trim(W))
      !
      case default
        write(fAll,'(A)') trim(LL)
      !
      case("ENDINPUT")
        exit Do1
        !
      case("INCLUDE") 
        ! Get Include File Name Level 1        
        call LinToWrd(L,W,EoL,"NO") 
        sIncludeFile1 = W
        
        if (File_Exist_System(sIncludeFile1)) then 
          
          if(iDebug>2) &
          & write(*,*) "INCLUDE LEVEL 1 : ", trim(sIncludeFile1)
          
          ! Open Include File Level 1
          call GetUnit(fAdd)
          open(fAdd,file=trim(sIncludeFile1))
          call Files_Index_Include(fHtm,sIncludeFile1)
          write(fAll,'(A,A)') "!!!!!!!!!!!!Insert File_begin!!!!!!!!!!!!"
          write(fAll,'(A,A)') "! include ", trim(sIncludeFile1)
          
          !=========== read FILE include LEVEL 1 =======================
          Do2: do
            !
            read(fAdd,'(A)',iostat=ios) LL
            if(ios/=0) exit Do2
            L=LL
            call LinToWrd(L,W1,EoL)
            call AppendToEnd(L,W1,EoL)
            !
            if(W1(1:1)=='!') cycle Do2 !skip comment lines
            if(trim(W1)=='/*') then !skip comment block
              do
                read(fInn,'(A)',iostat=ios) L
                if(ios/=0) then
                  print '(A)', &
                  & "!!!WARNING!!! COMMENT block UNTERMINATED !!!WARNING!!!"
                  exit Do2
                end if
                call LinToWrd(L,W1,EoL)
                if(trim(W1)=='*/') exit
              enddo
              cycle Do2
            end if
            !
            if(trim(W1)=="ENDINPUT")  exit Do2
            !
            if(trim(W1)=="INCLUDE") then 
              
              ! Get Include File Name Level 2
              call LinToWrd(L,W1,EoL,"NO") 
              sIncludeFile2 = W1
              
              if(File_Exist_System(sIncludeFile2)) then 
                
                if(iDebug>2) &
                & write(*,*) "INCLUDE LEVEL 2 : .. ", trim(sIncludeFile2)
                
                ! Open Include File Level 2
                call GetUnit(fAdd2)
                open(fAdd2,file=trim(sIncludeFile2))
                call Files_Index_Include(fHtm,sIncludeFile2)
                write(fAll,'(A)') "!!!!!!!!!!!!Insert File_begin!!!!!!!!!!!!"
                write(fAll,'(A,A)') "! include ", trim(sIncludeFile2)
                
                !=========== read FILE LEVEL 2 =========================
                do3: do 
                  read(fAdd2,'(A)',iostat=ios) LL
                  if(ios/=0) exit do3
                  L=trim(LL)
                  call LinToWrd(L,W2,EoL)
                  call AppendToEnd(L,W2,EoL)
                  !
                  if(W2(1:1)=='!') cycle do3 !skip comment lines
                  if(trim(W2)=='/*') then !skip comment block
                    do
                      read(fInn,'(A)',iostat=ios) L
                      if(ios/=0) then
                        print '(A)', &
                        & "!!!WARNING!!! COMMENT block UNTERMINATED !!!WARNING!!!"
                        exit Do3
                      end if
                      call LinToWrd(L,W2,EoL)
                      if(trim(W2)=='*/') exit
                    enddo
                    cycle do3
                  end if
                  !
                  if(trim(W2)=="ENDINPUT") exit do3
                  !
                  if(trim(W2)=="INCLUDE") then
                    print '(A)', "!!!WARNING!!! TOO Many Nested includes !!!WARNING!!!"
                    call LinToWrd(L,W2,EoL,"NO") 
                    call Stop_ ("INCLUDE LEVEL 3 :"//trim(W2)//" NOT IMPLEMENTED")
                  end if
                  !
                  write(fAll,'(A)') trim(LL)
                enddo do3
                
                ! Close Include File Level 2
                write(fAll,'(A)') "!!!!!!!!!!!!Insert File _end!!!!!!!!!!!!!"
                close(fAdd2)
              
              else
                
                call Stop_("File "//trim(sIncludeFile2)//" NOT FOUND -> Can Not include !!!")
              
              end if
              !
            end if ! nested include
            !
            write(fAll,'(A)') trim(LL)
          enddo Do2
          ! Close Include File Level 1
          write(fAll,'(A)') "!!!!!!!!!!!!Insert File _end!!!!!!!!!!!!!"
          close(fAdd)
        else
          call Stop_("FILE "//trim(sIncludeFile1)//" NOT FOUND -> Can Not include !!!")
        end if
      end select
    enddo Do1
    ! Close Root File
    close(fInn)
    ! Close Result Merged File
    !! write(fAll,'(A)') "ENDINPUT" !!! if(Mod="APPend") 
    close(fAll)

  end subroutine Files_BuildInput_IncludeFiles

  !---

  subroutine Files_Vars_Init(sFilNam)
    !==========================================================
    ! Init the Directories and Input Files Names
    !==========================================================
    implicit none
    character(*),intent(in):: sFilNam
    !---
    cTitle=""

    ! defaut directories
    DirOut="out_"
    DirLog="log_"

    DirDtbOut="dtb_out_"
    DirDtbLog="dtb_log_"

    ! make the input file the default file for all data
    NamFLogK=trim(sFilNam) ! file for logK
    NamFEle= trim(sFilNam) ! file for elements
    NamFSol= trim(sFilNam) ! file for solutions
    NamFKin= trim(sFilNam) ! file for kinetics
    NamFPtz= trim(sFilNam) ! file for pitzer parameters

  end subroutine Files_Vars_Init

end module M_Files

