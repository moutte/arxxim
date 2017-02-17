module M_IO_Menu
!.tools for input from KeyBoard / Screen

  implicit none
  private
  !
  integer :: ISelection = 0 ! selection number in file
  !
  public :: IO_Menu
  public :: IO_Options_Read
  !
  logical:: Ok_Path=      .false.
  logical:: Ok_Speciation=.false.
  logical:: Ok_TPpath=    .false.
  logical:: Ok_Simplex=   .false.
  logical:: Ok_Dynamic=   .false.
  logical:: Ok_Database=  .false.
  logical:: Ok_MixModel=  .false.
  
contains

!---

subroutine IO_Menu(S,OkCmd,DebugLevel)
  use M_FArgC,only: F_IARGC, F_GETARG
  use M_IoTools,only: WrdToInt,CarToInt,Str_Upper
  use M_VarStr
  !
  implicit none
  !---
  character(len=*),intent(out):: S
  logical,         intent(out):: OkCmd
  integer,         intent(out):: DebugLevel
  !---
  character(len=1) :: T
  integer          :: J
  logical          :: OkFInn
  
  DebugLevel= 0

  if(IARGC()>1) then
    
    ! 2 arguments on command line -> 2nd arg' is the debug level
    call GETARG(2,S)
    write(*,*) "READING ARG2 = ", S
    if(len_trim(S)==1) then
      J= CarToInt(S(1:1))
      if(J>=0) DebugLevel=J
      OkCmd=.false.
    else
      OkCmd=.true.
    end if

  else
    
    OkCmd=.false.
    
  end if

  if(.not. OkCmd) then
    
    !call Options_Read(DebugLevel)
    call IO_FInn_Selection(S, OkFInn)
    
    if (.not. OkFInn) then
      call IO_Menu_Selection(S)
    else
      write(*,'(A,A)') "TEST COMPUTE ", S
    end if
    
  end if

  call Str_Upper(S)
  
end subroutine IO_Menu

!---

subroutine IO_FInn_Selection(S, Ok)
  !===========================================================
  ! Select option number I in the TEST block from file FInn
  !===========================================================
  use M_Test_Read
  implicit none
  character(len=*),intent(out):: S
  logical,intent(out) :: Ok
  !----
  !! integer :: nRun = 0
  !---
  !S= ""
  if(ISelection==0) then
    call Test_Read_Init(Ok) ! should be called before !
    if(.not.Ok) return
  end if
  !
  ISelection = ISelection + 1
  Ok = Test_Read_OptCompute(ISelection, S)
  
  !~ print *,"IO_FInn_Selection=",S
  !~ pause
  !
end subroutine IO_FInn_Selection

!---

subroutine IO_Menu_Selection(S)
  !========================================================
  ! Print user menu selection and get user run option => S
  !========================================================
  use M_VarStr
  implicit none
  character(len=*),intent(out):: S
  !---
  type(T_VarStr)   :: Str
  !---
  print '(/,A,/)',"INTERACTIVE MENU selectION"

  if(Ok_Speciation) then
  
    print   '(A)',   "SPC:    Fluid Speciation"
    print   '(A)',   "EQU:    Equilibrium with other phases"
    
    if(Ok_TPpath) then
      print   '(A)',   "SPCTP:  Speciation  along T,P sequence"
      print   '(A)',   "EQUTP:  Equilibrium along T,P sequence"
    end if
    
    if(Ok_Path) then
      print '(A)',   "SPCPATH:  Speciation along Path"
      print '(A)',   "EQUPATH:  Speciation along Path"
    end if
    
    if(Ok_Dynamic) &
    & print '(A)',   "DYN:    Speciation & Dynamic"
    
    print   '(A)',   ""
    
  end if

  ! if(Ok_Database) then
  !   print '(A,/)',     "Database tests :"
  !   print '(A  )',     "  LOGK:      Build LogK database"
  !   print '(A  )',     "  DTBSHO:    Write Test Files"
  !   print '(A  )',     "  DTBFLUID:  Write G,H,S,Rho of H2O on T,P grid"
  !   print '(A  )',     "  DTBH2OHKF: Write solvent prop's of H2O on T,P grid"
  !   print '(A  )',     "  DTBAQU:    Write aqu'species prop's on T,P grid"
  ! end if

  ! if(Ok_MixModel) then
  !   print '(A)',     "  SOL:       Test mixing model"
  ! end if

  if(Ok_Simplex) then
     print '(/,A,/)', "Simplex tests :"
     print '(  A  )', "  SPLTP:   Simplex on TP path"
     print '(  A  )', "  SPLPATH: Simplex on composition path"
     print '(  A  )', "  SPLMIX:  Simplex with mixtures"
  end if
  !
  print   '(A)',     ""
  !
  print   '(A,/)',   "NEW:  New Input File"
  print   '(A,/)',   "REF:  Refresh"
  print   '(A,/)',   "Q:QUIT"
  !
  write(*,'(A)',advance="no") "Select : "
  call VarStr_Get(String=Str)

  S=VarStr_Char(Str)

end subroutine IO_Menu_Selection

!--

subroutine IO_Options_Read(DebugLevel)
  integer,intent(out):: DebugLevel
  call Options_Read(DebugLevel)
end subroutine

!--

subroutine Options_Read(DebugLevel)
  !========================================================
  ! Read available blocks in the input file
  ! Then set status of the different menu options
  ! Ok_Option <=> Option available in the menu
  !========================================================
  use M_IOTools !, only:dimV,LinToWrd,GetUnit
  use M_Files, only: NamFInn
  !use M_System_Vars,only: System_Type
  implicit none
  !
  integer,intent(out):: DebugLevel
  !
  character(len=512):: L,W
  logical:: EoL
  integer:: F,ios
  !
  Ok_Path=      .false.
  Ok_Speciation=.false.
  Ok_Simplex=   .false.
  Ok_Dynamic=   .false.
  Ok_MixModel=  .false.
  !
  call GetUnit(F)
  open(F,file=trim(NamFInn))

  DoFile: do

    read(F,'(A)',iostat=ios) L; if(ios/=0) exit DoFile
    call LinToWrd(L,W,EoL)
    if(W(1:1)=='!') cycle DoFile !skip comment lines
    call AppendToEnd(L,W,EoL)

    select case(trim(W))

    case("ENDINPUT")
      exit DoFile
    case("TP.TABLE")
      Ok_TPpath= .true.
    case("DEBUG")
      call LinToWrd(L,W,EoL)
      call WrdToInt(W,DebugLevel)
    case("ELEMENT")
      Ok_Database= .true.
    case("DYNAMIC")
      Ok_Dynamic= .true.
    case("PATH")
      Ok_Path=.true.
    case("SOLUTION.MODEL","MIXTURE.MODEL")
      Ok_MixModel=.true.
    case("SYSTEM","SYSTEM.AQUEOUS")
      Ok_Speciation=.true.
      !System_Type="AQUEOUS"
    case("SYSTEM.GLOBAL")
      Ok_Speciation=.true.
      !System_Type="GLOBAL" !not used yet
    case("SYSTEM.MOMAS")
      Ok_Speciation=.true.
      !System_Type="MOMAS"  !for "virtual" system, cf momas benches
    case("SYSTEM.SIMPLEX")
      Ok_Simplex= .true.
      !System_Type="SIMPLEX"

    end select

  end do DoFile

  close(F)
end subroutine Options_Read

end module M_IO_Menu
