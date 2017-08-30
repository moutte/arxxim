module M_Box_Read

  use M_Kinds
  use M_Trace

  use M_IOTools  
  use M_SetOption

  implicit none

  private

  !// Public Functions

  public:: Box_Read_Time
  public:: Box_Read_Reactor
  public:: Box_Read_Numeric
  public:: Box_Read_Rock
  public:: Box_Read_Debug
  public:: Box_Read_Lumping

contains

!---

subroutine Box_Read_Time(Ok)

!=========================================================
! Read Time Parameters
! BOX.TIME 
!=========================================================

  use M_Files, only: NamFInn
  use M_Box_TimeLoop_Vars

  implicit none

  !--- arguments
  logical, intent(out) :: Ok

  !--- local variables
  character(len=255):: L
  character(len=80) :: W,W1,W2,WW
  logical           :: EoL
  integer           :: F,ios

  !---------------------------------------------------------
  call Info_("Box_Read_Time")

  Ok=.false.

  !// Set Default Values
  call Info_("Set Default Values")

  call SetOption("TUNIT",  "SECOND", TUnit )
  call SetOption_Time("Tfinal", "0.    ", TFinal, TUnit )
  call SetOption_Time("DTINIT", "0.    ", dTInit, TUnit )
  call SetOption_Time("DTMIN",  "0.    ", dTMin,  TUnit )
  call SetOption_Time("DTMAX",  "0.    ", dTMax,  TUnit )

  call Info_("Read Bloc")

  !// Open Input File
  call GetUnit(F)
  call OpenFile(F,file=trim(NamFInn))
  !
  Do01: do 

    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit Do01
    call LinToWrd(L,W,EoL)

    if(W(1:1)=='!')   cycle Do01 !skip comment lines

    if(W=="END" .and. .not. Eol ) then
      call LinToWrd(L,WW,EoL)
      W=trim(W)//trim(WW)
    end if

    select case(trim(W))

    case("ENDINPUT")
      
      exit  Do01

    case("BOX.TIME")
      
      Ok=.true.
      Do02: do
      read(F,'(A)',iostat=ios) L
      if(ios/=0) exit Do01

      call LinToWrd(L,W1,EoL)
      call AppendToEnd(L,W1,EoL)

      if(W1(1:1)=='!') cycle Do02
      call LinToWrd(L,W2,EoL)

      !-> W2 is second word on line
      !-> either second keyword or numeric' param' (or 'empty', i.e. "!")

      select case(trim(W1)) 

      case("ENDINPUT")
        exit Do01

      case("END","ENDBOX.TIME")
        exit Do02

      case("TUNIT")
      select case(trim(W2)) 
      case("SECOND"); call SetOption("TUNIT", "SECOND", TUnit )
      case("HOUR");   call SetOption("TUNIT", "HOUR",   TUnit )
      case("DAY" );   call SetOption("TUNIT", "DAY",    TUnit )
      case("YEAR");   call SetOption("TUNIT", "YEAR",   TUnit )
      case default
      call FATAL_("BOX.TIME|TUNIT " //trim(W2)//" unknown keyword")
      
      end select

    case("Tfinal") ; call SetOption_Time("Tfinal", W2,  TFinal, TUnit )
    case("DTINIT") ; call SetOption_Time("DTINIT", W2,  dTInit, TUnit ) 
    case("DTMIN")  ; call SetOption_Time("DTMIN",  W2,  dTMin,  TUnit ) 
    case("DTMAX")  ; call SetOption_Time("DTMAX",  W2,  dTMax,  TUnit ) 

    case default
      call FATAL_("BOX.TIME = " //trim(W1)//" unknown keyword")

    end select
    
  enddo Do02

  end select
  
  enddo Do01

!// call closeFILE Input File
call closeFILE(F)

call Info_("Box_ReadTime Done")

end subroutine Box_Read_Time

!---

subroutine Box_Read_Reactor(Ok)

!=========================================================
! Read Box Reactor Parameters
! BOX.REACTOR
!=========================================================

  use M_Files,only: NamFInn

  use M_Box_Vars
  use M_Box_Thermo_Vars
  use M_Box_TimeLoop_Vars

  implicit none

  !-- arguments
  logical,       intent(out)  :: Ok

  !-- local variables
  character(len=255):: L
  character(len=80) :: W,W1,W2
  logical           :: EoL
  integer           :: F,ios

  !-------------------------------------------------------
  call Info_("Box_Read_Reactor")

  Ok=.false.

  !// Default Values
  call Info_("Set Default Values")

  call SetOption( "RHOF",    "1000.  ",   RhoF )
  call SetOption( "PHif",    "0.5    ",   PhiF )
  call SetOption( "VBOX",    "1.     ",   Vbox )
  call SetOption_FlowRate( "FOUT",    "0.     ",   Fout, TUnit )
  call SetOption_FlowRate( "FINJ",    "0.     ",   Finj, TUnit ) 
  call SetOption( "POROSITY","0.5    ",    PhiF )  

  call Info_("Read Bloc")
  !// Open Input File
  call GetUnit(F)
  call OpenFile(F,file=trim(NamFInn))

  Do01: do 
    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit Do01

    call LinToWrd(L,W,EoL)
    call AppendToEnd(L,W,EoL)

    if(W(1:1)=='!') cycle Do01 !skip comment lines

    select case(trim(W))

    case("ENDINPUT")
      exit  Do01

    case("BOX.REACTOR")
      Ok=.true.
      
      Do02: do
      
        read(F,'(A)',iostat=ios) L
        if(ios/=0) exit Do01

        call LinToWrd(L,W1,EoL)
        call AppendToEnd(L,W1,EoL)

        if(W1(1:1)=='!') cycle Do02

        call LinToWrd(L,W2,EoL)

        !-> W2 is second word on line
        !-> either second keyword or numeric' param' (or 'empty', i.e. "!")

        select case(trim(W1))

        case("ENDINPUT")
          exit Do01

        case("END","ENDBOX.REACTOR")
          exit Do02

        case("RHOF") ;     call SetOption( "RHOF",    W2, RhoF )
        case("PHif") ;     call SetOption( "PHif",    W2, PhiF )
        case("VBOX") ;     call SetOption( "VBOX",    W2, Vbox )
        case("FOUT") ;     call SetOption_FlowRate( "FOUT",    W2, Fout, TUnit )
        if (Fout < Zero) call Fatal_("Fout Should be positive !")
        case("FINJ")  ;    call SetOption_FlowRate( "FINJ",    W2, Finj, TUnit )
        if (Finj > Zero) call Fatal_("Finj Should be negative !")
        case("POROSITY") ; call SetOption( "POROSITY",W2, PhiF )  

        case default
        call Fatal_("BOX.REACTOR = "//trim(W1)//" unknown keyword")

        end select

      enddo Do02

    end select
  
  enddo Do01

  !// call closeFILE Input File
  call closeFILE(F)

  call Info_("Box_ReadReactor Done")

end subroutine Box_Read_Reactor

!---

subroutine Box_Read_Numeric(Ok)
!=========================================================
! Read Box Numeric Parameters
! BOX.NUMERIC
!=========================================================

  use M_Files,only: NamFInn
  use M_Box_Solver_Vars
  use M_Box_Newton_Vars
  use M_Box_Vars

  implicit none

  !-- arguments
  logical,intent(out)::Ok

  !-- local variables
  character(len=255):: L
  character(len=80) :: W,W1,W2
  logical           :: EoL
  integer           :: F,ios

  !----------------------------------------------------------
  call Info_("Box_Read_Numeric")
  Ok=.false.

  !// Default Values
  call Info_("Set Default Values")

  call SetOption("METHOD",   "NEWTLNSRCH", cMethod)
  call SetOption("MAXITER",  "100     ", NewtMaxIts)  
  call SetOption("ITERMAX",  "100     ", NewtIterMax) 
  call SetOption("ITERMIN",  "1       ", NewtIterMin) 
  call SetOption("NEWTTOLF", "1.D-6   ", NewtTolF) 
  call SetOption("NEWTTOLX", "1.D-10  ", NewtTolX) 
  call SetOption("FINDif",   .false.   , bFinDif) 
  call SetOption("TESTMAX",  .false.   , TestMax) 
  call SetOption("implicit", .true.    , Implicit_Surface)
  call SetOption("LOGFORMIN",   .false.    , LogForMin)
  call SetOption("LOGFORFLUID", .true.     , LogForFluid)
  call SetOption("EPSMINERAL", "1.D-09 "   , EpsMineral)

  call Info_("Read Bloc")

  !// Open Input File
  call GetUnit(F)
  call OpenFile(F,file=trim(NamFInn))

  Do01: do
  
    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit Do01

    call LinToWrd(L,W,EoL)
    call AppendToEnd(L,W,EoL)

    if(W(1:1)=='!')   cycle Do01 !skip comment lines

    SEL0: select case(trim(W))

    case("ENDINPUT"); exit  Do01

    case("BOX.NUMERIC")
      Ok=.true.
      
      Do02: do
      
      read(F,'(A)',iostat=ios) L
      if(ios/=0) exit Do01

      call LinToWrd(L,W1,EoL)
      call AppendToEnd(L,W1,EoL)

      if(W1(1:1)=='!') cycle Do02

      call LinToWrd(L,W2,EoL)

      !-> W2 is second word on line
      !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
      !

      SEL1: select case(trim(W1)) 

      case("ENDINPUT")
        exit Do01

      case("END","ENDBOX.NUMERIC")
        exit Do02 

      case("METHOD")  !!then
        select case(trim(W2))
        !
        case("NEWTON");       call SetOption("METHOD", "NEWTON",      cMethod) 
        case("NEWTLNSRCH");   call SetOption("METHOD", "NEWTLNSRCH",  cMethod)
        case("BROYDEN");      call SetOption("METHOD", "BROYDEN",     cMethod)
        case("NEWTONCHESS");  call SetOption("METHOD", "NEWTONCHESS", cMethod)
        case("TENSOLVE");     call SetOption("METHOD", "TENSOLVE",    cMethod)
        case("NEWTONWALKER"); call SetOption("METHOD", "NEWTONWALKER",cMethod)
        case("NEWTONKELLEY"); call SetOption("METHOD", "NEWTONKELLEY",cMethod)
        case default
        call Fatal_("BOX.NUMERIC|METHOD = "// trim(W2)//" unknown keyword")
        !
        end select
      
      case("MAXITER");  call SetOption("MAXITER",  W2, NewtMaxIts)  
      case("ITERMAX");  call SetOption("ITERMAX",  W2, NewtIterMax) 
      case("ITERMIN");  call SetOption("ITERMIN",  W2, NewtIterMin) 
      case("NEWTTOLF"); call SetOption("NEWTTOLF", W2, NewtTolF) 
      case("NEWTTOLX"); call SetOption("NEWTTOLX", W2, NewtTolX) 
      case("EPSMINERAL"); call SetOption("EPSMINERAL", W2, EpsMineral) 

      case("FINDif") ; call SetOption("FINDif",  .true.,  bFinDif) 
      case("TESTMAX"); call SetOption("TESTMAX", .false., TestMax) 

      case("LOGFORMIN")  ; call SetOption("LOGFORMIN",   .true.    , LogForMin)
      case("LOGFORFLUID"); call SetOption("LOGFORFLUID", .true.    , LogForFluid)

      case("SURFACE")
        select case(trim(W2))
        case("implicit"); call SetOption("implicit", .true.,  Implicit_Surface)
        case("EXPLICIT"); call SetOption("implicit", .false., Implicit_Surface)
        case default
        call Fatal_("BOX.NUMERIC|SURFACE = "// trim(W2)//" unknown keyword")
        end select

      case default
        call Fatal_("BOX.NUMERIC = "//trim(W1)//" unknown keyword")
      
      end select SEL1

      enddo Do02

    end select SEL0
  
  enddo Do01

  !// call closeFILE File
  call closeFILE(F)

  call Info_("Box_ReadNumeric Done")

end subroutine Box_Read_Numeric

!---

subroutine Box_Read_Rock(Ok)
!=========================================================
! Read Box Rock Parameters
! DYNAMIC.ROCK 
!=========================================================
  use M_T_KinModel
  use M_Global_Vars, only: vSpc,vFas,vKinFas,vKinModel 
  use M_KinFas_Read, only: T_LnkKin,KinFas_BuildLnk,KinFas_LnkToVec 
  use M_Dynam_Vars,  only: sModelSurf
  use M_KinModel_Read
  use M_Box_Debug_Vars, only : LDebug_KinRate

  logical,intent(out)::Ok
  type(T_LnkKin),pointer:: LnkKin
  real(dp):: x
  integer :: I,N
  !

  !---
  call Info_("Box_ReadRock")

  Ok=.false.

  !//======================================
  !// Read The Bloc KINETICS
  !//======================================
  call KinModel_Init(vSpc) 

  !//======================================
  !// Read The Bloc DYNAMIC ROCK
  !//======================================
  call KinFas_BuildLnk( &  
  & vFas, &           ! IN
  & vKinModel, &      ! IN
  & sModelSurf, &     ! IN
  !---
  & N, &              ! OUT
  & LnkKin )          ! OUT

  if (N>0) then 
  Ok=.true. 
  if(allocated(vKinFas)) deallocate(vKinFas)
  allocate(vKinFas(N))

  call KinFas_LnkToVec ( LnkKin, vFas, vKinFas) 

  else 
  Ok=.false.
  call FATAL_("NO Rock Found")
  return 
  end if

  !// Info
  call Info_("Dynamic Rocks")

  if (LDebug_KinRate )  then 
    print '(5A16)', & 
    & "Kinetic_Phase___", &
    & "Thermo_Model____", &
    & "Kinetic_Model___", &
    & "________Radius__", &
    & "________SurfKg__"


    do I=1,size(vKinFas) 
      if(vKinFas(I)%iKin>0) &
      print '(3(A15,1X),2(G15.3,1X))', & 
      & vKinFas(I)%NamKF, &
      & vFas(vKinFas(I)%iFas)%NamFs, &
      & vKinModel(vKinFas(I)%iKin)%Name, &
      & vKinFas(I)%Dat%Radius, &
      & vKinFas(I)%Dat%SurfKg 
    enddo
  end if
  !//======================================
  !// Normalize PhiM to SUM(VolFract)=1
  !//======================================

  x = SUM(vKinFas(:)%Dat%PhiM, MASK=vKinFas(:)%Dat%cSat /= "MINIMAL") 
  if (x>Zero) then 
    do I=1,size(vKinFas) 
      if (vKinFas(I)%Dat%cSat /= "MINIMAL") then 
        vKinFas(I)%Dat%PhiM= vKinFas(I)%Dat%PhiM/x 
      else 
        vKinFas(I)%Dat%PhiM= Zero 
      end if
    enddo
  else 
    call Fatal_("Normalize KinFas PhiM : Sum(PhiM) less than Zero") 
  end if

  call Info_("Box_ReadRock Done")

end subroutine Box_Read_Rock

!---

subroutine Box_Read_Debug(Ok)

!=========================================================
! Read Debug Parameters
! BOX.DEBUG
!=========================================================

  use M_Files, only: NamFInn
  use M_Box_Debug_Vars

  implicit none

  !--- arguments
  logical, intent(out) :: Ok

  !--- local variables
  character(len=255):: L
  character(len=80) :: W,W1,W2,WW
  logical           :: EoL
  integer           :: F,ios

  !---------------------------------------------------------
  call Info_("Box_Read_Debug")

  Ok=.false.

  !// Set Default Values
  call Info_("Set Default Values")

  call SetOption("SYSTEM",   .false., LDebug_System) 
  call SetOption("INIT",     .false., LDebug_Init)
  call SetOption("CALC",     .false., LDebug_Calc)
  call SetOption("KINTHERMO",.false., LDebug_KinThermo)
  call SetOption("KINRATE",    .false., LDebug_KinRate)
  call SetOption("GAMMA",    .false., LDebug_Gamma)
  call SetOption("REFSTATE", .false., LDebug_RefState)
  call SetOption("INJECT",   .false., LDebug_Inject)
  call SetOption("VOLUME",   .false., LDebug_Volume)
  call SetOption("ENTRIES",  .false., LDebug_Entries)
  call SetOption("VARIABLES",.false., LDebug_Variables)
  call SetOption("RESIDUAL", .false., LDebug_Residual)
  call SetOption("JACOBIAN", .false., LDebug_Jacobian)
  call SetOption("NEWTON",   .false., LDebug_Newton)
  call SetOption("INFO",     .false., LDebug_Info)
  call SetOption("WARNING",  .false., LDebug_Warning)
  call SetOption("TIMELOOP", .false., LDebug_TimeLoop)
  call SetOption("TIMESTEP", .false., LDebug_TimeStep)
  call SetOption("CHECK",    .false., LDebug_Check)
  call SetOption("TPUPDATE", .false., LDebug_TPUpdate)
  call SetOption("LUMPING" , .false., LDebug_Lumping)
  call Info_("Read Bloc")

  !// Open Input File
  call GetUnit(F)
  call OpenFile(F,file=trim(NamFInn))
  !
  Do01: do 

    read(F,'(A)',iostat=ios) L
    if(ios/=0) exit Do01
    call LinToWrd(L,W,EoL)

    if(W(1:1)=='!')   cycle Do01 !skip comment lines

    if(W=="END" .and. .not. Eol ) then
    call LinToWrd(L,WW,EoL)
    W=trim(W)//trim(WW)
    end if

    select case(trim(W))

    case("ENDINPUT")
    exit  Do01

    case("BOX.DEBUG")
    Ok=.true.
    
    Do02: do
    
      read(F,'(A)',iostat=ios) L
      if(ios/=0) exit Do01

      call LinToWrd(L,W1,EoL)
      call AppendToEnd(L,W1,EoL)

      if(W1(1:1)=='!') cycle Do02
        call LinToWrd(L,W2,EoL)

      !-> W1 is second word on line
      !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
      select case(trim(W1))

      case("ENDINPUT")
        exit Do01

      case("END","ENDBOX.DEBUG")
        exit Do02

      case("SYSTEM")    ; call SetOption("SYSTEM",   .true., LDebug_System) 
      case("INIT")      ; call SetOption("INIT",     .true., LDebug_Init)
      case("CHECK")     ; call SetOption("CHECK",    .true., LDebug_Check)
      case("CALC")      ; call SetOption("CALC",     .true., LDebug_Calc)
      case("KINTHERMO") ; call SetOption("KINTHERMO",.true., LDebug_KinThermo)
      case("KINRATE")  ;  call SetOption("KINRATE", .true., LDebug_KinRate)
      case("GAMMA")     ; call SetOption("GAMMA",    .true., LDebug_Gamma)
      case("REFSTATE")  ; call SetOption("REFSTATE", .true., LDebug_RefState)
      case("INJECT")    ; call SetOption("INJECT",   .true., LDebug_Inject)
      case("VOLUME")    ; call SetOption("VOLUME",   .true., LDebug_Volume)
      case("ENTRIES")   ; call SetOption("ENTRIES",  .true., LDebug_Entries)
      case("VARIABLES") ; call SetOption("VARIABLES",.true., LDebug_Variables)
      case("RESIDUAL")  ; call SetOption("RESIDUAL", .true., LDebug_Residual)
      case("JACOBIAN")  ; call SetOption("JACOBIAN", .true., LDebug_Jacobian)
      case("NEWTON")    ; call SetOption("NEWTON",   .true., LDebug_Newton)
      case("INFO")      ; call SetOption("INFO",     .true., LDebug_Info)
      case("WARNING")   ; call SetOption("WARNING",  .true., LDebug_Warning)
      case("TIMELOOP")  ; call SetOption("TIMELOOP", .true., LDebug_TimeLoop)
      case("TIMESTEP")  ; call SetOption("TIMESTEP", .true., LDebug_TimeStep)
      case("TPUPDATE")  ; call SetOption("TPUPDATE", .true., LDebug_TPUpdate)   
      case("LUMPING")   ; call SetOption("LUMPING",  .true., LDebug_Lumping)
      
      case default
      
      call Fatal_("BOX.DEBUG = "//trim(W1)//" unknown keyword")             
      
      end select
    
    enddo Do02

    end select
  end do Do01

  !// call closeFILE Input File
  call closeFILE(F)

  call Info_("Box_ReadTime Done")

end subroutine Box_Read_Debug

!---

subroutine Box_Read_Lumping(Ok)
  use M_Box_Lumping
  implicit none
  logical :: Ok
  !---
  call Box_Lumping_Read (Ok)   

end subroutine Box_Read_Lumping


end module M_Box_Read
