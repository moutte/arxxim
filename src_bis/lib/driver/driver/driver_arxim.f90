subroutine Driver_Arxim
  !---------------------------------------------------------------------
  ! Driver_Arxim
  !---------------------------------------------------------------------
  ! Interactive Program Arxim
  ! Get Input Filename and Compute Mode From [ User Command Line ]
  !---------------------------------------------------------------------
  use M_Trace,       only: iDebug, Stop_, Pause_, LWarning
  use M_Trace,       only: Trace_Close, Trace_Init
  use M_Files,       only: Files_PathRead, Files_BuildInput 
  use M_Global_Tools,only: Global_Zero, Global_Clean
  use M_System_Vars, only: System_Clean, vCpn
  use M_IO_Menu
  !
  implicit none
  !
  logical:: OkCmd,Ok
  logical:: OkSMode
  character(len=15):: S
  !
  call Trace_Init
  !
  call Files_PathRead(Ok)
  if(.not. Ok) call Stop_("ERRORS")
  ! 
  call Files_BuildInput
  DoMain: do
    !
    call Global_Zero
    !
    call IO_Menu(S,OkCmd,iDebug)
    LWarning= (iDebug>0)
    !
    !------------------------------------------------------
    !// SPECIAL MODES : Q, REF, NEW
    OkSMode = .false.
    select case(trim(S)) 
    case("Q") 
      OkSMode = .true.
      exit DoMain
    case("NEW")
      call Files_PathRead(Ok)
      call Files_BuildInput
      OkSMode = .true.
    case("REF") !to call when input file has been MODIFIED
      call Files_BuildInput !-> renew arxim_inn
      OkSMode = .true.
    case default
      OkSMode = .false.
    end select
    !------------------------------------------------------
    !// COMPUTE MODES : SPC, EQU, DYN, BOX, SPCPATH, etc ...
    if (.not. OkSMode) then
      call Driver_Arxim_ComputeSequence(S, OkCmd, OkSMode)
    end if
    !
    call System_Clean
    call Global_Clean
    !
    if(OkCmd) exit
  end do DoMain
  !
  call Trace_Close
  write(*,'(A)') "PERFECT"

end subroutine Driver_Arxim

!---

subroutine Driver_Arxim_Options(sFilename, sMode)
!---------------------------------------------------------------------
! Driver_Arxim_Options
!---------------------------------------------------------------------
! Batch Program Arxim 
! Get Input Filename and Compute Mode From [ Arguments ]
!---------------------------------------------------------------------
  use M_Trace,       only: iDebug, Stop_, Pause_, LWarning
  use M_Trace,       only: Trace_Close,Trace_Init
  use M_Files,       only: Files_PathRead_From_FileName, Files_BuildInput 
  use M_Global_Tools,only: Global_Zero, Global_Clean
  use M_System_Vars, only: System_Clean, vCpn
  use M_IO_Menu
  !
  implicit none
  !
  logical:: OkCmd,Ok, OkSMode
  character(len=*):: sMode
  character(len=*):: sFilename
  !
  call Trace_Init
  !
  call Files_PathRead_From_FileName(sFileName, Ok) 
  if(.not. Ok) call Stop_("ERRORS")
  ! 
  call Files_BuildInput
  !
  call Global_Zero
  !
  OkCmd=.true. ! Direct Command 
  call IO_Options_Read(iDebug)
  LWarning= (iDebug>0)
  !
  write(*,'(A,A)') "TEST COMPUTE ", sMode
  !
  call Driver_Arxim_ComputeSequence(sMode, OkCmd, OkSMode)
  !
  call System_Clean
  call Global_Clean
  !
  call Trace_Close
  !
  if (OkSMode)      write(*,'(A)') "PERFECT"
  if (.not.OkSMode) write(*,'(A)') "ERRORS"

end subroutine Driver_Arxim_Options

!---

subroutine Driver_Arxim_ComputeSequence(S, OkCmd, OkSMode)
  !
  use M_Trace,       only: Pause_
  use M_Basis_Vars,  only: Basis_CleanAll
  use M_Basis,       only: Basis_Change
  use M_System_Vars, only: System_Clean,vCpn
  use M_Global_Build,only: Global_Build
  use M_System_Tools,only: System_Build
  use M_Dtb_Test
  use M_Equil
  use M_Dynam
  use M_Dynam_Column
  use M_Path
  use M_GEM_Vars
  use M_GEM_Theriak
  use M_GEM_Path
  use M_GEM_Build
  use M_DiscretModel_Test
  use M_IO_Menu
  !
  implicit none
  !
  character(len=*), intent(in):: S
  logical :: OkCmd
  logical, intent(out) :: OkSMode

  logical :: Ok
  
  !--- COMPUTE Sequence 
  OkSMode = .true.
  
  select case(trim(S))
  !
  case default !-> in case S is unknown as a code ...
    !if(OkCmd) then
      print &
      & '(/,A)',trim(S)//"-> Unknown Code In command line...???" 
      OkSMode = .false.
    !else
    !  print &
    !  & '(/,A)',"WHAT do YOU SAY ?? RETRY ...???"
    !  OkSMode = .false.
    !end if     
  !------------------------------------------------------------------SPC
  case("SPC")
    call Global_Build
    call System_Build
    !
    call Equil_Calc("SPC")
    !
    call Basis_CleanAll
    call System_Clean
    !
  case("SPCTP","TPSPC")
    call Global_Build
    call System_Build
    !
    call Path_Execute("SPCTP")
    !
    call Basis_CleanAll
    call System_Clean
    
  case("SPCPATH","PATHSPC")
    call Global_Build
    call System_Build
    !
    call Path_Execute("SPC__")
    !
    call Basis_CleanAll
  !-----------------------------------------------------------------/SPC
  !
  !------------------------------------------------------------------EQU
  case("EQU","EQ2")
    call Global_Build
    call System_Build
    !
    !if(count(vCpn(:)%Statut=="MOBILE")>0) then
    if(count(vCpn(:)%Statut=="MOBILE")>0 .or. &
    &  count(vCpn(:)%Statut=="BUFFER")>0) then
      call Equil_Calc("SPC")
      call Basis_Change("EQU",vCpn)
    end if
    call Equil_Calc("EQ2")
    !
    call Basis_CleanAll
    call System_Clean
    !
  case("EQM")
    call Global_Build
    call System_Build
    !
    !if(count(vCpn(:)%Statut=="MOBILE")>0) then
    if(count(vCpn(:)%Statut=="MOBILE")>0 .or. &
    &  count(vCpn(:)%Statut=="BUFFER")>0) then
      call Equil_Calc("SPC")
      call Basis_Change("EQU",vCpn)
    end if
    call Equil_Calc("EQM")
    !
    call Basis_CleanAll
    call System_Clean
    !
  case("EQ1")
    call Global_Build
    call System_Build
    !
    if(count(vCpn(:)%Statut=="MOBILE")>0 .or. &
    &  count(vCpn(:)%Statut=="BUFFER")>0) then
    !if(count(vCpn(:)%Statut=="MOBILE")>0) then
      call Equil_Calc("SPC")
      call Basis_Change("EQU",vCpn)
    end if
    call Equil_Calc("EQ1")
    !
    call Basis_CleanAll
    call System_Clean

  case("EQUTP","TPEQU","EQ2TP")
    call Global_Build
    call System_Build
    !
    call Path_Execute("EQ2TP")
    !
    call Basis_CleanAll
    !
  case("EQ1TP")  
    call Global_Build
    call System_Build
    !
    call Path_Execute("EQ1TP")
    !
    call Basis_CleanAll
    !
  ! case("EQ0TP")  
  !   call Global_Build
  !   call System_Build
  !   !
  !   call Path_Execute("EQ0TP")
  !   !
  !   call Basis_CleanAll
    
  case("EQMTP")  
      call Global_Build
      call System_Build
      !
      call Path_Execute("EQMTP")
      !
      call Basis_CleanAll
      
  case("EQUPATH","PATHEQU","EQ2PATH")
    call Global_Build
    call System_Build
    !
    call Path_Execute("EQ2__")
    !
    call Basis_CleanAll
    !
  case("EQ1PATH")
    call Global_Build
    call System_Build
    !
    call Path_Execute("EQ1__")
    !
    call Basis_CleanAll
    !
  case("EQMPATH")
    call Global_Build
    call System_Build
    !
    call Path_Execute("EQM__")
    !
    call Basis_CleanAll
  !-----------------------------------------------------------------/EQU
  
  !------------------------------------------------------------------DYN
  case("BOXCRES")
    call Global_Build 
    call Driver_Box_Cres

  case("BOX")
    call Global_Build 
    call Driver_Box

  case("DYN")
    call Global_Build
    !
    call Dynam_Initialize
    !-> initial speciation for box fluid and inject fluid
    !
    call Dynam_Box
    !
    call Basis_CleanAll
    call System_Clean
    
  case("COLUMN")
    call Global_Build
    !! print *,"Global_Build_end"  ; call pause_
    !
    call Dynam_Initialize
    !-> initial speciation for box fluid and inject fluid
    !
    call Dynam_Column
    !
    call Basis_CleanAll
    call System_Clean
  !-----------------------------------------------------------------/DYN
  
  !------------------------------------------------------------------GEM
  ! case("SPLMIX")
  !   call Global_Build
  !   call GEM_Build
  !   call GEM_Theriak_Single(.false.)
  
  ! case("SPLMIXID")
  !   call Global_Build
  !   call GEM_Build
  !   call GEM_Theriak_Single(.true.)
  
  case("GEM")
    call Global_Build
    call GEM_Build
    call GEM_Theriak_Single
    call GEM_Vars_Clean
  
  ! case("GEMID")
  !   call Global_Build
  !   call GEM_Build
  !   call GEM_Theriak_Single(.true.)
  
  case("GEMPATH")
    call Global_Build
    call GEM_Build
    call GEM_Theriak_Path
    call GEM_Vars_Clean
    
  case("SPLTP")
    call Global_Build
    call GEM_Build
    call GEM_Path("TP")
  
  case("SPLPATH")
    call Global_Build
    call GEM_Build
    call GEM_Path("PATH")
  !-----------------------------------------------------------------/GEM

  !------------------------------------------------------------------DTB
  case("EOSIFP")
    call test_EosIfp
     
  case("DTBEQ36")
    call Global_Build
    call Dtb_Test_EQ36
    
  case("DTBSOL")
    call Global_Build
    !! call Dtb_Test_Mixture_1
    !call Dtb_Test_Mixture_2
    !
  case("DTBFLUID")
    call Global_Build
    call Dtb_Test_Fluid
  
  case("DTBH2OHKF")
    call Global_Build
    call Dtb_Test_H2OHkf
    !
  case("DTBSPC")
    call Global_Build
    call Dtb_Test_Species
    !
  case("DTBAQU")
    call Global_Build
    call System_Build
    !-> will work only on system's species,
    !   because could produce too many data on whole base
    call Dtb_Test_AquHkf
    !
  case("DTBLOGK")
    call Global_Build
    call Dtb_Tabulate("LOGK")
    
  case("DTBVOLUME")
    call Global_Build
    call Dtb_Tabulate("VOLUME")
    
  case("DTBGIBBS")
    call Global_Build
    call Dtb_Tabulate("GIBBS")
  
  case("DTBGIBBS2")
    call Global_Build
    call Dtb_Tabulate("GIBBS2")
    
  case("SSAS1","ssas1") !uses ssas01.inn as input !!!
    !
    call Global_Build
    call System_Build
    !!call Test_Ssas(1)
    !
  case("SSAS2","ssas2") !uses ssas02.inn as input !!!
    !
    call Global_Build
    call System_Build
    !!call Test_Ssas(2)
   ! !
  end select

end subroutine Driver_Arxim_ComputeSequence
