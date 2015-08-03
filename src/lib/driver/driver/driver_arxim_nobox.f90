SUBROUTINE Driver_Arxim_NoBox
  !---------------------------------------------------------------------
  ! Driver_Arxim
  !---------------------------------------------------------------------
  ! Interactive Program Arxim
  ! Get Input Filename and Compute Mode From [ User Command Line ]
  !---------------------------------------------------------------------
  USE M_Trace,       ONLY: iDebug, Stop_, Trace_Close, Trace_Init,Pause_,LWarning
  USE M_Files,       ONLY: Files_PathRead, Files_BuildInput 
  USE M_Global_Tools,ONLY: Global_Zero, Global_Clean
  USE M_System_Vars, ONLY: System_Clean, vCpn
  USE M_IO_Menu
  !
  IMPLICIT NONE
  !
  LOGICAL:: OkCmd,Ok
  LOGICAL:: OkSMode
  CHARACTER(LEN=15):: S
  !
  CALL Trace_Init
  !
  CALL Files_PathRead(Ok)
  IF(.NOT. Ok) CALL Stop_("ERRORS")
  ! 
  CALL Files_BuildInput
  DoMain: DO
    !
    CALL Global_Zero
    !
    CALL IO_Menu(S,OkCmd,iDebug)
    LWarning= (iDebug>0)
    !
    !------------------------------------------------------
    !// SPECIAL MODES : Q, REF, NEW
    OkSMode = .FALSE.
    SELECT CASE(TRIM(S)) 
    CASE("Q") 
      OkSMode = .TRUE.
      EXIT DoMain
    CASE("NEW")
      CALL Files_PathRead(Ok)
      CALL Files_BuildInput
      OkSMode = .TRUE.
    CASE("REF") !to call when input file has been modified
      CALL Files_BuildInput !-> renew arxim_inn
      OkSMode = .TRUE.
    CASE DEFAULT
      OkSMode = .FALSE.
    END SELECT
    !------------------------------------------------------
    !// COMPUTE MODES : SPC, EQU, DYN, BOX, SPCPATH, etc ...
    IF (.NOT. OkSMode) THEN
      CALL Driver_Arxim_ComputeSequence(S, OkCmd, OkSMode)
    END IF
    !
    CALL System_Clean
    CALL Global_Clean
    !
    IF(OkCmd) EXIT
  ENDDO DoMain
  !
  CALL Trace_Close
  WRITE(*,'(A)') "PERFECT"

END SUBROUTINE Driver_Arxim_NoBox

!---

SUBROUTINE Driver_Arxim_Options(sFilename, sMode)
!---------------------------------------------------------------------
! Driver_Arxim_Options
!---------------------------------------------------------------------
! Batch Program Arxim 
! Get Input Filename and Compute Mode From [ Arguments ]
!---------------------------------------------------------------------
  USE M_Trace,       ONLY: iDebug, Stop_,Trace_Close,Trace_Init,Pause_,LWarning
  USE M_Files,       ONLY: Files_PathRead_From_FileName, Files_BuildInput 
  USE M_Global_Tools,ONLY: Global_Zero, Global_Clean
  USE M_System_Vars, ONLY: System_Clean, vCpn
  USE M_IO_Menu
  !
  IMPLICIT NONE
  !
  LOGICAL:: OkCmd,Ok, OkSMode
  CHARACTER(LEN=*):: sMode
  CHARACTER(LEN=*):: sFilename
  !
  CALL Trace_Init
  !
  CALL Files_PathRead_From_FileName(sFileName, Ok) 
  IF(.NOT. Ok) CALL Stop_("ERRORS")
  ! 
  CALL Files_BuildInput
  !
  CALL Global_Zero
  !
  OkCmd=.TRUE. ! Direct Command 
  CALL IO_Options_Read(iDebug)
  LWarning= (iDebug>0)
  !
  WRITE(*,'(A,A)') "TEST COMPUTE ", sMode
  !
  CALL Driver_Arxim_ComputeSequence(sMode, OkCmd, OkSMode)
  !
  CALL System_Clean
  CALL Global_Clean
  !
  CALL Trace_Close
  !
  if (OkSMode)      WRITE(*,'(A)') "PERFECT"
  if (.not.OkSMode) WRITE(*,'(A)') "ERRORS"

END SUBROUTINE Driver_Arxim_Options

!---

SUBROUTINE Driver_Arxim_ComputeSequence(S, OkCmd, OkSMode)
  
  USE M_Basis_Vars,  ONLY: Basis_CleanAll
  USE M_Basis,       ONLY: Basis_Change
  USE M_System_Vars, ONLY: System_Clean,vCpn
  USE M_Global_Build,ONLY: Global_Build
  USE M_System_Tools,ONLY: System_Build
  USE M_Dtb_Test
  USE M_Equil
  USE M_Dynam
  USE M_Dynam_Column
  USE M_Path
  USE M_GEM_Vars
  USE M_Simplex_Theriak
  USE M_Simplex_Path
  USE M_Simplex_Build
  USE M_DiscretModel_Test
  USE M_IO_Menu

  IMPLICIT NONE

  CHARACTER(len=*) :: S
  LOGICAL :: OkCmd
  LOGICAL, intent(out) :: OkSMode

  LOGICAL :: Ok
  
  !--- COMPUTE Sequence 
  OkSMode = .true.
  
  SELECT CASE(TRIM(S))
    !
  CASE DEFAULT !-> in case S is unknown as a code ...
    !IF(OkCmd) THEN
      PRINT &
      & '(/,A)',TRIM(S)//"-> Unknown Code In command line...???" 
      OkSMode = .false.
    !ELSE
    !  PRINT &
    !  & '(/,A)',"WHAT DO YOU SAY ?? RETRY ...???"
    !  OkSMode = .false.
    !ENDIF
     
  CASE("SPC")
    CALL Global_Build
    CALL System_Build
    !
    CALL Equil_Calc("SPC")
    !
    CALL Basis_CleanAll
    CALL System_Clean
    !
  CASE("EQU","EQ2")
    CALL Global_Build
    CALL System_Build
    !
    !IF(COUNT(vCpn(:)%Statut=="MOBILE")>0) THEN
    IF(COUNT(vCpn(:)%Statut=="MOBILE")>0 .OR. &
    &  COUNT(vCpn(:)%Statut=="BUFFER")>0) THEN
      CALL Equil_Calc("SPC")
      CALL Basis_Change("EQU",vCpn)
    ENDIF
    CALL Equil_Calc("EQ2")
    !
    CALL Basis_CleanAll
    CALL System_Clean
    !
  CASE("EQM")
    CALL Global_Build
    CALL System_Build
    !
    !IF(COUNT(vCpn(:)%Statut=="MOBILE")>0) THEN
    IF(COUNT(vCpn(:)%Statut=="MOBILE")>0 .OR. &
    &  COUNT(vCpn(:)%Statut=="BUFFER")>0) THEN
      CALL Equil_Calc("SPC")
      CALL Basis_Change("EQU",vCpn)
    ENDIF
    CALL Equil_Calc("EQM")
    !
    CALL Basis_CleanAll
    CALL System_Clean
    !
  CASE("EQ1")
    CALL Global_Build
    CALL System_Build
    !
    IF(COUNT(vCpn(:)%Statut=="MOBILE")>0 .OR. &
    &  COUNT(vCpn(:)%Statut=="BUFFER")>0) THEN
    !IF(COUNT(vCpn(:)%Statut=="MOBILE")>0) THEN
      CALL Equil_Calc("SPC")
      CALL Basis_Change("EQU",vCpn)
    ENDIF
    CALL Equil_Calc("EQ1")
    !
    CALL Basis_CleanAll
    CALL System_Clean

  CASE("DYN")
    CALL Global_Build
    !
    CALL Dynam_Initialize
    !-> initial speciation for box fluid and inject fluid
    !
    CALL Dynam_Box
    !
    CALL Basis_CleanAll
    CALL System_Clean
    
  CASE("COLUMN")
    CALL Global_Build
    !
    CALL Dynam_Initialize
    !-> initial speciation for box fluid and inject fluid
    !
    CALL Dynam_Column
    !
    CALL Basis_CleanAll
    CALL System_Clean
  
  CASE("SPCTP","TPSPC")
    CALL Global_Build
    CALL System_Build
    !
    CALL Path_Execute("SPCTP")
    !
    CALL Basis_CleanAll
    CALL System_Clean
    
  CASE("EQUTP","TPEQU","EQ2TP")
    CALL Global_Build
    CALL System_Build
    !
    CALL Path_Execute("EQ2TP")
    !
    CALL Basis_CleanAll
  
  CASE("EQ1TP")  
    CALL Global_Build
    CALL System_Build
    !
    CALL Path_Execute("EQ1TP")
    !
    CALL Basis_CleanAll
    
  ! CASE("EQ0TP")  
  !   CALL Global_Build
  !   CALL System_Build
  !   !
  !   CALL Path_Execute("EQ0TP")
  !   !
  !   CALL Basis_CleanAll
    
  CASE("EQMTP")  
      CALL Global_Build
      CALL System_Build
      !
      CALL Path_Execute("EQMTP")
      !
      CALL Basis_CleanAll
      
  CASE("SPCPATH","PATHSPC")
    CALL Global_Build
    CALL System_Build
    !
    CALL Path_Execute("SPC__")
    !
    CALL Basis_CleanAll
  CASE("EQUPATH","PATHEQU","EQ2PATH")
    CALL Global_Build
    CALL System_Build
    !
    CALL Path_Execute("EQ2__")
    !
    CALL Basis_CleanAll
    !
  CASE("EQ1PATH")
    CALL Global_Build
    CALL System_Build
    !
    CALL Path_Execute("EQ1__")
    !
    CALL Basis_CleanAll
    !
  CASE("EQMPATH")
    CALL Global_Build
    CALL System_Build
    !
    CALL Path_Execute("EQM__")
    !
    CALL Basis_CleanAll
    !
  ! CASE("SPLMIX")
  !   CALL Global_Build
  !   CALL Simplex_Build
  !   CALL Simplex_Theriak(.FALSE.)
  
  ! CASE("SPLMIXID")
  !   CALL Global_Build
  !   CALL Simplex_Build
  !   CALL Simplex_Theriak(.TRUE.)
  
  CASE("GEM")
    CALL Global_Build
    CALL Simplex_Build
    CALL Simplex_Theriak
    CALL GEM_Vars_Clean
  
  ! CASE("GEMID")
  !   CALL Global_Build
  !   CALL Simplex_Build
  !   CALL Simplex_Theriak(.TRUE.)
  
  CASE("GEMPATH")
    CALL Global_Build
    CALL Simplex_Build
    CALL Simplex_Theriak_Path
    CALL GEM_Vars_Clean
    
  CASE("SPLTP")
    CALL Global_Build
    CALL Simplex_Build
    CALL Simplex_Path("TP")
  
  CASE("SPLPATH")
    CALL Global_Build
    CALL Simplex_Build
    CALL Simplex_Path("PATH")

  CASE("DTBEQ36")
    CALL Global_Build
    CALL Dtb_Test_EQ36
    
  CASE("DTBSOL")
    CALL Global_Build
    !! CALL Dtb_Test_Mixture_1
    !CALL Dtb_Test_Mixture_2
    !
  CASE("DTBFLUID")
    CALL Global_Build
    CALL Dtb_Test_Fluid
  
  CASE("DTBH2OHKF")
    CALL Global_Build
    CALL Dtb_Test_H2OHkf
    !
  CASE("DTBSPC")
    CALL Global_Build
    CALL Dtb_Test_Species
    !
  CASE("DTBAQU")
    CALL Global_Build
    CALL System_Build
    !-> will work only on system's species,
    !   because could produce too many data on whole base
    CALL Dtb_Test_AquHkf
    !
  CASE("DTBLOGK")
    CALL Global_Build
    CALL Dtb_Tabulate("LOGK")
    
  CASE("DTBGIBBS")
    CALL Global_Build
    CALL Dtb_Tabulate("GIBBS")
  
  CASE("DTBGIBBS2")
    CALL Global_Build
    CALL Dtb_Tabulate("GIBBS2")
    
  CASE("SSAS1","ssas1") !uses ssas01.inn as input !!!
    !
    CALL Global_Build
    CALL System_Build
    !!CALL Test_Ssas(1)
    !
  CASE("SSAS2","ssas2") !uses ssas02.inn as input !!!
    !
    CALL Global_Build
    CALL System_Build
    !!CALL Test_Ssas(2)
   ! !
  ENDSELECT

END SUBROUTINE Driver_Arxim_ComputeSequence
