MODULE M_IO_Menu
!.tools for input from KeyBoard / Screen

  IMPLICIT NONE
  PRIVATE
  !
  INTEGER :: ISelection = 0 ! selection number in file
  !
  PUBLIC :: IO_Menu
  PUBLIC :: IO_Options_Read
  !
  LOGICAL:: Ok_Path=      .FALSE.
  LOGICAL:: Ok_Speciation=.FALSE.
  LOGICAL:: Ok_TPpath=    .FALSE.
  LOGICAL:: Ok_Simplex=   .FALSE.
  LOGICAL:: Ok_Dynamic=   .FALSE.
  LOGICAL:: Ok_Database=  .FALSE.
  LOGICAL:: Ok_MixModel=  .FALSE.
  
CONTAINS

!---

SUBROUTINE IO_Menu(S,OkCmd,DebugLevel)
  USE M_FArgC,ONLY: F_IARGC, F_GETARG
  USE M_IoTools,ONLY: WrdToInt,CarToInt,Str_Upper
  USE M_VarStr
  !
  IMPLICIT NONE
  !---
  CHARACTER(LEN=*),INTENT(OUT):: S
  LOGICAL,         INTENT(OUT):: OkCmd
  INTEGER,         INTENT(OUT):: DebugLevel
  !---
  CHARACTER(LEN=1) :: T
  INTEGER          :: J
  LOGICAL          :: OkFInn
  
  DebugLevel= 0

  IF(IARGC()>1) THEN
    
    ! 2 arguments on command line -> 2nd arg' is the debug level
    CALL GETARG(2,S)
    WRITE(*,*) "READING ARG2 = ", S
    IF(LEN_TRIM(S)==1) THEN
      J= CarToInt(S(1:1))
      IF(J>=0) DebugLevel=J
      OkCmd=.FALSE.
    ELSE
      OkCmd=.TRUE.
    ENDIF

  ELSE
    
    OkCmd=.FALSE.
    
  ENDIF

  IF(.NOT. OkCmd) THEN
    
    !CALL Options_Read(DebugLevel)
    CALL IO_FInn_Selection(S, OkFInn)
    
    IF (.NOT. OkFInn) THEN
      CALL IO_Menu_Selection(S)
    ELSE
      WRITE(*,'(A,A)') "TEST COMPUTE ", S
    END IF
    
  END IF

  CALL Str_Upper(S)
  
ENDSUBROUTINE IO_Menu

!---

SUBROUTINE IO_FInn_Selection(S, Ok)
  !===========================================================
  ! Select option number I in the TEST block from file FInn
  !===========================================================
  USE M_Test_Read
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(OUT):: S
  LOGICAL,INTENT(OUT) :: Ok
  !----
  !! INTEGER :: nRun = 0
  !---
  !S= ""
  IF(ISelection==0) THEN
    CALL Test_Read_Init(Ok) ! should be called before !
    IF(.NOT.Ok) RETURN
  ENDIF
  !
  ISelection = ISelection + 1
  Ok = Test_Read_OptCompute(ISelection, S)
  
  !~ print *,"IO_FInn_Selection=",S
  !~ pause
  !
END SUBROUTINE IO_FInn_Selection

!---

SUBROUTINE IO_Menu_Selection(S)
  !========================================================
  ! Print user menu selection and get user run option => S
  !========================================================
  USE M_VarStr
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(OUT):: S
  !---
  TYPE(T_VarStr)   :: Str
  !---
  PRINT '(/,A,/)',"INTERACTIVE MENU SELECTION"

  IF(Ok_Speciation) THEN
  
    PRINT   '(A)',   "SPC:    Fluid Speciation"
    PRINT   '(A)',   "EQU:    Equilibrium with other phases"
    
    IF(Ok_TPpath) THEN
      PRINT   '(A)',   "SPCTP:  Speciation  along T,P sequence"
      PRINT   '(A)',   "EQUTP:  Equilibrium along T,P sequence"
    ENDIF
    
    IF(Ok_Path) THEN
      PRINT '(A)',   "SPCPATH:  Speciation along Path"
      PRINT '(A)',   "EQUPATH:  Speciation along Path"
    ENDIF
    
    IF(Ok_Dynamic) &
    & PRINT '(A)',   "DYN:    Speciation & Dynamic"
    
    PRINT   '(A)',   ""
    
  ENDIF

  !~ IF(Ok_Database) THEN
     !~ PRINT '(A,/)',     "Database tests :"
     !~ PRINT '(A  )',     "  LOGK:      Build LogK database"
     !~ PRINT '(A  )',     "  DTBSHO:    Write Test Files"
     !~ PRINT '(A  )',     "  DTBFLUID:  Write G,H,S,Rho of H2O on T,P grid"
     !~ PRINT '(A  )',     "  DTBH2OHKF: Write solvent prop's of H2O on T,P grid"
     !~ PRINT '(A  )',     "  DTBAQU:    Write aqu'species prop's on T,P grid"
  !~ ENDIF

  !~ IF(Ok_MixModel) THEN
     !~ PRINT '(A)',     "  SOL:       Test mixing model"
  !~ ENDIF

  IF(Ok_Simplex) THEN
     PRINT '(/,A,/)', "Simplex tests :"
     PRINT '(  A  )', "  SPLTP:   Simplex on TP path"
     PRINT '(  A  )', "  SPLPATH: Simplex on composition path"
     PRINT '(  A  )', "  SPLMIX:  Simplex with mixtures"
  ENDIF
  !
  PRINT   '(A)',     ""
  !
  PRINT   '(A,/)',   "NEW:  New Input File"
  PRINT   '(A,/)',   "REF:  Refresh"
  PRINT   '(A,/)',   "Q:QUIT"
  !
  WRITE(*,'(A)',ADVANCE='NO') "Select : "
  CALL VarStr_Get(String=Str)

  S=VarStr_Char(Str)

END SUBROUTINE IO_Menu_Selection

!--

SUBROUTINE IO_Options_Read(DebugLevel)
  INTEGER,INTENT(OUT):: DebugLevel
  CALL Options_Read(DebugLevel)
END SUBROUTINE

!--

SUBROUTINE Options_Read(DebugLevel)
  !========================================================
  ! Read available blocks in the input file
  ! Then set status of the different menu options
  ! Ok_Option <=> Option available in the menu
  !========================================================
  USE M_IOTools !, ONLY:dimV,LinToWrd,GetUnit
  USE M_Files, ONLY: NamFInn
  !USE M_System_Vars,ONLY: System_Type
  IMPLICIT NONE
  !
  INTEGER,INTENT(OUT):: DebugLevel
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL:: EoL
  INTEGER:: F,ios
  !
  Ok_Path=      .FALSE.
  Ok_Speciation=.FALSE.
  Ok_Simplex=   .FALSE.
  Ok_Dynamic=   .FALSE.
  Ok_MixModel=  .FALSE.
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFInn))

  DoFile: DO

    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)

    SELECT CASE(TRIM(W))

    CASE("ENDINPUT")
      EXIT DoFile
    CASE("TP.TABLE")
      Ok_TPpath= .TRUE.
    CASE("DEBUG")
      CALL LinToWrd(L,W,EoL)
      CALL WrdToInt(W,DebugLevel)
    CASE("ELEMENT")
      Ok_Database= .TRUE.
    CASE("DYNAMIC")
      Ok_Dynamic= .TRUE.
    CASE("PATH")
      Ok_Path=.TRUE.
    CASE("SOLUTION.MODEL","MIXTURE.MODEL")
      Ok_MixModel=.TRUE.
    CASE("SYSTEM","SYSTEM.AQUEOUS")
      Ok_Speciation=.TRUE.
      !System_Type="AQUEOUS"
    CASE("SYSTEM.GLOBAL")
      Ok_Speciation=.TRUE.
      !System_Type="GLOBAL" !not used yet
    CASE("SYSTEM.MOMAS")
      Ok_Speciation=.TRUE.
      !System_Type="MOMAS"  !for "virtual" system, cf momas benches
    CASE("SYSTEM.SIMPLEX")
      Ok_Simplex= .TRUE.
      !System_Type="SIMPLEX"

    ENDSELECT

  ENDDO DoFile

  CLOSE(F)
ENDSUBROUTINE Options_Read

ENDMODULE M_IO_Menu
