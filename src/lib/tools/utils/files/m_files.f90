MODULE M_Files

  USE M_Files_Vars
  USE M_Files_Index
  USE M_FileUtils
  IMPLICIT NONE
  !
  PRIVATE
  !
  !// Private variables
  CHARACTER(LEN=80) :: NamFIn0

  !// Exported M_Files_Vars Variables
  PUBLIC :: NamFInn
  PUBLIC :: cTitle
  PUBLIC :: NamFLogK,NamFEle,NamFKin,NamFSol,NamFPtz
  PUBLIC :: NamDtbAqu,NamDtbMin,NamDtbMlt
  PUBLIC :: DirOut,DirLog
  PUBLIC :: DirDtbOut,DirDtbLog 

  !// Public Methods
  PUBLIC :: Files_Vars_Init
  PUBLIC :: Files_BuildInput
  PUBLIC :: Files_PathRead
  PUBLIC :: Files_PathRead_From_FileName
  PUBLIC :: Info_PathRead
  
  ! Exported M_Files_Index methods
  PUBLIC :: Files_Index_Write
  PUBLIC :: Files_Index_Open
  PUBLIC :: Files_Index_Close

  ! Exported M_Files_Utils methods
  PUBLIC :: File_Exist
  PUBLIC :: File_Write_Date_Time 
  PUBLIC :: File_Path

  !// Private Methods
  PRIVATE :: Files_PathRead_From_CommandArgs
  PRIVATE :: Files_PathRead_From_UserInput
  
CONTAINS

  SUBROUTINE Info_PathRead
    !===================================================================
    ! Print the Input FileName 
    !===================================================================
    USE M_Trace,  ONLY: fHtm
    IMPLICIT NONE
    !--
    write(*,'(A)') "Arxim Script File : "// TRIM(NamFIn0)
    call Files_Index_InputFile(fHtm, NamFIn0)

  END SUBROUTINE Info_PathRead

  !---

  SUBROUTINE Files_PathRead(Ok)
    !===================================================================
    ! Read the Input FileName ( Default = "arxim.inn" )
    !===================================================================
    USE M_Trace, ONLY: Stop_
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: Ok
    !--
    Ok = .false.
    IF (.NOT.Ok) CALL Files_PathRead_From_CommandArgs(Ok)
    IF (.NOT.Ok) CALL Files_PathRead_From_FileName("arxim.inn",Ok)
    IF (.NOT.Ok) CALL Files_PathRead_From_UserInput(Ok)
    IF (.NOT.Ok) CALL Stop_("Error Reading Input File")

    CALL Info_PathRead
    
  END SUBROUTINE Files_PathRead

  !---

  SUBROUTINE Files_BuildInput
    !===================================================================
    ! Build the input data file from main file + include
    ! Init Files Variables
    !===================================================================
    IMPLICIT NONE
    !---
    ! Merge include files
    CALL Files_BuildInput_IncludeFiles

    ! Init files names variable
    CALL Files_Vars_Init(NamFInn)

  END SUBROUTINE Files_BuildInput

  !---

  SUBROUTINE Files_PathRead_From_FileName(sFileName, Ok) 
    !===================================================================
    ! Read the Input FileName from the command line argumants
    ! -> Default Extension = .inn
    ! Output in NamFIn0 
    !===================================================================
    IMPLICIT NONE
    LOGICAL,INTENT(OUT):: Ok
    CHARACTER(LEN=*), INTENT(IN):: sFileName
    !--
    Ok=.FALSE.
    NamFIn0="NOFILE"
    
    IF(.NOT. File_Exist(TRIM(sFileName))) THEN
       Ok=.FALSE.
       NamFIn0="NOFILE"
       PRINT '(A)',TRIM(sFileName)//"-> File Not Found ...!!!"
    ELSE
       Ok=.TRUE.
       NamFIn0=TRIM(sFileName)
    END IF
    
  END SUBROUTINE Files_PathRead_From_FileName

  !---

  SUBROUTINE Files_PathRead_From_CommandArgs(Ok) 
    !===================================================================
    ! Read the Input FileName from the command line argumants
    ! -> Default Extension = .inn
    ! Output in NamFIn0 
    !===================================================================
    USE M_FArgC, ONLY: F_GETARG, F_IARGC
    USE M_VarStr,ONLY: T_VarStr,VarStr_Get,VarStr_Char
    IMPLICIT NONE
    !--
    LOGICAL,INTENT(OUT):: Ok
    !--
    CHARACTER(LEN=255):: Str
    !--
    Ok=.FALSE.
    NamFIn0="NOFILE"

    IF(F_IARGC()>0) THEN
    
      CALL F_GETARG(1,Str)
      
      IF(INDEX(Str,".")<1) Str=TRIM(Str)//".inn"
      ! Default file name suffix is .inn
      
      IF(.NOT. File_Exist(TRIM(Str))) THEN
        Ok=.FALSE.
        PRINT '(A)',TRIM(Str)//"-> File Not Found ...!!!"
      ELSE
        Ok=.TRUE.
        NamFIn0=TRIM(Str)
      END IF
      
    END IF
    
    RETURN
  END SUBROUTINE Files_PathRead_From_CommandArgs

  !--

  SUBROUTINE Files_PathRead_From_UserInput(Ok) 
    !===================================================================
    ! Read the Input FileName by using interactive questions
    ! -> Default FileName  = arxim.inn
    ! -> Default Extension = .inn
    ! Output in NamFIn0 
    !===================================================================
    USE M_VarStr,ONLY: T_VarStr,VarStr_Get,VarStr_Char
    IMPLICIT NONE
    !--
    LOGICAL,INTENT(OUT):: Ok
    !--
    TYPE(T_VarStr):: VStr
    CHARACTER(LEN=255):: Str
    !--
    Ok=.FALSE.
    NamFIn0="NOFILE"
    
    DoNam: DO
       WRITE(*,'(A)') "Input File Name [ Return = arxim.inn ] ?"
       CALL VarStr_Get(String=VStr)
       IF(VarStr_Char(VStr)=="") THEN
          Str="arxim.inn"
       ELSE
          Str=VarStr_Char(VStr)
          ! Default file name suffix is .inn 
          IF(INDEX(Str,".")<1) Str=TRIM(Str)//".inn"        
       ENDIF
       IF(File_Exist(Str)) THEN
          Ok=.TRUE.
          EXIT DoNam
       ELSE                   
          Ok=.FALSE.
          PRINT '(A)',"File Not Found ...!!!"
       ENDIF
    ENDDO DoNam
    NamFIn0=TRIM(Str)
    
  END SUBROUTINE Files_PathRead_From_UserInput

  !---

  SUBROUTINE Files_BuildInput_IncludeFiles
    !==========================================================
    ! Build the Input file NamFinn = "arxim_inn.tmp"
    !
    !  -> Merge nested include files 
    !  -> Skip comments
    ! 
    ! Input  File  = NamFIn0 ( + INCLUDE Files )
    ! Output File  = NamFinn = "arxim_inn.tmp" 
    ! 
    !==========================================================
    USE M_Trace,  ONLY: iDebug,Pause_,fHtm,Stop_
    USE M_IOTools,ONLY: GetUnit,LinToWrd,AppendToEnd
    !
    INTEGER:: fAll,fInn,fAdd,fAdd2
    CHARACTER(LEN=512):: L, LL, W, W1, W2
    CHARACTER(LEN=512):: sIncludeFile1, sIncludeFile2
    LOGICAL:: EoL
    INTEGER:: ios     
    !
    NamFInn="arxim_inn.tmp"
    !
    CALL GetUnit(fAll); OPEN(fAll,FILE=TRIM(NamFInn))
    CALL GetUnit(fInn); OPEN(fInn,FILE=TRIM(NamFIn0))

    !=========== READ ROOT FILE  ================================
    Do1: DO
      !
      READ(fInn,'(A)',IOSTAT=ios) LL
      IF(ios/=0) EXIT Do1 
      L=TRIM(LL)
      CALL LinToWrd(L,W,EoL)
      CALL AppendToEnd(L,W,EoL)
      !
      IF(W(1:1)=='!') CYCLE Do1 ! skip comment lines
      IF(TRIM(W)=='/*') THEN !skip comment block
        DO
          READ(fInn,'(A)',IOSTAT=ios) L
          IF(ios/=0) THEN
            PRINT '(A)', "!!!WARNING!!! COMMENT BLOCK UNTERMINATED !!!WARNING!!!"
            EXIT Do1
          ENDIF
          CALL LinToWrd(L,W,EoL)
          IF(TRIM(W)=='*/') EXIT
        ENDDO
        CYCLE Do1
      ENDIF
      !
      SELECT CASE(TRIM(W))
      !
      CASE DEFAULT
        WRITE(fAll,'(A)') TRIM(LL)
      !
      CASE("ENDINPUT")
        EXIT Do1
        !
      CASE("INCLUDE") 
        ! Get Include File Name Level 1        
        CALL LinToWrd(L,W,EoL,"NO") 
        sIncludeFile1 = W
        
        IF (File_Exist_System(sIncludeFile1)) THEN 
          
          WRITE(*,*) "INCLUDE LEVEL 1 : ", TRIM(sIncludeFile1)
          
          ! Open Include File Level 1
          CALL GetUnit(fAdd)
          OPEN(fAdd,FILE=TRIM(sIncludeFile1))
          CALL Files_Index_Include(fHtm,sIncludeFile1)
          WRITE(fAll,'(A,A)') "!!!!!!!!!!!!Insert File_begin!!!!!!!!!!!!"
          WRITE(fAll,'(A,A)') "! INCLUDE ", TRIM(sIncludeFile1)
          
          !=========== READ FILE INCLUDE LEVEL 1 =======================
          Do2: DO
            !
            READ(fAdd,'(A)',IOSTAT=ios) LL
            IF(ios/=0) EXIT Do2
            L=LL
            CALL LinToWrd(L,W1,EoL)
            CALL AppendToEnd(L,W1,EoL)
            !
            IF(W1(1:1)=='!') CYCLE Do2 !skip comment lines
            IF(TRIM(W1)=='/*') THEN !skip comment block
              DO
                READ(fInn,'(A)',IOSTAT=ios) L
                IF(ios/=0) THEN
                  PRINT '(A)', &
                  & "!!!WARNING!!! COMMENT BLOCK UNTERMINATED !!!WARNING!!!"
                  EXIT Do2
                ENDIF
                CALL LinToWrd(L,W1,EoL)
                IF(TRIM(W1)=='*/') EXIT
              ENDDO
              CYCLE Do2
            ENDIF
            !
            IF(TRIM(W1)=="ENDINPUT")  EXIT Do2
            !
            IF(TRIM(W1)=="INCLUDE") THEN 
              
              ! Get Include File Name Level 2
              CALL LinToWrd(L,W1,EoL,"NO") 
              sIncludeFile2 = W1
              
              IF(File_Exist_System(sIncludeFile2)) THEN 
                
                WRITE(*,*) "INCLUDE LEVEL 2 : .. ", TRIM(sIncludeFile2)
                
                ! Open Include File Level 2
                CALL GetUnit(fAdd2)
                OPEN(fAdd2,FILE=TRIM(sIncludeFile2))
                CALL Files_Index_Include(fHtm,sIncludeFile2)
                WRITE(fAll,'(A)') "!!!!!!!!!!!!Insert File_begin!!!!!!!!!!!!"
                WRITE(fAll,'(A,A)') "! INCLUDE ", TRIM(sIncludeFile2)
                
                !=========== READ FILE LEVEL 2 =========================
                do3: DO 
                  READ(fAdd2,'(A)',IOSTAT=ios) LL
                  IF(ios/=0) EXIT do3
                  L=TRIM(LL)
                  CALL LinToWrd(L,W2,EoL)
                  CALL AppendToEnd(L,W2,EoL)
                  !
                  IF(W2(1:1)=='!') CYCLE do3 !skip comment lines
                  IF(TRIM(W2)=='/*') THEN !skip comment block
                    DO
                      READ(fInn,'(A)',IOSTAT=ios) L
                      IF(ios/=0) THEN
                        PRINT '(A)', &
                        & "!!!WARNING!!! COMMENT BLOCK UNTERMINATED !!!WARNING!!!"
                        EXIT Do3
                      ENDIF
                      CALL LinToWrd(L,W2,EoL)
                      IF(TRIM(W2)=='*/') EXIT
                    ENDDO
                    CYCLE do3
                  ENDIF
                  !
                  IF(TRIM(W2)=="ENDINPUT") EXIT do3
                  !
                  IF(TRIM(W2)=="INCLUDE") THEN
                    PRINT '(A)', "!!!WARNING!!! TOO Many Nested INCLUDEs !!!WARNING!!!"
                    CALL LinToWrd(L,W2,EoL,"NO") 
                    CALL Stop_ ("INCLUDE LEVEL 3 :"//TRIM(W2)//" NOT IMPLEMENTED")
                  END IF
                  !
                  WRITE(fAll,'(A)') TRIM(LL)
                ENDDO do3
                
                ! Close Include File Level 2
                WRITE(fAll,'(A)') "!!!!!!!!!!!!Insert File _end!!!!!!!!!!!!!"
                CLOSE(fAdd2)
              
              ELSE
                
                CALL Stop_("File "//TRIM(sIncludeFile2)//" NOT FOUND -> Can Not INCLUDE !!!")
              
              ENDIF
              !
            ENDIF ! nested INCLUDE
            !
            WRITE(fAll,'(A)') TRIM(LL)
          ENDDO Do2
          ! Close Include File Level 1
          WRITE(fAll,'(A)') "!!!!!!!!!!!!Insert File _end!!!!!!!!!!!!!"
          CLOSE(fAdd)
        ELSE
          CALL Stop_("FILE "//TRIM(sIncludeFile1)//" NOT FOUND -> Can Not INCLUDE !!!")
        ENDIF
      ENDSELECT
    ENDDO Do1
    ! Close Root File
    CLOSE(fInn)
    ! Close Result Merged File
    !! WRITE(fAll,'(A)') "ENDINPUT" !!! IF(Mod="APPEND") 
    CLOSE(fAll)

  END SUBROUTINE Files_BuildInput_IncludeFiles

  !---

  SUBROUTINE Files_Vars_Init(sFilNam)
    !==========================================================
    ! Init the Directories and Input Files Names
    !==========================================================
    IMPLICIT NONE
    CHARACTER(*),INTENT(IN):: sFilNam
    !---
    cTitle=""

    ! defaut directories
    DirOut="out_"
    DirLog="log_"

    DirDtbOut="dtb_out_"
    DirDtbLog="dtb_log_"

    ! make the input file the default file for all data
    NamFLogK=TRIM(sFilNam) ! file for logK
    NamFEle= TRIM(sFilNam) ! file for elements
    NamFSol= TRIM(sFilNam) ! file for solutions
    NamFKin= TRIM(sFilNam) ! file for kinetics
    NamFPtz= TRIM(sFilNam) ! file for pitzer parameters

  END SUBROUTINE Files_Vars_Init

END MODULE M_Files

