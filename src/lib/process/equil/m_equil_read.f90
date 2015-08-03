MODULE M_Equil_Read
  USE M_Trace,ONLY: iDebug,fTrc,T_,Stop_,Pause_
  USE M_Kinds
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Equil_Read_YesList
  PUBLIC:: Equil_Read_Debug
  PUBLIC:: Equil_Read_Conditions
  PUBLIC:: Equil_Read_Numeric
  PUBLIC:: Equil_Read_PhaseAdd
  
CONTAINS

SUBROUTINE Equil_Read_PhaseAdd(vFas)
!--
!-- read block SYSTEM.ROCK,
!-- -> read the non-fluid part of the system,
!-- -> update vFas%Mole
!--
  USE M_IOTools
  USE M_Files,  ONLY: NamFInn
  USE M_T_Phase,ONLY: T_Phase,Phase_Index
  !
  TYPE(T_Phase),DIMENSION(:),INTENT(INOUT):: vFas
  !
  CHARACTER(LEN=512):: L,W,W1
  LOGICAL :: EoL
  INTEGER :: F,ios,I
  REAL(dp):: X
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Read_PhaseAdd"
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFInn))
  !
  DoFile: DO
    ! 
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)
      !
      CASE("ENDINPUT"); EXIT DoFile
      !
      CASE("SYSTEM.ROCK") !
        DoBlock: DO
          !
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoBlock !skip comment lines
          CALL AppendToEnd(L,W,EoL)
          !
          SELECT CASE(W)
          CASE("ENDINPUT")              ;  EXIT DoFile
          CASE("END","ENDSYSTEM.ROCK")  ;  EXIT DoBlock
          CASE("MOLE","GRAM")           ;  CALL LinToWrd(L,W1,EoL)
          CASE DEFAULT
            CALL Stop_(TRIM(W)//" <- unknown keyword in SYSTEM.ROCK !!!") 
          ENDSELECT
          !
          I= Phase_Index(W1,vFas)
          !
          IF(I==0) &
          & CALL Stop_(TRIM(W1)//" <-phase unknown (SYSTEM.ROCK)")
          !
          IF(vFas(I)%iSpc==0) &
          & CALL Stop_(TRIM(W1)//" <-not valid, USE ONLY pure non'aqu'phases (SYSTEM.ROCK)")
          !
          CALL LinToWrd(L,W1,EoL)
          CALL WrdToReal(W1,X)
          !
          IF(X<=Zero) CALL Stop_(TRIM(W1)//" <-not valid as species amount  (SYSTEM.ROCK)")
          IF(W=="GRAM") X= X /vFas(I)%WeitKg /1.0D3
          !
          vFas(I)%Mole= X
          !
        ENDDO DoBlock
      !ENDCASE("SYSTEM.ROCK")
      !
    ENDSELECT
    !
  ENDDO DoFile
  CLOSE(F)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Read_PhaseAdd"
  !
ENDSUBROUTINE Equil_Read_PhaseAdd

SUBROUTINE Equil_Read_YesList(nFas,vFas,vYesList)
!--
!-- reads lists of included or excluded phases -> update vYesList
!--
  USE M_T_Phase,ONLY: T_Phase
  !
  INTEGER,      INTENT(IN) :: nFas
  TYPE(T_Phase),INTENT(IN) :: vFas(:)
  LOGICAL,      INTENT(OUT):: vYesList(:)
  !
  LOGICAL:: vInclud(nFas),vExclud(nFas)
  !vInclude= included phases; vExclude= excluded phases
  INTEGER:: I
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Read_YesList"
  !
  !-------------------------------- read list of "equilibrium phases" --
  vInclud=.FALSE.
  vExclud=.FALSE.
  !
  CALL Equil_Read_EquilPhase(vFas,vInclud,vExclud) !,InOk,ExOk)
  !
  vYesList=.TRUE. !-> default= all phases from DATAbase taken into account
  IF(COUNT(vInclud)>0) vYesList(:)= vInclud(:) !-> ONLY included phases
  IF(COUNT(vExclud)>0) vYesList(:)= vYesList(:) .AND. (.NOT. vExclud(:))
  !
  !-------------------------------/ read list of "equilibrium phases" --
  !
  !! vYesList(1)=.FALSE. !-> ??? to exclude species H2O as phase ???
  !
  !------------------------------------------------------------ trace --
  IF(iDebug>1) THEN
    PRINT '(/,A)',"< -- List of Phases --"
    DO I=1,nFas
      IF(vYesList(I)) PRINT '(A)',vFas(I)%NamFs
    ENDDO
    PRINT '(A,/)',"</-- List of Phases --"
    IF(iDebug==4) CALL Pause_
  ENDIF
  !
  IF(iDebug>0) THEN
    DO I=1,nFas
      IF(vYesList(I)) WRITE(fTrc,'(A)') vFas(I)%NamFs
    ENDDO
  ENDIF
  !-----------------------------------------------------------/ trace --
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Read_YesList"
ENDSUBROUTINE Equil_Read_YesList
  !  
SUBROUTINE Equil_Read_Debug(DebFormula,DebNewt,DebJacob)
!--
!-- read block CONDITIONS from arxim.inn to retrive DEBUG options
!--
  USE M_Trace,ONLY: LWarning
  USE M_IoTools
  USE M_Files,ONLY: NamFInn
  !
  LOGICAL,INTENT(OUT):: DebFormula,DebNewt,DebJacob
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W
  LOGICAL           :: EoL,Ok
  INTEGER           :: f,N,ios
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Read_Debug"
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  Ok=.FALSE.
  !
  DoFile: DO 
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    
    SELECT CASE(TRIM(W))
    !
    CASE("ENDINPUT"); EXIT  DoFile 
    !
    CASE("CONDITIONS") !LOGICAL global variables
      Ok=.TRUE.
      !
      DoDebug: DO
        
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=='!') CYCLE DoDebug
        CALL AppendToEnd(L,W,EoL)
        
        SELECT CASE(TRIM(W))
          !
          CASE("ENDINPUT"); EXIT DoFile
          
          CASE("END","ENDCONDITIONS"); EXIT DoDebug
          
          CASE("DEBUG")
            CALL LinToWrd(L,W,EoL); CALL WrdToInt(W,N)
            IF(N>0 .AND. iDebug==0) iDebug=N
            !iDebug is 0 IF it has not been already assigned in command line
          !
        ENDSELECT
      
      ENDDO DoDebug
      !
    ENDSELECT
  ENDDO DoFile
  CLOSE(f)
  !
  LWarning= (iDebug>0)
  !
  DebFormula= (iDebug==4)
  DebJacob=   (iDebug==4)
  DebNewt=    (iDebug==4 .OR. iDebug==9)
  !
  IF(.NOT.Ok) THEN 
    IF(iDebug>0) WRITE(fTrc,'(A)') &
    & "Block CONDITIONS not Found -> using default values for DEBUG"
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Read_Debug"
ENDSUBROUTINE Equil_Read_Debug

SUBROUTINE Equil_Read_OutputOptions(OutDistrib)
!--
!-- provisional, may be used for other output options ...
!-- read block CONDITIONS from arxim.inn for OUTPUT options
!--
  USE M_IoTools
  USE M_Files,ONLY: NamFInn
  !
  LOGICAL,INTENT(OUT):: OutDistrib
  !
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W
  LOGICAL           :: EoL,Ok
  INTEGER           :: f,ios
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Read_OutputOptions"
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  Ok=.FALSE.
  !
  OutDistrib= .FALSE.
  !
  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(TRIM(W))
      !
      CASE("ENDINPUT")
        EXIT  DoFile 
      !
      CASE("CONDITIONS")
        Ok=.TRUE.
        !
        DoDebug: DO
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoDebug
          CALL AppendToEnd(L,W,EoL)
          SELECT CASE(TRIM(W))
            CASE("ENDINPUT"); EXIT DoFile
            CASE("END","ENDCONDITIONS"); EXIT DoDebug
            !
            CASE("OUTPUT.DISTRIB"); OutDistrib=.TRUE. 
            !
          ENDSELECT
        ENDDO DoDebug
        !
    ENDSELECT
    !
  ENDDO DoFile
  !
  CLOSE(f)
  !
  IF(.NOT.Ok) THEN 
    IF(iDebug>0) WRITE(fTrc,'(A)') &
    & "Block CONDITIONS not Found -> using default values for OUTPUT options"
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Read_OutputOptions"
  !
ENDSUBROUTINE Equil_Read_OutputOptions

SUBROUTINE Equil_Read_EquilPhase(vFas,vInclude,vExclude) !,InOk,ExOk)
  USE M_IOTools !, ONLY:dimV,LinToWrd,GetUnit
  USE M_Files,ONLY: NamFInn
  USE M_T_Phase,ONLY: T_Phase,Phase_Index
  !
  TYPE(T_Phase),INTENT(IN) :: vFas(:)
  LOGICAL,      INTENT(OUT):: vInclude(:)
  LOGICAL,      INTENT(OUT):: vExclude(:)
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL:: EoL
  INTEGER:: F,I,ios
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Read_EquilPhase"
  !
  vInclude=.FALSE. !InOk=.FALSE.; 
  vExclude=.FALSE. !ExOk=.FALSE.; 
  !
  !!vInclude(:)= vFas(:)%iSpc/=0
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFInn))
  DoFile: DO
  
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    SELECT CASE(TRIM(W))
      
      CASE("ENDINPUT"); EXIT DoFile
      
      CASE("EQUIL.INCLUDE","EQU.INCLUDE")
        
        DoInclude: DO
          
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoInclude !skip comment lines
          CALL AppendToEnd(L,W,EoL)
          SELECT CASE(TRIM(W))
            CASE("ENDINPUT"); EXIT DoFile
            CASE("END","ENDEQU.INCLUDE","ENDEQUIL.INCLUDE"); EXIT DoInclude
          ENDSELECT
          
          I= Phase_Index(W,vFas)
          
          IF(I>0) THEN
            IF(vFas(I)%iSpc/=0) THEN !restriction to PURE phases, provisionally
              vInclude(I)=.TRUE.
              IF(iDebug>2) PRINT      '(A,A1,G15.6)', vFas(I)%NamFs,T_,vFas(I)%Grt
              IF(iDebug>0) WRITE(fTrc,'(A,A1,G15.6)') vFas(I)%NamFs,T_,vFas(I)%Grt
            ELSE
              !!vInclude(I)=.TRUE.
              !!IF(iDebug>2) PRINT      '(A,A1,G15.6)', vFas(I)%NamFs,T_,vFas(I)%Grt
              !!IF(iDebug>0) WRITE(fTrc,'(A,A1,G15.6)') vFas(I)%NamFs,T_,vFas(I)%Grt
              IF(iDebug>2) PRINT      '(2A)', TRIM(W)," is not PURE -> not included"
              IF(iDebug>0) WRITE(fTrc,'(2A)') TRIM(W)," is not PURE -> not included"
            ENDIF
          ELSE
            IF(iDebug>2) PRINT '(3A)',"WARNING, ",TRIM(W)," is not in current database"
            IF(iDebug>0) WRITE(fTrc,'(3A)') "WARNING !!! ",TRIM(W)," is not in current database"
          ENDIF
        
        ENDDO DoInclude
        
        !!InOk= COUNT(vInclude)>0
      !ENDCASE
      
      CASE("EQUIL.EXCLUDE","EQU.EXCLUDE")
        
        DoExclude: DO
          
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoExclude !skip comment lines
          CALL AppendToEnd(L,W,EoL)
          SELECT CASE(TRIM(W))
            CASE("ENDINPUT"); EXIT DoFile
            CASE("END","ENDEQUIL.EXCLUDE","ENDEQU.EXCLUDE"); EXIT DoExclude
          ENDSELECT
          
          I= Phase_Index(W,vFas)
          IF(I>0) vExclude(I)=.TRUE.
          
          IF(I>0) THEN !======================================< trace ==
            IF(iDebug>0) &
            & WRITE(fTrc,'(2A,A1,G15.6)') "EXCLUDED, Grt=",vFas(I)%NamFs,T_,vFas(I)%Grt
            IF(iDebug>1) &
            & PRINT '(2A,1X,G15.6)',"EXCLUDED, Grt=",vFas(I)%NamFs,vFas(I)%Grt
          ENDIF !============================================</ trace ==
          
        ENDDO DoExclude
        
        !!ExOk= COUNT(vExclude)>0
      !ENDCASE
    ENDSELECT
    
  ENDDO DoFile
  CLOSE(F)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Read_EquilPhase"
  !
ENDSUBROUTINE Equil_Read_EquilPhase

SUBROUTINE Equil_Read_Conditions(initMolNu,bOH2)
!--
!-- READ block CONDITIONS from arxim.inn
!--
  USE M_IoTools
  USE M_Files,ONLY: NamFInn
  USE M_Equil_Vars,ONLY: LogForAqu,DirectSub
  !
  REAL(dp), INTENT(OUT)  :: initMolNu
  LOGICAL,  INTENT(OUT)  :: bOH2
  !
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W
  LOGICAL           :: EoL,Ok
  INTEGER           :: f,ios
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_ReadConditions"
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  !
  Ok=.FALSE.
  !
  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)
      !
      CASE("ENDINPUT"); EXIT DoFile 
      !
      CASE("CONDITIONS")
        Ok=.TRUE.
        DoCond: DO
          !
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoCond
          CALL AppendToEnd(L,W,EoL)
          
          SELECT CASE(W)
          
            CASE("ENDINPUT")             ;  EXIT DoFile
            CASE("END","ENDCONDITIONS")  ;  EXIT DoCond
            
            CASE("H2OBALANCE")
              bOH2=.TRUE.
            
            CASE("INITIAL")
              CALL LinToWrd(L,W,EoL)
              CALL WrdToReal(W,initMolNu)
            
            !CASE("MOLEFORAQU")
            !  LogForAqu= .FALSE.
            !CASE("DIRECT")
            !  DirectSub= .TRUE.
          
          ENDSELECT
          !
        ENDDO DoCond
      !ENDCASE("CONDITIONS") 
    ENDSELECT
    !
  ENDDO DoFile
  !
  CLOSE(f)
  !
  IF(.NOT.Ok) THEN 
    IF(iDebug>0)  WRITE(fTrc,'(A)') "Block CONDITIONS not Found ...!!!"
    IF(iDebug>0)  PRINT '(A)',"WARNING: Block CONDITIONS not Found, using default values"
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_ReadConditions"
  !
ENDSUBROUTINE Equil_Read_Conditions

SUBROUTINE Equil_Read_Numeric( &
& initMolNu, &
& bOH2,      &
& LogForAqu, &
& DirectSub, &
& cMethod,   &
& bFinDif,   &
& Error,     &
& ErrorMsg   &
& )
!--
!-- READ block EQU.NUMERIC
!--
  USE M_IoTools
  USE M_Files,ONLY: NamFInn
  !
  REAL(dp),    INTENT(OUT)  :: initMolNu
  LOGICAL,     INTENT(INOUT):: bOH2
  LOGICAL,     INTENT(INOUT):: LogForAqu,DirectSub
  CHARACTER(*),INTENT(OUT)  :: cMethod
  LOGICAL,     INTENT(INOUT):: bFinDif
  LOGICAL,     INTENT(OUT)  :: Error
  CHARACTER(*),INTENT(OUT)  :: ErrorMsg
  !
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W
  LOGICAL           :: EoL,Ok
  INTEGER           :: f,ios
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Equil_Read_Numeric"
  !
  Error= .FALSE.
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  !
  Ok=.FALSE.
  
  DoFile: DO 
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    SELECT CASE(W)
      !
      CASE("ENDINPUT"); EXIT DoFile 
      !
      CASE("EQU.NUMERIC")
        Ok=.TRUE.
        DoCond: DO
        
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoCond
          CALL AppendToEnd(L,W,EoL)
          
          SELECT CASE(W)
          
            CASE("ENDINPUT")                ;  EXIT DoFile
            CASE("END","ENDEQU.NUMERIC")  ;  EXIT DoCond
            
            CASE("H2OBALANCE")
              bOH2=.TRUE.
            
            CASE("INITIAL")
              CALL LinToWrd(L,W,EoL)
              CALL WrdToReal(W,initMolNu)
            
            CASE("LOGFORAQU")
              CALL LinToWrd(L,W,EoL)
              SELECT CASE(TRIM(W))
                CASE("TRUE", "T", "YES")  ;  LogForAqu= .TRUE.
                CASE("FALSE", "F", "NO")  ;  LogForAqu= .FALSE.
                CASE DEFAULT
                  Error= .TRUE.
                  ErrorMsg= TRIM(W)//" is invalid keyword in EQU.NUMERIC"
              ENDSELECT
            
            !! obsolete : replace by LOGFORAQU= YES/NO !!
            !~ CASE("AQUEOUS")
              !~ CALL LinToWrd(L,W,EoL)
              !~ SELECT CASE(TRIM(W))
                !~ CASE("MOLE")  ;  LogForAqu= .FALSE.
                !~ CASE("LOG")   ;  LogForAqu= .TRUE.
                !~ CASE DEFAULT
                  !~ Error= .TRUE.
                  !~ ErrorMsg= TRIM(W)//" is invalid keyword in EQU.NUMERIC"
              !~ ENDSELECT
              
            CASE("DIRECT")
              CALL LinToWrd(L,W,EoL)
              SELECT CASE(TRIM(W))
                CASE("TRUE", "T", "YES")  ;  DirectSub= .TRUE.
                CASE("FALSE", "F", "NO")  ;  DirectSub= .FALSE.
                CASE DEFAULT
                  Error= .TRUE.
                  ErrorMsg= TRIM(W)//" is invalid keyword for DIRECT in EQU.NUMERIC"
              ENDSELECT
            
            CASE("JACOBIAN")
              CALL LinToWrd(L,W,EoL)
              SELECT CASE(TRIM(W))
                CASE("NUMERIC")   ;  bFinDif= .TRUE.
                CASE("ANALYTIC")  ;  bFinDif= .FALSE.
                CASE DEFAULT
                  Error= .TRUE.
                  ErrorMsg= TRIM(W)//" is invalid keyword in EQU.NUMERIC"
              ENDSELECT
            
            CASE("METHOD")
              CALL LinToWrd(L,W,EoL)
              SELECT CASE(TRIM(W))
                CASE("NEWTON")       ;   cMethod= TRIM(W)
                CASE("NEWTLNSRCH")   ;   cMethod= TRIM(W)
                CASE("NEWTONPRESS")  ;   cMethod= TRIM(W)
                CASE("BROYDEN")      ;   cMethod= TRIM(W)
                CASE("NEWTONCHESS")  ;   cMethod= TRIM(W)
                CASE("TENSOLVE_1")   ;   cMethod= TRIM(W)
                CASE("NEWTONKELLEY") ;   cMethod= TRIM(W)
                CASE("NEWTONWALKER") ;   cMethod= TRIM(W)
                CASE DEFAULT
                  Error= .TRUE.
                  ErrorMsg= TRIM(W)//" is invalid METHOD in EQU.NUMERIC"
              END SELECT
            !_
            CASE DEFAULT
              Error= .TRUE.
              ErrorMsg= TRIM(W)//" is invalid keyword in EQU.NUMERIC"
              !PRINT *,TRIM(W)//" = unknown keyword in EQU.NUMERIC !!"
            
          ENDSELECT
        
        ENDDO DoCond
      !ENDCASE("CONDITIONS") 
    ENDSELECT
    
  ENDDO DoFile
  CLOSE(f)
  !
  IF(.NOT.Ok) THEN 
    IF(iDebug>0)  WRITE(fTrc,'(A)') "Block EQU.NUMERIC not Found ...!!!"
    !IF(iDebug>0)  PRINT '(A)',"WARNING: Block CONDITIONS not Found, using default values"
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Equil_Read_Numeric"
  !
ENDSUBROUTINE Equil_Read_Numeric

ENDMODULE M_Equil_Read

