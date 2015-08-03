MODULE M_Dynam_Read
  USE M_Kinds
  USE M_Trace
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC:: Dynam_Read
  PUBLIC:: Dynam_ReadColumn

CONTAINS

SUBROUTINE Dynam_Read(DynTime,DynBox,DynBoxUser,Ok)
  USE M_Numeric_Const,ONLY: LN10
  USE M_KinModel_Read,ONLY: KinModel_Init
  USE M_Dynam_Tools,  ONLY: Dynam_TimeFactor
  !
  USE M_Global_Vars,  ONLY: vSpc,vKinModel
  USE M_Dynam_Vars,   ONLY: T_DynTime,T_DynBox,TimeFactor,TUnit
  USE M_Dynam_Vars,   ONLY: sModelSurf,Implicit_Surface,VBox0
  USE M_Dynam_Vars,   ONLY: LogForAqu,DirectSub,bFinDif
  USE M_Dynam_Vars,   ONLY: AdjustVolBox
  !
  TYPE(T_DynTime),INTENT(INOUT):: DynTime
  TYPE(T_DynBox), INTENT(INOUT):: DynBox
  TYPE(T_DynBox), INTENT(OUT)  :: DynBoxUser
  LOGICAL,        INTENT(OUT)  :: Ok
  !
  LOGICAL:: ReadOk
  CHARACTER(LEN=255):: Msg
  REAL(dp):: A
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dynam_Read"
  !
  Ok=.TRUE.
  !
  !-----------------------------------------------read time parameters--
  CALL Dynam_ReadTime(DynTime,ReadOk)
  !-> TUnit,TFinal,dTime,dTMin,dTMax,dTSav,TAgain
  IF (.NOT.ReadOk) THEN
     CALL Stop_("Block DYNAMIC.TIME not found")
     Ok = .FALSE.
     RETURN
  END IF
  DynTime%TimeFactor= Dynam_TimeFactor(DynTime%TUnit)
  !--------------------------------------------------------------------/
  !
  !------------------------------------------------read box parameters--
  CALL Dynam_ReadBox(DynBox,ReadOk,Msg)
  !-> VBox,dX,etc.
  IF(.NOT. ReadOk .AND. iDebug>2) PRINT '(A)',TRIM(Msg)
  !
  DynBoxUser= DynBox
  !
  IF(AdjustVolBox .AND. DynBox%VBox /= VBox0) THEN
    A= (DynBox%VBox/VBox0)**(1.0D0/3.0D0)
    DynBox%VBox=   VBox0
    DynBox%DX=     DynBox%DX     /A
    DynBox%UDarcy= DynBox%UDarcy /A
  ENDIF
  !
  IF(iDebug>2) THEN
    WRITE(fTrc,'(/,A)')     "<  BOX VOL/DX/UDARCY ADJUSTED==>"
    WRITE(fTrc,'(A,2G15.6)') "  VBox  =", DynBoxUser%VBox,   DynBox%VBox
    WRITE(fTrc,'(A,2G15.6)') "  dX    =", DynBoxUser%dX,     DynBox%dX
    WRITE(fTrc,'(A,2G15.6)') "  UDarcy=", DynBoxUser%UDarcy, DynBox%UDarcy
    WRITE(fTrc,'(A,/)')     "</ BOX VOL/DX/UDARCY ADJUSTED==>"
  ENDIF
  !---------------------------------------------------------------------
  !
  !------------------------------------------ read numeric parameters --
  CALL Dynam_ReadNumeric(ReadOk)
  !
  !IF(.NOT. LogForAqu) DirectSub= .TRUE.
  !IF(.NOT. LogForAqu) bFinDif=   .TRUE.
  IF(DirectSub) bFinDif=   .TRUE.
  !
  IF(sModelSurf=="CRUNCH") Implicit_Surface=.FALSE.
  !---------------------------------------------------------------------
  !
  !---------------------------------------------- read kinetic models --
  CALL KinModel_Init(vSpc)
  !! ReadOk= SIZE(vKinModel)>0
  !! IF(.NOT. ReadOk) CALL Stop_("NO Kinetic Data Found") !______STOP
  IF(SIZE(vKinModel)>0 .AND. iDebug>2) CALL Check_KinRate_CalcQsK
  !---------------------------------------------/ read kinetic models --
  !
  !------------------------------------- read box mineral composition --
  CALL Dynam_ReadRock(ReadOk)
  !! IF(.NOT.ReadOk)  CALL Stop_("NO Kinetic Minerals Found")
  !---------------------------------------------------------------------
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dynam_Read"
  !
  RETURN
ENDSUBROUTINE Dynam_Read

SUBROUTINE Dynam_ReadTime(T,Ok)
!--
!-- scan the input file, read block DYNAMIC.TIME
!-- -> temporal conditions of the run
!-- TYPE:: T_DynTime
!--   CHARACTER(LEN-6):: TUnit !-- "DAY","YEAR","HOUR"
!--   REAL(dp)::  &
!--   & Time,     & !-- current time
!--   & dTime,    & !-- current time step
!--   & dTmin,    & !-- minimal time step
!--   & dTMax,    & !-- maximal time step
!--   & dTSav       !-- time laps between two records on x_time.tab
!-- END TYPE T_DynTime
!--
  USE M_Files,ONLY: NamFInn
  USE M_IOTools
  USE M_Dynam_Vars, ONLY: SteadyState_Stop,TimeIsSeconds
  USE M_Dynam_Vars, ONLY: T_DynTime,TUnit,TFinal,dTime,dTMin,dTMax,dTSav,TimeFactor
  USE M_Dynam_Tools,ONLY: Dynam_TimeFactor
  !
  TYPE(T_DynTime),INTENT(INOUT):: T
  LOGICAL,        INTENT(OUT)  :: Ok
  !
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W,W1,W2,WW
  LOGICAL           :: EoL
  INTEGER           :: F,ios
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dynam_ReadTime"
  !
  CALL GetUnit(F)
  CALL OPENFILE(F,FILE=TRIM(NamFInn))
  Ok=.FALSE.
  !
  Do01: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT Do01
    CALL LinToWrd(L,W,EoL)
    !
    IF(W(1:1)=='!')   CYCLE Do01 !skip comment lines
    IF(W=="END" .AND. .NOT. Eol ) THEN !when END followed by keyword, appEND
      CALL LinToWrd(L,WW,EoL); W=TRIM(W)//TRIM(WW)
    ENDIF
    !
    SELECT CASE(TRIM(W))
    !
    CASE("ENDINPUT"); EXIT  Do01
    !
    CASE("DYNAMIC.TIME","DYNAMIC")
      !
      Ok=.TRUE.
      !
      Do02: DO
        !
        READ(F,'(A)',IOSTAT=ios) L  ;  IF(ios/=0) EXIT Do01
        CALL LinToWrd(L,W1,EoL)     ;  CALL AppendToEnd(L,W1,EoL)
        IF(W1(1:1)=='!') CYCLE Do02
        CALL LinToWrd(L,W2,EoL)
        !-> W2 is second word on line
        !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
        !
        SELECT CASE(TRIM(W1))
        !
        !----------------------------------------------timing parameters
        CASE("ENDINPUT")                                    ;  EXIT Do01
        CASE("END","ENDDYNAMIC.TIME","ENDDYNAMIC")          ;  EXIT Do02
        !
        CASE("STEADY")
          SteadyState_Stop= (TRIM(W2)=="STOP")
        CASE("TUNIT")
          SELECT CASE(TRIM(W2))
          CASE("SECOND")  ;  T%TUnit="SECOND"
          CASE("MINUTE")  ;  T%TUnit="MINUTE"
          CASE("HOUR")    ;  T%TUnit="HOUR"
          CASE("DAY" )    ;  T%TUnit="DAY"
          CASE("YEAR")    ;  T%TUnit="YEAR"
          CASE DEFAULT
            Ok= .FALSE.
            CALL Stop_ &
            & (TRIM(W1)//" <-unknown TUNIT in DYNAMIC.TIME / DYNAMIC")
          END SELECT
        !
        CASE("TFIN")   ;  CALL WrdToReal(W2,T%TFinal) !duration of simulation!
        CASE("TFINAL") ;  CALL WrdToReal(W2,T%TFinal) !duration of simulation!
        CASE("DTIME")  ;  CALL WrdToReal(W2,T%dTime)  !initial time step     !
        CASE("DTMIN")  ;  CALL WrdToReal(W2,T%dTMin)  !minimal time step     !
        CASE("DTMAX")  ;  CALL WrdToReal(W2,T%dTMax)  !maximal time step     !
        CASE("DTSAV")  ;  CALL WrdToReal(W2,T%dTSav)  !edition time step     !
        !! CASE DEFAULT
        !!   Ok= .FALSE.
        !!   CALL Stop_ &
        !!   & (TRIM(W1)//" <-unknown KeyWord in DYNAMIC.TIME or DYNAMIC")
        !
        !---------------------------------------------/timing parameters
        END SELECT
        !
      END DO Do02
    !_CASE("DYNAMIC")
    END SELECT
  END DO Do01
  CALL CLOSEFILE(F)
  !
  IF(T%dTSav==Zero) T%dTSav= T%TFinal/2.E3
  !
  !------------------------ if user given max time step is too small ---
  IF(T%dTMax>T%TFinal/100._dp) T%dTMax= T%TFinal/100._dp
  !-------------------- max time step set to 1/100 of total duration ---
  IF(T%dTMax==Zero) T%dTMax= T%TFinal/100._dp
  !
  T%TimeFactor= Dynam_TimeFactor(T%TUnit)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dynam_ReadTime"
  !
  RETURN
ENDSUBROUTINE Dynam_ReadTime

SUBROUTINE Dynam_ReadRock(Ok)
!--
!-- read block DYNAMIC.ROCK
!--
  USE M_Global_Vars, ONLY: vSpc,vFas,vKinFas,vKinModel
  USE M_KinFas_Read, ONLY: T_LnkKin,KinFas_BuildLnk,KinFas_LnkToVec
  USE M_Dynam_Vars,  ONLY: sModelSurf
  !
  LOGICAL,INTENT(OUT):: Ok
  !
  TYPE(T_LnkKin),POINTER:: LnkKin
  REAL(dp):: x
  INTEGER :: I,N
  !
  !------------------------------------------------------read ROCK block
  !
  CALL KinFas_BuildLnk( &  !
  & vFas,vKinModel,sModelSurf, & !IN
  & N,LnkKin)                    !OUT
  !
  IF(N>0) THEN !mod 12/06/2008 18:15
    Ok=.TRUE.
    !
    IF(ALLOCATED(vKinFas)) DEALLOCATE(vKinFas); ALLOCATE(vKinFas(N))
    !
    CALL KinFas_LnkToVec(LnkKin,vFas,vKinFas) !IN,OUT
    !
  ELSE
    Ok=.FALSE.
    RETURN !------------------------------------------------------RETURN
  ENDIF
  !
  IF(iDebug>2) THEN
    !
    PRINT '(/,A)',"Checking input from DYNAMIC.ROCK:"
    !
    IF(SIZE(vKinFas)>0) THEN
      
      PRINT '(5A16)', &
      & "Kinetic_Phase___", &
      & "Thermo_Model____", &
      & "Kinetic_Model___", &
      & "__-___Radius/m__", &
      & "________Surf/g__"
      
      DO I=1,SIZE(vKinFas)
        ! IF(vKinFas(I)%iFas<1) &
        ! & CALL Stop_("problem with "//TRIM(vKinFas(I)%Name))
        IF(vKinFas(I)%iKin>0) THEN
          PRINT '(3(A15,1X),2(G15.3,1X))', &
          & vKinFas(I)%NamKF, &
          & vFas(vKinFas(I)%iFas)%NamFs, &
          & vKinModel(vKinFas(I)%iKin)%Name, &
          & vKinFas(I)%Dat%Radius, &
          & vKinFas(I)%Dat%SurfKg
        ELSE
          PRINT '(3(A15,1X),2(G15.6,1X))', &
          & vKinFas(I)%NamKF, &
          & vFas(vKinFas(I)%iFas)%NamFs, &
          & "EQUILIBRIUM", &
          & vKinFas(I)%Dat%Radius, &
          & vKinFas(I)%Dat%SurfKg/1.0D3
        ENDIF
      ENDDO
    
    ELSE
      
      PRINT *,"EMPTY BOX ... !!"
    
    ENDIF
    !
    CALL Pause_
  ENDIF
  !
  !-----------------------------------------------------/read ROCK block
  !
  !------------------------------- normalize PhiM to SUM(VolFract)-1 ---
  x= SUM(vKinFas(:)%Dat%PhiM) !, MASK=vKinFas(:)%Dat%cSat /= "MINIMAL")
  IF(x>Zero) THEN
    DO I=1,SIZE(vKinFas)
      IF(vKinFas(I)%Dat%cSat /= "M") THEN           ! M(INIMAL
        vKinFas(I)%Dat%PhiM= vKinFas(I)%Dat%PhiM /x
      ELSE
        vKinFas(I)%Dat%PhiM= Zero
      ENDIF
    ENDDO
  ELSE
    CALL Stop_("SUM(vKinFas(1:nMk)%PhiM <0 ???") !------------------STOP
  ENDIF
  !---------------------------------------------------------------------
  !
  RETURN
ENDSUBROUTINE Dynam_ReadRock

SUBROUTINE Dynam_ReadBox(Box,Ok,Msg)
!--
!-- read block DYNAMIC.BOX -> spatial conditions
!--
  USE M_Files,ONLY: NamFInn
  USE M_IOTools
  USE M_Dynam_Vars,ONLY: T_DynBox,sModelSurf
  !
  TYPE(T_DynBox),INTENT(INOUT):: Box
  LOGICAL,       INTENT(OUT)  :: Ok
  CHARACTER(*),  INTENT(OUT)  :: Msg
  !
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W,W1,W2
  LOGICAL           :: EoL
  INTEGER           :: F,ios
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dynam_ReadBox"
  !
  CALL GetUnit(F)
  CALL OPENFILE(F,FILE=TRIM(NamFInn))
  Ok=.FALSE.
  Msg= "Ok"
  !
  Do01: DO
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT Do01
    CALL LinToWrd(L,W,EoL); CALL AppendToEnd(L,W,EoL)
    IF(W(1:1)=='!')   CYCLE Do01 !skip comment lines
    !
    SELECT CASE(TRIM(W))
    !
    CASE("ENDINPUT"); EXIT  Do01
    !
    CASE("DYNAMIC.BOX","DYNAMIC")
      Ok=.TRUE.
      Do02: DO
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT Do01
        CALL LinToWrd(L,W1,EoL); CALL AppENDToEnd(L,W1,EoL)
        IF(W1(1:1)=='!') CYCLE Do02
        CALL LinToWrd(L,W2,EoL)
        !-> W2 is second word on line
        !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
        !
        SELECT CASE(TRIM(W1))
        !
        CASE("ENDINPUT"); EXIT Do01
        CASE("END","ENDDYNAMIC.BOX","ENDDYNAMIC"); EXIT Do02
        !
        CASE("MODEL","SURFACE")
          SELECT CASE(TRIM(W2))
          CASE("SPHERE") ; sModelSurf= "SPHERE" ! (= the default model)
          CASE("CRUNCH") ; sModelSurf= "CRUNCH"
          CASE DEFAULT
            Ok= .FALSE.
            Msg="In DYNAMIC / DYNAMIC.BOX, "//TRIM(W2)//" is Unknown Keyword"
          END SELECT
          !
        CASE("VOLUME")
          SELECT CASE(TRIM(W2))
            CASE("FREE");  Box%VFixed= .FALSE. !UpdateMassFluid=.FALSE.
            CASE("FIXED"); Box%VFixed= .TRUE.  !UpdateMassFluid=.TRUE.
            CASE DEFAULT
              Ok= .FALSE.
              Msg= "In DYNAMIC / DYNAMIC.BOX, VOLUME is either FREE or FIXED"
          END SELECT
        !---------------------------------------------------------------------!
        !--------------------------------------------------spatial parameters.!
        CASE("DX");       CALL WrdToReal(W2,Box%dX)     !length of box........!
        CASE("VOLBOX");   CALL WrdToReal(W2,Box%VBox)   !volume of box........!
        CASE("UDARCY");   CALL WrdToReal(W2,Box%UDarcy) !flux rate,length/time!
        CASE("POROSITY"); CALL WrdToReal(W2,Box%PhiF)   !initial porosity.....!
        CASE("CELLS");    CALL WrdToInt (W2,Box%nCell)  !nr of boxes..........!
        !---------------------------------------------------------------------!
        !
        CASE DEFAULT
          Ok= .FALSE.
          Msg= "In DYNAMIC / DYNAMIC.BOX, "//TRIM(W1)//" is unknown KeyWord !!"
        !
        END SELECT
      ENDDO Do02
    !_CASE("DYNAMIC")
    END SELECT
  ENDDO Do01
  CALL CLOSEFILE(F)
  !
  IF(.NOT. Ok) Msg= "Block DYNAMIC.BOX not found, using default values !!!"
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dynam_ReadBox"
  !
ENDSUBROUTINE Dynam_ReadBox

SUBROUTINE Dynam_ReadColumn(Column,Ok,Msg)
!--
!-- read block DYNAMIC.1D
!--
  USE M_Files,ONLY: NamFInn
  USE M_IOTools
  USE M_Dynam_Vars,ONLY: T_DynColumn
  !
  TYPE(T_DynColumn),INTENT(INOUT):: Column
  LOGICAL,          INTENT(OUT)  :: Ok
  CHARACTER(*),     INTENT(OUT)  :: Msg
  !
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=80) :: W,W1,W2
  LOGICAL           :: EoL
  INTEGER           :: F,ios
  
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dynam_ReadColumn"
  !
  CALL GetUnit(F)
  CALL OPENFILE(F,FILE=TRIM(NamFInn))
  !
  Ok=.FALSE.
  Msg= "Ok"
  !
  Do01: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT Do01
    CALL LinToWrd(L,W,EoL); CALL AppendToEnd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE Do01 !skip comment lines
    !
    SELECT CASE(TRIM(W))
    !
    CASE("ENDINPUT"); EXIT  Do01
    !
    CASE("DYNAMIC.COLUMN")
      Ok=.TRUE.
      Do02: DO
        !
        READ(F,'(A)',IOSTAT=ios) L
        IF(ios/=0) EXIT Do01
        !
        CALL LinToWrd(L,W1,EoL)
        CALL AppendToEnd(L,W1,EoL)
        !
        IF(W1(1:1)=='!') CYCLE Do02
        !
        CALL LinToWrd(L,W2,EoL)
        !-> W2 is second word on line
        !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
        !
        SELECT CASE(TRIM(W1))
        !
        CASE("ENDINPUT"); EXIT Do01
        CASE("END","ENDDYNAMIC.COLUMN"); EXIT Do02
        !
        !-------------------------------------------spatial parameters--
        CASE("DX")      ;   CALL WrdToReal(W2,Column%dX)
        CASE("DT")      ;   CALL WrdToReal(W2,Column%dT)
        CASE("DISP")    ;   CALL WrdToReal(W2,Column%Disp)
        CASE("UDARCY")  ;   CALL WrdToReal(W2,Column%UDarcy)
        CASE("TFIN")    ;   CALL WrdToReal(W2,Column%Duration)
        CASE("TSAVE")   ;   CALL WrdToReal(W2,Column%Time_Save)
        CASE("NCELL")   ;   CALL WrdToInt (W2,Column%nCell)
        CASE("METHOD")  ;   CALL WrdToInt (W2,Column%Method)
        !---------------------------------------------------------------
        !
        CASE DEFAULT
          Ok= .FALSE.
          Msg= "In DYNAMIC.COLUMN, "//TRIM(W1)//" is unknown KeyWord !!"
        !
        END SELECT
      ENDDO Do02
    !
    END SELECT
    !
  ENDDO Do01
  !
  CALL CLOSEFILE(F)
  
  IF(.NOT. Ok) Msg= "Block DYNAMIC.COLUMN not found !!!"
  
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dynam_ReadColumn"
  
  RETURN
ENDSUBROUTINE Dynam_ReadColumn

SUBROUTINE Dynam_ReadNumeric(Ok)
!--
!-- read block DYNAMIC.NUMERIC -> numerical options
!--
  USE M_Files,ONLY: NamFInn
  USE M_IOTools
  USE M_Dynam_Vars,ONLY: &
  & iCtrlTime,cMethod,Implicit_Surface, &
  & LogForMin,LogForAqu,DirectSub, &
  & bFinDIF,DebNewt,TestJacob,TestMax, &
  & NewtMaxIts,NewtIterMax,NewtIterMin, &
  & NewtTOLF,NewtTOLX
  !
  LOGICAL,INTENT(OUT)::Ok
  !
  CHARACTER(LEN=255):: L
  CHARACTER(LEN=255) :: W0,W1,W2
  LOGICAL           :: EoL
  INTEGER           :: F,ios
  !
  CALL GetUnit(F)
  CALL OPENFILE (F,FILE=TRIM(NamFInn))
  Ok=.FALSE.
  
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< Dynam_ReadNumeric"
    
  Do01: DO
    !
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT Do01
    CALL LinToWrd(L,W0,EoL); CALL AppendToEnd(L,W0,EoL)
    IF(W0(1:1)=='!')   CYCLE Do01 !skip comment lines
    !
    SELECT CASE(TRIM(W0))
    !
    CASE("ENDINPUT")  ;  EXIT Do01
    !
    !------------------------------------- param's for numeric method --
    CASE("DYNAMIC.NUMERIC","DYN.NUMERIC")
      !
      Ok=.TRUE.
      !
      DoReadNum: DO
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT Do01
        !READ(F,'(A)') L
        CALL LinToWrd(L,W1,EoL)
        IF(W1(1:1)=='!') CYCLE DoReadNum
        CALL AppendToEnd(L,W1,EoL)
        !
        ! print *,'Dynam_ReadNumeric,W1=',TRIM(W1)   ;  pause
        
        SELECT CASE(TRIM(W1))
        CASE("ENDINPUT")
          EXIT Do01
        CASE("END","ENDDYNAMIC.NUMERIC","ENDDYN.NUMERIC")
          EXIT DoReadNum
        END SELECT
        !
        CALL LinToWrd(L,W2,EoL)
        !-> W2 is second word on line
        !-> either second keyword or numeric' param' (or 'empty', i.e. "!")
        !
        ! print *,'Dynam_ReadNumeric,W2=',TRIM(W2)  ;  pause
        
        SELECT CASE(TRIM(W1))
        !!CASE("STEADY.STOP");
        CASE("SURFACE")
          SELECT CASE(TRIM(W2))
            CASE("IMPLICIT")  ; Implicit_Surface= .TRUE.
            CASE("EXPLICIT")  ; Implicit_Surface= .FALSE.
            CASE DEFAULT
              CALL Stop_(TRIM(W2)//"-> unknown keyword (IMPLICIT,EXPLICIT)")
          END SELECT
        !
        CASE("MAXITER")  ; CALL WrdToInt(W2,NewtMaxIts)
        CASE("ITERMAX")  ; CALL WrdToInt(W2,NewtIterMax)
        CASE("ITERMIN")  ; CALL WrdToInt(W2,NewtIterMin)
        !
        !IF(iter>NewtIterMax) dT0=dT0/2.0D0
        !IF(iter>NewtIterMax) dT0=dT0/2.0D0
        !IF(iter<NewtIterMin) dT0=dT0*2.0D0
        !
        CASE("CONTROL")  ; CALL WrdToInt(W2,iCtrlTime)   !
        !
        CASE("METHOD")  !!THEN
          SELECT CASE(TRIM(W2))
          CASE("NEWTON")       ;   cMethod= TRIM(W2)  ! iMethod= 1
          CASE("NEWTLNSRCH")   ;   cMethod= TRIM(W2)  ! iMethod= 2
          CASE("NEWTONPRESS")  ;   cMethod= TRIM(W2)  ! iMethod= 2
          CASE("BROYDEN")      ;   cMethod= TRIM(W2)  ! iMethod= 3
          CASE("NEWTONCHESS")  ;   cMethod= TRIM(W2)  ! iMethod= 4
          CASE("TENSOLVE_1")   ;   cMethod= TRIM(W2)  ! iMethod= 5
          CASE("NEWTONKELLEY") ;   cMethod= TRIM(W2)  ! iMethod= 6
          CASE("NEWTONWALKER") ;   cMethod= TRIM(W2)  ! iMethod= 7 !-> default value
          CASE DEFAULT         ;   CALL Warning(W2)
          END SELECT
        !_
        CASE("NEWTTOLF")
          !convergence criterion on FUNCTION values
          CALL WrdToReal(W2,NewtTolF)
        !CASE("NEWTTOLMIN")
        !  !criterion for spurious convergence
        !  CALL WrdToReal(W2,NewtTolF)
        !CASE("NEWTTOLX")
        !  !convergence criterion on dx
        !  CALL WrdToReal(W2,NewtTolX)
        !
        !------------------------------------------------- obsolete --
        CASE("NOTLOG")  ;       LogForMin=.FALSE.
        CASE("LOGMIN")  ;       LogForMin=.TRUE.
        CASE("FINDIF")  ;       bFinDIF=  .TRUE.
        !
        CASE("TESTMAX")  ;      TestMax=  .TRUE.
        !------------------------------------------------/ obsolete --
        !
        CASE("LOGFORAQU")
          SELECT CASE(TRIM(W2))
            CASE("TRUE", "T", "YES")  ;  LogForAqu= .TRUE.
            CASE("FALSE", "F", "NO")  ;  LogForAqu= .FALSE.
          ENDSELECT

        CASE("LOGFORMIN")
          SELECT CASE(TRIM(W2))
          CASE("TRUE", "T", "YES")  ;  LogForMin= .TRUE.
          CASE("FALSE", "F", "NO")  ;  LogForMin= .FALSE.
          ENDSELECT

        CASE("DIRECT")
          SELECT CASE(TRIM(W2))
          CASE("TRUE", "T", "YES")  ;  DirectSub= .TRUE.
          CASE("FALSE", "F", "NO")  ;  DirectSub= .FALSE.
          ENDSELECT

        CASE("JACOBIAN")
          SELECT CASE(TRIM(W2))
          CASE("NUMERIC")   ;  bFinDif= .TRUE.
          CASE("ANALYTIC")  ;  bFinDif= .FALSE.
          ENDSELECT

        CASE DEFAULT
          CALL Warning(W1)

        END SELECT
        
        ! print *,"W1=",TRIM(W1)  ;  pause

      ENDDO DoReadNum
    !_CASE("DYNAMIC.NUMERIC")
    !---/
    END SELECT
    !
  ENDDO Do01
  !
  CALL CLOSEFILE(F)
  !
  DebNewt=   (iDebug>2)
  TestJacob= (iDebug>2)
  !
  IF(.NOT. &
  &       (NewtMaxIts>NewtIterMax &
  & .AND. NewtIterMax>NewtIterMin)) THEN
    IF(iDebug>0) WRITE(fTrc,'(A)') &
    & "Must Have MAXITER>ITERMAX .AND. ITERMAX>ITERMIN"
    CALL Stop_("Must Have MAXITER>ITERMAX .AND. ITERMAX>ITERMIN")
  ENDIF
  !
  IF(iCtrlTime<1 .OR. iCtrlTime>2) iCtrlTime=1
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dynam_ReadNumeric"
  !
  RETURN
ENDSUBROUTINE Dynam_ReadNumeric

SUBROUTINE Warning(W)
  CHARACTER(LEN=*),INTENT(IN):: W
  IF(iDebug>2) &
  & WRITE(fTrc,'(A)') "!!WARNING!! "//TRIM(W)//" <-unknown KeyWord !!"
  IF(iDebug>2) THEN
    PRINT '(A)',"!!WARNING!! "//TRIM(W)//" <-unknown KeyWord !!"
    CALL Pause_
  ENDIF
  RETURN
ENDSUBROUTINE Warning

SUBROUTINE Check_KinRate_CalcQsK
!--
!-- check kinetic models
!--
  USE M_IoTools,   ONLY: GetUnit, Closefile, OpenFile
  USE M_Files,     ONLY: DirLog,Files_Index_Write
  USE M_Global_Vars,ONLY: vKinModel, nAq
  USE M_KinRate,   ONLY: KinRate_CalcQsKFactor
  !
  REAL(dp):: Aff,Qsk,VmQsk
  INTEGER :: I,J,F
  REAL(dp),DIMENSION(1:nAq):: vDX,vX
  !
  CALL GetUnit(F)
  CALL OPENFILE(F,FILE=TRIM(DirLog)//"kinrate_vmqsk.log")
  !
  CALL Files_Index_Write(fHtm,&
  & TRIM(DirLog)//"kinrate_vmqsk.log",&
  & "DYNAMIC/LOG: check KinModels")
  !
  WRITE(F,'(2(A,A1))',ADVANCE="NO") "Ord",T_,"Aff",T_
  DO I=1,SIZE(vKinModel)
    WRITE(F,'(A,A1)',ADVANCE="NO") vKinModel(I)%Name,T_
  ENDDO
  WRITE(F,*)
  !
  Aff= Zero
  vX=  Zero
  J=1
  DO
    WRITE(F,'(I3,A1,G15.6,A1)',ADVANCE="NO") J,T_,Aff,T_
    QsK=EXP(-Aff)
    DO I=1,SIZE(vKinModel)
      !
      CALL KinRate_CalcQsKFactor(&
      & "D",          & !IN was "DISSOLU"
      & vKinModel(I), & !IN: Kinetic model
      & QsK,          & !IN
      & vX,           & !IN !dQsKdLnXi
      & VmQsK,        & !OUT
      & vDX)            !OUT dVmQdLnX_M
      !
      WRITE(F,'(G15.6,A1)',ADVANCE="NO") VmQsK,T_
    ENDDO
    WRITE(F,*)
    !
    Aff= Aff +0.1_dp; J=J+1
    IF(Aff>6.0_dp) EXIT
    !
  ENDDO
  !
  CALL CLOSEFILE(F)
  !
ENDSUBROUTINE Check_KinRate_CalcQsK

ENDMODULE M_Dynam_Read

