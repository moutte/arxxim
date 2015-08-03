MODULE M_Dtb_Read_DtbMinHkf
  USE M_Kinds
  USE M_IoTools
  USE M_Trace
  USE M_Dtb_Read_Vars
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: DtbMinHKF_Build
  PUBLIC:: DtbMinHKF_Read_Old
  PUBLIC:: DtbMinHKF_Read
  !
  TYPE(T_LisMinHkf),POINTER,SAVE:: LisCur, LisPrev

CONTAINS

SUBROUTINE DtbMinHKF_Build
!--
!-- transfer linked list to array of T_DtbMinHkf
!--
  USE M_Dtb_Vars,ONLY: vDtbMinHkf
  !
  INTEGER:: I
  !
  IF(iDebug>5) PRINT '(A)',"< DtbMinHKF_Build"
  !
  IF(ALLOCATED(vDtbMinHkf)) DEALLOCATE(vDtbMinHkf)
  ALLOCATE(vDtbMinHkf(nMinHkf))
  !
  LisCur=>LisMinHkf
  !
  I=0
  DO WHILE (ASSOCIATED(LisCur))
    I=I+1
    vDtbMinHkf(I)=LisCur%Value
    IF(iDebug>5) PRINT '(2A)', vDtbMinHkf(I)%Name,TRIM(vDtbMinHkf(I)%Num) !; CALL PAUSE_
    LisPrev=>LisCur
    LisCur=> LisCur%next
    DEALLOCATE(LisPrev)
  ENDDO
  !
  IF(iDebug>5) PRINT '(A)',"</ DtbMinHKF_Build"
  IF(iDebug>5) CALL Pause_
  RETURN
ENDSUBROUTINE DtbMinHKF_Build

SUBROUTINE DtbMinHKF_Read(F,vEle,N)
!--
!-- reads data for minerals & gases directly from (modified) slop*.dat, or obigt.dat
!-- "new" formats are without subblocks for minerals / gas:
!-- SPECIES HSV.HKF
!--   MIN M001 FORSTERITE ../..
!--   MIN ../..
!--   GAS G001 CO2,G ../..
!--   GAS ../..
!-- END SPECIES
!--
!-- "old" formats have separate blocks for MINERAL and GAS:
!-- SPECIES MIN.HKF
!-- MINERAL
!--   M001 FORSTERITE ../..
!--   ../..
!-- END MINERAL
!-- GAS
!--   G001 CO2,G ../..
!--   ../..
!-- END GAS
!-- END SPECIES
!--
!-- NB !!! parameters for phase transition are not read !!!
!--
  USE M_Files,      ONLY: DirDtbLog,Files_Index_Write
  USE M_Dtb_Const,  ONLY: CalToJoule,Tref,DtbConv_Benson
  USE M_T_Element,  ONLY: T_Element,Formula_Read,Formula_Build,Element_Index
  USE M_Dtb_Read_Vars,ONLY: LisMinHkf
  USE M_T_DtbMinHKF
  USE M_Dtb_Read_Tools
  !
  USE M_Dtb_Vars,   ONLY: DtbFormat
  !
  INTEGER,INTENT(IN):: F !input file
  TYPE(T_Element),INTENT(IN):: vEle(:)
  INTEGER,INTENT(INOUT):: N
  !
  CHARACTER(LEN=255):: L,W,sFormul
  CHARACTER(LEN=4)  :: ICode
  CHARACTER(LEN=10) :: sFormat
  TYPE(T_DtbMinHkf) :: M
  REAL(dp),DIMENSION(dimV)::vX
  LOGICAL :: EoL, fOk, SubBlock
  ! LOGICAL :: bMin, bGas
  INTEGER :: I,ios,K,ZSp,FF,Div,iEl,ieO_
  REAL(dp):: H_Calc
  !
  INTEGER,ALLOCATABLE::vStoik(:) ! for Formula_Read
  !
  !--- for header processing --
  LOGICAL :: IsHeader
  INTEGER,PARAMETER:: nField= 10
  CHARACTER(LEN=12):: vStrField(1:nField)
  ! vStrField contains names of all possible fields in a database
  INTEGER :: vIField(1:nField)
  ! CHARACTER(LEN=12):: vStrUnit(1:nField)
  !---/ for header processing --
  !
  !--- for formula translation --
  CHARACTER(LEN=6):: CodFormula ! SCFORM | ECFORM
  ! SCFORM = compact formula,  e.g. SiO2
  ! ECFORM = extended formula, e.g. SI(1)O(2)
  CHARACTER(LEN=2),ALLOCATABLE:: vElement(:)
  LOGICAL:: EcformIsOk
  INTEGER:: fFormula
  !---/
  !
  ! IF (N>1) LisCur => LisMinHkf
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DtbMinHKF_Read"
  !
  !----------------------------------------------------------------trace
  FF= 0
  IF(iDebug>2) THEN
    WRITE(fTrc,'(A)') &
    & "stoikio in "//TRIM(DirDtbLog)//"minhkf_stoik.log"
    !
    CALL GetUnit(FF)
    OPEN(FF,FILE=TRIM(DirDtbLog)//"minhkf_stoik.log")
    !
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirDtbLog)//"minhkf_stoik.log",&
    & "check species stoikio of min.species from supcrt-like database")
    !
    WRITE(FF,"(A15,A1)",ADVANCE="NO") "NAME",T_
    DO iEl=1,SIZE(vEle)
      WRITE(FF,"(A3,A1)",ADVANCE="NO") vEle(iEl)%NamEl,T_
    ENDDO
    WRITE(FF,"(A7,A1,A7)") "Div",T_,"FORMULA"
    !
  ENDIF
  !---------------------------------------------------------------/trace
  !
  SubBlock= .FALSE.
  sFormat= "HKF_CAL"
  !
  !--- for formula translation --
  CodFormula= "ECFORM"
  ALLOCATE(vElement(SIZE(vEle)+3))
  CALL DtbRead_Build_vElement(vEle,vElement)
  fFormula= 0
  !---/
  !
  !~ IF(iDebug>3) WRITE(fTrc,'(6(A,A1))') &
  !~ & "M%Name",T_, "M%G0R",T_, "M%S0Ele",T_, "M%H0R",T_, "H_Calc",T_, "M%H0R-H_Calc",T_
  !
  ALLOCATE(vStoik(1:SIZE(vEle)))
  ieO_= Element_Index("O__",vEle)
  !
  FilCode=TRIM(DtbFormat) !default, also read from file
  !
  !----------------------------------initialize vStrField(:), vIField(:)
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
  !
  !--- scan default field list
  !~ IF (FilCode(1:5)=="OBIGT") THEN
    !~ L= "SOURCE NAME SKIP SCFORM SKIP SKIP SKIP SKIP PARAMETERS"
  !~ ELSE ! case of SLOP98.DAT
    !~ L= "SOURCE NAME SKIP SKIP ECFORM SKIP SKIP PARAMETERS"
  !~ ENDIF
  !~ !
  !~ !IF(iDebug==4) PRINT *,"< Default values"
  !~ CALL FieldList_Read(L,vStrField,vIField)
  !~ !IF(iDebug==4) PRINT *,"</ Default values"
  !---/ scan default field list
  !
  IF(vIField(4)/=0 .AND. fFormula==0 .AND. iDebug>2) THEN
  ! -> files contains compact formulas
    CALL GetUnit(fFormula)
    OPEN(fFormula,file="debug_formula.log")
    WRITE(fFormula,'(A,/)') "resuts of formula conversion"
  ENDIF
  !---------------------------------/initialize vStrField(:), vIField(:)
  !
  DoFile: DO
    !
    CALL DtbMinHkf_Zero(M)
    !
    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=="!") CYCLE DoFile
    CALL AppendToEnd(L,W,EoL)
    !
    !------------------------------------------- read header, if present
    !--------- if the line begins with any member of array vStrField(:),
    !--------------------------------- then it contains the line headers
    IsHeader= .FALSE.
    DO I=1,nField
      IF(TRIM(W)==TRIM(vStrField(I))) THEN
        IsHeader= .TRUE.
        EXIT
      ENDIF
    END DO
    !
    IF(IsHeader) THEN
      L= TRIM(W)//" "//TRIM(L)
      !IF(iDebug==4) PRINT *,"< values from file"
      CALL FieldList_Read(L,vStrField,vIField)
      !IF(iDebug==4) PRINT *,"</ values from file"

      IF(vIField(4)/=0 .AND. fFormula==0 .AND. iDebug>2) THEN
        CALL GetUnit(fFormula)
        OPEN(fFormula,file="debug_formula.log")
        WRITE(fFormula,'(A,/)') "results of formula conversion"
      ENDIF

      CYCLE DoFile
    ENDIF
    !-------------------------------------------/read header, if present
    !
    !-----------------------------------------process first word of line
    SELECT CASE(W)
    !
    CASE("ENDINPUT")
      EXIT DoFile
    !
    CASE("ENDSPECIES")
      EXIT DoFile
    !
    CASE("END")
      EXIT DoFile
    !
    !----------------------------------------------------for old formats
    !old! CASE("END")
    !old!   IF(SubBlock) THEN
    !old!     SubBlock= .FALSE.
    !old!     CYCLE DoFile
    !old!   ELSE
    !old!     EXIT DoFile
    !old!   ENDIF
    !old! CASE("ENDMINERAL","ENDGAS")
    !old!   SubBlock= .FALSE.
    !old!   CYCLE DoFile
    !old! !
    !old! CASE("MINERAL")
    !old!   SubBlock= .TRUE.
    !old!   M%Typ= "MIN"
    !old!   CYCLE DoFile
    !old! CASE("GAS")
    !old!   SubBlock= .TRUE.
    !old!   M%Typ= "GAS"
    !old!   CYCLE DoFile
    !old! CASE("CODE")
    !old!   !-> can change the CODE, i.e. the DATA source inside the block
    !old!   CALL LinToWrd(L,W,EoL)
    !old!   FilCode= TRIM(W)
    !old!   !
    !old!   ! if FilCode changes, must re-initialize vIField ...!
    !old!   IF (FilCode(1:5)=="OBIGT") THEN
    !old!     L= "INDEX NAME SKIP SCFORM TYPE SKIP SKIP SKIP PARAMETERS"
    !old!   ELSE ! case of SLOP98.DAT
    !old!     L= "INDEX NAME SKIP SKIP ECFORM SKIP SKIP PARAMETERS"
    !old!   ENDIF
    !old!   CALL FieldList_Read(L,vStrField,vIField)
    !old!   !
    !old!   CYCLE DoFile
    !----------------------------------------------/for old formats --
    !
    CASE DEFAULT
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= TRIM(W)//" "//TRIM(L)
      !
    ENDSELECT
    !------------------------------------/ process first word of line --
    !
    M%Num=TRIM(FilCode)
    !
    !--------------------------------------------- scan the data line --
    !--------------------------------- up to column before PARAMETERS --
    DO I= 1,vIField(10)-1
      !
      !vStrField(1:nField)= &
      !!!!!"____________","____________","____________","____________","____________",
      !& (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
      !&   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
      !
      CALL LinToWrd(L,W,EoL,"NO")
      !
      IF(EoL) CYCLE DoFile ! data line, contains no sufficient data -> skip !!
      !
      IF(I==vIField(7)) THEN  !CODE / SOURCE
        CALL Str_Upper(W)  ;  FilCode= TRIM(W)
      END IF

      IF(I==vIField(8)) THEN  !FORMAT
        CALL Str_Upper(W)  ;  sFormat= TRIM(W)
      END IF

      IF(I==vIField(1)) THEN  !TYPE
        !old! CALL Str_Upper(W)  ;  M%Typ= TRIM(W)
        CALL Str_Upper(W)
        SELECT CASE(TRIM(W))
        ! CASE("AQU")  ;  M%Typ="AQU"
        CASE("MIN")  ;  M%Typ="MIN"
        CASE("GAS")  ;  M%Typ="GAS"
        CASE("LIQ")  ;  M%Typ="LIQ"
        CASE DEFAULT ;  CALL Stop_(TRIM(W)//"<< unknown TYPE in database !!...")          
        END SELECT
      END IF

      IF(I==vIField(2)) THEN  !INDEX
        CALL Str_Upper(W)  ;  M%Num= TRIM(W)
      END IF

      IF(I==vIField(3)) THEN  !NAME
        CALL Str_Upper(W)  ;  M%Name= TRIM(W)
      ENDIF

      IF(I==vIField(4)) THEN  !SCFORM
        !
        CALL DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
        IF(.NOT.EcformIsOk) CYCLE DoFile !-------------------------CYCLE
        !
        CALL Str_Upper(W)  ;  M%Formula=TRIM(W)
      ENDIF

      IF(I==vIField(5)) THEN  !ECFORM
        CALL Str_Upper(W)  ;  M%Formula=  TRIM(W)
      END IF

    END DO
    !------------------------ scan the data line (left to PARAMETERS) --
    !
    !--------------------------------------------------- read formula --
    CALL Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
    IF(.NOT. fOk) CYCLE DoFile !==============================< CYCLE ==
    !--------------------------------------------------/ read formula --
    !
    M%S0Ele=  DOT_PRODUCT(vStoik(:), vEle(:)%S0) /REAL(Div) !in Joule
    M%WeitKg= DOT_PRODUCT(vStoik(:), vEle(:)%WeitKg) /REAL(Div)
    M%Div= Div
    !
    CALL ReadRValsV(L,K,vX)
    !
    IF (FilCode(1:5)=="OBIGT") THEN
      ! in OBIGT database, Cp is tabulated fot min'species
      ! but not used as input parameters -> vX(4) ignored
      M%G0R= vX(1)
      M%H0R= vX(2)
      M%S0_= vX(3)
      M%V0R= vX(5)
      M%MK1(1)= vX(6)
      M%MK1(2)= vX(7) *1.0D-3
      M%MK1(3)= vX(8) *1.0D5
      M%NTran= 0
    ELSE
      ! default FilCode is SLOP
      M%G0R=    vX(1)
      M%H0R=    vX(2)
      M%S0_=    vX(3)
      M%V0R=    vX(4)
      M%MK1(1)= vX(5)
      M%MK1(2)= vX(6) *1.0D-3
      M%MK1(3)= vX(7) *1.0D5
      M%NTran= 0
    ENDIF
    !
    IF(sFormat(1:7)=="HKF_CAL") THEN
      M%G0R= M%G0R*CalToJoule
      M%H0R= M%H0R*CalToJoule
      M%S0_= M%S0_*CalToJoule
      M%MK1(1)= M%MK1(1)
      M%MK1(2)= M%MK1(2)*CalToJoule
      M%MK1(3)= M%MK1(3)*CalToJoule
    ENDIF
    !
    N=N+1
    !
    CALL IntToStr4(N,ICode)
    M%Num= TRIM(FilCode)//"_"//TRIM(ICode)
    !
    IF(.NOT. DtbConv_Benson) M%G0R= M%G0R - Tref *M%S0Ele !!!Berman Convention!!!
    !
    !----------------------------------------- "fill" the linked list --
    !CALL Save_Record(N,M)
    IF(N==1) THEN
      ALLOCATE(LisMinHkf)
      NULLIFY(LisMinHkf%next)
      LisMinHkf%Value= M
      LisCur=> LisMinHkf
    ELSE
      ALLOCATE(LisCur%next)
      NULLIFY(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=>LisCur%next
    ENDIF
    !-----------------------------------------/"fill" the linked list --
    !
    H_Calc= M%G0R + Tref *M%S0_ - Tref*M%S0Ele
    !~ IF(iDebug>3) WRITE(fTrc, '(A,A1,5(G15.6,A1))') &
    !~ & M%Name,T_, M%G0R,T_, M%S0Ele,T_, M%H0R,T_, H_Calc,T_, M%H0R-H_Calc,T_ ! M%H0R-H_Calc
    !
    IF(FF>0) THEN
      WRITE(FF,"(A15,A1)",ADVANCE="NO") M%Name,T_
      DO iEl=1,SIZE(vEle)
        WRITE(FF,"(I7,A1)",ADVANCE="NO") vStoik(iEl),T_
      ENDDO
      WRITE(FF,'(4(G15.6,A1))',ADVANCE="NO") &
      & M%G0R,T_,H_Calc,T_,M%H0R,T_,M%S0_,T_
      !
      ! build new formula, with fixed element order
      CALL Formula_Build(vEle,vStoik,Zsp,M%Div,sFormul)
      !!sFormul=""
      !!DO iEl=1,SIZE(vEle)
      !!  IF(vStoik(iEl)>0) sFormul=TRIM(sFormul)//TRIM(vEle(iEl)%NamEl)//""
      !!ENDDO
      WRITE(FF,"(I3,A1,A39,A1,A15,A1,A)") &
      & M%Div,T_,M%Formula,T_,M%Name,T_,TRIM(sFormul)
    ENDIF
    !
  END DO DoFile
  !
  DEALLOCATE(vStoik)
  !
  DEALLOCATE(vElement)
  IF(fFormula>0) CLOSE(fFormula)
  !
  IF(FF>0) CLOSE(FF)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DtbMinHKF_Read"
  !
  RETURN
ENDSUBROUTINE DtbMinHKF_Read

SUBROUTINE Save_Record(N,M)
  INTEGER,          INTENT(IN):: N
  TYPE(T_DtbMinHkf),INTENT(IN):: M
  !
    IF(N==1) THEN
      ALLOCATE(LisMinHkf)
      NULLIFY(LisMinHkf%next)
      LisMinHkf%Value= M
      LisCur=> LisMinHkf
    ELSE
      ALLOCATE(LisCur%next)
      NULLIFY(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=>LisCur%next
    ENDIF
  !
END SUBROUTINE Save_Record

SUBROUTINE DtbMinHKF_Read_Old(F,vEle,N)
!--
!-- reads data for minerals directly from (modified) slop*.dat, or obigt.dat
!-- this is for "old" formats, with separate blocks for MINERAL and GAS:
!-- SPECIES MIN.HKF
!-- MINERAL
!--   ..
!-- END MINERAL
!-- GAS
!--   ..
!-- END GAS
!-- END SPECIES
!-- NB !!! parameters for phase transition are not read !!!
!--
  USE M_Dtb_Const,  ONLY: CalToJoule,Tref,DtbConv_Benson
  USE M_T_Element,  ONLY: T_Element,Formula_Read,Formula_Build,Element_Index
  USE M_Files,      ONLY: DirDtbLog,Files_Index_Write
  USE M_Dtb_Vars,   ONLY: DtbFormat
  USE M_Dtb_Read_Tools
  !
  INTEGER,INTENT(IN):: F !input file
  TYPE(T_Element),INTENT(IN):: vEle(:)
  INTEGER,INTENT(INOUT):: N
  !
  !TYPE(T_LisMinHkf),POINTER,SAVE:: LisCur
  !
  CHARACTER(LEN=255):: L,W,sFormul
  CHARACTER(LEN=4)  :: ICode
  TYPE(T_DtbMinHkf) :: M
  REAL(dp),DIMENSION(dimV)::vX
  LOGICAL :: EoL, fOk
  LOGICAL :: bMin, bGas
  INTEGER :: ios,K,ZSp,FF,Div,iEl,ieO_
  REAL(dp):: H_Calc
  !
  INTEGER,ALLOCATABLE::vStoik(:) !for Formula_Read
  !
  !--- for formula translation --
  CHARACTER(LEN=6):: CodFormula ! SCFORM | ECFORM
  ! SCFORM = compact formula,  e.g. SiO2
  ! ECFORM = extended formula, e.g. SI(1)O(2)
  CHARACTER(LEN=2),ALLOCATABLE:: vElement(:)
  LOGICAL:: EcformIsOk
  INTEGER:: fFormula
  !---/
  !
  ! IF (N>1) LisCur => LisMinHkf
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DtbMinHKF_Read"
  !
  !----------------------------------------------------------------trace
  FF= 0
  IF(iDebug>2) THEN
    WRITE(fTrc,'(A)') &
    & "stoikio in "//TRIM(DirDtbLog)//"minhkf_stoik.log"
    !
    CALL GetUnit(FF)
    OPEN(FF,FILE=TRIM(DirDtbLog)//"minhkf_stoik.log")
    !
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirDtbLog)//"minhkf_stoik.log",&
    & "check species stoikio of min.species from supcrt-like DATAbase")
    !
    WRITE(FF,"(A15,A1)",ADVANCE="NO") "NAME",T_
    DO iEl=1,SIZE(vEle)
      WRITE(FF,"(A3,A1)",ADVANCE="NO") vEle(iEl)%NamEl,T_
    ENDDO
    WRITE(FF,"(A7,A1,A7)") "Div",T_,"FORMULA"
    !
  ENDIF
  !---------------------------------------------------------------/trace
  !
  !--- for formula translation --
  CodFormula= "ECFORM"
  ALLOCATE(vElement(SIZE(vEle)+3))
  CALL DtbRead_Build_vElement(vEle,vElement)
  fFormula= 0
  !---/
  !
  !~ IF(iDebug>3) WRITE(fTrc,'(6(A,A1))') &
  !~ & "M%Name",T_, "M%G0R",T_, "M%S0Ele",T_, "M%H0R",T_, "H_Calc",T_, "M%H0R-H_Calc",T_
  !
  ALLOCATE(vStoik(1:SIZE(vEle)))
  ieO_= Element_Index("O__",vEle)
  !
  FilCode=TRIM(DtbFormat) !default, also read from file
  !
  DoFile: DO

    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=="!") CYCLE DoFile
    CALL AppEndToEnd(L,W,EoL)
    !
    SELECT CASE(W)
      !
      CASE("END","ENDSPECIES","ENDINPUT")
        EXIT DoFile
        !
      CASE("CODE")
        CALL LinToWrd(L,W,EoL)
        FilCode=TRIM(W)
        CYCLE DoFile
        !
      CASE("FORMULA")
        CALL LinToWrd(L,W,EoL)
        CodFormula= TRIM(W)
        IF(iDebug>2 .and. CodFormula=="SCFORM") THEN
          CALL GetUnit(fFormula)
          OPEN(fFormula,file="debug_formula.log")
          WRITE(fFormula,'(A,/)') "resuts of formula conversion"
        ENDIF
        CYCLE DoFile
        !
      CASE("MINERAL","GAS")
        !
        SELECT CASE(W)
          CASE("MINERAL"); bGas=.FALSE.; bMin=.TRUE.
          CASE("GAS");     bGas=.TRUE.;  bMin=.FALSE.
        END SELECT
        !
        DoReadMin: DO
          !
          READ(F,'(A)',IOSTAT=ios) L
          IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=="!") CYCLE DoReadMin
          CALL AppendToEnd(L,W,EoL)
          !
          SELECT CASE(W)
            CASE("ENDINPUT","ENDSPECIES")
              EXIT DoFile
            CASE("END","ENDMINERAL","ENDGAS")
              CYCLE DoFile !EXIT DoReadMin
          END SELECT
          !
          M%Num=TRIM(W)
          IF(bGas) M%Typ="GAS"
          IF(bMin) M%Typ="MIN"
          !
          IF (FilCode(1:5)=="OBIGT") THEN
            !
            !L= "SOURCE NAME SKIP SCFORM SKIP SKIP SKIP SKIP PARAMETERS"
            CodFormula="SCFORM"
            !
            CALL LinToWrd(L,W,EoL); M%Name=TRIM(W) !IF(iDebug>0) WRITE(fTrc,"(A)") M%Name
            CALL LinToWrd(L,W,EoL) !skip ABBRV
            CALL LinToWrd(L,W,EoL,"NO") !-> compact formula, CHARACTER CASE is conserved :!!
            !
            IF(CodFormula=="SCFORM") THEN
              CALL DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
              IF(.not.EcformIsOk) THEN
                IF(iDebug>3) CALL Warning_("!!! Cannot translate "//TRIM(W))
                CYCLE DoReadMin
              ENDIF
            ENDIF
            CALL Str_Upper(W)
            !
            M%Formula=TRIM(W)
            !
            CALL LinToWrd(L,W,EoL) !skip STATE
            CALL LinToWrd(L,W,EoL) !skip SOURCE1
            CALL LinToWrd(L,W,EoL) !skip SOURCE2
            CALL LinToWrd(L,W,EoL) !DATE
            !
            CALL ReadRValsV(L,K,vX)
            !
          ELSE ! default FilCode is SLOP98
            !
            !L= "SOURCE SKIP NAME SKIP ECFORM SKIP SKIP PARAMETERS"
            !
            CALL LinToWrd(L,W,EoL) !skip abredged name ; M%Abbr=TRIM(W)
            CALL LinToWrd(L,W,EoL) ; M%Name=TRIM(W)
            CALL LinToWrd(L,W,EoL) !skip scform
            CALL LinToWrd(L,W,EoL) ; M%Formula=TRIM(W)
            CALL LinToWrd(L,W,EoL) !skip REF
            CALL LinToWrd(L,W,EoL) !skip DATE
            !
            CALL ReadRValsV(L,K,vX)
            !
          ENDIF
          !
          !-------------------------------------------------read formula
          CALL Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
          IF(.NOT. fOk) CYCLE DoReadMin !--------------------------CYCLE
          !------------------------------------------------/read formula
          CALL Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
          !
          IF (FilCode(1:5)=="OBIGT") THEN
            ! in OBIGT database, Cp is tabulated fot min'species
            ! but not used as input parameters -> vX(4) ignored
            M%G0R= vX(1) *CalToJoule
            M%H0R= vX(2) *CalToJoule
            M%S0_= vX(3) *CalToJoule
            M%V0R= vX(5)
            M%MK1(1)= vX(6) *CalToJoule
            M%MK1(2)= vX(7) *CalToJoule *1.0D-3
            M%MK1(3)= vX(8) *CalToJoule *1.0D5
            M%NTran= 0
          ELSE ! default FilCode is SLOP98
            M%G0R=    vX(1) *CalToJoule
            M%H0R=    vX(2) *CalToJoule
            M%S0_=    vX(3) *CalToJoule
            M%V0R=    vX(4)
            M%MK1(1)= vX(5) *CalToJoule
            M%MK1(2)= vX(6) *CalToJoule *1.0D-3
            M%MK1(3)= vX(7) *CalToJoule *1.0D5
            M%NTran= 0
          ENDIF

          !IF(fOk) THEN
          N=N+1
          !
          CALL IntToStr4(N,ICode)
          M%Num= TRIM(FilCode)//"_"//TRIM(ICode)
          !
          M%S0Ele=  DOT_PRODUCT(vStoik(:),vEle(:)%S0) /REAL(Div) ! is in Joule !!
          M%WeitKg= DOT_PRODUCT(vStoik(:),vEle(:)%WeitKg) /REAL(Div)
          M%Div=    Div
          !
          IF(.NOT. DtbConv_Benson) M%G0R= M%G0R - Tref *M%S0Ele !!!Berman Convention!!!
          !
          !--- "fill" the linked list
          IF(N==1) THEN
            ALLOCATE(LisMinHkf)
            NULLIFY(LisMinHkf%next)
            LisMinHkf%Value= M
            LisCur=> LisMinHkf
          ELSE
            ALLOCATE(LisCur%next)
            NULLIFY(LisCur%next%next)
            LisCur%next%Value=M
            LisCur=>LisCur%next
          ENDIF
          !---/
          !
          H_Calc= M%G0R + Tref *M%S0_ - Tref*M%S0Ele
          !~ IF(iDebug>3) WRITE(fTrc, '(A,A1,5(G15.6,A1))') &
          !~ & M%Name,T_, M%G0R,T_, M%S0Ele,T_, M%H0R,T_, H_Calc,T_, M%H0R-H_Calc,T_ ! M%H0R-H_Calc
          !
          IF(FF>0) THEN
            WRITE(FF,"(A15,A1)",ADVANCE="NO") M%Name,T_
            DO iEl=1,SIZE(vEle)
              WRITE(FF,"(I3,A1)",ADVANCE="NO") vStoik(iEl),T_
            ENDDO
            WRITE(FF,'(4(G15.6,A1))',ADVANCE="NO") &
            & M%G0R,T_,H_Calc,T_,M%H0R,T_,M%S0_,T_
            !
            ! build new formula, with fixed element order
            CALL Formula_Build(vEle,vStoik,Zsp,M%Div,sFormul)
            !!sFormul=""
            !!DO iEl=1,SIZE(vEle)
            !!  IF(vStoik(iEl)>0) sFormul=TRIM(sFormul)//TRIM(vEle(iEl)%NamEl)//""
            !!ENDDO
            WRITE(FF,"(I3,A1,A39,A1,A15,A1,A)") &
            & M%Div,T_,M%Formula,T_,M%Name,T_,TRIM(sFormul)
          ENDIF
          !  !
          !ENDIF
        END DO DoReadMin
      !ENDCASE("MINERAL")
    END SELECT !CASE(W0)
  END DO DoFile
  !
  DEALLOCATE(vStoik)
  !
  DEALLOCATE(vElement)
  IF(fFormula>0) CLOSE(fFormula)
  !
  IF(FF>0) CLOSE(FF)
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DtbMinHKF_Read"
  !
  RETURN
ENDSUBROUTINE DtbMinHKF_Read_Old

SUBROUTINE FieldList_Read( &
& L,          &
& vStrField, &
& vIField)
  USE M_IoTools
  !
  CHARACTER(LEN=*),INTENT(INOUT):: L
  CHARACTER(LEN=12),INTENT(IN):: vStrField(:)
  INTEGER,INTENT(OUT):: vIField(:)
  !
  CHARACTER(LEN=80):: W
  LOGICAL:: EoL
  INTEGER:: I,J
  !
  vIField(:)= 0
  I=0
  DO
    CALL LinToWrd(L,W,EoL)
    I=I+1
    DO J=1,SIZE(vStrField)
      IF( TRIM(W)==TRIM(vStrField(J)) ) THEN
        vIField(J)= I
        EXIT
      ENDIF
    END DO
    IF(EoL) EXIT
  END DO
  !
  !! IF(iDebug==4) THEN
  !!   DO I=1,SIZE(vStrField)
  !!     PRINT *,vIField(I),TRIM(vStrField(I))
  !!   ENDDO
  !! ENDIF
  !CALL PAUSE_
  !
  ! (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  !   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
  !
  IF(vIField(10)==0) & ! for "PARAMETERS"
  CALL Stop_( &
  & "in FieldList_Read: keyword not found for "//TRIM(vStrField(10)))

  IF(vIField(1)==0) & ! for "TYPE"  !!MIN/GAS/AQU
  CALL Stop_( &
  & "in FieldList_Read: keyword not found for "//TRIM(vStrField(1)))

  IF(vIField(3)==0) & ! for "NAME"
  CALL Stop_( &
  & "in FieldList_Read: keyword not found for "//TRIM(vStrField(3)))

  IF(vIField(4)==0 .AND. vIField(5)==0) & ! for ECFORM/SCFORM
  CALL Stop_( &
  & "in FieldList_Read: keyword not found for "//TRIM(vStrField(4))//"_"//TRIM(vStrField(5)))
  !
END SUBROUTINE FieldList_Read

ENDMODULE M_Dtb_Read_DtbMinHkf
