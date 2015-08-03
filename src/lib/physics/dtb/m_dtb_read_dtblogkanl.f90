MODULE M_Dtb_Read_DtbLogKAnl
  USE M_Kinds
  USE M_Trace
  USE M_Dtb_Read_Vars,ONLY: T_LisLogKAnl,LisLogKAnl,nLogKAnl
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: DtbLogKAnl_Read
  PUBLIC:: DtbLogKAnl_Build
  !
CONTAINS

SUBROUTINE DtbLogKAnl_Build
!--
!-- transfer linked list to array of T_DtbLogKAnl
!--
  USE M_Dtb_Vars,ONLY: vDtbLogKAnl
  !
  TYPE(T_LisLogKAnl),POINTER:: LisCur, LisPrev
  INTEGER:: I
  !
  IF(iDebug==5) PRINT '(A)',"< DtbLogKAnl_Build"
  !
  IF(ALLOCATED(vDtbLogKAnl)) DEALLOCATE(vDtbLogKAnl)
  ALLOCATE(vDtbLogKAnl(nLogKAnl))
  !
  I=0
  LisCur=> LisLogKAnl
  DO WHILE (ASSOCIATED(LisCur))
    I=I+1
    vDtbLogKAnl(I)= LisCur%Value
    IF(iDebug==5) &
    & PRINT '(2A)', vDtbLogKAnl(I)%Name,TRIM(vDtbLogKAnl(I)%Num)
    LisPrev=>LisCur
    LisCur=> LisCur%next
    DEALLOCATE(LisPrev)
  ENDDO
  !
  IF(iDebug==5) PRINT '(A)',"</ DtbLogKAnl_Build"
  IF(iDebug==5) CALL Pause_
  RETURN
ENDSUBROUTINE DtbLogKAnl_Build

LOGICAL FUNCTION LnkSpc_Found(L,W) !,M)
!--
!-- find index of string W in linked list
!-- return corresponding species Spc
!--
  USE M_T_DtbLogKAnl, ONLY: T_DtbLogKAnl
  !
  TYPE(T_LisLogKAnl),POINTER    :: L
  CHARACTER(LEN=*),  INTENT(IN) :: W
  ! TYPE(T_DtbLogKAnl),INTENT(OUT):: M
  !
  TYPE(T_LisLogKAnl),POINTER:: P,pPrev
  INTEGER::I
  !
  P=> L
  I=  0
  LnkSpc_Found=.FALSE.
  DO WHILE (ASSOCIATED(P))
    I=I+1
    IF(TRIM(W)==TRIM(P%Value%Name)) THEN
      LnkSpc_Found=.TRUE.
      !M= P%Value
      EXIT
    ENDIF
    pPrev=> P
    P=>     P%next
  ENDDO
  !
ENDFUNCTION LnkSpc_Found

SUBROUTINE DtbLogKAnl_Read(F,vEle,N)
!--
!-- build linked list from database in logK/analytic format
!--
  USE M_Dtb_Const,ONLY: T_CK
  USE M_T_Element,ONLY: T_Element,Formula_Read,Formula_Build,Element_Index
  USE M_Files,    ONLY: DirDtbLog,Files_Index_Write
  USE M_IoTools
  USE M_Dtb_Read_Tools
  !
  USE M_T_DtbLogKAnl, ONLY: T_DtbLogKAnl,DtbLogKAnl_New
  USE M_Dtb_Vars,     ONLY: DtbLogK_Dim,DtbLogK_vTPCond
  !
  INTEGER,        INTENT(IN)   :: F !input file
  TYPE(T_Element),INTENT(IN)   :: vEle(:)
  INTEGER,        INTENT(INOUT):: N
  !
  TYPE(T_LisLogKAnl),POINTER,SAVE:: LisCur
  !
  CHARACTER(LEN=512):: L
  CHARACTER(LEN=80):: W
  TYPE(T_DtbLogKAnl) :: M
  LOGICAL :: EoL, fOk
  LOGICAL :: EcformIsOk
  INTEGER :: ios,ZSp,Div,mDum,ff,I
  REAL(dp):: X,Rho
  REAL(dp):: vX(dimV)
  !
  !! CHARACTER(LEN=80):: vWord(1:12)
  !
  !--- for header processing --
  INTEGER,PARAMETER:: nField= 13
  CHARACTER(LEN=12):: vStrField(1:nField)
  ! vStrField contains names of all possible fields in a database
  INTEGER:: vIField(1:nField)
  LOGICAL:: IsHeader
  !
  ! CHARACTER(LEN=12):: vStrUnit(1:nField)
  !---/ for header processing --
  !
  !~ CHARACTER(LEN=512):: FieldList_Default= &
  !~ & "TYPE SOURCE NAME SCFORM SIZE PARAMETERS"
  !
  INTEGER,DIMENSION(:),ALLOCATABLE::vStoik !for Formula_Read
  !
  !--- for formula translation --
  CHARACTER(LEN=6):: CodFormula
  ! CodFormula is either "SCFORM" or "ECFORM" (names inherited from SUPCRT files):
  !   SCFORM = compact formula,  e.g. SiO2
  !   ECFORM = extended formula, e.g. SI(1)O(2)
  CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE:: vElement  !
  INTEGER:: fFormula
  !---/
  !
  CHARACTER(LEN=15):: CodFitting ! type of the fitting function, default: PHREEQ
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DtbLogKAnl_Read"
  !
  !--------------------------------------------------------------- trace
  IF(iDebug>2) THEN
    CALL GetUnit(ff)
    OPEN(ff,FILE="dtb_logKAnl_check.tab")
  ENDIF
  !------------------------------------------------------------- / trace
  !
  fFormula= 0
  CodFormula= "ECFORM"
  !
  !------------------------------ initialize vStrField(:), vIField(:) --
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  ", &
  &   "SIZE        ","VOLUME      ","DENSITY     " /)
  !
  L= "TYPE INDEX NAME SCFORM SIZE PARAMETERS" != default field list
  IF(iDebug==4) PRINT *,"< Default values"
  CALL FieldList_Read(L,vStrField,vIField)
  IF(iDebug==4) PRINT *,"</ Default values"
  !
  IF(vIField(4)/=0 .AND. fFormula==0 .AND. iDebug>2) THEN
  ! -> files contains compact formulas
    CALL GetUnit(fFormula)
    OPEN(fFormula,file="debug_formula.log")
    WRITE(fFormula,'(A,/)') "resuts of formula conversion"
  ENDIF
  !---------------------------------/initialize vStrField(:), vIField(:)
  !
  CodFitting= "PHREEQC"
  !
  !--for formula translation
  ALLOCATE(vElement(SIZE(vEle)+3))
  CALL DtbRead_Build_vElement(vEle,vElement)
  !--/for formula translation
  !
  ALLOCATE(vStoik(1:SIZE(vEle)))
  !
  !------------------------------- build a linked list of all species --
  !----------------------------- consistent with current element list --
  DoFile: DO

    READ(f,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile
    CALL AppendToEnd(L,W,EoL) ! if W=="END" append next word
    !
    CALL DtbLogKAnl_New(M)
    !
    IsHeader= .FALSE.
    !--------------------------------------------read header, if present 
    !--------- if the line begins with any member of array vStrField(:), 
    !--------------------------------- then it contains the line headers 
    DO I=1,nField
      IF(TRIM(W)==TRIM(vStrField(I))) THEN
        IsHeader= .TRUE.
        EXIT
      ENDIF
    END DO
    IF(IsHeader) THEN
      L= TRIM(W)//" "//TRIM(L)
      IF(iDebug>2) PRINT *,"< values from file"
      CALL FieldList_Read(L,vStrField,vIField)
      IF(iDebug>2) PRINT *,"</ values from file"

      IF(vIField(4)/=0 .AND. fFormula==0 .AND. iDebug>2) THEN
        CALL GetUnit(fFormula)
        OPEN(fFormula,file="debug_formula.log")
        WRITE(fFormula,'(A,/)') "resuts of formula conversion"
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
    CASE("END","ENDSPECIES")
      EXIT DoFile
    !
    !---------------------------------------------- for old formats --
    CASE("FITTING")
      CALL LinToWrd(L,W,EoL)
      CodFitting= TRIM(W)
      CYCLE DoFile !---------------------------------------------CYCLE
    !
    CASE("FORMULA")
      CALL LinToWrd(L,W,EoL)
      CodFormula= TRIM(W)
      !
      IF(CodFormula=="SCFORM") THEN
        IF(vIField(4)==0) THEN
          vIField(4)= vIField(5)
          vIField(5)= 0
        ENDIF
        IF(iDebug>2 .AND. fFormula==0) THEN
          CALL GetUnit(fFormula)
          OPEN(fFormula,FILE="debug_formula.log")
          WRITE(fFormula,'(A,/)') "resuts of formula conversion"
        ENDIF
      ENDIF
      !
      CYCLE DoFile !---------------------------------------------CYCLE
    !
    !old! !---------------------------------------------/ for old formats
    !old! CASE DEFAULT
    !old!   !--- ignore line beginning with unknown keyword !!! ??? --
    !old!   CYCLE DoFile !----------------------------------------- CYCLE
    !old! !
    !
    CASE DEFAULT
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= TRIM(W)//" "//TRIM(L)
    !
    ENDSELECT
    !----------------------------------------/process first word of line
    !
    M%Num=     "0"
    M%Name=    "0"
    M%Formula= "0"
    !---------------------------------------------------scan a data line
    DO I= 1,vIField(10)-1
    ! read words up to position before the parameters
      !
      CALL LinToWrd(L,W,EoL,"NO")
      !
      !vStrField(1:12)= &
      !!!!!"____________","____________","____________","____________","____________",
      !& (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
      !&   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  ", &
      !&   "SIZE        ","VOLUME      ","DENSITY     "/)
      !
      IF(I==vIField(1)) THEN  !TYPE
        ! CALL Str_Upper(W)  ;  M%Typ= TRIM(W)
        CALL Str_Upper(W)
        SELECT CASE(TRIM(W))
        CASE("AQU")  ;  M%Typ="AQU"
        CASE("MIN")  ;  M%Typ="MIN"
        CASE("GAS")  ;  M%Typ="GAS"
        CASE("LIQ")  ;  M%Typ="LIQ"
        CASE("XCH")  ;  M%Typ="XCH"
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

      IF(I==vIField(9)) THEN  !FITTING
        CALL Str_Upper(W)  ;  CodFitting= TRIM(W)
      END IF

      IF(I==vIField(11)) CALL WrdToReal(W,X) ! SIZE
      IF(I==vIField(12)) CALL WrdToReal(W,X) ! VOLUME
      IF(I==vIField(13)) CALL WrdToReal(W,X) ! DENSITY

    END DO
    !
    CALL ReadRValsV(L,mDum,vX)
    !-------------------------------------------------/ scan a data line
    !
    !~ IF(vIField(2)/=0)  PRINT *, "M%Num=      ",TRIM(M%Num)
    !~ IF(vIField(3)/=0)  PRINT *, "M%Name=     ",TRIM(M%Name)
    !~ IF(vIField(4)/=0)  PRINT *, "M%Formula=  ",TRIM(M%Formula)
    !~ IF(vIField(5)/=0)  PRINT *, "M%Formula=  ",TRIM(M%Formula)
    !~ IF(vIField(9)/=0)  PRINT *, "CodFitting= ",TRIM(CodFitting)
    !CALL PAUSE_
    !
    !---------------------------- check whether Name not already read --
    IF(LnkSpc_Found(LisLogKAnl,M%Name)) THEN
      IF(iDebug>0) WRITE(fTrc,'(A)') &
      & TRIM(M%Name)//"= name already used -> SKIPPED"
      !PRINT *, TRIM(M%Name)//" <-this name already used by another species ???"
      !CALL PAUSE_
      CYCLE DoFile
    ENDIF
    !---/
    !
    !------------------------------------------------------ read formula
    CALL Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
    IF(.NOT. fOk) CYCLE DoFile !==============================< CYCLE ==
    !-----------------------------------------------------/ read formula
    !
    !-- compute Molecular Weit
    M%WeitKg= DOT_PRODUCT(vStoik(:),vEle(:)%WeitKg) /REAL(Div)
    M%Div= Div
    M%Chg= Zsp
    !
    !------------size parameter for aqu'species, density for min'species
    IF(    vIField(11)/=0 & ! SIZE
    & .OR. vIField(12)/=0 & ! VOLUME, m^3 /mole ?
    & .OR. vIField(13)/=0 & ! DENSITY, Kg /m^3
    & ) THEN
      SELECT CASE(M%Typ)
        CASE("AQU")
          M%AquSize=X
        CASE("MIN","GAS")
          IF(vIField(12)/=0) THEN
            M%V0R= X *1.0D-6 ! for THERMODDEM database, CM^3 to M^3
            Rho= M%WeitKg /M%V0R
          ELSE
            Rho=X
            M%V0R= M%WeitKg /Rho
          ENDIF
      ENDSELECT
    ENDIF
    !-------- default value of size parameter for charged aqu'species --
    IF(M%Typ=="AQU" .and. M%Chg/=0 .and. M%AquSize<=Zero) M%AquSize= 3.72D0
    !-----------/size parameter for aqu'species, density for min'species
    !
    !----------------------------------------------read function coeff's
    SELECT CASE(TRIM(CodFitting))
      CASE("FIXED","NONE")      ;  M%iFitting= 0
      CASE("PHREEQC")           ;  M%iFitting= 1
      CASE("CHRISTOV")          ;  M%iFitting= 2
      CASE DEFAULT
        CALL Stop_(TRIM(W)//" invalid Fitting mode")
    END SELECT

    SELECT CASE(M%iFitting)
    CASE(0)
      M%vX(1)= vX(1)
    CASE(1)
      IF(M%Typ=="AQU") THEN  ;  M%vX(1:6)=  vX(1:6)
      ELSE                   ;  M%vX(1:6)= -vX(1:6)
      ! for PHREEQC or THERMODDEM analytical,
      ! logK of MINERALS is for DISSOCIATION
      ENDIF
    CASE(2)
      M%vX(1:6)= vX(1:6)
    END SELECT
    !---------------------------------------------/read function coeff's
    !
    !WRITE(fff,'(A,A1)') TRIM(L) !-> WRITE the logK's
    !
    !M%Fitting= TRIM(CodFitting)
    !
    IF(iDebug>2) CALL Test(ff,M) !---------------------------------trace
    !
    !---------------------------------------------"fill" the linked list
    N=N+1
    IF(N==1) THEN
      ALLOCATE(LisLogKAnl)
      NULLIFY(LisLogKAnl%next)
      LisLogKAnl%Value= M
      LisCur=> LisLogKAnl
    ELSE
      ALLOCATE(LisCur%next)
      NULLIFY(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=>LisCur%next
    ENDIF
    !--------------------------------------------/"fill" the linked list
    !
  ENDDO DoFile
  !
  DEALLOCATE(vStoik)
  DEALLOCATE(vElement)
  IF(fFormula>0) CLOSE(fFormula)
  !
  IF(iDebug>2) CLOSE(ff)!------------------------------------------trace
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DtbLogKAnl_Read"
  !
  IF(N==0) CALL Stop_("NO SPECIES FOUND ....")
  !
ENDSUBROUTINE DtbLogKAnl_Read

SUBROUTINE FieldList_Read( &
& L,          &
& vStrField, &
& vIField)
  
  USE M_IoTools
  
  !---------------------------------------------------------------------
  CHARACTER(LEN=*),INTENT(INOUT):: L
  CHARACTER(LEN=12),INTENT(IN):: vStrField(:)
  INTEGER,INTENT(OUT):: vIField(:)
  !---------------------------------------------------------------------
  CHARACTER(LEN=80):: W
  LOGICAL:: EoL
  INTEGER:: I,J
  !---------------------------------------------------------------------
  
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

  IF(iDebug==4) THEN
    DO I=1,SIZE(vStrField)
      PRINT *,vIField(I),TRIM(vStrField(I))
    ENDDO
  ENDIF
  !CALL PAUSE_
  !
  IF(vIField(10)==0) & ! for "PARAMETERS"
  CALL Stop_( &
  & "in FieldList_Read: keyword not found for "//TRIM(vStrField(10)))

  IF(vIField(3)==0) & ! for "NAME"
  CALL Stop_( &
  & "in FieldList_Read: keyword not found for "//TRIM(vStrField(3)))

  IF(vIField(4)==0 .AND. vIField(5)==0) & ! for ECFORM/SCFORM
  CALL Stop_( &
  & "in FieldList_Read: keyword not found for "//TRIM(vStrField(4))//"_"//TRIM(vStrField(5)))
  !
  !PRINT *,vIField(1),vIField(2),vIField(3),vIField(4),vIField(5),vIField(9),vIField(10)
  !CALL PAUSE_
END SUBROUTINE FieldList_Read

SUBROUTINE Test(ff,M)
  USE M_T_Species
  USE M_Numeric_Const,ONLY: Ln10
  USE M_T_DtbLogKAnl
  INTEGER,INTENT(IN):: ff
  TYPE(T_DtbLogKAnl),INTENT(IN) :: M
  !
  TYPE(T_Species):: S
  REAL(dp):: vT(1:5),vLogK(1:5),Pbar
  INTEGER :: I
  vT(1:5)= (/0.01,25.,50.,75.,100./)
  Pbar= 1.0D0
  DO I=1,5
    CALL DtbLogKAnl_Calc(M,vT(I)+273.15D0,Pbar,S)
    vLogK(I)= - S%G0RT /Ln10
  ENDDO
  WRITE(ff,'(2(A,1X))',ADVANCE="NO") M%Typ,M%Name
  WRITE(ff,'(5(G15.6,1X))',ADVANCE="NO") (vLogK(I),I=1,5)
  IF(M%V0R>Zero) THEN
    WRITE(ff,'(G15.6,1X)',   ADVANCE="NO") M%WeitKg /M%V0R
  ELSE
    WRITE(ff,'(A,1X)',   ADVANCE="NO") "_"
  ENDIF
  WRITE(ff,'(A,1X)') M%Formula
END SUBROUTINE Test

ENDMODULE M_Dtb_Read_DtbLogKAnl







