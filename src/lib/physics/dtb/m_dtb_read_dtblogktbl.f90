MODULE M_Dtb_Read_DtbLogKTbl
  USE M_Kinds
  USE M_Trace
  USE M_Dtb_Read_Vars,ONLY: T_LisLogKTbl,LisLogKTbl,nLogKTbl
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: DtbLogKTbl_Read
  PUBLIC:: DtbLogKTbl_Build
  PUBLIC:: DtbLogK_TPCond_Read
  !
CONTAINS

SUBROUTINE DtbLogKTbl_Build
!--
!-- transfer linked list to array of T_DtbLogKTbl
!--
  USE M_Dtb_Vars,ONLY: vDtbLogKTbl
  !
  TYPE(T_LisLogKTbl),POINTER:: LisCur, LisPrev
  INTEGER:: I
  !
  IF(iDebug==5) PRINT '(A)',"< DtbLogKTbl_Build"
  !
  IF(ALLOCATED(vDtbLogKTbl)) DEALLOCATE(vDtbLogKTbl)
  ALLOCATE(vDtbLogKTbl(nLogKTbl))
  !
  I=0
  LisCur=> LisLogKTbl
  DO WHILE (ASSOCIATED(LisCur))
    I=I+1
    vDtbLogKTbl(I)=LisCur%Value
    !
    IF(iDebug==5) PRINT '(2A)', vDtbLogKTbl(I)%Name,TRIM(vDtbLogKTbl(I)%Num)
    !
    LisPrev=>LisCur
    LisCur=> LisCur%next
    DEALLOCATE(LisPrev)
  ENDDO
  !
  IF(iDebug==5) PRINT '(A)',"< DtbLogKTbl_Build"
  IF(iDebug==5) CALL Pause_
  RETURN
ENDSUBROUTINE DtbLogKTbl_Build

LOGICAL FUNCTION LnkSpc_Found(L,W,M)
!--
!-- find index of string W in linked list
!-- return corresponding species Spc
!--
  USE M_T_DtbLogKTbl, ONLY: T_DtbLogKTbl
  !
  TYPE(T_LisLogKTbl),POINTER    :: L
  CHARACTER(LEN=*),  INTENT(IN) :: W
  TYPE(T_DtbLogKTbl),INTENT(OUT):: M
  !
  TYPE(T_LisLogKTbl),POINTER:: P,pPrev
  INTEGER::I
  !
  P=> L
  I=  0
  LnkSpc_Found=.FALSE.
  DO WHILE (ASSOCIATED(P))
    I=I+1
    IF(TRIM(w)==TRIM(P%Value%Name)) THEN
      LnkSpc_Found=.TRUE.
      M= P%Value
      EXIT;
    ENDIF
    pPrev=> P
    P=>     P%next
  ENDDO
ENDFUNCTION LnkSpc_Found

SUBROUTINE DtbLogKTbl_Read(F,vEle,N)
!--
!-- build linked list from database in logK format
!--
  USE M_Dtb_Const,ONLY: T_CK
  USE M_T_Element,ONLY: T_Element,Formula_Read,Formula_Build,Element_Index
  USE M_Files,    ONLY: DirDtbLog,Files_Index_Write
  USE M_IoTools
  USE M_Dtb_Read_Tools
  !
  USE M_T_DtbLogKTbl, ONLY: T_DtbLogKTbl,DimLogK_Max,DtbLogKTbl_Calc_Init
  USE M_Dtb_Vars,     ONLY: DtbLogK_Dim,DtbLogK_vTPCond
  !
  INTEGER,        INTENT(IN)   :: F !input file
  TYPE(T_Element),INTENT(IN)   :: vEle(:)
  INTEGER,        INTENT(INOUT):: N
  !
  TYPE(T_LisLogKTbl),POINTER,SAVE:: LisCur
  !
  CHARACTER(LEN=512):: L,W
  TYPE(T_DtbLogKTbl) :: M
  REAL(dp),DIMENSION(dimV):: vX
  REAL(dp):: X,Rho
  LOGICAL :: EoL, fOk, EcformIsOk
  INTEGER :: I,ff
  INTEGER :: ios,ZSp,Div,mDum
  !
  INTEGER,DIMENSION(:),ALLOCATABLE::vStoik !for Formula_Read
  !
  !--for formula translation
  CHARACTER(LEN=6):: CodFormula
  ! CodFormula is either "SCFORM" or "ECFORM" (names inherited from SUPCRT files):
  !   SCFORM = compact formula,  e.g. SiO2
  !   ECFORM = extended formula, e.g. SI(1)O(2)
  CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE:: vElement  !
  INTEGER:: fFormula
  !--/
  !
  !--------------------------------------------var for header processing
  LOGICAL :: IsHeader
  INTEGER,PARAMETER:: nField= 13
  CHARACTER(LEN=12):: vStrField(1:nField)
  ! vStrField contains names of all possible fields in a database
  INTEGER :: vIField(1:nField)
  !
  ! CHARACTER(LEN=12):: vStrUnit(1:nField)
  !-------------------------------------------/var for header processing
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DtbLogKTbl_Read"
  !
  !----------------------------------------------------------------trace
  IF(iDebug>2) THEN
    CALL GetUnit(ff)
    OPEN(ff,FILE="dtb_logktbl_check.tab")
  ENDIF
  !---------------------------------------------------------------/trace
  !
  fFormula= 0
  CodFormula= "ECFORM"
  ALLOCATE(vElement(SIZE(vEle)+3))
  CALL DtbRead_Build_vElement(vEle,vElement)
  !
  ALLOCATE(vStoik(1:SIZE(vEle)))
  !
  !--initialize vStrField(:), vIField(:)
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  ", &
  &   "SIZE        ","VOLUME      ","DENSITY     " /)
  !
  L= "TYPE INDEX NAME ECFORM SIZE PARAMETERS" != default field list
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
  !-----------------------------------build a linked list of all species
  !---------------------------------consistent with current element list
  DoFile: DO
    !
    READ(f,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile
    CALL AppendToEnd(L,W,EoL) ! if W=="END" append next word
    !
    !--------------------------------------------read header, if present
    !----------if the line begins with any member of array vStrField(:),
    !----------------------------------then it contains the line headers
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
      IF(iDebug>2) PRINT *,"< values from file"
      CALL FieldList_Read(L,vStrField,vIField)
      IF(iDebug>2) PRINT *,"</ values from file"

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
    CASE("END","ENDSPECIES")
      EXIT DoFile
    !
    !----------------------------------------------------for old formats
    CASE("FORMULA")
      CALL LinToWrd(L,W,EoL)
      CodFormula= TRIM(W)

      IF(CodFormula=="SCFORM") THEN
        IF(vIField(4)==0) THEN
          vIField(4)= vIField(5)
          vIField(5)= 0
        ENDIF
        IF(iDebug>2) THEN
          CALL GetUnit(fFormula)
          OPEN(fFormula,FILE="debug_formula.log")
          WRITE(fFormula,'(A,/)') "resuts of formula conversion"
        ENDIF
      ENDIF

      CYCLE DoFile !-----------------------------------------------CYCLE
    !
    !---------------------------------------------/ for old formats --
    !old! CASE DEFAULT
    !old!   !--- ignore line beginning with unknown keyword !!! ??? --
    !old!   CYCLE DoFile !-----------------------------------------CYCLE
    !old! !
    CASE DEFAULT
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= TRIM(W)//" "//TRIM(L)
    !
    ENDSELECT
    !----------------------------------------/process first word of line
    !
    M%Num= "0"
    !---------------------------------------------------scan a data line
    DO I= 1,vIField(10)-1
    ! read words up to position before the parameters
      !
      CALL LinToWrd(L,W,EoL,"NO")
      !
      !& (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
      !&   "SIZE        ","VOLUME      ","DENSITY     ","FITTING     ","PARAMETERS  ", &
      !&   "SKIP        ","SOURCE      "/)
      !
      IF(I==vIField(1)) THEN  ! TYPE
        !
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
        !
      END IF

      IF(I==vIField(2)) THEN  ! INDEX
        CALL Str_Upper(W)  ;  M%Num= TRIM(W)
      END IF

      IF(I==vIField(3)) THEN  ! NAME
        CALL Str_Upper(W)  ;  M%Name= TRIM(W)
      ENDIF

      IF(I==vIField(4)) THEN  ! SCFORM
        !
        CALL DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
        IF(.NOT.EcformIsOk) CYCLE DoFile !-------------------------CYCLE
        !
        CALL Str_Upper(W)  ;  M%Formula=TRIM(W)
      ENDIF

      IF(I==vIField(5)) THEN  ! ECFORM
        CALL Str_Upper(W)  ;  M%Formula=  TRIM(W)
      END IF

      !IF(I==vIField(9)) THEN
      !  CALL Str_Upper(W)  ;  CodFitting= TRIM(W)
      !END IF

      IF(I==vIField(11)) CALL WrdToReal(W,X) ! SIZE
      IF(I==vIField(12)) CALL WrdToReal(W,X) ! VOLUME
      IF(I==vIField(13)) CALL WrdToReal(W,X) ! DENSITY

    END DO
    !
    CALL ReadRValsV(L,mDum,vX)
    !--------------------------------------------------/scan a data line
    !
    IF(mDum<DtbLogK_Dim) &
    & CALL Stop_("Dimension of LogK array < Dimension of TP table !!")
    !
    !-------------------------------------------------------read formula
    CALL Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
    IF(.NOT. fOk) CYCLE DoFile !-----------------------------------CYCLE
    !------------------------------------------------------/read formula
    !
    !-- compute Molecular Weit
    M%WeitKg= DOT_PRODUCT(vStoik(:),vEle(:)%WeitKg) /REAL(Div)
    M%Div= Div
    M%Chg= Zsp
    !
    !---------size parameter for aqu'species, or density for min'species
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
    !--default value of size parameter for charged aqu'species
    IF(M%Typ=="AQU" .and. M%Chg/=0 .and. M%AquSize<=Zero) M%AquSize= 3.72D0
    !----/read size parameter for aqu'species, or volume for min'species
    !
    !WRITE(fff,'(A,A1)') TRIM(L) !-> WRITE the logK's
    !
    !--- read list of logK's --
    M%DimLogK=  DtbLogK_Dim
    M%vLogK(1:DimLogK_Max)= vX(1:DimLogK_Max)
    !---/ read list of logK's --
    !
    M%vTdgK(1:DtbLogK_Dim)= DtbLogK_vTPCond(1:DtbLogK_Dim)%TdgC + T_CK
    !
    CALL DtbLogKTbl_Calc_Init(M)
    !
    !! IF(iDebug==5) WRITE(ff,'(A3,A1,A15,A1,A39,A1,I3)') &
    !! & M%Typ,T_,M%Name,T_,M%Formula,T_,N
    !
    IF(iDebug>2) CALL Test(ff,M) !--trace
    !
    !---------------------------------------------"fill" the linked list
    N=N+1
    IF(N==1) THEN
      ALLOCATE(LisLogKTbl)
      NULLIFY(LisLogKTbl%next)
      LisLogKTbl%Value= M
      LisCur=> LisLogKTbl
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
  !
  DEALLOCATE(vElement)
  IF(fFormula>0) CLOSE(fFormula)
  !
  IF(iDebug>2) CLOSE(ff)!--trace
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "< DtbLogKTbl_Read"
  IF(N==0) CALL Stop_("NO SPECIES FOUND ....")
  !
ENDSUBROUTINE DtbLogKTbl_Read

SUBROUTINE FieldList_Read( &
& L,          &
& vStrField, &
& vIField)
  USE M_IoTools
  !
  CHARACTER(LEN=*),INTENT(INOUT):: L
  CHARACTER(LEN=*),INTENT(IN):: vStrField(:)
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

  IF(iDebug==4) THEN
    DO I=1,SIZE(vStrField)
      PRINT *,vIField(I),TRIM(vStrField(I))
    ENDDO
  ENDIF
  !CALL PAUSE_
  !
  !& (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  !&   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  ", &
  !&   "SIZE        ","VOLUME      ","DENSITY     " /)
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
  USE M_T_DtbLogKTbl
  
  INTEGER,INTENT(IN):: ff
  TYPE(T_DtbLogKtbl),INTENT(IN) :: M
  
  TYPE(T_Species):: S
  REAL(dp):: vT(1:5),vLogK(1:5),Pbar
  INTEGER :: I
  
  vT(1:5)= (/0.01,25.,50.,75.,100./)
  Pbar= 1.0D0
  DO I=1,5
    CALL DtbLogKtbl_Calc(M,vT(I)+273.15D0,Pbar,S)
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

SUBROUTINE DtbLogK_TPCond_Read(N)
!--
!-- reads block TP.TABLE, build DtbLogK_vTPCond
!--
  USE M_IoTools !,ONLY:GetUnit,dimV
  USE M_Files,     ONLY: NamFInn
  USE M_Dtb_Const, ONLY: T_CK
  USE M_T_Tpcond,  ONLY: T_TPCond
  USE M_Fluid_Calc,ONLY: Eos_H2O_psat
  !
  USE M_Dtb_Vars,    ONLY: DtbLogK_vTPCond,DtbLogK_Dim
  USE M_T_DtbLogKTbl,ONLY: DimLogK_Max
  !
  INTEGER,INTENT(OUT):: N
  !
  LOGICAL :: EoL,Ok
  INTEGER :: i,mDum,ios,f
  CHARACTER(LEN=512):: L,W,W1
  !
  REAL(dp):: TdgK,Pbar,PSatBar
  REAL(dp):: vX(dimV)
  TYPE(T_TPCond):: vCond(dimV)
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< DtbLogK_TPCond_Read"
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFInn))
  !
  N= 0
  Ok=.FALSE.
  !
  vCond(:)%TdgC= Zero
  vCond(:)%Pbar= Zero
  !
  DoFile: DO
    !
    READ(f,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!')   CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    !
    SELECT CASE(W)
    !
    CASE("ENDINPUT"); EXIT  DoFile
    !
    CASE("TP.TABLE","TPTABLE")
      Ok=  .TRUE.
      N= dimV
      !
      !---------------------------------------------- read table loop --
      DoTPtable: DO
        !
        READ(f,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,Eol)
        IF(W(1:1)=='!') CYCLE DoTPtable
        CALL AppendToEnd(L,W,EoL)
        !
        SELECT CASE(W)
        !
        CASE("ENDINPUT"); EXIT DoFile
        CASE("END","ENDTPTABLE","ENDTP.TABLE"); EXIT DoTPtable
        !
        CASE("TDGC","TDGK", "PBAR","PMPA","TDEGC","TDEGK")
          !
          !!CALL ReadRValsV(L,mDum,vX)
          !
          i=0
          DO
            CALL LinToWrd(L,W1,EoL)
            i=i+1
            IF(i>DimV) EXIT
            IF(TRIM(W1)=="PSAT") THEN
              IF(vCond(i)%TdgC>Zero) THEN
                CALL Eos_H2O_psat(vCond(i)%TdgC+T_CK,vX(i))
                !PRINT *,"TdgC= ", vCond(i)%TdgC
                !PRINT *,"psat= ", vX(i)  ;  CALL Pause_
              ELSE
                CALL Stop_("error in TP.TABLE: invalid TdgC for Psat computation")
              ENDIF
            ELSE
              CALL WrdToReal(W1,vX(i))
            ENDIF
            IF(EoL) EXIT
          ENDDO
          !
          mDum= i
          N=min(N,mDum)
          !
          IF(iDebug>0) WRITE(fTrc,'(A,A1,A,2I3)') TRIM(W),T_," DIM= ", mDum, N
          !
          SELECT CASE(TRIM(W))
          CASE("TDGC");  vCond(:)%TdgC=vX(:)
          CASE("TDGK");  vCond(:)%TdgC=vX(:)-T_CK
          !
          CASE("PBAR");  vCond(:)%Pbar=vX(:)       !pressure in Bar
          CASE("PMPA");  vCond(:)%Pbar=vX(:)*1.D-1 !MegaPascal to Bar
          !
          CASE("TDEGC"); vCond(:)%TdgC=vX(:)       ; CALL Warning_("Depreciated Keyword TDEGC")
          CASE("TDEGK"); vCond(:)%TdgC=vX(:)-T_CK  ; CALL Warning_("Depreciated Keyword TDEGK")
          !
          END SELECT
        !ENDCASE("TDGC",..
        !
        CASE("NAME")
          i=0
          DO
            CALL LinToWrd(L,W,EoL)
            i=i+1
            IF(I>DimV) EXIT
            vCond(i)%Name=TRIM(W)
            IF(iDebug>0) WRITE(fTrc,'(I3,1X,A15)') i,vCond(i)%Name
            IF(EoL) EXIT
          ENDDO
        !ENDCASE("NAME")
        !
        CASE DEFAULT; CALL Stop_(TRIM(W)//"<<unknown Keyword...")
        !
        END SELECT !ENDIF
        !ENDIF
      ENDDO DoTPtable
      !---------------------------------------------/ read table loop --
    !END CASE("ENDTP.TABLE")
    END SELECT
    !
  ENDDO DoFile
  !
  CLOSE(f)
  !
  IF(Ok) THEN
    !! IF(any(vCond(1:N)%TdgC<Zero)) CALL Stop_("all temperatures should be >0 !!!")
    IF(ANY(vCond(1:N)%Pbar<1.D-9)) CALL Stop_("all pressures should be >0 !!!")
    !
    DO I=1,N
      TdgK= vCond(I)%TdgC +T_CK
      Pbar= vCond(I)%Pbar
      IF(TdgK <Zero) TdgK= 200.0D0 ! Zero
      !
      !-- check that pressure is not below
      !-- the vapor saturation curve for pure water
      !-- critical point H2O- 647.25D0 TdgK / 220.55D0 Pbar
      IF (TdgK<=647.25D0) THEN; CALL Eos_H2O_psat(TdgK,PSatBar)
      ELSE                    ; PSatBar= 220.55D0
      ENDIF
      !-- IF pressure is too low, adjust to saturation pressure at TdgK
      IF(Pbar<PSatBar) Pbar=PSatBar
      !
      vCond(I)%TdgC= TdgK -T_CK
      vCond(I)%Pbar= Pbar
    ENDDO
    !
    IF(N>DimLogK_Max) N= DimLogK_Max
    !
    DtbLogK_Dim= N
    !
    IF(ALLOCATED(DtbLogK_vTPCond)) DEALLOCATE(DtbLogK_vTPCond)
    ALLOCATE(DtbLogK_vTPCond(N))
    DtbLogK_vTPCond(1:N)= vCond(1:N)
    !
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ DtbLogK_TPCond_Read"
  !
ENDSUBROUTINE DtbLogK_TPCond_Read

ENDMODULE M_Dtb_Read_DtbLogKTbl
