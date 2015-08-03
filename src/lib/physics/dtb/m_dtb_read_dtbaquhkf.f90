MODULE M_Dtb_Read_DtbAquHkf
  USE M_Kinds
  USE M_IoTools
  USE M_Trace
  USE M_Dtb_Read_Vars
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: DtbAquHKF_Read
  PUBLIC:: DtbAquHKF_Build
  !
CONTAINS

SUBROUTINE DtbAquHKF_Build
!--
!-- transfer linked list to array of T_DtbAquHkf
!--
  USE M_Dtb_Vars,ONLY: vDtbAquHkf
  !
  TYPE(T_LisAquHkf),POINTER:: LisAquCur, LisAquPrev
  INTEGER:: I
  !
  IF(iDebug>2) PRINT '(A)',"< DtbAquHKF_Build"
  !
  IF(ALLOCATED(vDtbAquHkf)) DEALLOCATE(vDtbAquHkf)
  ALLOCATE(vDtbAquHkf(nAquHkf))
  !
  I=0
  LisAquCur=>LisAquHkf
  DO WHILE (ASSOCIATED(LisAquCur))
    I=I+1
    vDtbAquHkf(I)=LisAquCur%Value
    IF(iDebug>5) PRINT '(2A)', vDtbAquHkf(I)%Name,TRIM(vDtbAquHkf(I)%Num)
    LisAquPrev=>LisAquCur
    LisAquCur=> LisAquCur%next
    DEALLOCATE(LisAquPrev)
  ENDDO
  !
  IF(iDebug>2) PRINT '(A)',"</ DtbAquHKF_Build"
  !
  IF(iDebug>5) CALL Pause_
  
  RETURN
ENDSUBROUTINE DtbAquHKF_Build

SUBROUTINE DtbAquHKF_Read(F,vEle,N)
!--
!-- reads data directly from modified SLOP98.DAT
!-- (the newest version of SPRONS.DAT,the text-file for SupCrt)
!-- ("normally", SUPCRT92.FOR reads from DPRONS,
!-- a "direct-access file" that is built from SLOP98.DAT using CPRONS92.FOR)
!--
  USE M_Dtb_Const,  ONLY: T_CK,Tref,Pref,S0_Hydrogen
  USE M_Files,      ONLY: DirDtbLog,Files_Index_Write
  USE M_T_Element,  ONLY: T_Element,Formula_Build,Formula_Read
  USE M_Dtb_Read_Tools
  !
  USE M_Dtb_Vars,   ONLY: DtbFormat
  
  !---------------------------------------------------------------------
  INTEGER,        INTENT(IN)   :: F !input file
  TYPE(T_Element),INTENT(IN)   :: vEle(:)
  INTEGER,        INTENT(INOUT):: N
  !---------------------------------------------------------------------
  TYPE(T_LisAquHkf),POINTER,SAVE:: LisCur
  TYPE(T_DtbAquHkf):: M
  !
  CHARACTER(LEN=255):: L,W,sFormul
  !
  CHARACTER(LEN=4)  :: ICode
  LOGICAL :: EoL,fOk,EcformIsOk
  INTEGER :: ios,K,ZSp,Div_,iEl
  INTEGER :: I,FF
  REAL(dp):: X
  REAL(dp),DIMENSION(dimV)::vX
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
  INTEGER,ALLOCATABLE::vStoik(:) !for Formula_Read
  !
  !--- for formula translation --
  CHARACTER(LEN=6):: CodFormula
  ! CodFormula is either "SCFORM" or "ECFORM" (names inherited from SUPCRT files):
  !   SCFORM = compact formula,  e.g. SiO2
  !   ECFORM = extended formula, e.g. SI(1)O(2)
  CHARACTER(LEN=2),ALLOCATABLE:: vElement(:)  !
  INTEGER:: fFormula
  !---/
  !---------------------------------------------------------------------
  ! IF (N>1) LisCur => LisAquHkf
  !
  IF(iDebug>0) WRITE(fTrc,"(/,A)") "< DtbReadAquHKF"
  !
  !--------------------------------------------------------------- trace
  FF= 0
  IF(iDebug>2) THEN
    WRITE(fTrc,'(A)') &
    & "stoikio in "//TRIM(DirDtbLog)//"aqu_hkf_stoik.log"
    !
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirDtbLog)//"aquhkf_stoik.log",&
    & "check species stoikio of aqu.species from hkf-type database")
    !
    CALL GetUnit(FF)
    OPEN(FF,FILE=TRIM(DirDtbLog)//"aquhkf_stoik.log")
    !
    WRITE(FF,"(A15,A1)",ADVANCE="NO") "NAME",T_
    DO iEl=1,SIZE(vEle)
      WRITE(FF,"(A3,A1)",ADVANCE="NO") vEle(iEl)%NamEl,T_
    ENDDO
    WRITE(FF,"(A7,A1,A7)") "Div",T_,"FORMULA"
    !
  ENDIF
  !--------------------------------------------------------------/ trace
  !
  !--- for formula translation --
  CodFormula= "ECFORM"
  ALLOCATE(vElement(SIZE(vEle)+3))
  CALL DtbRead_Build_vElement(vEle,vElement)
  fFormula= 0
  !---/
  !
  ALLOCATE(vStoik(1:SIZE(vEle)))
  !
  FilCode= TRIM(DtbFormat)
  !
  !--------------------------------- initialize vStrField(:), vIField(:)
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
  !
  !--- scan default field list
  IF (FilCode(1:5)=="OBIGT") THEN
    L= "INDEX NAME SKIP SCFORM TYPE SKIP SKIP SKIP PARAMETERS"
  ELSE ! case of SLOP98.DAT
    L= "INDEX NAME SKIP SKIP ECFORM SKIP SKIP PARAMETERS"
  ENDIF
  !
  IF(iDebug==4) PRINT *,"< Default values"
  CALL FieldList_Read(L,vStrField,vIField)
  IF(iDebug==4) PRINT *,"</ Default values"
  !---/ scan default field list
  !
  IF(vIField(4)/=0 .AND. fFormula==0 .AND. iDebug>2) THEN
  ! -> files contains compact formulas
    CALL GetUnit(fFormula)
    OPEN(fFormula,file="debug_formula.log")
    WRITE(fFormula,'(A,/)') "resuts of formula conversion"
  ENDIF
  !--------------------------------/ initialize vStrField(:), vIField(:)
  !
  !---------------------------------- build a linked list of all species
  !-------------------------------- consistent with current element list
  DoFile: DO
    !
    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=="!") CYCLE DoFile
    CALL AppendToEnd(L,W,EoL)
    !
    !---------------------------------------- read header, if present --
    !------ if the line begins with any member of array vStrField(:), --
    !------------------------------ then it contains the line headers --
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
      IF(iDebug==4) PRINT *,"< values from file"
      CALL FieldList_Read(L,vStrField,vIField)
      IF(iDebug==4) PRINT *,"</ values from file"

      IF(vIField(4)/=0 .AND. fFormula==0 .AND. iDebug>2) THEN
        CALL GetUnit(fFormula)
        OPEN(fFormula,file="debug_formula.log")
        WRITE(fFormula,'(A,/)') "results of formula conversion"
      ENDIF

      CYCLE DoFile
    ENDIF
    !------------------------------------------/ read header, if present
    !
    !---------------------------------------- process first word of line
    SELECT CASE(W)
    !
    CASE("ENDINPUT")
      EXIT DoFile
    !
    CASE("END","ENDSPECIES")
      EXIT DoFile
    !
    !------------------------------------------------ for old formats --
    CASE("CODE")
      !-> can change the CODE, i.e. the DATA source inside the block
      CALL LinToWrd(L,W,EoL)
      FilCode=TRIM(W)
      !
      ! if FilCode changes, must re-initialize vIField ...!
      IF (FilCode(1:5)=="OBIGT") THEN
        L= "INDEX NAME SKIP SCFORM TYPE SKIP SKIP SKIP PARAMETERS"
      ELSE ! case of SLOP98.DAT
        L= "INDEX NAME SKIP SKIP ECFORM SKIP SKIP PARAMETERS"
      ENDIF
      CALL FieldList_Read(L,vStrField,vIField)
      !
      CYCLE DoFile

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
    !--------------------------------------------------/ for old formats
    CASE DEFAULT
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= TRIM(W)//" "//TRIM(L)
      !
    ENDSELECT
    !---------------------------------------/ process first word of line
    !
    M%Num=     "0"
    M%Name=    "0"
    M%Formula= "0"
    !-------------------------------------------------- scan a data line
    DO I= 1,vIField(10)-1
    ! read words up to position before the parameters
      !
      CALL LinToWrd(L,W,EoL,"NO")
      !
      ! (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
      !   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
      !
      IF(I==vIField(2)) THEN  ! INDEX
        CALL Str_Upper(W)  ;  M%Num= TRIM(W)
      END IF

      IF(I==vIField(3)) THEN  ! NAME
        CALL Str_Upper(W)  ;  M%Name= TRIM(W)
      ENDIF

      IF(I==vIField(4)) THEN ! SCFORM
        !
        CALL DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
        IF(.NOT.EcformIsOk) CYCLE DoFile !-------------------------CYCLE
        !
        CALL Str_Upper(W)  ;  M%Formula=TRIM(W)
      ENDIF

      IF(I==vIField(5)) THEN ! ECFORM
        CALL Str_Upper(W)  ;  M%Formula=  TRIM(W)
      END IF

      !IF(I==vIField(9)) THEN
      !  CALL Str_Upper(W)  ;  CodFitting= TRIM(W)
      !END IF

      !~ IF(I==vIField(6)) CALL WrdToReal(W,X) ! SIZE
      !~ IF(I==vIField(7)) CALL WrdToReal(W,X) ! VOLUME
      !~ IF(I==vIField(8)) CALL WrdToReal(W,X) ! DENSITY

    END DO
    !
    CALL ReadRValsV(L,K,vX)
    !-------------------------------------------------/ scan a data line
    !!
    !! SELECT CASE(W)
    !!   CASE("ENDINPUT"); EXIT DoFile !DoFile
    !!   CASE("END","ENDSPECIES"); EXIT DoFile !DoFile
    !!   !!CASE("END"); EXIT DoFile
    !!   CASE("FORMULA")
    !!     CALL LinToWrd(L,W,EoL)
    !!     CodFormula= TRIM(W)
    !!     IF(iDebug>2 .and. CodFormula=="SCFORM") THEN
    !!       CALL GetUnit(fFormula)
    !!       OPEN(fFormula,file="debug_formula.log")
    !!       WRITE(fFormula,'(A,/)') "resuts of formula conversion"
    !!     ENDIF
    !!     CYCLE DoFile
    !! END SELECT
    !! !
    !! M%Num=TRIM(W)
    !! !
    !! IF (FilCode(1:5)=="OBIGT") THEN
    !!   CodFormula="SCFORM"
    !!   !
    !!   CALL LinToWrd(L,W,EoL); M%Name=TRIM(W) !IF(iDebug>0) WRITE(fTrc,"(A)") M%Name
    !!   CALL LinToWrd(L,W,EoL) !skip ABBRV
    !!   CALL LinToWrd(L,W,EoL,"NO") !-> compact formula, CHARACTER CASE is conserved :!!
    !!   !
    !!   IF(CodFormula=="SCFORM") THEN
    !!     CALL DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
    !!     IF(.not.EcformIsOk) CYCLE DoFile
    !!   ENDIF
    !!   !
    !!   CALL Str_Upper(W)
    !!   !
    !!   M%Formula=TRIM(W)
    !!   !
    !!   CALL LinToWrd(L,W,EoL) !skip STATE
    !!   CALL LinToWrd(L,W,EoL) !skip SOURCE1
    !!   CALL LinToWrd(L,W,EoL) !skip SOURCE2
    !!   CALL LinToWrd(L,W,EoL) !DATE
    !!   !
    !!   CALL ReadRValsV(L,K,vX)
    !!   !
    !! ELSE ! default FilCode is SLOP98
    !!   !
    !!   CALL LinToWrd(L,W,EoL) !skip ABREV
    !!   CALL LinToWrd(L,W,EoL) !skip NAME
    !!   !
    !!   CALL LinToWrd(L,W,EoL); M%Name=TRIM(W) !IF(iDebug>0) WRITE(fTrc,"(A)") M%Name
    !!   CALL LinToWrd(L,W,EoL); M%Formula=TRIM(W) !IF(iDebug>0) WRITE(fTrc,"(A)") M%Formula
    !!   CALL LinToWrd(L,W,EoL) !DATE
    !!   CALL LinToWrd(L,W,EoL) !REF
    !!   !
    !!   CALL ReadRValsV(L,K,vX)
    !!   !
    !! ENDIF
    !
    !------------------------------------------------------ read formula
    CALL Formula_Read(M%Formula,vEle,ZSp,Div_,fOk,vStoik)
    IF(.NOT. fOk) CYCLE DoFile !-----------------------------------CYCLE
    !-----------------------------------------------------/ read formula
    !
    IF (FilCode(1:5)=="OBIGT") THEN
      ! in OBIGT database, Cp and V are tabulated fot aqu'species
      ! but not used as input parameters -> vX(4:5) ignored

      M%G0R=vX(1);   M%H0R=vX(2); M%S0_= vX(3)
      M%A1= vX(6);   M%A2= vX(7);  M%A3=  vX(8); M%A4=vX(9);
      M%C1= vX(10);  M%C2= vX(11)
      M%wref=vX(12)
      M%Chg= INT(vX(13))

    ELSE

      M%G0R=vX(1);  M%H0R=vX(2); M%S0_= vX(3)
      M%A1= vX(4);  M%A2= vX(5); M%A3=  vX(6); M%A4=vX(7);
      M%C1= vX(8);  M%C2= vX(9)
      M%wref= vX(10)
      M%Chg=  INT(vX(11))

    ENDIF
    !
    M%A1=    M%A1*1.0D-1
    M%A2=    M%A2*1.0D02
    M%A4=    M%A4*1.0D04
    M%C2=    M%C2*1.0D04
    M%Wref=  M%wref*1.0D05
    !
    N=N+1
    CALL IntToStr4(N,ICode)
    M%Num= TRIM(FilCode)//"_"//TRIM(ICode)

    M%S0Ele= DOT_PRODUCT(vStoik(:),vEle(:)%S0) - ZSp *S0_Hydrogen ! is in Joule
    M%S0Ele= M%S0Ele /REAL(Div_)

    M%WeitKg= DOT_PRODUCT(vStoik(:), vEle(:)%WeitKg) /REAL(Div_) !-> in Kg
    !
    !-------------------------------------------- "fill" the linked list
    IF(N==1) THEN
      ALLOCATE(LisAquHkf)
      NULLIFY(LisAquHkf%next)
      LisAquHkf%Value= M
      LisCur=> LisAquHkf
    ELSE
      ALLOCATE(LisCur%next)
      NULLIFY(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=>LisCur%next
    ENDIF
    !-------------------------------------------/ "fill" the linked list
    !
    IF(FF>0) THEN
      WRITE(FF,"(A15,A1)",ADVANCE="NO") M%Name,T_
      DO iEl=1,SIZE(vStoik)
        WRITE(FF,"(I3,A1)",ADVANCE="NO") vStoik(iEl),T_
      ENDDO
      X= M%H0R -Tref*M%S0_
      WRITE(FF,'(3(G15.6,A1))',ADVANCE="NO") &
      & M%G0R,T_, X,T_, (M%G0R-X)/Tref,T_
      !
      ! build new formula, with fixed element order
      CALL Formula_Build(vEle,vStoik,Zsp,Div_,sFormul)
      !!!build new formula (without coeff"s !!), with fixed element order
      !!sFormul=""
      !!DO iEl=1,SIZE(vEle)
      !!  IF(vStoik(iEl)>0) sFormul=TRIM(sFormul)//TRIM(vEle(iEl)%NamEl)//""
      !!ENDDO
      WRITE(FF,"(I3,A1,3(A,A1))") &
      & M%Div,T_,TRIM(M%Formula),T_,TRIM(M%Name),T_,TRIM(sFormul)
      !
    ENDIF
    !
  ENDDO DoFile
  !
  DEALLOCATE(vStoik)
  !
  DEALLOCATE(vElement)
  IF(fFormula>0) CLOSE(fFormula)
  !
  IF(FF>0) CLOSE(FF)
  !
  IF(iDebug>0) WRITE(fTrc,"(A,/)") "</ DtbReadAquHKF"
  !
  RETURN
ENDSUBROUTINE DtbAquHKF_Read

SUBROUTINE FieldList_Read( &
& L,         &
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

  IF(iDebug==4) THEN
    DO I=1,SIZE(vStrField)
      PRINT *,vIField(I),TRIM(vStrField(I))
    ENDDO
  ENDIF
  !PAUSE_
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
  !PAUSE_
END SUBROUTINE FieldList_Read

ENDMODULE M_Dtb_Read_DtbAquHkf
