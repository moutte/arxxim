MODULE M_Dtb_Read_DtbMinThr
  USE M_Kinds
  USE M_IoTools
  USE M_Trace
  USE M_Dtb_Read_Vars
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: DtbMinThr_Build
  PUBLIC:: DtbMinThr_Read
  PUBLIC:: DtbMinThr_Read_Line
  PUBLIC:: DtbMinThr_ToLine
  !
  TYPE(T_LisMinThr),POINTER,SAVE:: LisCur, LisPrev

CONTAINS

SUBROUTINE DtbMinTHR_Build
!--
!-- -> build vDtbMinThr
!-- transform LinkedList LisMinThr to mineral table vDtbMinThr
!--
  USE M_T_Element,  ONLY: T_Element,Formula_Read
  USE M_Dtb_Vars,   ONLY: vDtbMinThr
  !
  INTEGER:: I
  !
  IF(iDebug==5) PRINT '(A)',"< DtbMinTHR_Build"
  !
  IF(ALLOCATED(vDtbMinThr)) DEALLOCATE(vDtbMinThr); ALLOCATE(vDtbMinThr(nMinThr))
  LisCur=>LisMinThr
  I=0
  DO WHILE (ASSOCIATED(LisCur))
    I=I+1
    vDtbMinThr(I)=LisCur%Value
    IF(iDebug==5) PRINT '(2A)', vDtbMinThr(I)%Name,TRIM(vDtbMinThr(I)%Num) !; PAUSE
    LisPrev=>LisCur
    LisCur=> LisCur%next
    DEALLOCATE(LisPrev)
  ENDDO
  !
  IF(iDebug==5) PRINT '(A)',"</ DtbMinTHR_Build"
  IF(iDebug==5) CALL Pause_
  !
  RETURN
ENDSUBROUTINE DtbMinTHR_Build

SUBROUTINE DtbMinThr_Read(F,vEle,N)
!--
!-- -> build or Append LisMinThr
!-- reads thermodyn. data from (modified) theriak-compatible datafile
!--
  USE M_T_Element,  ONLY: T_Element,Element_Index
  USE M_T_DtbMinThr,ONLY: DtbMinThr_Zero
  USE M_Files,      ONLY: DirDtbLog,Files_Index_Write
  USE M_Dtb_Read_Tools
  USE M_Dtb_Vars,   ONLY: DtbFormat
  !
  INTEGER,INTENT(IN):: F !input file
  TYPE(T_Element),INTENT(IN):: vEle(:)
  INTEGER,INTENT(INOUT):: N
  !
  CHARACTER(LEN=512):: L,W,Wb
  LOGICAL :: EoL
  LOGICAL :: bMin, bGas
  REAL(dp):: X
  INTEGER :: ios,I,iEl,FF,ieO_
  TYPE(T_DtbMinThr) :: M
  !
  CHARACTER(LEN=6):: CodFormula ! SCFORM | ECFORM
  ! SCFORM = compact formula,  e.g. SiO2
  ! ECFORM = extended formula, e.g. SI(1)O(2)
  CHARACTER(LEN=2),DIMENSION(:),ALLOCATABLE:: vElement
  INTEGER:: fFormula
  LOGICAL:: EcformIsOk
  !
  IF(iDebug>0) WRITE(fTrc,"(/,A)") "< DtbMinThr_Read"
  !
  CodFormula= "ECFORM" ! is default value
  ALLOCATE(vElement(SIZE(vEle)+3))
  CALL DtbRead_Build_vElement(vEle,vElement)
  !
  fFormula= 0
  !
  !------------------------------------------------------------ trace --
  IF(iDebug>2) THEN
    WRITE(fTrc,'(A)') &
    & "stoikio in "//TRIM(DirDtbLog)//"min_thr_stoik.log"
    CALL GetUnit(FF)
    OPEN(FF,FILE=TRIM(DirDtbLog)//"min_thr_stoik.log")
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirDtbLog)//"min_thr_stoik.log",&
    & "check species stoikio of min.species from theriak-type database")
    WRITE(FF,"(A15,A1)",ADVANCE="NO") "NAME",T_
    DO iEl=1,SIZE(vEle)
      WRITE(FF,"(A3,A1)",ADVANCE="NO") vEle(iEl)%NamEl,T_
    ENDDO
    WRITE(FF,"(A7,A1,A7)") "Div",T_,"FORMULA"
  ENDIF
  !-----------------------------------------------------------/ trace --
  !
  FilCode=TRIM(DtbFormat) !-> default value
  ieO_= Element_Index("O__",vEle)
  !
  M%Name= "ZZZ"
  M%Num=  "0"
  !
  DoFile: DO
    
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=="!") CYCLE DoFile
    CALL AppendToEnd(L,W,EoL)
    
    SELECT CASE(W)
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
    CASE("CODE")
      CALL LinToWrd(L,W,EoL)
      FilCode=TRIM(W)
      CYCLE DoFile
    !_
    CASE("ENDINPUT"); EXIT DoFile
    !
    CASE("END","ENDSPECIES"); EXIT DoFile
    !
    CASE("MINERAL","GAS")
    
      SELECT CASE(W)
      CASE("MINERAL") ; bGas=.FALSE.; bMin=.TRUE.
      CASE("GAS")     ; bGas=.TRUE.;  bMin=.FALSE.
      END SELECT
      
      DoReadMin: DO
        
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        
        IF(W(1:1)=="!") CYCLE DoReadMin !skip comment lines
        
        CALL AppendToEnd(L,W,EoL)
        
        SELECT CASE(W)
        CASE("ENDINPUT"); EXIT DoFile
        CASE("ENDSPECIES"); EXIT DoFile
        CASE("END","ENDMINERAL","ENDGAS"); CYCLE DoFile
        END SELECT
        
        IF(W(1:1)/="&") THEN
        !IF first char is not "&", the line may contain name, formula, etc
        
          IF(M%Name/="ZZZ") THEN !before proceeding,  SAVE current M in list
            CALL Save_Record
            CALL DtbMinThr_Zero(M)
          ENDIF
          IF(bGas) M%Typ="GAS"
          IF(bMin) M%Typ="MIN"
          M%Name=TRIM(W)
          
          !--------------------------------processing compact formulas--
          IF(CodFormula=="SCFORM") THEN
            CALL LinToWrd(L,W,EoL,"NO")
            !-> compact formula, CHARACTER CASE is conserved :!!
            CALL DtbRead_Build_ExtendedFormula(fFormula,vElement,W,EcformIsOk)
            IF(.not.EcformIsOk) THEN !CYCLE DoReadMin
              CALL Stop_("!!! Cannot translate "//TRIM(W))
            ENDIF
            CALL Str_Upper(W)
            M%Formula=TRIM(W)
          ELSE
            CALL LinToWrd(L,W,EoL)
            M%Formula= TRIM(W)
          ENDIF
          !--------------------------------------------------------end--
          
          CALL LinToWrd(L,W,EoL) !; M%Abbr=TRIM(W)
          !
          CYCLE DoReadMin
          
        ENDIF
        
        IF(W(1:1)=="&") THEN
        ! if first char is "&",
        ! then the line is the continuation of preceding, i.e. contains data
          
          CALL LinToWrd(L,W,EoL)
          
          SELECT CASE(TRIM(W))
          
          CASE("ST")
            CALL ReadRVals5(L,M%G0R,M%H0R,M%S0_,M%V0R,X)
          
          CASE("C1","CP1")
            CALL ReadRVals5(L,M%K1,M%K4,M%K3,M%K8,X)
            !K1,K4,K3,K8
          CASE("C2","CP2")
            CALL ReadRVals5(L,M%K6,M%K2,M%K5,M%K7,M%K9)
            !Cp= K6/T + K2*T + K5*TT + K7*SQT + K9*TT*T
          CASE("C3","CP3")
              CALL ReadRVals5(L,M%K1,M%K2,M%K3,M%K4,X)
          
          CASE("V1")
            M%codVol=1
            CALL ReadRVals5(L,M%VTA,M%VTB,M%VPA,M%VPB,X)
            M%VTA=M%VTA*M%V0R/1.0D5 ; M%VTB=M%VTB*M%V0R/1.0D5
            M%VPA=M%VPA*M%V0R/1.0D5 ; M%VPB=M%VPB*M%V0R/1.0D8
          CASE("VHP")
            M%codVol=4
            CALL ReadRVals5(L,M%VTA,M%VTB,M%TKRI,M%SMA,M%VPA)
          CASE("VH2")
            M%codVol=4
            CALL ReadRVals5(L,M%D1, M%D2, M%D3,X,X)
          
          CASE("D1")
            M%DIS=.TRUE.
            CALL ReadRVals5(L,M%D1, M%D4, M%D3, M%D8, M%D6)
          CASE("D2")
            CALL ReadRVals5(L,M%D2, M%D5, M%TD0,M%TDMAX,M%VADJ)
          
          CASE("TR1","T1","TR3","T3") !for Berman dtbase"s
            M%nLanda=M%nLanda+1; I=M%nLanda
            CALL ReadRVals5(L,M%TQ1B(I),M%TRE(I),M%ASPK(I),M%BSPK(I),M%DHTR(I))
          CASE("TR2","T2")
            IF (M%nLanda>0) THEN
              I=M%nLanda
              CALL ReadRVals5(L,M%TEQ(I),M%DVTR(I),M%DVDT(I),M%DVDP(I),X)
            ENDIF
          
          CASE("REDKWON","VDWAALS")
          ! Redlich-Kwong, VanDerWaals
            M%CodGas=TRIM(W)
            CALL ReadRVals5(L,M%AA0,M%AAT,M%BB0,M%BBT,X)
          CASE("SOAVE","PENGROB")
          ! Soave Cubic, Peng Robinson cubic
            M%CodGas=TRIM(W)
            CALL ReadRVals5(L,M%TCrit,M%PCrit,M%ACentric,X,X)
            !!WRITE(fTrc,'(A,3G15.6)') M%CodGas,M%TCrit,M%PCrit,M%ACentric
          CASE("PRSV")
          ! 
            M%CodGas="PRSV"
            CALL ReadRVals5(L,M%TCrit,M%PCrit,M%ACentric,X,X)
          
          CASE("SPECIAL")
            CALL LinToWrd(L,Wb,EoL); M%Special= TRIM(Wb)
          
          !CASE("V2"); M%codVol=2; CALL ReadRVals6(L,M%VAA,M%VAB,M%VB,            X,X,X)  !VO2=.TRUE.
          !CASE("V3"); M%codVol=3; CALL ReadRVals6(L,M%VL0,M%VLA,M%VLN, M%VL2,     X,X)   !VO3=.TRUE.
          !CASE("AQ1") K1,K2,K8,K9,K7  AQU=.TRUE. PROVID="AQU"
          !CASE("AQ2") D1,D2,D3,D4
          !CASE("AQP") D1,D2           AQU=.TRUE. PROVID="AQP"
          !CASE("SPC") CASE            SPC=.TRUE.
          !CASE("TL1") THEN; M%TL1=.TRUE.; CALL ReadRVals6(L,M%TKRI,M%SMA,X,X,X,X); ENDIF !!!obsolete ???
          
          END SELECT
          
        ENDIF !(W(1:1)=="&")
      
      ENDDO DoReadMin
      
    !ENDCASE("MINERAL","GAS")
      
    END SELECT
    
  ENDDO DoFile
  !
  IF(M%Name/="ZZZ") CALL Save_Record !Save last record
  !
  DEALLOCATE(vElement)
  IF(fFormula>0) CLOSE(fFormula)
  !
  IF(iDebug>2) CLOSE(FF)
  IF(iDebug>0) WRITE(fTrc,"(A,/)") "</ DtbMinThr_Read"
  
  RETURN
  !
CONTAINS

SUBROUTINE Save_Record
  USE M_T_Element,ONLY: Formula_Read,Formula_Build
  !
  CHARACTER(LEN=4)  :: ICode
  CHARACTER(LEN=512):: sFormul
  INTEGER,DIMENSION(1:SIZE(vEle))::vStoik !for Formula_Read
  LOGICAL :: fOk
  INTEGER :: ZSp,Div
  !
  CALL Formula_Read(M%Formula,vEle,ZSp,Div,fOk,vStoik)
  IF(fOk) THEN
    !
    N=N+1
    CALL IntToStr4(N,ICode)
    M%Num=TRIM(FilCode)//"_"//TRIM(ICode)
    !
    M%S0Ele= DOT_PRODUCT(vStoik(:), vEle(:)%S0) /REAL(Div) !in Joule
    M%WeitKg=DOT_PRODUCT(vStoik(:), vEle(:)%WeitKg) /REAL(Div)
    M%Div= Div
    !
    IF(N==1) THEN
      ALLOCATE(LisMinThr)
      NULLIFY(LisMinThr%next)
      LisMinThr%Value=M
      LisCur=> LisMinThr
    ELSE
      ALLOCATE(LisCur%next)
      NULLIFY(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=> LisCur%next
    ENDIF
    !
    IF(iDebug>2) THEN
      WRITE(FF,"(A15,A1)",ADVANCE="NO") M%Name,T_
      DO iEl=1,SIZE(vEle)
        WRITE(FF,"(I5,A1)",ADVANCE="NO") vStoik(iEl),T_
      ENDDO
      !build new formula, with fixed element order
      CALL Formula_Build(vEle,vStoik,Zsp,M%Div,sFormul)
      !!sFormul=""
      !!DO iEl=1,SIZE(vEle)
      !!  IF(vStoik(iEl)>0) sFormul=TRIM(sFormul)//TRIM(vEle(iEl)%NamEl)//""
      !!ENDDO
      WRITE(FF,"(I5,A1,A39,A1,A15,A1,A)") &
      & M%Div,T_,M%Formula,T_,M%Name,T_,TRIM(sFormul)
    ENDIF
    !
  ENDIF
ENDSUBROUTINE Save_Record

ENDSUBROUTINE DtbMinThr_Read

SUBROUTINE DtbMinThr_Read_Line(F,vEle,N)
!--
!---> build or Append LisMinThr
!--
!-- reads thermodyn. data from theriak-compatible single line datafile
!-- this single line reading is called when DtbFormat is "HSV.THR"
!-- (multiline is called for "MIN.THR" or "GAS.THR")
!-- differences from the multiline format:
!--   1/ whether the species is MIN or GAS is indicated for each species,
!--   and not for a whole block
!--   2/ columns on left are free format
!-----------------------------------------------------------------------
!-- first column- AQU/MIN/GAS
!-- then columns for species description (Name, Formula, Data source, ...),
!-- then one or more blocks of 6 columns,
!-- each with a code (ST,C1,...) on column 1 and data on columns 2..6
!-----------------------------------------------------------------------

  USE M_T_Element,  ONLY: T_Element,Element_Index
  USE M_T_Element,  ONLY: Formula_Read
  USE M_T_DtbMinThr,ONLY: DtbMinThr_Zero
  USE M_Files,      ONLY: DirDtbLog,Files_Index_Write
  USE M_Dtb_Read_Tools
  USE M_Dtb_Vars,   ONLY: DtbFormat
  !
  INTEGER,        INTENT(IN):: F !input file
  TYPE(T_Element),INTENT(IN):: vEle(:)
  !
  INTEGER,        INTENT(INOUT):: N
  !
  CHARACTER(LEN=512):: L,W,Wb
  LOGICAL :: EoL
  REAL(dp):: X
  INTEGER :: ios,I,iEl,FF,ieO_,iLine
  TYPE(T_DtbMinThr) :: M
  !TYPE(T_LisMinThr),POINTER,SAVE:: LisCur
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
  !--- for Formula_Read
  INTEGER :: vStoik(1:SIZE(vEle))
  LOGICAL :: fOk
  INTEGER :: ZSp,Div
  !---/
  !
  !--- for formula translation --
  CHARACTER(LEN=6):: CodFormula
  ! CodFormula is either "SCFORM" or "ECFORM" (names inherited from SUPCRT files):
  !   SCFORM = compact formula,  e.g. SiO2
  !   ECFORM = extended formula, e.g. SI(1)O(2)
  CHARACTER(LEN=2),ALLOCATABLE:: vElement(:)  !
  INTEGER:: fFormula
  LOGICAL:: EcformIsOk
  !---/
  !
  IF(iDebug>0) WRITE(fTrc,"(/,A)") "< DtbMinThr_Read_line"
  !
  !--- for formula translation --
  CodFormula= "ECFORM" ! is default value
  ALLOCATE(vElement(SIZE(vEle)+3))
  CALL DtbRead_Build_vElement(vEle,vElement)
  fFormula= 0
  !---/
  !
  ! IF(N>0) LisCur=> LisMinThr
  !
  !------------------------------------------------------------ trace --
  IF(iDebug<3) THEN
    FF= 0
  ELSE
    !
    WRITE(fTrc,'(A)') &
    & "stoikio in "//TRIM(DirDtbLog)//"min_thr_stoik.log"
    !
    CALL GetUnit(FF)
    OPEN(FF,FILE=TRIM(DirDtbLog)//"min_thr_stoik.log")
    !
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirDtbLog)//"min_thr_stoik.log",&
    & "check species stoikio of min.species from theriak-type database")
    !
    WRITE(FF,"(A15,A1)",ADVANCE="NO") "NAME",T_
    DO iEl=1,SIZE(vEle)
      WRITE(FF,"(A3,A1)",ADVANCE="NO") vEle(iEl)%NamEl,T_
    ENDDO
    WRITE(FF,"(A7,A1,A7)") "Div",T_,"FORMULA"
    !
  ENDIF
  !-----------------------------------------------------------/ trace --
  !
  FilCode= TRIM(DtbFormat) !-> default value
  iLine= 0
  ieO_= Element_Index("O__",vEle)
  !
  !------------------------------ initialize vStrField(:), vIField(:) --
  vStrField(1:nField)= &
  !!!!"____________","____________","____________","____________","____________",
  & (/"TYPE        ","INDEX       ","NAME        ","SCFORM      ","ECFORM      ", &
  &   "SKIP        ","SOURCE      ","FORMAT      ","FITTING     ","PARAMETERS  " /)
  !
  !--- scan default field list
  L= "TYPE NAME ECFORM SKIP SOURCE SKIP PARAMETERS"
  !
  !IF(iDebug==4) PRINT *,"< Default values"
  CALL FieldList_Read(L,vStrField,vIField)
  !IF(iDebug==4) PRINT *,"</ Default values"
  !---/ scan default field list
  !
  IF(vIField(4)/=0 .AND. fFormula==0 .AND. iDebug>2) THEN
  ! -> files contains compact formulas
    CALL GetUnit(fFormula)
    OPEN(fFormula,file="debug_formula.log")
    WRITE(fFormula,'(A,/)') "resuts of formula conversion"
  ENDIF
  !-----------------------------/ initialize vStrField(:), vIField(:) --
  !
  !------------------------------- build a linked list of all species --
  !----------------------------- consistent with current element list --
  DoFile: DO

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
    !---------------------------------------/ read header, if present --
    !
    !------------------------------------- process first word of line --
    SELECT CASE(W)
    !
    CASE("ENDINPUT")
      EXIT DoFile
    !
    CASE("END","ENDSPECIES")
      EXIT DoFile
    !
    CASE DEFAULT
      ! in other cases, the line (may) contain data
      ! -> re-assemble W at beginning of L for further processing
      L= TRIM(W)//" "//TRIM(L)
    !
    ENDSELECT
    !------------------------------------/ process first word of line --
    !
    CALL DtbMinThr_Zero(M)
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
        IF(.NOT.EcformIsOk) CYCLE DoFile !====================< CYCLE ==
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
    !----------------------------------------- loop on 6-cells groups --
    DO
    ! read species data
    ! (= organized as groups of 6 cells:
    !  1 group= 1 code (ST,C1,V1,...) followed by 5 numeric data cells)
    !
      CALL LinToWrd(L,W,EoL)
      IF(EoL) EXIT

      SELECT CASE(TRIM(W))

      CASE("SPECIAL")
        CALL LinToWrd(L,Wb,EoL)
        M%Special= TRIM(Wb)
        
      CASE("ST")
        CALL ReadRVals5(L,M%G0R,M%H0R,M%S0_,M%V0R,X)

      CASE("C1","CP1")
        CALL ReadRVals5(L,M%K1,M%K4,M%K3,M%K8,X)
        ! K1,K4,K3,K8

      CASE("C2","CP2")
        CALL ReadRVals5(L,M%K6,M%K2,M%K5,M%K7,M%K9)
        ! Cp= K6/T + K2*T + K5*TT + K7*SQT + K9*TT*T

      CASE("C3","CP3")
        CALL ReadRVals5(L,M%K1,M%K2,M%K3,M%K4,X)
        ! Cp= K1 + K2*T + K3/TT + K4/SQT -> Maier-Kelley

      CASE("V1")
        M%codVol=1
        CALL ReadRVals5(L,M%VTA,M%VTB,M%VPA,M%VPB,X)
        !
        M%VTA=M%VTA*M%V0R/1.0D5  ;  M%VTB=M%VTB*M%V0R/1.0D5
        M%VPA=M%VPA*M%V0R/1.0D5  ;  M%VPB=M%VPB*M%V0R/1.0D8

      CASE("VHP")
        M%codVol=4
        CALL ReadRVals5(L,M%VTA,M%VTB,M%TKRI,M%SMA,M%VPA)

      CASE("VH2")
        M%codVol=4
        CALL ReadRVals5(L,M%D1, M%D2, M%D3,X,X)

      CASE("D1")
        M%DIS=.TRUE.
        CALL ReadRVals5(L,M%D1, M%D4, M%D3, M%D8,   M%D6)

      CASE("D2")
        CALL ReadRVals5(L,M%D2, M%D5, M%TD0,M%TDMAX,M%VADJ)

      CASE("TR1","T1","TR3","T3") !for Berman database"s
        M%nLanda= M%nLanda +1
        I= M%nLanda
        CALL ReadRVals5(L,M%TQ1B(I),M%TRE(I),M%ASPK(I),M%BSPK(I),M%DHTR(I))

      CASE("TR2","T2")
        IF (M%nLanda>0) THEN
          I=M%nLanda
          CALL ReadRVals5(L,M%TEQ(I),M%DVTR(I),M%DVDT(I),M%DVDP(I),X)
        ENDIF

      CASE("REDKWON","VDWAALS")
        M%CodGas=TRIM(W)
        CALL ReadRVals5(L,M%AA0,M%AAT,M%BB0,M%BBT,X)

      CASE("SOAVE","PENGROB")
        M%CodGas= TRIM(W)
        CALL ReadRVals5(L,M%TCrit,M%PCrit,M%ACentric,X,X)

      CASE("PRSV")
        M%CodGas= "PRSV"
        CALL ReadRVals5(L,M%TCrit,M%PCrit,M%ACentric,X,X)

      CASE DEFAULT
        CALL Stop_(TRIM(W)//"= unknown code in HSV/Theriak format")

      END SELECT

    ENDDO
    !--------------------------------------/ loop on 6-cells groups --
    !
    CALL Save_Record
    !
  ENDDO DoFile
  !
  !! PRINT *, "nMinThr=", N  ;  PAUSE
  !
  DEALLOCATE(vElement)
  IF(fFormula>0) CLOSE(fFormula)
  !
  IF(FF>0) CLOSE(FF)
  IF(iDebug>0) WRITE(fTrc,"(A,/)") "</ DtbMinThr_Read"
  !
  RETURN
  !
CONTAINS

  SUBROUTINE Save_Record
  !--
  !-- save record in the linked list LisMinThr
  !--
    USE M_T_Element,ONLY: Formula_Build
    !
    CHARACTER(LEN=4)  :: ICode
    CHARACTER(LEN=512):: sFormul
    !
    N=N+1
    !
    iLine= iLine +1
    CALL IntToStr4(iLine,ICode)
    M%Num= TRIM(M%Num)//"_"//TRIM(ICode)
    !
    IF(N==1) THEN
      ALLOCATE(LisMinThr)
      NULLIFY(LisMinThr%next)
      LisMinThr%Value=M
      LisCur=> LisMinThr
    ELSE
      ALLOCATE(LisCur%next)
      NULLIFY(LisCur%next%next)
      LisCur%next%Value=M
      LisCur=> LisCur%next
    ENDIF
    !
    IF(FF>0) THEN !===========================================< trace ==
      WRITE(FF,"(A15,A1)",ADVANCE="NO") M%Name,T_
      DO iEl=1,SIZE(vEle)
        WRITE(FF,"(I5,A1)",ADVANCE="NO") vStoik(iEl),T_
      ENDDO
      ! build new formula, with fixed element order
      CALL Formula_Build(vEle,vStoik,Zsp,M%Div,sFormul)
      !!sFormul=""
      !!DO iEl=1,SIZE(vEle)
      !!  IF(vStoik(iEl)>0) sFormul=TRIM(sFormul)//TRIM(vEle(iEl)%NamEl)//""
      !!ENDDO
      WRITE(FF,"(I5,A1,A39,A1,A15,A1,A)") &
      & M%Div,T_,M%Formula,T_,M%Name,T_,TRIM(sFormul)
    ENDIF !==================================================</ trace ==
    !
  ENDSUBROUTINE Save_Record

ENDSUBROUTINE DtbMinThr_Read_Line

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

  !old! IF(iDebug==4) THEN
  !old!   DO I=1,SIZE(vStrField)
  !old!     PRINT *,vIField(I),TRIM(vStrField(I))
  !old!   ENDDO
  !old! ENDIF
  !PAUSE
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

SUBROUTINE ReadRVals5(Line,x1,x2,x3,x4,x5)
  USE M_IOTools,ONLY: LinToWrd,WrdToReal
  !
  CHARACTER(LEN=*),INTENT(INOUT):: Line
  REAL(dp),        INTENT(OUT)  ::x1,x2,x3,x4,x5
  !
  CHARACTER(255)   Word
  INTEGER  ::i
  LOGICAL  ::EoL
  REAL(dp),DIMENSION(1:5)::vX
  !
  vX(:)= Zero
  EoL=  .TRUE.
  DO i=1,5
    CALL LinToWrd(Line,Word,EoL)
    CALL WrdToReal(Word,vX(i))
    IF(EoL) EXIT
  ENDDO
  !
  x1=vX(1)
  x2=vX(2)
  x3=vX(3)
  x4=vX(4)
  x5=vX(5)
  !
ENDSUBROUTINE ReadRVals5

SUBROUTINE DtbMinThr_ToLine(f)
!--
!-- write a multiline theriak-compatible datafile to a monoline equivalent
!--
  USE M_Files,      ONLY: DirDtbLog,Files_Index_Write
  USE M_Dtb_Vars,   ONLY: DtbFormat
  !
  INTEGER,INTENT(IN):: f
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL :: EoL,bGas,bMin
  INTEGER :: ios,I
  TYPE(T_DtbMinThr) :: M
  !
  IF(iDebug>0) WRITE(fTrc,"(/,A)") "< DtbMinThr_ToLine"
  !
  WRITE(fTrc,'(A)') &
  & "stoikio in "//TRIM(DirDtbLog)//"min_thr_line.log"
  IF(fLin==0) THEN
    CALL GetUnit(fLin)
    OPEN(fLin,FILE=TRIM(DirDtbLog)//"min_thr_line.log")
    CALL Files_Index_Write(fHtm,&
    & TRIM(DirDtbLog)//"min_thr_line.log",&
    & "multiline DATAfile from therak converted to monoline")
  ENDIF
  !
  FilCode=TRIM(DtbFormat) !-> default value
  M%Name="Z"
  !
  DoFile: DO

    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=="!")   CYCLE DoFile
    CALL AppendToEnd(L,W,EoL)

    SELECT CASE(W)
      !
      CASE("CODE")
        CALL LinToWrd(L,W,EoL)
        FilCode=TRIM(W)
        CYCLE DoFile

      CASE("ENDINPUT")
        EXIT DoFile

      CASE("ENDSPECIES")
        EXIT DoFile

      CASE("MINERAL","GAS")

        SELECT CASE(W)
          CASE("MINERAL"); bGas=.FALSE.; bMin=.TRUE.
          CASE("GAS");     bGas=.TRUE.;  bMin=.FALSE.
        END SELECT

        DoReadMin: DO
          !
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          IF(W(1:1)=="!") CYCLE DoReadMin !skip comment lines
          CALL AppendToEnd(L,W,EoL)
          !
          SELECT CASE(W)
            CASE("ENDINPUT"); EXIT DoFile
            CASE("ENDSPECIES"); EXIT DoFile
            CASE("END","ENDMINERAL","ENDGAS"); CYCLE DoFile !EXIT DoReadMin
          END SELECT

          IF(W(1:1)/="&") THEN
          !IF first char is not "&", the line may contain name, formula, etc
            WRITE(fLin,*)
            M%Name=TRIM(W)
            IF(bGas) WRITE(fLin,'(A,A1)',ADVANCE="NO") "GAS",T_
            IF(bMin) WRITE(fLin,'(A,A1)',ADVANCE="NO") "MIN",T_
            CALL LinToWrd(L,W,EoL); M%Formula=TRIM(W)
            CALL LinToWrd(L,W,EoL) !; M%Abbr=TRIM(W)
            WRITE(fLin,'(3(A,A1))',ADVANCE="NO") &
            & TRIM(FilCode),T_,M%Name,T_,M%Formula,T_
            !
            CYCLE DoReadMin
          ENDIF

          IF(W(1:1)=="&") THEN
            ! if first char is "&", 
            ! then the line is the continuation of preceding, i.e. contains data
            DO I=1,7
              CALL LinToWrd(L,W,EoL)
              IF(.NOT. EoL) THEN
                WRITE(fLin,'(A,A1)',ADVANCE="NO") TRIM(W),T_
              ELSE
                WRITE(fLin,'(A,A1)',ADVANCE="NO") "_",T_
              ENDIF
            ENDDO
          ENDIF !(W(1:1)=="&")

        ENDDO DoReadMin
      !ENDCASE("MINERAL")
    END SELECT
  ENDDO DoFile
  !
  WRITE(fLin,*)
  !
  IF(iDebug>0) WRITE(fTrc,"(A,/)") "</ DtbMinThr_ToLine"
  !
  RETURN
ENDSUBROUTINE DtbMinThr_ToLine

ENDMODULE M_Dtb_Read_DtbMinThr
