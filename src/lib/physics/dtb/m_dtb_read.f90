MODULE M_Dtb_Read
!--
!-- reading databases
!--
  USE M_Kinds
  USE M_IoTools
  USE M_Trace,        ONLY: iDebug,T_,fTrc,fHtm,Stop_,Pause_, Warning_
  !
  USE M_T_DtbAquHkf,  ONLY: T_DtbAquHkf
  USE M_T_DtbMinHkf,  ONLY: T_DtbMinHkf
  USE M_T_DtbMinThr,  ONLY: T_DtbMinThr
  USE M_T_DtbLogKTbl, ONLY: T_DtbLogKTbl
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC:: Dtb_Read_Species
  PUBLIC:: Dtb_Read_Format
  !
CONTAINS

SUBROUTINE Dtb_Read_Species(vEle)
!--
!-- read species block
!-- -> IF species logk EXIT
!--    ELSE build vDtbAquHkf, vDtbMinThr, etc.
!--
  USE M_T_Element,ONLY: T_Element
  USE M_Files,    ONLY: NamFInn
  !
  USE M_Dtb_Read_Vars
  USE M_Dtb_Read_DtbMinHkf
  USE M_Dtb_Read_DtbAquHkf
  USE M_Dtb_Read_DtbMinThr
  USE M_Dtb_Read_DtbLogKTbl
  USE M_Dtb_Read_DtbLogKAnl
  !
  USE M_Dtb_Vars, ONLY: Dtb_Vars_Zero, DtbFormat
  !
  TYPE(T_Element),DIMENSION(:),INTENT(IN):: vEle
  !
  CHARACTER(LEN=255):: L,W
  LOGICAL           :: EoL
  INTEGER:: F,ios,N
  !
  IF(iDebug>0) WRITE(fTrc,"(/,A)") "< Dtb_Read_Species"
  !
  CALL Dtb_Vars_Zero
  !
  CALL Dtb_Read_Format
  !
  !! PRINT *,"DtbFormat",DtbFormat ; pause_
  !
  IF(DtbFormat=="LOGKTBL") THEN
    CALL DtbLogK_TPCond_Read(N)
    IF(N==0) CALL Stop_("TP table not found for LogKTbl DATAbase")
    !! PRINT *,"DtbLogK_Dim=",N ; pause_
  ENDIF
  !
  CALL GetUnit(F); OPEN(F,FILE=TRIM(NamFInn))
  !
  !-- dimensions and lined lists initialized outside loop
  !-- -> should allow reading several successive SPECIES blocks
  !
  nAquHkf=  0 ; NULLIFY(LisAquHkf)
  nMinHkf=  0 ; NULLIFY(LisMinHkf)
  nMinThr=  0 ; NULLIFY(LisMinThr)
  nLogKTbl= 0 ; NULLIFY(LisLogKTbl)
  nLogKAnl= 0 ; NULLIFY(LisLogKAnl)
  !
  DoFile: DO 
    
    READ(F,'(A)',IOSTAT=ios) L
    IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    
    SELECT CASE(W)
      
      CASE("SPECIES")
        
        IF(EoL) THEN
          !LogK TP.table is considered the default database format
          DtbFormat= "LOGKTBL"
        ELSE
          !
          CALL LinToWrd(L,W,EoL)
          !
          IF(TRIM(W)=="LOGK" .or. TRIM(W)=="LOGK.TABLE" ) W= "LOGKTBL"
          !
          IF(TRIM(W)=="LOGK.ANALYTIC") W= "LOGKANL"
          !
          DtbFormat= TRIM(W)
          !
        ENDIF
        
        SELECT CASE(DtbFormat) !each of them will END on "ENDSPECIES"
          
        CASE("MIN.THR")
          CALL DtbMinThr_Read(f,vEle,nMinThr)      !-> LisMinThr
        CASE("GAS.THR")
          CALL DtbMinThr_Read(f,vEle,nMinThr)      !-> LisMinThr
        CASE("HSV.THR")
          CALL DtbMinThr_Read_Line(f,vEle,nMinThr) !-> LisMinThr
        
        !---------------------------------------------------------------
        CASE("AQU.HKF")
          CALL DtbAquHKF_Read(f,vEle,nAquHkf) !-> LisAquHkf
        
        !---------------------------------------------------------------
        CASE("MIN.HKF")
          CALL DtbMinHKF_Read_Old(f,vEle,nMinHkf) !-> LisMinHkf
        CASE("GAS.HKF")
          CALL DtbMinHKF_Read_Old(f,vEle,nMinHkf) !-> LisMinHkf
        CASE("HSV.HKF")
          CALL DtbMinHKF_Read(f,vEle,nMinHkf) !-> LisMinHkf
        
        !---------------------------------------------------------------
        CASE("LOGKTBL")
          CALL DtbLogKTbl_Read(f,vEle,nLogKTbl) !-> LisLogKTbl
        CASE("LOGKANL")
          CALL DtbLogKAnl_Read(f,vEle,nLogKAnl) !-> LisLogKTbl
        
        !---------------------------------------------------------------
        CASE DEFAULT
          CALL Stop_(TRIM(DtbFormat)//" <-unknown format !!??")
            
        END SELECT
        
    END SELECT
  
  ENDDO DoFile
  CLOSE(F)
  !
  IF(fLin/=0) CLOSE(fLin)
  !
  !! IF(nAquHkf==0) CALL Stop_("Found NO Aqu.Species ...") !________STOP
  !
  IF(nAquHkf>0)  CALL DtbAquHKF_Build  !-> vDtbAquHkf
  IF(nMinHkf>0)  CALL DtbMinHKF_Build  !-> vDtbMinHkf
  IF(nMinThr>0)  CALL DtbMinThr_Build  !-> vDtbMinThr
  IF(nLogKTbl>0) CALL DtbLogKTbl_Build !-> vDtbLogKTbl
  IF(nLogKAnl>0) CALL DtbLogKAnl_Build !-> vDtbLogKAnl
  !
  !IF(iDebug>0) THEN
  !  IF(nAquHkf>0) CALL DtbAquHKF_ShoStoik(vEle)
  !  IF(nMinHkf>0) CALL DtbMinHKF_ShoStoik(vEle)
  !  IF(nMinThr>0) CALL DtbMinTHR_ShoStoik(vEle)
  !ENDIF
  !
  IF(iDebug>2) CALL Dtb_Read_Check
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ Dtb_Read_Species"
  
  RETURN
ENDSUBROUTINE Dtb_Read_Species

SUBROUTINE Dtb_Read_Format
!--
!-- scan the input file for all SPECIES block
!-- IF at least one SPECIES block is of log format,
!-- THEN DtbFormat is logK,
!-- and any SPECIES block of different format is discarded
!--

  USE M_Files,ONLY: NamFInn
  USE M_Dtb_Vars,ONLY: DtbFormat
  
  CHARACTER(LEN=255):: L,W
  LOGICAL           :: EoL
  INTEGER:: ios, f
  
  CALL GetUnit(F); OPEN(F,FILE=TRIM(NamFInn))
  !
  DoFile: DO 
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppendToEnd(L,W,EoL)
    
    SELECT CASE(W)
    !
    CASE("ENDINPUT")
      EXIT DoFile
    !
    CASE("SPECIES")
      IF(EoL) THEN
        DtbFormat= "LOGKTBL"
        ! <-LogK TP.table is considered the default DATAbase format
      ELSE
        CALL LinToWrd(L,W,EoL)
        DtbFormat= TRIM(W)
        IF(TRIM(DtbFormat)=="LOGK") DtbFormat= "LOGKTBL"
      ENDIF
      IF(DtbFormat=="LOGKTBL") THEN
        CLOSE(F)
        RETURN
      ENDIF
    !
    END SELECT
  ENDDO DoFile
  !
  CLOSE(F)
  !
  RETURN
ENDSUBROUTINE Dtb_Read_Format

SUBROUTINE Dtb_Read_Check
  USE M_Dtb_Vars
  USE M_Dtb_Write
  !
  IF(SIZE(vDtbAquHkf)>0) CALL DtbAquHkf_Write !check results DtbAquHKF_Read
  IF(SIZE(vDtbMinHkf)>0) CALL DtbMinHkf_Write !check results DtbMinHkf_Read
  IF(SIZE(vDtbMinThr)>0) CALL DtbMinThr_Write !check results DtbMinThr_Read
  !
ENDSUBROUTINE Dtb_Read_Check

ENDMODULE M_Dtb_Read

!comments from deCapitani files
!
!THE NAMES OF THE GASES CAN NOT BE CHANGED
!AS THEY ARE USED TO RECOGNIZE WHICH GAS ROUTINE SHOULD BE USED.
!ST=  code for standard STATE PROPERTIES - G, H, S, V
!   UNITS: J, J, J/DEG, J/bar
!   NOTE: G=  H - T*S  WHERE S is third law entropy,
!   not entropy of formation; i.e. it"s an apparent G of formation
!CP calculated with eqn:
!   Cp=  K0 + K1/SQRT(T) + K2/T/T + K3/T/T/T + K4/T + K5*T + K6*T*T
!C1=  CODE FOR THE FIRST LINE OF CP TERMS - K0, K1, K2, K3
!C2=  CODE FOR THE 2ND LINE OF CP TERMS - K4, K5, K6
!C3=  CODE FOR MAIER-KELLY CP EQUATION - K0, K5, K2
!   UNITS: J AND DEGREES
!CP(disorder) calculated by:
!   Cp(dis)=  D0 + D1/SQRT(T) + D2/T/T + D3/T + D4*T + D5*T*T
!   To get H, S (disorder), this eqn integrated by TMIN and TMAX
!D1=  CODE FOR THE FIRST LINE OF DISORDER TERMS - D0, D1, D2, D3
!   UNITS: J  AND DEGREES
!D2=  CODE FOR THE SECOND LINE OF DISORDER TERMS - D4, D5, TMIN, TMAX
!   UNITS: J AND DEGREES
!LAMBDA HEAT CAPACITY EQUATION IS:
!   CP(LAMBDA)=  T * (L1 + L2 * T)**2
!T1 (T3,T5 FOR MORE TRANSITIONS)=  CODE FOR THE FIRST LINE OF TRANSITION TERMS
!   T TRANS(K), REFERENCE T FOR INTEGRATION, LAMBDA CP TERM L1, LAMBDA CP TERM L2, HEAT OF TRANSITION
!   UNITS: J AND DEGREES 
!   (SEE BERMAN AND BROWN, 1985, CMP, 89, 168-183)
!T2 (T4,T6 FOR MORE TRANSITIONS)=  CODE FOR THE SECOND LINE OF TRANS TERMS 
!   DT/DP SLOPE OF TRANSITION, PARAMETER FOR HIGHER ORDER FIT OF TRANSITION T AS FUNCTION OF P, 
!   DV/DT FOR LOW T POLYMORPH, DV/DP FOR LOW T POLYMORPH
!   UNITS: JOULES, BARS, DEGREES
!VOLUME EQUATION:
!   V(P,T)/V(1,298)=  1 + V1(T-298) + V2(T-298)**2 + V3(P-1) + V4(P-1)**2
!V1=  CODE FOR VOLUME PARAMETERS - V1, V2,V3,V4
!   NOTE: V1,V2,V3 terms need to be divided by 10**5, V4 by 10**8
