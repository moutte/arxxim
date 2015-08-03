 MODULE M_Formula_Arxim_Standard
   
  !===============================================================
  ! Purpose :  Arxim Formulas Utils Standard
  !---------------------------------------------------------------
  ! Arxim Formula Syntax 
  !---------------------------------------------------------------
  ! CaCO3        = Ca(3)C(1)O(3) 
  ! HCO3-        = H(1)C(1)O(3)-(1)
  ! Ca+2         = Ca(1)+(2)
  ! CO3(MgCa)0.5 = C(2)O(6)MG(1)CA(1)/(2) 
  !             or C(2)O(6)MG(1)CA(1)DIV(2)
  !===============================================================
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Formula_Arxim_Build_Standard
  PUBLIC :: Formula_Arxim_Read_Standard

  LOGICAL, PARAMETER :: DebFormula = .false.

  !---

CONTAINS
  
  
  SUBROUTINE Formula_Arxim_Read_Standard(sFormul,vNameEle,Z,nDiv,fOk,vStoik) 
  !=============================================================================
  !== READ a chemical formula sFormul
  !== to a stoikiometric vector (vStoik)/Div and the charge Z
  !== according to the element list vEle
  !==
  !==-> used to read the elemental decompositon of a substance from its formula
  !==
  !== CAVEAT
  !==   vStoik is an array of integers; nDiv the formula divisor
  !==   thus, the "REAL" stoichio coeff for element vEle(i) is vstoik(i)/Div !!!
  !==
  !=============================================================================

    USE M_Formula_Utils
    USE M_IOTools,ONLY: Str_Append, WrdToInt
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),INTENT(IN) :: sFormul       !the stoichio'formula of a substance
    CHARACTER(LEN=*),INTENT(IN) :: vNameEle(:)   !the current element name list
    INTEGER,         INTENT(OUT):: Z,nDiv        !the charge and the "divisor"
    LOGICAL,         INTENT(OUT):: fOk           !the formula is consistent (with elemental basis)
    INTEGER,         INTENT(OUT):: vStoik(:)     !the "formula vector"
    !   
    CHARACTER(LEN=80):: sRest
    CHARACTER(LEN=3) :: sElem
    CHARACTER(LEN=7) :: sCoeff
    INTEGER          :: Coeff,I,LenRest,iEl
    !-------
    fOk= .TRUE.
    
    sRest= TRIM(sFormul)//"_"
    
    IF(TRIM(sRest)=="_") THEN
      fOk= .FALSE.
      RETURN
    ENDIF
    
    vStoik= 0
    Coeff=  0
    Z=      0
    nDiv=   1
    
    DO
      IF (sRest(1:1)=='_') EXIT               !END of formula
      
      LenRest= LEN_TRIM(sRest)
      I= SCAN(sRest,'(')                     !find 1st occurrence of '('
      
      IF(I==0) EXIT                          !no '(' found --> end of string
      
      IF(I>1) THEN
      
        sElem=sRest(1:I-1)                   !ElementName
        sRest=sRest(I+1:LenRest)             !rest of string
        I=SCAN(sRest,')')
        IF(I>1) THEN
          sCoeff=sRest(1:I-1)                !Element Coeff
          sRest=sRest(I+1:LenRest)           !rest of string
        ENDIF
      
      ENDIF
       
      CALL Str_Append(sElem,3)
      CALL WrdToInt(sCoeff,Coeff)
       
      SELECT CASE(TRIM(sElem))
        CASE("___"); EXIT                    !end of formula
        CASE("DIV"); nDiv=Coeff; CYCLE
        CASE("/__"); nDiv=Coeff; CYCLE
        CASE("E__"); Z= Coeff;   CYCLE
        CASE("+__"); Z= Coeff;   CYCLE
        CASE("-__"); Z=-Coeff;   CYCLE
      ENDSELECT
       
      !IF(DebFormula) WRITE(fTrc,'(A)') TRIM(sElem)
      iEl=Name_Index(sElem,vNameEle)
       
      IF(iEl==0) THEN
      
        fOk=.FALSE.  !element not found in input element list
        RETURN
       
      ELSE
      
        vStoik(iEl)=Coeff
        ! IF(iFirst==0) iFirst=iEl !iFirst==0 checks whether iFirst already allocated
        ! IF(DebFormula) WRITE(fTrc,'(A13,A3,A1,I3,A1,I3)') &
        ! & "El iEl Coeff ",sElem,T_,iEl,T_,Coeff  !,T_,Z
        
      ENDIF
       
    ENDDO
    
    RETURN
  ENDSUBROUTINE Formula_Arxim_Read_Standard

  !---
  
  SUBROUTINE Formula_Arxim_Build_Standard(vNameEle,vStoik,Zsp,nDiv,S) !,Discret_,nDiv1,nDiv2,I_,J_,S)
    !=================================================
    ! Purpose : build a formula with a "static" element order,
    !           -> to make species sorting easier ...
    !=================================================
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: vNameEle(:)
    INTEGER,          INTENT(IN) :: vStoik(:)
    INTEGER,          INTENT(IN) :: Zsp,nDiv
    CHARACTER(LEN=*), INTENT(OUT):: S
    !
    CHARACTER(LEN=4):: Str
    INTEGER:: iEl
    S=""
    DO iEl=SIZE(vNameEle),1,-1 !-> will have H and O as last elements 
      IF(vStoik(iEl)>0) THEN
        WRITE(Str,'(I4)') vStoik(iEl) !; IF(vStoik(iEl)<10) Str(1:1)='0'
        S= TRIM(S)//vNameEle(iEl)(1:2)//"("//TRIM(ADJUSTL(Str))//")"
      ENDIF
    ENDDO
    IF(Zsp>0) THEN
      WRITE(Str,'(I4)') Zsp  ; S= TRIM(S)//"+("//TRIM(ADJUSTL(Str))//")"
    ENDIF
    IF(Zsp<0) THEN
      WRITE(Str,'(I4)') -Zsp  ; S= TRIM(S)//"-("//TRIM(ADJUSTL(Str))//")"
    ENDIF
    IF(nDiv>1) THEN
      WRITE(Str,'(I4)') nDiv  ; S= TRIM(S)//"/("//TRIM(ADJUSTL(Str))//")"
    ENDIF
  ENDSUBROUTINE Formula_Arxim_Build_Standard

END MODULE M_Formula_Arxim_Standard
