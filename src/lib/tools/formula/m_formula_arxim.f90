MODULE M_Formula_Arxim
!===============================================================
! Purpose :  Arxim Formulas Utils 
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

  !// PUBLIC FUNCTIONs
  PUBLIC :: Formula_Arxim_Build
  PUBLIC :: Formula_Arxim_Read
  
  INTERFACE Formula_Arxim_Build
    MODULE PROCEDURE Formula_Arxim_Build_Real
    MODULE PROCEDURE Formula_Arxim_Build_Int
  END INTERFACE

  LOGICAL, PARAMETER :: DebFormula = .false.

CONTAINS

  SUBROUTINE Formula_Arxim_Build_Int ( vName, vCoef, OutFormula )
  !=====================================================================
  ! Purpose : Create an Arxim Formula from a table of Integer Coefficients
  ! Remark. Name can contain the pseudo-elements (+) and (-) 
  !=====================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: vName(:) 
    CHARACTER(LEN=*), INTENT(OUT) :: OutFormula
    !---
    INTEGER, INTENT(IN) :: vCoef(:)
    INTEGER :: n, i
    CHARACTER(LEN=20) :: scoef
    !---
    OutFormula = ''
    n = SIZE(vName)
    DO i =1,n
      IF (vCoef(i)==0) THEN 
      !// nothing to DO
      ELSE
        scoef=''
        WRITE(scoef,'(I6)') vCoef(i)
        scoef = ADJUSTL(scoef)
        OutFormula = TRIM(OutFormula)//TRIM(vName(i))//'('//TRIM(scoef)//')'
      END IF
    END DO
    !
  END SUBROUTINE Formula_Arxim_Build_Int

  !---
  
  SUBROUTINE Formula_Arxim_Build_Real ( vName, vCoef, OutFormula )
  !=====================================================================
  ! Purpose : Create an Arxim Formula from a table of Real Coefficients$
  ! Remark. Name can contain the pseudo-elements (+) and (-) 
  !=====================================================================
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: vName(:) 
    CHARACTER(LEN=*), INTENT(OUT) :: OutFormula
    !---
    REAL(kind=8), INTENT(IN) :: vCoef(:)
    INTEGER :: n, i
    CHARACTER(LEN=20) :: scoef
    REAL(kind=8) :: epsilon = 1.d-20
    !---
    OutFormula = ''
    n = SIZE(vName)
    DO i =1,n
      IF (ABS(vCoef(i))<epsilon) THEN 
      !// nothing to DO
      ELSE
        !// WRITE(*,*) TRIM(name(i)), ':', vCoef(i)
        scoef=''
        WRITE(scoef,'(F10.6)') vCoef(i)
        scoef = ADJUSTL(scoef)
        OutFormula = TRIM(OutFormula)//TRIM(vName(i))//'('//TRIM(scoef)//')'
      END IF
    END DO

  END SUBROUTINE Formula_Arxim_Build_Real

 
  !---

  SUBROUTINE Formula_Arxim_Read(sFormul,vNameEle,Z,nDiv,fOk,vStoik) 
  !=====================================================================
  !.read a chemical formula sFormul to a stoikiometric vector (vStoik)/Div 
  !.and the charge Z according to the element list vEle
  !.
  !.CAVEAT
  !.vStoik is an array of integers; nDiv the formula divisor
  !.thus, the "REAL" stoichio coeff for element vEle(i) is vstoik(i)/Div !!!
  !.
  !.-> used to read the elemental decompositon of a substance from its formula
  !=====================================================================
    USE M_Formula_Utils
    USE M_IOTools,ONLY: WrdToInt
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
    vStoik= 0
    Coeff=  0; Z=0; nDiv=1
    
    DO
      
      IF (sRest(1:1)=='_') EXIT      !end of formula
      
      LenRest=LEN_TRIM(sRest)
      I=SCAN(sRest,'(')              !find 1st occurrence of '('
      
      IF(I==0) EXIT                  !no '(' found --> END of string
      
      IF(I>1) THEN
        sElem=sRest(1:I-1)           !ElementName
        sRest=sRest(I+1:LenRest)     !rest of string
        I=SCAN(sRest,')')
        IF(I>1) THEN
          sCoeff=sRest(1:I-1)        !Element Coeff
          sRest=sRest(I+1:LenRest)   !rest of string
        ENDIF
      ENDIF
      
      !CALL Str_AppEND(sElem,3)
      CALL WrdToInt(sCoeff,Coeff)
      
      SELECT CASE(TRIM(sElem))
      CASE(""); EXIT                 !end of formula
      CASE("DIV") ; nDiv=Coeff ; CYCLE
      CASE("/")   ; nDiv=Coeff ; CYCLE
      CASE("E")   ; Z= Coeff   ; CYCLE
      CASE("+")   ; Z= Coeff   ; CYCLE
      CASE("-")   ; Z=-Coeff   ; CYCLE
      END SELECT
      
      !IF(DebFormula) WRITE(fTrc,'(A)') TRIM(sElem)
      iEl=Name_Index(sElem,vNameEle)
      IF(iEl==0) THEN
        fOk=.FALSE.
        RETURN    !element not found in input element list
      ELSE
        vStoik(iEl)=Coeff
        !IF(iFirst==0) iFirst=iEl !iFirst==0 checks whether iFirst already allocated
        !IF(DebFormula) WRITE(fTrc,'(A13,A3,A1,I3,A1,I3)') "El iEl Coeff ",sElem,T_,iEl,T_,Coeff  !,T_,Z
      ENDIF
    
    ENDDO
    
    RETURN
  ENDSUBROUTINE Formula_Arxim_Read

END MODULE M_Formula_Arxim
