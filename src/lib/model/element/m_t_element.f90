MODULE M_T_Element
  USE M_Kinds
  !
  IMPLICIT NONE
  !
  PRIVATE
  !-- note on structured DATA TYPEs USEd in the project --
  !
  !-- at the moment, all data types implemented are of fixed size, 
  !-- they do not contain fields of ALLOCATABLE or POINTER TYPE,
  !
  !-- links between objects are realized using indexes,
  !-- e.g. a component c comprises an integer index C%iEle 
  !-- that points to the element vEle(C%iEle) of the element list vEle
  !
  !-- CHARACTER strings are all of fixed size
  !-- (generally - n*8 -1, for better output on text editor using tabs of length 8 ...)
  !
  PUBLIC:: T_Element
  PUBLIC:: Element_Index
  PUBLIC:: Element_Zero
  PUBLIC:: Formula_Read
  PUBLIC:: Formula_Build
  !
  TYPE:: T_Element
  !-- container for storing a CHEMICAL ELEMENT
  !-- currently not used, because of near equivalence element - component in closed systems
  !-- should be used in future developement (element-chemistry, component-termodynamics)
    CHARACTER(LEN=3):: NamEl
    CHARACTER(LEN=3):: Redox !either "FIX", or "VAR" (for Cl, Fe, Cr, As, ...)
    !=> whether the default oxydation state of the element is fixed or variable !!
    REAL(dp)     :: &
    & WeitKg, & ! atomic weight (kg)
    & S0        ! entropy (J) at 25°C/1atm
    !             (used for conversions between thermodyn. conventions)
    INTEGER      :: Z 
    ! Z= nominal valency of element, used for computing species charge, redox state, etc.
  ENDTYPE T_Element
  TYPE(T_Element):: Element_Zero= T_Element("ZZZ","VAR",One,Zero,0)
  !
CONTAINS  
 
!! TYPE(T_Element) FUNCTION Element_Zero
!!   Element_Zero%Name=   "Z"
!!   Element_Zero%ReDOx=  "VAR"
!!   Element_Zero%WeitKg= One
!!   Element_Zero%S0=     Zero
!!   Element_Zero%Z=      0
!! ENDFUNCTION Element_Zero

INTEGER FUNCTION Element_Index(Str,V)
!--
!-- -> position of T_element with %Name--Str in V
!--
  CHARACTER(LEN=*),INTENT(IN):: Str
  TYPE(T_Element), INTENT(IN):: V(:)
  INTEGER     ::I
  !
  Element_Index=0
  IF(SIZE(V)==0) RETURN
  !
  I=0
  DO
    I=I+1 
    IF(TRIM(Str)==TRIM(V(I)%NamEl)) THEN
      Element_Index=I
      EXIT
    ENDIF
    IF(I==SIZE(V)) EXIT
  ENDDO !IF Str not found -> I=0
  !
  RETURN
ENDFUNCTION Element_Index

!---

SUBROUTINE Formula_Read( &
& sFormul,vEle, &
& Z,nDiv,fOk,vStoik)

  USE M_Formula_Parser
  
  IMPLICIT NONE
  !
  CHARACTER(LEN=*),INTENT(IN) :: sFormul !the stoichio'formula of a substance
  TYPE(T_Element), INTENT(IN) :: vEle(:) !the current element list
  !
  INTEGER,INTENT(OUT):: Z,nDiv ! the charge and the "divisor"
  LOGICAL,INTENT(OUT):: fOk    ! the formula is consistent (with elemental basis)
  INTEGER,INTENT(OUT):: vStoik(SIZE(vEle)) !the "formula vector"
  
  CALL Formula_Arxim_Read_Standard (sFormul,vEle(:)%NamEl,Z,nDiv,fOk,vStoik)
  
  RETURN
END SUBROUTINE Formula_Read

!---

SUBROUTINE Formula_Build(vEle,vStoik,Zsp,nDiv,S) !,Discret_,nDiv1,nDiv2,I_,J_,S)
  USE M_Formula_Parser
  IMPLICIT NONE
  !build a formula with a "static" element order, -> to make species sorting easier ...
  TYPE(T_Element), INTENT(IN) :: vEle(:)
  INTEGER,         INTENT(IN) :: vStoik(:)
  INTEGER,         INTENT(IN) :: Zsp,nDiv
  CHARACTER(LEN=*),INTENT(OUT):: S

  CALL Formula_Arxim_Build_Standard(vEle(:)%NamEl,vStoik,Zsp,nDiv,S)
  
END SUBROUTINE Formula_Build

ENDMODULE M_T_Element

