MODULE M_Element_Read
  USE M_Kinds
  USE M_Trace,    ONLY: iDebug,T_,fTrc,Stop_,Warning_
  USE M_T_Element,ONLY: T_Element
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC:: Elements_BuildLnk
  PUBLIC:: Elements_LnkToVec
  PUBLIC:: LnkEle_Build
  PUBLIC:: T_LnkEle
  PUBLIC:: Element_Read_Redox
  PUBLIC:: Element_Read_Entropy
  !
  TYPE:: T_LnkEle
    TYPE(T_Element)         ::Value
    TYPE(T_LnkEle), POINTER ::Next
  ENDTYPE T_LnkEle
  !
CONTAINS

SUBROUTINE Element_Read_Entropy(vEle)
!--
!-- read entropy of elements at standard state
!--
!-- entropy data for elements
!-- are needed for conversion from Berman-Brown
!-- to Benson-Helgeson convention
!--
  USE M_IOTools !, ONLY:dimV,LinToWrd,GetUnit
  USE M_Files,    ONLY: NamFInn
  USE M_T_Element,ONLY: T_Element,Element_Index
  !
  TYPE(T_Element),INTENT(INOUT):: vEle(:)
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL:: EoL
  INTEGER:: F,ios,i
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "< Elements_Read_Entropy"
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFInn))
  !
  DoFile: DO 
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppENDToEND(L,W,EoL)
    SELECT CASE(W)
      CASE("ELEMENTS.ENTROPY") !!!!!!canevas for READing one "block"
        DoBlock: DO
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          CALL AppendToEND(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoBlock !skip comment lines
          SELECT CASE(W)
            CASE("END","ENDELEMENTS.ENTROPY"); EXIT DoFile
          ENDSELECT
          CALL Str_AppEND(W,3)
          i= Element_Index(W(1:3),vEle)
          IF(i>0) THEN
            CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,vEle(i)%S0)
          ELSE
            CALL Warning_("Element "//TRIM(W)//" not in base")
          ENDIF
        ENDDO DoBlock
    ENDSELECT
  ENDDO DoFile
  CLOSE(F)
  !
  IF(iDebug>0) THEN
    WRITE(fTrc,'(A,/)') "entropy of element at standard conditions:"
    DO i=1,SIZE(vEle)
      WRITE(fTrc,"(A3,A1,G12.3)") &
      & vEle(I)%NamEl, T_, vEle(I)%S0
    ENDDO
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "</ Elements_Read_Entropy"
ENDSUBROUTINE Element_Read_Entropy
  
SUBROUTINE Element_Read_Redox(vEle)
!--
!-- read redox state of (some) elements
!--
!-- format of block:
!--   ELEMENTS.REDOX
!--     FE 2
!--     ../..
!--   END
!--
  USE M_IOTools !, ONLY:dimV,LinToWrd,GetUnit
  USE M_Files,    ONLY: NamFInn
  USE M_T_Element,ONLY: T_Element,Element_Index
  !
  TYPE(T_Element),INTENT(INOUT):: vEle(:)
  !
  CHARACTER(LEN=512):: L,W
  LOGICAL:: EoL
  INTEGER:: F,ios,i
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "< Elements_Read_Redox"
  !
  CALL GetUnit(F)
  OPEN(F,FILE=TRIM(NamFInn))
  !
  DoFile: DO 
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=='!') CYCLE DoFile !skip comment lines
    CALL AppENDToEND(L,W,EoL)
    SELECT CASE(W)
      CASE("ELEMENTS.REDOX") !!!!!!canevas for READing one "block"
        DoBlock: DO
          !
          READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
          CALL LinToWrd(L,W,EoL)
          CALL AppendToEND(L,W,EoL)
          IF(W(1:1)=='!') CYCLE DoBlock !skip comment lines
          SELECT CASE(W)
            CASE("END","ENDELEMENTS.REDOX"); EXIT DoFile
          ENDSELECT
          !
          CALL Str_Append(W,3)
          i= Element_Index(W(1:3),vEle)
          IF(i>0) THEN
            CALL LinToWrd(L,W,EoL)
            CALL WrdToInt(W,vEle(i)%Z)
          ELSE
            CALL Warning_("Element "//TRIM(W)//" not in base")
          ENDIF
          !
        ENDDO DoBlock
    ENDSELECT
  ENDDO DoFile
  CLOSE(F)
  !
  IF(iDebug>0) THEN
    WRITE(fTrc,'(A,/)') "oxydation states, default values:"
    DO i=1,SIZE(vEle)
      WRITE(fTrc,"(A3,A1,I3)") &
      & vEle(I)%NamEl, T_, vEle(I)%Z
    ENDDO
  ENDIF
  !
  IF(iDebug>0) WRITE(fTrc,'(/,A,/)') "</ Elements_Read_Redox"
ENDSUBROUTINE Element_Read_ReDOx
  
SUBROUTINE LnkEle_Build(B,E,L,P)
  LOGICAL               ::b
  TYPE(T_Element)     ::E
  TYPE(T_LnkEle),POINTER::L,P
  IF(B) NULLIFY(L)
  IF(B) THEN; ALLOCATE(L);     NULLIFY(L%next);      L%Value=E;       P=>L
  ELSE;       ALLOCATE(P%next);NULLIFY(P%next%next); P%next%Value=E;  P=>P%next
  ENDIF
ENDSUBROUTINE LnkEle_Build

SUBROUTINE Elements_LnkToVec(LnkEle,vEle)
  USE M_T_Element,ONLY: T_Element
  !
  TYPE(T_LnkEle), POINTER    :: LnkEle
  TYPE(T_Element),INTENT(OUT):: vEle(:)
  !
  TYPE(T_LnkEle),POINTER:: pCur, pPrev
  INTEGER:: I
  !
  I=0
  pCur=>LnkEle
  DO WHILE (ASSOCIATED(pCur))
    I= I+1
    vEle(I)=pCur%Value
    pPrev=>pCur; pCur=> pCur%next; DEALLOCATE(pPrev)
  ENDDO
  !
ENDSUBROUTINE Elements_LnkToVec

SUBROUTINE Elements_BuildLnk(LnkEle,nEle)
  USE M_IOTools !,ONLY: LinToWrd
  USE M_Files,    ONLY: NamFEle
  USE M_T_Element,ONLY: T_Element
  !
  TYPE(T_LnkEle),POINTER    :: LnkEle
  INTEGER,       INTENT(OUT):: nEle
  !
  LOGICAL:: L_Aqueous=.TRUE.
  !
  TYPE(T_Element)   :: Ele
  CHARACTER(LEN=512):: L,W,sListElem
  LOGICAL           :: EoL
  INTEGER           :: f,ios
  TYPE(T_LnkEle),POINTER::pCur
  !
  IF(iDebug>0) WRITE(fTrc,'(A)') "< Elements_BuildLnk"
  !
  ! nEle-> total nr in database, not the nCp of the run
  !
  CALL GetUnit(f)
  OPEN(f,FILE=TRIM(NamFele))
  !
  !-------------------------------------------------------- build linked list --
  DoFile: DO 
    READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
    CALL LinToWrd(L,W,EoL)
    IF(W(1:1)=="!") CYCLE DoFile !skip comment lines
    CALL AppENDToEnd(L,W,EoL)
    SELECT CASE(W)
    !
    CASE("ENDINPUT"); EXIT DoFile
    !
    CASE("ELEMENTS")
      nEle=0; sListElem=""
      !
      IF(L_Aqueous) THEN
        Ele%NamEl="O__"  ;  Ele%WeitKg=0.015999D0  ;  Ele%Z=-2  ;  Ele%S0=102.57D0
        nEle=nEle+1
        CALL LnkEle_Build(nEle==1,Ele,LnkEle,pCur)
        sListElem=TRIM(sListElem)//TRIM(Ele%NamEl)
        !
        Ele%NamEl="H__"  ;  Ele%WeitKg=0.001008D0  ;  Ele%Z=1   ;  Ele%S0=65.34D0
        nEle=nEle+1
        CALL LnkEle_Build(nEle==1,Ele,LnkEle,pCur)
        sListElem=TRIM(sListElem)//TRIM(Ele%NamEl)
        !
        Ele%NamEl="OX_"  ;  Ele%WeitKg=Zero        ;  Ele%Z=1   ;  Ele%S0=Zero
        nEle=nEle+1
        CALL LnkEle_Build(nEle==1,Ele,LnkEle,pCur)
        sListElem=TRIM(sListElem)//TRIM(Ele%NamEl)
        !
        !-> sListElem="O__H__OX_"
        IF(iDebug==5) PRINT '(A)',TRIM(sListElem)
      ENDIF
      !
      LoopReadElem: DO
        !
        READ(F,'(A)',IOSTAT=ios) L; IF(ios/=0) EXIT DoFile
        CALL LinToWrd(L,W,EoL)
        IF(W(1:1)=="!") CYCLE LoopReadElem !skip comment lines
        CALL AppendToEnd(L,W,EoL)
        SELECT CASE(W)
          CASE("ENDINPUT"); EXIT DoFile
          CASE("END","ENDELEMENTS"); EXIT LoopReadElem
        END SELECT
        
        CALL Str_Append(W,3)
        IF(INDEX(sListElem,TRIM(W))<1) THEN
          !
          CALL Str_Append(W,3)
          Ele%NamEl=TRIM(W)
          IF(iDebug>0) WRITE(fTrc,"(A3)") Ele%NamEl
          !
          sListElem=TRIM(sListElem)//TRIM(W)
          IF(iDebug==5) PRINT '(A)',TRIM(sListElem)
          !
          CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,Ele%WeitKg)
          !
          CALL LinToWrd(L,W,EoL); CALL WrdToInt (W,Ele%Z)
          IF(Ele%Z==0) THEN  ; Ele%ReDOx="VAR"
          ELSE               ; Ele%ReDOx="FIX"
          !
          ENDIF
          CALL LinToWrd(L,W,EoL)
          !
          CALL LinToWrd(L,W,EoL); CALL WrdToReal(W,Ele%S0)
          !
          nEle=nEle+1
          !!IF(nEle <= nElMax) THEN
          CALL LnkEle_Build(nEle==1,Ele,LnkEle,pCur)
          !!ELSE
          !!  CALL Warning_("Too many elements >> element Skipped")
          !!ENDIF
        ENDIF
        !
      ENDDO LoopReadElem
      
      IF(iDebug>0) WRITE(fTrc,"(A,I3)") "nEle=",nEle
      
    !ENDCASE("ELEMENT")
    END SELECT
  ENDDO DoFile
  CLOSE(f)
  !-------------------------------------------------------/ build linked list --
  !
  IF(nEle==0) CALL Stop_("Found NO Elements ... (missing path ??)")
  !
  IF(iDebug>0) WRITE(fTrc,'(A)') "</ Elements_BuildLnk"
ENDSUBROUTINE Elements_BuildLnk

ENDMODULE M_Element_Read
