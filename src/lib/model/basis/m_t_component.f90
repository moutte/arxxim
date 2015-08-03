MODULE M_T_Component
  USE M_Kinds
  USE M_T_Species,ONLY: nElMax
  ! 
  IMPLICIT NONE 
  PRIVATE 
  PUBLIC:: T_Component 
  PUBLIC:: Component_Zero
  PUBLIC:: Component_Index
  PUBLIC:: Component_Print
  PUBLIC:: Component_Stoikio
  PUBLIC:: CpnMolMinim
  !
  TYPE:: T_Component !container for storing a thermodynamic component
    CHARACTER(LEN=23):: NamCp
    CHARACTER(LEN=71):: Formula !-> currently used only for simplex !!!
    CHARACTER(LEN=7) :: Statut  !status of component; INERT,MOBILE,BALANCE,BUFFER,...
    INTEGER :: vStoikCp(0:nElMax)
    INTEGER :: iEle
    INTEGER :: iSpc
    REAL(dp):: Factor
    !
    CHARACTER(LEN=23):: namSol
    ! namSol= name of phase controlling this component (used for mobile cpn') 
    !   for an aqueous species or a pure species, namSol="Z" 
    !   for a non aqueous solution, %namSol is the name of the solid or gas solution phase 
    !iMix, iPol = used for multi-solution system (e.g. SS-AS)
    !-> iFas is calculated from %namSol, c%iFas=MixPhase_Index(IN=vMixFas,IN=c%namSol) 
    INTEGER :: iMix !index of "controlling" mixture phase in vMixPhase (iFas=0 for aqu. or pure species)
    INTEGER :: iPol !index of the species in the end-member list of that phase 
    ! 
    REAL(dp):: Mole   !stores the current mole number (for inert component) 
    REAL(dp):: LnAct  !stores the current activity (for mobile component) 
    ! 
  ENDTYPE T_Component
  !
  REAL(dp):: CpnMolMinim= 1.D-16
 
CONTAINS 

SUBROUTINE Component_Zero(C)
  TYPE(T_Component),INTENT(OUT):: C
  !
  C%NamCp=   "Z"
  C%Formula= "Z"
  C%Statut=  "INERT"
  C%vStoikCp(:)= 0
  C%vStoikCp(0)= 1 ! divider
  C%iEle=     0
  C%iSpc=     0
  C%Factor=   1.0_dp
  C%namSol=  "Z"
  C%iMix=     0
  C%iPol=     0
  C%Mole=     Zero
  C%LnAct=    Zero
  !
END SUBROUTINE Component_Zero
 
INTEGER FUNCTION Component_Index(Str,V)
!--
!-- finds component index according to its %NamCp
!--
  CHARACTER(LEN=*), INTENT(IN):: Str
  TYPE(T_Component),INTENT(IN):: V(:)
  !
  INTEGER::I 
  
  Component_Index=0
  IF(SIZE(V)==0) RETURN
  !
  I=0
  DO 
    I=I+1 
    IF(TRIM(Str)==TRIM(V(I)%NamCp)) THEN
      Component_Index=I  ;  EXIT
    ENDIF 
    IF(I==SIZE(V)) EXIT 
  ENDDO !IF Str not found -> I=0
  
  RETURN
END FUNCTION Component_Index
 
SUBROUTINE Component_Print(f,vEle,vSpc,C) 
  USE M_T_Element,ONLY: T_Element 
  USE M_T_Species,ONLY: T_Species 
  ! 
  INTEGER,          INTENT(IN):: f 
  TYPE(T_Element),  INTENT(IN):: vEle(:) 
  TYPE(T_Species),  INTENT(IN):: vSpc(:) 
  TYPE(T_Component),INTENT(IN):: C 
  
  WRITE(f,'(3(A,1X),2(A15,1X),2G12.3)') & 
  &  C%NamCp, vEle(C%iEle)%NamEl,C%Statut, &
  &  vSpc(C%iSpc)%NamSp, &
  &  C%NamSol, &
  &  C%Mole,C%LnAct
  
  RETURN
END SUBROUTINE Component_Print 
 
SUBROUTINE Component_Stoikio(vEle,ieOx,Cpn,fOk)
!--
!-- !!! currently used only in simplex routines !!!
!-- build stoichiometry vector of component Cpn according to element base vEle
!-- NEW: stoichiometry saved as S%vStoikCp
!--
  USE M_T_Element,ONLY: T_Element,Formula_Read
  !
  TYPE(T_Element),  INTENT(IN)::    vEle(:)
  INTEGER,          INTENT(IN)::    ieOx
  TYPE(T_Component),INTENT(INOUT):: Cpn
  LOGICAL,          INTENT(OUT)::   fOk
  !
  INTEGER:: vStoik(1:SIZE(vEle))
  INTEGER:: ZSp,nDiv,nEl,ZExc
  !
  ! print *,"=Component_Stoikio=",TRIM(Cpn%Formula)  ;  pause
  !
  nEl=  SIZE(vEle)
  !
  CALL Formula_Read(  &
  & Cpn%Formula,vEle, &
  & ZSp,nDiv,fOk,vStoik)
  !
  Cpn%vStoikCp(:)=     0
  Cpn%vStoikCp(1:nEl)= vStoik(1:nEl)
  Cpn%vStoikCp(0)=     nDiv
  Cpn%vStoikCp(nEl+1)= ZSp  !-> component charge
  !
  RETURN
  !
  ZExc=DOT_PRODUCT(vStoik(1:nEl),vEle(1:nEl)%Z) !-> oxidation state
  !
  IF(ieOx/=0) THEN
    Cpn%vStoikCp(ieOx)= ZSp - ZExc
    fOk=.TRUE.
  ENDIF
  !
  RETURN
END SUBROUTINE Component_Stoikio

END MODULE M_T_Component
 
