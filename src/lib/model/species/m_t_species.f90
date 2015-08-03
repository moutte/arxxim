MODULE M_T_Species
  USE M_Kinds
  !
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC:: T_Species, T_SpcData, T_SpeciesDtb
  PUBLIC:: Species_Zero
  PUBLIC:: Species_Index
  PUBLIC:: Species_Rename
  PUBLIC:: Species_Stoikio
  PUBLIC:: Species_Stoikio_Calc
  PUBLIC:: Species_EntropyZero
  !
  PUBLIC:: nElMax
  !
  TYPE:: T_SpcData !container for storing a species' data
    REAL(dp):: & !store current values
    & Mole, &    !mole number in reference system
    !            (either in 1 kg H2O, either in the box volume, etc.)
    & LAct, &    !ln(activity)
    & LGam       !ln(activ'coeff)
  ENDTYPE T_SpcData
  !
  TYPE:: T_SpeciesThermo
  ! container for storing current thermo data of a SPECIES
  ! for future use ??
    REAL(dp):: G0rt,V0,H0,S0,Cp0,LnPhi,logK
  ENDTYPE T_SpeciesThermo
  !
  INTEGER,PARAMETER:: nElMax= 64
  !
  TYPE T_SpeciesDtb
    ! pointer to a species in a database-specific array
    CHARACTER(LEN=7) :: DtbModel !selector: "AQU_HKF","MIN_HKF","LOGKTBL","LOGKANL",etc
    INTEGER          :: Indx     !index of species in the corresponding vDtbEoS
    CHARACTER(LEN=15):: DtbTrace !retrace origin of data, herited from vDtb*%Num
  ENDTYPE T_SpeciesDtb
  !
  TYPE:: T_Species !container for storing a SPECIES' structure
    CHARACTER(LEN=3 ):: Typ   !AQU/MIN/GAS/XCH
    CHARACTER(LEN=23):: NamSp !
    !
    INTEGER:: iDtb
    INTEGER:: iDiscret
    ! iDtb/=0
    ! -> species pointing to a pure species in the thermodynamic database
    ! iDtb gives index of species in vSpcDtb
    ! iDiscret/= 0 -> species generated by discretization of a mixture
    ! iDiscret gives index of species in vDiscretParam
    !
    CHARACTER(LEN=71):: Formula
    INTEGER:: vStoikio(0:nElMax)
    INTEGER:: Div !"formula factor", i.e. COMMON denominator of all element numbers
    !
    INTEGER:: Z   !charge
    ! Note: actual charge is Z/REAL(Div) ???
    !
    REAL(dp):: &
    & WeitKg,  &  !molecule weight (formula weight), in Kg
    & AquSize, &  !for aqu.species, SIZE PARAMETER for Debye-Huckel relations
    !! & Vs,Ss,Cps,Hs,Gs, &    !solvation terms, for aqu'species, at current (T,P) 
    & G0rt,V0,H0,S0,Cp0,LnFug
    !-> values at current (T,P) !! (NOT those at ref'conditions)
    !G0rt=G/RT of FORMATION -> will produce deltaG's of ASSOCIATION/PRECIPITATION
    !caveat:
    !logK database contain, in lieu of G0, the log10K -> G0/RT= - log10K *ln10
    !
    TYPE(T_SpcData):: Dat !store values of vMolF, vLnAct, vLnGam, ..
    !
  ENDTYPE T_Species
  !
  ! the field "phase" could be added in order to deal with vMulti-solution systems, 
  ! e.g. equilibrium of aqueous phase with a gas mixture of a solid solution
  !
  ! currently, the phase is associated through the iPhase index
  !
  ! in simple cases, originally considered, of equil. between aqu. solution and pure minerals or gases,
  ! the field "typ" is sufficient
  !
  ! NB1: what is referred to as "Gibbs energies" is actually G/RT, i.e. -LnK
  ! NB2: if the database is in "logK", it is generally log10, then LnK=Ln10*logK
  ! deltaG's for sec'species and minerals are for FORMATION reactions:
  !   Prim'Species -> Second'Species
  ! -> DG0(iSc)=  vSpc(iSc)%G0 - DOT_PRODUCT(tStNu(1:nCp,iSc),vSpc(1:nCp)%G0)
  ! each time basis is changed, DG0 must be re-calculated
  ! -> thus, DG0 not stored in vSpc
  !
  LOGICAL,PARAMETER:: Use_Synonym= .TRUE.
  
CONTAINS

INTEGER FUNCTION Species_Index(Str,V) !,nSpc_)
!--
!-- position of species with %NamSp--Str in vSpc_(1:nSpc_) --
!--
  CHARACTER(LEN=*),INTENT(IN):: Str
  TYPE(T_Species), INTENT(IN):: V(:)
  !
  INTEGER::I
  !
  Species_Index=0
  IF(SIZE(V)==0) RETURN
  I=0
  DO
    I=I+1
    !IF(TRIM(Str)==TRIM(V(I)%NamSp)) THEN
    IF(Synonym(TRIM(Str),TRIM(V(I)%NamSp))) THEN
      Species_Index=I
      EXIT
    ENDIF
    IF(I==SIZE(V)) EXIT
  ENDDO
  !->> if Str not found, I=0
  !
  RETURN
END FUNCTION Species_Index

LOGICAL FUNCTION Synonym(S1,S2)
  CHARACTER(LEN=*),INTENT(IN):: S1,S2
  !
  INTEGER:: J1,J2
  CHARACTER(LEN=4):: S
  CHARACTER(LEN=23):: SS1,SS2
  !
  Synonym= (TRIM(S1)==TRIM(S2))
  IF(Synonym) RETURN
  IF(.NOT. Use_Synonym) RETURN
  !
  SS1= TRIM(S1)
  SS2= TRIM(S2)
  !
  !--- to consider ",AQ" and "(AQ)" as synonyms
  S= ",AQ"
  J1= INDEX(SS1,TRIM(S))  ;  IF(J1>1) SS1= SS1(1:J1-1)//"_AQ"
  J2= INDEX(SS2,TRIM(S))  ;  IF(J2>1) SS2= SS2(1:J2-1)//"_AQ"
  S= "(AQ)"
  J1= INDEX(SS1,TRIM(S))  ;  IF(J1>1) SS1= SS1(1:J1-1)//"_AQ"
  J2= INDEX(SS2,TRIM(S))  ;  IF(J2>1) SS2= SS2(1:J2-1)//"_AQ"
  IF(TRIM(SS1)==TRIM(SS2)) THEN
    Synonym= .TRUE.
    RETURN
  ENDIF
  !
  !--- to consider ",G" and "(G)" as synonyms
  S= ",G"
  J1= INDEX(SS1,TRIM(S))  ;  IF(J1>1) SS1= SS1(1:J1-1)//"_G"
  J2= INDEX(SS2,TRIM(S))  ;  IF(J2>1) SS2= SS2(1:J2-1)//"_G"
  S= "(G)"
  J1= INDEX(SS1,TRIM(S))  ;  IF(J1>1) SS1= SS1(1:J1-1)//"_G"
  J2= INDEX(SS2,TRIM(S))  ;  IF(J2>1) SS2= SS2(1:J2-1)//"_G"
  IF(TRIM(SS1)==TRIM(SS2)) THEN
    Synonym= .TRUE.
    RETURN
  ENDIF
  !
  RETURN
END FUNCTION Synonym

SUBROUTINE Species_Rename(sName)
  CHARACTER(LEN=*),INTENT(INOUT):: sName
  !
  INTEGER:: J
  CHARACTER(LEN=4):: S
  
  IF(Use_Synonym) RETURN
  !
  !--- to consider ",AQ" and "(AQ)" as synonyms
  S= ",AQ"
  J= INDEX(sName,TRIM(S))
  IF(J>1) sName= sName(1:J-1)//"_AQ"
  S= "(AQ)"
  J= INDEX(sName,TRIM(S))
  IF(J>1) sName= sName(1:J-1)//"_AQ"
  
  !--- to consider ",G" and "(G)" as synonyms
  S= ",G"
  J= INDEX(sName,TRIM(S))
  IF(J>1) sName= sName(1:J-1)//"_G"
  S= "(G)"
  J= INDEX(sName,TRIM(S))
  IF(J>1) sName= sName(1:J-1)//"_G"
  
END SUBROUTINE Species_Rename

SUBROUTINE SpeciesDtb_Zero(S)
  TYPE(T_SpeciesDtb),INTENT(OUT):: S
  !
  S%DtbModel=  "NONE"
  S%Indx=      0
  S%DtbTrace=  "NONE"
  !
  RETURN
END SUBROUTINE SpeciesDtb_Zero

SUBROUTINE Species_Zero(S)
  TYPE(T_Species),INTENT(OUT):: S
  !
  S%Typ=     "Z"
  S%NamSp=   "Z"
  S%Formula= "Z"
  !! S%DtbTrace="Z"
  !! S%DtbModel="Z"
  S%iDtb=     0
  S%iDiscret= 0
  S%Z=        0
  S%Div=      1
  S%WeitKg=   Zero
  S%AquSize=  3.72D0 !! default aqu'SIZE
  !
  S%vStoikio(:)= 0 != S%vStoikio(1:nEl), stoichiometry in INTEGER values
  S%vStoikio(1)= 1 != the COMMON divisor of the chemical formula
  !
  S%G0rt=  Zero
  S%V0=    One
  S%H0=    Zero
  S%S0=    Zero
  S%Cp0=   Zero
  S%LnFug= Zero
  !! S%Rho=      One
  !!S%S0Ele= Zero
  !S%Vs,S%Ss,S%Cps,S%Hs,S%Gs
  !S%G0rt,S%V0,S%H0,S%S0,S%Cp0,S%logK 
END SUBROUTINE Species_Zero

SUBROUTINE Species_Stoikio_Calc(vEle,vSpc,Ok)
  USE M_T_Element,  ONLY: T_Element,Element_Index
  !
  TYPE(T_Element), INTENT(IN)   :: vEle(:)
  TYPE(T_Species), INTENT(INOUT):: vSpc(:) !modifs in nDiv,Z,Oxy
  LOGICAL,         INTENT(OUT)  :: Ok
  !
  INTEGER:: I,ieOx
  LOGICAL:: fOk
  !
  ieOx= Element_Index("OX_",vEle)
  Ok= .true.
  !
  DO I=1,SIZE(vSpc)
    !IF(vSpc(i)%iDiscret==0) THEN
    IF(vSpc(i)%iDtb>0) THEN
      ! for discretized species,
      ! use DiscretSpecies_Stoikio_Calc
      CALL Species_Stoikio(vEle,ieOx,vSpc(I),fOk)
      IF(.not. fOk) THEN
        Ok= .false.
        !! Msg= "Stoikiometry problem in species "//TRIM(vSpc(I)%NamSp)
        RETURN
      ENDIF
    ENDIF
  ENDDO
  !
  RETURN
END SUBROUTINE Species_Stoikio_Calc

SUBROUTINE Species_Stoikio(vEle,ieOx,S,fOk)
!--
!-- build stoichiometry vector of species S according to element base vEle
!-- NEW: stoichiometry saved as S%vStoikio
!--
  USE M_T_Element,ONLY: T_Element,Formula_Read
  !
  TYPE(T_Element),INTENT(IN)::    vEle(:)
  INTEGER,        INTENT(IN)::    ieOx
  TYPE(T_Species),INTENT(INOUT):: S
  LOGICAL,        INTENT(OUT)::   fOk
  !
  INTEGER:: vStoik(1:SIZE(vEle))
  INTEGER:: ZSp,nDiv,nEl,ZExc
  !
  nEl= SIZE(vEle)
  S%vStoikio(:)= 0
  !
  !write(71,'(2A)') S%NamSp,S%Formula
  CALL Formula_Read(S%Formula,vEle,ZSp,nDiv,fOk,vStoik)
  !
  IF(.NOT. fOk) RETURN
  !
  S%vStoikio(1:nEl)= vStoik(1:nEl)
  S%vStoikio(0)=     nDiv
  S%vStoikio(nEl+1)= ZSp
  !
  S%Div= nDiv !-> "formula factor"
  S%Z=   ZSp  !-> species charge
  !
  ZExc=DOT_PRODUCT(vStoik(1:nEl),vEle(1:nEl)%Z) !-> oxidation state
  !
  IF(ieOx/=0) THEN
    S%vStoikio(ieOx)= ZSp - ZExc
    fOk=.TRUE.
  !<added 2010May
  ELSE
    fOk= (Zsp==ZExc)
  !</
  ENDIF
  !
END SUBROUTINE Species_Stoikio

REAL(dp) FUNCTION Species_EntropyZero(vEle,S)
  USE M_Dtb_Const,ONLY: S0_Hydrogen
  USE M_T_Element,ONLY: T_Element
  !
  TYPE(T_Element),INTENT(IN):: vEle(:)
  TYPE(T_Species),INTENT(IN):: S
  !
  INTEGER:: N
  !
  N= SIZE(vEle)
  !
  Species_EntropyZero= &
  & ( DOT_PRODUCT(S%vStoikio(1:N),vEle(1:N)%S0)  &
  &   - S%vStoikio(N+1) *S0_Hydrogen           ) &
  & /REAL(S%vStoikio(0))
  !
END FUNCTION Species_EntropyZero

REAL(dp) FUNCTION Species_WeitKg(vEle,S)
  USE M_T_Element,ONLY: T_Element
  TYPE(T_Element),INTENT(IN):: vEle(:)
  TYPE(T_Species),INTENT(IN):: S
  !
  INTEGER:: nEl
  !
  nEl= SIZE(vEle)
  Species_WeitKg= &
  & sum(vEle(1:nEl)%WeitKg * S%vStoikio(1:nEl)) /REAL(S%vStoikio(0))
  !
END FUNCTION Species_WeitKg

ENDMODULE M_T_Species
