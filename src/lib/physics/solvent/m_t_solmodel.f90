MODULE M_T_SolModel
!--implementation of asymmetric mixing models,
!--water-like, with a solvent and solutes, 
!--and molality-based concentrations
!--(as opposed to symmetric mixing models, implemented with a T_MixModel)

  USE M_Kinds
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC:: T_SolModel
  PUBLIC:: T_SolModelDat
  PUBLIC:: vSolModelAct
  
  PUBLIC:: SolModel_Spc_Init
  
  TYPE:: T_SolModelDat
    ! container for current values of some properties of a solution phase
    ! density, dielectric, parameters for the activity model
    ! here, Rho is the density of pure solvent, not the solution !!
    REAL(dp):: Rho, Eps, dhA, dhB, Bdot
  ENDTYPE T_SolModelDat
  
  TYPE:: T_SolModel  
    !--- describes an assymetric mixing model (SOLUTION model),
    !--- with a solvent and solutes
    !--- (for symmetric mixing models (MIXTURE model), use a T_MixModel)
    !
    ! CAVEAT PROGRAMMATOR:
    ! whereas T_MixModel is currently a fixed size structure,
    ! with fixed maximum numbers of end members, margules, etc
    ! T_SolModel is a 'dynamic' structure, with the number of solute 
    ! species determined at run time
    !
    CHARACTER(LEN=23):: Name  !
    CHARACTER(LEN=3) :: Typ   !"LIQ"|"MIN"|"GAS"|...
    !
    INTEGER:: iActModel
    != index of the activity model in vSolModelAct
    !
    INTEGER :: iSolvent
    != index of the solvent species in current species list
    REAL(dp):: MolWeitSv
    !
    INTEGER:: nSpecies
    != number of species, including solvent = dim' of vISpecies
    INTEGER,POINTER:: vISpecies(:)
    != indices of species involved,
    ! vISpecies(i) is index of species i in current species list
    ! (= vSpc in M_Global_Vars)
    
    LOGICAL,POINTER:: vIsNonPolar(:)
    ! for EQ3 option: apply activ'coeff' on neutral species
    ! cf EQ36, v7.0, user's manual (Wolery, 1992)
    ! activ'coeff' of CO2(aq) in a (otherwise pure) NaCl solution
    ! as a function of ionic strength
    ! uses an expression proposed by Drummond, 1981
    ! this should be applied ONLY to non-polar species
    ! (O2(aq), H2(aq), N2(aq), ...),
    ! for which a salting out effect is expected
    !
    TYPE(T_SolModelDat):: Dat
    ! current values of some properties of the solution phase,
    ! density, dielectric, DH parameters ...
    
  ENDTYPE T_SolModel
  
  ! vSolModelAct= 
  ! list of name codes for the different activity models
  ! implemented for assymetric (molality-based) solutions
  INTEGER,PARAMETER:: nSolModelAct= 12
  CHARACTER(7),PARAMETER:: vSolModelAct(1:nSolModelAct) = &
  & (/ &
  & "IDEAL  ",   & ! 1
  & "DH1    ",   & ! 2
  & "DH1EQ3 ",   & ! 3
  & "DH2    ",   & ! 4
  & "DH2EQ3 ",   & ! 5
  & "DAV_1  ",   & ! 6
  & "DAV_2  ",   & ! 7
  & "PITZER ",   & ! 8
  & "SAMSON ",   & ! 9
  & "HKF81  ",   & ! 10
  & "SIT    ",   & ! 11
  & "name12 "    & ! 12
  & /)
  
CONTAINS

SUBROUTINE SolModel_Spc_Init(Str,vSpc,S,Ok,OkMsg)
!-- initialize S%iSolvent, S%nSpecies, S%vISpecies, S%MolWeitSv
!-- according to a species database vSpc,
!-- and given the name, Str, of the solvent species

  USE M_Trace
  USE M_T_Species,ONLY: T_Species

  CHARACTER(*),    INTENT(IN)   :: Str !name of solvent species
  TYPE(T_Species), INTENT(IN)   :: vSpc(:)
  TYPE(T_SolModel),INTENT(OUT)  :: S
  LOGICAL,         INTENT(OUT)  :: Ok
  CHARACTER(*),    INTENT(OUT)  :: OkMsg

  INTEGER:: N,I

  IF(iDebug>0) WRITE(fTrc,'(/,A)') "< SolModel_Spc_Init"

  Ok= .TRUE.
  OkMsg= ""
  S%iSolvent= 0

  N= 0
  DO I=1,SIZE(vSpc)
    IF(TRIM(Str)==TRIM(vSpc(I)%NamSp)) THEN
      ! solvent
      S%iSolvent=  I
      S%MolWeitSv= vSpc(S%iSolvent)%WeitKg
    ELSE
      ! solutes
      IF(vSpc(I)%Typ=="AQU") N= N+1
    ENDIF
  ENDDO

  IF(N==0) THEN
    Ok= .FALSE.
    OkMsg=                              "!! Solute species not found !!"
    RETURN !------------------------------------------------------------
  ENDIF
  !
  IF(S%iSolvent==0) THEN
    Ok= .FALSE.
    OkMsg=                             "!! Solvent species not found !!"
    RETURN !------------------------------------------------------------
  ENDIF
  !
  S%nSpecies= N
  ALLOCATE(S%vISpecies(N))
  !
  N= 0
  DO I=1,SIZE(vSpc)
    IF(vSpc(I)%Typ=="AQU" .AND. TRIM(Str)/=TRIM(vSpc(I)%NamSp)) THEN
      N= N+1
      S%vISpecies(N)= I
    ENDIF
  ENDDO
  !
  ! print *,"<SolModel_Spc_Init"
  ! print *,S%iSolvent,vSpc(S%iSolvent)%Name
  ! print *,S%vISpecies(1),"...",S%vISpecies(S%nSpecies)
  ! print *,"</SolModel_Spc_Init"
  ! pause
  !
  IF(iDebug>0) WRITE(fTrc,'(A,/)') "</ SolModel_Spc_Init"
  !
END SUBROUTINE SolModel_Spc_Init

ENDMODULE M_T_SolModel
