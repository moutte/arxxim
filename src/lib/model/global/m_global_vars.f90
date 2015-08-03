MODULE M_Global_Vars
!--
!-- MODULE for definition of structures, constants and global variables
!--
  USE M_Kinds
  !
  USE M_T_Element,     ONLY: T_Element
  USE M_T_Species,     ONLY: T_Species,T_SpeciesDtb
  USE M_T_MixModel,    ONLY: T_MixModel
  USE M_T_SolModel,    ONLY: T_SolModel
  USE M_T_MixPhase,    ONLY: T_MixPhase
  USE M_T_SolPhase,    ONLY: T_SolPhase
  USE M_T_Phase,       ONLY: T_Phase
  USE M_T_DiscretModel,ONLY: T_DiscretModel,T_DiscretParam
  !
  USE M_T_KinFas,   ONLY: T_KinFas
  USE M_T_Kinmodel, ONLY: T_KinModel
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  TYPE(T_Element),   DIMENSION(:),ALLOCATABLE,PUBLIC:: vEle
  TYPE(T_SpeciesDtb),DIMENSION(:),ALLOCATABLE,PUBLIC:: vSpcDtb
  TYPE(T_Species),   DIMENSION(:),ALLOCATABLE,PUBLIC:: vSpc
  !! TYPE(T_SpcData),  DIMENSION(:),ALLOCATABLE,PUBLIC:: vSpcData
  !
  TYPE(T_MixModel), DIMENSION(:),ALLOCATABLE,PUBLIC:: vMixModel
  ! vMixModel stores all mixing models available for the current run
  TYPE(T_MixPhase), DIMENSION(:),ALLOCATABLE,PUBLIC:: vMixFas
  ! vMixFas stores all phase data for the current run
  !
  TYPE(T_SolModel), DIMENSION(:),ALLOCATABLE,PUBLIC:: vSolModel
  TYPE(T_SolPhase), DIMENSION(:),ALLOCATABLE,PUBLIC:: vSolFas
  TYPE(T_SolModel),PUBLIC:: SolModel
  !
  TYPE(T_Phase),    DIMENSION(:),ALLOCATABLE,PUBLIC:: vFas
  !
  TYPE(T_DiscretParam),DIMENSION(:),ALLOCATABLE,PUBLIC:: vDiscretParam
  TYPE(T_DiscretModel),DIMENSION(:),ALLOCATABLE,PUBLIC:: vDiscretModel
  !
  REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC:: tFormula
  !1:nCp,1:nSp= Stoikio of All Species in Terms of Elements, in Columns
  !
  !!LOGICAL, DIMENSION(:),ALLOCATABLE,PUBLIC:: vLAqu,vLMin,vLGas 
  !a collection of masks, all of DIMENSION 1:nSp, to flag the nature of the species
  !
  TYPE(T_KinFas),   DIMENSION(:),ALLOCATABLE,PUBLIC:: vKinFas !-> DATA on kinetic minerals
  TYPE(T_KinModel), DIMENSION(:),ALLOCATABLE,PUBLIC:: vKinModel !-> DATA on kinetic models
  !
  INTEGER,PUBLIC:: &
  & nAq, &  !nAq= NrOfAqueousSpecies=  Prim'Species(1:nCp) U Second'Species(nCp+1:nAq)
  & nMn, &  !nMn= NumberOfMinerals !at the moment, nMn includes all pure phases, including gases
  & nGs !, &  !nGs= NumberOfGases !for future use 
  !& nCp, &  !nCp= NrOfPrimarySpecies (idem NrOfComponents) = redundant with SIZE(vCpn)
  !& nSp, &  !nSp= TotalNrOfSpecies = redundant with SIZE(vSpc)
  !& nFs, &  !number of non aqueous Multi-component phases
  !& nMk
  !
  !TYPE:: T_GlobalState
  !  !CONTAINS a "state" of the chemical space
  !  !i.e. abundances of phases,
  !  !     abundances of species and their activities in their respective phases
  !  REAL(dp), ALLOCATABLE:: vFasMole(:)
  !  REAL(dp), ALLOCATABLE:: vSpcMole(:),vSpcLnAct(:)
  !  
  !ENDTYPE T_State
  !
ENDMODULE M_Global_Vars



