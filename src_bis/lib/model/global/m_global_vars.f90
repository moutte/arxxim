module M_Global_Vars
!--
!-- module for definition of structures, constants and global variables
!--
  use M_Kinds
  !
  use M_T_Element,     only: T_Element
  use M_T_Species,     only: T_Species,T_SpeciesDtb
  use M_T_MixModel,    only: T_MixModel
  use M_T_SolModel,    only: T_SolModel
  use M_T_MixPhase,    only: T_MixPhase
  use M_T_SolPhase,    only: T_SolPhase
  use M_T_Phase,       only: T_Phase
  use M_T_DiscretModel,only: T_DiscretModel,T_DiscretParam
  !
  use M_T_KinFas,   only: T_KinFas
  use M_T_Kinmodel, only: T_KinModel
  !
  implicit none
  !
  private
  !
  type(T_Element),   dimension(:),allocatable,public:: vEle
  type(T_SpeciesDtb),dimension(:),allocatable,public:: vSpcDtb
  type(T_Species),   dimension(:),allocatable,public:: vSpc
  !! type(T_SpcData),  dimension(:),allocatable,public:: vSpcData
  !
  type(T_MixModel), dimension(:),allocatable,public:: vMixModel
  ! vMixModel stores all mixing models available for the current run
  type(T_MixPhase), dimension(:),allocatable,public:: vMixFas
  ! vMixFas stores all phase data for the current run
  !
  type(T_SolModel), dimension(:),allocatable,public:: vSolModel
  type(T_SolPhase), dimension(:),allocatable,public:: vSolFas
  type(T_SolModel),public:: SolModel
  !
  type(T_Phase),    dimension(:),allocatable,public:: vFas
  !
  type(T_DiscretParam),dimension(:),allocatable,public:: vDiscretParam
  type(T_DiscretModel),dimension(:),allocatable,public:: vDiscretModel
  !
  real(dp),dimension(:,:),allocatable,public:: tFormula
  !1:nCp,1:nSp= Stoikio of All Species in Terms of Elements, in Columns
  !
  !!logical, dimension(:),allocatable,public:: vLAqu,vLMin,vLGas 
  !a collection of masks, all of dimension 1:nSp, to flag the nature of the species
  !
  type(T_KinFas),   dimension(:),allocatable,public:: vKinFas !-> data on kinetic minerals
  type(T_KinModel), dimension(:),allocatable,public:: vKinModel !-> data on kinetic models
  !
  integer,public:: &
  & nAq, &  !nAq= NrOfAqueousSpecies=  Prim'Species(1:nCp) U Second'Species(nCp+1:nAq)
  & nMn, &  !nMn= NumberOfMinerals !at the moment, nMn includes all pure phases, including gases
  & nGs !, &  !nGs= NumberOfGases !for future use 
  !& nCp, &  !nCp= NrOfPrimarySpecies (idem NrOfComponents) = redundant with size(vCpn)
  !& nSp, &  !nSp= TotalNrOfSpecies = redundant with size(vSpc)
  !& nFs, &  !number of non aqueous Multi-component phases
  !& nMk
  !
  !type:: T_GlobalState
  !  !contains a "state" of the chemical space
  !  !i.e. abundances of phases,
  !  !     abundances of species and their activities in their respective phases
  !  real(dp), allocatable:: vFasMole(:)
  !  real(dp), allocatable:: vSpcMole(:),vSpcLnAct(:)
  !  
  !end type T_State
  !
end module M_Global_Vars



