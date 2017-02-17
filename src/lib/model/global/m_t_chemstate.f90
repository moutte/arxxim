module M_T_ChemState
!--
!-- describe the state of a chemical system
!--

  use M_Kinds

  use M_T_Element,     only: T_Element
  use M_T_Species,     only: T_Species,T_SpeciesDtb
  use M_T_MixModel,    only: T_MixModel
  use M_T_SolModel,    only: T_SolModel
  
  ! use M_T_DiscretModel,only: T_DiscretModel,T_DiscretParam
  
  ! use M_T_MixPhase,    only: T_MixPhase
  ! use M_T_SolPhase,    only: T_SolPhase
  ! use M_T_Phase,       only: T_Phase

  ! use M_T_KinFas,   only: T_KinFas
  ! use M_T_Kinmodel, only: T_KinModel

  implicit none
  
  private
  
  type,public:: T_ChemSpace
    
    type(T_Element),    allocatable:: vEle         (:)
    type(T_SpeciesDtb), allocatable:: vSpcDtb      (:)
    type(T_Species),    allocatable:: vSpc         (:)
    type(T_MixModel),   allocatable:: vMixModel    (:)
    type(T_SolModel),   allocatable:: vSolModel    (:)

    ! type(T_DiscretParam),allocatable:: vDiscretParam(:)
    ! type(T_DiscretModel),allocatable:: vDiscretModel(:)
      
    ! type(T_MixPhase),   allocatable:: vMixFas      (:)
    ! type(T_SolPhase),   allocatable:: vSolFas      (:)
    ! type(T_Phase),      allocatable:: vFas         (:)
    
    contains
      procedure :: New_
      procedure :: Free_
      procedure :: Set_
      procedure :: Get_
  
  end type T_ChemSpace
  
  type:: T_ChemicalState
    
    real(dp):: TdgK,Pbar
    
    real(dp),allocatable:: vMolFracSpc(:)
    real(dp),allocatable:: vMolFas(:)
  
  end type T_ChemicalState
  !
  ! real(dp),allocatable,public:: tFormula(:,:)
  ! !1:nCp,1:nSp= Stoikio of All Species in Terms of Elements, in Columns
  ! !
  ! !!logical, dimension(:),allocatable,public:: vLAqu,vLMin,vLGas 
  ! !a collection of masks, all of dimension 1:nSp, to flag the nature of the species
  ! !
  ! type(T_KinFas),   allocatable,public:: vKinFas(:) !-> data on kinetic minerals
  ! type(T_KinModel), allocatable,public:: vKinModel(:) !-> data on kinetic models
  ! !
  ! integer,public:: &
  ! & nAq, &  !nAq= NrOfAqueousSpecies=  Prim'Species(1:nCp) U Second'Species(nCp+1:nAq)
  ! & nMn, &  !nMn= NumberOfMinerals !at the moment, nMn includes all pure phases, including gases
  ! & nGs !, &  !nGs= NumberOfGases !for future use 
  
  !& nCp, &  !nCp= NrOfPrimarySpecies (idem NrOfComponents) = redundant with size(vCpn)
  !& nSp, &  !nSp= TotalNrOfSpecies = redundant with size(vSpc)
  !& nFs, &  !number of non aqueous Multi-component phases
  !& nMk

  !type:: T_GlobalState
  !  !contains a "state" of the chemical space
  !  !i.e. abundances of phases,
  !  !     abundances of species and their activities in their respective phases
  !  real(dp), allocatable:: vFasMole(:)
  !  real(dp), allocatable:: vSpcMole(:),vSpcLnAct(:)
  !  
  !end type T_GlobalState
  !
contains

subroutine New_( & !
& self,          & !
& nEle, nSpcDtb,nSpc,nMixModel,nSolModel &
! ,nMixFas,nSolFas,nFas    & !
& )
  !---------------------------------------------------------------------
  class(T_ChemSpace),intent(inout):: Self
  integer:: nEle, nSpcDtb,nSpc,nMixModel,nSolModel
  ! integer:: nMixFas,nSolFas,nFas
  !---------------------------------------------------------------------
  ! self%nEle      = nEle     
  ! self%nSpcDtb   = nSpcDtb  
  ! self%nSpc      = nSpc     
  ! self%nMixModel = nMixModel
  ! self%nSolModel = nSolModel
  ! self%nMixFas   = nMixFas  
  ! self%nSolFas   = nSolFas  
  ! self%nFas      = nFas     
  
  ! nDiscretParam
  ! nDiscretModel
  
  allocate (self%vEle         (nEle      ))
  allocate (self%vSpcDtb      (nSpcDtb   ))
  allocate (self%vSpc         (nSpc      ))
  allocate (self%vMixModel    (nMixModel ))
  allocate (self%vSolModel    (nSolModel ))
  ! allocate (self%vMixFas      (nMixFas   ))
  ! allocate (self%vSolFas      (nSolFas   ))
  ! allocate (self%vFas         (nFas      ))

  ! allocate (vDiscretParam(:))
  ! allocate (vDiscretModel(:))
  
  return
end subroutine New_

subroutine Free_(self)
  !---------------------------------------------------------------------
  class(T_ChemSpace),intent(inout):: Self
  !---------------------------------------------------------------------
  ! self%nEle      = 0
  ! self%nSpcDtb   = 0
  ! self%nSpc      = 0
  ! self%nMixModel = 0
  ! self%nSolModel = 0
  ! self%nMixFas   = 0
  ! self%nSolFas   = 0
  ! self%nFas      = 0
  ! nDiscretParam
  ! nDiscretModel
  !
  deallocate (self%vEle     )
  deallocate (self%vSpcDtb  )
  deallocate (self%vSpc     )
  deallocate (self%vMixModel)
  deallocate (self%vSolModel)

  ! deallocate (self%vMixFas  )
  ! deallocate (self%vSolFas  )
  ! deallocate (self%vFas     )
  !
  ! allocate (vDiscretParam(:))
  ! allocate (vDiscretModel(:))
  
  return
end subroutine Free_

subroutine Set_( & !
& self,          & !
& vEle     ,     & !
& vSpcDtb  ,     & !
& vSpc     ,     & !
& vMixModel,     & !
& vSolModel      & !
! & vMixFas  ,     & !
! & vSolFas  ,     & !
! & vFas           & !
& )
  !---------------------------------------------------------------------
  class(T_ChemSpace),intent(inout):: Self
  type(T_Element),   intent(in):: vEle     (:)
  type(T_SpeciesDtb),intent(in):: vSpcDtb  (:)
  type(T_Species),   intent(in):: vSpc     (:)
  type(T_MixModel),  intent(in):: vMixModel(:)
  type(T_SolModel),  intent(in):: vSolModel(:)
  ! type(T_MixPhase),  intent(in):: vMixFas  (:)
  ! type(T_SolPhase),  intent(in):: vSolFas  (:)
  ! type(T_Phase),     intent(in):: vFas     (:)
  !---------------------------------------------------------------------
  self%vEle     (:)= vEle     (:)
  self%vSpcDtb  (:)= vSpcDtb  (:)
  self%vSpc     (:)= vSpc     (:)
  self%vMixModel(:)= vMixModel(:)
  self%vSolModel(:)= vSolModel(:)
  ! self%vMixFas  (:)= vMixFas  (:)
  ! self%vSolFas  (:)= vSolFas  (:)
  ! self%vFas     (:)= vFas     (:)
  
  return
end subroutine Set_
 
subroutine Get_( & !
& self     ,     & ! 
& vEle     ,     & !
& vSpcDtb  ,     & !
& vSpc     ,     & !
& vMixModel,     & !
& vSolModel      & !
! & vMixFas  ,     & !
! & vSolFas  ,     & !
! & vFas,          & !
& )
  !---------------------------------------------------------------------
  class(T_ChemSpace),intent(in):: Self
  type(T_Element),   intent(out):: vEle     (:)
  type(T_SpeciesDtb),intent(out):: vSpcDtb  (:)
  type(T_Species),   intent(out):: vSpc     (:)
  type(T_MixModel),  intent(out):: vMixModel(:)
  type(T_SolModel),  intent(out):: vSolModel(:)
  ! type(T_MixPhase),  intent(out):: vMixFas  (:)
  ! type(T_SolPhase),  intent(out):: vSolFas  (:)
  ! type(T_Phase),     intent(out):: vFas     (:)
  !---------------------------------------------------------------------
  
  vEle     (:)= self%vEle     (:)
  vSpcDtb  (:)= self%vSpcDtb  (:)
  vSpc     (:)= self%vSpc     (:)
  vMixModel(:)= self%vMixModel(:)
  vSolModel(:)= self%vSolModel(:)
  ! vMixFas  (:)= self%vMixFas  (:)
  ! vSolFas  (:)= self%vSolFas  (:)
  ! vFas     (:)= self%vFas     (:)
  
  return
end subroutine Get_

end module M_T_ChemState


